# Copyright 2015 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
This module contains algorithms based on linear programming techniques, including mixed-integer linear programming
"""

from __future__ import print_function

import logging
import warnings

import numpy
from IProgress.progressbar import ProgressBar
from IProgress.widgets import Bar, Percentage
from pandas import DataFrame
from sympy import Add

import cobra
from cobra.util import fix_objective_as_constraint
from cobra.exceptions import OptimizationError
from cobra.flux_analysis import find_essential_reactions

import optlang
from optlang.duality import convert_linear_problem_to_dual

import cameo
from cameo import config
from cameo import ui
from cameo.core.model_dual import convert_to_dual
from cameo.core.strain_design import StrainDesignMethodResult, StrainDesignMethod, StrainDesign
from cameo.core.target import ReactionKnockoutTarget
from cameo.core.utils import get_reaction_for
from cameo.flux_analysis.analysis import phenotypic_phase_plane, flux_variability_analysis
from cameo.flux_analysis.simulation import fba
from cameo.flux_analysis.structural import find_coupled_reactions_nullspace
from cameo.util import reduce_reaction_set, decompose_reaction_groups

logger = logging.getLogger(__name__)

__all__ = ["OptKnock", "GrowthCouplingPotential"]


class OptKnock(StrainDesignMethod):
    """
    OptKnock.

    OptKnock solves a bi-level optimization problem, finding the set of knockouts that allows maximal
    target production under optimal growth.

    Parameters
    ----------
    model : cobra.Model
        A model to be used for finding optimal knockouts. Always set a non-zero lower bound on
        biomass reaction before using OptKnock.
    exclude_reactions : iterable of str or Reaction objects
        Reactions that will not be knocked out. Excluding reactions can give more realistic results
        and decrease running time. Essential reactions and exchanges are always excluded.
    remove_blocked : boolean (default True)
        If True, reactions that cannot carry flux (determined by FVA) will be removed from the model.
        This reduces running time significantly.
    fraction_of_optimum : If not None, this value will be used to constrain the inner objective (e.g. growth) to
        a fraction of the optimal inner objective value. If inner objective is not constrained manually
        this argument should be used. (Default: None)
    exclude_non_gene_reactions : If True (default), reactions that are not associated with genes will not be
        knocked out. This results in more practically relevant solutions as well as shorter running times.
    use_nullspace_simplification: Boolean (default True)
        Use a basis for the nullspace to find groups of reactions whose fluxes are multiples of each other. From
        each of these groups only 1 reaction will be included as a possible knockout

    Examples
    --------
    >>> from cameo import models
    >>> from cameo.strain_design.deterministic import OptKnock
    >>> model = models.bigg.e_coli_core
    >>> model.reactions.Biomass_Ecoli_core_w_GAM.lower_bound = 0.1
    >>> model.solver = "gurobi" # Using gurobi or cplex is recommended
    >>> optknock = OptKnock(model)
    >>> result = optknock.run(k=2, target="EX_ac_e", max_results=3)
    """

    def __init__(self, model, exclude_reactions=None, remove_blocked=True, fraction_of_optimum=0.1,
                 exclude_non_gene_reactions=True, use_nullspace_simplification=True, *args, **kwargs):
        super(OptKnock, self).__init__(*args, **kwargs)
        self._model = model.copy()
        self._original_model = model

        if "gurobi" in config.solvers:
            logger.info("Changing solver to Gurobi and tweaking some parameters.")
            if "gurobi_interface" not in model.solver.interface.__name__:
                model.solver = "gurobi"
            # The tolerances are set to the minimum value. This gives maximum precision.
            problem = model.solver.problem
            problem.params.NodeMethod = 1  # primal simplex node relaxation
            problem.params.FeasibilityTol = 1e-9
            problem.params.OptimalityTol = 1e-3
            problem.params.IntFeasTol = 1e-9
            problem.params.MIPgapAbs = 1e-9
            problem.params.MIPgap = 1e-9
        elif "cplex" in config.solvers:
            logger.debug("Changing solver to cplex and tweaking some parameters.")
            if "cplex_interface" not in self._model.solver.interface.__name__:
                self._model.solver = "cplex"
            problem = self._model.solver.problem
            problem.parameters.mip.strategy.startalgorithm.set(1)
            problem.parameters.simplex.tolerances.feasibility.set(1e-8)
            problem.parameters.simplex.tolerances.optimality.set(1e-8)
            problem.parameters.mip.tolerances.integrality.set(1e-8)
            problem.parameters.mip.tolerances.absmipgap.set(1e-8)
            problem.parameters.mip.tolerances.mipgap.set(1e-8)
        else:
            warnings.warn("You are trying to run OptKnock with %s. This might not end well." %
                          self._model.solver.interface.__name__.split(".")[-1])

        if fraction_of_optimum is not None:
            fix_objective_as_constraint(self._model, fraction=fraction_of_optimum)
        blocked = self._remove_blocked_reactions() if remove_blocked else []
        if exclude_reactions:
            # Convert exclude_reactions to reaction ID's
            exclude_reactions = [
                r.id if isinstance(r, cobra.Reaction) else r for r in exclude_reactions
            ]
            # if a blocked reaction were in exclude reactions, it would raise
            # because blocked reactions are removed from the model
            exclude_reactions = [r_id for r_id in exclude_reactions if r_id not in blocked]
            for r_id in exclude_reactions:
                if r_id not in self._model.reactions:
                    raise ValueError("Excluded reaction {} is not in the model".format(r_id))
        else:
            exclude_reactions = []
        if exclude_non_gene_reactions:
            exclude_reactions += [r.id for r in self._model.reactions if not r.genes]

        self._build_problem(exclude_reactions, use_nullspace_simplification)

    def _remove_blocked_reactions(self):
        fva_res = flux_variability_analysis(self._model, fraction_of_optimum=0)
        blocked = [
            reaction for reaction, row in fva_res.data_frame.iterrows()
            if (round(row["lower_bound"], config.ndecimals) == round(
                row["upper_bound"], config.ndecimals) == 0)
        ]
        self._model.remove_reactions(blocked)
        return blocked

    def _reduce_to_nullspace(self, reactions):
        self.reaction_groups = find_coupled_reactions_nullspace(self._model)
        reaction_groups_keys = [set(group) for group in self.reaction_groups]
        reduced_reactions = reduce_reaction_set(reactions, reaction_groups_keys)
        return reduced_reactions

    def _build_problem(self, exclude_reactions, use_nullspace_simplification):
        logger.debug("Starting to formulate OptKnock problem")

        self.essential_reactions = find_essential_reactions(self._model, processes=1).union(self._model.boundary)
        if exclude_reactions:
            self.exclude_reactions = set.union(
                self.essential_reactions,
                set(self._model.reactions.get_by_id(r) for r in exclude_reactions)
            )

        reactions = set(self._model.reactions) - self.exclude_reactions
        if use_nullspace_simplification:
            reactions = self._reduce_to_nullspace(reactions)
        else:
            self.reaction_groups = None

        self._make_dual()

        self._combine_primal_and_dual()
        logger.debug("Primal and dual successfully combined")

        y_vars = {}
        constrained_dual_vars = set()
        for reaction in reactions:
            if reaction not in self.exclude_reactions and reaction.lower_bound <= 0 <= reaction.upper_bound:
                y_var, constrained_vars = self._add_knockout_constraints(reaction)
                y_vars[y_var] = reaction
                constrained_dual_vars.update(constrained_vars)
        self._y_vars = y_vars

        primal_objective = self._model.solver.objective
        dual_objective = self._model.solver.interface.Objective.clone(
            self._dual_problem.objective, model=self._model.solver)

        reduced_expression = Add(*((c * v) for v, c in dual_objective.expression.as_coefficients_dict().items()
                                   if v not in constrained_dual_vars))
        dual_objective = self._model.solver.interface.Objective(reduced_expression, direction=dual_objective.direction)

        optimality_constraint = self._model.solver.interface.Constraint(
            primal_objective.expression - dual_objective.expression,
            lb=0, ub=0, name="inner_optimality")
        self._model.solver.add(optimality_constraint)
        logger.debug("Inner optimality constrained")

        logger.debug("Adding constraint for number of knockouts")
        knockout_number_constraint = self._model.solver.interface.Constraint(
            Add(*y_vars), lb=len(y_vars), ub=len(y_vars)
        )
        self._model.solver.add(knockout_number_constraint)
        self._number_of_knockouts_constraint = knockout_number_constraint

    def _make_dual(self):
        dual_problem = convert_to_dual(self._model.solver)
        self._dual_problem = dual_problem
        logger.debug("Dual problem successfully created")

    def _combine_primal_and_dual(self):
        primal_problem = self._model.solver
        dual_problem = self._dual_problem

        for var in dual_problem.variables:
            var = primal_problem.interface.Variable.clone(var)
            primal_problem.add(var)
        for const in dual_problem.constraints:
            const = primal_problem.interface.Constraint.clone(const, model=primal_problem)
            primal_problem.add(const)

    def _add_knockout_constraints(self, reaction):
        interface = self._model.solver.interface
        y_var = interface.Variable("y_" + reaction.id, type="binary")

        self._model.solver.add(interface.Constraint(reaction.flux_expression - 1000 * y_var, ub=0))
        self._model.solver.add(interface.Constraint(reaction.flux_expression + 1000 * y_var, lb=0))

        constrained_vars = []

        if reaction.upper_bound != 0:
            dual_forward_ub = self._model.solver.variables["dual_" + reaction.forward_variable.name + "_ub"]
            self._model.solver.add(interface.Constraint(dual_forward_ub - 1000 * (1 - y_var), ub=0))
            constrained_vars.append(dual_forward_ub)
        if reaction.lower_bound != 0:
            dual_reverse_ub = self._model.solver.variables["dual_" + reaction.reverse_variable.name + "_ub"]
            self._model.solver.add(interface.Constraint(dual_reverse_ub - 1000 * (1 - y_var), ub=0))
            constrained_vars.append(dual_reverse_ub)

        return y_var, constrained_vars

    def run(self, max_knockouts=5, biomass=None, target=None, max_results=1, *args, **kwargs):
        """
        Perform the OptKnock simulation

        Parameters
        ----------
        target: str, Metabolite or Reaction
            The design target
        biomass: str, Metabolite or Reaction
            The biomass definition in the model
        max_knockouts: int
            Max number of knockouts allowed
        max_results: int
            Max number of different designs to return if found

        Returns
        -------
        OptKnockResult
        """
        # TODO: why not required arguments?
        if biomass is None or target is None:
            raise ValueError('missing biomass and/or target reaction')
        target = get_reaction_for(self._model, target, add=False)
        biomass = get_reaction_for(self._model, biomass, add=False)

        knockout_list = []
        fluxes_list = []
        production_list = []
        biomass_list = []
        loader_id = ui.loading()
        with self._model:
            self._model.objective = target.id
            self._number_of_knockouts_constraint.lb = self._number_of_knockouts_constraint.ub - max_knockouts
            count = 0
            while count < max_results:
                try:
                    solution = self._model.optimize(raise_error=True)
                except OptimizationError as e:
                    logger.debug("Problem could not be solved. Terminating and returning " + str(count) + " solutions")
                    logger.debug(str(e))
                    break

                knockouts = tuple(reaction for y, reaction in self._y_vars.items() if round(y.primal, 3) == 0)
                assert len(knockouts) <= max_knockouts

                if self.reaction_groups:
                    combinations = decompose_reaction_groups(self.reaction_groups, knockouts)
                    for kos in combinations:
                        knockout_list.append({r.id for r in kos})
                        fluxes_list.append(solution.fluxes)
                        production_list.append(solution.objective_value)
                        biomass_list.append(solution.fluxes[biomass.id])
                else:
                    knockout_list.append({r.id for r in knockouts})
                    fluxes_list.append(solution.fluxes)
                    production_list.append(solution.objective_value)
                    biomass_list.append(solution.fluxes[biomass.id])

                # Add an integer cut
                y_vars_to_cut = [y for y in self._y_vars if round(y.primal, 3) == 0]
                integer_cut = self._model.solver.interface.Constraint(Add(*y_vars_to_cut),
                                                                      lb=1,
                                                                      name="integer_cut_" + str(count))

                if len(knockouts) < max_knockouts:
                    self._number_of_knockouts_constraint.lb = self._number_of_knockouts_constraint.ub - len(knockouts)
                self._model.add_cons_vars(integer_cut)
                count += 1

            ui.stop_loader(loader_id)

            return OptKnockResult(self._original_model, knockout_list, fluxes_list,
                                  production_list, biomass_list, target.id, biomass)


class RobustKnock(StrainDesignMethod):
    pass


class OptKnockResult(StrainDesignMethodResult):
    __method_name__ = "OptKnock"

    def __init__(self, model, knockouts, fluxes, production_fluxes, biomass_fluxes, target, biomass, *args, **kwargs):
        super(OptKnockResult, self).__init__(self._generate_designs(knockouts), *args, **kwargs)
        self._model = model
        self._knockouts = knockouts
        self._fluxes = fluxes
        self._production_fluxes = production_fluxes
        self._biomass_fluxes = biomass_fluxes
        self._target = target
        self._biomass = biomass
        self._processed_knockouts = None

    @staticmethod
    def _generate_designs(knockouts):
        designs = []
        for knockout_design in knockouts:
            designs.append(StrainDesign([ReactionKnockoutTarget(ko for ko in knockout_design)]))

        return designs

    def _process_knockouts(self):
        progress = ProgressBar(maxval=len(self._knockouts), widgets=["Processing solutions: ", Bar(), Percentage()])

        self._processed_knockouts = DataFrame(columns=["reactions", "size", self._target,
                                                       "biomass", "fva_min", "fva_max"])

        for i, knockouts in progress(enumerate(self._knockouts)):
            try:
                with self._model:
                    [self._model.reactions.get_by_id(ko).knock_out() for ko in knockouts]
                    fva = flux_variability_analysis(self._model, fraction_of_optimum=0.99, reactions=[self.target])
                self._processed_knockouts.loc[i] = [knockouts, len(knockouts), self.production[i], self.biomass[i],
                                                    fva.lower_bound(self.target), fva.upper_bound(self.target)]
            except OptimizationError:
                self._processed_knockouts.loc[i] = [numpy.nan for _ in self._processed_knockouts.columns]

    @property
    def knockouts(self):
        return self._knockouts

    @property
    def fluxes(self):
        return self._fluxes

    @property
    def production(self):
        return self._production_fluxes

    @property
    def biomass(self):
        return self._biomass_fluxes

    @property
    def target(self):
        return self._target

    def display_on_map(self, index=0, map_name=None, palette="YlGnBu"):
        with self._model:
            for ko in self.data_frame.loc[index, "reactions"]:
                self._model.reactions.get_by_id(ko).knock_out()
            fluxes = fba(self._model)
            fluxes.display_on_map(map_name=map_name, palette=palette)

    def plot(self, plotter, index=0, grid=None, width=None, height=None, title=None, palette=None, **kwargs):
        wt_production = phenotypic_phase_plane(self._model, objective=self._target, variables=[self._biomass.id])
        with self._model:
            for ko in self.data_frame.loc[index, "reactions"]:
                self._model.reactions.get_by_id(ko).knock_out()

            mt_production = phenotypic_phase_plane(self._model, objective=self._target, variables=[self._biomass.id])
        if title is None:
            title = "Production Envelope"

        dataframe = DataFrame(columns=["ub", "lb", "value", "strain"])
        for _, row in wt_production.iterrows():
            _df = DataFrame([[row['objective_upper_bound'], row['objective_lower_bound'], row[self._biomass.id], "WT"]],
                            columns=dataframe.columns)
            dataframe = dataframe.append(_df)
        for _, row in mt_production.iterrows():
            _df = DataFrame([[row['objective_upper_bound'], row['objective_lower_bound'], row[self._biomass.id], "MT"]],
                            columns=dataframe.columns)
            dataframe = dataframe.append(_df)

        plot = plotter.production_envelope(dataframe, grid=grid, width=width, height=height, title=title,
                                           x_axis_label=self._biomass.id, y_axis_label=self._target, palette=palette)
        plotter.display(plot)

    @property
    def data_frame(self):
        if self._processed_knockouts is None:
            self._process_knockouts()
        data_frame = DataFrame(self._processed_knockouts)
        data_frame.sort_values("size", inplace=True)
        data_frame.index = [i for i in range(len(data_frame))]
        return data_frame

    def _repr_html_(self):
        html_string = """
        <h3>OptKnock:</h3>
        <ul>
            <li>Target: %s</li>
        </ul>
        %s""" % (self._target, self.data_frame._repr_html_())
        return html_string


class GrowthCouplingPotential(StrainDesignMethod):
    """
    Optimize the Growth Coupling Potential through a combination of knockouts, knockins and medium additions

    Use the 'run' method to optimize the MILP problem and find the best solution.

    In order to find multiple sub-optimal solutions through the use of solution pools, please refer to the
    relevant solvers documentation. The underlying solver object can be accessed through
        GrowthCouplingPotential(...).model.solver.problem

    The method is described in the paper "OptCouple: Joint simulation of gene knockouts,
    insertions and medium modifications for prediction of growth-coupled strain designs"
    https://www.sciencedirect.com/science/article/pii/S2214030118300373

    Attributes
    ----------
    model : cobra.Model
        A cobra model object containing the target reaction as well as all native and heterologous reactions.
    target : cobra.Reaction or string
        The reaction (or reaction ID) that is to be growth-coupled
    knockout_reactions : list
        A list of ID's for the reactions that the algorithm should attempt to knock out
    knockin_reactions : list
        A List of ID's for the reactions that are heterologous and that the algorithm should attempt to add
    medium_additions : list
        A list of metabolite ID's for the the metabolites that should potentially be added to the medium
    n_knockouts : int
        The maximum allowed number of knockouts
    n_knockins : int
        The maximum allowed number of knockins
    n_medium : int
        The maximum allowed number of medium additions
    use_solution_pool : Bool
        Set to True when using the Gurobi solver to identify multiple alternative solutions at once
    """
    def __init__(
            self, model, target, biomass_id, knockout_reactions, knockin_reactions, medium_additions,
            n_knockouts, n_knockin, n_medium, remove_blocked=True, use_solution_pool=False, **kwargs
    ):
        super(GrowthCouplingPotential, self).__init__(**kwargs)

        # Check model for forced fluxes, which should not be present. The zero solution must be allowed.
        for r in model.reactions:
            if r.lower_bound > 0 or r.upper_bound < 0:
                raise ValueError(
                    "%s has a forced flux. Remove the reaction or change its bounds for the function to work." % r.id)

        # If a cobra.Reaction is given as target, convert to its ID string.
        if isinstance(target, cobra.Reaction):
            target = target.id
        self.target = target

        if target not in knockout_reactions:
            knockout_reactions.append(target)

        self.original_model = model
        self.biomass_id = biomass_id
        self.model = model = model.copy()

        # Add boundary reactions for these metabolites named "ADD_<metabolitename>" to the model, with bounds (-10,0).
        medium_addition_reactions = []  # List of the addition reactions added to the model
        for m in medium_additions:
            met = model.metabolites.get_by_id(m)
            r = model.add_boundary(met, type="medium addition", reaction_id=("ADD_" + m), lb=-10, ub=0)
            # r.type = "medium"
            medium_addition_reactions.append(r.id)

        # Set the solver to Gurobi for the fastest result. Set to CPLEX if Gurobi is not available.
        if "gurobi" in cobra.util.solver.solvers.keys():
            logger.info("Changing solver to Gurobi and tweaking some parameters.")
            if "gurobi_interface" not in model.solver.interface.__name__:
                model.solver = "gurobi"
            # The tolerances are set to the minimum value. This gives maximum precision.
            problem = model.solver.problem
            problem.params.NodeMethod = 1  # primal simplex node relaxation
            problem.params.FeasibilityTol = 1e-9
            problem.params.OptimalityTol = 1e-3
            problem.params.IntFeasTol = 1e-9
            problem.params.MIPgapAbs = 1e-9
            problem.params.MIPgap = 1e-9

        elif "cplex" in cobra.util.solver.solvers.keys():
            logger.warning(
                "Changing solver to CPLEX, as Gurobi is not available."
                "This may cause a big slowdown and limit options afterwards.")
            if "cplex_interface" not in model.solver.interface.__name__:
                model.solver = "cplex"
            # The tolerances are set to the minimum value. This gives maximum precision.
            problem = model.solver.problem
            problem.parameters.mip.strategy.startalgorithm.set(1)  # primal simplex node relaxation
            problem.parameters.simplex.tolerances.feasibility.set(1e-9)
            problem.parameters.simplex.tolerances.optimality.set(1e-3)
            problem.parameters.mip.tolerances.integrality.set(1e-9)
            problem.parameters.mip.tolerances.absmipgap.set(1e-9)
            problem.parameters.mip.tolerances.mipgap.set(1e-9)

        else:
            logger.warning("You are trying to run 'GrowthCouplingPotential' with %s. This might not end well." %
                           model.solver.interface.__name__.split(".")[-1])

        # Remove reactions that are blocked: no flux through these reactions possible.
        # This will reduce the search space for the solver, if not done already.
        if remove_blocked:
            blocked_reactions = cameo.flux_analysis.analysis.find_blocked_reactions(model)
            model.remove_reactions(blocked_reactions)
            blocked_reaction_ids = {r.id for r in blocked_reactions}
            knockout_reactions = {r for r in knockout_reactions if r not in blocked_reaction_ids}
            knockin_reactions = {r for r in knockin_reactions if r not in blocked_reaction_ids}
            medium_additions = {r for r in medium_additions if r not in blocked_reaction_ids}
            logger.debug("Removed " + str(len(blocked_reactions)) + " reactions that were blocked")

        # Make dual
        copied_model = model.copy()
        copied_model.optimize()
        dual_problem = convert_linear_problem_to_dual(copied_model.solver)
        logger.debug("Dual problem successfully created")

        # Combine primal and dual
        primal_problem = model.solver

        for var in dual_problem.variables:  # All variables in the dual are copied to the primal
            var = primal_problem.interface.Variable.clone(var)
            primal_problem.add(var)
        for const in dual_problem.constraints:  # All constraints in the dual are copied to the primal
            const = primal_problem.interface.Constraint.clone(const, model=primal_problem)
            primal_problem.add(const)
        logger.debug("Dual and primal combined")

        dual_problem.optimize()

        # Dictionaries to hold the binary control variables:
        native_y_vars = {}  # 1 for knockout, 0 for active
        heterologous_y_vars = {}  # 1 for knockin, 0 for inactive
        medium_y_vars = {}  # 1 for medium addition (up to -10), 0 for no addition

        # Now the fun stuff
        constrained_dual_vars = set()
        M = max(-cobra.Configuration().lower_bound, cobra.Configuration().upper_bound)  # Large value for big-M method

        # Fill the dictionaries with binary variables
        # For the knockouts: only for the native reactions which are not exchanges
        for reac_id in knockout_reactions:
            reaction = model.reactions.get_by_id(reac_id)

            # Add constraint variables
            interface = model.solver.interface
            y_var = interface.Variable("y_" + reaction.id, type="binary")

            # Constrain the primal: flux through reactions within (-1000, 1000)
            model.solver.add(interface.Constraint(reaction.flux_expression - M * (1 - y_var), ub=0,
                                                  name="primal_y_const_" + reaction.id + "_ub"))
            model.solver.add(interface.Constraint(reaction.flux_expression + M * (1 - y_var), lb=0,
                                                  name="primal_y_const_" + reaction.id + "_lb"))

            # Constrain the dual
            constrained_vars = []

            if reaction.upper_bound != 0:
                dual_forward_ub = model.solver.variables["dual_" + reaction.forward_variable.name + "_ub"]
                model.solver.add(interface.Constraint(dual_forward_ub - M * y_var, ub=0))
                constrained_vars.append(dual_forward_ub)
            if reaction.lower_bound != 0:
                dual_reverse_ub = model.solver.variables["dual_" + reaction.reverse_variable.name + "_ub"]
                model.solver.add(interface.Constraint(dual_reverse_ub - M * y_var, ub=0))
                constrained_vars.append(dual_reverse_ub)
            constrained_dual_vars.update(constrained_vars)

            # Add to dictionary with binary variables for native reactions:
            native_y_vars[y_var] = reaction

        # For the knockins and medium additions:
        for reac_id in list(knockin_reactions) + list(medium_addition_reactions):
            reaction = model.reactions.get_by_id(reac_id)
            # Add constraint variables
            interface = model.solver.interface
            y_var = interface.Variable("y_" + reaction.id, type="binary")

            # Constrain the primal: flux through reaction within (-1000, 1000)
            model.solver.add(interface.Constraint(reaction.flux_expression - M * y_var, ub=0,
                                                  name="primal_y_const_" + reaction.id + "_ub"))
            model.solver.add(interface.Constraint(reaction.flux_expression + M * y_var, lb=0,
                                                  name="primal_y_const_" + reaction.id + "_lb"))

            # Constrain the dual
            constrained_vars = []

            if reaction.upper_bound != 0:
                dual_forward_ub = model.solver.variables["dual_" + reaction.forward_variable.name + "_ub"]
                model.solver.add(interface.Constraint(dual_forward_ub - M * (1 - y_var), ub=0))
                constrained_vars.append(dual_forward_ub)
            if reaction.lower_bound != 0:
                dual_reverse_ub = model.solver.variables["dual_" + reaction.reverse_variable.name + "_ub"]
                model.solver.add(interface.Constraint(dual_reverse_ub - M * (1 - y_var), ub=0))
                constrained_vars.append(dual_reverse_ub)
            constrained_dual_vars.update(constrained_vars)

            # Add y variable to the corresponding modifications dictionary
            if reac_id in medium_addition_reactions:  # reaction.type == "medium":
                medium_y_vars[y_var] = reaction
            elif reac_id in knockin_reactions:  # reaction.type == "heterologous":
                heterologous_y_vars[y_var] = reaction
            else:
                raise RuntimeError("Something is wrong")

        logger.debug("Control variables created")

        # Set the objective
        primal_objective = model.solver.objective
        dual_objective = interface.Objective.clone(
            dual_problem.objective, model=model.solver
        )

        reduced_expression = optlang.symbolics.Add(
            *((c * v) for v, c in dual_objective.expression.as_coefficients_dict().items()
              if v not in constrained_dual_vars)
        )
        dual_objective = interface.Objective(reduced_expression, direction=dual_objective.direction)

        full_objective = interface.Objective(primal_objective.expression - dual_objective.expression, direction="max")
        model.objective = full_objective
        logger.debug("Objective created")

        # Add number of knockouts constraint. ub=K+1 as the target will be forced to be 1 as well
        knockout_number_constraint = model.solver.interface.Constraint(
            optlang.symbolics.Add(*native_y_vars), lb=0, ub=n_knockouts + 1, name="number_of_knockouts_constraint"
        )
        model.solver.add(knockout_number_constraint)

        # Add number of knockins constraint
        knockin_number_constraint = model.solver.interface.Constraint(
            optlang.symbolics.Add(*heterologous_y_vars), lb=0, ub=n_knockin, name="number_of_knockins_constraint"
        )
        model.solver.add(knockin_number_constraint)

        # Add number of medium additions constraint
        medium_additions_number_constraint = model.solver.interface.Constraint(
            optlang.symbolics.Add(*medium_y_vars), lb=0, ub=n_medium, name="number_of_medium_additions_constraint"
        )
        model.solver.add(medium_additions_number_constraint)

        logger.debug("Added constraint for number of knockouts, knockins and medium additions")

        # Force target control variable to be 1, but allow flux in primal
        model.solver.constraints["primal_y_const_" + target + "_lb"].lb = -2000
        model.solver.constraints["primal_y_const_" + target + "_ub"].ub = 2000
        model.variables["y_" + target].lb = 1

        self.model = model
        self.native_y_vars = native_y_vars
        self.heterologous_y_vars = heterologous_y_vars
        self.medium_y_vars = medium_y_vars

        self.knockout_reactions = knockout_reactions
        self.knockin_reactions = knockin_reactions
        self.medium_additions = medium_additions

        self.use_solution_pool = use_solution_pool

        self.model.solver.configuration.presolve = True

        if use_solution_pool:
            if "gurobi_interface" not in model.solver.interface.__name__:
                raise NotImplementedError("Solution pools are currently only available with the Gurobi solver")

            # Find the N best solutions (N is specified as pool_size in run())
            self.model.solver.problem.params.PoolSearchMode = 2

            # Don't store solutions with objective = 0
            self.model.solver.problem.params.PoolGap = 0.99

    def run(self, pool_size=None, pool_time_limit=None):
        """
        Run the GrowthCouplingPotential method

        pool_size: None or int
            The number of alternative solutions to find (0 - 2,000,000,000)
        pool_time_limit: None or int
            The maximum time (in seconds) to spend finding solutions. Will terminate with a time_limit status.
        """
        if self.use_solution_pool:
            if pool_size is not None:
                self.model.solver.problem.params.PoolSolutions = pool_size

            if pool_time_limit is not None:
                self.model.solver.problem.params.TimeLimit = pool_time_limit

        sol = self.model.optimize()

        if self.use_solution_pool:
            solution_pool = self.extract_gurobi_solution_pool()
            return solution_pool
        else:
            objective_value = sol.objective_value
            knockouts = [
                reac.id for var, reac in self.native_y_vars.items() if var.primal > 0.9 and reac.id != self.target
            ]
            knockins = [reac.id for var, reac in self.heterologous_y_vars.items() if var.primal > 0.9]
            medium_additions = [reac.id for var, reac in self.medium_y_vars.items() if var.primal > 0.9]
            return {
                "knockouts": knockouts,
                "knockins": knockins,
                "medium": medium_additions,
                "obj_val": objective_value
            }

    def extract_gurobi_solution_pool(self):
        """
        Extracts the individual solutions from an MILP solution pool created by the Gurobi solver.
        Returns the solution pool and a header with information to accompany the list of dicts.
        """

        target_id = self.target
        biomass_id = self.biomass_id

        milp = self.model
        if "gurobi_interface" not in milp.solver.interface.__name__:
            raise ValueError(
                "The solver for the MILP is not Gurobi. This function only works for Gurobi created solution pools.")
        if target_id not in milp.solver.problem.getVars():
            raise ValueError("Specified target reaction ID %s not found in the variable list." % target_id)
        if biomass_id not in milp.solver.problem.getVars():
            raise ValueError("Specified biomass reaction ID %s not found in the variable list." % biomass_id)

        solution_pool = []
        for solution_number in range(milp.solver.problem.solcount):
            solution = {}
            milp.solver.problem.params.SolutionNumber = solution_number
            solution["number"] = solution_number
            # Note that the objective is the growth coupling potential, not the growth or production rate
            solution["objective"] = milp.solver.problem.PoolObjVal

            # Prepare dict for looping over variables
            target_flux_name = "target_flux_(%s)" % target_id
            biomass_id_reverse = biomass_id + "_reverse_"
            target_id_reverse = target_id + "_reverse_"

            solution["growth_coupling_slope"] = None
            solution[target_flux_name] = 0
            solution["growth_rate"] = 0
            solution["growth_without_target"] = None
            solution["knockout"] = []
            solution["knockin"] = []
            solution["medium_addition"] = []

            for v in milp.solver.problem.getVars():
                if v.VarName.startswith("y_") and v.xn > 0.99:
                    reaction_id = v.VarName[2:]  # removing "y_"
                    if reaction_id == self.target:
                        continue
                    if reaction_id in self.knockout_reactions:
                        solution["knockout"].append(reaction_id)
                    elif reaction_id in self.knockin_reactions:
                        solution["knockin"].append(reaction_id)
                    elif reaction_id in self.medium_additions:
                        solution["medium_addition"].append(reaction_id)
                    else:
                        raise RuntimeError(reaction_id)

                # Extract target flux and growth rate. Subtract the forward and reverse fluxes to get the net flux.
                if v.VarName == target_id:
                    solution[target_flux_name] += v.xn
                if v.VarName.startswith(target_id_reverse):
                    solution[target_flux_name] -= v.xn
                if v.VarName == biomass_id:
                    solution["growth_rate"] += v.xn
                if v.VarName.startswith(biomass_id_reverse):
                    solution["growth_rate"] -= v.xn


            solution["growth_without_target"] = solution["growth_rate"] - solution["objective"]
            if round(solution["objective"], 7) > 0:
                solution["growth_coupling_slope"] = solution[target_flux_name] / solution["objective"]

            solution_pool.append(solution)

        return solution_pool

    @staticmethod
    def write_solutions_to_file(solutions, filename):
        import pandas as pd
        pd.DataFrame(solutions).to_csv(filename)
