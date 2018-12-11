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

from cobra.util import fix_objective_as_constraint
from cobra.exceptions import OptimizationError
from cobra.flux_analysis import find_essential_reactions

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
from cameo.visualization.plotting import plotter

logger = logging.getLogger(__name__)

__all__ = ["OptKnock"]


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
    >>> model.solver = "cplex" # Using cplex is recommended
    >>> optknock = OptKnock(model)
    >>> result = optknock.run(k=2, target="EX_ac_e", max_results=3)
    """

    def __init__(self, model, exclude_reactions=None, remove_blocked=True, fraction_of_optimum=0.1,
                 exclude_non_gene_reactions=True, use_nullspace_simplification=True, *args, **kwargs):
        super(OptKnock, self).__init__(*args, **kwargs)
        self._model = model.copy()
        self._original_model = model

        if "cplex" in config.solvers:
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
        if remove_blocked:
            self._remove_blocked_reactions()
        if not exclude_reactions:
            exclude_reactions = []
        if exclude_non_gene_reactions:
            exclude_reactions += [r for r in self._model.reactions if not r.genes]

        self._build_problem(exclude_reactions, use_nullspace_simplification)

    def _remove_blocked_reactions(self):
        fva_res = flux_variability_analysis(self._model, fraction_of_optimum=0)
        # FIXME: Iterate over the index only (reaction identifiers).
        blocked = [
            self._model.reactions.get_by_id(reaction) for reaction, row in fva_res.data_frame.iterrows()
            if (round(row["lower_bound"], config.ndecimals) == round(
                row["upper_bound"], config.ndecimals) == 0)
        ]
        self._model.remove_reactions(blocked)

    def _reduce_to_nullspace(self, reactions):
        self.reaction_groups = find_coupled_reactions_nullspace(self._model)
        reaction_groups_keys = [set(group) for group in self.reaction_groups]
        reduced_reactions = reduce_reaction_set(reactions, reaction_groups_keys)
        return reduced_reactions

    def _build_problem(self, essential_reactions, use_nullspace_simplification):
        logger.debug("Starting to formulate OptKnock problem")

        self.essential_reactions = find_essential_reactions(self._model, processes=1).union(self._model.exchanges)
        if essential_reactions:
            self.essential_reactions.update(set(get_reaction_for(self._model, r) for r in essential_reactions))

        reactions = set(self._model.reactions) - self.essential_reactions
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
            if reaction not in self.essential_reactions and reaction.lower_bound <= 0 <= reaction.upper_bound:
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

    def plot(self, index=0, grid=None, width=None, height=None, title=None, palette=None, **kwargs):
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
