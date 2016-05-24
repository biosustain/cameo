# -*- coding: utf-8 -*-
# Copyright 2013 Novo Nordisk Foundation Center for Biosustainability, DTU.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

# See the License for the specific language governing permissions and
# limitations under the License.

"""A solver-based model class and other extensions of cobrapy objects.
"""

from __future__ import absolute_import, print_function
from functools import partial

import six

import time
import csv
from copy import copy, deepcopy

import types

import cobra
import sympy
from sympy import Add
from sympy import Mul

from sympy.core.singleton import S
import optlang
from pandas import DataFrame, pandas

from cameo.util import TimeMachine, inheritdocstring
from cameo import config
from cameo import exceptions

from cameo.exceptions import SolveError, Infeasible, UndefinedSolution
from .reaction import Reaction
from .solution import LazySolution, Solution
from cameo.core.metabolite import Metabolite
from cameo.core.gene import Gene

import logging

__all__ = ['to_solver_based_model', 'SolverBasedModel']

logger = logging.getLogger(__name__)

add = Add._from_args
mul = Mul._from_args


def to_solver_based_model(cobrapy_model, solver_interface=optlang):
    """Convert a cobrapy model into a solver-based model.

    Parameters
    ----------
    cobrapy_model : cobra.core.Model
    solver_interface : solver_interface, optional
        For example, optlang.glpk_interface or any other optlang interface (the default is optlang.interface).
    """

    solver_interface = config.solvers.get(solver_interface, solver_interface)
    solver_based_model = SolverBasedModel(
        solver_interface=solver_interface, description=cobrapy_model)
    return solver_based_model


@six.add_metaclass(inheritdocstring)
class SolverBasedModel(cobra.core.Model):
    """Implements a model with an attached optlang solver instance.

    Every model manipulation is immediately reflected in the solver instance.
    """

    def __init__(self, description=None, solver_interface=optlang, **kwargs):
        super(SolverBasedModel, self).__init__(description, **kwargs)
        cleaned_reactions = cobra.core.DictList()
        for reaction in self.reactions:
            if isinstance(reaction, Reaction):
                cleaned_reactions.append(reaction)
            else:
                cleaned_reactions.append(Reaction.clone(reaction, model=self))
        self.reactions = cleaned_reactions

        cleaned_genes = cobra.core.DictList()
        for gene in self.genes:
            if isinstance(gene, Gene):
                cleaned_genes.append(gene)
            else:
                cleaned_genes.append(Gene.clone(gene, model=self))
        self.genes = cleaned_genes

        cleaned_metabolites = cobra.core.DictList()
        for metabolite in self.metabolites:
            if isinstance(metabolite, Metabolite):
                cleaned_metabolites.append(metabolite)
            else:
                cleaned_metabolites.append(Metabolite.clone(metabolite, model=self))
        self.metabolites = cleaned_metabolites

        for metabolite in self.metabolites:
            metabolite._model = self
            metabolite._reaction = {self.reactions.get_by_id(re.id) for re in metabolite.reactions}

        for gene in self.genes:
            gene._model = self
            gene._reaction = {self.reactions.get_by_id(re.id) for re in gene.reactions}

        for reaction in self.reactions:
            reaction._genes = {self.genes.get_by_id(gene.id) for gene in reaction.genes}
            reaction._metabolites = {
                self.metabolites.get_by_id(met.id): coef for met, coef in reaction.metabolites.items()
            }

        self._solver = solver_interface.Model()
        self._populate_solver(self.reactions)
        self._timestamp_last_optimization = None
        self.solution = LazySolution(self)

    def __copy__(self):
        return self.__deepcopy__()

    def __deepcopy__(self):
        return self.copy()

    def copy(self):
        """Needed for compatibility with cobrapy."""
        model_copy = super(SolverBasedModel, self).copy()
        try:
            model_copy._solver = deepcopy(self.solver)
        except Exception:  # pragma: no cover # Cplex has an issue with deep copies
            model_copy._solver = copy(self.solver)  # pragma: no cover
        return model_copy

    def _repr_html_(self):  # pragma: no cover
        template = """<table>
<tr>
<td>Name</td>
<td>%(name)s</td>
</tr>
<tr>
<td>Number of metabolites</td>
<td>%(num_metabolites)s</td>
</tr>
<tr>
<td>Number of reactions</td>
<td>%(num_reactions)s</td>
</tr>
<tr>
<td>Reactions</td>
<td><div style="width:100%%; max-height:300px; overflow:auto">%(reactions)s</div></td>
</tr>
</table>"""
        return template % {'name': self.id, 'num_metabolites': len(self.metabolites),
                           'num_reactions': len(self.reactions),
                           'reactions': '<br>'.join([r.build_reaction_string() for r in self.reactions])}

    @property
    def objective(self):
        """The model objective."""
        return self.solver.objective

    @objective.setter
    def objective(self, value):
        if isinstance(value, six.string_types):
            try:
                value = self.reactions.get_by_id(value)
            except KeyError:
                raise ValueError("No reaction with the id %s in the model" % value)
        if isinstance(value, Reaction):
            if value.model is not self:
                raise ValueError("%r does not belong to the model" % value)
            self.solver.objective = self.solver.interface.Objective(value.flux_expression, sloppy=True)
        elif isinstance(value, self.solver.interface.Objective):
            self.solver.objective = value
        # TODO: maybe the following should be allowed
        # elif isinstance(value, optlang.interface.Objective):
        # self.solver.objective = self.solver.interface.Objective.clone(value)
        elif isinstance(value, sympy.Basic):
            self.solver.objective = self.solver.interface.Objective(value, sloppy=False)
        else:
            raise TypeError('%r is not a valid objective for %r.' % (value, self.solver))

    @property
    def solver(self):
        """Attached solver instance.

        Very useful for accessing the optimization problem directly. Furthermore, can be used to define additional
        non-metabolic constraints.

        Examples
        --------
        >>> new_constraint_from_objective = model.solver.interface.Constraint(model.objective.expression, lb=0.99)
        >>> model.solver.add(new_constraint)

        """
        return self._solver

    @solver.setter
    def solver(self, value):
        not_valid_interface = ValueError(
            '%s is not a valid solver interface. Pick from %s, or specify an optlang interface (e.g. optlang.glpk_interface).' % (
                value, list(config.solvers.keys())))
        if isinstance(value, six.string_types):
            try:
                interface = config.solvers[value]
            except KeyError:
                raise not_valid_interface
        elif isinstance(value, types.ModuleType) and hasattr(value, 'Model'):
            interface = value
        else:
            raise not_valid_interface
        self._solver = interface.Model.clone(self._solver)

    @property
    def exchanges(self):
        """Exchange reactions in model.

        Reactions that either don't have products or substrates.
        """
        return [reaction for reaction in self.reactions if len(reaction.reactants) == 0 or len(reaction.products) == 0]

    def add_metabolites(self, metabolite_list):
        super(SolverBasedModel, self).add_metabolites(metabolite_list)
        for met in metabolite_list:
            if met.id not in self.solver.constraints:
                self.solver.add(self.solver.interface.Constraint(S.Zero, name=met.id, lb=0, ub=0))

    def _populate_solver(self, reaction_list):
        """Populate attached solver with constraints and variables that model the provided reactions."""
        constr_terms = dict()
        objective_terms = list()
        metabolites = {}
        for reaction in reaction_list:

            if reaction.reversibility:
                forward_variable = self.solver.interface.Variable(reaction._get_forward_id(), lb=0,
                                                                  ub=reaction._upper_bound)
                reverse_variable = self.solver.interface.Variable(reaction._get_reverse_id(), lb=0,
                                                                  ub=-1 * reaction._lower_bound)
            elif 0 == reaction.lower_bound and reaction.upper_bound == 0:
                forward_variable = self.solver.interface.Variable(reaction._get_forward_id(), lb=0, ub=0)
                reverse_variable = self.solver.interface.Variable(reaction._get_reverse_id(), lb=0, ub=0)
            elif reaction.lower_bound >= 0:
                forward_variable = self.solver.interface.Variable(reaction.id, lb=reaction._lower_bound,
                                                                  ub=reaction._upper_bound)
                reverse_variable = self.solver.interface.Variable(reaction._get_reverse_id(), lb=0, ub=0)
            elif reaction.upper_bound <= 0:
                forward_variable = self.solver.interface.Variable(reaction.id, lb=0, ub=0)
                reverse_variable = self.solver.interface.Variable(reaction._get_reverse_id(),
                                                                  lb=-1 * reaction._upper_bound,
                                                                  ub=-1 * reaction._lower_bound)
            self.solver.add(forward_variable)
            self.solver.add(reverse_variable)

            for metabolite, coeff in six.iteritems(reaction.metabolites):
                if metabolite.id in constr_terms:
                    constr_terms[metabolite.id].append(
                        sympy.Mul._from_args([sympy.RealNumber(coeff), forward_variable]))
                else:
                    constr_terms[metabolite.id] = [sympy.Mul._from_args([sympy.RealNumber(coeff), forward_variable])]
                    metabolites[metabolite.id] = metabolite
                constr_terms[metabolite.id].append(
                    sympy.Mul._from_args([sympy.RealNumber(-1 * coeff), reverse_variable]))

            if reaction._objective_coefficient != 0.:
                objective_terms.append(reaction._objective_coefficient * reaction.flux_expression)

        new_constraints = list()
        for met_id, terms in six.iteritems(constr_terms):
            expr = sympy.Add._from_args(terms)
            try:
                self.solver.constraints[met_id] += expr
            except KeyError:
                new_constraints.append(self.solver.interface.Constraint(expr, name=met_id, lb=0, ub=0, sloppy=True))
        self.solver.add(new_constraints, sloppy=True)

        objective_expression = sympy.Add(*objective_terms)
        if self.solver.objective is None:
            self.solver.objective = self.solver.interface.Objective(objective_expression, name='obj', direction='max')
        else:
            # TODO: remove this weird hack. Looks like some weird issue with lazy objective expressions in CPLEX and GLPK interface in optlang.
            self.solver.objective.variables
            self.solver.objective += objective_expression

    def add_reactions(self, reaction_list):
        cloned_reaction_list = list()
        for reaction in reaction_list:  # this is necessary for cobrapy compatibility
            if not isinstance(reaction, Reaction):
                cloned_reaction_list.append(Reaction.clone(reaction, model=self))
            else:
                cloned_reaction_list.append(reaction)

        # cobrapy will raise an exceptions if one of the reactions already exists in the model (before adding any reactions)
        super(SolverBasedModel, self).add_reactions(cloned_reaction_list)
        self._populate_solver(cloned_reaction_list)

    def remove_reactions(self, the_reactions, delete=True, remove_orphans=False):
        super(SolverBasedModel, self).remove_reactions(the_reactions, delete=delete, remove_orphans=remove_orphans)

    def add_demand(self, metabolite, prefix="DM_", time_machine=None):
        from warnings import warn
        warn('"add_demand" function is replaced with "add_exchange".', PendingDeprecationWarning)
        return self.add_exchange(metabolite, prefix=prefix, time_machine=time_machine)

    def add_exchange(self, metabolite, demand=True, prefix='DM_', bound=1000.0, time_machine=None):
        """Add an exchange reaction for a metabolite (demand=TRUE: metabolite --> Ø or demand=False: 0 --> metabolite )

        Parameters
        ----------
        metabolite : Metabolite
        demand : bool, optional
            True for sink type exchange, False for uptake type exchange
        prefix : str, optional
            A prefix that will be added to the metabolite ID to be used as the demand reaction's ID (defaults to 'DM_').
        bound : float, optional
            Upper bound for sink reaction / lower bound for uptake (multiplied by -1)
        time_machine : TimeMachine, optional
            A TimeMachine instance that enables undoing.

        Returns
        -------
        Reaction
            The created demand reaction.
        """
        id = str(prefix + metabolite.id)
        name = "Exchange %s" % metabolite.name if prefix != "DM_" else "Demand %s" % metabolite.name
        if id in self.reactions:
            raise ValueError("The metabolite already has a demand reaction.")

        reaction = Reaction()
        reaction.id = id
        reaction.name = name

        reaction.add_metabolites({metabolite: -1})
        if demand:
            reaction.upper_bound = bound
            reaction.lower_bound = 0
        else:
            reaction.upper_bound = 0
            reaction.lower_bound = -bound



        if time_machine is not None:
            time_machine(do=partial(self.add_reactions, [reaction]),
                         undo=partial(self.remove_reactions, [reaction], delete=False))
        else:
            self.add_reactions([reaction])
        return reaction

    def fix_objective_as_constraint(self, time_machine=None, fraction=1):
        """Fix current objective as an additional constraint (e.g., ..math`c^T v >= max c^T v`).

        Parameters
        ----------
        time_machine : TimeMachine, optional
            A TimeMachine instance can be provided, making it easy to undo this modification.

        Returns
        -------
        None
        """
        fix_objective_name = 'Fixed_objective_{}'.format(self.objective.name)
        if fix_objective_name in self.solver.constraints:
            self.solver.remove(fix_objective_name)
        objective_value = self.solve().objective_value * fraction
        constraint = self.solver.interface.Constraint(self.objective.expression,
                                                      name=fix_objective_name)
        if self.objective.direction == 'max':
            constraint.lb = objective_value
        else:
            constraint.ub = objective_value
        if time_machine is None:
            self.solver._add_constraint(constraint, sloppy=True)
        else:
            time_machine(do=partial(self.solver._add_constraint, constraint, sloppy=True),
                         undo=partial(self.solver.remove, constraint))

    def add_ratio_constraint(self, reaction1, reaction2, ratio, prefix='ratio_constraint_'):
        """Adds a ratio constraint (reaction1/reaction2 = ratio) to the model.

        Parameters
        ----------
        reaction1 : str or Reaction
            A reaction or a reaction ID.
        reaction2 : str or Reaction
            A reaction or a reaction ID.
        ratio : float
            The ratio in reaction1/reaction2 = ratio
        prefix : str
            The prefix that will be added to the constraint ID (defaults to 'ratio_constraint_').

        Returns
        -------
        optlang.Constraint
            The constraint name will be composed of `prefix`
            and the two reaction IDs (e.g. 'ratio_constraint_reaction1_reaction2').

        Examples
        --------
        model.add_ratio_constraint('r1', 'r2', 0.5)
        print model.solver.constraints['ratio_constraint_r1_r2']
        > ratio_constraint: ratio_constraint_r1_r2: 0 <= -0.5*r1 + 1.0*PGI <= 0
        """
        if isinstance(reaction1, bytes):
            reaction1 = self.reactions.get_by_id(reaction1)
        if isinstance(reaction2, bytes):
            reaction2 = self.reactions.get_by_id(reaction2)
        ratio_constraint = self.solver.interface.Constraint(
            reaction1.flux_expression - ratio * reaction2.flux_expression, lb=0, ub=0,
            name=prefix + reaction1.id + '_' + reaction2.id)
        self.solver._add_constraint(ratio_constraint, sloppy=True)
        return ratio_constraint

    def optimize(self, objective_sense=None, solution_type=Solution, **kwargs):
        """OptlangBasedModel implementation of optimize.

        Exists only for compatibility reasons. Uses model.solve() instead.
        """
        self._timestamp_last_optimization = time.time()
        if objective_sense is not None:
            original_direction = self.objective.direction
            self.objective.direction = {'minimize': 'min', 'maximize': 'max'}[objective_sense]
            self.solver.optimize()
            self.objective.direction = original_direction
        else:
            self.solver.optimize()
        solution = solution_type(self)
        self.solution = solution
        return solution

    def solve(self, solution_type=LazySolution, *args, **kwargs):
        """Optimize model.

        Parameters
        ----------
        solution_type : Solution or LazySolution, optional
            The type of solution that should be returned (defaults to LazySolution).

        Returns
        -------
        Solution or LazySolution
        """
        solution = self.optimize(solution_type=solution_type, *args, **kwargs)
        if solution.status is not 'optimal':
            logger.debug('No optimal solution found. Turning on presolve.')
            self.solver.configuration.presolve = True
            solution = self.optimize(solution_type=solution_type, *args, **kwargs)
            logger.debug('hey, %s, %s' % (solution.status, self.solver.status))
            self.solver.configuration.presolve = False
            if solution.status is not 'optimal':
                status = solution.status
                # GLPK 4.45 hack http://lists.gnu.org/archive/html/help-glpk/2013-09/msg00015.html
                if status == 'undefined' and self.solver.interface.__name__ == 'optlang.glpk_interface':  # and self.solver.interface.glp_version() == '4.45':
                    status = 'infeasible'
                raise exceptions._OPTLANG_TO_EXCEPTIONS_DICT.get(status, SolveError)(
                    'Solving model %s did not return an optimal solution. The returned solution status is "%s"' % (
                        self, status))
            return solution
        else:
            return solution

    def __dir__(self):
        # Hide 'optimize' from user.
        fields = sorted(dir(type(self)) + list(self.__dict__.keys()))
        fields.remove('optimize')
        return fields

    def essential_reactions(self, threshold=1e-6):
        """Return a list of essential reactions.

        Parameters
        ----------
        threshold : float (default 1e-6)
            Minimal objective flux to be considered viable.

        Returns
        -------
        list
            List of essential reactions
        """
        essential = []
        try:
            solution = self.solve()

            for reaction_id, flux in six.iteritems(solution.fluxes):
                if abs(flux) > 0:
                    reaction = self.reactions.get_by_id(reaction_id)
                    with TimeMachine() as tm:
                        reaction.knock_out(time_machine=tm)
                        try:
                            sol = self.solve()
                        except Infeasible:
                            essential.append(reaction)
                        else:
                            if sol.f < threshold:
                                essential.append(reaction)

        except SolveError as e:
            logger.error('Cannot determine essential reactions for un-optimal model.')
            raise e

        return essential

    def essential_genes(self, threshold=1e-6):
        """Return a list of essential genes.

        Parameters
        ----------
        threshold : float (default 1e-6)
            Minimal objective flux to be considered viable.

        Returns
        -------
        list
            List of essential genes
        """
        essential = []
        try:
            solution = self.solve()
            genes_to_check = set()
            for reaction_id, flux in six.iteritems(solution.fluxes):
                if abs(flux) > 0:
                    genes_to_check.update(self.reactions.get_by_id(reaction_id).genes)
            for gene in genes_to_check:
                with TimeMachine() as tm:
                    gene.knock_out(time_machine=tm)
                    try:
                        sol = self.solve()
                    except Infeasible:
                        essential.append(gene)
                    else:
                        if sol.f < threshold:
                            essential.append(gene)

        except SolveError as e:
            logger.error('Cannot determine essential genes for un-optimal model.')
            raise e

        return essential

    @property
    def medium(self):
        """Current medium."""
        reaction_ids = []
        reaction_names = []
        lower_bounds = []
        upper_bounds = []
        for ex in self.exchanges:
            metabolite = list(ex.metabolites.keys())[0]
            coeff = ex.metabolites[metabolite]
            if coeff * ex.lower_bound > 0:
                reaction_ids.append(ex.id)
                reaction_names.append(ex.name)
                lower_bounds.append(ex.lower_bound)
                upper_bounds.append(ex.upper_bound)

        return DataFrame({'reaction_id': reaction_ids,
                          'reaction_name': reaction_names,
                          'lower_bound': lower_bounds,
                          'upper_bound': upper_bounds},
                         index=None, columns=['reaction_id', 'reaction_name', 'lower_bound', 'upper_bound'])

    # TODO: describe the formats in doc
    def load_medium(self, medium, copy=False):
        """
        Loads a medium into the model. If copy is true it will return
        a copy of the model. Otherwise it applies the medium to itself.
        Supported formats
        TODO

        :param medium: can be a file, a pandas DataFrame or dictionary.
        :param copy: boolean, optional
        :return:
        """

        if copy:
            model = self.copy()
        else:
            model = self
        if isinstance(medium, dict):
            model._load_medium_from_dict(model, medium)
        elif isinstance(medium, pandas.DataFrame):
            model._load_medium_from_dataframe(model, medium)
        elif isinstance(medium, six.string_types):
            model._load_medium_from_file(model, medium)
        else:
            raise AssertionError("input type (%s) is not valid" % type(medium))

        return model

    def _ids_to_reactions(self, reactions):
        """Translate reaction IDs into reactions (skips reactions)."""
        clean_reactions = list()
        for reaction in reactions:
            if isinstance(reaction, six.string_types):
                clean_reactions.append(self.reactions.get_by_id(reaction))
            elif isinstance(reaction, Reaction):
                clean_reactions.append(reaction)
            else:
                raise Exception('%s is not a reaction or reaction ID.' % reaction)
        return clean_reactions

    def change_objective(self, value, time_machine=None):
        """
        Changes the objective of the model to the given value. Allows passing a time machine to
        revert the change later
        """
        if time_machine is None:
            self.objective = value
        else:
            time_machine(do=partial(setattr, self, "objective", value),
                         undo=partial(setattr, self, "objective", self.objective))

    def reaction_for(self, value, time_machine=None, add=True):
        """
        Converts an object into a reaction.

        If a Metabolite or a Metabolite id is given, it will return an exchange or demand reaction if it exists.
        If *add* is true, it adds a demand reaction if it does not exist.

        Parameters
        ----------
        value: str, Reaction or Metabolite
            An object that can be converted to a reaction
        time_machine: TimeMachine
            Can be used when *add* is True to revert the model
        add: bool
            Adds a demand reaction for a metabolite if a metabolite is found for *value*

        Returns
        -------
        Reaction

        Raises
        ------
        KeyError
            If *value* does not match any Reaction or Metabolite

        """

        if isinstance(value, Reaction):
            value = self.reactions.get_by_id(value.id)

        if isinstance(value, six.string_types):
            try:
                value = self.reactions.get_by_id(value)
            except KeyError:
                try:
                    value = self.metabolites.get_by_id(value)
                except KeyError:
                    raise KeyError("Invalid target %s." % value)

        if isinstance(value, cobra.core.Metabolite):
            try:
                value = self.reactions.get_by_id("EX_%s" % value.id)
            except KeyError:
                try:
                    value = self.reactions.get_by_id("DM_%s" % value.id)
                except KeyError as e:
                    if add:
                        value = self.add_demand(value, time_machine=time_machine)
                    else:
                        raise e

        if value is None:
            raise KeyError(None)

        return value

    @staticmethod
    def _load_medium_from_dict(model, medium):
        for rid, values in six.iteritems(medium):
            if model.reactions.has_id(rid):
                model.reactions.get_by_id(rid).lower_bound = values[0]
                model.reactions.get_by_id(rid).upper_bound = values[1]

    @staticmethod
    def _load_medium_from_file(model, file_path, delimiter="\t"):
        with open(file_path, "rb") as csv_file:
            reader = csv.reader(csv_file, delimiter=delimiter)
            for row in reader:
                if model.reactions.has_id(row[0]):
                    model.reactions.get_by_id(row[0]).lower_bound = float(row[1])
                    model.reactions.get_by_id(row[0]).upper_bound = float(row[2])

    @staticmethod
    def _load_medium_from_dataframe(model, medium):
        for i in six.moves.range(len(medium) - 1):
            rid = medium['reaction_id'][i]
            if model.reactions.has_id(rid):
                model.reactions.get_by_id(rid).lower_bound = medium['lower_bound'][i]
                model.reactions.get_by_id(rid).upper_bound = medium['upper_bound'][i]
