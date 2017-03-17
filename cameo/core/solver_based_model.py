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

import csv
import logging
import time
import types
from copy import copy, deepcopy
from functools import partial

import cobra
import numpy as np
import optlang
import six
import sympy
from pandas import DataFrame, pandas
from sympy import Add
from sympy import Mul
from sympy.core.singleton import S

from cameo import config
from cameo import exceptions
from cameo.core.gene import Gene
from cameo.core.metabolite import Metabolite
from cameo.exceptions import SolveError, Infeasible
from cameo.util import TimeMachine, inheritdocstring, AutoVivification
from .reaction import Reaction
from .solution import LazySolution, Solution

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
                reaction.model = self
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
                metabolite._model = self
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
                self.metabolites.get_by_id(met.id): coef for met, coef in six.iteritems(reaction.metabolites)
            }

        self._solver = solver_interface.Model()
        self._solver.objective = solver_interface.Objective(S.Zero)
        self._populate_solver(self.reactions, self.metabolites)
        self._timestamp_last_optimization = None
        self.solution = LazySolution(self)

    @property
    def non_functional_genes(self):
        """All non-functional genes in this model
        Returns
        -------
        frozenset
            set with the genes that are marked as non-functional
        """
        return frozenset(gene for gene in self.genes if not gene.functional)

    def __copy__(self):
        return self.__deepcopy__()

    def __deepcopy__(self):
        return self.copy()

    def copy(self):
        """Needed for compatibility with cobrapy."""
        model_copy = super(SolverBasedModel, self).copy()
        for reac in model_copy.reactions:
            reac._reset_var_cache()
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
            '%s is not a valid solver interface. '
            'Pick from %s, or specify an optlang interface (e.g. optlang.glpk_interface).' % (
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
        for reaction in self.reactions:
            reaction._reset_var_cache()
        self._solver = interface.Model.clone(self._solver)

    @property
    def exchanges(self):
        """Exchange reactions in model.

        Reactions that either don't have products or substrates.
        """
        return [reaction for reaction in self.reactions if reaction.is_exchange]

    def add_metabolites(self, metabolite_list):
        super(SolverBasedModel, self).add_metabolites(metabolite_list)
        for met in metabolite_list:
            if met.id not in self.solver.constraints:
                constraint = self.solver.interface.Constraint(S.Zero, name=met.id, lb=0, ub=0)
                self.solver.add(constraint)

    def add_metabolite(self, metabolite):
        self.add_metabolites([metabolite])

    def _populate_solver(self, reaction_list, metabolite_list=None):
        """Populate attached solver with constraints and variables that model the provided reactions."""
        constraint_terms = AutoVivification()
        if metabolite_list is not None:
            for met in metabolite_list:
                constraint = self.solver.interface.Constraint(S.Zero, name=met.id, lb=0, ub=0)
                self.solver.add(constraint)

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
            self.solver.update()

            for metabolite, coeff in six.iteritems(reaction.metabolites):
                if metabolite.id in self.solver.constraints:
                    constraint = self.solver.constraints[metabolite.id]
                else:
                    constraint = self.solver.interface.Constraint(S.Zero, name=metabolite.id, lb=0, ub=0)
                    self.solver.add(constraint, sloppy=True)

                constraint_terms[constraint][forward_variable] = coeff
                constraint_terms[constraint][reverse_variable] = -coeff

            objective_coeff = reaction._objective_coefficient
            if objective_coeff != 0.:
                if self.solver.objective is None:
                    self.solver.objective = self.solver.interface.Objective(0, direction='max')
                if self.solver.objective.direction == 'min':
                    self.solver.objective.direction = 'max'
                self.solver.objective.set_linear_coefficients(
                    {forward_variable: objective_coeff, reverse_variable: -objective_coeff})

        self.solver.update()
        for constraint, terms in six.iteritems(constraint_terms):
            constraint.set_linear_coefficients(terms)

    def add_reactions(self, reaction_list):
        cloned_reaction_list = list()
        for reaction in reaction_list:  # this is necessary for cobrapy compatibility
            if not isinstance(reaction, Reaction):
                cloned_reaction_list.append(Reaction.clone(reaction))
            else:
                cloned_reaction_list.append(reaction)

        # cobrapy will raise an exceptions if one of the reactions already exists in the model (before adding any
        # reactions)
        super(SolverBasedModel, self).add_reactions(cloned_reaction_list)
        for reac in cloned_reaction_list:
            reac.model = self
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

    def add_ratio_constraint(self, expr1, expr2, ratio, prefix='ratio_constraint_'):
        """Adds a ratio constraint (expr1/expr2 = ratio) to the model.

        Parameters
        ----------
        expr1 : str, Reaction, list or sympy.Expression
            A reaction, a reaction ID or a linear expression.
        expr2 : str, Reaction, list or sympy.Expression
            A reaction, a reaction ID or a linear expression.
        ratio : float
            The ratio in expr1/expr2 = ratio
        prefix : str
            The prefix that will be added to the constraint ID (defaults to 'ratio_constraint_').

        Returns
        -------
        optlang.Constraint
            The constraint name will be composed of `prefix`
            and the two reaction IDs (e.g. 'ratio_constraint_reaction1_reaction2').

        Examples
        --------
        >>> model.add_ratio_constraint('r1', 'r2', 0.5)
        >>> print(model.solver.constraints['ratio_constraint_r1_r2'])
        ratio_constraint: ratio_constraint_r1_r2: 0 <= -0.5 * r1 + 1.0 * PGI <= 0
        """
        if isinstance(expr1, six.string_types):
            prefix_1 = expr1
            expr1 = self.reactions.get_by_id(expr1).flux_expression
        elif isinstance(expr1, Reaction):
            prefix_1 = expr1.id
            expr1 = expr1.flux_expression
        elif isinstance(expr1, list):
            prefix_1 = "+".join(r.id for r in expr1)
            expr1 = sum([r.flux_expression for r in expr1], S.Zero)
        elif not isinstance(expr1, sympy.Expr):
            raise ValueError("'expr1' is not a valid expression")
        else:
            prefix_1 = str(expr1)

        if isinstance(expr2, six.string_types):
            prefix_2 = expr2
            expr2 = self.reactions.get_by_id(expr2).flux_expression
        elif isinstance(expr2, Reaction):
            prefix_2 = expr2.id
            expr2 = expr2.flux_expression
        elif isinstance(expr2, list):
            prefix_2 = "+".join(r.id for r in expr2)
            expr2 = sum([r.flux_expression for r in expr2], S.Zero)
        elif not isinstance(expr2, sympy.Expr):
            raise ValueError("'expr2' is not a valid expression")
        else:
            prefix_2 = str(expr2)

        ratio_constraint = self.solver.interface.Constraint(expr1 - ratio * expr2,
                                                            lb=0,
                                                            ub=0,
                                                            name=prefix + prefix_1 + '_' + prefix_2)

        self.solver.add(ratio_constraint, sloppy=True)
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
            raise exceptions._OPTLANG_TO_EXCEPTIONS_DICT.get(solution.status, SolveError)(
                'Solving model %s did not return an optimal solution. The returned solution status is "%s"' % (
                    self, solution.status))
        else:
            return solution

    def __dir__(self):
        # Hide 'optimize' from user.
        fields = sorted(dir(type(self)) + list(self.__dict__.keys()))
        fields.remove('optimize')
        return fields

    def essential_metabolites(self, threshold=1e-6, force_steady_state=False):
        """Return a list of essential metabolites.

        This can be done in 2 ways:

        1. Implementation follows the description in [1]:
            "All fluxes around the metabolite M should be restricted to only produce the metabolite,
             for which balancing constraint of mass conservation is relaxed to allow nonzero values
             of the incoming fluxes whereas all outgoing fluxes are limited to zero."

        2. Force Steady State approach:
            All reactions consuming the metabolite are restricted to only produce the metabolite. A demand
            reaction is added to sink the metabolite produced to keep the problem feasible under
            the S.v = 0 constraint.

        Briefly, for each metabolite, all reactions that consume that metabolite are blocked and if that makes the
        model either infeasible or results in near-zero flux in the model objective, then the metabolite is
        considered essential.

        Parameters
        ----------
        threshold : float (default 1e-6)
            Minimal objective flux to be considered viable.
        force_steady_state: bool
            If True, uses approach 2.

        References
        ----------
        .. [1] Kim, P.-J., Lee, D.-Y., Kim, T. Y., Lee, K. H., Jeong, H., Lee, S. Y., & Park, S. (2007).
         Metabolite essentiality elucidates robustness of Escherichia coli metabolism. PNAS, 104(34), 13638–13642
        """

        essential_metabolites = []

        # Essential metabolites are only in reactions that carry flux.
        metabolites = set()
        solution = self.solve()

        for reaction_id, flux in six.iteritems(solution.fluxes):
            if abs(flux) > 0:
                reaction = self.reactions.get_by_id(reaction_id)
                metabolites.update(reaction.metabolites.keys())

        for metabolite in metabolites:
            with TimeMachine() as tm:
                metabolite.knock_out(time_machine=tm, force_steady_state=force_steady_state)
                try:
                    solution = self.solve()
                    if solution.f < threshold:
                        essential_metabolites.append(metabolite)
                except Infeasible:
                    essential_metabolites.append(metabolite)

        return essential_metabolites

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
    def S(self):
        metabolite_index = {metabolite.id: index for index, metabolite in enumerate(self.metabolites)}
        stoichiometric_matrix = np.zeros((len(self.metabolites), len(self.reactions)))

        for i, reaction in enumerate(self.reactions):
            for metabolite, coefficient in six.iteritems(reaction.metabolites):
                j = metabolite_index[metabolite.id]
                stoichiometric_matrix[j, i] = coefficient

        return stoichiometric_matrix


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
    def load_medium(self, medium, copy=False, delimiter="\t"):
        """
        Loads a medium into the model. If copy is true it will return
        a copy of the model. Otherwise it applies the medium to itself.
        Supported formats
        TODO

        Parameters
        ----------
        medium: str, pandas.DataFrame, dict.

        copy: boolean, optional
            If True copies the model, otherwise the changes will happen inplace.
        delimiter: str
            Only if loading the medium from a file.

        Returns
        -------
        SolverBasedModel
            If copy=True, returns a copy of the model.

        """

        if copy:
            model = self.copy()
        else:
            model = self
        if isinstance(medium, dict):
            model._load_medium_from_dict(medium)
        elif isinstance(medium, pandas.DataFrame):
            model._load_medium_from_dataframe(medium)
        elif isinstance(medium, six.string_types):
            model._load_medium_from_file(medium, delimiter=delimiter)
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

    def _reaction_for(self, value, time_machine=None, add=True):
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

    def _load_medium_from_dict(self, medium):
        assert isinstance(medium, dict)
        for ex_reaction in self.exchanges:
            ex_reaction.lower_bound = medium.get(ex_reaction.id, 0)

    def _load_medium_from_file(self, file_path, delimiter="\t"):
        medium = {}

        with open(file_path, "rb") as csv_file:
            reader = csv.reader(csv_file, delimiter=delimiter)

            for row in reader:
                self.reactions.get_by_id(row[0])
                medium[row[0]] = row[1]

        self._load_medium_from_dict(medium)

    def _load_medium_from_dataframe(self, medium):
        assert isinstance(medium, DataFrame)
        for ex_reaction in self.exchanges:
            if ex_reaction.id in medium.reaction_id.values:
                medium_row = medium[medium.reaction_id == ex_reaction.id]
                ex_reaction.lower_bound = medium_row.lower_bound.values[0]
            else:
                ex_reaction.lower_bound = 0
