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

import time
import datetime
import csv
from copy import copy, deepcopy
from functools import partial
import types

import cobra as _cobrapy
import sympy
from sympy import Add
from sympy import Mul
from sympy.core.singleton import S
import optlang
from pandas import DataFrame, pandas

from cameo.util import TimeMachine
from cameo import exceptions
from cameo import config
from cameo.exceptions import SolveError, Infeasible, UndefinedSolution
from .reaction import Reaction
from .solution import LazySolution, Solution

import logging
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


class SolverBasedModel(_cobrapy.core.Model):
    """Implements a model with an attached optlang solver instance.

    Every model manipulation is immediately reflected in the solver instance.
    """

    def __init__(self, description=None, solver_interface=optlang, **kwargs):
        super(SolverBasedModel, self).__init__(description, **kwargs)
        self._reversible_encoding = 'split'
        cleaned_reactions = _cobrapy.core.DictList()
        for reaction in self.reactions:
            if isinstance(reaction, Reaction):
                cleaned_reactions.append(reaction)
            else:
                cleaned_reactions.append(Reaction.clone(reaction, model=self))
        self.reactions = cleaned_reactions
        for gene in self.genes:
            gene._model = self
            gene_reactions = list()
            for reaction in gene.reactions:
                model_reaction = self.reactions.get_by_id(reaction.id)
                gene_reactions.append(model_reaction)
            gene._reaction = set(gene_reactions)
        for metabolite in self.metabolites:
            metabolite._model = self
            metabolite_reactions = list()
            for reaction in metabolite.reactions:
                model_reaction = self.reactions.get_by_id(reaction.id)
                metabolite_reactions.append(model_reaction)
            metabolite._reaction = set(metabolite_reactions)
        self._solver = solver_interface.Model()
        self._populate_solver_from_scratch()
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
        except:  # pragma: no cover # Cplex has an issue with deep copies
            model_copy._solver = copy(self.solver) # pragma: no cover
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
        return self.solver.objective

    @objective.setter
    def objective(self, value):
        if isinstance(value, str):
            value = self.reactions.get_by_id(value)
        if isinstance(value, Reaction):
            if self.reversible_encoding == 'split' and value.reverse_variable is not None:
                obj_expression = Add._from_args((Mul._from_args((S.One, value.variable)), Mul._from_args((S.NegativeOne, value.reverse_variable))))
                self.solver.objective = self.solver.interface.Objective(obj_expression, sloppy=True)
            else:
                obj_expression = Mul._from_args((S.One, value.variable))
                self.solver.objective = self.solver.interface.Objective(obj_expression, sloppy=True)
        elif isinstance(value, self.solver.interface.Objective):
            self.solver.objective = value
        # TODO: maybe the following should be allowed
        # elif isinstance(value, optlang.interface.Objective):
        # self.solver.objective = self.solver.interface.Objective.clone(value)
        elif isinstance(value, sympy.Basic):
            self.solver.objective = self.solver.interface.Objective(value, sloppy=False)
        else:
            raise Exception('%s is not a valid objective.' % value)

    @property
    def solver(self):
        return self._solver

    @solver.setter
    def solver(self, value):
        not_valid_interface = ValueError('%s is not a valid solver interface. Pick from %s, or specify an optlang interface (e.g. optlang.glpk_interface).' % (value, config.solvers.keys()))
        if isinstance(value, types.StringType):
            try:
                interface = config.solvers[value]
            except KeyError:
                raise not_valid_interface
        elif isinstance(value, types.ModuleType) and hasattr(value, 'Model'):
            interface = value
        else:
            raise not_valid_interface
        self._solver = interface.Model()
        self._populate_solver_from_scratch()  #FIXME: This ignores non-reaction variables and constraints

    @property
    def exchanges(self):
        return [reaction for reaction in self.reactions if len(reaction.reactants) == 0 or len(reaction.products) == 0]

    def _populate_solver_from_scratch(self):
        objective_terms = list()
        constr_terms = dict()
        for rxn in self.reactions:
            lower_bound, upper_bound = rxn._lower_bound, rxn._upper_bound
            if lower_bound < 0 and upper_bound > 0:  # i.e. lower_bound < 0 and upper_bound > 0
                var = self.solver._add_variable(self.solver.interface.Variable(rxn.id, lb=0, ub=upper_bound))
                aux_var = self.solver._add_variable(
                    self.solver.interface.Variable(rxn._get_reverse_id(), lb=0, ub=-1 * lower_bound))
            else:
                var = self.solver._add_variable(
                    self.solver.interface.Variable(rxn.id, lb=lower_bound, ub=upper_bound))
            if rxn.objective_coefficient != 0.:
                objective_terms.append(sympy.Mul._from_args((sympy.RealNumber(rxn.objective_coefficient), var)))
            for met, coeff in rxn._metabolites.iteritems():
                if constr_terms.has_key(met.id):
                    constr_terms[met.id] += [(sympy.RealNumber(coeff), var)]
                else:
                    constr_terms[met.id] = [(sympy.RealNumber(coeff), var)]
                if lower_bound < 0 and upper_bound > 0:
                    constr_terms[met.id] += [(-1 * sympy.RealNumber(coeff), aux_var)]

        for met_id, terms in constr_terms.iteritems():
            expr = sympy.Add._from_args([sympy.Mul._from_args((coeff, var))
                                         for coeff, var in terms])
            constr = self.solver.interface.Constraint(expr, lb=0, ub=0, name=met_id)
            try:
                self.solver._add_constraint(constr, sloppy=False)  # TODO: should be True
            except Exception, e:
                print(e)
                raise
        objective_expression = sympy.Add._from_args(objective_terms)
        self.solver.objective = self.solver.interface.Objective(objective_expression, name='obj', direction='max')

    @property
    def reversible_encoding(self):
        return self._reversible_encoding

    @reversible_encoding.setter
    def reversible_encoding(self, value):
        if self._reversible_encoding == value:
            pass
        else:
            if value == 'unsplit':
                for reaction in self.reactions:
                    if reaction.reversibility:
                        reaction.variable.lb = -1 * reaction.reverse_variable.ub
                        reaction.reverse_variable.ub = 0
            elif value == 'split':
                for reaction in self.reactions:
                    if reaction.reversibility:
                        reaction.reverse_variable.ub = -1 * reaction.variable.lb
                        reaction.variable.lb = 0
            else:
                raise ValueError('%s is not a valid encoding. Try one of %s instead.' % (value, ('unsplit', 'split')))
            self._reversible_encoding = value

    def add_metabolites(self, metabolite_list):
        super(SolverBasedModel, self).add_metabolites(metabolite_list)
        for met in metabolite_list:
            if not self.solver.constraints.has_key(met.id):
                self.solver.add(self.solver.interface.Constraint(S.Zero, name=met.id, lb=0, ub=0))

    def add_reactions(self, reaction_list):
        cloned_reaction_list = list()
        for reaction in reaction_list:  # this is necessary for cobrapy compatibility
            if not isinstance(reaction, Reaction):
                cloned_reaction_list.append(Reaction.clone(reaction, model=self))
            else:
                cloned_reaction_list.append(reaction)

        # cobrapy will raise an exceptions if one of the reactions already exists in the model (before adding any reactions)
        super(SolverBasedModel, self).add_reactions(cloned_reaction_list)

        constr_terms = dict()
        metabolites = {}
        for reaction in cloned_reaction_list:
            if reaction.reversibility and self._reversible_encoding == "split":
                reaction_variable = self.solver.interface.Variable(reaction.id, lb=0, ub=reaction._upper_bound)
                aux_var = self.solver.interface.Variable(reaction._get_reverse_id(), lb=0, ub=-reaction._lower_bound)
                self.solver._add_variable(aux_var)
            else:
                reaction_variable = self.solver.interface.Variable(reaction.id, lb=reaction._lower_bound,
                                                                   ub=reaction._upper_bound)
            self.solver._add_variable(reaction_variable)

            for metabolite, coeff in reaction.metabolites.iteritems():
                if metabolite.id in constr_terms:
                    constr_terms[metabolite.id].append(
                        sympy.Mul._from_args([sympy.RealNumber(coeff), reaction_variable]))
                else:
                    constr_terms[metabolite.id] = [sympy.Mul._from_args([sympy.RealNumber(coeff), reaction_variable])]
                    metabolites[metabolite.id] = metabolite
                if reaction.reversibility and self._reversible_encoding == "split":
                    constr_terms[metabolite.id].append(sympy.Mul._from_args([sympy.RealNumber(-1*coeff), aux_var]))

        for met_id, terms in constr_terms.iteritems():
            expr = sympy.Add._from_args(terms)
            try:
                self.solver.constraints[met_id] += expr
            except KeyError:
                self.solver._add_constraint(self.solver.interface.Constraint(expr, name=met_id, lb=0, ub=0))

    def remove_reactions(self, the_reactions):
        for reaction in the_reactions:
            self.solver.remove(reaction.id)
        super(SolverBasedModel, self).remove_reactions(the_reactions)

    def add_demand(self, metabolite, prefix="DM_"):
        demand_reaction = Reaction(prefix + metabolite.id)
        demand_reaction.add_metabolites({metabolite: -1})
        demand_reaction.lower_bound = 0
        demand_reaction.upper_bound = 1000
        self.add_reactions([demand_reaction])
        return demand_reaction

    def add_ratio_constraint(self, reaction1, reaction2, ratio, prefix='ratio_constraint'):
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
            The prefix that will be added to the constraint ID.

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
        if isinstance(reaction1, types.StringType):
            reaction1 = self.reactions.get_by_id(reaction1)
        if isinstance(reaction2, types.StringType):
            reaction2 = self.reactions.get_by_id(reaction2)

        if reaction1.reverse_variable is not None:
            term1 = reaction1.variable - reaction1.reverse_variable
        else:
            term1 = reaction1.variable

        if reaction2.reverse_variable is not None:
            term2 = reaction2.variable - reaction2.reverse_variable
        else:
            term2 = reaction2.variable

        ratio_constraint = self.solver.interface.Constraint(term1 - ratio * term2, lb=0, ub=0, name='ratio_constraint_'+reaction1.id+'_'+reaction2.id)
        self.solver._add_constraint(ratio_constraint, sloppy=True)
        return ratio_constraint

    def optimize(self, new_objective=None, objective_sense=None, solution_type=Solution, **kwargs):
        """OptlangBasedModel implementation of optimize. Returns lazy solution object. Exists for compatibility reasons. Uses model.solve() instead."""
        if new_objective is None or new_objective == 0:
            pass
        else:
            # TODO: This i going to be deprecated soon anyway ...
            objective_formula = sympy.Add()
            [setattr(x, 'objective_coefficient', 0.) for x in self.reactions]
            if isinstance(new_objective, dict):
                for the_reaction, the_coefficient in new_objective.iteritems():
                    if isinstance(the_reaction, int):
                        the_reaction = self.reactions[the_reaction]
                    else:
                        if hasattr(the_reaction, 'id'):
                            the_reaction = the_reaction.id
                        the_reaction = self.reactions.get_by_id(the_reaction)
                    the_reaction.objective_coefficient = the_coefficient
                    objective_formula += the_coefficient * \
                                         self.solver.variables[the_reaction.id]
            else:
                # Allow for objectives to be constructed from multiple reactions
                if not isinstance(new_objective, list) and \
                        not isinstance(new_objective, tuple):
                    new_objective = [new_objective]
                for the_reaction in new_objective:
                    if isinstance(the_reaction, int):
                        the_reaction = self.reactions[the_reaction]
                    else:
                        if hasattr(the_reaction, 'id'):
                            the_reaction = the_reaction.id
                        the_reaction = self.reactions.get_by_id(the_reaction)
                    the_reaction.objective_coefficient = 1.
                    objective_formula += 1. * \
                                         self.solver.variables[the_reaction.id]

            if objective_formula != 0:
                self.solver.objective = self.solver.interface.Objective(
                    objective_formula, direction={'minimize': 'min', 'maximize': 'max'}[objective_sense])
        timestamp_formatter = lambda timestamp: datetime.datetime.fromtimestamp(timestamp).strftime(
            "%Y-%m-%d %H:%M:%S:%f")
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
        """Optimize model."""
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
                if status == 'undefined' and self.solver.interface.__name__ == 'optlang.glpk_interface' and self.solver.interface.glp_version() == '4.45':
                    status = 'infeasible'
                raise exceptions._OPTLANG_TO_EXCEPTIONS_DICT.get(status, SolveError)(
                    'Solving model %s did not return an optimal solution. The returned solution status is "%s"' % (
                        self, status))
            return solution
        else:
            return solution

    def __dir__(self):
        # Hide 'optimize' from user.
        fields = sorted(dir(type(self)) + self.__dict__.keys())
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
        except SolveError as e:
            print('Cannot determine essential reactions for un-optimal model.')
            raise e
        for reaction_id, flux in solution.fluxes.iteritems():
            if abs(flux) > 0:
                reaction = self.reactions.get_by_id(reaction_id)
                with TimeMachine() as tm:
                    reaction.knock_out(time_machine=tm)
                    try:
                        sol = self.solve()
                    except (Infeasible, UndefinedSolution):
                        essential.append(reaction)
                    else:
                        if sol.f < threshold:
                            essential.append(reaction)
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
        except SolveError as e:
            print('Cannot determine essential genes for un-optimal model.')
            raise e
        genes_to_check = set()
        for reaction_id, flux in solution.fluxes.iteritems():
            if abs(flux) > 0:
                genes_to_check.update(self.reactions.get_by_id(reaction_id).genes)
        for gene in genes_to_check:
            reactions = _cobrapy.manipulation.delete.find_gene_knockout_reactions(self, [gene])
            with TimeMachine() as tm:
                for reaction in reactions:
                    reaction.knock_out(time_machine=tm)
                try:
                    sol = self.solve()
                except (Infeasible, UndefinedSolution):
                    essential.append(gene)
                else:
                    if sol.f < threshold:
                        essential.append(gene)
        return essential

    def medium(self):
        reaction_ids = []
        reaction_names = []
        lower_bounds = []
        upper_bounds = []
        for ex in self.exchanges:
            metabolite = ex.metabolites.keys()[0]
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
        elif isinstance(medium, str):
            model._load_medium_from_file(model, medium)
        else:
            raise AssertionError("input type (%s) is not valid" % type(medium))

        return model

    def _ids_to_reactions(self, reactions):
        """Translate reaction IDs into reactions (skips reactions)."""
        clean_reactions = list()
        for reaction in reactions:
            if isinstance(reaction, str):
                clean_reactions.append(self.reactions.get_by_id(reaction))
            elif isinstance(reaction, Reaction):
                clean_reactions.append(reaction)
            else:
                raise Exception('%s is not a reaction or reaction ID.' % reaction)
        return clean_reactions

    @staticmethod
    def _load_medium_from_dict(model, medium):
        for rid, values in medium.iteritems():
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
        for i in xrange(len(medium) - 1):
            rid = medium['reaction_id'][i]
            if model.reactions.has_id(rid):
                model.reactions.get_by_id(rid).lower_bound = medium['lower_bound'][i]
                model.reactions.get_by_id(rid).upper_bound = medium['upper_bound'][i]