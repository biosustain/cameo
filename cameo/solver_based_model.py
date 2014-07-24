# Copyright 2013 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import csv

import time
from copy import deepcopy, copy
from cobra.core import Solution

import optlang

from cameo.util import TimeMachine
from cameo import exceptions
from cameo.exceptions import SolveError, Infeasible, Unbounded, FeasibleButNotOptimal, \
    UndefinedSolution
from cameo.flux_analysis.analysis import _flux_variability_analysis

from cobra.core.Reaction import Reaction as OriginalReaction
from cobra.core.Model import Model
from cobra.core.DictList import DictList
from cobra.manipulation.delete import find_gene_knockout_reactions

import sympy
from sympy.core.add import _unevaluated_Add
from sympy.core.mul import _unevaluated_Mul
from sympy import Mul
from sympy.core.singleton import S

from pandas import Series, DataFrame, pandas

from functools import partial

_SOLVER_INTERFACES = {}

try:
    from optlang import glpk_interface

    _SOLVER_INTERFACES['glpk'] = optlang.glpk_interface
except ImportError:
    pass
try:
    from optlang import cplex_interface

    _SOLVER_INTERFACES['cplex'] = optlang.cplex_interface
except ImportError:
    pass


def to_solver_based_model(cobrapy_model, solver_interface=optlang, deepcopy_model=True):
    """Convert a core model into a solver-based model."""

    solver_interface = _SOLVER_INTERFACES.get(solver_interface, solver_interface)
    solver_based_model = SolverBasedModel(
        solver_interface=solver_interface, description=cobrapy_model, deepcopy_model=deepcopy_model)
    for y in ['reactions', 'genes', 'metabolites']:
        for x in solver_based_model.__dict__[y]:
            setattr(x, '_model', solver_based_model)
    new_reactions = DictList()
    for reaction in solver_based_model.reactions:
        new_reaction = Reaction(name=reaction.name)
        new_reaction.__dict__.update(reaction.__dict__)
        new_reaction._objective_coefficient = reaction.objective_coefficient
        new_reactions.append(new_reaction)
    solver_based_model.reactions = new_reactions
    for metabolite in solver_based_model.metabolites:
        metabolite._reaction = set([solver_based_model.reactions.get_by_id(reaction.id)
                                    for reaction in metabolite._reaction])
    for gene in solver_based_model.genes:
        gene._reaction = set([solver_based_model.reactions.get_by_id(reaction.id)
                              for reaction in gene._reaction])
    return solver_based_model


class LazySolution(object):
    """This class implements a lazy evaluating version of the original cobrapy Solution class."""

    def __init__(self, model):
        self.model = model
        self._time_stamp = self.model._timestamp_last_optimization
        self._f = None

    def __str__(self):
        return str(DataFrame({'primal': Series(self.x_dict), 'dual': Series(self.y_dict)}))

    def _check_freshness(self):
        if self._time_stamp != self.model._timestamp_last_optimization:
            raise UndefinedSolution(
                'The solution (capture around %s) has become invalid as the model has been re-optimized recently (%s).' % (
                    time.ctime(self._time_stamp), time.ctime(self.model._timestamp_last_optimization)))

    def as_cobrapy_solution(self):
        return Solution(self.f, x=self.x,
                        x_dict=self.x_dict, y=self.y, y_dict=self.y_dict,
                        the_solver=None, the_time=0, status=self.status)

    @property
    def f(self):
        self._check_freshness()
        if self._f is None:
            return self.model.solver.objective.value
        else:
            return self._f

    @f.setter
    def f(self, value):
        self._f = value

    @property
    def x(self):
        self._check_freshness()
        return [variable.primal for variable in self.model.solver.variables.values()]

    @property
    def x_dict(self):
        self._check_freshness()
        return dict([(variable.name, variable.primal) for variable in self.model.solver.variables.values()])

    @property
    def y(self):
        self._check_freshness()
        return [variable.dual for variable in self.model.solver.variables.values()]

    @property
    def y_dict(self):
        self._check_freshness()
        return dict([(variable.name, variable.dual) for variable in self.model.solver.variables.values()])

    @property
    def status(self):
        self._check_freshness()
        return self.model.solver.status

    @property
    def primal_dict(self):
        return self.x_dict

    @property
    def dual_dict(self):
        return self.y_dict

    def get_primal_by_id(self, reaction_id):
        self._check_freshness()
        return self.model.reactions.get_by_id(reaction_id).variable.primal


class Reaction(OriginalReaction):
    """docstring for Reaction"""

    @classmethod
    def clone(cls, reaction, model=None):
        new_reaction = cls(name=reaction.name)
        for attribute, value in reaction.__dict__.iteritems():
            setattr(new_reaction, attribute, value)
        if not isinstance(reaction.get_model(), SolverBasedModel):
            new_reaction._model = None
        if model is not None:
            new_reaction._model = model
        return new_reaction

    def __init__(self, name=None):
        super(Reaction, self).__init__(name=name)
        self._lower_bound = 0
        self._upper_bound = 1000.
        self._objective_coefficient = 0.

    def __str__(self):
        return self.build_reaction_string()

    @property
    def variable(self):
        model = self.get_model()
        if model is not None:
            return model.solver.variables[self.id]
        else:
            return None

    @property
    def lower_bound(self):
        model = self.get_model()
        if model is not None:
            return model.solver.variables[self.id].lb
        else:
            return self._lower_bound

    @lower_bound.setter
    def lower_bound(self, value):
        if self.get_model() is not None:
            try:
                self.get_model().solver.variables[self.id].lb = value
            except ValueError:
                self.get_model().solver.variables[self.id].ub = value
                self.get_model().solver.variables[self.id].lb = value
            self._lower_bound = value
        else:
            self._lower_bound = value

    @property
    def upper_bound(self):
        if self.get_model() is not None:
            return self.get_model().solver.variables[self.id].ub
        else:
            return self._upper_bound

    @upper_bound.setter
    def upper_bound(self, value):
        if self.get_model() is not None:
            try:
                self.get_model().solver.variables[self.id].ub = value
            except ValueError:
                self.get_model().solver.variables[self.id].lb = value
                self.get_model().solver.variables[self.id].ub = value
            self._upper_bound = value
        else:
            self._upper_bound = value

    @property
    def objective_coefficient(self):
        return self._objective_coefficient

    @objective_coefficient.setter
    def objective_coefficient(self, value):
        model = self.get_model()
        if model is not None:
            model.solver._set_linear_objective_term(self.variable, value)

        self._objective_coefficient = value

    @property
    def effective_lower_bound(self):
        model = self.get_model()
        return _flux_variability_analysis(model, reactions=[self])[self.id]['minimum']

    @property
    def effective_upper_bound(self):
        model = self.get_model()
        return _flux_variability_analysis(model, reactions=[self])[self.id]['maximum']


class SolverBasedModel(Model):
    """Implements a model with an attached optlang solver instance.

    Every model manipulation is immediately reflected in the solver instance.
    """

    def __init__(self, solver_interface=optlang, description=None, deepcopy_model=False, **kwargs):
        if deepcopy_model and isinstance(description, Model):
            description = description.copy()
        super(SolverBasedModel, self).__init__(description, **kwargs)
        self._solver = solver_interface.Model()
        self._populate_solver_from_scratch()

    def __copy__(self):
        return self.__deepcopy__()

    def __deepcopy__(self):
        return self.copy()

    def copy(self):
        """Needed for compatibility with cobrapy."""
        model_copy = super(SolverBasedModel, self).copy()
        model_copy._solver = deepcopy(model_copy.solver)
        return model_copy

    def _repr_html_(self):
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
            self.solver.objective = self.solver.interface.Objective(
                Mul._from_args([S.One, self.solver.variables[value]]), sloppy=True)
        elif isinstance(value, Reaction):
            self.solver.objective = self.solver.interface.Objective(
                Mul._from_args([S.One, self.solver.variables[value.id]]), sloppy=True)
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
        interface = _SOLVER_INTERFACES.get(value, value)
        if self._solver is None:
            self._solver = interface.Model()
            self._populate_solver_from_scratch()
        else:
            self._solver = interface.Model.clone(self._solver)

    @property
    def exchanges(self):
        return [reaction for reaction in self.reactions if len(reaction.reactants) == 0 or len(reaction.products) == 0]

    def _populate_solver_from_scratch(self):
        objective = self.solver.interface.Objective(0)
        constr_terms = dict()
        for rxn in self.reactions:
            var = self.solver._add_variable(
                self.solver.interface.Variable(rxn.id, lb=rxn.lower_bound, ub=rxn.upper_bound))
            if rxn.objective_coefficient != 0.:
                objective.expression += rxn.objective_coefficient * var
            for met, coeff in rxn._metabolites.iteritems():
                if constr_terms.has_key(met.id):
                    constr_terms[met.id] += [(sympy.RealNumber(coeff), var)]
                else:
                    constr_terms[met.id] = [(sympy.RealNumber(coeff), var)]

        for met_id, terms in constr_terms.iteritems():
            expr = _unevaluated_Add(*[_unevaluated_Mul(coeff, var)
                                      for coeff, var in terms])
            constr = self.solver.interface.Constraint(expr, lb=0, ub=0, name=met_id)
            try:
                self.solver._add_constraint(constr, sloppy=True)
            except Exception, e:
                print e
                raise
        self.solver.objective = objective

    def add_metabolites(self, metabolite_list):
        super(SolverBasedModel, self).add_metabolites(metabolite_list)
        for met in metabolite_list:
            self.solver.add(self.solver.interface.Constraint(S.Zero, name=met.name, lb=0, ub=0))

    def add_reactions(self, reaction_list):
        cloned_reaction_list = list()
        for reaction in reaction_list:  # this is necessary for cobrapy compatibility
            if not isinstance(reaction, Reaction):
                cloned_reaction_list.append(Reaction.clone(reaction, model=self))
            else:
                cloned_reaction_list.append(reaction)
        for reaction in cloned_reaction_list:
            try:
                reaction_variable = self.solver.variables[reaction.id]
            except KeyError:
                self.solver.add(
                    self.solver.interface.Variable(reaction.id, lb=reaction._lower_bound, ub=reaction._upper_bound))
                reaction_variable = self.solver.variables[reaction.id]

            metabolite_coeff_dict = reaction.metabolites
            for metabolite, coeff in metabolite_coeff_dict.iteritems():
                try:
                    metabolite_constraint = self.solver.constraints[metabolite.id]
                except KeyError:  # cannot override add_metabolites here as it is not used by cobrapy in add_reactions
                    self.solver.add(self.solver.interface.Constraint(S.Zero, lb=0, ub=0,
                                                                     name=metabolite.id))  # TODO: 1 will not work ...
                    metabolite_constraint = self.solver.constraints[metabolite.id]
                metabolite_constraint += coeff * reaction_variable

        super(SolverBasedModel, self).add_reactions(cloned_reaction_list)

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

    def optimize(self, new_objective=None, objective_sense='maximize', solution_type=LazySolution, **kwargs):
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
        self._timestamp_last_optimization = time.time()
        self.solver.optimize()
        solution = solution_type(self)
        self.solution = solution
        return solution

    def solve(self, *args, **kwargs):
        """Optimize model."""
        solution = self.optimize(*args, **kwargs)
        if solution.status is not 'optimal':
            self.solver.configuration.presolve = True
            solution = self.optimize(*args, **kwargs)
            self.solver.configuration.presolve = False
            if solution.status is not 'optimal':
                raise exceptions._OPTLANG_TO_EXCEPTIONS_DICT.get(solution.status, SolveError)(
                    'Solving model %s did not return an optimal solution. The returned solution status is "%s"' % (
                        self, solution.status))
        else:
            return solution

    def __dir__(self):
        # Hide 'optimize' from user.
        fields = sorted(dir(type(self)) + self.__dict__.keys())
        fields.remove('optimize')
        return fields

    def essential_reactions(self, threshold=1e-6):
        essential = []
        time_machine = TimeMachine()
        try:
            solution = self.solve()
        except SolveError as e:
            print 'Cannot determine essential reactions for un-optimal model.'
            raise e
        for reaction_id, flux in solution.x_dict.iteritems():
            if abs(flux) > 0:
                reaction = self.reactions.get_by_id(reaction_id)
                time_machine(do=partial(setattr, reaction, 'lower_bound', 0),
                             undo=partial(setattr, reaction, 'lower_bound', reaction.lower_bound))
                time_machine(do=partial(setattr, reaction, 'upper_bound', 0),
                             undo=partial(setattr, reaction, 'upper_bound', reaction.upper_bound))
                try:
                    sol = self.solve()
                except (Infeasible, UndefinedSolution):
                    essential.append(reaction)
                else:
                    if sol.f < threshold:
                        essential.append(reaction)
                finally:
                    time_machine.reset()
        return essential

    def essential_genes(self, threshold=1e-6):
        essential = []
        time_machine = TimeMachine()
        try:
            solution = self.solve()
        except SolveError as e:
            print 'Cannot determine essential genes for unoptimal model.'
            raise e
        genes_to_check = set()
        for reaction_id, flux in solution.x_dict.iteritems():
            if abs(flux) > 0:
                genes_to_check.update(self.reactions.get_by_id(reaction_id).genes)
        for gene in genes_to_check:
            reactions = find_gene_knockout_reactions(self, [gene])
            for reaction in reactions:
                time_machine(do=partial(setattr, reaction, 'lower_bound', 0),
                             undo=partial(setattr, reaction, 'lower_bound', reaction.lower_bound))
                time_machine(do=partial(setattr, reaction, 'upper_bound', 0),
                             undo=partial(setattr, reaction, 'upper_bound', reaction.upper_bound))
            try:
                sol = self.solve()
            except (Infeasible, UndefinedSolution):
                essential.append(gene)
            else:
                if sol.f < threshold:
                    essential.append(gene)
                time_machine.reset()
            finally:
                time_machine.reset()
        return essential

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
            self._load_medium_from_dict(model, medium)
        elif isinstance(medium, pandas.DataFrame):
            self._load_medium_from_dataframe(model, medium)
        elif isinstance(medium, str):
            self._load_medium_from_file(model, medium)
        else:
            raise AssertionError("input type (%s) is not valid" % type(medium))

        return model


    @staticmethod
    def _load_medium_from_dict(model, medium):
        for rid, values in medium.iteritens():
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

if __name__ == '__main__':
    from cameo import load_model

    model = load_model('../../test/data/EcoliCore.xml')
