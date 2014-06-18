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

import time
from copy import deepcopy, copy

import optlang

from cobra.core.Reaction import Reaction as OriginalReaction
from cobra.core.Model import Model
from cobra.core.DictList import DictList

import sympy
from sympy.core.add import _unevaluated_Add
from sympy.core.mul import _unevaluated_Mul
from sympy import Mul
from sympy.core.singleton import S

from pandas import Series, DataFrame


def to_solver_based_model(model, solver=None, deepcopy_model=True):
    """Convert a core model into a solver-based model."""

    if solver is None:
        # Generate the default optlang solver instance
        solver = optlang.Model()
    sbmodel = OptlangBasedModel(
        solver=solver, description=model, deepcopy_model=deepcopy_model)
    for y in ['reactions', 'genes', 'metabolites']:
        for x in sbmodel.__dict__[y]:
            setattr(x, '_model', sbmodel)
    new_reactions = DictList()
    for reaction in sbmodel.reactions:
        new_reaction = Reaction(name=reaction.name)
        new_reaction.__dict__.update(reaction.__dict__)
        new_reactions.append(new_reaction)
    sbmodel.reactions = new_reactions
    for metabolite in sbmodel.metabolites:
        metabolite._reaction = set([sbmodel.reactions.get_by_id(reaction.id)
                                    for reaction in metabolite._reaction])
    for gene in sbmodel.genes:
        gene._reaction = set([sbmodel.reactions.get_by_id(reaction.id)
                              for reaction in gene._reaction])
    return sbmodel


class LazySolution(object):
    """This class implements a lazily evaluating version of the original cobrapy Solution class."""

    def __init__(self, model):
        self.model = model
        self._time_stamp = self.model._timestamp_last_optimization
        self._f = None

    def __str__(self):
        return str(DataFrame({'primal': Series(self.x_dict), 'dual': Series(self.y_dict)}))

    def _check_freshness(self):
        if self._time_stamp != self.model._timestamp_last_optimization:
            print self._time_stamp
            print self.model._timestamp_last_optimization
            raise Exception(
                'The solution (capture around %s) has become invalid as the model has been re-optimized recently (%s).' % (
                    time.ctime(self._time_stamp), time.ctime(self.model._timestamp_last_optimization)))

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


class Reaction(OriginalReaction):
    """docstring for Reaction"""

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
            objective = model.objective
            # print model.solver.variables[self.id]
            # print value * model.solver.variables[self.id]
            # print model.objective.expression
            # model.objective.expression += value * model.solver.variables[self.id]
            expression = copy(model.objective.expression)
            model.objective = optlang.Objective(expression + value * model.solver.variables[self.id])

        self._objective_coefficient = value

        # def add_metabolites(self, *args, **kwargs):
        #     super(Reaction, self).add_metabolites(*args, **kwargs)
        #     return self


class OptlangBasedModel(Model):
    """Implements a model with an attached optlang solver instance.

    Every model manipulation is immediately reflected in the solver instance.
    """

    def __init__(self, solver=None, description=None, deepcopy_model=False, **kwargs):
        _timestamp_last_optimization = None
        if deepcopy_model and isinstance(description, Model):
            description = description.copy()
        super(OptlangBasedModel, self).__init__(description, **kwargs)
        if isinstance(solver, optlang.interface.Model):
            self._solver = solver
            self._populate_solver_from_scratch()
        else:
            self._solver = solver

    def __copy__(self):
        return self.__deepcopy__()

    def __deepcopy__(self):
        return self.copy()

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
            self.solver.objective = optlang.Objective(
                Mul._from_args([S.One, self.solver.variables[value]]), sloppy=True)
        elif isinstance(value, Reaction):
            self.solver.objective = optlang.Objective(
                Mul._from_args([S.One, self.solver.variables[value.id]]), sloppy=True)
        elif isinstance(value, optlang.Objective):
            self.solver.objective = value
        elif isinstance(value, sympy.Basic):
            self.solver.objective = optlang.Objective(value, sloppy=True)
        else:
            raise Exception('%s is not a valid objective.' % value)

    @property
    def solver(self):
        return self._solver

    @solver.setter
    def solver(self, value):
        if self._solver is None:
            self._solver = value
            self._populate_solver_from_scratch()
        else:
            previous_solver = self._solver
            self._solver = value
            self._populate_solver_from_other_solver(previous_solver)

            # TODO: something like this ...
            # self._populate_solver()
            # should follow ...
            # or if previous solver

    @property
    def exchanges(self):
        return [reaction for reaction in self.reactions if len(reaction.reactants) == 0 or len(reaction.products) == 0]

    def copy(self):
        """Needed for compatibility with cobrapy."""
        model_copy = super(OptlangBasedModel, self).copy()
        model_copy.solver = deepcopy(model_copy.solver)
        return model_copy

    def _populate_solver_from_scratch(self):
        """

        """
        objective = optlang.Objective(0)
        constr_terms = dict()
        for rxn in self.reactions:
            var = self.solver._add_variable(
                optlang.Variable(rxn.id, lb=rxn.lower_bound, ub=rxn.upper_bound))
            if rxn.objective_coefficient != 0.:
                objective.expression += rxn.objective_coefficient * var
            for met, coeff in rxn._metabolites.iteritems():
                if constr_terms.has_key(met.id):
                    constr_terms[met.id] += [(sympy.Real(coeff), var)]
                else:
                    constr_terms[met.id] = [(sympy.Real(coeff), var)]

        for met_id, terms in constr_terms.iteritems():
            expr = _unevaluated_Add(*[_unevaluated_Mul(coeff, var)
                                      for coeff, var in terms])
            constr = optlang.Constraint(expr, lb=0, ub=0, name=met_id)
            try:
                self.solver._add_constraint(constr, sloppy=True)
            except Exception, e:
                print e
                raise
        self.solver.objective = objective

    def _populate_solver_from_other_solver(self, other_solver):
        """

        """
        pass

    def add_metabolites(self, metabolite_list):
        super(OptlangBasedModel, self).add_metabolites(metabolite_list)
        for met in metabolite_list:
            self.solver.add(optlang.Constraint(name=met.name))

    def add_reactions(self, reaction_list):
        for reaction in reaction_list:
            try:
                reaction_variable = self.solver.variables[reaction.id]
            except KeyError:
                self.solver.add(optlang.Variable(reaction.id, lb=reaction.lower_bound, ub=reaction.upper_bound))
                reaction_variable = self.solver.variables[reaction.id]

            metabolite_coeff_dict = reaction.metabolites
            for metabolite, coeff in metabolite_coeff_dict.iteritems():
                try:
                    metabolite_constraint = self.solver.constraints[metabolite.id]
                except KeyError:
                    self.solver.add(optlang.Constraint(1, lb=0, ub=0, name=metabolite.id))
                    metabolite_constraint = self.solver.constraints[metabolite.id]
                metabolite_constraint += coeff * reaction_variable

        super(OptlangBasedModel, self).add_reactions(reaction_list)

    def remove_reactions(self, the_reactions):
        for reaction in the_reactions:
            self.solver.remove(reaction.id)
        super(OptlangBasedModel, self).remove_reactions(the_reactions)

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
                self.solver.objective = optlang.Objective(
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
            raise Exception(solution.status)
            return solution
        else:
            return solution

    def __dir__(self):
        # Hide 'optimize' from user.
        fields = sorted(dir(type(self)) + self.__dict__.keys())
        fields.remove('optimize')
        return fields


if __name__ == '__main__':
    import pickle
    import time
    from optlang.glpk_interface import Model as GlpkModel
    from cobra.manipulation.delete import delete_model_genes

    with open("../tests/data/iJO1366.pickle") as fhandle:
        model = pickle.load(fhandle)
        sbmodel = to_solver_based_model(model)

    choice = -444
    if choice == 1:
        solver = GlpkModel()
        sbmodel = to_solver_based_model(model, solver)
    elif choice == 2:
        from cobra.flux_analysis.variability import flux_variability_analysis

        from cobra.flux_analysis.variability import flux_variability_analysis

        with open("../tests/data/iJO1366.pickle") as fhandle:
            model = pickle.load(fhandle)
        model.reactions.get_by_id(
            'Ec_biomass_iJO1366_core_53p95M').lower_bound = 0.9823718127269797

        # Run standard cobrapy fva
        if False:
            t = time.time()
            fva_sol_model = flux_variability_analysis(
                model, the_reactions=model.reactions,
                fraction_of_optimum=1., solver="glpk")
            elapsed = time.time() - t
            print 'Cobrapy FVA time elapsed: ', elapsed
            # 473.228835106 last time I checked

        # Run solver-based-model fva
        solver = GlpkModel()
        sbmodel = to_solver_based_model(model, solver)
        sbmodel.solver.configuration.verbosity = 3
        # t = time.time()
        # fva_sol_sbmodel = flux_variability_analysis(
        #     sbmodel, the_reactions=sbmodel.reactions[0:100],
        #     fraction_of_optimum=1.
        # )
        # elapsed = time.time() - t
        # print 'optlang based FVA time elapsed: ', elapsed

        # Run fva manually
        fva_sol_sbmodel2 = dict()
        for rxn in sbmodel.reactions:
            fva_sol_sbmodel2[rxn.id] = dict()
        t = time.time()
        for rxn in sbmodel.reactions:
            # print "Minimize", rxn
            sbmodel.solver.objective = optlang.Objective(
                1. * sbmodel.solver.variables[rxn.id], direction='min')
            status = sbmodel.solver.optimize()

            if status is 'optimal':
                fva_sol_sbmodel2[rxn.id][
                    'minimum'] = sbmodel.solver.objective.value
            else:
                fva_sol_sbmodel2[rxn.id]['minimum'] = None
                # lb = sbmodel.solver.objective.value
                # if abs(lb) < 10**-6:
                #     lb = 0.
                # sbmodel.solver.variables[rxn.id].lb = lb
        for rxn in sbmodel.reactions:
            # print "Maximize", rxn
            sbmodel.solver.objective = optlang.Objective(
                1. * sbmodel.solver.variables[rxn.id], direction='max')
            status = sbmodel.solver.optimize()
            if status is 'optimal':
                fva_sol_sbmodel2[rxn.id][
                    'maximum'] = sbmodel.solver.objective.value
                # ub = sbmodel.solver.objective.value
                # if abs(ub) < 10**-6:
                #     ub = 0.
                # sbmodel.solver.variables[rxn.id].ub = ub
            else:
                fva_sol_sbmodel2[rxn.id][
                    'maximum'] = sbmodel.solver.objective.value
        elapsed = time.time() - t
        print 'optlang based FVA time elapsed: ', elapsed
        # 15.3940598965 last time checked
        # open('test.lp', 'w').write(str(sbmodel.solver))

        # stuff = sbmodel._populate_solver_from_scratch()
        # sbmodel = OptlangBasedModel(solver=solver)
        # sbmodel.add_metabolites(model.metabolites[0:3])
        # sbmodel.add_reactions(model.reactions)
        # print sbmodel.reactions.EX_12ppd__R_e
        # print sbmodel.solver
        # print sbmodel.solver
        # sbmodel.remove_reactions(['DM_AACALD'])
        # print sbmodel.solver
        # for key in fva_sol_model.keys():
        #     print key
        #     print fva_sol_model[key]
        #     print fva_sol_sbmodel[key]
        #     print fva_sol_sbmodel2[key]

    elif choice == 3:
        solver = GlpkModel()
        with open("../../test/data/salmonella.pickle") as fhandle:
            model = pickle.load(fhandle)

        cobra_model = to_solver_based_model(model, solver)
        print "move to test"
        # TODO: Add in tests for each function
        cumulative_deletions = False
        disable_orphans = False
        gene_list = ['STM1067', 'STM0227']
        # The following reactions are trimmed when STM1332 and STM1101 are
        # deleted
        dependent_reactions = set(['3HAD121',
                                   '3HAD160',
                                   '3HAD80',
                                   '3HAD140',
                                   '3HAD180',
                                   '3HAD100',
                                   '3HAD181',
                                   '3HAD120',
                                   '3HAD60',
                                   '3HAD141',
                                   '3HAD161',
                                   'T2DECAI',
                                   '3HAD40'])
        delete_model_genes(cobra_model, gene_list)
        symmetric_difference = dependent_reactions.symmetric_difference(
            [x.id for x in cobra_model._trimmed_reactions])
        if len(symmetric_difference) == 0:
            print 'Successful deletion of %s' % repr(gene_list)
        else:
            print 'Failed deletion of %s\n%s reactions did not match' % (repr(gene_list),
                                                                         repr(symmetric_difference))
    elif choice == 4:
        from functools import partial
        from time_machine import TimeMachine

        with open("../tests/data/salmonella.pickle") as fhandle:
            model = pickle.load(fhandle)
        model = to_solver_based_model(model)
        tm = TimeMachine()
        rxn = model.reactions.get_by_id('ENO')
        tm(do=partial(setattr, rxn, 'lower_bound', -666), undo=partial(setattr,
                                                                       rxn, 'lower_bound',
                                                                       model.reactions.get_by_id('ENO').lower_bound))
        print model.reactions.get_by_id('ENO').lower_bound
        tm.undo()
        print model.reactions.get_by_id('ENO').lower_bound
    elif choice == 5:
        with open("../../test/data/salmonella.pickle") as fhandle:
            model = to_solver_based_model(pickle.load(fhandle))
        print model.solver.constraints['atp_c']
        demand_rxn = Reaction(name='DM_atp')
        demand_rxn.add_metabolites({model.metabolites.get_by_id('atp_c'): -1})
        print demand_rxn.build_reaction_string()
        print "Is there already a model defined?", demand_rxn.get_model()
        model.add_reaction(demand_rxn)
        print model.reactions.get_by_id('DM_atp').build_reaction_string()
        print model.solver.constraints['atp_c']
