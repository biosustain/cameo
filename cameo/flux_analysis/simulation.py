# Copyright 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from functools import partial
import sympy
from cameo.util import TimeMachine
from cameo.exceptions import SolveError


def fba(model, objective=None):
    """Perform flux balance analysis."""
    tm = TimeMachine()
    if objective is not None:
        tm(do=partial(setattr, model, 'objective', objective),
           undo=partial(setattr, model, 'objective', model.objective))
    try:
        solution = model.solve()
    except SolveError:
        pass
    finally:
        tm.reset()
    return solution


def pfba(model, objective=None):
    tm = TimeMachine()
    try:
        if objective is not None:
            tm(do=partial(setattr, model, 'objective', objective),
               undo=partial(setattr, model, 'objective', model.objective))
        try:
            obj_val = model.solve().f
        except SolveError as e:
            print "pfba could not determine maximum objective value for\n%s." % model.objective
            raise e
        if model.objective.direction == 'max':
            fix_obj_constraint = model.solver.interface.Constraint(model.objective.expression, lb=obj_val)
        else:
            fix_obj_constraint = model.solver.interface.Constraint(model.objective.expression, ub=obj_val)
        tm(do=partial(model.solver._add_constraint, fix_obj_constraint),
           undo=partial(model.solver._remove_constraint, fix_obj_constraint))
        obj_terms = list()
        aux_constraint_terms = dict()
        for reaction in model.reactions:
            if reaction.reversibility:
                aux_variable = model.solver.interface.Variable(reaction.id + '_aux', lb=0, ub=-1 * reaction.lower_bound)
                tm(do=partial(setattr, reaction, 'lower_bound', 0),
                   undo=partial(setattr, reaction, 'lower_bound', reaction.lower_bound))
                obj_terms.append(aux_variable)
                tm(do=partial(model.solver._add_variable, aux_variable),
                   undo=partial(model.solver._remove_variable, aux_variable))
                for metabolite, coefficient in reaction.metabolites.iteritems():
                    try:
                        aux_constraint_terms[metabolite].append(-1 * coefficient * aux_variable)
                    except KeyError:
                        aux_constraint_terms[metabolite] = list()
                        aux_constraint_terms[metabolite].append(-1 * coefficient * aux_variable)
        pfba_obj = model.solver.interface.Objective(sympy.Add._from_args(obj_terms), direction='min')
        tm(do=partial(setattr, model, 'objective', pfba_obj),
           undo=partial(setattr, model, 'objective', model.objective))
        for metabolite, terms in aux_constraint_terms.iteritems():
            tm(do=partial(model.solver.constraints[metabolite.id].__iadd__, sympy.Add._from_args(terms)),
               undo=partial(model.solver.constraints[metabolite.id].__isub__, sympy.Add._from_args(terms)))
        try:
            solution = model.solve()
            result = dict()
            result['flux_sum'] = solution.f
            result['fluxes'] = solution.x_dict
        except SolveError as e:
            print "gimme could not determine an optimal solution for objective %s" % model.objective
            raise e
    finally:
        tm.reset()
    return result


def moma(model, objective=None):
    pass


def lmoma(model, objective=None, wt_reference=None):
    if wt_reference is None:
        wt_reference = pfba(model, objective=objective)

    tm = TimeMachine()
    try:
        obj_terms = list()
        aux_constraint_terms = dict()
        for reaction in model.reactions:
            obj_terms.append(reaction.variable)
            for metabolite, coefficient in reaction.metabolites.iteritems():
                try:
                    aux_constraint_terms[metabolite].append(reaction - reaction.x_dict[reaction.id])
                except KeyError:
                    aux_constraint_terms[metabolite] = list()
                    aux_constraint_terms[metabolite].append(reaction - reaction.x_dict[reaction.id])

        lmoma_obj = model.solver.interface.Objective(sympy.Add._from_args(obj_terms), direction='min')

        tm(do=partial(setattr, model, 'objective', lmoma_obj),
           undo=partial(setattr, model, 'objective', model.objective))
        for metabolite, terms in aux_constraint_terms.iteritems():
            tm(do=partial(model.solver.constraints[metabolite.id].__iadd__, sympy.Add._from_args(terms)),
               undo=partial(model.solver.constraints[metabolite.id].__isub__, sympy.Add._from_args(terms)))
        try:
            solution = model.solve()
        except SolveError as e:
            print "gimme could not determine an optimal solution for objective %s" % model.objective
            print model.solver
            raise e

    finally:
        # tic = time.time()
        tm.reset()
    # print time.time() - tic
    return solution

def _cycle_free_flux(model, fluxes, fix=[]):
    """Remove cycles from a flux-distribution (http://cran.r-project.org/web/packages/sybilcycleFreeFlux/index.html)."""
    tm = TimeMachine()
    exchange_reactions = model.exchanges
    exchange_ids = [exchange.id for exchange in exchange_reactions]
    internal_reactions = [reaction for reaction in model.reactions if reaction.id not in exchange_ids]
    for exchange in exchange_reactions:
        exchange_flux = fluxes[exchange.id]
        tm(do=partial(setattr, exchange, 'lower_bound', exchange_flux),
           undo=partial(setattr, exchange, 'lower_bound', exchange.lower_bound))
        tm(do=partial(setattr, exchange, 'upper_bound', exchange_flux),
           undo=partial(setattr, exchange, 'upper_bound', exchange.upper_bound))
    obj_terms = list()
    for internal_reaction in internal_reactions:
        internal_flux = fluxes[internal_reaction.id]
        if internal_flux >= 0:
            obj_terms.append(sympy.Mul._from_args([sympy.S.One, internal_reaction.variable]))
            tm(do=partial(setattr, internal_reaction, 'lower_bound', 0),
               undo=partial(setattr, internal_reaction, 'lower_bound', internal_reaction.lower_bound))
            tm(do=partial(setattr, internal_reaction, 'upper_bound', internal_flux),
               undo=partial(setattr, internal_reaction, 'upper_bound', internal_reaction.upper_bound))
        elif internal_flux < 0:
            obj_terms.append(sympy.Mul._from_args([sympy.S.NegativeOne, internal_reaction.variable]))
            tm(do=partial(setattr, internal_reaction, 'lower_bound', internal_flux),
               undo=partial(setattr, internal_reaction, 'lower_bound', internal_reaction.lower_bound))
            tm(do=partial(setattr, internal_reaction, 'upper_bound', 0),
               undo=partial(setattr, internal_reaction, 'upper_bound', internal_reaction.upper_bound))
        else:
            pass
    for reaction_id in fix:
        reaction_to_fix = model.reactions.get_by_id(reaction_id)
        tm(do=partial(setattr, reaction_to_fix, 'lower_bound', fluxes[reaction_id]),
           undo=partial(setattr, reaction_to_fix, 'lower_bound', reaction_to_fix.lower_bound))
        tm(do=partial(setattr, reaction_to_fix, 'upper_bound', fluxes[reaction_id]),
           undo=partial(setattr, reaction_to_fix, 'upper_bound', reaction_to_fix.upper_bound))
    tm(do=partial(setattr, model, 'objective',
                  model.solver.interface.Objective(sympy.Add._from_args(obj_terms), name='Flux minimization',
                                                   direction='min', sloppy=True)),
       undo=partial(setattr, model, 'objective', model.objective))
    solution = model.optimize()
    tm.reset()
    return solution.x_dict


if __name__ == '__main__':
    import time
    from cobra.io import read_sbml_model
    from cobra.flux_analysis.parsimonious import optimize_minimal_flux
    from cameo import load_model

    # sbml_path = '../../tests/data/EcoliCore.xml'
    sbml_path = '../../tests/data/iJO1366.xml'

    cb_model = read_sbml_model(sbml_path)
    cb_model.optimize()
    print sum([abs(val) for val in cb_model.solution.x_dict.values()])
    model = load_model(sbml_path)

    model.solver = 'glpk'

    print "cobra pfba"
    tic = time.time()
    optimize_minimal_flux(cb_model)
    print "flux sum:", sum([abs(val) for val in cb_model.solution.x_dict.values()])
    print "cobra pfba runtime:", time.time() - tic

    print "pfba"
    tic = time.time()
    solution = pfba(model)
    print "flux sum:",
    print sum([abs(val) for val in solution['fluxes'].values()])
    print "cameo pfba runtime:", time.time() - tic

    # print model.solver
