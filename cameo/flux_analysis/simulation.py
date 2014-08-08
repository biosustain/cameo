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
import copy

import time
from functools import partial
import sympy
from cameo.util import TimeMachine
from cameo.exceptions import SolveError
from sympy import Add
from sympy import Mul

add = Add._from_args
mul = Mul._from_args

def fba(model, objective=None, *args, **kwargs):
    """Perform flux balance analysis."""
    tm = TimeMachine()
    if objective is not None:
        tm(do=partial(setattr, model, 'objective', objective),
           undo=partial(setattr, model, 'objective', model.objective))
    try:
        solution = model.solve()
        tm.reset()
        return solution
    except SolveError as e:
        tm.reset()
        raise e


def pfba(model, objective=None, *args, **kwargs):
    tm = TimeMachine()
    tm(do=partial(setattr, model, 'reversible_encoding', 'split'),
       undo=partial(setattr, model, 'reversible_encoding', model.reversible_encoding))
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

        pfba_obj = model.solver.interface.Objective(add(
            [mul((sympy.singleton.S.One, variable)) for variable in model.solver.variables.values()]),
                                                    direction='min', sloppy=True)
        # tic = time.time()
        tm(do=partial(setattr, model, 'objective', pfba_obj),
           undo=partial(setattr, model, 'objective', model.objective))
        # print "obj: ", time.time() - tic
        try:
            solution = model.solve()
            tm.reset()
            return solution
        except SolveError as e:
            tm.reset()
            print "pfba could not determine an optimal solution for objective %s" % model.objective
            raise e
    except Exception as e:
        tm.reset()
        raise e


def moma(model, objective=None, *args, **kwargs):
    pass


def lmoma(model, wt_reference=None):
    tm = TimeMachine()
    original_objective = copy.copy(model.objective)
    try:
        obj_terms = list()
        constraints = list()
        variables = list()
        for rid, flux_value in wt_reference.iteritems():
            reaction = model.reactions.get_by_id(rid)
            pos_var = model.solver.interface.Variable("u_%s_pos" % rid, lb=0)
            neg_var = model.solver.interface.Variable("u_%s_neg" % rid, lb=0)

            model.solver._add_variable(pos_var)
            model.solver._add_variable(neg_var)
            variables.extend([pos_var, neg_var])

            obj_terms.append(pos_var)
            obj_terms.append(neg_var)

            # ui = vi - wt
            reaction_variable = reaction.variable
            constraint_a = model.solver.interface.Constraint(
                add([pos_var, mul([sympy.singleton.S.NegativeOne, reaction_variable])]), lb=-flux_value, sloppy=True)

            constraint_b = model.solver.interface.Constraint(add([neg_var, reaction.variable]), lb=flux_value,
                                                             sloppy=True)
            constraints.extend([constraint_a, constraint_b])

        for constraint in constraints:
            model.solver._add_constraint(constraint, sloppy=True)

        for term in obj_terms:
            model.solver._set_linear_objective_term(term, 1.)
        model.solver.objective.direction = 'min'

        try:
            solution = model.solve()
            return solution
        except SolveError as e:
            print "lmoma could not determine an optimal solution for objective %s" % model.objective
            raise e
    finally:
        tic = time.time()
        model.solver._remove_variables(variables)
        model.objective = original_objective
        model.solver._remove_constraints(constraints)
        tm.reset()
        print 'lmoma reset:', time.time() - tic


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
            obj_terms.append(mul([sympy.S.One, internal_reaction.variable]))
            tm(do=partial(setattr, internal_reaction, 'lower_bound', 0),
               undo=partial(setattr, internal_reaction, 'lower_bound', internal_reaction.lower_bound))
            tm(do=partial(setattr, internal_reaction, 'upper_bound', internal_flux),
               undo=partial(setattr, internal_reaction, 'upper_bound', internal_reaction.upper_bound))
        elif internal_flux < 0:
            obj_terms.append(mul([sympy.S.NegativeOne, internal_reaction.variable]))
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
                  model.solver.interface.Objective(add(obj_terms), name='Flux minimization',
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
    from cameo.solver_based_model import to_solver_based_model

    # sbml_path = '../../tests/data/EcoliCore.xml'
    sbml_path = '../../tests/data/iJO1366.xml'

    cb_model = read_sbml_model(sbml_path)
    model = load_model(sbml_path)

    # model.solver = 'glpk'

    print "cobra fba"
    tic = time.time()
    cb_model.optimize(solver='cglpk')
    print "flux sum:", sum([abs(val) for val in cb_model.solution.x_dict.values()])
    print "cobra fba runtime:", time.time() - tic

    print "cobra pfba"
    tic = time.time()
    optimize_minimal_flux(cb_model, solver='cglpk')
    print "flux sum:", sum([abs(val) for val in cb_model.solution.x_dict.values()])
    print "cobra pfba runtime:", time.time() - tic

    print "pfba"
    tic = time.time()
    solution = pfba(model)
    print "flux sum:",
    print sum([abs(val) for val in solution.x_dict.values()])
    print "cameo pfba runtime:", time.time() - tic

    print "lmoma"
    ref = solution.x_dict
    tic = time.time()
    solution = lmoma(model, wt_reference=ref)
    res = solution.x_dict
    print "flux distance:",
    print sum([abs(res[v] - ref[v]) for v in res.keys()])
    print "cameo lmoma runtime:", time.time() - tic

    print "lmoma w/ ko"
    tic = time.time()
    model.reactions.PGI.lower_bound = 0
    model.reactions.PGI.upper_bound = 0
    solution = lmoma(model, wt_reference=ref)
    res = solution.x_dict
    print "flux distance:",
    print sum([abs(res[v] - ref[v]) for v in res.keys()])
    print "cameo lmoma runtime:", time.time() - tic

    # print model.solver