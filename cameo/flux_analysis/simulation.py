# Copyright 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.
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
"""
Methods for simulating flux distributions.
Currently implements:
* fba - Flux Balance Analysis
* pfba - Parsimonious Flux Balance Analysis
* lmoma - (Linear) Minimization of Metabolic Adjustment
* room - Regulatory On/Off Minimization

"""




from __future__ import absolute_import, print_function
import cameo
from cameo.flux_analysis.distance import ManhattanDistance, RegulatoryOnOffDistance

__all__ = ['fba', 'pfba', 'moma', 'lmoma', 'room']

import six

from functools import partial

import sympy
from sympy import Add
from sympy import Mul

from cameo.util import TimeMachine
from cameo.exceptions import SolveError
from cameo.core.solution import Solution
from cameo.core.result import FluxDistributionResult

import logging
logger = logging.getLogger(__name__)

add = Add._from_args
mul = Mul._from_args
NegativeOne = sympy.singleton.S.NegativeOne
One = sympy.singleton.S.One
RealNumber = sympy.RealNumber


def fba(model, objective=None, *args, **kwargs):
    """Flux Balance Analysis.

    Parameters
    ----------
    model: SolverBasedModel
    objective: a valid objective - see SolverBaseModel.objective (optional)

    Returns
    -------
    FluxDistributionResult
        Contains the result of the linear solver.

    """
    with TimeMachine() as tm:
        if objective is not None:
            tm(do=partial(setattr, model, 'objective', objective),
               undo=partial(setattr, model, 'objective', model.objective))
        solution = model.solve()
        result = FluxDistributionResult(solution)
    return result

def pfba(model, objective=None, *args, **kwargs):
    """Parsimonious Flux Balance Analysis.

    Parameters
    ----------
    model: SolverBasedModel
    objective: str or reaction or optlang.Objective
        An objective to be minimized/maximized for

    Returns
    -------
    FluxDistributionResult
        Contains the result of the linear solver.

    """
    with TimeMachine() as tm:
        original_objective = model.objective.expression
        try:
            if objective is not None:
                tm(do=partial(setattr, model, 'objective', objective),
                   undo=partial(setattr, model, 'objective', original_objective))
            try:
                obj_val = fba(model)[model.objective]
            except SolveError as e:
                logger.debug("pfba could not determine maximum objective value for\n%s." % model.objective)
                raise e
            if model.objective.direction == 'max':
                fix_obj_constraint = model.solver.interface.Constraint(model.objective.expression, lb=obj_val)
            else:
                fix_obj_constraint = model.solver.interface.Constraint(model.objective.expression, ub=obj_val)
            tm(do=partial(model.solver._add_constraint, fix_obj_constraint),
               undo=partial(model.solver._remove_constraint, fix_obj_constraint))
            pfba_obj = model.solver.interface.Objective(
                add([mul([One, var]) for r in model.reactions for var in [r.forward_variable, r.reverse_variable]]),
                direction='min', sloppy=True)
            tm(do=partial(setattr, model, 'objective', pfba_obj),
               undo=partial(setattr, model, 'objective', original_objective))
            try:
                solution = model.solve()
                result = FluxDistributionResult(solution)
                tm.reset()
                return result
            except SolveError as e:
                logger.error("pfba could not determine an optimal solution for objective %s" % model.objective)
                raise e
        except Exception as e:
            raise e

def moma(model, reference=None, *args, **kwargs):
    raise NotImplementedError('Quadratic MOMA not yet implemented.')


def lmoma(model, reference=None, *args, **kwargs):
    """Linear Minimization Of Metabolic Adjustment.

    Parameters
    ----------
    model: SolverBasedModel
    reference: dict

    Returns
    -------
    FluxDistributionResult
        Contains the result of the linear solver.

    """

    if isinstance(model, cameo.Model):
        model = ManhattanDistance(model, reference)

    assert isinstance(model, ManhattanDistance)

    return model.minimize(*args, **kwargs)


def room(model, reference=None, delta=0.03, epsilon=0.001, *args, **kwargs):
    """Regulatory On/Off Minimization [1].

    Parameters
    ----------
    model: SolverBasedModel
    reference: dict or FluxDistributionResult
    objective: str or reaction or optlang.Objective
        An objective to be minimized/maximized for

    Returns
    -------
    FluxDistributionResult
        Contains the result of the linear solver.

    References
    ----------
    [1] Tomer Shlomi, Omer Berkman and Eytan Ruppin, "Regulatory on/off minimization of metabolic
    flux changes after genetic perturbations", PNAS 2005 102 (21) 7695-7700; doi:10.1073/pnas.0406346102
    """

    if isinstance(model, cameo.Model):
        model = RegulatoryOnOffDistance(model, reference, epsilon=epsilon, delta=delta)

    assert isinstance(model, RegulatoryOnOffDistance)

    return model.minimize(*args, **kwargs)

def _cycle_free_flux(model, fluxes, fix=[]):
    """Remove cycles from a flux-distribution (http://cran.r-project.org/web/packages/sybilcycleFreeFlux/index.html)."""
    # import time
    tm = TimeMachine()
    try:
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
        # tic = time.time()
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
        # print 'bounds', time.time() - tic
        for reaction_id in fix:
            reaction_to_fix = model.reactions.get_by_id(reaction_id)
            tm(do=partial(setattr, reaction_to_fix, 'lower_bound', fluxes[reaction_id]),
               undo=partial(setattr, reaction_to_fix, 'lower_bound', reaction_to_fix.lower_bound))
            tm(do=partial(setattr, reaction_to_fix, 'upper_bound', fluxes[reaction_id]),
               undo=partial(setattr, reaction_to_fix, 'upper_bound', reaction_to_fix.upper_bound))
        # tic = time.time()
        tm(do=partial(setattr, model, 'objective',
                      model.solver.interface.Objective(add(obj_terms), name='Flux minimization',
                                                       direction='min', sloppy=True)),
           undo=partial(setattr, model, 'objective', model.objective))
        # print 'blub', time.time() - tic
        try:
            # model.solver.configuration.verbosity = 3
            solution = model.solve()
            # model.solver.configuration.verbosity = 0
        except SolveError as e:
            logger.warning("Couldn't remove cycles from reference flux distribution.")
            raise e
        # print 'returning'
        return solution.x_dict
    finally:
        # tic = time.time()
        # print 'resetting'
        tm.reset()
        # print 'reset', time.time() - tic


def reset_model(model, cache):
    """
    When the simulation was not volatile, uses the cache to
    revert the model to it's original state.

    Parameters
    ----------
    model: SolverBasedModel
    cache: dict
        The cache must contain the added variables, constrains and
        the original objective

    """
    model.solver._remove_variables(list(cache['variables'].values()))
    model.objective = cache['original_objective']
    model.solver._remove_constraints(list(cache['constraints'].values()))


if __name__ == '__main__':
    import time
    from cobra.io import read_sbml_model
    from cobra.flux_analysis.parsimonious import optimize_minimal_flux
    from cameo import load_model

    # sbml_path = '../../tests/data/EcoliCore.xml'
    sbml_path = '../../tests/data/iJO1366.xml'

    cb_model = read_sbml_model(sbml_path)
    model = load_model(sbml_path)

    # model.solver = 'glpk'

    print("cobra fba")
    tic = time.time()
    cb_model.optimize(solver='cglpk')
    print("flux sum:", sum([abs(val) for val in list(cb_model.solution.x_dict.values())]))
    print("cobra fba runtime:", time.time() - tic)

    print("cobra pfba")
    tic = time.time()
    optimize_minimal_flux(cb_model, solver='cglpk')
    print("flux sum:", sum([abs(val) for val in list(cb_model.solution.x_dict.values())]))
    print("cobra pfba runtime:", time.time() - tic)

    print("pfba")
    tic = time.time()
    solution = pfba(model)
    print("flux sum:", end=' ')
    print(sum([abs(val) for val in list(solution.x_dict.values())]))
    print("cameo pfba runtime:", time.time() - tic)

    print("lmoma")
    ref = solution.x_dict
    tic = time.time()
    solution = lmoma(model, reference=ref)
    res = solution.x_dict
    print("flux distance:", end=' ')
    print(sum([abs(res[v] - ref[v]) for v in list(res.keys())]))
    print("cameo lmoma runtime:", time.time() - tic)

    print("room")
    tic = time.time()
    solution = room(model, reference=ref)
    res = solution.x_dict
    print(sum([abs(res[v] - ref[v]) for v in list(res.keys())]))
    print("cameo room runtime:", time.time() - tic)

    print("flux distance:", end=' ')
    print(sum([abs(res[v] - ref[v]) for v in list(res.keys())]))
    print("sum yi:", solution.f)
    print("cameo room runtime:", time.time() - tic)

    print("lmoma w/ ko")
    tic = time.time()
    model.reactions.PGI.lower_bound = 0
    model.reactions.PGI.upper_bound = 0
    solution = lmoma(model, reference=ref)
    print("flux distance:", end=' ')
    print(sum([abs(solution[v] - ref[v]) for v in list(solution.keys())]))
    print("cameo lmoma runtime:", time.time() - tic)

    # print model.solver