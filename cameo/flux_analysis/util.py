# -*- coding: utf-8 -*-
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

from functools import partial

from cameo.exceptions import SolveError
from cameo.util import TimeMachine

import sympy
from sympy import Add, Mul
from cameo.flux_analysis.simulation import add_pfba

import logging

__all__ = ['remove_infeasible_cycles', 'fix_pfba_as_constraint']

FloatOne = sympy.Float(1)
logger = logging.getLogger(__name__)


def remove_infeasible_cycles(model, fluxes, fix=()):
    """Remove thermodynamically infeasible cycles from a flux distribution.

    Arguments
    ---------
    model : SolverBasedModel
        The model that generated the flux distribution.
    fluxes : dict
        The flux distribution containing infeasible loops.

    Returns
    -------
    dict
        A cycle free flux distribution.

    References
    ----------
    .. [1]	A. A. Desouki, F. Jarre, G. Gelius-Dietrich, and M. J. Lercher, “CycleFreeFlux: efficient removal of
            thermodynamically infeasible loops from flux distributions.”
    """
    with TimeMachine() as tm:
        # make sure the original object is restored
        tm(do=int, undo=partial(setattr, model, 'objective', model.objective))
        exchange_reactions = model.exchanges
        exchange_ids = [exchange.id for exchange in exchange_reactions]
        internal_reactions = [reaction for reaction in model.reactions if reaction.id not in exchange_ids]
        for exchange in exchange_reactions:
            exchange_flux = fluxes[exchange.id]
            tm(do=partial(setattr, exchange, 'lower_bound', exchange_flux),
               undo=partial(setattr, exchange, 'lower_bound', exchange.lower_bound))
            tm(do=partial(setattr, exchange, 'upper_bound', exchange_flux),
               undo=partial(setattr, exchange, 'upper_bound', exchange.upper_bound))
        cycle_free_objective_list = []
        for internal_reaction in internal_reactions:
            internal_flux = fluxes[internal_reaction.id]
            if internal_flux >= 0:
                cycle_free_objective_list.append(Mul._from_args((FloatOne, internal_reaction.forward_variable)))
                tm(do=partial(setattr, internal_reaction, 'lower_bound', 0),
                   undo=partial(setattr, internal_reaction, 'lower_bound', internal_reaction.lower_bound))
                tm(do=partial(setattr, internal_reaction, 'upper_bound', internal_flux),
                   undo=partial(setattr, internal_reaction, 'upper_bound', internal_reaction.upper_bound))
            else:  # internal_flux < 0:
                cycle_free_objective_list.append(Mul._from_args((FloatOne, internal_reaction.reverse_variable)))
                tm(do=partial(setattr, internal_reaction, 'lower_bound', internal_flux),
                   undo=partial(setattr, internal_reaction, 'lower_bound', internal_reaction.lower_bound))
                tm(do=partial(setattr, internal_reaction, 'upper_bound', 0),
                   undo=partial(setattr, internal_reaction, 'upper_bound', internal_reaction.upper_bound))

        cycle_free_objective = model.solver.interface.Objective(
            Add._from_args(cycle_free_objective_list), direction="min", sloppy=True
        )
        model.objective = cycle_free_objective

        for reaction_id in fix:
            reaction_to_fix = model.reactions.get_by_id(reaction_id)
            tm(do=partial(setattr, reaction_to_fix, 'lower_bound', fluxes[reaction_id]),
               undo=partial(setattr, reaction_to_fix, 'lower_bound', reaction_to_fix.lower_bound))
            tm(do=partial(setattr, reaction_to_fix, 'upper_bound', fluxes[reaction_id]),
               undo=partial(setattr, reaction_to_fix, 'upper_bound', reaction_to_fix.upper_bound))

        try:
            solution = model.solve()
        except SolveError as e:
            logger.warning("Couldn't remove cycles from reference flux distribution.")
            raise e
        result = solution.x_dict
        return result


def fix_pfba_as_constraint(model, multiplier=1, fraction_of_optimum=1, time_machine=None):
    """Fix the pFBA optimum as a constraint

    Useful when setting other objectives, like the maximum flux through given reaction may be more realistic if all
    other fluxes are not allowed to reach their full upper bounds, but collectively constrained to max sum.

    Parameters
    ----------
    model : cameo.core.SolverBasedModel
        The model to add the pfba constraint to
    multiplier : float
        The multiplier of the minimal sum of all reaction fluxes to use as the constraint.
    fraction_of_optimum : float
        The fraction of the objective value's optimum to use as constraint when getting the pFBA objective's minimum
    time_machine : TimeMachine, optional
        A TimeMachine instance can be provided, making it easy to undo this modification.
    """

    fix_constraint_name = '_fixed_pfba_constraint'
    if fix_constraint_name in model.solver.constraints:
        model.solver.remove(fix_constraint_name)
    with TimeMachine() as tm:
        add_pfba(model, time_machine=tm, fraction_of_optimum=fraction_of_optimum)
        pfba_objective_value = model.optimize().objective_value * multiplier
        constraint = model.solver.interface.Constraint(model.objective.expression,
                                                       name=fix_constraint_name,
                                                       ub=pfba_objective_value)
    if time_machine is None:
        model.solver._add_constraint(constraint, sloppy=True)
    else:
        time_machine(do=partial(model.solver._add_constraint, constraint, sloppy=True),
                     undo=partial(model.solver.remove, constraint))
