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

from cobra.exceptions import OptimizationError

import sympy
from sympy import Add, Mul
from cobra.flux_analysis.parsimonious import add_pfba

import logging

__all__ = ['remove_infeasible_cycles', 'fix_pfba_as_constraint']

FloatOne = sympy.Float(1)
logger = logging.getLogger(__name__)


def remove_infeasible_cycles(model, fluxes, fix=()):
    """Remove thermodynamically infeasible cycles from a flux distribution.

    Parameters
    ---------
    model : cobra.Model
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
    with model:
        # make sure the original object is restored
        exchange_reactions = model.boundary
        exchange_ids = [exchange.id for exchange in exchange_reactions]
        internal_reactions = [reaction for reaction in model.reactions if reaction.id not in exchange_ids]
        for exchange in exchange_reactions:
            exchange_flux = fluxes[exchange.id]
            exchange.bounds = (exchange_flux, exchange_flux)
        cycle_free_objective_list = []
        for internal_reaction in internal_reactions:
            internal_flux = fluxes[internal_reaction.id]
            if internal_flux >= 0:
                cycle_free_objective_list.append(Mul._from_args((FloatOne, internal_reaction.forward_variable)))
                internal_reaction.bounds = (0, internal_flux)
            else:  # internal_flux < 0:
                cycle_free_objective_list.append(Mul._from_args((FloatOne, internal_reaction.reverse_variable)))
                internal_reaction.bounds = (internal_flux, 0)
        cycle_free_objective = model.solver.interface.Objective(
            Add._from_args(cycle_free_objective_list), direction="min", sloppy=True
        )
        model.objective = cycle_free_objective

        for reaction_id in fix:
            reaction_to_fix = model.reactions.get_by_id(reaction_id)
            reaction_to_fix.bounds = (fluxes[reaction_id], fluxes[reaction_id])
        try:
            solution = model.optimize(raise_error=True)
        except OptimizationError as e:
            logger.warning("Couldn't remove cycles from reference flux distribution.")
            raise e
        result = solution.fluxes
        return result


def fix_pfba_as_constraint(model, multiplier=1, fraction_of_optimum=1):
    """Fix the pFBA optimum as a constraint

    Useful when setting other objectives, like the maximum flux through given reaction may be more realistic if all
    other fluxes are not allowed to reach their full upper bounds, but collectively constrained to max sum.

    Parameters
    ----------
    model : cobra.Model
        The model to add the pfba constraint to
    multiplier : float
        The multiplier of the minimal sum of all reaction fluxes to use as the constraint.
    fraction_of_optimum : float
        The fraction of the objective value's optimum to use as constraint when getting the pFBA objective's minimum
    """

    fix_constraint_name = '_fixed_pfba_constraint'
    if fix_constraint_name in model.solver.constraints:
        model.solver.remove(fix_constraint_name)
    with model:
        add_pfba(model, fraction_of_optimum=fraction_of_optimum)
        pfba_objective_value = model.slim_optimize(error_value=None) * multiplier
        constraint = model.solver.interface.Constraint(model.objective.expression,
                                                       name=fix_constraint_name,
                                                       ub=pfba_objective_value)
    model.add_cons_vars(constraint, sloppy=True)
