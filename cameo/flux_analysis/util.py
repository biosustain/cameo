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

__all__ = ['remove_infeasible_cycles']

import copy
from functools import partial

from cameo.exceptions import SolveError
from cameo.util import TimeMachine

import logging

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
        # make sure the orignal object is restored
        tm(do=int, undo=partial(setattr, model, 'objective', copy.copy(model.objective)))
        exchange_reactions = model.exchanges
        exchange_ids = [exchange.id for exchange in exchange_reactions]
        internal_reactions = [reaction for reaction in model.reactions if reaction.id not in exchange_ids]
        for exchange in exchange_reactions:
            exchange_flux = fluxes[exchange.id]
            tm(do=partial(setattr, exchange, 'lower_bound', exchange_flux),
               undo=partial(setattr, exchange, 'lower_bound', exchange.lower_bound))
            tm(do=partial(setattr, exchange, 'upper_bound', exchange_flux),
               undo=partial(setattr, exchange, 'upper_bound', exchange.upper_bound))
        for internal_reaction in internal_reactions:
            internal_flux = fluxes[internal_reaction.id]
            if internal_flux >= 0:
                model.solver._set_linear_objective_term(internal_reaction.forward_variable, 1.)
                tm(do=partial(setattr, internal_reaction, 'lower_bound', 0),
                   undo=partial(setattr, internal_reaction, 'lower_bound', internal_reaction.lower_bound))
                tm(do=partial(setattr, internal_reaction, 'upper_bound', internal_flux),
                   undo=partial(setattr, internal_reaction, 'upper_bound', internal_reaction.upper_bound))
            elif internal_flux < 0:
                model.solver._set_linear_objective_term(internal_reaction.reverse_variable, -1.)
                tm(do=partial(setattr, internal_reaction, 'lower_bound', internal_flux),
                   undo=partial(setattr, internal_reaction, 'lower_bound', internal_reaction.lower_bound))
                tm(do=partial(setattr, internal_reaction, 'upper_bound', 0),
                   undo=partial(setattr, internal_reaction, 'upper_bound', internal_reaction.upper_bound))
            else:
                pass
        model.objective.direction = 'min'

        for reaction_id in fix:
            reaction_to_fix = model.reactions.get_by_id(reaction_id)
            tm(do=partial(setattr, reaction_to_fix, 'lower_bound', fluxes[reaction_id]),
               undo=partial(setattr, reaction_to_fix, 'lower_bound', reaction_to_fix.lower_bound))
            tm(do=partial(setattr, reaction_to_fix, 'upper_bound', fluxes[reaction_id]),
               undo=partial(setattr, reaction_to_fix, 'upper_bound', reaction_to_fix.upper_bound))

        try:
            solution = model.solve()
            result = solution.x_dict
        except SolveError as e:
            logger.warning("Couldn't remove cycles from reference flux distribution.")
            raise e

    return result
