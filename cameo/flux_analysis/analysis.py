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


from copy import copy
import itertools
from collections import OrderedDict

import numpy
from cobra.core import Reaction

from cameo import config
from cameo.util import TimeMachine, partition
from functools import partial

from sympy import Add, Mul
from sympy.core.singleton import S

from optlang import Objective

from pandas import DataFrame


def _ids_to_reactions(model, reactions):
    """translates reaction IDs into reactions (skips reactions)."""
    clean_reactions = list()
    for reaction in reactions:
        if isinstance(reaction, str):
            clean_reactions.append(model.reactions.get_by_id(reaction))
        elif isinstance(reaction, Reaction):
            clean_reactions.append(reaction)
        else:
            raise Exception('%s is not a reaction or reaction ID.' % reaction)
    return clean_reactions


def _flux_variability_analysis(model, reactions=None):
    """Flux-variability analysis."""
    original_objective = copy(model.objective)
    if reactions is None:
        reactions = model.reactions
    else:
        reactions = _ids_to_reactions(model, reactions)
    fva_sol = OrderedDict()
    for reaction in reactions:
        fva_sol[reaction.id] = dict()
        model.objective = reaction
        model.objective.direction = 'min'
        solution = model.optimize()
        if solution.status == 'optimal':
            fva_sol[reaction.id]['minimum'] = model.solution.f
        else:
            fva_sol[reaction.id]['minimum'] = model.solution.status
    for reaction in reactions:
        model.objective = reaction
        model.objective.direction = 'max'
        solution = model.optimize()
        if solution.status == 'optimal':
            fva_sol[reaction.id]['maximum'] = model.solution.f
        else:
            fva_sol[reaction.id]['maximum'] = model.solution.status
    model.objective = original_objective
    return fva_sol


class FvaFunctionObject(object):
    """"""
    def __init__(self, model, fva):
        self.model = model
        self.fva = fva

    def __call__(self, reactions):
        return self.fva(self.model, reactions)


def flux_variability_analysis(model, reactions=None, view=config.default_view):
    """Flux-variability analysis."""
    if reactions is None:
        reactions = model.reactions
    reaction_chunks = (chunk for chunk in partition(reactions, len(view)))
    func_obj = FvaFunctionObject(model, _flux_variability_analysis)
    chunky_results = view.map(func_obj, reaction_chunks)
    solution = {}
    for dictionary in chunky_results:
        solution.update(dictionary)
    return solution


def _production_envelope_inner(model, point, variable_reactions):
    for (reaction, coordinate) in zip(variable_reactions, point):
        reaction.lower_bound = coordinate
        reaction.upper_bound = coordinate
    interval = []
    model.objective.direction = 'min'
    solution = model.optimize()
    if solution.status == 'optimal':
        interval.append(solution.f)
    else:
        model.solver.configuration.presolver = True
        model.solver.configuration.verbosity = 3
        solution = model.optimize()
        if solution.status == 'optimal':
            interval.append(solution.f)
        else:
            interval.append(solution.status)
    model.objective.direction = 'max'
    solution = model.optimize()
    if solution.status == 'optimal':
        interval.append(solution.f)
    else:
        model.solver.configuration.presolver = True
        model.solver.configuration.verbosity = 3
        solution = model.optimize()
        if solution.status == 'optimal':
            interval.append(solution.f)
        else:
            interval.append(solution.status)
    return interval


def production_envelope(model, target, variables=[], points=20, view=None):
    """Calculate a production envelope ..."""
    tm = TimeMachine()
    original_objective = copy(model.objective)
    init_bookmark = tm(do=partial(setattr, model, 'objective', target),
                       undo=partial(setattr, model, 'objective', original_objective))
    variable_reactions = _ids_to_reactions(model, variables)
    variables_min_max = flux_variability_analysis(
        model, reactions=variable_reactions)
    grid = [numpy.linspace(val['minimum'], val['maximum'], points, endpoint=True)
            for key, val in variables_min_max.iteritems()]
    generator = itertools.product(*grid)
    original_bounds = dict([(reaction, (reaction.lower_bound, reaction.upper_bound))
                            for reaction in variable_reactions])
    envelope = OrderedDict()
    for point in generator:
        envelope[point] = _production_envelope_inner(
            model, point, variable_reactions)

    for reaction, bounds in original_bounds.iteritems():
        reaction.lower_bound = bounds[0]
        reaction.upper_bound = bounds[1]
    tm.undo(init_bookmark)
    final_envelope = OrderedDict()
    for i, reaction in enumerate(variable_reactions):
        final_envelope[reaction.id] = [elem[i] for elem in envelope.keys()]
    final_envelope[target] = envelope.values()
    return final_envelope


def phenotypic_phase_plane(model, target, variables=[], points=20, view=None):
    """Calculate a production envelope ..."""
    tm = TimeMachine()
    original_objective = copy(model.objective)
    init_bookmark = tm(do=partial(setattr, model, 'objective', target),
                       undo=partial(setattr, model, 'objective', original_objective))
    variable_reactions = _ids_to_reactions(model, variables)
    variables_min_max = flux_variability_analysis(
        model, reactions=variable_reactions)
    grid = [numpy.linspace(val['minimum'], val['maximum'], points, endpoint=True)
            for key, val in variables_min_max.iteritems()]
    generator = itertools.product(*grid)
    original_bounds = dict([(reaction, (reaction.lower_bound, reaction.upper_bound))
                            for reaction in variable_reactions])
    envelope = OrderedDict()
    for point in generator:
        envelope[point] = _production_envelope_inner(
            model, point, variable_reactions)

    for reaction, bounds in original_bounds.iteritems():
        reaction.lower_bound = bounds[0]
        reaction.upper_bound = bounds[1]
    tm.undo(init_bookmark)
    final_envelope = OrderedDict()
    for i, reaction in enumerate(variable_reactions):
        final_envelope[reaction.id] = [elem[i] for elem in envelope.keys()]
    final_envelope[target] = envelope.values()
    return final_envelope


def cycle_free_fva(model, reactions=None, sloppy=False):
    """Flux-variability analysis."""
    original_objective = copy(model.objective)
    if reactions is None:
        reactions = model.reactions
    else:
        reactions = _ids_to_reactions(model, reactions)
    fva_sol = OrderedDict()
    for reaction in reactions:
        fva_sol[reaction.id] = dict()
        model.objective = reaction
        model.objective.direction = 'min'
        solution = model.optimize()
        if solution.status == 'optimal':
            bound = solution.f
            if sloppy and bound > -100:
                fva_sol[reaction.id]['minimum'] = bound
            else:
                v0_fluxes = solution.x_dict
                v1_cycle_free_fluxes = cycle_free_flux(model, v0_fluxes)
                if abs(v1_cycle_free_fluxes[reaction.id] - bound) < 10 ** -6:
                    fva_sol[reaction.id]['minimum'] = bound
                else:
                    v2_one_cycle_fluxes = cycle_free_flux(model, v0_fluxes, fix=[reaction.id])
                    zero_in_v1_but_not_in_v2 = list()
                    model_broken_cyle = copy(model)
                    for key, v1_flux in v1_cycle_free_fluxes.iteritems():
                        if v1_flux == 0 and v2_one_cycle_fluxes[key] != 0:
                            knockout_reaction = model_broken_cyle.reactions.get_by_id(key)
                            knockout_reaction.lower_bound = 0.
                            knockout_reaction.upper_bound = 0.
                    model_broken_cyle.objective.direction = 'min'
                    solution = model_broken_cyle.optimize()
                    if solution.status == 'optimal':
                        fva_sol[reaction.id]['minimum'] = solution.f
        else:
            fva_sol[reaction.id]['minimum'] = model.solution.status
    for reaction in reactions:
        model.objective = reaction
        model.objective.direction = 'max'
        solution = model.optimize()
        if solution.status == 'optimal':
            bound = solution.f
            if sloppy and bound < 100:
                fva_sol[reaction.id]['maximum'] = bound
            else:
                v0_fluxes = solution.x_dict
                v1_cycle_free_fluxes = cycle_free_flux(model, v0_fluxes)

                if abs(v1_cycle_free_fluxes[reaction.id] - bound) < 10 ** -6:
                    fva_sol[reaction.id]['maximum'] = v0_fluxes[reaction.id]
                else:
                    v2_one_cycle_fluxes = cycle_free_flux(model, v0_fluxes, fix=[reaction.id])
                    model_broken_cyle = copy(model)
                    for key, v1_flux in v1_cycle_free_fluxes.iteritems():
                        if v1_flux == 0 and v2_one_cycle_fluxes[key] != 0:
                            knockout_reaction = model_broken_cyle.reactions.get_by_id(key)
                            knockout_reaction.lower_bound = 0.
                            knockout_reaction.upper_bound = 0.
                    model_broken_cyle.objective.direction = 'max'
                    solution = model_broken_cyle.optimize()
                    if solution.status == 'optimal':
                        fva_sol[reaction.id]['maximum'] = solution.f
        else:
            fva_sol[reaction.id]['maximum'] = model.solution.status
    model.objective = original_objective
    return fva_sol


def cycle_free_flux(model, fluxes, fix=[]):
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
            obj_terms.append(Mul._from_args([S.One, internal_reaction.variable]))
            tm(do=partial(setattr, internal_reaction, 'lower_bound', 0),
               undo=partial(setattr, internal_reaction, 'lower_bound', internal_reaction.lower_bound))
            tm(do=partial(setattr, internal_reaction, 'upper_bound', internal_flux),
               undo=partial(setattr, internal_reaction, 'upper_bound', internal_reaction.upper_bound))
        elif internal_flux < 0:
            obj_terms.append(Mul._from_args([S.NegativeOne, internal_reaction.variable]))
            tm(do=partial(setattr, internal_reaction, 'lower_bound', internal_flux),
               undo=partial(setattr, internal_reaction, 'lower_bound', internal_reaction.lower_bound))
            tm(do=partial(setattr, internal_reaction, 'upper_bound', 0),
               undo=partial(setattr, internal_reaction, 'upper_bound', internal_reaction.upper_bound))
        else:
            pass
            # print internal_flux, internal_reaction
    for reaction_id in fix:
        reaction_to_fix = model.reactions.get_by_id(reaction_id)
        tm(do=partial(setattr, reaction_to_fix, 'lower_bound', fluxes[reaction_id]),
           undo=partial(setattr, reaction_to_fix, 'lower_bound', reaction_to_fix.lower_bound))
        tm(do=partial(setattr, reaction_to_fix, 'upper_bound', fluxes[reaction_id]),
           undo=partial(setattr, reaction_to_fix, 'upper_bound', reaction_to_fix.upper_bound))
    tm(do=partial(setattr, model, 'objective',
                  Objective(Add._from_args(obj_terms), name='Flux minimization', direction='min', sloppy=True)),
       undo=partial(setattr, model, 'objective', model.objective))
    solution = model.optimize()
    tm.reset()
    return solution.x_dict