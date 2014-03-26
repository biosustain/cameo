# Copyright 2013 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


from copy import copy
import itertools
from collections import OrderedDict
from functools import partial
import numpy
from cobra.core import Reaction
from cameo.config import default_view
from cameo.util import partition, TimeMachine
from IPython.parallel import interactive


def trace_dead_routes(model, dead_reaction):
    """Explain a dead reaction by recursively tracing dead paths through the system."""
    for metabolite in dead_reaction.metabolites:
        demand_reaction = model.add_demand(metabolite)
        model.objective = demand_reaction
        solution = model.optimize()
        if solution.status == 'optimal':
            print solution.f
        model.remove_reactions([demand_reaction])


def _ids_to_reactions(model, reactions):
    """translates reaction ids into reactions (skips reactions)."""
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
    # tm = TimeMachine()
    original_objective = copy(model.objective)
    # tm(do=, undo=partial(setattr, model, 'objective', original_objective))
    if reactions is None:
        reactions = model.reactions
    else:
        reactions = _ids_to_reactions(model, reactions)
    fva_sol = OrderedDict()
    for reaction in reactions:
        fva_sol[reaction.id] = dict()
        # tm(do=partial(setattr, model, 'objective', reaction), undo=partial(setattr, model, 'objective', original_objective))
        model.objective = reaction
        model.objective.direction = 'min'
        solution = model.optimize()
        if solution.status == 'optimal':
            fva_sol[reaction.id]['minimum'] = model.solution.f
        else:
            fva_sol[reaction.id]['minimum'] = model.solution.status
        # tm.reset()
    for reaction in reactions:
        # tm(do=partial(setattr, model, 'objective', reaction), undo=partial(setattr, model, 'objective', original_objective))
        model.objective = reaction
        model.objective.direction = 'max'
        solution = model.optimize()
        if solution.status == 'optimal':
            fva_sol[reaction.id]['maximum'] = model.solution.f
        else:
            fva_sol[reaction.id]['maximum'] = model.solution.status
        # tm.reset()
    model.objective = original_objective
    return fva_sol


class FvaFunctionObject(object):

    def __init__(self, model):
        self.model = model

    def __call__(self, reactions):
        return _flux_variability_analysis(self.model, reactions)


def flux_variability_analysis(model, reactions=None, view=default_view):
    """Flux-variability analysis."""
    if reactions is None:
        reactions = model.reactions
    reaction_chunks = (chunk for chunk in partition(reactions, len(view)))

    func_obj = FvaFunctionObject(model)
    chunky_results = view.map(func_obj, reaction_chunks)
    # chunky_results = view.map(lambda reactions, model=model:_flux_variability_analysis(model, reactions), reaction_chunks)
    solution = {}
    for dictionary in chunky_results:
        solution.update(dictionary)
    return solution


def _production_envelope_inner(model, point, variable_reactions):
    envelope = OrderedDict()
    for (reaction, coordinate) in zip(variable_reactions, point):
        if reaction.lower_bound > coordinate:
            reaction.lower_bound = coordinate
            reaction.upper_bound = coordinate
        else:
            reaction.upper_bound = coordinate
            reaction.lower_bound = coordinate
        solution = model.optimize()
        if solution.status == 'optimal':
            return solution.f
        else:
            # print zip(variable_reactions, point), solution.status, solution.f
            model.solver.configuration.presolver = True
            model.solver.configuration.verbosity = 3
            solution = model.optimize()
            if solution.status == 'optimal':
                return solution.f
            else:
                return 0


def production_envelope2(model, target, variables=[], points=20, view=None):
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


def production_envelope(model, target, variables=[], points=20):
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
        for (reaction, coordinate) in zip(variable_reactions, point):
            if reaction.lower_bound > coordinate:
                reaction.lower_bound = coordinate
                reaction.upper_bound = coordinate
            else:
                reaction.upper_bound = coordinate
                reaction.lower_bound = coordinate

        solution = model.optimize()
        if solution.status == 'optimal':
            envelope[point] = solution.f
        else:
            # print zip(variable_reactions, point), solution.status, solution.f
            model.solver.configuration.presolver = True
            model.solver.configuration.verbosity = 3
            solution = model.optimize()
            if solution.status == 'optimal':
                envelope[point] = solution.f
            else:
                envelope[point] = 0
            model.solver.configuration.presolver = False
            model.solver.configuration.verbosity = 0
            # tm.undo(bookmark)
    for reaction, bounds in original_bounds.iteritems():
        reaction.lower_bound = bounds[0]
        reaction.upper_bound = bounds[1]
    tm.undo(init_bookmark)
    final_envelope = OrderedDict()
    for i, reaction in enumerate(variable_reactions):
        final_envelope[reaction.id] = [elem[i] for elem in envelope.keys()]
    final_envelope[target] = envelope.values()
    return final_envelope


if __name__ == '__main__':

    # repickled = pickle.loads(pickle.dumps(model))
    # print "#"*80
    # print model.solver.variables.keys()
    # print "#"*80
    # print repickled.solver.variables.keys()
    # print "#"*80

    client = parallel.Client()
    client.block = True
    view = client.direct_view()

    t1 = time.time()
    fva_sol_p = flux_variability_analysis_parallel(model, view=view)
    t2 = time.time()
    print "Execution time: %s" % (t2 - t1)
    # print fva_sol_p

    print "#" * 80
    # fva_sol = flux_variability_analysis(model, ['Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2',
    #     'EX_o2_LPAREN_e_RPAREN_'])
    t1 = time.time()
    fva_sol = flux_variability_analysis(model)
    t2 = time.time()
    print "Execution time: %s" % (t2 - t1)
    # print fva_sol

    for key, val in fva_sol.iteritems():
        val2 = fva_sol_p[key]
        # print key, val, val2
        print key, abs(val['minimum'] - val['minimum']), abs(val['maximum'] - val['maximum'])


    # envelope = production_envelope(model,
    # 'EX_ac_LPAREN_e_RPAREN_',
    # variables=['Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2',
    # 'EX_o2_LPAREN_e_RPAREN_'],
    # points=10
    # )
