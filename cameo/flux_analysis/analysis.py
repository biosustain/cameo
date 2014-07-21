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


import itertools
from copy import copy
from collections import OrderedDict
from functools import partial
import numpy
from cobra.core import Reaction
from cameo import config
from cameo.exceptions import UndefinedSolution, Infeasible
from cameo.util import TimeMachine, partition
from cameo.parallel import SequentialView
from cameo.flux_analysis.simulation import _cycle_free_flux
import pandas


def flux_variability_analysis(model, reactions=None, fraction_of_optimum=0., remove_cycles=True, view=None):
    """Flux variability analysis.

    Parameters
    ----------
    model: SolverBasedModel
    reactions: None or iterable
        The list of reaction whose lower and upper bounds should be determined.
        If `None`, all reactions in `model` will be assessed.
    view: SequentialView or MultiprocessingView or ipython.cluster.DirectView
        A parallelization view.

    Returns
    -------
    pandas.DataFrame
        Pandas DataFrame containing the results of the flux variability analysis.

    """
    tm = TimeMachine()
    if view is None:
        view = config.default_view
    if reactions is None:
        reactions = model.reactions
    try:
        if fraction_of_optimum > 0.:
            if model.objective.direction == 'max':
                fix_obj_constraint = model.solver.interface.Constraint(model.objective.expression,
                                                                       lb=fraction_of_optimum * obj_val)
            else:
                fix_obj_constraint = model.solver.interface.Constraint(model.objective.expression,
                                                                       ub=fraction_of_optimum * obj_val)
            tm(do=partial(model.solver._add_constraint, fix_obj_constraint),
               undo=partial(model.solver._remove_constraint, fix_obj_constraint))
        reaction_chunks = (chunk for chunk in partition(reactions, len(view)))
        if remove_cycles == True:
            func_obj = _FvaFunctionObject(model, _cycle_free_fva)
        else:
            func_obj = _FvaFunctionObject(model, _flux_variability_analysis)
        chunky_results = view.map(func_obj, reaction_chunks)
        solution = pandas.concat(chunky_results)
    finally:
        tm.reset()
    return solution


def phenotypic_phase_plane(model, variables=[], objective=None, points=20, view=None):
    """Phenotypic phase plane analysis.

    Parameters
    ----------
    model: SolverBasedModel
    variables: str or reaction or iterable
        A reaction ID, reaction, or list of reactions to be varied.
    objective: str or reaction or optlang.Objective
        An objective to be minimized/maximized for
    points: int or iterable
        Number of points to be interspersed between the variable bounds.
        A list of same same dimensions as `variables` can be used to specify
        variable specific numbers of points.
    view: SequentialView or MultiprocessingView or ipython.cluster.DirectView
        A parallelization view.

    Returns
    -------
    pandas.DataFrame
        Pandas DataFrame containing the phenotypic phase plane.

    """
    if view is None:
        view = SequentialView()
    tm = TimeMachine()
    original_objective = copy(model.objective)
    if objective is not None:
        tm(do=partial(setattr, 'objective', objective), undo=partial(setattr, 'objective', model.objective))

    variable_reactions = _ids_to_reactions(model, variables)
    variables_min_max = flux_variability_analysis(model, reactions=variable_reactions)
    grid = [numpy.linspace(lower_bound, upper_bound, points, endpoint=True) for reaction_id, lower_bound, upper_bound in
            variables_min_max.itertuples()]
    grid_generator = itertools.product(*grid)
    original_bounds = dict([(reaction, (reaction.lower_bound, reaction.upper_bound))
                            for reaction in variable_reactions])

    chunks_of_points = partition(list(grid_generator), len(view))
    evaluator = _PhenotypicPhasePlaneChunkEvaluator(model, variable_reactions)
    chunk_results = view.map(evaluator, chunks_of_points)
    envelope = reduce(list.__add__, chunk_results)

    # for point in grid_generator:
    # envelope.append(point + _production_envelope_inner(model, point, variable_reactions))

    for reaction, bounds in original_bounds.iteritems():
        reaction.lower_bound = bounds[0]
        reaction.upper_bound = bounds[1]

    variable_reactions_ids = [reaction.id for reaction in variable_reactions]
    return pandas.DataFrame(envelope,
                            columns=(variable_reactions_ids + ['objective_lower_bound', 'objective_upper_bound']))


def _ids_to_reactions(model, reactions):
    """Translate reaction IDs into reactions (skips reactions)."""
    clean_reactions = list()
    for reaction in reactions:
        if isinstance(reaction, str):
            clean_reactions.append(model.reactions.get_by_id(reaction))
        elif isinstance(reaction, Reaction):
            clean_reactions.append(reaction)
        else:
            raise Exception('%s is not a reaction or reaction ID.' % reaction)
    return clean_reactions


class _FvaFunctionObject(object):
    def __init__(self, model, fva):
        self.model = model
        self.fva = fva

    def __call__(self, reactions):
        return self.fva(self.model, reactions)

def _flux_variability_analysis(model, reactions=None):
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
            fva_sol[reaction.id]['lower_bound'] = model.solution.f
        else:
            fva_sol[reaction.id]['lower_bound'] = model.solution.status
    for reaction in reactions:
        model.objective = reaction
        model.objective.direction = 'max'
        solution = model.optimize()
        if solution.status == 'optimal':
            fva_sol[reaction.id]['upper_bound'] = model.solution.f
        else:
            fva_sol[reaction.id]['upper_bound'] = model.solution.status
    model.objective = original_objective
    return pandas.DataFrame.from_dict(fva_sol, orient='index')


def _cycle_free_fva(model, reactions=None, sloppy=True):
    """Cycle free flux-variability analysis. (http://cran.r-project.org/web/packages/sybilcycleFreeFlux/index.html)

    Parameters
    ----------
    model: SolverBasedModel
    reactions: list
        List of reactions whose flux-ranges should be determined.
    sloppy: boolean
        If true, only abs(v) > 100 are checked to be futile cycles.

    Rer
    """
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
                fva_sol[reaction.id]['lower_bound'] = bound
            else:
                v0_fluxes = solution.x_dict
                v1_cycle_free_fluxes = _cycle_free_flux(model, v0_fluxes)
                if abs(v1_cycle_free_fluxes[reaction.id] - bound) < 10 ** -6:
                    fva_sol[reaction.id]['lower_bound'] = bound
                else:
                    v2_one_cycle_fluxes = _cycle_free_flux(model, v0_fluxes, fix=[reaction.id])
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
                        fva_sol[reaction.id]['lower_bound'] = solution.f
        else:
            fva_sol[reaction.id]['lower_bound'] = model.solution.status
    for reaction in reactions:
        model.objective = reaction
        model.objective.direction = 'max'
        solution = model.optimize()
        if solution.status == 'optimal':
            bound = solution.f
            if sloppy and bound < 100:
                fva_sol[reaction.id]['upper_bound'] = bound
            else:
                v0_fluxes = solution.x_dict
                v1_cycle_free_fluxes = _cycle_free_flux(model, v0_fluxes)

                if abs(v1_cycle_free_fluxes[reaction.id] - bound) < 10 ** -6:
                    fva_sol[reaction.id]['upper_bound'] = v0_fluxes[reaction.id]
                else:
                    v2_one_cycle_fluxes = _cycle_free_flux(model, v0_fluxes, fix=[reaction.id])
                    model_broken_cyle = copy(model)
                    for key, v1_flux in v1_cycle_free_fluxes.iteritems():
                        if v1_flux == 0 and v2_one_cycle_fluxes[key] != 0:
                            knockout_reaction = model_broken_cyle.reactions.get_by_id(key)
                            knockout_reaction.lower_bound = 0.
                            knockout_reaction.upper_bound = 0.
                    model_broken_cyle.objective.direction = 'max'
                    solution = model_broken_cyle.optimize()
                    if solution.status == 'optimal':
                        fva_sol[reaction.id]['upper_bound'] = solution.f
        else:
            fva_sol[reaction.id]['upper_bound'] = model.solution.status
    model.objective = original_objective
    fva_sol = pandas.DataFrame.from_dict(fva_sol, orient='index')
    return fva_sol


class _PhenotypicPhasePlaneChunkEvaluator(object):
    def __init__(self, model, variable_reactions):
        self.model = model
        self.variable_reactions = variable_reactions

    def __call__(self, points):
        return [self._production_envelope_inner(point) for point in points]

    def _production_envelope_inner(self, point):
        for (reaction, coordinate) in zip(self.variable_reactions, point):
            reaction.lower_bound, reaction.upper_bound = coordinate, coordinate
        interval = []
        self.model.objective.direction = 'min'
        try:
            solution = self.model.solve().f
        except (Infeasible,
                UndefinedSolution):  # TODO: should really just be Infeasible (see http://lists.gnu.org/archive/html/help-glpk/2013-09/msg00015.html)
            solution = 0
        interval.append(solution)
        self.model.objective.direction = 'max'
        try:
            solution = self.model.optimize().f
        except (Infeasible, UndefinedSolution):
            solution = 0
        interval.append(solution)
        return point + tuple(interval)


if __name__ == '__main__':
    import time
    from cameo import load_model
    from cameo.parallel import MultiprocessingView

    model = load_model('../../tests/data/EcoliCore.xml')
    # model.solver = 'cplex'
    view = MultiprocessingView()
    tic = time.time()
    ppp = phenotypic_phase_plane(model,
                                 ['EX_o2_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_', 'EX_nh4_LPAREN_e_RPAREN_'],
                                 view=view, points=30)
    # print ppp
    # print ppp.describe()
    print time.time() - tic

    view = SequentialView()
    tic = time.time()
    ppp = phenotypic_phase_plane(model,
                                 ['EX_o2_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_', 'EX_nh4_LPAREN_e_RPAREN_'],
                                 view=view, points=30)
    print time.time() - tic