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

from __future__ import absolute_import, print_function

__all__ = ['find_blocked_reactions', 'flux_variability_analysis', 'phenotypic_phase_plane', 'fbid']

from cobra.core import Reaction

import six
from six.moves import zip

import itertools
from copy import copy
from collections import OrderedDict
from functools import partial, reduce
import numpy
import pandas

import cameo
from cameo import config
from cameo.exceptions import Infeasible, Unbounded, SolveError
from cameo.util import TimeMachine, partition
from cameo.parallel import SequentialView
from cameo.flux_analysis.simulation import _cycle_free_flux
from cameo.core.result import Result
from cameo.ui import notice
from cameo.visualization import plotting

import logging

logger = logging.getLogger(__name__)


def find_blocked_reactions(model):
    """Determine reactions that cannot carry steady-state flux.

    Parameters
    ----------
    model: SolverBasedModel

    Returns
    -------
    list
        A list of reactions.

    """
    with TimeMachine() as tm:
        for exchange in model.exchanges:
            tm(do=partial(setattr, exchange, 'lower_bound', -999999),
               undo=partial(setattr, exchange, 'lower_bound', exchange.lower_bound))
            tm(do=partial(setattr, exchange, 'upper_bound', 999999),
               undo=partial(setattr, exchange, 'upper_bound', exchange.upper_bound))
        fva_solution = flux_variability_analysis(model)
    return [model.reactions.get_by_id(id) for id in
            fva_solution.data_frame.query('upper_bound == lower_bound == 0').index]


def flux_variability_analysis(model, reactions=None, fraction_of_optimum=0., remove_cycles=False, view=None):
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
    if view is None:
        view = config.default_view
    if reactions is None:
        reactions = model.reactions
    with TimeMachine() as tm:
        if fraction_of_optimum > 0.:
            try:
                obj_val = model.solve().f
            except SolveError as e:
                logger.debug(
                    "flux_variability_analyis was not able to determine an optimal solution for objective %s" % model.objective)
                raise e
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
    return FluxVariabilityResult(solution)


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
    PhenotypicPhasePlaneResult
        The phenotypic phase plane.

    """
    if isinstance(variables, str):
        variables = [variables]
    elif isinstance(variables, cameo.core.reaction.Reaction):
        variables = [variables]

    if view is None:
        view = config.default_view
    with TimeMachine() as tm:
        if objective is not None:
            tm(do=partial(setattr, model, 'objective', objective),
               undo=partial(setattr, model, 'objective', model.objective))

        variable_reactions = model._ids_to_reactions(variables)
        variables_min_max = flux_variability_analysis(model, reactions=variable_reactions, view=SequentialView())
        grid = [numpy.linspace(lower_bound, upper_bound, points, endpoint=True) for
                reaction_id, lower_bound, upper_bound in
                variables_min_max.data_frame.itertuples()]
        grid_generator = itertools.product(*grid)
        original_bounds = dict([(reaction, (reaction.lower_bound, reaction.upper_bound))
                                for reaction in variable_reactions])

        chunks_of_points = partition(list(grid_generator), len(view))
        evaluator = _PhenotypicPhasePlaneChunkEvaluator(model, variable_reactions)
        chunk_results = view.map(evaluator, chunks_of_points)
        envelope = reduce(list.__add__, chunk_results)

        for reaction, bounds in six.iteritems(original_bounds):
            reaction.lower_bound = bounds[0]
            reaction.upper_bound = bounds[1]

        variable_reactions_ids = [reaction.id for reaction in variable_reactions]
        phase_plane = pandas.DataFrame(envelope, columns=(variable_reactions_ids +
                                                          ['objective_lower_bound', 'objective_upper_bound']))

        if objective is None:
            objective = model.objective

        if isinstance(objective, Reaction):
            objective = objective.id
        else:
            objective = str(objective)

        return PhenotypicPhasePlaneResult(phase_plane, variable_reactions_ids, objective)


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
        reactions = model._ids_to_reactions(reactions)
    fva_sol = OrderedDict()
    [lb_flag, ub_flag] = [False, False]
    for reaction in reactions:
        fva_sol[reaction.id] = dict()
        model.objective = reaction
        model.objective.direction = 'min'
        try:
            solution = model.solve()
            fva_sol[reaction.id]['lower_bound'] = solution.f
        except Unbounded:
            fva_sol[reaction.id]['lower_bound'] = -numpy.inf
        except Infeasible:
            lb_flag = True

        model.objective.direction = 'max'
        try:
            solution = model.solve()
            fva_sol[reaction.id]['upper_bound'] = solution.f
        except Unbounded:
            fva_sol[reaction.id]['upper_bound'] = numpy.inf
        except Infeasible:
            ub_flag = True

        if lb_flag is True and ub_flag is True:
            fva_sol[reaction.id]['lower_bound'] = 0
            fva_sol[reaction.id]['upper_bound'] = 0
            [lb_flag, ub_flag] = [False, False]
        elif lb_flag is True and ub_flag is False:
            fva_sol[reaction.id]['lower_bound'] = fva_sol[reaction.id]['upper_bound']
            lb_flag = False
        elif lb_flag is False and ub_flag is True:
            fva_sol[reaction.id]['upper_bound'] = fva_sol[reaction.id]['lower_bound']
            ub_flag = False

    model.objective = original_objective
    df = pandas.DataFrame.from_dict(fva_sol, orient='index')
    lb_higher_ub = df[df.lower_bound > df.upper_bound]
    try:  # this is an alternative solution to what I did above with flags
        assert ((
                lb_higher_ub.lower_bound - lb_higher_ub.upper_bound) < 1e-6).all()  # Assert that these cases really only numerical artifacts
    except AssertionError as e:
        logger.debug(list(zip(model.reactions, (lb_higher_ub.lower_bound - lb_higher_ub.upper_bound) < 1e-6)))
    df.lower_bound[lb_higher_ub.index] = df.upper_bound[lb_higher_ub.index]
    return df


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
    cycle_count = 0
    try:
        original_objective = copy(model.objective)
        if reactions is None:
            reactions = model.reactions
        else:
            reactions = model._ids_to_reactions(reactions)
        fva_sol = OrderedDict()
        for reaction in reactions:
            fva_sol[reaction.id] = dict()
            model.objective = reaction
            model.objective.direction = 'min'
            try:
                solution = model.solve()
            except Unbounded:
                fva_sol[reaction.id]['lower_bound'] = -numpy.inf
                continue
            except Infeasible:
                fva_sol[reaction.id]['lower_bound'] = 0
                continue
            except Exception as e:
                raise e
            bound = solution.f
            if sloppy and bound > -100:
                fva_sol[reaction.id]['lower_bound'] = bound
            else:
                v0_fluxes = solution.x_dict
                v1_cycle_free_fluxes = _cycle_free_flux(model, v0_fluxes)
                if abs(v1_cycle_free_fluxes[reaction.id] - bound) < 10 ** -6:
                    fva_sol[reaction.id]['lower_bound'] = bound
                else:
                    cycle_count += 1
                    v2_one_cycle_fluxes = _cycle_free_flux(model, v0_fluxes, fix=[reaction.id])
                    tm = TimeMachine()
                    for key, v1_flux in six.iteritems(v1_cycle_free_fluxes):
                        if v1_flux == 0 and v2_one_cycle_fluxes[key] != 0:
                            knockout_reaction = model.reactions.get_by_id(key)
                            tm(do=partial(setattr, knockout_reaction, 'lower_bound', 0.),
                               undo=partial(setattr, knockout_reaction, 'lower_bound', knockout_reaction.lower_bound))
                            tm(do=partial(setattr, knockout_reaction, 'upper_bound', 0.),
                               undo=partial(setattr, knockout_reaction, 'upper_bound', knockout_reaction.lower_bound))
                    model.objective.direction = 'min'
                    try:
                        solution = model.solve()
                    except Unbounded:
                        fva_sol[reaction.id]['lower_bound'] = -numpy.inf
                    except Infeasible:
                        fva_sol[reaction.id]['lower_bound'] = 0
                    else:
                        fva_sol[reaction.id]['lower_bound'] = solution.f
                    finally:
                        tm.reset()
        for reaction in reactions:
            model.objective = reaction
            model.objective.direction = 'max'
            try:
                solution = model.solve()
            except Unbounded:
                fva_sol[reaction.id]['upper_bound'] = numpy.inf
                continue
            except Infeasible:
                fva_sol[reaction.id]['upper_bound'] = 0
                continue
            except Exception as e:
                raise e
            else:
                bound = solution.f
                if sloppy and bound < 100:
                    fva_sol[reaction.id]['upper_bound'] = bound
                else:
                    v0_fluxes = solution.x_dict
                    v1_cycle_free_fluxes = _cycle_free_flux(model, v0_fluxes)

                    if abs(v1_cycle_free_fluxes[reaction.id] - bound) < 10 ** -6:
                        fva_sol[reaction.id]['upper_bound'] = v0_fluxes[reaction.id]
                    else:
                        cycle_count += 1
                        v2_one_cycle_fluxes = _cycle_free_flux(model, v0_fluxes, fix=[reaction.id])
                        tm = TimeMachine()
                        for key, v1_flux in six.iteritems(v1_cycle_free_fluxes):
                            if v1_flux == 0 and v2_one_cycle_fluxes[key] != 0:
                                knockout_reaction = model.reactions.get_by_id(key)
                                tm(do=partial(setattr, knockout_reaction, 'lower_bound', 0.),
                                   undo=partial(setattr, knockout_reaction, 'lower_bound',
                                                knockout_reaction.lower_bound))
                                tm(do=partial(setattr, knockout_reaction, 'upper_bound', 0.),
                                   undo=partial(setattr, knockout_reaction, 'upper_bound',
                                                knockout_reaction.lower_bound))
                        model.objective.direction = 'max'
                        try:
                            solution = model.solve()
                        except Unbounded:
                            fva_sol[reaction.id]['upper_bound'] = numpy.inf
                        except Infeasible:
                            fva_sol[reaction.id]['upper_bound'] = 0
                        else:
                            fva_sol[reaction.id]['upper_bound'] = solution.f
                        finally:
                            tm.reset()
        df = pandas.DataFrame.from_dict(fva_sol, orient='index')
        lb_higher_ub = df[df.lower_bound > df.upper_bound]
        assert ((
                lb_higher_ub.lower_bound - lb_higher_ub.upper_bound) < 1e-6).all()  # Assert that these cases really only numerical artifacts
        df.lower_bound[lb_higher_ub.index] = df.upper_bound[lb_higher_ub.index]
    finally:
        model.objective = original_objective
    return df


class _PhenotypicPhasePlaneChunkEvaluator(object):
    def __init__(self, model, variable_reactions):
        self.model = model
        self.variable_reactions = variable_reactions

    def __call__(self, points):
        return [self._production_envelope_inner(point) for point in points]

    def _production_envelope_inner(self, point):
        tm = TimeMachine()
        try:
            for (reaction, coordinate) in zip(self.variable_reactions, point):
                tm(do=partial(setattr, reaction, 'lower_bound', coordinate),
                   undo=partial(setattr, reaction, 'lower_bound', reaction.lower_bound))
                tm(do=partial(setattr, reaction, 'upper_bound', coordinate),
                   undo=partial(setattr, reaction, 'upper_bound', reaction.upper_bound))
            interval = []
            tm(do=int, undo=partial(setattr, self.model.objective, 'direction', self.model.objective.direction))
            self.model.objective.direction = 'min'
            try:
                solution = self.model.solve().f
            except Infeasible:
                solution = 0
            interval.append(solution)
            self.model.objective.direction = 'max'
            try:
                solution = self.model.solve().f
            except Infeasible:
                solution = 0
            interval.append(solution)
        finally:
            tm.reset()
        return point + tuple(interval)


def fbid(model, knockouts, view=config.default_view, method="fva"):
    """
    Flux balance impact degree by Zhao et al 2013

    :param model: wild-type model
    :param knockouts: list of reaction knockouts
    :param method: the method to compute the perturbation. default is "fva" - Flux Variability Analysis.
        It can also be computed with "em" - Elementary modes
    :return: perturbation
    """

    if method == "fva":
        _fbid_fva(model, knockouts, view)
    elif method == "em":
        raise NotImplementedError("Elementary modes approach is not implemented")
    else:
        raise ValueError("%s method is not valid to compute FBIP" % method)


def _fbid_fva(model, knockouts, view):
    tm = TimeMachine()
    for reaction in model.reactions:
        if reaction.reversibility:
            tm(do=partial(setattr, reaction, 'lower_bound', -1),
               undo=partial(setattr, reaction, 'lower_bound', reaction.lower_bound))
            tm(do=partial(setattr, reaction, 'upper_bound', 1),
               undo=partial(setattr, reaction, 'upper_bound', reaction.upper_bound))
        else:
            tm(do=partial(setattr, reaction, 'lower_bound', 0),
               undo=partial(setattr, reaction, 'lower_bound', reaction.lower_bound))
            tm(do=partial(setattr, reaction, 'upper_bound', 1),
               undo=partial(setattr, reaction, 'upper_bound', reaction.upper_bound))

    wt_fva = flux_variability_analysis(model, view)
    for reaction in knockouts:
        tm(do=partial(setattr, reaction, 'upper_bound', 0),
           undo=partial(setattr, reaction, 'upper_bound', reaction.upper_bound))
        tm(do=partial(setattr, reaction, 'lower_bound', 0),
           undo=partial(setattr, reaction, 'lower_bound', reaction.lower_bound))

    mt_fva = flux_variability_analysis(model, view)

    perturbation = 0
    for reaction in model.reactions:
        if wt_fva[reaction.id] != 0 and mt_fva[reaction.id] == 0:
            perturbation += 1

    tm.reset()
    return perturbation


class PhenotypicPhasePlaneResult(Result):
    def __init__(self, phase_plane, variable_ids, objective, *args, **kwargs):
        super(PhenotypicPhasePlaneResult, self).__init__(*args, **kwargs)
        self._phase_plane = phase_plane
        self.variable_ids = variable_ids
        self.objective = objective

    @property
    def data_frame(self):
        return pandas.DataFrame(self._phase_plane)

    def plot(self, grid=None, width=None, height=None, title=None, axis_font_size=None):
        if len(self.variable_ids) > 1:
            notice("Multi-dimensional plotting is not supported")
            return
        plotting.plot_production_envelope(self._phase_plane, objective=self.objective, key=self.variable_ids[0],
                                          grid=grid, width=width, height=height, title=title,
                                          axis_font_size=axis_font_size)

    def __getitem__(self, item):
        return self._phase_plane[item]

    def iterrows(self):
        return self._phase_plane.iterrows()


class FluxVariabilityResult(Result):
    def __init__(self, data_frame, *args, **kwargs):
        super(FluxVariabilityResult, self).__init__(*args, **kwargs)
        self._data_frame = data_frame

    @property
    def data_frame(self):
        return self._data_frame

    def plot(self, grid=None, width=None, height=None, title=None, axis_font_size=None):
        raise NotImplementedError('Plotting of flux variability results has not been implemented yet.')

    def __getitem__(self, item):
        return self._data_frame[item]

    def iterrows(self):
        return self._data_frame.iterrows()
