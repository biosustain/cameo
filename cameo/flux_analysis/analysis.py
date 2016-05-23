# coding=utf-8
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

from cobra.core import Reaction, Metabolite

import six
from six.moves import zip
from numpy import trapz

import itertools
from collections import OrderedDict
from functools import partial, reduce
import numpy
import pandas

import cameo
from cameo import config
from cameo.exceptions import Infeasible, Unbounded, SolveError, UndefinedSolution
from cameo.util import TimeMachine, partition
from cameo.parallel import SequentialView
from cameo.core.result import Result
from cameo.ui import notice
from cameo.visualization.plotting import plotter
from cameo.flux_analysis.util import remove_infeasible_cycles

import logging

logger = logging.getLogger(__name__)

__all__ = ['find_blocked_reactions', 'flux_variability_analysis', 'phenotypic_phase_plane',
           'flux_balance_impact_degree']


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
            tm(do=partial(setattr, exchange, 'lower_bound', -9999),
               undo=partial(setattr, exchange, 'lower_bound', exchange.lower_bound))
            tm(do=partial(setattr, exchange, 'upper_bound', 9999),
               undo=partial(setattr, exchange, 'upper_bound', exchange.upper_bound))
        fva_solution = flux_variability_analysis(model)
    return [reaction for reaction in model.reactions
            if round(fva_solution.data_frame.loc[reaction.id, "lower_bound"], config.ndecimals) == 0 and
            round(fva_solution.data_frame.loc[reaction.id, "upper_bound"], config.ndecimals) == 0]


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
            model.fix_objective_as_constraint(fraction=fraction_of_optimum, time_machine=tm)
        tm(do=int, undo=partial(setattr, model, "objective", model.objective))
        reaction_chunks = (chunk for chunk in partition(reactions, len(view)))
        if remove_cycles:
            func_obj = _FvaFunctionObject(model, _cycle_free_fva)
        else:
            func_obj = _FvaFunctionObject(model, _flux_variability_analysis)
        chunky_results = view.map(func_obj, reaction_chunks)
        solution = pandas.concat(chunky_results)
    return FluxVariabilityResult(solution)


def phenotypic_phase_plane(model, variables=[], objective=None, points=20, view=None):
    """Phenotypic phase plane analysis [1].

    Parameters
    ----------
    model: SolverBasedModel
    variables: str or reaction or iterable
        A reaction ID, reaction, or list of reactions to be varied.
    objective: str or reaction or optlang.Objective or Metabolite, optional
        An objective, a reaction's flux, or a metabolite's production to be minimized/maximized
        (defaults to the current model objective).
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

    References
    ----------
    [1] Edwards, J. S., Ramakrishna, R. and Palsson, B. O. (2002). Characterizing the metabolic phenotype: a phenotype
        phase plane analysis. Biotechnology and Bioengineering, 77(1), 27–36. doi:10.1002/bit.10047
    """
    if isinstance(variables, str):
        variables = [variables]
    elif isinstance(variables, cameo.core.reaction.Reaction):
        variables = [variables]

    if view is None:
        view = config.default_view
    with TimeMachine() as tm:
        if objective is not None:
            try:
                objective = model.reaction_for(objective, time_machine=tm)
            except KeyError:
                pass

            tm(do=partial(setattr, model, 'objective', objective),
               undo=partial(setattr, model, 'objective', model.objective))

        variable_reactions = model._ids_to_reactions(variables)
        variables_min_max = flux_variability_analysis(model, reactions=variable_reactions, view=SequentialView())
        grid = [numpy.linspace(lower_bound, upper_bound, points, endpoint=True) for
                reaction_id, lower_bound, upper_bound in
                variables_min_max.data_frame.itertuples()]
        grid_generator = itertools.product(*grid)

        chunks_of_points = partition(list(grid_generator), len(view))
        evaluator = _PhenotypicPhasePlaneChunkEvaluator(model, variable_reactions)
        chunk_results = view.map(evaluator, chunks_of_points)
        envelope = reduce(list.__add__, chunk_results)

    variable_reactions_ids = []
    nice_variable_ids = []
    for reaction in variable_reactions:
        if hasattr(reaction, "nice_id"):
            variable_reactions_ids.append(reaction.id)
            nice_variable_ids.append(reaction.nice_id)
        else:
            variable_reactions_ids.append(reaction.id)
            nice_variable_ids.append(reaction.id)

    phase_plane = pandas.DataFrame(envelope, columns=(variable_reactions_ids +
                                                      ['objective_lower_bound', 'objective_upper_bound']))

    if objective is None:
        objective = model.objective

    if isinstance(objective, Reaction):
        if hasattr(objective, 'nice_id'):
            nice_objective_id = objective.nice_id
            objective = objective.id
        else:
            objective = objective.id
            nice_objective_id = objective
    else:
        objective = str(objective)
        nice_objective_id = str(objective)

    return PhenotypicPhasePlaneResult(phase_plane, variable_reactions_ids, objective,
                                      nice_variable_ids=nice_variable_ids, nice_objective_id=nice_objective_id)


class _FvaFunctionObject(object):
    def __init__(self, model, fva):
        self.model = model
        self.fva = fva

    def __call__(self, reactions):
        return self.fva(self.model, reactions)


def _flux_variability_analysis(model, reactions=None):
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

    df = pandas.DataFrame.from_dict(fva_sol, orient='index')
    lb_higher_ub = df[df.lower_bound > df.upper_bound]
    try:  # this is an alternative solution to what I did above with flags
        assert ((
                lb_higher_ub.lower_bound - lb_higher_ub.upper_bound) < 1e-6).all()  # Assert that these cases really only numerical artifacts
    except AssertionError:
        logger.debug(list(zip(model.reactions, (lb_higher_ub.lower_bound - lb_higher_ub.upper_bound) < 1e-6)))
    df.lower_bound[lb_higher_ub.index] = df.upper_bound[lb_higher_ub.index]
    return df


def _cycle_free_fva(model, reactions=None, sloppy=True, sloppy_bound=666):
    """Cycle free flux-variability analysis. (http://cran.r-project.org/web/packages/sybilcycleFreeFlux/index.html)

    Parameters
    ----------
    model : SolverBasedModel
    reactions : list
        List of reactions whose flux-ranges should be determined.
    sloppy : boolean, optional
        If true, only fluxes v with abs(v) > sloppy_bound are checked to be futile cycles (defaults to True).
    sloppy_bound : int, optional
        The threshold bound used by sloppy (defaults to the number of the beast).
    """
    cycle_count = 0
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
        bound = solution.f
        if sloppy and abs(bound) < sloppy_bound:
            fva_sol[reaction.id]['lower_bound'] = bound
        else:
            logger.debug('Determine if {} with bound {} is a cycle'.format(reaction.id, bound))
            v0_fluxes = solution.x_dict
            v1_cycle_free_fluxes = remove_infeasible_cycles(model, v0_fluxes)
            if abs(v1_cycle_free_fluxes[reaction.id] - bound) < 10 ** -6:
                fva_sol[reaction.id]['lower_bound'] = bound
            else:
                logger.debug('Cycle detected: {}'.format(reaction.id))
                cycle_count += 1
                v2_one_cycle_fluxes = remove_infeasible_cycles(model, v0_fluxes, fix=[reaction.id])
                with TimeMachine() as tm:
                    for key, v1_flux in six.iteritems(v1_cycle_free_fluxes):
                        if round(v1_flux, config.ndecimals) == 0 and round(v2_one_cycle_fluxes[key],
                                                                           config.ndecimals) != 0:
                            knockout_reaction = model.reactions.get_by_id(key)
                            tm(do=partial(setattr, knockout_reaction, 'lower_bound', 0.),
                               undo=partial(setattr, knockout_reaction, 'lower_bound', knockout_reaction.lower_bound))
                            tm(do=partial(setattr, knockout_reaction, 'upper_bound', 0.),
                               undo=partial(setattr, knockout_reaction, 'upper_bound', knockout_reaction.upper_bound))
                    model.objective.direction = 'min'
                    try:
                        solution = model.solve()
                    except Unbounded:
                        fva_sol[reaction.id]['lower_bound'] = -numpy.inf
                    except Infeasible:
                        fva_sol[reaction.id]['lower_bound'] = 0
                    else:
                        fva_sol[reaction.id]['lower_bound'] = solution.f

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
        bound = solution.f
        if sloppy and abs(bound) < sloppy_bound:
            fva_sol[reaction.id]['upper_bound'] = bound
        else:
            logger.debug('Determine if {} with bound {} is a cycle'.format(reaction.id, bound))
            v0_fluxes = solution.x_dict
            v1_cycle_free_fluxes = remove_infeasible_cycles(model, v0_fluxes)
            if abs(v1_cycle_free_fluxes[reaction.id] - bound) < 1e-6:
                fva_sol[reaction.id]['upper_bound'] = v0_fluxes[reaction.id]
            else:
                logger.debug('Cycle detected: {}'.format(reaction.id))
                cycle_count += 1
                v2_one_cycle_fluxes = remove_infeasible_cycles(model, v0_fluxes, fix=[reaction.id])
                with TimeMachine() as tm:
                    for key, v1_flux in six.iteritems(v1_cycle_free_fluxes):
                        if round(v1_flux, config.ndecimals) == 0 and round(v2_one_cycle_fluxes[key],
                                                                           config.ndecimals) != 0:
                            knockout_reaction = model.reactions.get_by_id(key)
                            tm(do=partial(setattr, knockout_reaction, 'lower_bound', 0.),
                               undo=partial(setattr, knockout_reaction, 'lower_bound',
                                            knockout_reaction.lower_bound))
                            tm(do=partial(setattr, knockout_reaction, 'upper_bound', 0.),
                               undo=partial(setattr, knockout_reaction, 'upper_bound',
                                            knockout_reaction.upper_bound))
                    model.objective.direction = 'max'
                    try:
                        solution = model.solve()
                    except Unbounded:
                        fva_sol[reaction.id]['upper_bound'] = numpy.inf
                    except Infeasible:
                        fva_sol[reaction.id]['upper_bound'] = 0
                    else:
                        fva_sol[reaction.id]['upper_bound'] = solution.f

    df = pandas.DataFrame.from_dict(fva_sol, orient='index')
    lb_higher_ub = df[df.lower_bound > df.upper_bound]
    # Assert that these cases really only numerical artifacts
    assert ((lb_higher_ub.lower_bound - lb_higher_ub.upper_bound) < 1e-6).all()
    df.lower_bound[lb_higher_ub.index] = df.upper_bound[lb_higher_ub.index]

    return df


class _PhenotypicPhasePlaneChunkEvaluator(object):
    def __init__(self, model, variable_reactions):
        self.model = model
        self.variable_reactions = variable_reactions

    def __call__(self, points):
        return [self._production_envelope_inner(point) for point in points]

    def _production_envelope_inner(self, point):
        with TimeMachine() as tm:
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
            except (Infeasible, UndefinedSolution):  # Hack to handle GLPK bug
                solution = 0
            interval.append(solution)
            self.model.objective.direction = 'max'
            try:
                solution = self.model.solve().f
            except (Infeasible, UndefinedSolution):
                solution = 0
            interval.append(solution)
        return point + tuple(interval)


def flux_balance_impact_degree(model, knockouts, view=config.default_view, method="fva"):
    """
    Flux balance impact degree by Zhao et al 2013

    Parameters
    ----------
    model: SolverBasedModel
        Wild-type model
    knockouts: list
        Reactions to knockout
    method: str
        The method to compute the perturbation. default is "fva" - Flux Variability Analysis.
        It can also be computed with "em" - Elementary modes (not implemented)

    Returns
    -------
    FluxBalanceImpactDegreeResult: perturbation
        The changes in reachable reactions (reactions that can carry flux)
    """

    if method == "fva":
        reachable_reactions, perturbed_reactions = _fbid_fva(model, knockouts, view)
    elif method == "em":
        raise NotImplementedError("Elementary modes approach is not implemented")
    else:
        raise ValueError("%s method is not valid to compute Flux Balance Impact Degree" % method)

    return FluxBalanceImpactDegreeResult(reachable_reactions, perturbed_reactions, method)


def _fbid_fva(model, knockouts, view):
    with TimeMachine() as tm:

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

        wt_fva = flux_variability_analysis(model, view=view, remove_cycles=False)
        wt_fva._data_frame['upper_bound'] = wt_fva._data_frame.upper_bound.apply(numpy.round)
        wt_fva._data_frame['lower_bound'] = wt_fva._data_frame.lower_bound.apply(numpy.round)

        reachable_reactions = wt_fva.data_frame.query("lower_bound != 0 | upper_bound != 0")

        for reaction in model._ids_to_reactions(knockouts):
            reaction.knock_out(tm)

        mt_fva = flux_variability_analysis(model, reactions=reachable_reactions.index, view=view, remove_cycles=False)
        mt_fva._data_frame['upper_bound'] = mt_fva._data_frame.upper_bound.apply(numpy.round)
        mt_fva._data_frame['lower_bound'] = mt_fva._data_frame.lower_bound.apply(numpy.round)

        perturbed_reactions = []
        for reaction in reachable_reactions.index:
            if wt_fva.upper_bound(reaction) != mt_fva.upper_bound(reaction) or wt_fva.lower_bound(
                    reaction) != wt_fva.lower_bound(reaction):
                perturbed_reactions.append(reaction)

        return list(reachable_reactions.index), perturbed_reactions


class PhenotypicPhasePlaneResult(Result):
    def __init__(self, phase_plane, variable_ids, objective,
                 nice_variable_ids=None, nice_objective_id=None, *args, **kwargs):
        super(PhenotypicPhasePlaneResult, self).__init__(*args, **kwargs)
        self._phase_plane = phase_plane
        self.variable_ids = variable_ids
        self.nice_variable_ids = nice_variable_ids
        self.objective = objective
        self.nice_objective_id = nice_objective_id

    @property
    def data_frame(self):
        return pandas.DataFrame(self._phase_plane)

    def plot(self, grid=None, width=None, height=None, title=None, axis_font_size=None, palette=None,
             points=None, points_colors=None, **kwargs):
        if len(self.variable_ids) > 1:
            notice("Multi-dimensional plotting is not supported")
            return

        if title is None:
            title = "Phenotypic Phase Plane"
        variable = self.variable_ids[0]
        x_axis_label = self.nice_variable_ids[0]
        y_axis_label = self.nice_objective_id

        dataframe = pandas.DataFrame(columns=["ub", "lb", "value", "strain"])
        for _, row in self.iterrows():
            _df = pandas.DataFrame([[row['objective_upper_bound'], row['objective_lower_bound'], row[variable], "WT"]],
                                   columns=dataframe.columns)
            dataframe = dataframe.append(_df)

        plot = plotter.production_envelope(dataframe, grid=grid, width=width, height=height,
                                           title=title, y_axis_label=y_axis_label, x_axis_label=x_axis_label,
                                           palette=palette, points=points, points_colors=points_colors)
        if grid is None:
            plotter.display(plot)

    def __getitem__(self, item):
        return self._phase_plane[item]

    def iterrows(self):
        return self._phase_plane.iterrows()

    @property
    def area(self):
        area = 0
        for variable_id in self.variable_ids:
            area += self.area_for(variable_id)
        return area

    def area_for(self, variable_id):
        auc_max = trapz(self._phase_plane.objective_upper_bound.values, x=self._phase_plane[variable_id])
        auc_min = trapz(self._phase_plane.objective_lower_bound.values, x=self._phase_plane[variable_id])
        return auc_max - auc_min


class FluxVariabilityResult(Result):
    def __init__(self, data_frame, *args, **kwargs):
        super(FluxVariabilityResult, self).__init__(*args, **kwargs)
        self._data_frame = data_frame

    @property
    def data_frame(self):
        return self._data_frame

    def plot(self, index=None, grid=None, width=None, height=None, title=None, palette=None, **kwargs):
        if index is None:
            index = self.data_frame.index[0:10]
        fva_result = self.data_frame.loc[index]
        if title is None:
            title = "Flux Variability Analysis"

        dataframe = pandas.DataFrame(columns=["lb", "ub", "strain", "reaction"])
        for reaction_id, row in fva_result.iterrows():
            _df = pandas.DataFrame([[row['lower_bound'], row['upper_bound'], "WT", reaction_id]],
                                   columns=dataframe.columns)
            dataframe = dataframe.append(_df)

        plot = plotter.flux_variability_analysis(dataframe, grid=grid, width=width, height=height,
                                                 title=title, x_axis_label="Reactions", y_axis_label="Flux limits",
                                                 palette=palette)
        if grid is None:
            plotter.display(plot)

    def __getitem__(self, item):
        return self._data_frame[item]

    def upper_bound(self, item):
        if isinstance(item, Reaction):
            item = item.id
        return self['upper_bound'][item]

    def lower_bound(self, item):
        if isinstance(item, Reaction):
            item = item.id
        return self['lower_bound'][item]

    def iterrows(self):
        return self._data_frame.iterrows()


class FluxBalanceImpactDegreeResult(Result):
    def __init__(self, reachable_reactions, perturbed_reactions, method, *args, **kwargs):
        super(FluxBalanceImpactDegreeResult, self).__init__(*args, **kwargs)
        self._method = method
        self._reachable_reactions = reachable_reactions
        self._perturbed_reactions = perturbed_reactions

    def __contains__(self, item):
        if isinstance(item, Reaction):
            item = item.id
        return item in self._reachable_reactions

    def _repr_html_(self):
        return """
        <h3>Flux Balance Impact Degree</h3>
        <ul>
            <li> Degree: %i</li>
            <li> Reactions: %i</li>
        </ul>
        %s
        """ % (self.degree, len(self._reachable_reactions), self.data_frame._repr_html_())

    @property
    def degree(self):
        return len(self._perturbed_reactions)

    @property
    def data_frame(self):
        data_frame = pandas.DataFrame(columns=["perturbed"])
        for reaction in self._reachable_reactions:
            data_frame.loc[reaction] = [reaction in self._perturbed_reactions]
        return data_frame

    def plot(self, grid=None, width=None, height=None, title=None):
        pass
