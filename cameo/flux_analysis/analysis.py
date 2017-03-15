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

import itertools
import logging
import re
from collections import OrderedDict
from functools import partial, reduce

import numpy
import pandas
import six
from cobra.core import Reaction, Metabolite
from numpy import trapz
from six.moves import zip
from sympy import S

import cameo
from cameo import config
from cameo.core.result import Result
from cameo.exceptions import Infeasible, Unbounded
from cameo.flux_analysis.util import remove_infeasible_cycles, fix_pfba_as_constraint
from cameo.parallel import SequentialView
from cameo.ui import notice
from cameo.util import TimeMachine, partition, _BIOMASS_RE_
from cameo.visualization.plotting import plotter

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
            exchange.change_bounds(-9999, 9999, tm)
        fva_solution = flux_variability_analysis(model)
    return frozenset(reaction for reaction in model.reactions
                     if round(fva_solution.lower_bound(reaction.id), config.ndecimals) == 0 and
                     round(fva_solution.upper_bound(reaction.id), config.ndecimals) == 0)


def flux_variability_analysis(model, reactions=None, fraction_of_optimum=0., pfba_factor=None,
                              remove_cycles=False, view=None):
    """Flux variability analysis.

    Parameters
    ----------
    model : cameo.core.SolverBasedModel
    reactions: None or iterable
        The list of reaction whose lower and upper bounds should be determined.
        If `None`, all reactions in `model` will be assessed.
    fraction_of_optimum : float
        Fix the objective of the model to a fraction of it's max. Expected to be within [0;1]. Lower values increase
        variability.
    pfba_factor : float
        If not None, fix the total sum of reaction fluxes to its minimum times a factor. Expected to be within [
        1;inf]. Higher factors increase flux variability to a certain point since the bound for the objective is
        still fixed.
    remove_cycles : bool
        If true, apply the CycleFreeFlux algorithm to remove loops from each simulated flux distribution.
    view: cameo.parallel.SequentialView or cameo.parallel.MultiprocessingView or ipython.cluster.DirectView
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
        if pfba_factor is not None:
            # don't add the objective-constraint again so fraction_of_optimum=0
            fix_pfba_as_constraint(model, multiplier=pfba_factor, time_machine=tm, fraction_of_optimum=0)
        tm(do=int, undo=partial(setattr, model, "objective", model.objective))
        reaction_chunks = (chunk for chunk in partition(reactions, len(view)))
        if remove_cycles:
            func_obj = _FvaFunctionObject(model, _cycle_free_fva)
        else:
            func_obj = _FvaFunctionObject(model, _flux_variability_analysis)
        chunky_results = view.map(func_obj, reaction_chunks)
        solution = pandas.concat(chunky_results)

    return FluxVariabilityResult(solution)


def phenotypic_phase_plane(model, variables=[], objective=None, source=None, points=20, view=None):
    """Phenotypic phase plane analysis [1].

    Implements a phenotypic phase plan analysis with interpretation same as
    that presented in [1] but calculated by optimizing the model for all
    steps of the indicated variables (instead of using shadow prices).

    Parameters
    ----------
    model: SolverBasedModel
    variables: str or reaction or iterable
        A reaction ID, reaction, or list of reactions to be varied.
    objective: str or reaction or optlang.Objective or Metabolite, optional
        An objective, a reaction's flux, or a metabolite's production to be minimized/maximized
        (defaults to the current model objective).
    source: Reaction or reaction identifier
        The reaction to use as the source when calculating mass and carbon yield. Set to the medium reaction with the
        highest input carbon flux if left none.
    points: int or iterable
        Number of points to be interspersed between the variable bounds.
        A list of same same dimensions as `variables` can be used to specify
        variable specific numbers of points.
    view: SequentialView or MultiprocessingView or ipython.cluster.DirectView
        A parallelization view.

    Returns
    -------
    PhenotypicPhasePlaneResult
        The phenotypic phase plane with flux, mol carbon input / mol carbon
        output, gram input / gram output

    References
    ----------
    [1] Edwards, J. S., Ramakrishna, R. and Palsson, B. O. (2002). Characterizing the metabolic phenotype: a phenotype
        phase plane analysis. Biotechnology and Bioengineering, 77(1), 27–36. doi:10.1002/bit.10047
    """

    if isinstance(variables, str):
        variables = [variables]
    elif isinstance(variables, cameo.core.reaction.Reaction):
        variables = [variables]
    variable_ids = [var if isinstance(var, str) else var.id for var in variables]

    if view is None:
        view = config.default_view
    with TimeMachine() as tm:
        if objective is not None:
            if isinstance(objective, Metabolite):
                try:
                    objective = model.reactions.get_by_id("DM_%s" % objective.id)
                except KeyError:
                    objective = model.add_demand(objective, time_machine=tm)
            # try:
            #     objective = model.reaction_for(objective, time_machine=tm)
            # except KeyError:
            #     pass

            model.change_objective(objective, time_machine=tm)

        if source:
            source_reaction = model._reaction_for(source)
        else:
            source_reaction = _get_c_source_reaction(model)

        variable_reactions = model._ids_to_reactions(variables)
        variables_min_max = flux_variability_analysis(model, reactions=variable_reactions, view=SequentialView())
        grid = [numpy.linspace(lower_bound, upper_bound, points, endpoint=True) for
                reaction_id, lower_bound, upper_bound in
                variables_min_max.data_frame.loc[variable_ids].itertuples()]
        grid_generator = itertools.product(*grid)
        chunks_of_points = partition(list(grid_generator), len(view))
        evaluator = _PhenotypicPhasePlaneChunkEvaluator(model, variable_reactions, source_reaction)
        chunk_results = view.map(evaluator, chunks_of_points)
        envelope = reduce(list.__add__, chunk_results)

    nice_variable_ids = [_nice_id(reaction) for reaction in variable_reactions]
    variable_reactions_ids = [reaction.id for reaction in variable_reactions]
    phase_plane = pandas.DataFrame(envelope,
                                   columns=(variable_reactions_ids +
                                            ['objective_lower_bound',
                                             'objective_upper_bound',
                                             'c_yield_lower_bound',
                                             'c_yield_upper_bound',
                                             'mass_yield_lower_bound',
                                             'mass_yield_upper_bound']))

    if objective is None:
        objective = model.objective
    nice_objective_id = _nice_id(objective)

    objective = objective.id if isinstance(objective, Reaction) else str(objective)
    return PhenotypicPhasePlaneResult(phase_plane, variable_reactions_ids, objective,
                                      nice_variable_ids=nice_variable_ids,
                                      source_reaction=_nice_id(source_reaction),
                                      nice_objective_id=nice_objective_id)


def _nice_id(reaction):
    if isinstance(reaction, Reaction):
        if hasattr(reaction, 'nice_id'):
            nice_id = reaction.nice_id
        else:
            nice_id = reaction
    else:
        nice_id = str(reaction)
    return nice_id


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
    lb_flags = dict()
    with TimeMachine() as tm:
        model.change_objective(S.Zero, time_machine=tm)

        model.objective.direction = 'min'
        for reaction in reactions:
            lb_flags[reaction.id] = False
            fva_sol[reaction.id] = dict()
            model.solver.objective.set_linear_coefficients({reaction.forward_variable: 1.,
                                                            reaction.reverse_variable: -1.})
            try:
                solution = model.solve()
                fva_sol[reaction.id]['lower_bound'] = solution.f
            except Unbounded:
                fva_sol[reaction.id]['lower_bound'] = -numpy.inf
            except Infeasible:
                lb_flags[reaction.id] = True
            model.solver.objective.set_linear_coefficients({reaction.forward_variable: 0.,
                                                            reaction.reverse_variable: 0.})

            assert model.objective.expression == 0, model.objective.expression

        model.objective.direction = 'max'
        for reaction in reactions:
            ub_flag = False
            model.solver.objective.set_linear_coefficients({reaction.forward_variable: 1.,
                                                            reaction.reverse_variable: -1.})

            try:
                solution = model.solve()
                fva_sol[reaction.id]['upper_bound'] = solution.f
            except Unbounded:
                fva_sol[reaction.id]['upper_bound'] = numpy.inf
            except Infeasible:
                ub_flag = True

            if lb_flags[reaction.id] is True and ub_flag is True:
                fva_sol[reaction.id]['lower_bound'] = 0
                fva_sol[reaction.id]['upper_bound'] = 0
                [lb_flags[reaction.id], ub_flag] = [False, False]
            elif lb_flags[reaction.id] is True and ub_flag is False:
                fva_sol[reaction.id]['lower_bound'] = fva_sol[reaction.id]['upper_bound']
                lb_flags[reaction.id] = False
            elif lb_flags[reaction.id] is False and ub_flag is True:
                fva_sol[reaction.id]['upper_bound'] = fva_sol[reaction.id]['lower_bound']
                ub_flag = False

            model.solver.objective.set_linear_coefficients({reaction.forward_variable: 0.,
                                                            reaction.reverse_variable: 0.})

            assert model.objective.expression == 0, model.objective.expression

            assert lb_flags[reaction.id] is False and ub_flag is False, "Something is wrong with FVA (%s)" % reaction.id

    df = pandas.DataFrame.from_dict(fva_sol, orient='index')
    lb_higher_ub = df[df.lower_bound > df.upper_bound]
    # this is an alternative solution to what I did above with flags
    # Assert that these cases really only numerical artifacts
    try:
        assert ((lb_higher_ub.lower_bound - lb_higher_ub.upper_bound) < 1e-6).all()
    except AssertionError:
        logger.debug(list(zip(model.reactions, (lb_higher_ub.lower_bound - lb_higher_ub.upper_bound) < 1e-6)))
    df.lower_bound[lb_higher_ub.index] = df.upper_bound[lb_higher_ub.index]
    return df


def _get_c_source_reaction(model):
    """ carbon source reactions

    Returns
    -------
    Reaction
       The medium reaction with highest input carbon flux
    """
    medium_reactions = [model.reactions.get_by_id(reaction) for reaction in model.medium.reaction_id]
    try:
        model.solve()
    except (Infeasible, AssertionError):
        return None
    source_reactions = [(reaction, reaction.flux * reaction.n_carbon) for reaction in medium_reactions if
                        reaction.flux < 0]
    sorted_sources = sorted(source_reactions, key=lambda reaction_tuple: reaction_tuple[1])
    return sorted_sources[0][0]


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
                            knockout_reaction.knock_out(time_machine=tm)
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
                            knockout_reaction.knock_out(time_machine=tm)
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
    def __init__(self, model, variable_reactions, source):
        self.model = model
        self.source = source
        self.variable_reactions = variable_reactions
        objective_reactions = [reaction for reaction in model.reactions
                               if reaction.objective_coefficient != 0]
        if len(objective_reactions) != 1:
            raise Exception('multiple objectives not supported')
        self.product_reaction = objective_reactions[0]

    @classmethod
    def carbon_flux(cls, reaction):
        """ carbon flux for reactions

        Parameters
        ----------
        reaction : Reaction
            the reaction to carbon return flux for

        Returns
        -------
        float
            reaction flux multiplied by number of carbon in reactants"""
        carbon = sum(metabolite.n_carbon for metabolite in reaction.reactants)
        try:
            return reaction.flux * carbon
        except AssertionError:
            # optimized flux can be out of bounds, in that case return as missing value
            return numpy.nan

    def carbon_yield(self):
        """ mol product per mol carbon input

        Returns
        -------
        float
            the mol carbon atoms in the product (as defined by the model objective) divided by the mol carbon in the
            input reactions (as defined by the model medium) or zero in case of division by zero arises"""
        carbon_input_flux = self.carbon_flux(self.source)
        carbon_output_flux = self.carbon_flux(self.product_reaction)
        try:
            return carbon_output_flux / (carbon_input_flux * -1)
        except ZeroDivisionError:
            return 0

    def mass_yield(self):
        """ gram product divided by gram of feeding source

        only defined when we have only one product (as defined by the model
        objective) and only one compound as carbon source (as defined by the
        model medium).

        Returns
        -------
        float
            gram product per 1 g of feeding source or nan if more than one product or feeding source
        """
        too_long = (len(self.source.metabolites) > 1,
                    len(self.product_reaction.metabolites) > 1)
        if any(too_long):
            return numpy.nan
        try:
            source_flux = self.source.flux
            product_flux = self.product_reaction.flux
        except AssertionError:
            return numpy.nan
        source_metabolite = list(self.source.metabolites)[0]
        product_metabolite = list(self.product_reaction.metabolites)[0]
        product_mass = product_metabolite.formula_weight
        source_mass = source_metabolite.formula_weight
        try:
            mol_prod_mol_src = product_flux / (source_flux * -1)
            return (mol_prod_mol_src * product_mass) / source_mass
        except ZeroDivisionError:
            return 0

    def __call__(self, points):
        return [self._production_envelope_inner(point) for point in points]

    def _interval_estimates(self):
        try:
            flux = self.model.solve().f
            carbon_yield = self.carbon_yield()
            mass_yield = self.mass_yield()
        except Infeasible:
            flux = 0
            carbon_yield = 0
            mass_yield = 0
        return flux, carbon_yield, mass_yield

    def _production_envelope_inner(self, point):
        with TimeMachine() as tm:
            for (reaction, coordinate) in zip(self.variable_reactions, point):
                reaction.change_bounds(coordinate, coordinate, time_machine=tm)
            interval = []
            interval_carbon_yield = []
            interval_mass_yield = []
            tm(do=int, undo=partial(setattr, self.model.objective, 'direction', self.model.objective.direction))

            self.model.objective.direction = 'min'
            flux, carbon_yield, mass_yield = self._interval_estimates()
            interval.append(flux)
            interval_carbon_yield.append(carbon_yield)
            interval_mass_yield.append(mass_yield)

            self.model.objective.direction = 'max'
            flux, carbon_yield, mass_yield = self._interval_estimates()
            interval.append(flux)
            interval_carbon_yield.append(carbon_yield)
            interval_mass_yield.append(mass_yield)

        intervals = tuple(interval) + tuple(interval_carbon_yield) + tuple(interval_mass_yield)
        return point + intervals


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
                reaction.change_bounds(-1, 1, time_machine=tm)
            else:
                reaction.change_bounds(0, 1, time_machine=tm)

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
                 nice_variable_ids=None, nice_objective_id=None,
                 source_reaction=None, *args, **kwargs):
        super(PhenotypicPhasePlaneResult, self).__init__(*args, **kwargs)
        self._phase_plane = phase_plane
        self.variable_ids = variable_ids
        self.nice_variable_ids = nice_variable_ids
        self.objective = objective
        self.nice_objective_id = nice_objective_id
        self.source_reaction = source_reaction

    @property
    def data_frame(self):
        return pandas.DataFrame(self._phase_plane)

    def plot(self, grid=None, width=None, height=None, title=None, axis_font_size=None, palette=None,
             points=None, points_colors=None, estimate='flux', **kwargs):
        """plot phenotypic phase plane result

        create a plot of a phenotypic phase plane analysis

        Parameters
        ----------
        grid: plotting grid
            the grid for plotting
        width: int
            the width of the plot
        height: int
            the height of the plot
        title: string
            the height of the plot
        axis_font_size: int
            the font sizes for the axis
        palette: string
            name of color palette to use, e.g. RdYlBlu
        points: iterable of points
            additional points to plot as x, y iterable
        points_colors: iterable of strings
            iterable with colors for the points
        estimate: string
            either flux, mass_yield (g output / g output) or c_yield (mol carbon output / mol carbon input)
        """
        possible_estimates = {'flux': ('objective_upper_bound',
                                       'objective_lower_bound',
                                       'flux',
                                       '[mmol gDW^-1 h^-1]'),
                              'mass_yield': ('mass_yield_upper_bound',
                                             'mass_yield_lower_bound',
                                             'mass yield, src={}'.format(self.source_reaction),
                                             '[g/g(src) h^-1]'),
                              'c_yield': ('c_yield_upper_bound',
                                          'c_yield_lower_bound',
                                          'carbon yield, src={}'.format(self.source_reaction),
                                          '[mmol(C)/mmol(C(src)) h^-1]')}
        if estimate not in possible_estimates:
            raise Exception('estimate must be one of %s' %
                            ', '.join(possible_estimates.keys()))
        upper, lower, description, unit = possible_estimates[estimate]
        if title is None:
            title = "Phenotypic Phase Plane ({})".format(description)
        if len(self.variable_ids) == 1:

            variable = self.variable_ids[0]
            y_axis_label = self._axis_label(self.objective, self.nice_objective_id, unit)
            x_axis_label = self._axis_label(variable, self.nice_variable_ids[0], '[mmol gDW^-1 h^-1]')

            dataframe = pandas.DataFrame(columns=["ub", "lb", "value", "strain"])
            for _, row in self.iterrows():
                _df = pandas.DataFrame([[row[upper], row[lower], row[variable], "WT"]],
                                       columns=dataframe.columns)
                dataframe = dataframe.append(_df)

            plot = plotter.production_envelope(dataframe, grid=grid, width=width, height=height,
                                               title=title, y_axis_label=y_axis_label, x_axis_label=x_axis_label,
                                               palette=palette, points=points, points_colors=points_colors)

        elif len(self.variable_ids) == 2:
            var_1 = self.variable_ids[0]
            var_2 = self.variable_ids[1]
            x_axis_label = self._axis_label(var_1, self.nice_variable_ids[0], '[mmol gDW^-1 h^-1]')
            y_axis_label = self._axis_label(var_2, self.nice_variable_ids[1], '[mmol gDW^-1 h^-1]')
            z_axis_label = self._axis_label(self.objective, self.nice_objective_id, unit)

            dataframe = pandas.DataFrame(columns=["ub", "lb", "value1", "value2", "strain"])
            for _, row in self.iterrows():
                _df = pandas.DataFrame([[row[upper], row[lower],
                                         row[var_1], row[var_2], "WT"]],
                                       columns=dataframe.columns)
                dataframe = dataframe.append(_df)

            plot = plotter.production_envelope_3d(dataframe, grid=grid, width=width, height=height,
                                                  title=title, y_axis_label=y_axis_label, x_axis_label=x_axis_label,
                                                  z_axis_label=z_axis_label, palette=palette, points=points,
                                                  points_colors=points_colors)

        else:
            notice("Multi-dimensional plotting is not supported")
            return

        if grid is None:
            plotter.display(plot)

    @staticmethod
    def _axis_label(variable_id, nice_variable_id, unit):
        if re.search(_BIOMASS_RE_, variable_id):
            return '{} [h^-1]'.format(nice_variable_id)
        else:
            return '{} {}'.format(nice_variable_id, unit)

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
        data_frame = self._phase_plane.sort_values(by=variable_id, ascending=True)
        auc_max = trapz(data_frame.objective_upper_bound.values, x=data_frame[variable_id])
        auc_min = trapz(data_frame.objective_lower_bound.values, x=data_frame[variable_id])
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
                                                 title=title, y_axis_label="Reactions", x_axis_label="Flux limits",
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
