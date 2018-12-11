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

from __future__ import absolute_import, print_function, division

import itertools
import logging
import re
from collections import OrderedDict
from functools import reduce
from operator import itemgetter

import numpy
import pandas
import six
from cobra import Reaction, Metabolite
from cobra.core import get_solution
from cobra.util import fix_objective_as_constraint, get_context
from cobra.exceptions import OptimizationError
from numpy import trapz
from six.moves import zip
from sympy import S
from optlang.interface import UNBOUNDED, OPTIMAL

import cameo
from cameo import config
from cameo.core.result import Result
from cameo.flux_analysis.util import remove_infeasible_cycles, fix_pfba_as_constraint
from cameo.parallel import SequentialView
from cameo.ui import notice
from cameo.util import partition, _BIOMASS_RE_
from cameo.core.utils import get_reaction_for
from cameo.visualization.plotting import plotter

logger = logging.getLogger(__name__)

__all__ = ['find_blocked_reactions', 'flux_variability_analysis', 'phenotypic_phase_plane',
           'flux_balance_impact_degree']


def knock_out_metabolite(metabolite, force_steady_state=False):
    """'Knockout' a metabolite. This can be done in 2 ways:

    1. Implementation follows the description in [1] "All fluxes around
    the metabolite M should be restricted to only produce the
    metabolite, for which balancing constraint of mass conservation is
    relaxed to allow nonzero values of the incoming fluxes whereas all
    outgoing fluxes are limited to zero."

    2. Force steady state All reactions consuming the metabolite are
    restricted to only produce the metabolite. A demand reaction is
    added to sink the metabolite produced to keep the problem feasible
    under the S.v = 0 constraint.


    Knocking out a metabolite overrules the constraints set on the
    reactions producing the metabolite.

    Parameters
    ----------
    force_steady_state: bool
        If True, uses approach 2.

    References
    ----------
    .. [1] Kim, P.-J., Lee, D.-Y., Kim, T. Y., Lee, K. H., Jeong, H.,
    Lee, S. Y., & Park, S. (2007). Metabolite essentiality elucidates
    robustness of Escherichia coli metabolism. PNAS, 104(34), 13638-13642

    """
    # restrict reactions to produce metabolite
    for rxn in metabolite.reactions:
        if rxn.metabolites[metabolite] > 0:
            rxn.bounds = (0, 0) if rxn.upper_bound < 0 \
                else (0, rxn.upper_bound)
        elif rxn.metabolites[metabolite] < 0:
            rxn.bounds = (0, 0) if rxn.lower_bound > 0 \
                else (rxn.lower_bound, 0)
    if force_steady_state:
        metabolite._model.add_boundary(metabolite, type="knock-out",
                                       lb=0, ub=1000,
                                       reaction_id="KO_{}".format(metabolite.id))
    else:
        previous_bounds = metabolite.constraint.lb, metabolite.constraint.ub
        metabolite.constraint.lb, metabolite.constraint.ub = -1000, 1000
        context = get_context(metabolite)
        if context:
            def reset():
                metabolite.constraint.lb, metabolite.constraint.ub = previous_bounds

            context(reset)


def find_essential_metabolites(model, threshold=1e-6, force_steady_state=False):
    """Return a list of essential metabolites.

    This can be done in 2 ways:

    1. Implementation follows the description in [1]:
        "All fluxes around the metabolite M should be restricted to only produce the metabolite,
         for which balancing constraint of mass conservation is relaxed to allow nonzero values
         of the incoming fluxes whereas all outgoing fluxes are limited to zero."

    2. Force Steady State approach:
        All reactions consuming the metabolite are restricted to only produce the metabolite. A demand
        reaction is added to sink the metabolite produced to keep the problem feasible under
        the S.v = 0 constraint.

    Briefly, for each metabolite, all reactions that consume that metabolite are blocked and if that makes the
    model either infeasible or results in near-zero flux in the model objective, then the metabolite is
    considered essential.

    Parameters
    ----------
    model : cobra.Model
        The model to find the essential metabolites for.
    threshold : float (default 1e-6)
        Minimal objective flux to be considered viable.
    force_steady_state: bool
        If True, uses approach 2.

    References
    ----------
    .. [1] Kim, P.-J., Lee, D.-Y., Kim, T. Y., Lee, K. H., Jeong, H., Lee, S. Y., & Park, S. (2007).
     Metabolite essentiality elucidates robustness of Escherichia coli metabolism. PNAS, 104(34), 13638–13642
    """

    essential = []
    # Essential metabolites are only in reactions that carry flux.
    metabolites = set()
    solution = model.optimize(raise_error=True)
    for reaction_id, flux in six.iteritems(solution.fluxes):
        if abs(flux) > 0:
            reaction = model.reactions.get_by_id(reaction_id)
            metabolites.update(reaction.metabolites.keys())

    for metabolite in metabolites:
        with model:
            knock_out_metabolite(metabolite, force_steady_state=force_steady_state)
            model.solver.optimize()
            if model.solver.status != OPTIMAL or model.objective.value < threshold:
                essential.append(metabolite)
    return essential


def find_blocked_reactions(model):
    """Determine reactions that cannot carry steady-state flux.

    Parameters
    ----------
    model: cobra.Model

    Returns
    -------
    list
        A list of reactions.

    """
    with model:
        for exchange in model.exchanges:
            exchange.bounds = (-9999, 9999)
        fva_solution = flux_variability_analysis(model)
    return frozenset(
        reaction for reaction in model.reactions
        if round(
            fva_solution.lower_bound(reaction.id),
            config.ndecimals) == 0 and round(
            fva_solution.upper_bound(reaction.id), config.ndecimals) == 0
    )


def flux_variability_analysis(model, reactions=None, fraction_of_optimum=0., pfba_factor=None,
                              remove_cycles=False, view=None):
    """Flux variability analysis.

    Parameters
    ----------
    model : cobra.Model
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
    with model:
        if fraction_of_optimum > 0.:
            fix_objective_as_constraint(model, fraction=fraction_of_optimum)
        if pfba_factor is not None:
            # don't add the objective-constraint again so fraction_of_optimum=0
            fix_pfba_as_constraint(model, multiplier=pfba_factor, fraction_of_optimum=0)
        reaction_chunks = (chunk for chunk in partition(reactions, len(view)))
        if remove_cycles:
            func_obj = _FvaFunctionObject(model, _cycle_free_fva)
        else:
            func_obj = _FvaFunctionObject(model, _flux_variability_analysis)
        chunky_results = view.map(func_obj, reaction_chunks)
        solution = pandas.concat(chunky_results)

    return FluxVariabilityResult(solution)


def phenotypic_phase_plane(model, variables, objective=None, source=None, points=20, view=None):
    """Phenotypic phase plane analysis [1].

    Implements a phenotypic phase plan analysis with interpretation same as
    that presented in [1] but calculated by optimizing the model for all
    steps of the indicated variables (instead of using shadow prices).

    Parameters
    ----------
    model: cobra.Model
    variables: str or reaction or iterable
        A reaction ID, reaction, or list of reactions to be varied.
    objective: str or reaction or optlang.Objective or Metabolite, optional
        An objective, a reaction's flux, or a metabolite's production to be minimized/maximized
        (defaults to the current model objective).
    source: Reaction or reaction identifier
        The reaction to use as the source when calculating mass and carbon yield. Set to the medium reaction with the
        highest input carbon flux if left None.
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

    if isinstance(variables, six.string_types):
        variables = [variables]
    elif isinstance(variables, Reaction):
        variables = [variables]
    variable_ids = [var if isinstance(var, six.string_types) else var.id for var in variables]

    if view is None:
        view = config.default_view
    with model:
        if objective is not None:
            if isinstance(objective, Metabolite):
                try:
                    objective = model.reactions.get_by_id("DM_%s" % objective.id)
                except KeyError:
                    objective = model.add_boundary(objective, type='demand')
            # try:
            #     objective = model.reaction_for(objective, time_machine=tm)
            # except KeyError:
            #     pass

            model.objective = objective

        if source is not None:
            source_reaction = get_reaction_for(model, source)
        else:
            source_reaction = get_c_source_reaction(model)

        variable_reactions = model.reactions.get_by_any(variables)
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
    phase_plane = pandas.DataFrame(
        envelope, columns=(variable_reactions_ids + [
            'objective_lower_bound',
            'objective_upper_bound',
            'c_yield_lower_bound',
            'c_yield_upper_bound',
            'mass_yield_lower_bound',
            'mass_yield_upper_bound'
        ])
    )

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
        elif len(reaction.metabolites) < 5:
            nice_id = reaction
        else:
            nice_id = reaction.id
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
        reactions = model.reactions.get_by_any(reactions)
    fva_sol = OrderedDict()
    lb_flags = dict()
    with model:
        model.objective = S.Zero

        model.objective.direction = 'min'
        for reaction in reactions:
            lb_flags[reaction.id] = False
            fva_sol[reaction.id] = dict()
            model.solver.objective.set_linear_coefficients({reaction.forward_variable: 1.,
                                                            reaction.reverse_variable: -1.})
            model.solver.optimize()
            if model.solver.status == OPTIMAL:
                fva_sol[reaction.id]['lower_bound'] = model.objective.value
            elif model.solver.status == UNBOUNDED:
                fva_sol[reaction.id]['lower_bound'] = -numpy.inf
            else:
                lb_flags[reaction.id] = True
            model.solver.objective.set_linear_coefficients({reaction.forward_variable: 0.,
                                                            reaction.reverse_variable: 0.})

            assert model.objective.expression == 0, model.objective.expression

        model.objective.direction = 'max'
        for reaction in reactions:
            ub_flag = False
            model.solver.objective.set_linear_coefficients({reaction.forward_variable: 1.,
                                                            reaction.reverse_variable: -1.})

            model.solver.optimize()
            if model.solver.status == OPTIMAL:
                fva_sol[reaction.id]['upper_bound'] = model.objective.value
            elif model.solver.status == UNBOUNDED:
                fva_sol[reaction.id]['upper_bound'] = numpy.inf
            else:
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


def get_c_source_reaction(model):
    """ carbon source reactions

    Returns
    -------
    Reaction
       The medium reaction with highest input carbon flux
    """
    try:
        model.slim_optimize(error_value=None)
    except OptimizationError:
        return None
    medium_reactions = model.reactions.get_by_any(list(model.medium))
    medium_reactions_fluxes = [(rxn, total_carbon_flux(rxn, consumption=True)) for rxn in medium_reactions]
    source_reactions = [(rxn, c_flux) for rxn, c_flux in medium_reactions_fluxes if c_flux > 0]
    try:
        return max(source_reactions, key=itemgetter(1))[0]
    except ValueError:
        return None


def total_carbon_flux(reaction, consumption=True):
    """summed product carbon flux for a reaction

    Parameters
    ----------
    reaction : Reaction
        the reaction to carbon return flux for
    consumption : bool
        flux for consumed metabolite, else produced

    Returns
    -------
    float
        reaction flux multiplied by number of carbon for the products of the reaction"""
    direction = 1 if consumption else -1
    c_flux = [reaction.flux * coeff * met.elements.get('C', 0) * direction
              for met, coeff in reaction.metabolites.items()]
    return sum([flux for flux in c_flux if flux > 0])


def single_flux(reaction, consumption=True):
    """flux into single product for a reaction

    only defined for reactions with single products

    Parameters
    ----------
    reaction : Reaction
        the reaction to product flux for
    consumption : bool
        flux for consumed metabolite, else produced

    Returns
    -------
    tuple
        metabolite, flux for the metabolite"""
    if len(list(reaction.metabolites)) != 1:
        raise ValueError('product flux only defined for single metabolite reactions')
    met, coeff = next(six.iteritems(reaction.metabolites))
    direction = 1 if consumption else -1
    return met, reaction.flux * coeff * direction


def _cycle_free_fva(model, reactions=None, sloppy=True, sloppy_bound=666):
    """Cycle free flux-variability analysis. (http://cran.r-project.org/web/packages/sybilcycleFreeFlux/index.html)

    Parameters
    ----------
    model : cobra.Model
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
        reactions = model.reactions.get_by_any(reactions)
    fva_sol = OrderedDict()
    for reaction in reactions:
        fva_sol[reaction.id] = dict()
        model.objective = reaction
        model.objective.direction = 'min'
        model.solver.optimize()
        if model.solver.status == UNBOUNDED:
            fva_sol[reaction.id]['lower_bound'] = -numpy.inf
            continue
        elif model.solver.status != OPTIMAL:
            fva_sol[reaction.id]['lower_bound'] = 0
            continue
        bound = model.objective.value
        if sloppy and abs(bound) < sloppy_bound:
            fva_sol[reaction.id]['lower_bound'] = bound
        else:
            logger.debug('Determine if {} with bound {} is a cycle'.format(reaction.id, bound))
            solution = get_solution(model)
            v0_fluxes = solution.fluxes
            v1_cycle_free_fluxes = remove_infeasible_cycles(model, v0_fluxes)
            if abs(v1_cycle_free_fluxes[reaction.id] - bound) < 10 ** -6:
                fva_sol[reaction.id]['lower_bound'] = bound
            else:
                logger.debug('Cycle detected: {}'.format(reaction.id))
                cycle_count += 1
                v2_one_cycle_fluxes = remove_infeasible_cycles(model, v0_fluxes, fix=[reaction.id])
                with model:
                    for key, v1_flux in six.iteritems(v1_cycle_free_fluxes):
                        if round(v1_flux, config.ndecimals) == 0 and round(v2_one_cycle_fluxes[key],
                                                                           config.ndecimals) != 0:
                            knockout_reaction = model.reactions.get_by_id(key)
                            knockout_reaction.knock_out()
                    model.objective.direction = 'min'
                    model.solver.optimize()
                    if model.solver.status == OPTIMAL:
                        fva_sol[reaction.id]['lower_bound'] = model.objective.value
                    elif model.solver.status == UNBOUNDED:
                        fva_sol[reaction.id]['lower_bound'] = -numpy.inf
                    else:
                        fva_sol[reaction.id]['lower_bound'] = 0

    for reaction in reactions:
        model.objective = reaction
        model.objective.direction = 'max'
        model.solver.optimize()
        if model.solver.status == UNBOUNDED:
            fva_sol[reaction.id]['upper_bound'] = numpy.inf
            continue
        elif model.solver.status != OPTIMAL:
            fva_sol[reaction.id]['upper_bound'] = 0
            continue
        bound = model.objective.value
        if sloppy and abs(bound) < sloppy_bound:
            fva_sol[reaction.id]['upper_bound'] = bound
        else:
            logger.debug('Determine if {} with bound {} is a cycle'.format(reaction.id, bound))
            solution = get_solution(model)
            v0_fluxes = solution.fluxes
            v1_cycle_free_fluxes = remove_infeasible_cycles(model, v0_fluxes)
            if abs(v1_cycle_free_fluxes[reaction.id] - bound) < 1e-6:
                fva_sol[reaction.id]['upper_bound'] = v0_fluxes[reaction.id]
            else:
                logger.debug('Cycle detected: {}'.format(reaction.id))
                cycle_count += 1
                v2_one_cycle_fluxes = remove_infeasible_cycles(model, v0_fluxes, fix=[reaction.id])
                with model:
                    for key, v1_flux in six.iteritems(v1_cycle_free_fluxes):
                        if round(v1_flux, config.ndecimals) == 0 and round(v2_one_cycle_fluxes[key],
                                                                           config.ndecimals) != 0:
                            knockout_reaction = model.reactions.get_by_id(key)
                            knockout_reaction.knock_out()
                    model.objective.direction = 'max'
                    model.solver.optimize()
                    if model.solver.status == OPTIMAL:
                        fva_sol[reaction.id]['upper_bound'] = model.objective.value
                    elif model.solver.status == UNBOUNDED:
                        fva_sol[reaction.id]['upper_bound'] = numpy.inf
                    else:
                        fva_sol[reaction.id]['upper_bound'] = 0

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
            raise NotImplementedError('complex objectives not supported')
        self.product_reaction = objective_reactions[0]

    def carbon_yield(self):
        """ mol product per mol carbon input

        Returns
        -------
        float
            the mol carbon atoms in the product (as defined by the model objective) divided by the mol carbon in the
            input reactions (as defined by the model medium) or zero in case of division by zero arises"""
        if self.source is None:
            return numpy.nan
        carbon_input_flux = total_carbon_flux(self.source, consumption=True)
        carbon_output_flux = total_carbon_flux(self.product_reaction, consumption=False)
        try:
            return carbon_output_flux / carbon_input_flux
        except ZeroDivisionError:
            return numpy.nan

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
        if self.source is None:
            return numpy.nan
        try:
            source, source_flux = single_flux(self.source, consumption=True)
            product, product_flux = single_flux(self.product_reaction, consumption=False)
        except ValueError:
            return numpy.nan
        mol_prod_mol_src = product_flux / source_flux
        return (mol_prod_mol_src * product.formula_weight) / source.formula_weight

    def __call__(self, points):
        return [self._production_envelope_inner(point) for point in points]

    def _interval_estimates(self):
        self.model.solver.optimize()
        if self.model.solver.status == OPTIMAL:
            flux = self.model.solver.objective.value
            carbon_yield = self.carbon_yield()
            mass_yield = self.mass_yield()
        else:
            flux = numpy.nan
            carbon_yield = numpy.nan
            mass_yield = numpy.nan
        return flux, carbon_yield, mass_yield

    def _production_envelope_inner(self, point):
        original_direction = self.model.objective.direction
        with self.model:
            for (reaction, coordinate) in zip(self.variable_reactions, point):
                reaction.bounds = (coordinate, coordinate)
            interval = []
            interval_carbon_yield = []
            interval_mass_yield = []

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
        self.model.objective.direction = original_direction

        intervals = tuple(interval) + tuple(interval_carbon_yield) + tuple(interval_mass_yield)
        return point + intervals


def flux_balance_impact_degree(model, knockouts, view=config.default_view, method="fva"):
    """
    Flux balance impact degree by Zhao et al 2013

    Parameters
    ----------
    model: cobra.Model
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
    with model:

        for reaction in model.reactions:
            if reaction.reversibility:
                reaction.bounds = (-1, 1)
            else:
                reaction.bounds = (0, 1)

        wt_fva = flux_variability_analysis(model, view=view, remove_cycles=False)
        wt_fva._data_frame['upper_bound'] = wt_fva._data_frame.upper_bound.apply(numpy.round)
        wt_fva._data_frame['lower_bound'] = wt_fva._data_frame.lower_bound.apply(numpy.round)

        reachable_reactions = wt_fva.data_frame.query("lower_bound != 0 | upper_bound != 0")

        for reaction in model.reactions.get_by_any(knockouts):
            reaction.knock_out()

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
            raise ValueError('estimate must be one of %s' % ', '.join(possible_estimates.keys()))
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


def n_carbon(reaction):
    return sum(metabolite.elements.get('C', 0) for metabolite in reaction.metabolites)
