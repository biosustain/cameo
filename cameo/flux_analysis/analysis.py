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

import itertools
from copy import copy
from collections import OrderedDict
from functools import partial
import numpy
from cobra.core import Reaction
from cameo import config
from cameo.exceptions import UndefinedSolution, Infeasible, Unbounded, SolveError
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
        if model.reversible_encoding == 'split':
            tm(do=partial(setattr, model, 'reversible_encoding', 'unsplit'),
               undo=partial(setattr, model, 'reversible_encoding', 'split'))
        if fraction_of_optimum > 0.:
            try:
                obj_val = model.solve()
            except SolveError as e:
                print "flux_variability_analyis was not able to determine an optimal solution for objective %s" % model.objective
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
    if not hasattr(variables, '__iter__'):
        variables = [variables]
    if view is None:
        view = config.default_view
    tm = TimeMachine()
    if model.reversible_encoding == 'split':
        tm(do=partial(setattr, model, 'reversible_encoding', 'unsplit'),
           undo=partial(setattr, model, 'reversible_encoding', 'split'))
    if objective is not None:
        tm(do=partial(setattr, model, 'objective', objective),
           undo=partial(setattr, model, 'objective', model.objective))

    variable_reactions = _ids_to_reactions(model, variables)
    variables_min_max = flux_variability_analysis(model, reactions=variable_reactions, view=view)
    grid = [numpy.linspace(lower_bound, upper_bound, points, endpoint=True) for reaction_id, lower_bound, upper_bound in
            variables_min_max.itertuples()]
    grid_generator = itertools.product(*grid)
    original_bounds = dict([(reaction, (reaction.lower_bound, reaction.upper_bound))
                            for reaction in variable_reactions])

    chunks_of_points = partition(list(grid_generator), len(view))
    evaluator = _PhenotypicPhasePlaneChunkEvaluator(model, variable_reactions)
    chunk_results = view.map(evaluator, chunks_of_points)
    envelope = reduce(list.__add__, chunk_results)

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
        try:
            solution = model.solve()
            fva_sol[reaction.id]['lower_bound'] = solution.f
        except Unbounded:
            fva_sol[reaction.id]['lower_bound'] = -numpy.inf
        except Infeasible:
            fva_sol[reaction.id]['lower_bound'] = 0
    for reaction in reactions:
        model.objective = reaction
        model.objective.direction = 'max'
        try:
            solution = model.solve()
            fva_sol[reaction.id]['upper_bound'] = solution.f
        except Unbounded:
            fva_sol[reaction.id]['upper_bound'] = numpy.inf
        except Infeasible:
            fva_sol[reaction.id]['lower_bound'] = 0
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
        try:
            model.objective = reaction
        except:
            pass
        model.objective.direction = 'min'
        try:
            solution = model.solve()
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
                    try:
                        solution = model_broken_cyle.optimize()
                    except Unbounded:
                        fva_sol[reaction.id]['lower_bound'] = -numpy.inf
        except Unbounded:
            fva_sol[reaction.id]['lower_bound'] = -numpy.inf
        except Infeasible:
            fva_sol[reaction.id]['lower_bound'] = 0
        except Exception as e:
            print reaction.id
            raise e
    for reaction in reactions:
        model.objective = reaction
        model.objective.direction = 'max'
        try:
            solution = model.solve()
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
                    try:
                        solution = model_broken_cyle.solve()
                        fva_sol[reaction.id]['upper_bound'] = solution.f
                    except Unbounded:
                        fva_sol[reaction.id]['upper_bound'] = numpy.inf
        except Unbounded:
            fva_sol[reaction.id]['upper_bound'] = numpy.inf
        except Infeasible:
            fva_sol[reaction.id]['upper_bound'] = 0
        except Exception as e:
            raise e
            # print model.reactions.Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2.lower_bound
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


def fbip(model, knockouts, view=config.default_view, method="fva"):
    """
    Flux balance impact degree by Zhao et al 2013

    :param model: wild-type model
    :param knockouts: list of reaction knockouts
    :param method: the method to compute the perturbation. default is "fva" - Flux Variability Analysis.
        It can also be computed with "em" - Elementary modes
    :return: perturbation
    """

    if method == "fva":
        _fbip_fva(model, knockouts, view)
    elif method == "em":
        raise NotImplementedError("Elementary modes approach is not implemented")
    else:
        raise ValueError("%s method is not valid to compute FBIP" % method)


def _fbip_fva(model, knockouts, view):
    tm = TimeMachine()
    for reaction in model.reactions:
        if reaction.reversibility:
            tm(do=partial(setattr, reaction, 'upper_bound', 1),
               undo=partial(setattr, reaction, 'upper_bound', reaction.upper_bound))
            tm(do=partial(setattr, reaction, 'lower_bound', -1),
               undo=partial(setattr, reaction, 'lower_bound', reaction.upper_bound))
        else:
            tm(do=partial(setattr, reaction, 'upper_bound', 1),
               undo=partial(setattr, reaction, 'upper_bound', reaction.upper_bound))
            tm(do=partial(setattr, reaction, 'lower_bound', 0),
               undo=partial(setattr, reaction, 'lower_bound', reaction.upper_bound))

    wt_fva = flux_variability_analysis(model, view)
    for reaction in knockouts:
        tm(do=partial(setattr, reaction, 'upper_bound', 0),
           undo=partial(setattr, reaction, 'upper_bound', reaction.upper_bound))
        tm(do=partial(setattr, reaction, 'lower_bound', 0),
           undo=partial(setattr, reaction, 'lower_bound', reaction.upper_bound))

    mt_fva = flux_variability_analysis(model, view)

    perturbation = 0
    for reaction in model.reactions:
        if wt_fva[reaction.id] != 0 and mt_fva[reaction.id] == 0:
            perturbation += 1

    tm.reset()
    return perturbation


def reaction_component_production(model, reaction):
    tm = TimeMachine()
    for metabolite in reaction.metabolites:
        test = Reaction("EX_%s_temp" % metabolite.id)
        test._metabolites[metabolite] = -1
        #hack frozen set from cobrapy to be able to add a reaction
        metabolite._reaction = set(metabolite._reaction)
        tm(do=partial(model.add_reactions, [test]), undo=partial(model.remove_reactions, [test]))
        tm(do=partial(setattr, model, 'objective', test.id), undo=partial(setattr, model, 'objective', model.objective))
        try:
            print metabolite.id, "= ", model.solve().f
        except SolveError:
            print metabolite, " cannot be produced (reactions: %s)" % metabolite.reactions
        finally:
            tm.reset()



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