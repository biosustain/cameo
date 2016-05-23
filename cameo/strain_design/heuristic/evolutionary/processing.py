# Copyright 2016 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
from functools import partial

from cameo import flux_variability_analysis
from cameo.util import TimeMachine


def process_knockout_solution(model, solution, simulation_method, simulation_kwargs,
                              biomass, target, substrate, objective_functions, cache=None):
    """

    Arguments
    ---------

    model: SolverBasedModel
        A constraint-based model
    solution: tuple - (reactions, knockouts)
        The output of a decoder
    simulation_method: function
        See see cameo.flux_analysis.simulation
    simulation_kwargs: dict
        Keyword arguments to run the simulation method
    biomass: Reaction
        Cellular biomass reaction
    target: Reaction
        The strain design target
    substrate: Reaction
        The main carbon source uptake rate
    objective_functions: list
        A list of cameo.strain_design.heuristic.evolutionary.objective_functions.ObjectiveFunction
    cache: ProblemCache
        A problem cache for performance improvment

    Returns
    -------

    list
        A list with: reactions, knockouts, size, fva_min, fva_max, target flux, biomass flux, yield, fitness,
        [fitness, [fitness]]
    """

    with TimeMachine() as tm:
        for ko in solution[0]:
            model.reactions.get_by_id(ko).knock_out(tm)

        reactions = reactions2filter(objective_functions)
        flux_dist = simulation_method(model, cache=cache, reactions=reactions, objective=biomass, **simulation_kwargs)
        tm(do=partial(setattr, model, "objective", biomass),
           undo=partial(setattr, model, "objective", model.objective))

        fva = flux_variability_analysis(model, fraction_of_optimum=0.99, reactions=[target])
        target_yield = flux_dist[target] / abs(flux_dist[substrate])
        return [solution[0], solution[1], len(solution[1]), fva.lower_bound(target),
                fva.upper_bound(target), flux_dist[target], flux_dist[biomass],
                target_yield] + [of(model, flux_dist, solution) for of in objective_functions]


def reactions2filter(objective_function):
    """
    Retrieve from the solvers memory only the reactions required by the objective functions. This is faster then
    reading all reactions from the solver output.
    """
    if isinstance(objective_function, list):
        reactions = []
        [reactions.extend(of.reactions) for of in objective_function]
    else:
        reactions = objective_function.reactions
    return set(reactions)
