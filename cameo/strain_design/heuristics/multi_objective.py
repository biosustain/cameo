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
from functools import partial
from random import Random
import inspyred
from cameo.flux_analysis.simulation import pfba


class HeuristicOptimization(object):
    def __init__(self, model=None, heuristic_method=inspyred.ec.emo.NSGA, objective_functions=None, random=None,
                 termination=inspyred.ec.terminators.evaluation_termination, *args, **kwargs):
        super(HeuristicOptimization, self).__init__(*args, **kwargs)
        if random is None:
            random = Random()
        self.random = random
        self.model = model
        self.termination = termination
        self.heuristic_method = heuristic_method
        self.objective_functions = objective_functions

    @property
    def heuristic_method(self):
        return self._heuristic_method

    @heuristic_method.setter
    def heuristic_method(self, heuristic_method):
        self._heuristic_method = heuristic_method(self.random)

    def run(self, **kwargs):
        return self.heuristic_method.evolve(**kwargs)


class KnockoutOptimization(HeuristicOptimization):
    def __init__(self, simulation_method=pfba, max_size=100, *args, **kwargs):
        super(KnockoutOptimization, self).__init__(*args, **kwargs)
        self.simulation_method = simulation_method
        self.solution_pool = None
        self.observer = None
        self.representation = None
        self.ko_type = None

    def _evaluate_individual(self, individual, tm, args):
        reactions = self._decode_individual(individual)[0]

        for reaction in reactions:
            tm(do=partial(setattr, reaction, 'lower_bound', 0),
               undo=partial(setattr, reaction, 'lower_bound', reaction.lower_bound))
            tm(do=partial(setattr, reaction, 'upper_bound', 0),
               undo=partial(setattr, reaction, 'upper_bound', reaction.upper_bound))

        try:
            solution = self.simulation_method(self.model)
            fitness = [of(solution) for of in self.objective_functions]
        except Exception:
            fitness = [0 for of in self.objective_functions]

        tm.reset()

        return fitness

        return inspyred.ec.emo.Pareto(fitness)