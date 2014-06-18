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
from random import sample
import inspyred
from cameo import config

from cameo.flux_analysis.simulation import pfba
from cameo.strain_design.heuristics import BestSolutionPool
from cameo.util import partition, TimeMachine

from cobra.manipulation.delete import find_gene_knockout_reactions


class HeuristicOptimization(object):
    def __init__(self, model=None, heuristic_method=inspyred.ec.GA, objective_function=None, *args, **kwargs):
        super(HeuristicOptimization, self).__init__(*args, **kwargs)
        self.model = model
        self.heuristic_method = heuristic_method
        self.objective_function = objective_function

    def run(self, **kwargs):
        return self.heuristic_method.evolve(**kwargs)


class KnockoutOptimization(HeuristicOptimization):

    class _ChunkEval(object):

        def __init__(self, optimization):
            self.optimization = optimization
            self.time_machine = TimeMachine()

        def __call__(self, population):
            return [self.optimization._evaluate_individual(i, self.time_machine) for i in population]

    def __init__(self, simulation_method=pfba, max_size=100, *args, **kwargs):
        super(HeuristicOptimization, self).__init__(*args, **kwargs)
        self.simulation_method = simulation_method
        self.solution_pool = BestSolutionPool(max_size)
       #self.observer = BestSolutionObserver()

    def _decode_individual(self, individual):
        raise NotImplementedError

    def _generate_population(self, args):
        raise NotImplementedError

    @inspyred.ec.variators.mutator
    def _mutation(self, random, individual, args):
        raise NotImplementedError

    def _evaluate_individual(self, individual, tm):

        reactions = self._decode_individual(individual)

        for reaction in reactions:
            tm(do=partial(setattr, reaction, 'lower_bound', 0),
               undo=partial(setattr, reaction, 'lower_bound', reaction.lower_bound))
            tm(do=partial(setattr, reaction, 'upper_bound', 0),
               undo=partial(setattr, reaction, 'upper_bound', reaction.lower_bound))

        try:
            solution = self.simulation_method(self.model)
            fitness = self.objective_function(solution)
        except Exception:
            fitness = 0

        tm.reset()

        return fitness

    def _evaluate_population(self, population, args):
        view = args.get('view')
        population_chunks = (chunk for chunk in partition(population, len(view)))
        func_obj = self._ChunkEval(self)
        results = view.map(func_obj, population_chunks)
        return reduce(list.__add__, results)

    def run(self, view=config.default_view, **kwargs):
        self.heuristic_method.generator = self._generate_population
        self.heuristic_method.evaluator = self._evaluate_population
        self.heuristic_method.observer = [self.observer]
        self.heuristic_method.variator = [inspyred.ec.variators.n_point_crossover, self._mutation]
        final_population = self.heuristic_method.evolve(view=view, **kwargs)
        return [final_population]


class ReactionKnockoutOptimization(KnockoutOptimization):
    def __init__(self, *args, **kwargs):
        super(HeuristicOptimization, self).__init__(*args, **kwargs)

    def _decode_individual(self, individual):
        return [self.model.reactions[index] for index in individual]

    def _generate_population(self, args):
        max_size = args.get('max_size', 9)
        individual = sample(xrange(len(self.model.reactions)), self.heuristic_method.random.randint(1, max_size))
        return individual

    def custom_mutation(self, random, individual, args):
        new_individual = list()
        for index in individual:
            if random.random() < args.get('mutation_rate', .1):
                new_individual.append(random.randint(0, len(self.model.reactions)-1))
            else:
                new_individual.append(index)

        return list(set(new_individual))


class GeneKnockoutOptimization(KnockoutOptimization):
    def __init__(self, *args, **kwargs):
        super(HeuristicOptimization, self).__init__(*args, **kwargs)

    def _decode_individual(self, individual):
        genes = [self.model.genes[index] for index in individual]
        reactions = [find_gene_knockout_reactions(gene) for gene in genes]
        return reduce(list.__add__, reactions)

    def _generate_population(self, args):
        max_size = args.get('max_size', 9)
        individual = sample(xrange(len(self.model.genes)), self.heuristic_method.random.randint(1, max_size))
        return individual

    def custom_mutation(self, random, individual, args):
        new_individual = list()
        for index in individual:
            if random.random() < args.get('mutation_rate', .1):
                new_individual.append(random.randint(0, len(self.model.genes)-1))
            else:
                new_individual.append(index)

            return list(set(new_individual))