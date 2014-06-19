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
import logging

from cameo import config
from cameo.flux_analysis.simulation import pfba
from cameo.strain_design.heuristics import BestSolutionPool, PlotObserver
from cameo.util import partition, TimeMachine

from cobra.manipulation.delete import find_gene_knockout_reactions

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


class HeuristicOptimization(object):
    def __init__(self, model=None, heuristic_method=inspyred.ec.GA, objective_function=None, random=None, *args,
                 **kwargs):
        super(HeuristicOptimization, self).__init__(*args, **kwargs)
        if random is None:
            random = Random()
        self.model = model
        self.heuristic_method = heuristic_method(random)
        self.objective_function = objective_function

    def run(self, **kwargs):
        return self.heuristic_method.evolve(**kwargs)


class _ChunkEvaluation(object):
    def __init__(self, optimization):
        self.optimization = optimization
        self.time_machine = TimeMachine()

    def __call__(self, population):
        return [self.optimization._evaluate_individual(i, self.time_machine) for i in population]


class KnockoutOptimization(HeuristicOptimization):
    def __init__(self, simulation_method=pfba, max_size=100, *args, **kwargs):
        super(KnockoutOptimization, self).__init__(*args, **kwargs)
        self.simulation_method = simulation_method
        self.solution_pool = BestSolutionPool(max_size)
        self.observer = PlotObserver()

    def _decode_individual(self, individual):
        raise NotImplementedError

    def _generate_individual(self, random, args):
        raise NotImplementedError

    def _mutator(self, random, candidates, args):
        return [self._mutation(random, candidate, args) for candidate in candidates]

    def _mutation(self, random, individual, args):
        raise NotImplementedError

    def _evaluate_individual(self, individual, tm):

        reactions = self._decode_individual(individual)

        for reaction in reactions:
            tm(do=partial(setattr, reaction, 'lower_bound', 0),
               undo=partial(setattr, reaction, 'lower_bound', reaction.lower_bound))
            tm(do=partial(setattr, reaction, 'upper_bound', 0),
               undo=partial(setattr, reaction, 'upper_bound', reaction.upper_bound))

        try:
            solution = self.simulation_method(self.model)
            fitness = self.objective_function(solution)
        except Exception, e:
            fitness = 0

        tm.reset()

        return fitness

    def _evaluate_population(self, candidates, args):
        view = args.get('view')
        population_chunks = (chunk for chunk in partition(candidates, len(view)))
        func_obj = _ChunkEvaluation(self)
        results = view.map(func_obj, population_chunks)
        fitnesses = reduce(list.__add__, results)
        for individual, fitness in zip(candidates, fitnesses):
            self.solution_pool.add(individual, fitness)

        return fitnesses

    def run(self, view=config.default_view, **kwargs):
        self.heuristic_method.observer = self.observer
        self.heuristic_method.variator = [inspyred.ec.variators.n_point_crossover, self._mutator]
        self.heuristic_method.terminator = inspyred.ec.terminators.evaluation_termination
        print kwargs
        final_population = self.heuristic_method.evolve(
            generator=self._generate_individual,
            maximize=True,
            view=view,
            evaluator=self._evaluate_population,
            **kwargs)
        return final_population


class ReactionKnockoutOptimization(KnockoutOptimization):
    def __init__(self, reactions=None, *args, **kwargs):
        super(ReactionKnockoutOptimization, self).__init__(*args, **kwargs)
        if reactions is None:
            reactions = set([r.id for r in self.model.reactions])

        essential_reactions = set([r.id for r in self.model.essential_reactions()])
        exchange_reactions = set([r.id for r in self.model.exchanges])
        self.reactions = list(reactions.difference(essential_reactions).difference(exchange_reactions))

    def _decode_individual(self, individual):
        return [self.model.reactions.get_by_id(self.reactions[index]) for index in individual]

    def _generate_individual(self, random, args):
        max_size = args.get('max_size', 9)
        individual = random.sample(xrange(len(self.reactions)), random.randint(1, max_size))
        return individual

    def _mutation(self, random, individual, args):
        new_individual = list()
        for index in individual:
            if random.random() < args.get('mutation_rate', .1):
                new_individual.append(random.randint(0, len(self.reactions) - 1))
            else:
                new_individual.append(index)

        return list(set(new_individual))


class GeneKnockoutOptimization(KnockoutOptimization):
    def __init__(self, genes=None, *args, **kwargs):
        super(GeneKnockoutOptimization, self).__init__(*args, **kwargs)
        if genes is None:
            genes = set([g.id for g in self.model.genes])

        essential_genes = set([g.id for g in self.model.essential_genes()])
        self.genes = list(genes.difference(essential_genes))

    def _decode_individual(self, individual):
        genes = [self.model.genes[index] for index in individual]
        reactions = [find_gene_knockout_reactions([gene]) for gene in genes]
        return reduce(list.__add__, reactions)

    def _generate_individual(self, random, args):
        max_size = args.get('max_size', 9)
        individual = random.sample(xrange(len(self.model.genes)), random.randint(1, max_size))
        return individual

    def _mutation(self, random, individual, args):
        new_individual = list()
        for index in individual:
            if random.random() < args.get('mutation_rate', .1):
                new_individual.append(random.randint(0, len(self.model.genes) - 1))
            else:
                new_individual.append(index)

            return list(set(new_individual))