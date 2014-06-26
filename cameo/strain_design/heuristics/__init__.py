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

from cameo.strain_design.heuristics import archivers
from cameo.strain_design.heuristics import plotters
from cameo import config
from cameo.flux_analysis.simulation import pfba
from cameo.util import partition, TimeMachine


import inspyred
import logging
from ordered_set import OrderedSet

from functools import partial
from random import Random



from cobra.manipulation.delete import find_gene_knockout_reactions

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


def _setup(algorithm, variator, selector, replacer, archiver, terminator):
    logger.debug("Setting up algorithm: %s" % algorithm.heuristic_method)
    algorithm.heuristic_method.variator = variator
    algorithm.heuristic_method.selector = selector
    algorithm.heuristic_method.replacer = replacer
    algorithm.heuristic_method.archiver = archiver
    algorithm.heuristic_method.terminator = terminator


PRE_CONFIGURED = {
    inspyred.ec.GA: lambda algorithm:
    _setup(
        algorithm,
        [
            inspyred.ec.variators.crossovers.n_point_crossover,
            algorithm._mutator
        ],
        inspyred.ec.selectors.rank_selection,
        inspyred.ec.replacers.generational_replacement,
        archivers.BestSolutionArchiver(),
        algorithm.termination
    ),

    inspyred.ec.SA: lambda algorithm:
    _setup(
        algorithm,
        [algorithm._mutator],
        inspyred.ec.selectors.default_selection,
        inspyred.ec.replacers.simulated_annealing_replacement,
        archivers.BestSolutionArchiver(),
        algorithm.termination
    ),
    inspyred.ec.emo.NSGA2: lambda algorithm:
        _setup(
            algorithm,
            [
                algorithm._mutator,
                inspyred.ec.variators.crossovers.n_point_crossover],
            inspyred.ec.selectors.tournament_selection,
            inspyred.ec.replacers.nsga_replacement,
            inspyred.ec.archivers.population_archiver,
            algorithm.termination
        ),
    inspyred.ec.emo.PAES: lambda algorithm:
        _setup(
            algorithm,
            [
                algorithm._mutator,
                inspyred.ec.variators.crossovers.n_point_crossover
            ],
            inspyred.ec.selectors.default_selection,
            inspyred.ec.replacers.paes_replacement,
            inspyred.ec.archivers.adaptive_grid_archiver,
            algorithm.termination
        )
}


class HeuristicOptimization(object):
    def __init__(self, model=None, heuristic_method=inspyred.ec.GA, objective_function=None, random=None,
                 termination=inspyred.ec.terminators.evaluation_termination, *args, **kwargs):
        super(HeuristicOptimization, self).__init__(*args, **kwargs)
        if random is None:
            random = Random()
        self.random = random
        self.model = model
        self.termination = termination
        self._objective_function = objective_function
        self._heuristic_method = heuristic_method
        self.heuristic_method = heuristic_method

    @property
    def objective_function(self):
        return self._objective_function

    @objective_function.setter
    def objective_function(self, objective_function):
        if self._heuristic_method.__module__ == inspyred.ec.ec.__name__ and isinstance(objective_function, list):
            if len(objective_function) == 1:
                self._objective_function = objective_function[0]
            else:
                raise TypeError("single objective heuristics do not support multiple objective functions")
        elif self._heuristic_method.__module__ == inspyred.ec.emo.__name__ and not isinstance(objective_function, list):
            self._objective_function = [objective_function]
        else:
            self._objective_function = objective_function

    @property
    def heuristic_method(self):
        return self._heuristic_method

    @heuristic_method.setter
    def heuristic_method(self, heuristic_method):
        if heuristic_method.__module__ == inspyred.ec.emo.__name__ and not self.is_mo():
            self._objective_function = [self.objective_function]
        elif heuristic_method.__module__ == inspyred.ec.ec.__name__ and self.is_mo():
            if len(self.objective_function) == 1:
                self._objective_function = self.objective_function[0]
            else:
                raise TypeError("single objective heuristics do not support multiple objective functions")
        self._heuristic_method = heuristic_method(self.random)

    def run(self, view=config.default_view, **kwargs):
        return self.heuristic_method.evolve(
            generator=self._generate_individual,
            maximize=True,
            view=view,
            evaluator=self._evaluate_population,
            **kwargs)

    def is_mo(self):
        return isinstance(self.objective_function, list)


class _ChunkEvaluation(object):
    def __init__(self, optimization):
        self.optimization = optimization
        self.time_machine = TimeMachine()

    def __call__(self, population):
        return [self.optimization.evaluate_individual(i, self.time_machine) for i in population]


class KnockoutOptimization(HeuristicOptimization):
    def __init__(self, simulation_method=pfba, *args, **kwargs):
        super(KnockoutOptimization, self).__init__(*args, **kwargs)
        self.simulation_method = simulation_method
        self.representation = None
        self.ko_type = None

    def _distance_function(self, candidate1, candidate2):
        return len(set(candidate1).symmetric_difference(set(candidate2)))

    def _decoder(self, individual):
        raise NotImplementedError

    def _generate_individual(self, random, args):
        max_size = args.get('max_size', 9)
        individual = random.sample(xrange(len(self.representation)), random.randint(1, max_size))
        return individual

    def _mutator(self, random, candidates, args):
        candidates = [self._mutation(random, candidate, args) for candidate in candidates]
        return [self._indel(random, candidate, args) for candidate in candidates]

    def _mutation(self, random, individual, args):
        new_individual = []
        for index in individual:
            if random.random() < args.get('mutation_rate', .1):
                new_individual.append(random.randint(0, len(self.representation) - 1))
            else:
                new_individual.append(index)

        return new_individual

    def _indel(self, random, individual, args):
        if random.random() < args.get('mutation_rate', .1):
            if random.random() > 0.5:
                if len(individual) > 1:
                    individual.pop(random.randint(0, len(individual) - 1))
            else:
                individual.append(random.sample(xrange(len(self.representation)), 1)[0])

        return list(OrderedSet(individual))

    def evaluate_individual(self, individual, tm):
        decoded = self._decoder(individual)
        reactions = decoded[0]

        for reaction in reactions:
            tm(do=partial(setattr, reaction, 'lower_bound', 0),
               undo=partial(setattr, reaction, 'lower_bound', reaction.lower_bound))
            tm(do=partial(setattr, reaction, 'upper_bound', 0),
               undo=partial(setattr, reaction, 'upper_bound', reaction.upper_bound))

        try:
            solution = self.simulation_method(self.model)
            fitness = self._calculate_fitness(solution, decoded)
        except Exception, e:
            logger.exception(e)
            if isinstance(self.objective_function, list):
                fitness = inspyred.ec.emo.Pareto(values=[0 for of in self.objective_function])
            else:
                fitness = 0

        tm.reset()

        return fitness

    def _calculate_fitness(self, solution, decoded):
        if self.is_mo():
            logger.debug("evaluate multi objective")
            return inspyred.ec.emo.Pareto(values=[of(self.model, solution, decoded) for of in self.objective_function])
        else:
            logger.debug("evaluate single objective")
            return self.objective_function(self.model, solution, decoded)

    def _evaluate_population(self, candidates, args):
        view = args.get('view')
        population_chunks = (chunk for chunk in partition(candidates, len(view)))
        func_obj = _ChunkEvaluation(self)
        results = view.map(func_obj, population_chunks)
        fitness = reduce(list.__add__, results)

        return fitness

    @HeuristicOptimization.heuristic_method.setter
    def heuristic_method(self, heuristic_method):
        HeuristicOptimization.heuristic_method.fset(self, heuristic_method)
        self._set_observer()
        try:
            PRE_CONFIGURED[heuristic_method](self)
        except KeyError:
            logger.warning("Please verify the variator is compatible with set representation")

    @HeuristicOptimization.objective_function.setter
    def objective_function(self, objective_function):
        HeuristicOptimization.objective_function.fset(self, objective_function)
        self._set_observer()

    def _set_observer(self):
        if self.is_mo():
            self.observer = plotters.ParetoPlotObserver(ofs=self.objective_function)
        else:
            self.observer = plotters.PlotObserver()

    def run(self, **kwargs):
        self.observer.reset()
        self.heuristic_method.observer = self.observer
        super(KnockoutOptimization, self).run(distance_function=self._distance_function, **kwargs)
        return KnockoutOptimizationResult(self.model, self.heuristic_method, self.simulation_method,
                                          self.heuristic_method.archive, self.objective_function,
                                          self.ko_type, self._decoder)


class KnockoutOptimizationResult(object):
    class KnockoutOptimizationSolution(object):
        def __init__(self, solution, model, simulation_method, decoder, *args, **kwargs):
            super(KnockoutOptimizationResult.KnockoutOptimizationSolution, self).__init__(*args, **kwargs)
            decoded_solution = decoder(solution.candidate)
            self.knockout_list = decoded_solution[1]
            result = self._simulate(decoded_solution[0], simulation_method, model)
            self.biomass = result.f
            self.fitness = solution.fitness

        def _simulate(self, reactions, method, model):
            tm = TimeMachine()
            for reaction in reactions:
                tm(do=partial(setattr, reaction, 'lower_bound', 0),
                   undo=partial(setattr, reaction, 'lower_bound', reaction.lower_bound))
                tm(do=partial(setattr, reaction, 'upper_bound', 0),
                   undo=partial(setattr, reaction, 'upper_bound', reaction.upper_bound))

            try:
                solution = method(model)
            except Exception, e:
                logger.exception(e)

            tm.reset()
            return solution

        def html_row(self):
            return "<tr>" \
                   "    <td>" + ", ".join([ko.id for ko in self.knockout_list]) + "</td>" \
                   "    <td>" + str(self.fitness) + "</td>" \
                   "    <td>" + str(self.biomass) + "</td>" \
                   "</tr>"

    def __init__(self, model, heuristic_method, simulation_method, solutions, objective_function,
                 ko_type, decoder, *args, **kwargs):
        super(KnockoutOptimizationResult, self).__init__(*args, **kwargs)
        self.model = model
        self.heuristic_method = heuristic_method
        self.simulation_method = simulation_method
        if isinstance(objective_function, list):
            self.objective_functions = objective_function
        else:
            self.objective_functions = [objective_function]
        self.ko_type = ko_type
        self.decoder = decoder
        self.solutions = [self._build_solution(s, model, simulation_method, decoder) for s in solutions]

    def _build_solution(self, solution, model, simulation_method, decoder):
        return KnockoutOptimizationResult.KnockoutOptimizationSolution(
            solution,
            model,
            simulation_method,
            decoder
        )

    def _repr_html_(self):
        results = "<h4>Result:</h4>" \
                  "<ul>" \
                  "    <li>model: " + self.model.id + "</li>" \
                  "    <li>heuristic: " + self.heuristic_method.__class__.__name__ + "</li>" \
                  "    <li>objective function: " + "|".join([o.__name__ for o in self.objective_functions]) + "</li>" \
                  "    <li>simulation method: " + self.simulation_method.__name__ + "</li>" \
                  "    <li>type: " + self.ko_type + "</li>" \
                  "</ul>" \
                  "<h4>Solutions:</h4>"
        solutions = "<table>" \
                    "   <thead>" \
                    "       <tr>" \
                    "           <th>Knockouts</th>" \
                    "           <th>Fitness</th>" \
                    "           <th>Biomass</th>" \
                    "       </tr>" \
                    "   </thead>" \
                    "   <tbody>" \
                    "       " + "\n".join([solution.html_row() for solution in self.solutions]) + "" \
                    "   </tbody>" \
                    "<table>" \

        return results + solutions


class ReactionKnockoutOptimization(KnockoutOptimization):
    def __init__(self, reactions=None, *args, **kwargs):
        super(ReactionKnockoutOptimization, self).__init__(*args, **kwargs)
        if reactions is None:
            reactions = set([r.id for r in self.model.reactions])

        essential_reactions = set([r.id for r in self.model.essential_reactions()])
        exchange_reactions = set([r.id for r in self.model.exchanges])
        self.representation = list(reactions.difference(essential_reactions).difference(exchange_reactions))
        self.ko_type = 'reaction'

    def _decoder(self, individual):
        reactions = [self.model.reactions.get_by_id(self.representation[index]) for index in individual]
        return [reactions, reactions]


class GeneKnockoutOptimization(KnockoutOptimization):
    def __init__(self, genes=None, *args, **kwargs):
        super(GeneKnockoutOptimization, self).__init__(*args, **kwargs)
        if genes is None:
            genes = set([g.id for g in self.model.genes])

        essential_genes = set([g.id for g in self.model.essential_genes()])
        self.representation = list(genes.difference(essential_genes))
        self.ko_type = 'gene'


    def _decoder(self, individual):
        genes = [self.model.genes.get_by_id(self.representation[index]) for index in individual]
        reactions = [find_gene_knockout_reactions(self.model, [gene]) for gene in genes]
        reactions = reduce(list.__add__, reactions)
        return [reactions, genes]