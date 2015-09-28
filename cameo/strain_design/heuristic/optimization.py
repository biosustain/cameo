# Copyright 2015 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import absolute_import, print_function
from cameo.strain_design import StrainDesignResult, StrainDesignMethod, StrainDesign

__all__ = ['ReactionKnockoutOptimization', 'GeneKnockoutOptimization']

from six.moves import range

import time
import numpy as np

from functools import reduce
from sympy import Symbol

from cobra.manipulation.delete import find_gene_knockout_reactions
from inspyred.ec.emo import Pareto

from cameo.exceptions import SolveError
from cameo.visualization.sympy_ext import And, Or
from cameo.strain_design.heuristic import archivers
from cameo.strain_design.heuristic import plotters
from cameo.strain_design.heuristic import observers
from cameo.strain_design.heuristic import variators
from cameo.strain_design.heuristic import generators
from cameo.strain_design.heuristic import decoders
from cameo.strain_design.heuristic import stats
from cameo import config
from cameo.flux_analysis.simulation import pfba, lmoma, moma, room
from cameo.util import partition, TimeMachine, memoize, ProblemCache
from pandas import DataFrame

import inspyred
import logging

from cameo.util import RandomGenerator as Random
from cameo.util import in_ipnb
from cameo.visualization import draw_knockout_result

REACTION_KNOCKOUT_TYPE = "reaction"
GENE_KNOCKOUT_TYPE = "gene"

SIZE = 'Size'
BIOMASS = 'Biomass'
KNOCKOUTS = 'Knockouts'
REACTIONS = 'Reactions'

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

PRE_CONFIGURED = {
    inspyred.ec.GA: [
        [
            variators.set_mutation,
            variators.set_indel,
            variators.set_n_point_crossover
        ],
        inspyred.ec.selectors.tournament_selection,
        inspyred.ec.replacers.generational_replacement,
        archivers.BestSolutionArchiver(),
    ],
    inspyred.ec.SA: [

        [
            variators.set_mutation,
            variators.set_indel
        ],
        inspyred.ec.selectors.default_selection,
        inspyred.ec.replacers.simulated_annealing_replacement,
        archivers.BestSolutionArchiver()
    ],
    inspyred.ec.emo.NSGA2: [
        [
            variators.set_mutation,
            variators.set_indel,
            variators.set_n_point_crossover
        ],
        inspyred.ec.selectors.tournament_selection,
        inspyred.ec.replacers.nsga_replacement,
        inspyred.ec.archivers.population_archiver
    ],
    inspyred.ec.emo.PAES: [
        [
            variators.set_mutation,
            variators.set_indel,
            variators.set_n_point_crossover
        ],
        inspyred.ec.selectors.default_selection,
        inspyred.ec.replacers.paes_replacement,
        inspyred.ec.archivers.adaptive_grid_archiver
    ]
}


def set_distance_function(candidate1, candidate2):
    return len(set(candidate1).symmetric_difference(set(candidate2)))


class HeuristicOptimization(StrainDesignMethod):
    """
    Blueprint for any model optimization based on heuristic methods.


    Attributes
    ----------
    model : SolverBasedModel
        A constraint-based model.
    heuristic_method : inspyred.ec.EvolutionaryComputation
        An evolutionary algorithm.
    objective_function : objective function or list(objective function)
        The objectives for the algorithm to maximize.
    seed : int
        A seed for random. It is auto-generated if None is given.
    termination : inspyred.ec.terminators
        A termination criteria for the algorithm. The default is inspyred.ec.terminators.evaluation_termination


    Methods
    -------
    run(view=config.default_view, maximize=True, **kwargs)


    See Also
    --------
    *inspyred.ec
    *cameo.config.default_view

    """

    def __init__(self, model=None, heuristic_method=inspyred.ec.GA, objective_function=None, seed=None,
                 termination=inspyred.ec.terminators.evaluation_termination, plot=True, progress=True,
                 *args, **kwargs):
        super(HeuristicOptimization, self).__init__(*args, **kwargs)
        logger.debug("Seed: %s" % seed)
        self.plot = plot
        self.progress = progress
        if seed is None:
            seed = int(time.time())
        self.seed = seed
        self.observers = []
        self.random = Random(seed=seed)
        self.model = model
        self.termination = termination
        self._objective_function = objective_function
        self._heuristic_method = None
        self.heuristic_method = heuristic_method
        self.heuristic_method.terminator = termination
        self._generator = None

    @property
    def objective_function(self):
        return self._objective_function

    @objective_function.setter
    def objective_function(self, objective_function):
        if self._heuristic_method.__module__ == inspyred.ec.ec.__name__ and isinstance(objective_function, list):
            if len(objective_function) == 1:
                self._objective_function = objective_function[0]
            else:
                raise TypeError("single objective heuristic do not support multiple objective functions")
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

    def _evaluator(self, candidates, args):
        raise NotImplementedError

    def run(self, view=config.default_view, maximize=True, **kwargs):
        for observer in self.observers:
            observer.reset()
        t = time.time()
        print(time.strftime("Starting optimization at %a, %d %b %Y %H:%M:%S", time.localtime(t)))
        res = self.heuristic_method.evolve(generator=self._generator,
                                           maximize=maximize,
                                           view=view,
                                           evaluator=self._evaluator,
                                           **kwargs)
        for observer in self.observers:
            observer.end()
        runtime = time.time() - t
        print(time.strftime("Finished after %H:%M:%S", time.gmtime(runtime)))

        return res

    def is_mo(self):
        return isinstance(self.objective_function, list)


class KnockoutEvaluator(object):
    """
    Evaluator for knockouts. It gets a representation, decodes into knockouts and simulates
    the model with a given method.

    Attributes
    ----------
    model : SolverBasedModel
        A constraint-based model
    decoder : KnockoutDecoder
        A decoder to convert the representation into knockouts
    objective_function : objective_function or list(objective_function)
        The objectives of the algorithm
    simulation_method : see flux_analysis.simulation
        The method use to simulate the knockouts
    simulation_kwargs : dict
        The extra parameters used by the simulation method

    See Also
    --------
    cameo.strain_design.heuristic.objective_function

    Methods
    -------
    __call__(population)
        calling the object will evaluate a population (see inspyred)

    """

    def __init__(self, model, decoder, objective_function, simulation_method, simulation_kwargs):
        self.model = model
        self.decoder = decoder
        self.objective_function = objective_function
        self.simulation_method = simulation_method
        self.simulation_kwargs = simulation_kwargs
        self.cache = ProblemCache(model)

    def __call__(self, population):
        res = [self.evaluate_individual(frozenset(i)) for i in population]
        self.cache.reset()
        return res

    @memoize
    def evaluate_individual(self, individual):
        decoded = self.decoder(individual)
        reactions = decoded[0]
        with TimeMachine() as tm:
            for reaction in reactions:
                reaction.knock_out(time_machine=tm)
            try:
                solution = self.simulation_method(self.model,
                                                  cache=self.cache,
                                                  volatile=False,
                                                  raw=True,
                                                  reactions=self._reactions2filter(),
                                                  **self.simulation_kwargs)
                fitness = self._calculate_fitness(solution, decoded)
            except SolveError as e:
                logger.debug(e)
                if isinstance(self.objective_function, list):
                    fitness = inspyred.ec.emo.Pareto(values=[0 for _ in self.objective_function])
                else:
                    fitness = 0

            return fitness

    def _reactions2filter(self):
        if isinstance(self.objective_function, list):
            reactions = []
            [reactions.extend(of.reactions) for of in self.objective_function]
        else:
            reactions = self.objective_function.reactions
        return set(reactions)

    def _calculate_fitness(self, solution, decoded):
        if isinstance(self.objective_function, list):
            logger.debug("evaluate multi objective")
            return inspyred.ec.emo.Pareto(values=[of(self.model, solution, decoded) for of in self.objective_function])
        else:
            logger.debug("evaluate single objective")
            return self.objective_function(self.model, solution, decoded)


class KnockoutOptimization(HeuristicOptimization):
    """
    Abstract class for knockout optimization.
    """

    def __init__(self, simulation_method=pfba, wt_reference=None, *args, **kwargs):
        """
         Attributes
        ----------
        same as HeuristicOptimization
        simulation_method: see flux_analysis.simulation
        wt_reference: dict
        """
        super(KnockoutOptimization, self).__init__(*args, **kwargs)
        self.wt_reference = wt_reference
        self._simulation_method = None
        self.simulation_method = simulation_method
        self.representation = None
        self._ko_type = None
        self._decoder = None
        self._generator = generators.set_generator

    @property
    def simulation_method(self):
        return self._simulation_method

    @simulation_method.setter
    def simulation_method(self, simulation_method):
        if simulation_method in [lmoma, moma, room] and self.wt_reference is None:
            logger.info("No WT reference found, generating using pfba.")
            self.wt_reference = pfba(self.model).fluxes
        self._simulation_method = simulation_method

    def _evaluator(self, candidates, args):
        view = args.get('view')
        population_chunks = (chunk for chunk in partition(candidates, len(view)))
        kwargs = {'reference': self.wt_reference}

        func_obj = KnockoutEvaluator(self.model, self._decoder, self.objective_function, self.simulation_method, kwargs)
        try:
            results = view.map(func_obj, population_chunks)
        except KeyboardInterrupt as e:
            view.shutdown()
            raise e

        fitness = reduce(list.__add__, results)

        return fitness

    @HeuristicOptimization.heuristic_method.setter
    def heuristic_method(self, heuristic_method):
        HeuristicOptimization.heuristic_method.fset(self, heuristic_method)
        self._set_observer()
        try:
            configuration = PRE_CONFIGURED[heuristic_method]
            self._setup(*configuration)
        except KeyError:
            logger.warning("Please verify the variator is compatible with set representation")

    def _setup(self, variator, selector, replacer, archiver):
        logger.debug("Setting up algorithm: %s" % self.heuristic_method)
        self.heuristic_method.variator = variator
        self.heuristic_method.selector = selector
        self.heuristic_method.replacer = replacer
        self.heuristic_method.archiver = archiver
        self.heuristic_method.terminator = self.termination

    @HeuristicOptimization.objective_function.setter
    def objective_function(self, objective_function):
        HeuristicOptimization.objective_function.fset(self, objective_function)
        self._set_observer()

    def _set_observer(self):
        self.observers = []

        if in_ipnb() and self.plot:
            if config.use_bokeh:
                if self.is_mo():
                    self.observers.append(plotters.IPythonBokehParetoPlotter(self.objective_function))
                else:
                    self.observers.append(plotters.IPythonBokehFitnessPlotter())
            elif config.use_matplotlib:
                pass
            else:
                pass

        else:
            if config.use_bokeh:
                pass
            else:
                pass
        if self.progress:
            self.observers.append(observers.ProgressObserver())

    def run(self, **kwargs):
        """
        Parameters
        ----------
        max_size: int
            Maximum size of a solution.
        variable_size: boolean
            If true, the solution size can change meaning that the combination of knockouts can have different sizes up to
            max_size. Otherwise it only produces knockout solutions with a fixed number of knockouts.

        """
        self.heuristic_method.observer = self.observers
        super(KnockoutOptimization, self).run(
            distance_function=set_distance_function,
            representation=self.representation,
            **kwargs)
        return KnockoutOptimizationResult(model=self.model,
                                          heuristic_method=self.heuristic_method,
                                          simulation_method=self.simulation_method,
                                          solutions=self.heuristic_method.archive,
                                          objective_function=self.objective_function,
                                          ko_type=self._ko_type,
                                          decoder=self._decoder,
                                          product=kwargs.get('product', None),
                                          biomass=kwargs.get('biomass', None),
                                          seed=self.seed,
                                          reference=self.wt_reference)


# TODO: Figure out a way to provide generic parameters for different simulation methods
class KnockoutOptimizationResult(StrainDesignResult):
    @staticmethod
    def merge(a, b):
        return a._merge(b)

    def __init__(self, model=None, heuristic_method=None, simulation_method=None, solutions=None,
                 objective_function=None, ko_type=None, decoder=None, product=None, biomass=None,
                 seed=None, reference=None, *args, **kwargs):
        super(KnockoutOptimizationResult, self).__init__(*args, **kwargs)
        self.biomass = biomass
        self.seed = seed
        self.reference = reference
        if product is None:
            self.products = []
        elif isinstance(product, str):
            self.products = [product]
        else:
            self.products = product

        self.model = model
        self.heuristic_method = heuristic_method
        self.simulation_method = simulation_method
        if isinstance(objective_function, list):
            self.objective_functions = objective_function
        else:
            self.objective_functions = [objective_function]
        self.ko_type = ko_type
        self.decoder = decoder
        self.solutions = self._build_solutions(solutions)

    def apply(self, column, function, *args, **kwargs):
        self.solutions[column].apply(function, *args, **kwargs)

    def __getstate__(self):
        return {'product': self.products,
                'model': self.model,
                'biomass': self.biomass,
                'reference': self.reference,
                'simulation_method': self.simulation_method,
                'heuristic_method.__class__': self.heuristic_method.__class__,
                'heuristic_method.maximize': self.heuristic_method.maximize,
                'heuristic_method.variator': self.heuristic_method.variator,
                'heuristic_method.terminator': self.heuristic_method.terminator,
                'heuristic_method.archiver': self.heuristic_method.archiver,
                'heuristic_method.termination_cause': self.heuristic_method.termination_cause,
                'heuristic_method._random': self.heuristic_method._random,
                'heuristic_method.generator': self.heuristic_method.generator,
                'heuristic_method._kwargs.representation': self.heuristic_method._kwargs.get('representation'),
                'heuristic_method._kwargs.max_size': self.heuristic_method._kwargs.get('max_size'),
                'heuristic_method._kwargs.variable_size': self.heuristic_method._kwargs.get('variable_size'),
                'heuristic_method._kwargs.pop_size': self.heuristic_method._kwargs.get('pop_size'),
                'heuristic_method._kwargs.mutation_rate': self.heuristic_method._kwargs.get('mutation_rate'),
                'heuristic_method._kwargs.crossover_rate': self.heuristic_method._kwargs.get('crossover_rate'),
                'heuristic_method._kwargs.num_elites': self.heuristic_method._kwargs.get('num_elites'),
                'objective_functions': self.objective_functions,
                'ko_type': self.ko_type,
                'solutions': self.solutions,
                'seed': self.seed}

    def __setstate__(self, d):
        self.products = d['product']
        self.model = d['model']
        self.biomass = d['biomass']
        self.simulation_method = d['simulation_method']
        self.seed = d['seed']
        self.reference = d['reference']
        random = d['heuristic_method._random']
        self.heuristic_method = d['heuristic_method.__class__'](random)
        self.heuristic_method.maximize = d['heuristic_method.maximize']
        self.heuristic_method.terminator = d['heuristic_method.terminator']
        self.heuristic_method.termination_cause = d['heuristic_method.termination_cause']
        self.heuristic_method.archiver = d['heuristic_method.archiver']
        self.heuristic_method._kwargs['representation'] = d['heuristic_method._kwargs.representation']
        self.heuristic_method._kwargs['max_size'] = d['heuristic_method._kwargs.max_size']
        self.heuristic_method._kwargs['variable_size'] = d['heuristic_method._kwargs.variable_size']
        self.heuristic_method._kwargs['pop_size'] = d['heuristic_method._kwargs.pop_size']
        self.heuristic_method._kwargs['mutation_rate'] = d['heuristic_method._kwargs.mutation_rate']
        self.heuristic_method._kwargs['crossover_rate'] = d['heuristic_method._kwargs.crossover_rate']
        self.heuristic_method._kwargs['num_elites'] = d['heuristic_method._kwargs.num_elites']
        self.objective_functions = d['objective_functions']
        self.ko_type = d['ko_type']
        self.solutions = d['solutions']

    def _build_solutions(self, solutions):
        aggregate = False
        aggregation_functions = {}
        knockouts = []
        biomass = []
        sizes = []
        reactions = []
        fitness = np.zeros((len(solutions), len(self.objective_functions)))
        products = np.zeros((len(solutions), len(self.products)))
        for i, solution in enumerate(solutions):
            mo = isinstance(solution.fitness, Pareto)
            proceed = True if mo else solution.fitness > 0
            if proceed:
                decoded_solution = self.decoder(solution.candidate)
                try:
                    simulation_result = self._simulate(decoded_solution[0])
                except SolveError as e:
                    logger.debug(e)
                    products[i] = [np.nan for _ in self.products]
                    fitness[i] = [np.nan for _ in self.objective_functions]
                    continue
                size = len(decoded_solution[1])

                if self.biomass:
                    biomass.append(simulation_result[self.biomass])

                if not mo:
                    fitness[i] = [solution.fitness]
                else:
                    fitness[i] = solution.fitness

                knockouts.append(And(*[Symbol(v.id) for v in decoded_solution[1]]))
                reactions.append(And(*[Symbol(v.id) for v in decoded_solution[0]]))
                sizes.append(size)
                products[i] = [simulation_result[p] for p in self.products]

        products = products[~np.isnan(products).any(axis=1)]
        fitness = fitness[~np.isnan(fitness).any(axis=1)]

        if self.ko_type == REACTION_KNOCKOUT_TYPE:
            data_frame = DataFrame({KNOCKOUTS: knockouts, SIZE: sizes})
        else:
            data_frame = DataFrame({KNOCKOUTS: knockouts, REACTIONS: reactions, SIZE: sizes})
            aggregation_functions[KNOCKOUTS] = lambda x: Or(*x.values)
            aggregate = True

        if self.biomass is not None:
            data_frame[BIOMASS] = biomass
            aggregation_functions[BIOMASS] = lambda x: x.values[0]

        for j, product in enumerate(self.products):
            data_frame[product] = products[:, j]
            aggregation_functions[product] = lambda x: x.values[0]

        for j, of in enumerate(self.objective_functions):
            data_frame["Fitness %i" % (j+1)] = fitness[:, j]
            aggregation_functions["Fitness %i" % (j+1)] = lambda x: x.values[0]

        if aggregate:
            columns = data_frame.columns
            data_frame = data_frame.groupby([REACTIONS, SIZE], as_index=False).aggregate(aggregation_functions)
            data_frame = data_frame[columns]

        return data_frame

    def _simulate(self, reactions):
        with TimeMachine() as tm:
            for reaction in reactions:
                reaction.knock_out(time_machine=tm)
            solution = self.simulation_method(self.model, reference=self.reference)
            return solution

    def _repr_html_(self):

        template = """
        <h4>Result:</h4>
        <ul>
            <li>model: %s</li>
            <li>heuristic: %s</li>
            <li>objective function: %s</li>
            <li>simulation method: %s</li>
            <li>type: %s</li>
        <ul>
        """
        model_id = self.model.id
        heuristic = self.heuristic_method.__class__.__name__
        of_string = "<br/>".join([o._repr_latex_() for o in self.objective_functions])
        simulation = self.simulation_method.__name__
        solutions = self.solutions._repr_html_()

        results = template % (model_id, heuristic, of_string, simulation, self.ko_type)
        return results + solutions

    def _merge(self, other_result):
        assert isinstance(other_result, self.__class__), "Cannot merge result with %s" % type(other_result)
        assert self.model.id == other_result.model.id, "Cannot merge results from different models"
        # assert self.objective_functions == other_result.objective_functions, \
        # "Cannot merge results with different objective functions"
        assert self.ko_type == other_result.ko_type, "Cannot merge results with resulting from different strategies"
        assert self.heuristic_method.__class__.__name__ == other_result.heuristic_method.__class__.__name__, \
            "Cannot merge results from different heuristic methods"

        self.solutions = self.solutions.append(other_result.solutions, ignore_index=True)
        self.solutions.drop_duplicates(subset=KNOCKOUTS, take_last=True, inplace=False)

        return self

    # TODO: find out how to plot an histogram (?) in bokeh
    def stats(self):
        stats_data = None
        if in_ipnb():
            if config.use_bokeh:
                stats_data = stats.BokehStatsData(self)
            elif config.use_matplotlib:
                pass
        else:
            stats_data = stats.CLIStatsData(self)

        stats_data.display()

    def visualize(self, index, map_name):
        if type == REACTION_KNOCKOUT_TYPE:
            knockouts = self.solutions[KNOCKOUTS][index]
        else:
            genes = [self.model.genes.get_by_id(g) for g in self.solutions[KNOCKOUTS][index]]
            knockouts = find_gene_knockout_reactions(self.model, genes)

        builder = draw_knockout_result(self.model, map_name, self.simulation_method, knockouts)
        return builder.display_in_notebook()

    def individuals(self):
        for index, row in self.solutions.iterrows():
            if len(self.objective_functions) == 1:
                fitness = row["Fitness 1"]
            else:
                fitness = Pareto(values=[row["Fitness %i"] % (i+1) for i in range(len(self.objective_functions))])
            if self.ko_type == GENE_KNOCKOUT_TYPE:
                for knockout in row["Knockouts"].args:
                    if isinstance(knockout, Symbol):
                        yield [str(knockout)], fitness
                    else:
                        yield [str(s) for s in knockout.args], fitness
            else:
                if isinstance(row["Knockouts"], Symbol):
                    yield [str(row["Knockouts"])], fitness
                else:
                    yield [str(s) for s in row["Knockouts"].args], fitness

    def plot(self, grid=None, width=None, height=None, title=None):
        pass

    def __iter__(self):
        for index, row in self.solutions.iterrows():
            if self.ko_type == GENE_KNOCKOUT_TYPE:
                yield StrainDesign(knockouts=row["Reactions"].args)
            else:
                yield StrainDesign(knockouts=row["Knockouts"].args)

    def data_frame(self):
        return DataFrame(self.solutions)


class ReactionKnockoutOptimization(KnockoutOptimization):
    """
    Knockout optimization using reactions.

    Attributes
    ----------
    model : SolverBasedModel
        A constraint-based model.
    heuristic_method : inspyred.ec.EvolutionaryComputation
        An evolutionary algorithm.
    objective_function : objective function or list(objective function)
        The objectives for the algorithm to maximize.
    seed : int
        A seed for random. It is auto-generated if None is given.
    termination : inspyred.ec.terminators
        A termination criteria for the algorithm. The default is inspyred.ec.terminators.evaluation_termination.
    simulation_method: flux_analysis.simulation
        The method used to simulate the model.
    wt_reference: dict
        A reference initial state for the optimization. It is required for flux_analysis.simulation.lmoma and
        flux_analysis.simulation.room. If not given, it will be computed using flux_analysis.simulation.pfba
    reactions: list
        A list of valid reactions to knockout. If None, then all reactions in the model will be knockout candidates
        except the ones defined in essential_reactions
    essential_reactions: list
        A list of reactions that cannot be knocked out. If None, then all essential reactions will be removed from
        the valid reactions set.

    Methods
    -------
    run(view=config.default_view, maximize=True, **kwargs)


    See Also
    --------
    *inspyred.ec
    *cameo.config.default_view

    Examples
    --------
    >>> from cameo import models
    >>> model = models.bigg.iJO1366
    >>> from cameo.strain_design.heuristic.objective_functions import biomass_product_coupled_yield
    >>> bpcy = biomass_product_coupled_yield(model.reactions.Ec_biomass_iJO1366_core_53p95,
    >>>                                      model.reactions.EX_succ_e),
    >>>                                      model.reactions.EX_glc__D_e)
    >>> knockout_optimization = ReactionKnockoutOptimization(model=model, objective_function=bpcy,
    >>>                                                      essential_reactions=["ATPM"])
    >>> knockout_optimization.run(max_evaluations=50000)


    """

    def __init__(self, reactions=None, essential_reactions=None, *args, **kwargs):
        super(ReactionKnockoutOptimization, self).__init__(*args, **kwargs)
        if reactions is None:
            self.reactions = set([r.id for r in self.model.reactions])
        else:
            self.reactions = reactions

        if essential_reactions is None:
            self.essential_reactions = set([r.id for r in self.model.essential_reactions()])
        else:
            self.essential_reactions = set([r.id for r in self.model.essential_reactions()] + essential_reactions)

        exchange_reactions = set([r.id for r in self.model.exchanges])
        self.representation = list(self.reactions.difference(self.essential_reactions).difference(exchange_reactions))
        self._ko_type = REACTION_KNOCKOUT_TYPE
        self._decoder = decoders.ReactionKnockoutDecoder(self.representation, self.model)


class GeneKnockoutOptimization(KnockoutOptimization):
    """
    Knockout optimization using genes.

    Attributes
    ----------
    model : SolverBasedModel
        A constraint-based model.
    heuristic_method : inspyred.ec.EvolutionaryComputation
        An evolutionary algorithm.
    objective_function : objective function or list(objective function)
        The objectives for the algorithm to maximize.
    seed : int
        A seed for random. It is auto-generated if None is given.
    termination : inspyred.ec.terminators
        A termination criteria for the algorithm. The default is inspyred.ec.terminators.evaluation_termination.
    simulation_method: flux_analysis.simulation
        The method used to simulate the model.
    wt_reference: dict
        A reference initial state for the optimization. It is required for flux_analysis.simulation.lmoma and
        flux_analysis.simulation.room. If not given, it will be computed using flux_analysis.simulation.pfba
    genes: list
        A list of valid genes to knockout. If None, then all genes in the model will be knockout candidates except the
         ones defined in essential_genes
    essential_genes: list
        A list of genes that cannot be knocked out. If None, then all essential genes will be removed from the valid
        genes set.

    Methods
    -------
    run(view=config.default_view, maximize=True, **kwargs)


    See Also
    --------
    *inspyred.ec
    *cameo.config.default_view

    Examples
    --------
    >>> from cameo import models
    >>> model = models.bigg.iJO1366
    >>> from cameo.strain_design.heuristic.objective_functions import biomass_product_coupled_yield
    >>> bpcy = biomass_product_coupled_yield(model.reactions.Ec_biomass_iJO1366_core_53p95,
    >>>                                      model.reactions.EX_succ_e),
    >>>                                      model.reactions.EX_glc__D_e)
    >>> knockout_optimization = GeneKnockoutOptimization(model=model, objective_function=bpcy)
    >>> knockout_optimization.run(max_evaluations=50000)

    """

    def __init__(self, genes=None, essential_genes=None, *args, **kwargs):
        super(GeneKnockoutOptimization, self).__init__(*args, **kwargs)
        if genes is None:
            self.genes = set([g.id for g in self.model.genes])
        else:
            self.genes = genes

        if essential_genes is None:
            self.essential_genes = set([g.id for g in self.model.essential_genes()])
        else:
            self.essential_genes = set([g.id for g in self.model.essential_genes()] + essential_genes)

        self.representation = list(self.genes.difference(self.essential_genes))
        self._ko_type = GENE_KNOCKOUT_TYPE
        self._decoder = decoders.GeneKnockoutDecoder(self.representation, self.model)


class KnockinKnockoutEvaluator(KnockoutEvaluator):
    def __call__(self, *args, **kwargs):
        pass


class KnockinKnockoutOptimizationResult:
    pass


class KnockinKnockoutOptimization(KnockoutOptimization):
    def __init__(self, *args, **kwargs):
        super(KnockinKnockoutOptimization, self).__init__(*args, **kwargs)

    def _evaluator(self, candidates, args):
        view = args.get('view')
        population_chunks = (chunk for chunk in partition(candidates, len(view)))
        func_obj = KnockinKnockoutEvaluator(self.model, self._decoder, self.objective_function, self.simulation_method)
        results = view.map(func_obj, population_chunks)
        fitness = reduce(list.__add__, results)

        return fitness

    def run(self, **kwargs):
        for observer in self.observers:
            observer.reset()
        self.heuristic_method.observer = self.observers
        super(KnockoutOptimization, self).run(
            keys=['knockout', 'knockin'],
            knockout_representation=self.representation,
            knockin_representation=self.knockin_representaion,
            candidate_size=self.max_size,
            variable_candidate_size=self.variable_size,
            **kwargs)
        return KnockinKnockoutOptimizationResult(model=self.model,
                                                 heuristic_method=self.heuristic_method,
                                                 simulation_method=self.simulation_method,
                                                 solutions=self.heuristic_method.archive,
                                                 objective_function=self.objective_function,
                                                 ko_type=self._ko_type,
                                                 decoder=self._decoder,
                                                 product=kwargs.get('product', None))


# TODO: implement a knockout knockin approach using a reaction db
class ReactionKnockinKnockoutOptimization(ReactionKnockoutOptimization, KnockinKnockoutOptimization):
    def __init__(self, reaction_db=None, *args, **kwargs):
        KnockinKnockoutOptimization.__init__(self, *args, **kwargs)
        ReactionKnockoutOptimization.__init__(self, *args, **kwargs)

    def run(self, **kwargs):
        KnockinKnockoutOptimization.run(self, **kwargs)
