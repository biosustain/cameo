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

import logging
import time
from functools import reduce

import inspyred
from pandas import DataFrame

from cameo import config
from cameo.core.result import Result
from cameo.exceptions import SolveError
from cameo.flux_analysis.simulation import pfba, lmoma, moma, room
from cameo.strain_design.heuristic.evolutionary import archives
from cameo.strain_design.heuristic.evolutionary import decoders
from cameo.strain_design.heuristic.evolutionary import generators
from cameo.strain_design.heuristic.evolutionary import observers
from cameo.strain_design.heuristic.evolutionary import plotters
from cameo.strain_design.heuristic.evolutionary import variators
from cameo.strain_design.heuristic.evolutionary import stats
from cameo.strain_design.heuristic.evolutionary.processing import reactions2filter
from cameo.util import RandomGenerator as Random
from cameo.util import in_ipnb
from cameo.util import partition, TimeMachine, memoize, ProblemCache

__all__ = ['ReactionKnockoutOptimization', 'GeneKnockoutOptimization']

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
        archives.BestSolutionArchive(),
    ],
    inspyred.ec.SA: [

        [
            variators.set_mutation,
            variators.set_indel
        ],
        inspyred.ec.selectors.default_selection,
        inspyred.ec.replacers.simulated_annealing_replacement,
        archives.BestSolutionArchive()
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


class HeuristicOptimization(object):
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
        self._seed = seed
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

    @property
    def seed(self):
        return self._seed

    @seed.setter
    def seed(self, seed):
        if seed is None:
            seed = int(time.time())
        self._seed = seed
        self.random = Random(seed=seed)

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
        try:
            res = [self._evaluate_individual(frozenset(i)) for i in population]
        except Exception as e:
            logger.error(e)
            res = None
        finally:
            self.cache.reset()

        return res

    @memoize
    def _evaluate_individual(self, individual):
        """
        Evaluates a single individual.

        Arguments
        ---------

        individual: set
            The encoded representation of a single individual.


        Returns
        -------
        fitness
            A single real value or a Pareto, depending on the number of objectives.
        """
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
                                                  reactions=reactions2filter(self.objective_function),
                                                  **self.simulation_kwargs)
                fitness = self._calculate_fitness(solution, decoded)
            except SolveError as e:
                logger.debug(e)
                if isinstance(self.objective_function, list):
                    fitness = inspyred.ec.emo.Pareto(values=[0 for _ in self.objective_function])
                else:
                    fitness = 0

            return fitness

    def _calculate_fitness(self, solution, decoded):
        if isinstance(self.objective_function, list):
            logger.debug("evaluate multiobjective solution")
            return inspyred.ec.emo.Pareto(values=[of(self.model, solution, decoded) for of in self.objective_function])
        else:
            logger.debug("evaluate single objective solution")
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
        self._simulation_kwargs = dict()
        self._simulation_kwargs['reference'] = wt_reference
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
        if simulation_method in [lmoma, moma, room] and self._simulation_kwargs.get("reference", None) is None:
            logger.warning("No WT reference found, generating using pfba.")
            self._simulation_kwargs['reference'] = pfba(self.model).fluxes
            logger.warning("Reference successfully computed.")
        self._simulation_method = simulation_method

    @property
    def simulation_kwargs(self):
        return self._simulation_kwargs

    @simulation_kwargs.setter
    def simulation_kwargs(self, simulation_kwargs):
        if self.simulation_method in [lmoma, moma, room] and simulation_kwargs.get("reference", None) is None:
            logger.warning("No WT reference found, generating using pfba.")
            simulation_kwargs['reference'] = pfba(self.model).fluxes
            logger.warning("Reference successfully computed.")
        self._simulation_kwargs = simulation_kwargs

    def _evaluator(self, candidates, args):
        view = args.get('view')
        population_chunks = (chunk for chunk in partition(candidates, len(view)))
        func_obj = KnockoutEvaluator(self.model, self._decoder, self.objective_function,
                                     self.simulation_method, self._simulation_kwargs)
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

    def run(self, max_size=10, variable_size=True, **kwargs):
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
                                          simulation_kwargs=self._simulation_kwargs,
                                          solutions=self.heuristic_method.archive,
                                          objective_function=self.objective_function,
                                          ko_type=self._ko_type,
                                          decoder=self._decoder,
                                          seed=self.seed)


class KnockoutOptimizationResult(Result):
    def __init__(self, model=None, heuristic_method=None, simulation_method=None, simulation_kwargs=None,
                 solutions=None, objective_function=None, ko_type=None, decoder=None, seed=None, *args, **kwargs):
        super(KnockoutOptimizationResult, self).__init__(*args, **kwargs)
        self.seed = seed
        self.model = model
        self.heuristic_method = heuristic_method
        self.simulation_method = simulation_method
        self.simulation_kwargs = simulation_kwargs or {}
        if isinstance(objective_function, list):
            self.objective_functions = objective_function
        else:
            self.objective_functions = [objective_function]
        self.ko_type = ko_type
        self._decoder = decoder
        self._solutions = self._decode_solutions(solutions)

    def _decode_solutions(self, solutions):
        decoded_solutions = DataFrame(columns=["reactions", "knockouts", "fitness"])
        for index, solution in enumerate(solutions):
            decoded_solutions.loc[index] = self._decoder(solution, flat=True) + [solution.fitness]
        return decoded_solutions

    def __len__(self):
        return len(self._solutions)

    def __getstate__(self):
        return {'model': self.model,
                'decoder': self._decoder,
                'simulation_method': self.simulation_method,
                'simulation_kwargs': self.simulation_kwargs,
                'heuristic_method.__class__': self.heuristic_method.__class__,
                'heuristic_method.maximize': self.heuristic_method.maximize,
                'heuristic_method.variator': self.heuristic_method.variator,
                'heuristic_method.terminator': self.heuristic_method.terminator,
                'heuristic_method.archiver': self.heuristic_method.archiver,
                'heuristic_method.termination_cause': self.heuristic_method.termination_cause,
                'heuristic_method._random': self.heuristic_method._random,
                'heuristic_method.generator': self.heuristic_method.generator,
                'heuristic_method._kwargs': self.heuristic_method._kwargs,
                'objective_functions': self.objective_functions,
                'ko_type': self.ko_type,
                'solutions': self._solutions,
                'seed': self.seed}

    def __setstate__(self, d):
        self.model = d['model']
        self.simulation_method = d['simulation_method']
        self.simulation_kwargs = d['simulation_kwargs']
        self.seed = d['seed']
        random = d['heuristic_method._random']
        self.heuristic_method = d['heuristic_method.__class__'](random)
        self.heuristic_method.maximize = d['heuristic_method.maximize']
        self.heuristic_method.terminator = d['heuristic_method.terminator']
        self.heuristic_method.termination_cause = d['heuristic_method.termination_cause']
        self.heuristic_method.archive = d['heuristic_method.archive']
        self.heuristic_method._kwargs = d['heuristic_method._kwargs']
        self.objective_functions = d['objective_functions']
        self.ko_type = d['ko_type']
        self._solutions = d['solutions']

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
        solutions = self.data_frame._repr_html_()

        results = template % (model_id, heuristic, of_string, simulation, self.ko_type)
        return results + solutions

    def __iter__(self):
        for _, row in self._solutions.iterrows():
            yield [row['reactions'], row['knockouts'], row['fitness']]

    def __iadd__(self, other):
        if not isinstance(other, self.__class__):
            raise AssertionError("Cannot merge result with %s" % type(other))
        if self.model.id != other.model.id:
            raise AssertionError("Cannot merge results from different models")
        if self.ko_type != other.ko_type:
            raise AssertionError("Cannot merge results with resulting from different strategies")
        if self.heuristic_method.__class__.__name__ != other.heuristic_method.__class__.__name__:
            raise AssertionError("Cannot merge results from different heuristic methods")

        self._solutions = self._solutions.append(other._solutions, ignore_index=True)
        self._solutions.drop_duplicates(subset="knockouts", take_last=True, inplace=True)

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

    def plot(self, grid=None, width=None, height=None, title=None):
        pass

    @property
    def data_frame(self):
        return DataFrame(self._solutions)


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
    >>> from cameo.strain_design.heuristic.evolutionary.objective_functions import biomass_product_coupled_yield
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
        logger.debug("Computing essential reactions...")
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
    >>> from cameo.strain_design.heuristic.evolutionary.objective_functions import biomass_product_coupled_yield
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
        logger.debug("Computing essential genes...")
        if essential_genes is None:
            self.essential_genes = set([g.id for g in self.model.essential_genes()])
        else:
            self.essential_genes = set([g.id for g in self.model.essential_genes()] + essential_genes)

        self.representation = list(self.genes.difference(self.essential_genes))
        self._ko_type = GENE_KNOCKOUT_TYPE
        self._decoder = decoders.GeneKnockoutDecoder(self.representation, self.model)


class KnockinKnockoutEvaluator(KnockoutEvaluator):
    def __init__(self, database, *args, **kwargs):
        super(KnockinKnockoutEvaluator, self).__init__(*args, **kwargs)
        self._database = database

    def __call__(self, *args, **kwargs):
        pass


class KnockinKnockoutOptimizationResult:
    pass


class KnockinKnockoutOptimization(KnockoutOptimization):
    def __init__(self, database=None, *args, **kwargs):
        super(KnockinKnockoutOptimization, self).__init__(*args, **kwargs)
        self._database = database

    def _evaluator(self, candidates, args):
        view = args.get('view')
        population_chunks = (chunk for chunk in partition(candidates, len(view)))
        func_obj = KnockinKnockoutEvaluator(self.model, self._decoder,
                                            self.objective_function, self.simulation_method,
                                            self.simulation_kwargs)
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
            knockin_representation=self._database,
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
