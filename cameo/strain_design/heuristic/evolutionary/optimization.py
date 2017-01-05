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
import types
from functools import reduce

import inspyred
import time
from pandas import DataFrame

from cameo import config
from cameo.core.result import Result
from cameo.core.solver_based_model import SolverBasedModel
from cameo.flux_analysis.simulation import pfba, lmoma, moma, room
from cameo.strain_design.heuristic.evolutionary import archives
from cameo.strain_design.heuristic.evolutionary import decoders
from cameo.strain_design.heuristic.evolutionary import evaluators
from cameo.strain_design.heuristic.evolutionary import generators
from cameo.strain_design.heuristic.evolutionary import observers
from cameo.strain_design.heuristic.evolutionary import plotters
from cameo.strain_design.heuristic.evolutionary import stats
from cameo.strain_design.heuristic.evolutionary import variators
from cameo.strain_design.heuristic.evolutionary.objective_functions import MultiObjectiveFunction, ObjectiveFunction
from cameo.util import RandomGenerator as Random
from cameo.util import in_ipnb
from cameo.util import partition


__all__ = ['ReactionKnockoutOptimization', 'GeneKnockoutOptimization', 'CofactorSwapOptimization']

REACTION_KNOCKOUT_TYPE = "reaction"
SWAP_TYPE = "cofactor-swap"
GENE_KNOCKOUT_TYPE = "gene"
NADH_NADPH = (['nad_c', 'nadh_c'], ['nadp_c', 'nadph_c'])

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
        self.observers = []
        self.model = model
        self._random = None
        self.termination = termination
        self._objective_function = objective_function
        self._heuristic_method = None
        self.heuristic_method = heuristic_method
        self.heuristic_method.terminator = termination

    @property
    def archiver(self):
        return self._heuristic_method.archiver

    @archiver.setter
    def archiver(self, archiver):
        self._heuristic_method.archiver = archiver

    @property
    def objective_function(self):
        return self._objective_function

    @objective_function.setter
    def objective_function(self, objective_function):
        if not isinstance(objective_function, ObjectiveFunction):
            raise TypeError("objective function is not instance of ObjectiveFunction")
        elif self._heuristic_method.__module__ == inspyred.ec.ec.__name__ and isinstance(objective_function,
                                                                                         MultiObjectiveFunction):
            raise TypeError("single objective heuristic do not support multiple objective functions")
        else:
            self._objective_function = objective_function

    @property
    def random(self):
        if self._random is None:
            self._random = Random()
        return self._random

    @property
    def heuristic_method(self):
        return self._heuristic_method

    @heuristic_method.setter
    def heuristic_method(self, heuristic_method):
        if heuristic_method.__module__ == inspyred.ec.ec.__name__ and isinstance(self.objective_function,
                                                                                 MultiObjectiveFunction):
            raise TypeError("single objective heuristics do not support multiple objective functions")
        self._heuristic_method = heuristic_method(self.random)

    def run(self, evaluator=None, generator=None, view=config.default_view, maximize=True, **kwargs):
        if kwargs.get('seed', None) is None:
            kwargs['seed'] = int(time.time())

        self._heuristic_method._random.seed(kwargs['seed'])

        for observer in self.observers:
            observer.reset()
        t = time.time()
        print(time.strftime("Starting optimization at %a, %d %b %Y %H:%M:%S", time.localtime(t)))
        res = self.heuristic_method.evolve(generator=generator,
                                           maximize=maximize,
                                           evaluator=evaluator,
                                           **kwargs)
        for observer in self.observers:
            observer.end()
        runtime = time.time() - t
        print(time.strftime("Finished after %H:%M:%S", time.gmtime(runtime)))

        return res


class EvaluatorWrapper(object):
    def __init__(self, view, evaluator):
        if not hasattr(view, 'map'):
            raise ValueError("View %s does not contain the required map function")
        if not (hasattr(evaluator, '__call__') or isinstance(evaluator, types.FunctionType)):
            raise ValueError("evaluator %s must be a function or callable")
        self.view = view
        self.evaluator = evaluator
        self.__name__ = "Wrapped %s" % EvaluatorWrapper.__class__.__name__

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.evaluator.reset()

    def __call__(self, candidates, args):
        population_chunks = (chunk for chunk in partition(candidates, len(self.view)))
        try:
            chunked_results = self.view.map(self.evaluator, population_chunks)
        except KeyboardInterrupt as e:
            self.view.shutdown()
            raise e

        fitness = reduce(list.__add__, chunked_results)

        return fitness


class TargetOptimization(HeuristicOptimization):
    """
    Abstract class for target optimization.
    """

    def __init__(self, simulation_method=pfba, wt_reference=None, *args, **kwargs):
        """
        Class for generic optimization algorithms for knockout (or similar) strain design methods

        Attributes
        ----------
        simulation_method: see flux_analysis.simulation
        wt_reference: dict
        simulation_method: method
           the simulation method to use for evaluating results
        evaluator: TargetEvaluator
           the class used to evaluate results
        """
        super(TargetOptimization, self).__init__(*args, **kwargs)
        self._simulation_kwargs = dict()
        self._simulation_kwargs['reference'] = wt_reference
        self._simulation_method = None
        self.simulation_method = simulation_method
        self.representation = None
        self._evaluator = None
        self._target_type = None
        self._decoder = None
        self._metadata = {}

    @property
    def metadata(self):
        return self._metadata

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
                if len(self.objective_function) > 1:
                    self.observers.append(plotters.IPythonBokehParetoPlotter(self.objective_function))
                else:
                    self.observers.append(plotters.IPythonBokehFitnessPlotter())

        else:
            if config.use_bokeh:
                pass
            else:
                pass
        if self.progress:
            self.observers.append(observers.ProgressObserver())

    def run(self, view=config.default_view, max_size=10, variable_size=True, **kwargs):
        """
        Parameters
        ----------
        max_size: int
            Maximum size of a solution, e.g., the maximum number of reactions or genes to knock-out or swap
        variable_size: boolean
            If true, the solution size can change meaning that the combination of knockouts can have different sizes up to
            max_size. Otherwise it only produces knockout solutions with a fixed number of knockouts.

        """

        if kwargs.get('seed', None) is None:
            kwargs['seed'] = int(time.time())

        self.heuristic_method.observer = self.observers

        with EvaluatorWrapper(view, self._evaluator) as evaluator:
            super(TargetOptimization, self).run(distance_function=set_distance_function,
                                                representation=self.representation,
                                                evaluator=evaluator,
                                                generator=generators.set_generator,
                                                max_size=max_size,
                                                **kwargs)

            return TargetOptimizationResult(model=self.model,
                                            heuristic_method=self.heuristic_method,
                                            simulation_method=self.simulation_method,
                                            simulation_kwargs=self._simulation_kwargs,
                                            solutions=self.heuristic_method.archive,
                                            objective_function=self.objective_function,
                                            target_type=self._target_type,
                                            decoder=self._decoder,
                                            seed=kwargs['seed'],
                                            metadata=self.metadata)


class KnockoutOptimization(TargetOptimization):
    """
    Abstract knockout optimization class.
    """
    def __init__(self, simulation_method=pfba, wt_reference=None, *args, **kwargs):
        super(KnockoutOptimization, self).__init__(simulation_method=simulation_method,
                                                   wt_reference=wt_reference,
                                                   *args, **kwargs)


class TargetOptimizationResult(Result):
    def __init__(self, model=None, heuristic_method=None, simulation_method=None,
                 simulation_kwargs=None, solutions=None, objective_function=None,
                 target_type=None, decoder=None, seed=None, metadata=None, *args, **kwargs):
        super(TargetOptimizationResult, self).__init__(*args, **kwargs)
        self.seed = seed
        self.model = model
        self.heuristic_method = heuristic_method
        self.simulation_method = simulation_method
        self.simulation_kwargs = simulation_kwargs or {}
        if isinstance(objective_function, list):
            self.objective_functions = objective_function
        else:
            self.objective_functions = [objective_function]
        self.target_type = target_type
        self._decoder = decoder
        self._metadata = metadata
        self._solutions = self._decode_solutions(solutions)

    def __len__(self):
        return len(self._solutions)

    def __getstate__(self):
        state = super(TargetOptimizationResult, self).__getstate__()
        state.update({'model': self.model,
                      'decoder': self._decoder,
                      'simulation_method': self.simulation_method,
                      'simulation_kwargs': self.simulation_kwargs,
                      'heuristic_method.__class__': self.heuristic_method.__class__,
                      'heuristic_method.maximize': self.heuristic_method.maximize,
                      'heuristic_method.variator': self.heuristic_method.variator,
                      'heuristic_method.terminator': self.heuristic_method.terminator,
                      'heuristic_method.archiver': self.heuristic_method.archiver,
                      'heuristic_method.archive': self.heuristic_method.archive,
                      'heuristic_method.termination_cause': self.heuristic_method.termination_cause,
                      'heuristic_method._random': self.heuristic_method._random,
                      'heuristic_method.generator': self.heuristic_method.generator,
                      'heuristic_method._kwargs': self.heuristic_method._kwargs,
                      'objective_functions': self.objective_functions,
                      'target_type': self.target_type,
                      'solutions': self._solutions,
                      'seed': self.seed,
                      'metadata': self._metadata})

    def __setstate__(self, state):
        super(TargetOptimizationResult, self).__setstate__(state)
        self.model = state['model']
        self.simulation_method = state['simulation_method']
        self.simulation_kwargs = state['simulation_kwargs']
        self.seed = state['seed']
        random = state['heuristic_method._random']
        self.heuristic_method = state['heuristic_method.__class__'](random)
        self.heuristic_method.maximize = state['heuristic_method.maximize']
        self.heuristic_method.terminator = state['heuristic_method.terminator']
        self.heuristic_method.termination_cause = state['heuristic_method.termination_cause']
        self.heuristic_method.archiver = state['heuristic_method.archiver']
        self.heuristic_method.archive = state['heuristic_method.archive']
        self.heuristic_method._kwargs = state['heuristic_method._kwargs']
        self.objective_functions = state['objective_functions']
        self.target_type = state['target_type']
        self._solutions = state['solutions']
        self._metadata = state['metadata']

    def _repr_html_(self):

        template = """
        <h4>Result:</h4>
        <ul>
            <li>model: %s</li>
            <li>heuristic: %s</li>
            <li>objective function: %s</li>
            <li>simulation method: %s</li>
            <li>target type: %s</li>
        <ul>
        """

        model_id = self.model.id
        heuristic = self.heuristic_method.__class__.__name__
        of_string = "<br/>".join([o._repr_latex_() for o in self.objective_functions])
        simulation = self.simulation_method.__name__
        solutions = self.data_frame._repr_html_()

        results = template % (model_id, heuristic, of_string, simulation, self.target_type)
        return results + solutions

    def __iter__(self):
        for _, row in self._solutions.iterrows():
            yield [row['targets'], row['fitness']]

    def __iadd__(self, other):
        if not isinstance(other, self.__class__):
            raise AssertionError("Cannot merge result with %s" % type(other))
        if self.model.id != other.model.id:
            raise AssertionError("Cannot merge results from different models")
        if self.target_type != other.target_type:
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

    @property
    def data_frame(self):
        return DataFrame(self._solutions)

    def _decode_solutions(self, solutions):
        decoded_solutions = DataFrame(columns=["targets", "fitness"])
        for index, solution in enumerate(solutions):
            knockouts = self._decoder(solution.candidate, flat=True)
            if len(knockouts) > 0:
                decoded_solutions.loc[index] = [knockouts, solution.fitness]

        return decoded_solutions


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
        self._target_type = REACTION_KNOCKOUT_TYPE
        self._decoder = decoders.ReactionSetDecoder(self.representation, self.model)
        self._evaluator = evaluators.KnockoutEvaluator(model=self.model,
                                                       decoder=self._decoder,
                                                       objective_function=self.objective_function,
                                                       simulation_method=self._simulation_method,
                                                       simulation_kwargs=self._simulation_kwargs)


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
        self._target_type = GENE_KNOCKOUT_TYPE
        self._decoder = decoders.GeneSetDecoder(self.representation, self.model)
        self._evaluator = evaluators.KnockoutEvaluator(model=self.model,
                                                       decoder=self._decoder,
                                                       objective_function=self.objective_function,
                                                       simulation_method=self._simulation_method,
                                                       simulation_kwargs=self._simulation_kwargs)


class CofactorSwapOptimization(TargetOptimization):
    """
    Optimize co-factor swapping

    As suggested in [1]_, flux through a given reaction can sometimes be optimized by swapping complementary
    co-factor. This class implements a search for reactions when swapped improve the given objective. Briefly,
    the approach is to

    - find reactions that have all the targeted co-factor pairs e.g. (nad_c -> nadp_c, nadh_c -> nadph_c)

    - add reactions that have the co-factors swapped and then by a search algorithm switching one off in favor of the
      other

    The implementation here differs from that in [1]_ in that we use a general purpose search algorithm rather than
    formulating the search as a mixed integer linear programming problem.

    References
    ----------
    .. [1] King, Zachary A., and Adam M. Feist. "Optimizing Cofactor Specificity of Oxidoreductase Enzymes for the
      Generation of Microbial Production Strains - OptSwap." Industrial Biotechnology 9, no. 4 (August 1,
      2013): 236-46. - doi:10.1089/ind.2013.0005.

    Parameters
    ----------
    model : SolverBasedModel
       the model to operator on
    cofactor_id_swaps : tuple
       a tuple of length 2 that defines two lists of metabolite identifiers that should be interchanged during the
       swap optimization see e.g. `NADH_NADPH` which is also the default.
    candidate_reactions : list
       reactions to consider for co-factor swap - if not given then search for all reactions that include the given
       cofactors
    skip_reactions : list
       reactions to not consider for co-factor swap, defaults to the objective function if not provided
    args, kwargs : keyword arguments
       passed on to super-classes, see in particular `objective_function`, `heuristic_method`, `termination`,
       `simulation_method`, `wt_reference`, of `HeuristicOptimization` and `max_size` of `HeuristicOptimization.run`

    Examples
    --------
    >>> from cameo import models
    >>> from cameo.strain_design.heuristic.evolutionary.objective_functions import product_yield
    >>> model = models.bigg.iJO1366
    >>> model.objective = model.reactions.EX_thr__L_e
    >>> model.reactions.BIOMASS_Ec_iJO1366_core_53p95M.lower_bound = 0.1
    >>> py = product_yield(model.reactions.EX_thr__L_e, model.reactions.EX_glc__D_e)
    >>> swap_optimization = CofactorSwapOptimization(model=model, objective_function=py)
    >>> swap_optimization.run(max_evaluations=2000, max_size=2)
    """

    def __init__(self, cofactor_id_swaps=NADH_NADPH, candidate_reactions=None, skip_reactions=None, *args, **kwargs):
        super(CofactorSwapOptimization, self).__init__(*args, **kwargs)
        self._target_type = SWAP_TYPE
        swap_pairs = ([self.model.metabolites.get_by_id(m) for m in cofactor_id_swaps[0]],
                      [self.model.metabolites.get_by_id(m) for m in cofactor_id_swaps[1]])
        self.metadata['swap_pairs'] = swap_pairs
        self.representation = candidate_reactions or self.find_swappable_reactions(self.model, swap_pairs)
        if skip_reactions:
            self.representation -= skip_reactions
        self._decoder = decoders.ReactionSetDecoder(self.representation, self.model)

        self._evaluator = evaluators.SwapEvaluator(model=self.model,
                                                   decoder=self._decoder,
                                                   objective_function=self.objective_function,
                                                   simulation_method=self._simulation_method,
                                                   simulation_kwargs=self._simulation_kwargs,
                                                   swap_pair=swap_pairs)

    @staticmethod
    def find_swappable_reactions(model, swaps):
        """
        Get all reactions that can undergo co-factor swapping

        Find reactions that have one set of the cofactors targeted for swapping and are mass balanced and updates the
        `candidate_reactions` attribute

        Arguments
        ---------
        model: SolverBasedModel
            A model with reactions to search on.
        swaps: tuple
            Pair of cofactors to swap.
        """

        def swap_search(mets):
            has_pairs = all(mets.get(m, False) for m in swaps[0]) or all(mets.get(m, False) for m in swaps[1])

            contains_all = all(mets.get(m, False) for m in swaps[0]) and all(mets.get(m, False) for m in swaps[1])

            return has_pairs and not contains_all

        candidates = model.reactions.query(search_function=swap_search, attribute='metabolites')
        return [r.id for r in candidates]
