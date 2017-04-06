# Copyright 2016 The Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import logging

from cobra.exceptions import OptimizationError

from cameo.core.manipulation import swap_cofactors
from cameo.strain_design.heuristic.evolutionary.decoders import SetDecoder
from cameo.strain_design.heuristic.evolutionary.objective_functions import ObjectiveFunction
from cameo.util import ProblemCache, memoize

logger = logging.getLogger(__name__)

__all__ = ['KnockoutEvaluator', 'SwapEvaluator']


class Evaluator(object):
    """
    Any evaluator that takes a population and returns a fitness.
    """
    def __call__(self, population):
        raise NotImplementedError

    def evaluate_individual(self, individual):
        raise NotImplementedError

    def reset(self):
        raise NotImplementedError


class TargetEvaluator(Evaluator):
    """
    Evaluator for targets in a model. It gets a representation, decodes into targets and simulates
    the model with a given method.

    Attributes
    ----------
    model : cobra.Model
        A constraint-based model
    decoder : SetDecoder
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
        if not isinstance(decoder, SetDecoder):
            raise ValueError("Invalid decoder %s" % decoder)
        self.decoder = decoder

        if not isinstance(objective_function, ObjectiveFunction):
            raise ValueError("'objective_function' must be instance of ObjectiveFunction (%s)" % objective_function)

        self.objective_function = objective_function
        self.simulation_method = simulation_method
        self.simulation_kwargs = simulation_kwargs
        self.cache = ProblemCache(model)

    def __call__(self, population):
        return [self.evaluate_individual(tuple(i)) for i in population]

    def reset(self):
        self.cache.reset()


class KnockoutEvaluator(TargetEvaluator):
    """
    Knockout evaluator for genes or reactions.
    """

    @memoize
    def evaluate_individual(self, individual):
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
        targets = self.decoder(individual)[0]
        with self.model:
            for target in targets:
                target.knock_out()
            try:
                solution = self.simulation_method(self.model,
                                                  cache=self.cache,
                                                  volatile=False,
                                                  raw=True,
                                                  reactions=self.objective_function.reactions,
                                                  **self.simulation_kwargs)
                fitness = self.objective_function(self.model, solution, targets)
            except OptimizationError as e:
                logger.debug(e)
                fitness = self.objective_function.worst_fitness()
            return fitness


class SwapEvaluator(TargetEvaluator):
    """ Evaluate reaction swaps where we knock one reaction in favor of another """

    def __init__(self, swap_pair=None, *args, **kwargs):
        super(SwapEvaluator, self).__init__(*args, **kwargs)
        self.swap_pair = swap_pair

    @memoize
    def evaluate_individual(self, individual):
        swap_reactions = self.decoder(individual)[0]
        with self.model:
            for reaction in swap_reactions:
                swap_cofactors(reaction, self.model, self.swap_pair, inplace=True)
            try:
                solution = self.simulation_method(self.model,
                                                  cache=self.cache,
                                                  volatile=False,
                                                  raw=True,
                                                  reactions=self.objective_function.reactions,
                                                  **self.simulation_kwargs)
                fitness = self.objective_function(self.model, solution, swap_reactions)
            except OptimizationError as e:
                logger.debug(e)
                fitness = self.objective_function.worst_fitness()

            return fitness


class KnockinKnockoutEvaluator(KnockoutEvaluator):
    def __init__(self, database, *args, **kwargs):
        super(KnockinKnockoutEvaluator, self).__init__(*args, **kwargs)
        self._database = database
