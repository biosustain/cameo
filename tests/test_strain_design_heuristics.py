# Copyright 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import absolute_import, print_function

import os
import pickle
import time
from collections import namedtuple
from math import sqrt
from tempfile import mkstemp

import inspyred
import numpy
import pytest
import six
from inspyred.ec import Bounder
from inspyred.ec.emo import Pareto
from ordered_set import OrderedSet
from six.moves import range

from cameo import config, fba
from cameo.core.manipulation import swap_cofactors
from cameo.parallel import SequentialView
from cameo.strain_design import OptGene
from cameo.strain_design.heuristic.evolutionary.archives import (BestSolutionArchive,
                                                                 Individual)
from cameo.strain_design.heuristic.evolutionary.decoders import (GeneSetDecoder,
                                                                 ReactionSetDecoder,
                                                                 SetDecoder)
from cameo.strain_design.heuristic.evolutionary.evaluators import KnockoutEvaluator
from cameo.strain_design.heuristic.evolutionary.generators import (linear_set_generator,
                                                                   multiple_chromosome_set_generator,
                                                                   set_generator)
from cameo.strain_design.heuristic.evolutionary.genomes import MultipleChromosomeGenome
from cameo.strain_design.heuristic.evolutionary.metrics import (euclidean_distance,
                                                                manhattan_distance)
from cameo.strain_design.heuristic.evolutionary.multiprocess.migrators import MultiprocessingMigrator
from cameo.strain_design.heuristic.evolutionary.objective_functions import (MultiObjectiveFunction,
                                                                            YieldFunction,
                                                                            biomass_product_coupled_min_yield,
                                                                            biomass_product_coupled_yield,
                                                                            number_of_knockouts,
                                                                            product_yield)
from cameo.strain_design.heuristic.evolutionary.optimization import (NADH_NADPH,
                                                                     CofactorSwapOptimization,
                                                                     EvaluatorWrapper,
                                                                     HeuristicOptimization,
                                                                     ReactionKnockoutOptimization,
                                                                     SolutionSimplification,
                                                                     TargetOptimizationResult,
                                                                     set_distance_function,
                                                                     GeneKnockoutOptimization)
from cameo.strain_design.heuristic.evolutionary.variators import (_do_set_n_point_crossover,
                                                                  multiple_chromosome_set_indel,
                                                                  multiple_chromosome_set_mutation,
                                                                  set_indel,
                                                                  set_mutation,
                                                                  set_n_point_crossover)
from cobra.flux_analysis import find_essential_genes, find_essential_reactions
from cameo.util import RandomGenerator as Random


try:
    from cameo.parallel import RedisQueue
except ImportError:
    RedisQueue = None

TRAVIS = bool(os.getenv('TRAVIS', False))

if os.getenv('REDIS_PORT_6379_TCP_ADDR'):
    REDIS_HOST = os.getenv('REDIS_PORT_6379_TCP_ADDR')  # wercker
else:
    REDIS_HOST = 'localhost'

SEED = 1234
CURRENT_PATH = os.path.dirname(__file__)

SOLUTIONS = [
    [[1, 2, 3], 0.1],
    [[1, 3, 2, 4], 0.1],
    [[2, 3, 4], 0.45],
    [[62, 51, 4], 0.2],
    [[5, 3, 4, 51], 0.9],
    [[5, 23, 41, 51], 0.9],
    [[5, 3, 4, 51, 31], 0.9],
    [[5, 3, 4, 51], 0.9],
    [[44, 12, 42, 51], 0.0],
    [[52, 22, 4, 11], 0.0]
]


@pytest.fixture(scope='function')
def generators(model):
    mockup_evolutionary_algorithm = namedtuple("EA", ["bounder"])
    args = {}
    args.setdefault('representation', [r.id for r in model.reactions])
    random = Random()
    return args, random, mockup_evolutionary_algorithm


@pytest.fixture(scope="module")
def objectives():
    single_objective_function = product_yield('product', 'substrate')
    multi_objective_function = MultiObjectiveFunction([
        product_yield('product', 'substrate'),
        number_of_knockouts()
    ])
    return single_objective_function, multi_objective_function


class TestMetrics:
    def test_euclidean_distance(self):
        distance = euclidean_distance({'a': 9}, {'a': 3})
        assert distance == sqrt((9 - 3) ** 2)

    def test_manhattan_distance(self):
        distance = manhattan_distance({'a': 9}, {'a': 3})
        assert distance == abs(9 - 3)


class TestBestSolutionArchive:
    def test_solution_string(self):
        sol1 = Individual(SOLUTIONS[0][0], SOLUTIONS[0][1])
        sol2 = Individual(SOLUTIONS[1][0], SOLUTIONS[1][1])
        sol3 = Individual(SOLUTIONS[2][0], SOLUTIONS[2][1])
        assert sol1.__str__() == "[1, 2, 3] - 0.1 sense: max"
        assert sol2.__str__() == "[1, 2, 3, 4] - 0.1 sense: max"
        assert sol3.__str__() == "[2, 3, 4] - 0.45 sense: max"

    def test_solution_comparison_maximization(self):
        sol1 = Individual(SOLUTIONS[0][0], SOLUTIONS[0][1])
        sol2 = Individual(SOLUTIONS[1][0], SOLUTIONS[1][1])
        sol3 = Individual(SOLUTIONS[2][0], SOLUTIONS[2][1])

        # test ordering
        assert sol1.__cmp__(sol2) == -1
        assert sol1.__cmp__(sol1) == 0
        assert sol1.__cmp__(sol3) == 1

        assert sol1 < sol2
        assert sol1 == sol1
        assert sol1 > sol3

        # test gt and lt
        assert sol1.__lt__(sol2)
        assert sol1.__gt__(sol3)
        assert not sol1.__lt__(sol1)
        assert not sol1.__gt__(sol1)
        assert not sol2.__lt__(sol1)
        assert not sol3.__gt__(sol1)

        # testing issubset
        assert sol1.issubset(sol2), "Solution 1 is subset of Solution 2"
        assert not sol2.issubset(sol1), "Solution 2 is not subset of Solution 1"
        assert sol3.issubset(sol2), "Solution 3 is subset of Solution 2"
        assert not sol2.issubset(sol3), "Solution 2 is not subset of Solution 3"
        assert not sol1.issubset(sol3), "Solution 1 is subset of Solution 3"
        assert not sol2.issubset(sol3), "Solution 3 is not subset of Solution 1"

        # test difference
        l = len(sol2.symmetric_difference(sol1))
        assert l == 1, "Difference between Solution 2 and 1 is (%s)" % sol2.symmetric_difference(sol1)
        l = len(sol3.symmetric_difference(sol2))
        assert l == 1, "Difference between Solution 3 and 1 is (%s)" % sol3.symmetric_difference(sol2)
        l = len(sol3.symmetric_difference(sol1))
        assert l == 2, "Difference between Solution 1 and 3 is (%s)" % sol3.symmetric_difference(sol1)

        assert sol1.improves(sol2), "Solution 1 is better than Solution 2"
        assert sol3.improves(sol2), "Solution 3 is better than Solution 2"
        assert not sol3.improves(sol1), "Solution 3 does not improve Solution 1"
        assert not sol2.improves(sol1), "Solution 2 does not improve Solution 1"
        assert not sol2.improves(sol3), "Solution 2 does not improve Solution 3"

    def test_solution_comparison_minimization(self):
        sol1 = Individual(SOLUTIONS[0][0], SOLUTIONS[0][1], maximize=False)
        sol2 = Individual(SOLUTIONS[1][0], SOLUTIONS[1][1], maximize=False)
        sol3 = Individual(SOLUTIONS[2][0], SOLUTIONS[2][1], maximize=False)

        # test ordering
        assert sol1.__cmp__(sol2) == -1
        assert sol1.__cmp__(sol1) == 0
        assert sol1.__cmp__(sol3) == -1
        assert sol3.__cmp__(sol1) == 1

        assert sol1 < sol2
        assert sol1 == sol1
        assert sol1 < sol3

        # test gt and lt
        assert sol1.__lt__(sol2)
        assert sol1.__lt__(sol3)
        assert not sol1.__gt__(sol1)
        assert not sol1.__lt__(sol1)
        assert sol2.__gt__(sol1)
        assert not sol3.__lt__(sol1)

        # testing issubset
        assert sol1.issubset(sol2), "Solution 1 is subset of Solution 2"
        assert not sol2.issubset(sol1), "Solution 2 is not subset of Solution 1"
        assert sol3.issubset(sol2), "Solution 3 is subset of Solution 2"
        assert not sol2.issubset(sol3), "Solution 2 is not subset of Solution 3"
        assert not sol1.issubset(sol3), "Solution 1 is subset of Solution 3"
        assert not sol2.issubset(sol3), "Solution 3 is not subset of Solution 1"

        # test difference
        l = len(sol2.symmetric_difference(sol1))
        assert l == 1, "Difference between Solution 2 and 1 is (%s)" % sol2.symmetric_difference(sol1)
        l = len(sol3.symmetric_difference(sol2))
        assert l == 1, "Difference between Solution 3 and 1 is (%s)" % sol3.symmetric_difference(sol2)
        l = len(sol3.symmetric_difference(sol1))
        assert l == 2, "Difference between Solution 1 and 3 is (%s)" % sol3.symmetric_difference(sol1)

        assert sol1.improves(sol2), "Solution 1 is better than Solution 2"
        assert not sol3.improves(sol2), "Solution 3 is not better than Solution 2"
        assert not sol3.improves(sol1), "Solution 3 does not improve Solution 1"
        assert not sol2.improves(sol1), "Solution 2 does not improve Solution 1"
        assert not sol2.improves(sol3), "Solution 2 does not improve Solution 3"

    def test_add_greater_solution_with_same_fitness(self):
        size = 1
        pool = BestSolutionArchive()
        pool.add(SOLUTIONS[0][0], SOLUTIONS[0][1], None, True, size)
        pool.add(SOLUTIONS[1][0], SOLUTIONS[1][1], None, True, size)
        assert pool.length() == 1, "Pool must keep one solution (length=%s)" % pool.length()
        best_solution = set(SOLUTIONS[0][0])
        best_fitness = SOLUTIONS[0][1]
        sol = pool.get(0)
        assert sol.candidate == best_solution, "Best solution set must be the first"
        assert sol.fitness == best_fitness, "Best solution fitness must be the first"

    def test_add_smaller_solution_with_same_fitness(self):
        size = 1
        pool = BestSolutionArchive()
        pool.add(SOLUTIONS[1][0], SOLUTIONS[1][1], None, True, size)
        pool.add(SOLUTIONS[0][0], SOLUTIONS[0][1], None, True, size)
        assert pool.length() == 1, "Pool must keep one solution (length=%s)" % pool.length()
        solution = set(SOLUTIONS[0][0])
        fitness = SOLUTIONS[0][1]
        sol = pool.get(0)
        assert sol.candidate == solution, "Best solution must be the first (%s)" % sol.candidate
        assert sol.fitness == fitness, "Best fitness must be the first (%s)" % sol.fitness

    def test_uniqueness_of_solutions(self):
        size = 2
        pool = BestSolutionArchive()
        pool.add(SOLUTIONS[1][0], SOLUTIONS[1][1], None, True, size)
        pool.add(SOLUTIONS[1][0], SOLUTIONS[1][1], None, True, size)

        assert pool.length() == 1, "Added repeated solution"

    def test_pool_size_limit(self):
        size = 1
        pool = BestSolutionArchive()
        pool.add(SOLUTIONS[0][0], SOLUTIONS[0][1], None, True, size)
        pool.add(SOLUTIONS[1][0], SOLUTIONS[1][1], None, True, size)
        pool.add(SOLUTIONS[2][0], SOLUTIONS[2][1], None, True, size)
        pool.add(SOLUTIONS[3][0], SOLUTIONS[3][1], None, True, size)
        pool.add(SOLUTIONS[4][0], SOLUTIONS[4][1], None, True, size)
        pool.add(SOLUTIONS[5][0], SOLUTIONS[5][1], None, True, size)
        pool.add(SOLUTIONS[6][0], SOLUTIONS[6][1], None, True, size)
        pool.add(SOLUTIONS[7][0], SOLUTIONS[7][1], None, True, size)
        pool.add(SOLUTIONS[8][0], SOLUTIONS[8][1], None, True, size)
        pool.add(SOLUTIONS[9][0], SOLUTIONS[9][1], None, True, size)
        assert pool.length() <= 1, "Pool must keep one solution (length=%s)" % pool.length()
        size = 2
        pool = BestSolutionArchive()
        pool.add(SOLUTIONS[0][0], SOLUTIONS[0][1], None, True, size)
        pool.add(SOLUTIONS[1][0], SOLUTIONS[1][1], None, True, size)
        pool.add(SOLUTIONS[2][0], SOLUTIONS[2][1], None, True, size)
        pool.add(SOLUTIONS[3][0], SOLUTIONS[3][1], None, True, size)
        pool.add(SOLUTIONS[4][0], SOLUTIONS[4][1], None, True, size)
        pool.add(SOLUTIONS[5][0], SOLUTIONS[5][1], None, True, size)
        pool.add(SOLUTIONS[6][0], SOLUTIONS[6][1], None, True, size)
        pool.add(SOLUTIONS[7][0], SOLUTIONS[7][1], None, True, size)
        pool.add(SOLUTIONS[8][0], SOLUTIONS[8][1], None, True, size)
        pool.add(SOLUTIONS[9][0], SOLUTIONS[9][1], None, True, size)
        assert pool.length() <= 2, "Pool must keep one solution (length=%s)" % pool.length()
        size = 3
        pool = BestSolutionArchive()
        pool.add(SOLUTIONS[0][0], SOLUTIONS[0][1], None, True, size)
        pool.add(SOLUTIONS[1][0], SOLUTIONS[1][1], None, True, size)
        pool.add(SOLUTIONS[2][0], SOLUTIONS[2][1], None, True, size)
        pool.add(SOLUTIONS[3][0], SOLUTIONS[3][1], None, True, size)
        pool.add(SOLUTIONS[4][0], SOLUTIONS[4][1], None, True, size)
        pool.add(SOLUTIONS[5][0], SOLUTIONS[5][1], None, True, size)
        pool.add(SOLUTIONS[6][0], SOLUTIONS[6][1], None, True, size)
        pool.add(SOLUTIONS[7][0], SOLUTIONS[7][1], None, True, size)
        pool.add(SOLUTIONS[8][0], SOLUTIONS[8][1], None, True, size)
        pool.add(SOLUTIONS[9][0], SOLUTIONS[9][1], None, True, size)
        assert pool.length() <= 3, "Pool must keep one solution (length=%s)" % pool.length()
        size = 4
        pool = BestSolutionArchive()
        pool.add(SOLUTIONS[0][0], SOLUTIONS[0][1], None, True, size)
        pool.add(SOLUTIONS[1][0], SOLUTIONS[1][1], None, True, size)
        pool.add(SOLUTIONS[2][0], SOLUTIONS[2][1], None, True, size)
        pool.add(SOLUTIONS[3][0], SOLUTIONS[3][1], None, True, size)
        pool.add(SOLUTIONS[4][0], SOLUTIONS[4][1], None, True, size)
        pool.add(SOLUTIONS[5][0], SOLUTIONS[5][1], None, True, size)
        pool.add(SOLUTIONS[6][0], SOLUTIONS[6][1], None, True, size)
        pool.add(SOLUTIONS[7][0], SOLUTIONS[7][1], None, True, size)
        pool.add(SOLUTIONS[8][0], SOLUTIONS[8][1], None, True, size)
        pool.add(SOLUTIONS[9][0], SOLUTIONS[9][1], None, True, size)
        assert pool.length() <= 4, "Pool must keep one solution (length=%s)" % pool.length()

    def test_callable_pool(self):
        pool = BestSolutionArchive()
        size = 3
        args = {'max_archive_size': 3}
        population = [Individual(SOLUTIONS[0][0], SOLUTIONS[0][1]),
                      Individual(SOLUTIONS[1][0], SOLUTIONS[1][1]),
                      Individual(SOLUTIONS[2][0], SOLUTIONS[2][1]),
                      Individual(SOLUTIONS[3][0], SOLUTIONS[3][1]),
                      Individual(SOLUTIONS[4][0], SOLUTIONS[4][1]),
                      Individual(SOLUTIONS[5][0], SOLUTIONS[5][1]),
                      Individual(SOLUTIONS[6][0], SOLUTIONS[6][1])]
        archive = pool(None, population, [], args)
        assert pool.length() == size

        for sol in pool:
            assert sol in archive


class TestObjectiveFunctions:
    class _MockupSolution:
        def __init__(self):
            self._primal = {}

        def set_primal(self, k, v):
            self._primal[k] = v

        def get_primal_by_id(self, k):
            return self._primal[k]

        @property
        def fluxes(self):
            return self._primal

    def _assert_is_pickable(self, of):
        assert isinstance(pickle.dumps(of), bytes)

    def test_base_yield_function(self, model):
        solution = self._MockupSolution()
        solution.set_primal('EX_ac_lp_e_rp_', 2)
        solution.set_primal('EX_glc_lp_e_rp_', -10)
        of = YieldFunction(model.reactions.EX_ac_lp_e_rp_, model.reactions.EX_glc_lp_e_rp_)
        self._assert_is_pickable(of)

        with pytest.raises(ValueError):
            YieldFunction({}, model.reactions.EX_glc_lp_e_rp_)
        with pytest.raises(ValueError):
            YieldFunction(None, model.reactions.EX_glc_lp_e_rp_)
        with pytest.raises(ValueError):
            YieldFunction([], model.reactions.EX_glc_lp_e_rp_)
        with pytest.raises(ValueError):
            YieldFunction(1, model.reactions.EX_glc_lp_e_rp_)
        with pytest.raises(ValueError):
            YieldFunction(model.reactions.EX_ac_lp_e_rp_, [])
        with pytest.raises(ValueError):
            YieldFunction(model.reactions.EX_ac_lp_e_rp_, 1)
        with pytest.raises(ValueError):
            YieldFunction(model.reactions.EX_ac_lp_e_rp_, {})
        with pytest.raises(ValueError):
            YieldFunction(model.reactions.EX_ac_lp_e_rp_, [])

    def test_biomass_product_coupled_yield(self):
        solution = self._MockupSolution()
        solution.set_primal('biomass', 0.6)
        solution.set_primal('product', 2)
        solution.set_primal('substrate', -10)

        of = biomass_product_coupled_yield("biomass", "product", "substrate")
        assert of.name == "bpcy = (biomass * product) / substrate"
        self._assert_is_pickable(of)
        fitness = of(None, solution, None)
        assert round(abs((0.6 * 2) / 10 - fitness), 7) == 0

        solution.set_primal('substrate', 0)

        fitness = of(None, solution, None)
        assert 0 == fitness

        solution.set_primal('substrate2', -5)
        solution.set_primal('substrate', -5)

        of2 = biomass_product_coupled_yield("biomass", "product", ["substrate", "substrate2"])
        assert of2.name == "bpcy = (biomass * product) / (substrate + substrate2)"
        self._assert_is_pickable(of2)
        fitness = of2(None, solution, None)
        assert round(abs((0.6 * 2) / 10 - fitness), 7) == 0

    def test_biomass_product_coupled_min_yield(self, model):
        biomass = "Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2"
        product = "EX_ac_lp_e_rp_"
        substrate = "EX_glc_lp_e_rp_"
        solution = self._MockupSolution()
        solution.set_primal(biomass, 0.263136)
        solution.set_primal(product, 16.000731)
        solution.set_primal(substrate, -10)

        of = biomass_product_coupled_min_yield(biomass, product, substrate)
        self._assert_is_pickable(of)
        assert of.name == "bpcy = (%s * min(%s)) / %s" % (biomass, product, substrate)
        reactions = [model.reactions.get_by_id(r) for r in ['ATPS4r', 'CO2t', 'GLUDy', 'PPS', 'PYK']]
        with model:
            for r in reactions:
                r.knock_out()
            fitness = of(model, solution, reactions)
        assert round(abs(0.414851 - fitness), 5) == 0

    def test_product_yield(self):
        solution = self._MockupSolution()
        solution.set_primal('biomass', 0.6)
        solution.set_primal('product', 2)
        solution.set_primal('substrate', -10)

        of = product_yield("product", "substrate", carbon_yield=False)
        assert of.name == "yield = (product / substrate)"
        self._assert_is_pickable(of)
        fitness = of(None, solution, None)
        assert round(abs(2.0 / 10.0 - fitness), 7) == 0

        solution.set_primal('substrate', 0)
        fitness = of(None, solution, None)
        assert 0 == fitness

        solution.set_primal('substrate', -5)
        solution.set_primal('substrate2', -5)

        of2 = product_yield('product', ['substrate', 'substrate2'], carbon_yield=False)
        assert of2.name == "yield = (product / (substrate + substrate2))"
        self._assert_is_pickable(of2)
        fitness = of2(None, solution, None)
        assert round(abs(2.0 / 10.0 - fitness), 7) == 0

    def test_number_of_knockouts(self):
        of_max = number_of_knockouts(sense='max')
        assert of_max.name == "max knockouts"
        of_min = number_of_knockouts(sense='min')
        assert of_min.name == "min knockouts"

        f1 = of_max(None, None, ['a', 'b'])
        f2 = of_max(None, None, ['a', 'b', 'c'])
        assert f2 > f1

        f1 = of_min(None, None, ['a', 'b'])
        f2 = of_min(None, None, ['a', 'b', 'c'])
        assert f1 > f2


class TestKnockoutEvaluator:
    def test_initializer(self, model):
        objective1 = biomass_product_coupled_yield(
            "Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2",
            "EX_ac_lp_e_rp_",
            "EX_glc_lp_e_rp_")
        decoder = ReactionSetDecoder(["PGI", "PDH", "FUM", "FBA", "G6PDH2r", "FRD7", "PGL", "PPC"], model)
        evaluator = KnockoutEvaluator(model, decoder, objective1, fba, {})
        assert evaluator.decoder == decoder
        assert evaluator.objective_function == objective1
        assert hasattr(evaluator, "__call__")

        objective2 = product_yield("EX_ac_lp_e_rp_", "EX_glc_lp_e_rp_")
        evaluator = KnockoutEvaluator(model, decoder, MultiObjectiveFunction([objective1, objective2]), fba, {})
        assert evaluator.objective_function.objectives == [objective1, objective2]

    def test_invalid_initializers(self, model):
        objective1 = biomass_product_coupled_yield(
            "Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2",
            "EX_ac_lp_e_rp_",
            "EX_glc_lp_e_rp_")
        decoder = ReactionSetDecoder(["PGI", "PDH", "FUM", "FBA", "G6PDH2r", "FRD7", "PGL", "PPC"], model)
        with pytest.raises(ValueError):
            KnockoutEvaluator(model, decoder, 1, fba, {})
        with pytest.raises(ValueError):
            KnockoutEvaluator(model, decoder, None, fba, {})
        with pytest.raises(ValueError):
            KnockoutEvaluator(model, decoder, [], fba, {})
        with pytest.raises(ValueError):
            KnockoutEvaluator(model, decoder, [2, 3], fba, {})
        with pytest.raises(ValueError):
            KnockoutEvaluator(model, decoder, [objective1], fba, {})
        with pytest.raises(ValueError):
            KnockoutEvaluator(model, None, [], fba, {})
        with pytest.raises(ValueError):
            KnockoutEvaluator(model, True, [], fba, {})

    def test_evaluate_single_objective(self, model):
        representation = ["ATPS4r", "PYK", "GLUDy", "PPS", "CO2t", "PDH",
                          "FUM", "FBA", "G6PDH2r", "FRD7", "PGL", "PPC"]
        decoder = ReactionSetDecoder(representation, model)
        objective1 = biomass_product_coupled_yield(
            "Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2",
            "EX_ac_lp_e_rp_",
            "EX_glc_lp_e_rp_")
        evaluator = KnockoutEvaluator(model, decoder, objective1, fba, {})
        fitness = evaluator([[0, 1, 2, 3, 4]])[0]

        assert abs(fitness - 0.41) < 0.02

    def test_ko_evaluate_single_objective_benchmark(self, benchmark, model):
        benchmark(self.test_evaluate_single_objective, model)

    def test_evaluate_multi_objective(self, model):
        representation = ["ATPS4r", "PYK", "GLUDy", "PPS", "CO2t", "PDH",
                          "FUM", "FBA", "G6PDH2r", "FRD7", "PGL", "PPC"]
        decoder = ReactionSetDecoder(representation, model)
        objective1 = biomass_product_coupled_yield(
            "Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2",
            "EX_ac_lp_e_rp_",
            "EX_glc_lp_e_rp_")
        objective2 = product_yield("EX_ac_lp_e_rp_", "EX_glc_lp_e_rp_", carbon_yield=False)
        objective = MultiObjectiveFunction([objective1, objective2])
        evaluator = KnockoutEvaluator(model, decoder, objective, fba, {})
        fitness = evaluator([[0, 1, 2, 3, 4]])[0]

        assert isinstance(fitness, Pareto)
        assert abs(fitness[0] - 0.41) < 0.02
        assert abs(fitness[1] - 1.57) < 0.035

    def test_ko_evaluate_multi_objective_benchmark(self, benchmark, model):
        benchmark(self.test_evaluate_multi_objective, model)

    def test_evaluate_infeasible_solution(self, model):
        representation = ["ENO", "ATPS4r", "PYK", "GLUDy", "PPS", "CO2t", "PDH",
                          "FUM", "FBA", "G6PDH2r", "FRD7", "PGL", "PPC"]

        decoder = ReactionSetDecoder(representation, model)
        objective1 = biomass_product_coupled_yield(
            "Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2",
            "EX_ac_lp_e_rp_",
            "EX_glc_lp_e_rp_")
        evaluator = KnockoutEvaluator(model, decoder, objective1, fba, {})
        fitness = evaluator([[0]])[0]
        assert fitness == 0


class TestWrappedEvaluator:
    def test_initializer(self):
        def evaluation_function(x):
            return 1

        evaluator = EvaluatorWrapper(config.default_view, evaluation_function)
        assert hasattr(evaluator, '__call__')
        assert hasattr(evaluator, 'view')
        assert hasattr(evaluator, 'evaluator')
        assert evaluator.view == config.default_view
        assert evaluator.evaluator == evaluation_function

    def test_invalid_initializer(self):
        with pytest.raises(ValueError):
            EvaluatorWrapper(config.default_view, None)
        with pytest.raises(ValueError):
            EvaluatorWrapper(config.default_view, 1)
        with pytest.raises(ValueError):
            EvaluatorWrapper(config.default_view, [1, 2, 3])
        with pytest.raises(ValueError):
            EvaluatorWrapper(lambda x: 1, config.default_view)
        with pytest.raises(ValueError):
            EvaluatorWrapper(None, lambda x: 1)
        with pytest.raises(ValueError):
            EvaluatorWrapper(123, lambda x: 1)


class TestSwapOptimization:
    def test_swap_reaction_identification(self, model):
        expected_reactions = ['ACALD', 'AKGDH', 'ALCD2x', 'G6PDH2r', 'GAPD', 'GLUDy', 'GLUSy', 'GND', 'ICDHyr',
                              'LDH_D', 'MDH', 'ME1', 'ME2', 'NADH16', 'PDH']

        swap_pairs = ([model.metabolites.get_by_id(m) for m in NADH_NADPH[0]],
                      [model.metabolites.get_by_id(m) for m in NADH_NADPH[1]])

        representation = CofactorSwapOptimization.find_swappable_reactions(model, swap_pairs)

        assert expected_reactions == representation
        assert 'PGI' not in representation

    def test_evaluate_swap(self, model):
        cofactors = ((model.metabolites.nad_c, model.metabolites.nadh_c),
                     (model.metabolites.nadp_c, model.metabolites.nadph_c))
        with model:
            swap_cofactors(model.reactions.ALCD2x, model, cofactors, inplace=True)
            assert model.metabolites.nadp_c in model.reactions.ALCD2x.metabolites
        assert model.metabolites.nad_c in model.reactions.ALCD2x.metabolites

        with model:
            model.reactions.Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2.lower_bound = 0.5
            py = product_yield(model.reactions.EX_etoh_lp_e_rp_, model.reactions.EX_glc_lp_e_rp_)
            model.objective = model.reactions.EX_etoh_lp_e_rp_
            swap_cofactors(model.reactions.ALCD2x, model, cofactors, inplace=True)
            reactions = ['GAPD', 'AKGDH', 'PDH', 'GLUDy', 'MDH']
            optimization = CofactorSwapOptimization(model=model, objective_function=py,
                                                    candidate_reactions=reactions)
            optimization_result = optimization.run(max_evaluations=10000, max_size=1, pop_size=100,
                                                   variable_size=False,
                                                   mutation_rate=0.5, seed=1485441961)
            # above should not have added anything to the history
            fitness = optimization_result.data_frame.fitness.max()
            assert round(abs(fitness - 0.322085), 3) == 0


class TestDecoders:
    def test_set_decoder(self, model):
        representation = [1, 2, 'a', 'b', None, '0']
        decoder = SetDecoder(representation, model)
        assert decoder([])[0] == []
        for i in range(len(representation)):
            assert decoder([i])[0] == [representation[i]]

    def test_reaction_set_decoder(self, model):
        decoder = ReactionSetDecoder([r.id for r in model.reactions], model)
        reactions = decoder([1, 2, 3, 4])[0]
        for i in range(1, 5):
            assert model.reactions[i] == reactions[i - 1]

    def test_reaction_set_decoder_with_groups(self, model):
        groups = [{model.reactions[1]: 1, model.reactions[11]: 1, model.reactions[12]: 5},
                  {model.reactions[2]: 1, model.reactions[13]: 1, model.reactions[14]: 5}]

        decoder = ReactionSetDecoder([r.id for r in model.reactions[0:10]], model, groups=groups)
        combinations = decoder([1, 2, 3, 4])
        for reactions in combinations:
            for i in range(1, 5):
                reaction = reactions[i - 1]
                group = next((g for g in groups if reaction in g), {reaction: 1})
                assert model.reactions[i] in group

    def test_gene_set_decoder(self, model):
        decoder = GeneSetDecoder([g.id for g in model.genes], model)
        genes = decoder([1, 2, 3, 4])[0]
        for i in range(1, 5):
            assert model.genes[i] == genes[i - 1]


class TestGenerators:
    def test_set_generator(self):
        random = Random(SEED)
        representation = ["a", "b", "c", "d", "e", "f"]
        max_size = 5
        variable_size = False
        expected = [[0, 1, 2, 4, 5],
                    [0, 2, 3, 4, 5],
                    [0, 1, 2, 3, 5],
                    [1, 2, 3, 4, 5],
                    [0, 2, 3, 4, 5]]

        for i in range(len(expected)):
            candidate = set_generator(random, dict(representation=representation,
                                                   max_size=max_size,
                                                   variable_size=variable_size))
            assert candidate == expected[i]

    def test_multiple_chromosome_set_generator(self):
        random = Random(SEED)
        args = dict(keys=["test_key_1", "test_key_2"],
                    test_key_1_representation=["a1", "a2", "a3", "a4", "a5"],
                    test_key_2_representation=["b1", "b2", "b3", "b4", "b5", "b6", "b7"],
                    test_key_1_max_size=3,
                    test_key_2_max_size=5,
                    variable_size=False)
        candidate = multiple_chromosome_set_generator(random, args)
        assert len(candidate['test_key_1']) == 3
        assert len(candidate['test_key_2']) == 5

    def test_fixed_size_set_generator(self, generators):
        args, random, _ = generators
        candidates_file = os.path.join(CURRENT_PATH, "data", "fix_size_candidates.pkl")
        random.seed(SEED)
        args.setdefault('variable_size', False)

        candidates = []

        args['max_size'] = 10
        for _ in range(1000):
            candidate = set_generator(random, args)
            assert len(candidate) == 10
            candidates.append(candidate)

        # with open(candidates_file, 'wb') as out_file:
        #     pickle.dump(candidates, out_file, protocol=2)

        with open(candidates_file, 'rb') as in_file:
            if six.PY3:
                expected_candidates = pickle.load(in_file, encoding="latin1")
            else:
                expected_candidates = pickle.load(in_file)

        assert candidates == expected_candidates

        args['max_size'] = 20
        for _ in range(1000):
            candidate = set_generator(random, args)
            assert len(candidate) == 20

    def test_variable_size_set_generator(self, generators):
        args, random, _ = generators
        candidates_file = os.path.join(CURRENT_PATH, "data", "variable_size_candidates.pkl")
        args.setdefault('variable_size', True)
        random.seed(SEED)
        candidates = []
        args['max_size'] = 10
        for _ in range(1000):
            candidate = set_generator(random, args)
            assert len(candidate) <= 10
            candidates.append(candidate)

        with open(candidates_file, 'rb') as in_file:
            if six.PY3:
                expected_candidates = pickle.load(in_file, encoding="latin1")
            else:
                expected_candidates = pickle.load(in_file)

        assert candidates == expected_candidates

        args['max_size'] = 20
        for _ in range(1000):
            candidate = set_generator(random, args)
            assert len(candidate) <= 20

    def test_fixed_size_linear_set_generator(self, generators):
        args, random, mockup_evolutionary_algorithm = generators
        ec = mockup_evolutionary_algorithm(Bounder(-10, 10))
        args.setdefault('variable_size', False)
        args['max_size'] = 10
        args['_ec'] = ec
        for _ in range(1000):
            candidate = linear_set_generator(random, args)
            for i, v in six.iteritems(candidate):
                assert isinstance(i, (int, numpy.int64, numpy.int32))
                assert isinstance(v, float)

            assert len(candidate) <= 10


class TestHeuristicOptimization:
    def test_default_initializer(self, model, objectives):
        single_objective_function, multi_objective_function = objectives
        heuristic_optimization = HeuristicOptimization(
            model=model,
            objective_function=single_objective_function
        )

        assert heuristic_optimization.model == model
        assert heuristic_optimization.objective_function == single_objective_function

        heuristic_optimization = HeuristicOptimization(
            model=model,
            objective_function=single_objective_function,
        )

        assert heuristic_optimization.model == model
        assert heuristic_optimization.objective_function == single_objective_function

    def test_multi_objective_initializer(self, model, objectives):
        single_objective_function, multi_objective_function = objectives
        heuristic_optimization = HeuristicOptimization(
            model=model,
            objective_function=multi_objective_function,
            heuristic_method=inspyred.ec.emo.NSGA2
        )

        assert heuristic_optimization.model == model
        assert len(heuristic_optimization.objective_function) == 2

        heuristic_optimization = HeuristicOptimization(
            model=model,
            objective_function=multi_objective_function,
            heuristic_method=inspyred.ec.emo.NSGA2,
        )

        assert heuristic_optimization.model == model
        assert len(heuristic_optimization.objective_function) == 2

    def test_invalid_initializer(self, model, objectives):
        single_objective_function, multi_objective_function = objectives
        with pytest.raises(TypeError):
            HeuristicOptimization(model=model,
                                  objective_function=multi_objective_function,
                                  heuristic_method=inspyred.ec.GA)

    def test_single_objective_function_with_multiobjective_initializer(self, model, objectives):
        single_objective_function, multi_objective_function = objectives
        heuristic_optimization = HeuristicOptimization(
            model=model,
            objective_function=single_objective_function,
            heuristic_method=inspyred.ec.emo.NSGA2
        )

        assert len(heuristic_optimization.objective_function) == 1

    def test_change_objective_function(self, model, objectives):
        single_objective_function, multi_objective_function = objectives
        single_objective_heuristic = HeuristicOptimization(
            model=model,
            objective_function=single_objective_function,
        )

        nok = number_of_knockouts()

        single_objective_heuristic.objective_function = nok
        assert nok == single_objective_heuristic.objective_function
        with pytest.raises(TypeError):
            single_objective_heuristic.objective_function(multi_objective_function)

        with pytest.raises(TypeError):
            single_objective_heuristic.objective_function(multi_objective_function)

        multiobjective_heuristic = HeuristicOptimization(
            model=model,
            objective_function=multi_objective_function,
            heuristic_method=inspyred.ec.emo.NSGA2
        )

        multiobjective_heuristic.objective_function = nok
        assert len(multiobjective_heuristic.objective_function) == 1
        assert multiobjective_heuristic.objective_function == nok

    def test_change_heuristic_method(self, model, objectives):
        single_objective_function, multi_objective_function = objectives
        single_objective_heuristic = HeuristicOptimization(
            model=model,
            objective_function=single_objective_function,
        )

        single_objective_heuristic.heuristic_method = inspyred.ec.emo.NSGA2
        assert len(single_objective_heuristic.objective_function) == 1

        multiobjective_heuristic = HeuristicOptimization(
            model=model,
            objective_function=multi_objective_function,
            heuristic_method=inspyred.ec.emo.NSGA2
        )

        with pytest.raises(TypeError):
            multiobjective_heuristic.heuristic_method(inspyred.ec.GA)
        multiobjective_heuristic.objective_function = single_objective_function
        multiobjective_heuristic.heuristic_method = inspyred.ec.GA

    def test_set_distance_function(self):
        s1 = {1, 2, 3}
        s2 = {1, 2, 3, 4}
        d = set_distance_function(s1, s2)
        assert d == 1
        s3 = {2, 3, 4}
        d = set_distance_function(s1, s3)
        assert d == 2
        d = set_distance_function(s3, s2)
        assert d == 1


class TestMigrators:
    @pytest.mark.skipif(RedisQueue is None, reason='redis not available')
    def test_migrator_constructor(self):
        migrator = MultiprocessingMigrator(max_migrants=1, host=REDIS_HOST)
        assert isinstance(migrator.migrants, RedisQueue)
        assert migrator.max_migrants == 1

        migrator = MultiprocessingMigrator(max_migrants=2, host=REDIS_HOST)
        assert isinstance(migrator.migrants, RedisQueue)
        assert migrator.max_migrants == 2

        migrator = MultiprocessingMigrator(max_migrants=3, host=REDIS_HOST)
        assert isinstance(migrator.migrants, RedisQueue)
        assert migrator.max_migrants == 3

    @pytest.mark.skipif(RedisQueue is None, reason='redis not available')
    def test_migrate_individuals_without_evaluation(self):
        population = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        random = Random(SEED)
        migrator = MultiprocessingMigrator(max_migrants=1, host=REDIS_HOST)
        assert isinstance(migrator.migrants, RedisQueue)
        assert migrator.max_migrants == 1

        migrator(random, population, {})
        assert len(migrator.migrants) == 1

        migrator(random, population, {})
        assert len(migrator.migrants) == 1


class TestOptimizationResult:
    def test_reaction_result(self, model):
        representation = [r.id for r in model.reactions]
        random = Random(SEED)
        args = {"representation": representation}

        solutions = BestSolutionArchive()
        for _ in range(10000):
            solutions.add(set_generator(random, args), random.random(), None, True, 100)

        decoder = ReactionSetDecoder(representation, model)

        result = TargetOptimizationResult(
            model=model,
            heuristic_method=None,
            simulation_method=fba,
            simulation_kwargs=None,
            solutions=solutions,
            objective_function=None,
            target_type="reaction",
            decoder=decoder,
            seed=SEED,
            simplify=False)

        assert result.target_type == "reaction"

        individuals = []
        for row in result:
            encoded = set(representation.index(v) for v in row[0])
            individual = Individual(encoded, row[1])
            assert individual not in individuals, "%s is repeated on result"
            individuals.append(individual)
            assert individual in solutions.archive
            assert solutions.archive.count(individual) == 1, "%s is unique in archive" % individual


@pytest.fixture(scope="function")
def reaction_ko_single_objective(model):
    objective = biomass_product_coupled_yield(
        "Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2", "EX_ac_lp_e_rp_", "EX_glc_lp_e_rp_")
    return ReactionKnockoutOptimization(model=model, simulation_method=fba, objective_function=objective)


@pytest.fixture(scope="function")
def reaction_ko_multi_objective(model):
    objective1 = biomass_product_coupled_yield(
        "Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2", "EX_ac_lp_e_rp_", "EX_glc_lp_e_rp_")
    objective2 = number_of_knockouts()
    objective = MultiObjectiveFunction([objective1, objective2])
    return ReactionKnockoutOptimization(model=model, simulation_method=fba, objective_function=objective,
                                        heuristic_method=inspyred.ec.emo.NSGA2)


class TestReactionKnockoutOptimization:
    def test_initializer(self, model):
        essential_reactions = set([r.id for r in find_essential_reactions(model)])
        objective = biomass_product_coupled_yield(
            "Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2", "EX_ac_lp_e_rp_", "EX_glc_lp_e_rp_")
        rko = ReactionKnockoutOptimization(model=model,
                                           simulation_method=fba,
                                           objective_function=objective)

        assert sorted(essential_reactions) == sorted(rko.essential_reactions)
        assert rko._target_type == "reaction"
        assert isinstance(rko._decoder, ReactionSetDecoder)

    def test_run_single_objective(self, reaction_ko_single_objective):
        # TODO: make optlang deterministic so this results can be permanently stored.
        _, result_file = mkstemp('.pkl')
        results = reaction_ko_single_objective.run(max_evaluations=3000, pop_size=10, view=SequentialView(), seed=SEED)
        assert len(results.data_frame.targets) > 0
        assert len(results.data_frame.targets) == len(results.data_frame.targets.apply(tuple).unique())

        with open(result_file, 'wb') as in_file:
            pickle.dump(results, in_file)

        with open(result_file, 'rb') as in_file:
            if six.PY3:
                expected_results = pickle.load(in_file, encoding="latin1")
            else:
                expected_results = pickle.load(in_file)

        assert results.seed == expected_results.seed

    def test_run_reaction_single_ko_objective_benchmark(self, benchmark, reaction_ko_single_objective):
        benchmark(reaction_ko_single_objective.run, max_evaluations=3000, pop_size=10, view=SequentialView(), seed=SEED)

    def test_run_with_time_limit(self, model):
        # TODO: make optlang deterministic so this results can be permanently stored.
        objective = biomass_product_coupled_yield(
            "Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2", "EX_ac_lp_e_rp_", "EX_glc_lp_e_rp_")

        rko = ReactionKnockoutOptimization(model=model,
                                           simulation_method=fba,
                                           objective_function=objective)

        start_time = time.time()
        rko.run(max_evaluations=3000000, pop_size=10, view=SequentialView(), seed=SEED, max_time=(1, 0))
        elapsed_time = time.time() - start_time

        assert elapsed_time < 1.25 * 60

    def test_optgene_with_time_limit(self, model):
        ko = OptGene(model)
        start_time = time.time()
        ko.run(target="EX_ac_lp_e_rp_",
               biomass="Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2",
               substrate="EX_glc_lp_e_rp_", max_evaluations=3000000, seed=SEED, max_time=(1, 0))
        elapsed_time = time.time() - start_time

        # assert elapsed_time < 1.25 * 60
        print(elapsed_time)
        assert elapsed_time < 2 * 60

    def test_run_multi_objective(self, model, reaction_ko_multi_objective):
        # TODO: make optlang deterministic so this results can be permanently stored.
        _, result_file = mkstemp('.pkl')
        results = reaction_ko_multi_objective.run(max_evaluations=3000, pop_size=10, view=SequentialView(), seed=SEED)

        assert len(results.data_frame.targets) == len(results.data_frame.targets.apply(tuple).unique())

        with open(result_file, 'wb') as in_file:
            pickle.dump(results, in_file)

        with open(result_file, 'rb') as in_file:
            if six.PY3:
                expected_results = pickle.load(in_file, encoding="latin1")
            else:
                expected_results = pickle.load(in_file)

        assert results.seed == expected_results.seed

    def test_run_reaction_ko_multi_objective_benchmark(self, benchmark, reaction_ko_multi_objective):
        benchmark(reaction_ko_multi_objective.run, max_evaluations=3000, pop_size=10, view=SequentialView(), seed=SEED)


@pytest.fixture(scope="function")
def gene_ko_single_objective(model):
    objective = biomass_product_coupled_yield(
        "Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2", "EX_ac_lp_e_rp_", "EX_glc_lp_e_rp_")
    return GeneKnockoutOptimization(model=model, simulation_method=fba, objective_function=objective)


@pytest.fixture(scope="function")
def gene_ko_multi_objective(model):
    objective1 = biomass_product_coupled_yield(
        "Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2", "EX_ac_lp_e_rp_", "EX_glc_lp_e_rp_")
    objective2 = number_of_knockouts()
    objective = MultiObjectiveFunction([objective1, objective2])
    return GeneKnockoutOptimization(model=model, simulation_method=fba, objective_function=objective,
                                    heuristic_method=inspyred.ec.emo.NSGA2)


class TestGeneKnockoutOptimization:
    def test_initializer(self, model):
        essential_genes = set([r.id for r in find_essential_genes(model)])
        objective = biomass_product_coupled_yield(
            "Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2", "EX_ac_lp_e_rp_", "EX_glc_lp_e_rp_")
        rko = GeneKnockoutOptimization(model=model,
                                       simulation_method=fba,
                                       objective_function=objective)

        assert sorted(essential_genes) == sorted(rko.essential_genes)
        assert rko._target_type == "gene"
        assert isinstance(rko._decoder, GeneSetDecoder)

    def test_run_single_objective(self, model, gene_ko_single_objective):
        # TODO: make optlang deterministic so this results can be permanently stored.
        _, result_file = mkstemp('.pkl')
        results = gene_ko_single_objective.run(max_evaluations=3000, pop_size=10, view=SequentialView(), seed=SEED)

        assert len(results.data_frame.targets) == len(results.data_frame.targets.apply(tuple).unique())

        with open(result_file, 'wb') as in_file:
            pickle.dump(results, in_file)

        with open(result_file, 'rb') as in_file:
            if six.PY3:
                expected_results = pickle.load(in_file, encoding="latin1")
            else:
                expected_results = pickle.load(in_file)

        assert results.seed == expected_results.seed

    def test_run_gene_ko_single_objective_benchmark(self, gene_ko_single_objective, benchmark):
        benchmark(gene_ko_single_objective.run, max_evaluations=3000, pop_size=10, view=SequentialView(), seed=SEED)

    def test_run_multi_objective(self, model, gene_ko_multi_objective):
        # TODO: make optlang deterministic so this results can be permanently stored.
        _, result_file = mkstemp('.pkl')

        results = gene_ko_multi_objective.run(max_evaluations=3000, pop_size=10, view=SequentialView(), seed=SEED)

        assert len(results.data_frame.targets) == len(results.data_frame.targets.apply(tuple).unique())

        with open(result_file, 'wb') as in_file:
            pickle.dump(results, in_file)

        with open(result_file, 'rb') as in_file:
            if six.PY3:
                expected_results = pickle.load(in_file, encoding="latin1")
            else:
                expected_results = pickle.load(in_file)

        assert results.seed == expected_results.seed

    def test_run_gene_ko_multi_objective_benchmark(self, gene_ko_multi_objective, benchmark):
        benchmark(gene_ko_multi_objective.run, max_evaluations=3000, pop_size=10, view=SequentialView(), seed=SEED)

    def test_run_with_time_limit(self, model):
        # TODO: make optlang deterministic so this results can be permanently stored.
        objective = biomass_product_coupled_yield(
            "Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2", "EX_ac_lp_e_rp_", "EX_glc_lp_e_rp_")

        rko = ReactionKnockoutOptimization(model=model,
                                           simulation_method=fba,
                                           objective_function=objective)

        start_time = time.time()
        rko.run(max_evaluations=3000000, pop_size=10, view=SequentialView(), seed=SEED, max_time=(1, 0))
        elapsed_time = time.time() - start_time

        assert elapsed_time < 1.25 * 60


class TestVariator:
    def test_set_n_point_crossover(self):
        mom = OrderedSet([1, 3, 5, 9, 10])
        dad = OrderedSet([2, 3, 7, 8])
        args = {
            "crossover_rate": 1.0,
            "num_crossover_points": 1,
            "candidate_size": 10
        }
        children = set_n_point_crossover(Random(SEED), [mom, dad], args)
        bro = OrderedSet([1, 3, 5, 8])
        sis = OrderedSet([2, 3, 7, 9, 10])
        assert bro == children[0]
        assert sis == children[1]

    def test_do_not_set_n_point_crossover(self):
        mom = OrderedSet([1, 3, 5, 9, 10])
        dad = OrderedSet([2, 3, 7, 8])
        args = {
            "crossover_rate": 0.0,
            "num_crossover_points": 1,
            "candidate_size": 10
        }
        children = set_n_point_crossover(Random(SEED), [mom, dad], args)
        assert mom == children[0]
        assert dad == children[1]

    def test_set_mutation(self):
        individual = OrderedSet([1, 3, 5, 9, 10])
        representation = list(range(10))
        args = {
            "representation": representation,
            "mutation_rate": 1.0
        }
        new_individuals = set_mutation(Random(SEED), [individual], args)
        assert len(new_individuals[0]) == len(individual)
        assert new_individuals[0] != individual
        assert new_individuals[0] == [0, 2, 4, 6, 7]

    def test_do_not_set_mutation(self):
        individual = OrderedSet([1, 3, 5, 9, 10])
        representation = list(range(10))
        args = {
            "representation": representation,
            "mutation_rate": 0.0
        }
        new_individuals = set_mutation(Random(SEED), [individual], args)
        assert len(new_individuals[0]) == len(individual)
        assert new_individuals[0] == individual

    def test_set_indel(self):
        individual = [1, 3, 5, 9, 10]
        representation = list(range(10))
        args = {
            "representation": representation,
            "indel_rate": 1.0
        }
        new_individuals = set_indel(Random(SEED), [individual], args)
        assert len(new_individuals[0]) != len(individual)
        assert new_individuals[0] == [1, 3, 5, 6, 9, 10]

    def test_do_not_set_indel(self):
        individual = [1, 3, 5, 9, 10]
        representation = list(range(10))
        args = {
            "representation": representation,
            "indel_rate": 0.0
        }
        new_individuals = set_indel(Random(SEED), [individual], args)
        assert len(new_individuals[0]) == len(individual)
        assert new_individuals[0] == individual

        args = {
            "representation": representation,
            "indel_rate": 1.0,
            "variable_size": False
        }
        new_individuals = set_indel(Random(SEED), [individual], args)
        assert len(new_individuals[0]) == len(individual)
        assert new_individuals[0] == individual

    def test_do_set_n_point_crossover(self):
        representation = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]
        int_representation = [representation.index(v) for v in representation]
        mom = OrderedSet([representation.index(v) for v in ["A", "B", "E", "K", "L", "M"]])
        dad = OrderedSet([representation.index(v) for v in ["A", "C", "I", "J", "K", "L"]])
        points = [4]
        children = _do_set_n_point_crossover(int_representation, mom, dad, points, Random(), len(mom))
        bro = OrderedSet([0, 1, 8, 9, 10, 11])
        sis = OrderedSet([0, 2, 4, 10, 11, 12])
        assert children[0] == bro
        assert children[1] == sis

    def test_multiple_chromosome_set_mutation(self):
        genome = MultipleChromosomeGenome(["A", "B"])
        genome["A"] = [1, 2, 3, 4]
        genome["B"] = [1, 5, 7, 10]
        representation = list(range(10))
        args = {
            "A_representation": representation,
            "B_representation": representation,
            "A_mutation_rate": 1,
            "B_mutation_rate": 1
        }

        new_individuals = multiple_chromosome_set_mutation(Random(SEED), [genome], args)
        assert new_individuals[0]["A"] == OrderedSet([0, 6, 7, 8])
        assert new_individuals[0]["B"] == OrderedSet([0, 6, 8, 9])

    def test_multiple_chromosome_set_indel(self):
        genome = MultipleChromosomeGenome(["A", "B"])
        genome["A"] = [1, 2, 3, 4]
        genome["B"] = [1, 5, 7, 10]
        representation = list(range(10))
        args = {
            "A_representation": representation,
            "B_representation": representation,
            "A_indel_rate": 1,
            "B_indel_rate": 1
        }

        random = Random(SEED)
        new_individuals = multiple_chromosome_set_indel(random, [genome for _ in range(5)], args)
        assert new_individuals[0]["A"] == OrderedSet([1, 2, 3, 4, 7])
        assert new_individuals[0]["B"] == OrderedSet([1, 5, 10])
        assert new_individuals[1]["A"] == OrderedSet([2, 3, 4])
        assert new_individuals[1]["B"] == OrderedSet([1, 5, 7, 8, 10])
        assert new_individuals[2]["A"] == OrderedSet([1, 2, 3, 4, 6])
        assert new_individuals[2]["B"] == OrderedSet([1, 5, 7])
        assert new_individuals[3]["A"] == OrderedSet([1, 2, 3, 4, 8])
        assert new_individuals[3]["B"] == OrderedSet([0, 1, 5, 7, 10])
        assert new_individuals[4]["A"] == OrderedSet([1, 2, 3, 4, 7])
        assert new_individuals[4]["B"] == OrderedSet([1, 5, 7, 8, 10])


class TestGenomes:
    def test_two_chromosomes(self):
        genome = MultipleChromosomeGenome(["A", "B"])
        assert isinstance(genome["A"], list)
        assert isinstance(genome["B"], list)
        genome["A"] = [1, 2, 3, 4]
        genome["B"] = ["A", "B", "C"]

        assert genome["A"] == OrderedSet([1, 2, 3, 4])
        assert genome["B"] == OrderedSet(["A", "B", "C"])

        del genome["A"]
        with pytest.raises(KeyError):
            genome.__getitem__("A")


def simplify_knockout_solutions_for_succ(iaf1260):
    representation = ["FUM", "SFGTHi", "DHACOAH", "ASPTRS"]
    solution = [0, 1, 2, 3]

    bpcy = biomass_product_coupled_min_yield("Ec_biomass_iAF1260_core_59p81M",
                                             "EX_succ_lp_e_rp_",
                                             "EX_glc_lp_e_rp_")

    decoder = ReactionSetDecoder(representation, iaf1260)
    evaluator = KnockoutEvaluator(iaf1260, decoder, bpcy, fba, {})
    simplification = SolutionSimplification(evaluator)

    new_solution = simplification(solution)
    assert [0] == new_solution
