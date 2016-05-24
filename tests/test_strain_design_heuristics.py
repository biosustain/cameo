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
import unittest
from collections import namedtuple
from math import sqrt

import inspyred
import numpy
import six
from cobra.manipulation.delete import find_gene_knockout_reactions
from inspyred.ec import Bounder
from ordered_set import OrderedSet
from pandas.util.testing import assert_frame_equal
from six.moves import range

from cameo import load_model, fba
from cameo.parallel import SequentialView, RedisQueue
from cameo.strain_design.heuristic.evolutionary.archives import Individual, BestSolutionArchive
from cameo.strain_design.heuristic.evolutionary.decoders import ReactionKnockoutDecoder, KnockoutDecoder, \
    GeneKnockoutDecoder
from cameo.strain_design.heuristic.evolutionary.generators import set_generator, unique_set_generator, \
    multiple_chromosome_set_generator, linear_set_generator
from cameo.strain_design.heuristic.evolutionary.genomes import MultipleChromosomeGenome
from cameo.strain_design.heuristic.evolutionary.metrics import euclidean_distance
from cameo.strain_design.heuristic.evolutionary.metrics import manhattan_distance
from cameo.strain_design.heuristic.evolutionary.multiprocess.migrators import MultiprocessingMigrator
from cameo.strain_design.heuristic.evolutionary.objective_functions import biomass_product_coupled_yield, \
    product_yield, number_of_knockouts
from cameo.strain_design.heuristic.evolutionary.optimization import HeuristicOptimization, \
    ReactionKnockoutOptimization, set_distance_function, KnockoutOptimizationResult
from cameo.strain_design.heuristic.evolutionary.variators import _do_set_n_point_crossover, set_n_point_crossover, \
    set_mutation, set_indel, multiple_chromosome_set_mutation, multiple_chromosome_set_indel
from cameo.util import RandomGenerator as Random

TRAVIS = os.getenv('TRAVIS', False)

if os.getenv('REDIS_PORT_6379_TCP_ADDR'):
    REDIS_HOST = os.getenv('REDIS_PORT_6379_TCP_ADDR')  # wercker
else:
    REDIS_HOST = 'localhost'

SEED = 1234

CURRENT_PATH = os.path.dirname(__file__)
MODEL_PATH = os.path.join(CURRENT_PATH, "data/EcoliCore.xml")

TEST_MODEL = load_model(MODEL_PATH, sanitize=False)

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


class TestMetrics(unittest.TestCase):
    def test_euclidean_distance(self):
        distance = euclidean_distance({'a': 9}, {'a': 3})
        self.assertEqual(distance, sqrt((9 - 3) ** 2))

    def test_manhattan_distance(self):
        distance = manhattan_distance({'a': 9}, {'a': 3})
        self.assertEqual(distance, abs(9 - 3))


class TestBestSolutionArchive(unittest.TestCase):
    def test_solution_string(self):
        sol1 = Individual(SOLUTIONS[0][0], SOLUTIONS[0][1])
        sol2 = Individual(SOLUTIONS[1][0], SOLUTIONS[1][1])
        sol3 = Individual(SOLUTIONS[2][0], SOLUTIONS[2][1])
        self.assertEqual(sol1.__str__(), "[1, 2, 3] - 0.1 sense: max")
        self.assertEqual(sol2.__str__(), "[1, 2, 3, 4] - 0.1 sense: max")
        self.assertEqual(sol3.__str__(), "[2, 3, 4] - 0.45 sense: max")

    def test_solution_comparison_maximization(self):
        sol1 = Individual(SOLUTIONS[0][0], SOLUTIONS[0][1])
        sol2 = Individual(SOLUTIONS[1][0], SOLUTIONS[1][1])
        sol3 = Individual(SOLUTIONS[2][0], SOLUTIONS[2][1])

        # test ordering
        self.assertEqual(sol1.__cmp__(sol2), -1)
        self.assertEqual(sol1.__cmp__(sol1), 0)
        self.assertEqual(sol1.__cmp__(sol3), 1)

        self.assertTrue(sol1 < sol2)
        self.assertTrue(sol1 == sol1)
        self.assertTrue(sol1 > sol3)

        # test gt and lt
        self.assertTrue(sol1.__lt__(sol2))
        self.assertTrue(sol1.__gt__(sol3))
        self.assertFalse(sol1.__lt__(sol1))
        self.assertFalse(sol1.__gt__(sol1))
        self.assertFalse(sol2.__lt__(sol1))
        self.assertFalse(sol3.__gt__(sol1))

        # testing issubset
        self.assertTrue(sol1.issubset(sol2), msg="Solution 1 is subset of Solution 2")
        self.assertFalse(sol2.issubset(sol1), msg="Solution 2 is not subset of Solution 1")
        self.assertTrue(sol3.issubset(sol2), msg="Solution 3 is subset of Solution 2")
        self.assertFalse(sol2.issubset(sol3), msg="Solution 2 is not subset of Solution 3")
        self.assertFalse(sol1.issubset(sol3), msg="Solution 1 is subset of Solution 3")
        self.assertFalse(sol2.issubset(sol3), msg="Solution 3 is not subset of Solution 1")

        # test difference
        l = len(sol2.symmetric_difference(sol1))
        self.assertEqual(l, 1, msg="Difference between Solution 2 and 1 is (%s)" % sol2.symmetric_difference(sol1))
        l = len(sol3.symmetric_difference(sol2))
        self.assertEqual(l, 1, msg="Difference between Solution 3 and 1 is (%s)" % sol3.symmetric_difference(sol2))
        l = len(sol3.symmetric_difference(sol1))
        self.assertEqual(l, 2, msg="Difference between Solution 1 and 3 is (%s)" % sol3.symmetric_difference(sol1))

        self.assertTrue(sol1.improves(sol2), msg="Solution 1 is better than Solution 2")
        self.assertTrue(sol3.improves(sol2), msg="Solution 3 is better than Solution 2")
        self.assertFalse(sol3.improves(sol1), msg="Solution 3 does not improve Solution 1")
        self.assertFalse(sol2.improves(sol1), msg="Solution 2 does not improve Solution 1")
        self.assertFalse(sol2.improves(sol3), msg="Solution 2 does not improve Solution 3")

    def test_solution_comparison_minimization(self):
        sol1 = Individual(SOLUTIONS[0][0], SOLUTIONS[0][1], maximize=False)
        sol2 = Individual(SOLUTIONS[1][0], SOLUTIONS[1][1], maximize=False)
        sol3 = Individual(SOLUTIONS[2][0], SOLUTIONS[2][1], maximize=False)

        # test ordering
        self.assertEqual(sol1.__cmp__(sol2), -1)
        self.assertEqual(sol1.__cmp__(sol1), 0)
        self.assertEqual(sol1.__cmp__(sol3), -1)
        self.assertEqual(sol3.__cmp__(sol1), 1)

        self.assertTrue(sol1 < sol2)
        self.assertTrue(sol1 == sol1)
        self.assertTrue(sol1 < sol3)

        # test gt and lt
        self.assertTrue(sol1.__lt__(sol2))
        self.assertTrue(sol1.__lt__(sol3))
        self.assertFalse(sol1.__gt__(sol1))
        self.assertFalse(sol1.__lt__(sol1))
        self.assertTrue(sol2.__gt__(sol1))
        self.assertFalse(sol3.__lt__(sol1))

        # testing issubset
        self.assertTrue(sol1.issubset(sol2), msg="Solution 1 is subset of Solution 2")
        self.assertFalse(sol2.issubset(sol1), msg="Solution 2 is not subset of Solution 1")
        self.assertTrue(sol3.issubset(sol2), msg="Solution 3 is subset of Solution 2")
        self.assertFalse(sol2.issubset(sol3), msg="Solution 2 is not subset of Solution 3")
        self.assertFalse(sol1.issubset(sol3), msg="Solution 1 is subset of Solution 3")
        self.assertFalse(sol2.issubset(sol3), msg="Solution 3 is not subset of Solution 1")

        # test difference
        l = len(sol2.symmetric_difference(sol1))
        self.assertEqual(l, 1, msg="Difference between Solution 2 and 1 is (%s)" % sol2.symmetric_difference(sol1))
        l = len(sol3.symmetric_difference(sol2))
        self.assertEqual(l, 1, msg="Difference between Solution 3 and 1 is (%s)" % sol3.symmetric_difference(sol2))
        l = len(sol3.symmetric_difference(sol1))
        self.assertEqual(l, 2, msg="Difference between Solution 1 and 3 is (%s)" % sol3.symmetric_difference(sol1))

        self.assertTrue(sol1.improves(sol2), msg="Solution 1 is better than Solution 2")
        self.assertFalse(sol3.improves(sol2), msg="Solution 3 is not better than Solution 2")
        self.assertFalse(sol3.improves(sol1), msg="Solution 3 does not improve Solution 1")
        self.assertFalse(sol2.improves(sol1), msg="Solution 2 does not improve Solution 1")
        self.assertFalse(sol2.improves(sol3), msg="Solution 2 does not improve Solution 3")

    def test_add_greater_solution_with_same_fitness(self):
        size = 1
        pool = BestSolutionArchive()
        pool.add(SOLUTIONS[0][0], SOLUTIONS[0][1], None, True, size)
        pool.add(SOLUTIONS[1][0], SOLUTIONS[1][1], None, True, size)
        self.assertEqual(pool.length(), 1, msg="Pool must keep one solution (length=%s)" % pool.length())
        best_solution = set(SOLUTIONS[0][0])
        best_fitness = SOLUTIONS[0][1]
        sol = pool.get(0)
        self.assertEqual(sol.candidate, best_solution, msg="Best solution set must be the first")
        self.assertEqual(sol.fitness, best_fitness, msg="Best solution fitness must be the first")

    def test_add_smaller_solution_with_same_fitness(self):
        size = 1
        pool = BestSolutionArchive()
        pool.add(SOLUTIONS[1][0], SOLUTIONS[1][1], None, True, size)
        pool.add(SOLUTIONS[0][0], SOLUTIONS[0][1], None, True, size)
        self.assertEqual(pool.length(), 1, msg="Pool must keep one solution (length=%s)" % pool.length())
        solution = set(SOLUTIONS[0][0])
        fitness = SOLUTIONS[0][1]
        sol = pool.get(0)
        self.assertEqual(sol.candidate, solution, msg="Best solution must be the first (%s)" % sol.candidate)
        self.assertEqual(sol.fitness, fitness, msg="Best fitness must be the first (%s)" % sol.fitness)

    def test_uniqueness_of_solutions(self):
        size = 2
        pool = BestSolutionArchive()
        pool.add(SOLUTIONS[1][0], SOLUTIONS[1][1], None, True, size)
        pool.add(SOLUTIONS[1][0], SOLUTIONS[1][1], None, True, size)

        self.assertEqual(pool.length(), 1, "Added repeated solution")

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
        self.assertLessEqual(pool.length(), 1, msg="Pool must keep one solution (length=%s)" % pool.length())
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
        self.assertLessEqual(pool.length(), 2, msg="Pool must keep one solution (length=%s)" % pool.length())
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
        self.assertLessEqual(pool.length(), 3, msg="Pool must keep one solution (length=%s)" % pool.length())
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
        self.assertLessEqual(pool.length(), 4, msg="Pool must keep one solution (length=%s)" % pool.length())

    def test_callable_pool(self):
        pool = BestSolutionArchive()
        size = 3
        args = {}
        args['max_archive_size'] = size
        population = [Individual(SOLUTIONS[0][0], SOLUTIONS[0][1]),
                      Individual(SOLUTIONS[1][0], SOLUTIONS[1][1]),
                      Individual(SOLUTIONS[2][0], SOLUTIONS[2][1]),
                      Individual(SOLUTIONS[3][0], SOLUTIONS[3][1]),
                      Individual(SOLUTIONS[4][0], SOLUTIONS[4][1]),
                      Individual(SOLUTIONS[5][0], SOLUTIONS[5][1]),
                      Individual(SOLUTIONS[6][0], SOLUTIONS[6][1])]
        archive = pool(None, population, [], args)
        self.assertEqual(pool.length(), size)

        for sol in pool:
            self.assertTrue(sol in archive)


class TestObjectiveFunctions(unittest.TestCase):
    class _MockupSolution():
        def __init__(self):
            self._primal = {}

        def set_primal(self, k, v):
            self._primal[k] = v

        def get_primal_by_id(self, k):
            return self._primal[k]

        @property
        def fluxes(self):
            return self._primal

    def test_biomass_product_coupled_yield(self):
        solution = self._MockupSolution()
        solution.set_primal('biomass', 0.6)
        solution.set_primal('product', 2)
        solution.set_primal('substrate', -10)

        of = biomass_product_coupled_yield("biomass", "product", "substrate")
        self.assertEqual(of.name, "bpcy = (biomass * product) / substrate")

        fitness = of(None, solution, None)
        self.assertAlmostEqual((0.6 * 2) / 10, fitness)

        solution.set_primal('substrate', 0)

        fitness = of(None, solution, None)
        self.assertEquals(0, fitness)

    def test_yield(self):
        solution = self._MockupSolution()
        solution.set_primal('biomass', 0.6)
        solution.set_primal('product', 2)
        solution.set_primal('substrate', -10)

        of = product_yield("product", "substrate")
        self.assertEqual(of.name, "yield = (product / substrate)")
        fitness = of(None, solution, None)
        self.assertAlmostEqual(2.0 / 10.0, fitness)

        solution.set_primal('substrate', 0)
        fitness = of(None, solution, None)
        self.assertEquals(0, fitness)

    def test_number_of_knockouts(self):
        of_max = number_of_knockouts(sense='max')
        self.assertEqual(of_max.name, "max knockouts")
        of_min = number_of_knockouts(sense='min')
        self.assertEqual(of_min.name, "min knockouts")

        f1 = of_max(None, None, [['a', 'b'], ['a', 'b']])
        f2 = of_max(None, None, [['a', 'b'], ['a', 'b', 'c']])
        self.assertGreater(f2, f1)

        f1 = of_min(None, None, [['a', 'b'], ['a', 'b']])
        f2 = of_min(None, None, [['a', 'b'], ['a', 'b', 'c']])
        self.assertGreater(f1, f2)


class TestDecoders(unittest.TestCase):
    def setUp(self):
        self.model = TEST_MODEL

    def test_abstract_decoder(self):
        decoder = KnockoutDecoder(None, self.model)
        self.assertRaises(NotImplementedError, decoder, [])

    def test_reaction_knockout_decoder(self):
        decoder = ReactionKnockoutDecoder([r.id for r in self.model.reactions], self.model)
        reactions1, reactions2 = decoder([1, 2, 3, 4])
        self.assertTrue(sorted(reactions1, key=lambda x: x.id) == sorted(reactions2, key=lambda x: x.id))

    def test_gene_knockout_decoder(self):
        decoder = GeneKnockoutDecoder([g.id for g in self.model.genes], self.model)
        reactions1, genes = decoder([1, 2, 3, 4])
        reactions2 = find_gene_knockout_reactions(self.model, genes)
        self.assertTrue(sorted(reactions1, key=lambda x: x.id) == sorted(reactions2, key=lambda x: x.id))


class TestGenerators(unittest.TestCase):
    mockup_evolutionary_algorithm = namedtuple("EA", ["bounder"])

    def setUp(self):
        self.model = TEST_MODEL
        self.args = {}
        self.args.setdefault('representation', [r.id for r in self.model.reactions])
        self.random = Random()

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
            self.assertEqual(candidate, expected[i])

    def test_multiple_chromossome_set_generator(self):
        random = Random(SEED)
        args = dict(keys=["test_key_1", "test_key_2"],
                    test_key_1_representation=["a1", "a2", "a3", "a4", "a5"],
                    test_key_2_representation=["b1", "b2", "b3", "b4", "b5", "b6", "b7"],
                    test_key_1_max_size=3,
                    test_key_2_max_size=5,
                    variable_size=False)
        candidate = multiple_chromosome_set_generator(random, args)
        self.assertEqual(len(candidate['test_key_1']), 3)
        self.assertEqual(len(candidate['test_key_2']), 5)

    def test_fixed_size_set_generator(self):
        self.args.setdefault('variable_size', False)

        self.args['max_size'] = 10
        for _ in range(1000):
            candidate = set_generator(self.random, self.args)
            self.assertEqual(len(candidate), 10)
            candidate = unique_set_generator(self.random, self.args)
            self.assertEqual(len(candidate), 10)

        self.args['max_size'] = 20
        for _ in range(1000):
            candidate = set_generator(self.random, self.args)
            self.assertEqual(len(candidate), 20)
            candidate = unique_set_generator(self.random, self.args)
            self.assertEqual(len(candidate), 20)

    def test_variable_size_set_generator(self):
        self.args.setdefault('variable_size', True)

        self.args['max_size'] = 10
        for _ in range(1000):
            candidate = set_generator(self.random, self.args)
            self.assertLessEqual(len(candidate), 10)
            candidate = unique_set_generator(self.random, self.args)
            self.assertLessEqual(len(candidate), 10)

        self.args['max_size'] = 20
        for _ in range(1000):
            candidate = set_generator(self.random, self.args)
            self.assertLessEqual(len(candidate), 20)
            candidate = unique_set_generator(self.random, self.args)
            self.assertLessEqual(len(candidate), 20)

    def test_fixed_size_linear_set_generator(self):
        ec = self.mockup_evolutionary_algorithm(Bounder(-10, 10))
        self.args.setdefault('variable_size', True)
        self.args['max_size'] = 10
        self.args['_ec'] = ec
        for _ in range(1000):
            candidate = linear_set_generator(self.random, self.args)
            for i, v in six.iteritems(candidate):
                self.assertIsInstance(i, (int, numpy.int64, numpy.int32))
                self.assertIsInstance(v, float)

            self.assertLessEqual(len(candidate), 10)


class TestHeuristicOptimization(unittest.TestCase):
    def setUp(self):
        self.model = TEST_MODEL
        self.single_objective_function = product_yield('product', 'substrate')
        self.multiobjective_function = [
            product_yield('product', 'substrate'),
            number_of_knockouts()
        ]

    def test_default_initializer(self):
        heuristic_optimization = HeuristicOptimization(
            model=self.model,
            objective_function=self.single_objective_function
        )

        self.assertIsNone(heuristic_optimization._generator)
        self.assertIsNot(heuristic_optimization.seed, SEED)
        self.assertEqual(heuristic_optimization.model, self.model)
        self.assertEqual(heuristic_optimization.objective_function, self.single_objective_function)

        heuristic_optimization = HeuristicOptimization(
            model=self.model,
            objective_function=self.single_objective_function,
            seed=SEED
        )

        self.assertIsNone(heuristic_optimization._generator)
        self.assertEqual(heuristic_optimization.seed, SEED)
        self.assertEqual(heuristic_optimization.model, self.model)
        self.assertEqual(heuristic_optimization.objective_function, self.single_objective_function)

    def test_multiobjective_initializer(self):
        heuristic_optimization = HeuristicOptimization(
            model=self.model,
            objective_function=self.multiobjective_function,
            heuristic_method=inspyred.ec.emo.NSGA2
        )

        self.assertIsNone(heuristic_optimization._generator)
        self.assertIsNot(heuristic_optimization.seed, SEED)
        self.assertEqual(heuristic_optimization.model, self.model)
        self.assertEqual(len(heuristic_optimization.objective_function), 2)

        heuristic_optimization = HeuristicOptimization(
            model=self.model,
            objective_function=self.multiobjective_function,
            heuristic_method=inspyred.ec.emo.NSGA2,
            seed=SEED
        )

        self.assertIsNone(heuristic_optimization._generator)
        self.assertEqual(heuristic_optimization.seed, SEED)
        self.assertEqual(heuristic_optimization.model, self.model)
        self.assertEqual(len(heuristic_optimization.objective_function), 2)

    def test_invalid_initializer(self):
        self.assertRaises(TypeError, HeuristicOptimization,
                          model=self.model,
                          objective_function=self.multiobjective_function,
                          heuristic_method=inspyred.ec.GA)

    def test_single_objective_function_with_multiobjective_initializer(self):
        heuristic_optimization = HeuristicOptimization(
            model=self.model,
            objective_function=self.single_objective_function,
            heuristic_method=inspyred.ec.emo.NSGA2
        )

        self.assertEqual(len(heuristic_optimization.objective_function), 1)

    def test_change_objective_function(self):
        single_objective_heuristic = HeuristicOptimization(
            model=self.model,
            objective_function=self.single_objective_function,
        )

        nok = number_of_knockouts()

        single_objective_heuristic.objective_function = nok
        self.assertEqual(nok, single_objective_heuristic.objective_function)
        self.assertFalse(single_objective_heuristic.is_mo())
        self.assertRaises(TypeError,
                          single_objective_heuristic.objective_function,
                          self.multiobjective_function)

        single_objective_heuristic.objective_function = [nok]
        self.assertEqual(nok, single_objective_heuristic.objective_function)
        self.assertFalse(single_objective_heuristic.is_mo())
        self.assertRaises(TypeError, single_objective_heuristic.objective_function, self.multiobjective_function)

        multiobjective_heuristic = HeuristicOptimization(
            model=self.model,
            objective_function=self.multiobjective_function,
            heuristic_method=inspyred.ec.emo.NSGA2
        )

        multiobjective_heuristic.objective_function = nok
        self.assertEqual(len(multiobjective_heuristic.objective_function), 1)
        self.assertEqual(multiobjective_heuristic.objective_function[0], nok)
        self.assertTrue(multiobjective_heuristic.is_mo())

    def test_change_heuristic_method(self):
        single_objective_heuristic = HeuristicOptimization(
            model=self.model,
            objective_function=self.single_objective_function,
        )

        single_objective_heuristic.heuristic_method = inspyred.ec.emo.NSGA2
        self.assertTrue(single_objective_heuristic.is_mo())
        self.assertEqual(len(single_objective_heuristic.objective_function), 1)

        multiobjective_heuristic = HeuristicOptimization(
            model=self.model,
            objective_function=self.multiobjective_function,
            heuristic_method=inspyred.ec.emo.NSGA2
        )

        self.assertRaises(TypeError, multiobjective_heuristic.heuristic_method, inspyred.ec.GA)
        multiobjective_heuristic.objective_function = self.single_objective_function
        multiobjective_heuristic.heuristic_method = inspyred.ec.GA
        self.assertFalse(multiobjective_heuristic.is_mo())

    def test_set_distance_function(self):
        s1 = {1, 2, 3}
        s2 = {1, 2, 3, 4}
        d = set_distance_function(s1, s2)
        self.assertEqual(d, 1)
        s3 = {2, 3, 4}
        d = set_distance_function(s1, s3)
        self.assertEqual(d, 2)
        d = set_distance_function(s3, s2)
        self.assertEqual(d, 1)


class TestMigrators(unittest.TestCase):
    def setUp(self):
        self.population = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        self.random = Random(SEED)

    # unittest.skipIf(os.getenv('WERCKER', False), 'Currently not working on wercker as redis is not running on localhost')
    def test_migrator_constructor(self):
        migrator = MultiprocessingMigrator(max_migrants=1, host=REDIS_HOST)
        self.assertIsInstance(migrator.migrants, RedisQueue)
        self.assertEqual(migrator.max_migrants, 1)

        migrator = MultiprocessingMigrator(max_migrants=2, host=REDIS_HOST)
        self.assertIsInstance(migrator.migrants, RedisQueue)
        self.assertEqual(migrator.max_migrants, 2)

        migrator = MultiprocessingMigrator(max_migrants=3, host=REDIS_HOST)
        self.assertIsInstance(migrator.migrants, RedisQueue)
        self.assertEqual(migrator.max_migrants, 3)

    # unittest.skipIf(os.getenv('WERCKER', False), 'Currently not working on wercker as redis is not running on localhost')
    def test_migrate_individuals_without_evaluation(self):
        migrator = MultiprocessingMigrator(max_migrants=1, host=REDIS_HOST)
        self.assertIsInstance(migrator.migrants, RedisQueue)
        self.assertEqual(migrator.max_migrants, 1)

        migrator(self.random, self.population, {})
        self.assertEqual(len(migrator.migrants), 1)

        migrator(self.random, self.population, {})
        self.assertEqual(len(migrator.migrants), 1)


class TestKnockoutOptimizationResult(unittest.TestCase):
    def setUp(self):
        self.model = TEST_MODEL
        self.representation = [r.id for r in self.model.reactions]
        random = Random(SEED)
        args = {"representation": self.representation}
        self.solutions = BestSolutionArchive()
        for _ in range(10000):
            self.solutions.add(set_generator(random, args), random.random(), None, True, 100)
        self.decoder = ReactionKnockoutDecoder(self.representation, self.model)

    def test_reaction_result(self):
        result = KnockoutOptimizationResult(
            model=self.model,
            heuristic_method=None,
            simulation_method=fba,
            simulation_kwargs=None,
            solutions=self.solutions,
            objective_function=None,
            ko_type="reaction",
            decoder=self.decoder,
            seed=SEED)

        self.assertEqual(result.ko_type, "reaction")

        individuals = []
        for row in result:
            encoded = set(self.representation.index(v) for v in row[0])
            individual = Individual(encoded, row[2])
            self.assertNotIn(individual, individuals, msg="%s is repeated on result")
            individuals.append(individual)
            self.assertIn(individual, self.solutions.archive)
            self.assertEqual(self.solutions.archive.count(individual), 1, msg="%s is unique in archive" % individual)


class TestReactionKnockoutOptimization(unittest.TestCase):
    def setUp(self):
        self.model = TEST_MODEL
        self.essential_reactions = set([r.id for r in self.model.essential_reactions()])

    def test_initialize(self):
        rko = ReactionKnockoutOptimization(model=self.model,
                                           simulation_method=fba,
                                           seed=SEED)

        self.assertTrue(sorted(self.essential_reactions) == sorted(rko.essential_reactions))
        self.assertEqual(rko._ko_type, "reaction")
        self.assertTrue(isinstance(rko._decoder, ReactionKnockoutDecoder))

    @unittest.skipIf(True, 'Broken ..')
    def test_run_single_objective(self):
        result_file = os.path.join(CURRENT_PATH, "data", "reaction_knockout_single_objective.pkl")
        objective = biomass_product_coupled_yield(
            "Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2",
            "EX_ac_LPAREN_e_RPAREN_",
            "EX_glc_LPAREN_e_RPAREN_")

        rko = ReactionKnockoutOptimization(model=self.model,
                                           simulation_method=fba,
                                           objective_function=objective,
                                           seed=SEED)

        # self.assertEqual(rko.random.random(), 0.1915194503788923)

        results = rko.run(max_evaluations=3000, pop_size=10, view=SequentialView())

        # self.assertEqual(rko.random.random(), 0.9268454219291495)

        with open(result_file, 'rb') as in_file:
            if six.PY3:
                expected_results = pickle.load(in_file, encoding="latin1")
            else:
                expected_results = pickle.load(in_file)

        assert_frame_equal(results.solutions, expected_results.solutions)

    @unittest.skipIf(True, 'Broken ..')
    def test_run_multiobjective(self):
        result_file = os.path.join(CURRENT_PATH, "data", "reaction_knockout_multi_objective.pkl")
        objective1 = biomass_product_coupled_yield(
            "Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2",
            "EX_ac_LPAREN_e_RPAREN_",
            "EX_glc_LPAREN_e_RPAREN_")

        objective2 = number_of_knockouts()
        objective = [objective1, objective2]

        rko = ReactionKnockoutOptimization(model=self.model,
                                           simulation_method=fba,
                                           objective_function=objective,
                                           heuristic_method=inspyred.ec.emo.NSGA2,
                                           seed=SEED)

        # self.assertEqual(rko.random.random(), 0.1915194503788923)

        results = rko.run(max_evaluations=3000, pop_size=10, view=SequentialView())

        # print(rko.random.random(), 0.545818634701)

        with open(result_file, 'rb') as in_file:
            if six.PY3:
                expected_results = pickle.load(in_file, encoding="latin1")
            else:
                expected_results = pickle.load(in_file)

        assert_frame_equal(results.solutions, expected_results.solutions)

    def test_evaluator(self):
        pass


class VariatorTestCase(unittest.TestCase):
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
        self.assertEqual(bro, children[0])
        self.assertEqual(sis, children[1])

    def test_set_mutation(self):
        individual = OrderedSet([1, 3, 5, 9, 10])
        representation = list(range(10))
        args = {
            "representation": representation,
            "mutation_rate": 1
        }
        new_individuals = set_mutation(Random(SEED), [individual], args)
        self.assertEqual(len(new_individuals[0]), len(individual))
        self.assertNotEqual(new_individuals[0], individual)
        self.assertEqual(new_individuals[0], [0, 2, 4, 6, 7])

    def test_set_indel(self):
        individual = [1, 3, 5, 9, 10]
        representation = list(range(10))
        args = {
            "representation": representation,
            "indel_rate": 1
        }
        new_individuals = set_indel(Random(SEED), [individual], args)
        self.assertNotEqual(len(new_individuals[0]), len(individual))
        self.assertEqual(new_individuals[0], [1, 3, 5, 6, 9, 10])

    def test_do_set_n_point_crossover(self):
        representation = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]
        int_representation = [representation.index(v) for v in representation]
        mom = OrderedSet([representation.index(v) for v in ["A", "B", "E", "K", "L", "M"]])
        dad = OrderedSet([representation.index(v) for v in ["A", "C", "I", "J", "K", "L"]])
        points = [4]
        children = _do_set_n_point_crossover(int_representation, mom, dad, points, Random(), len(mom))
        bro = OrderedSet([0, 1, 8, 9, 10, 11])
        sis = OrderedSet([0, 2, 4, 10, 11, 12])
        self.assertEqual(children[0], bro)
        self.assertEqual(children[1], sis)

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
        self.assertEqual(new_individuals[0]["A"], OrderedSet([0, 6, 7, 8]))
        self.assertEqual(new_individuals[0]["B"], OrderedSet([0, 6, 8, 9]))

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
        self.assertEqual(new_individuals[0]["A"], OrderedSet([1, 2, 3, 4, 7]))
        self.assertEqual(new_individuals[0]["B"], OrderedSet([1, 5, 10]))
        self.assertEqual(new_individuals[1]["A"], OrderedSet([2, 3, 4]))
        self.assertEqual(new_individuals[1]["B"], OrderedSet([1, 5, 7, 8, 10]))
        self.assertEqual(new_individuals[2]["A"], OrderedSet([1, 2, 3, 4, 6]))
        self.assertEqual(new_individuals[2]["B"], OrderedSet([1, 5, 7]))
        self.assertEqual(new_individuals[3]["A"], OrderedSet([1, 2, 3, 4, 8]))
        self.assertEqual(new_individuals[3]["B"], OrderedSet([0, 1, 5, 7, 10]))
        self.assertEqual(new_individuals[4]["A"], OrderedSet([1, 2, 3, 4, 7]))
        self.assertEqual(new_individuals[4]["B"], OrderedSet([1, 5, 7, 8, 10]))


class GenomesTestCase(unittest.TestCase):
    def test_two_chromosomes(self):
        genome = MultipleChromosomeGenome(["A", "B"])
        self.assertIsInstance(genome["A"], OrderedSet)
        self.assertIsInstance(genome["B"], OrderedSet)
        genome["A"] = [1, 2, 3, 4]
        genome["B"] = ["A", "B", "C"]

        self.assertEqual(genome["A"], OrderedSet([1, 2, 3, 4]))
        self.assertEqual(genome["B"], OrderedSet(["A", "B", "C"]))

        del genome["A"]
        self.assertRaises(KeyError, genome.__getitem__, "A")
