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

import os
import unittest
import inspyred
import pickle

from pandas.util.testing import assert_frame_equal

from cameo import load_model, fba, config
from cameo.util import RandomGenerator as Random
from cameo.strain_design.heuristic.optimization import HeuristicOptimization, ReactionKnockoutOptimization, \
    set_distance_function
from cameo.strain_design.heuristic.archivers import SolutionTuple, BestSolutionArchiver
from cameo.strain_design.heuristic.decoders import ReactionKnockoutDecoder, KnockoutDecoder, GeneKnockoutDecoder
from cameo.strain_design.heuristic.generators import set_generator, unique_set_generator
from cameo.strain_design.heuristic.objective_functions import biomass_product_coupled_yield, product_yield, \
    number_of_knockouts
from cobra.manipulation.delete import find_gene_knockout_reactions
from cameo.parallel import SequentialView, MultiprocessingView

config.default_view = MultiprocessingView(processes=2)

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


class TestBestSolutionArchiver(unittest.TestCase):
    def test_solution_string(self):
        sol1 = SolutionTuple(SOLUTIONS[0][0], SOLUTIONS[0][1])
        sol2 = SolutionTuple(SOLUTIONS[1][0], SOLUTIONS[1][1])
        sol3 = SolutionTuple(SOLUTIONS[2][0], SOLUTIONS[2][1])
        self.assertEqual(sol1.__str__(), "[1, 2, 3] - 0.1")
        self.assertEqual(sol2.__str__(), "[1, 2, 3, 4] - 0.1")
        self.assertEqual(sol3.__str__(), "[2, 3, 4] - 0.45")

    def test_solution_comparison(self):
        sol1 = SolutionTuple(SOLUTIONS[0][0], SOLUTIONS[0][1])
        sol2 = SolutionTuple(SOLUTIONS[1][0], SOLUTIONS[1][1])
        sol3 = SolutionTuple(SOLUTIONS[2][0], SOLUTIONS[2][1])

        #test ordering
        self.assertEqual(sol1.__cmp__(sol2), -1)
        self.assertEqual(sol1.__cmp__(sol1), 0)
        self.assertEqual(sol1.__cmp__(sol3), 1)

        #testing issubset
        self.assertTrue(sol1.issubset(sol2), msg="Solution 1 is subset of Solution 2")
        self.assertFalse(sol2.issubset(sol1), msg="Solution 2 is not subset of Solution 1")
        self.assertTrue(sol3.issubset(sol2), msg="Solution 3 is subset of Solution 2")
        self.assertFalse(sol2.issubset(sol3), msg="Solution 2 is not subset of Solution 3")
        self.assertFalse(sol1.issubset(sol3), msg="Solution 1 is subset of Solution 3")
        self.assertFalse(sol2.issubset(sol3), msg="Solution 3 is not subset of Solution 1")

        #test difference
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

    def test_add_greater_solution_with_same_fitness(self):
        size = 1
        pool = BestSolutionArchiver()
        pool.add(SOLUTIONS[0][0], SOLUTIONS[0][1], size)
        pool.add(SOLUTIONS[1][0], SOLUTIONS[1][1], size)
        self.assertEqual(pool.length(), 1, msg="Pool must keep one solution (length=%s)" % pool.length())
        best_solution = set(SOLUTIONS[0][0])
        best_fitness = SOLUTIONS[0][1]
        sol = pool.get(0)
        self.assertEqual(sol.candidate, best_solution, msg="Best solution set must be the first")
        self.assertEqual(sol.fitness, best_fitness, msg="Best solution fitness must be the first")

    def test_add_smaller_solution_with_same_fitness(self):
        size = 1
        pool = BestSolutionArchiver()
        pool.add(SOLUTIONS[1][0], SOLUTIONS[1][1], size)
        pool.add(SOLUTIONS[0][0], SOLUTIONS[0][1], size)
        self.assertEqual(pool.length(), 1, msg="Pool must keep one solution (length=%s)" % pool.length())
        solution = set(SOLUTIONS[0][0])
        fitness = SOLUTIONS[0][1]
        sol = pool.get(0)
        self.assertEqual(sol.candidate, solution, msg="Best solution must be the first (%s)" % sol.candidate)
        self.assertEqual(sol.fitness, fitness, msg="Best fitness must be the first (%s)" % sol.fitness)

    def test_pool_size_limit(self):
        size = 1
        pool = BestSolutionArchiver()
        pool.add(SOLUTIONS[0][0], SOLUTIONS[0][1], size)
        pool.add(SOLUTIONS[1][0], SOLUTIONS[1][1], size)
        pool.add(SOLUTIONS[2][0], SOLUTIONS[2][1], size)
        pool.add(SOLUTIONS[3][0], SOLUTIONS[3][1], size)
        pool.add(SOLUTIONS[4][0], SOLUTIONS[4][1], size)
        pool.add(SOLUTIONS[5][0], SOLUTIONS[5][1], size)
        pool.add(SOLUTIONS[6][0], SOLUTIONS[6][1], size)
        pool.add(SOLUTIONS[7][0], SOLUTIONS[7][1], size)
        pool.add(SOLUTIONS[8][0], SOLUTIONS[8][1], size)
        pool.add(SOLUTIONS[9][0], SOLUTIONS[9][1], size)
        self.assertLessEqual(pool.length(), 1, msg="Pool must keep one solution (length=%s)" % pool.length())
        size = 2
        pool = BestSolutionArchiver()
        pool.add(SOLUTIONS[0][0], SOLUTIONS[0][1], size)
        pool.add(SOLUTIONS[1][0], SOLUTIONS[1][1], size)
        pool.add(SOLUTIONS[2][0], SOLUTIONS[2][1], size)
        pool.add(SOLUTIONS[3][0], SOLUTIONS[3][1], size)
        pool.add(SOLUTIONS[4][0], SOLUTIONS[4][1], size)
        pool.add(SOLUTIONS[5][0], SOLUTIONS[5][1], size)
        pool.add(SOLUTIONS[6][0], SOLUTIONS[6][1], size)
        pool.add(SOLUTIONS[7][0], SOLUTIONS[7][1], size)
        pool.add(SOLUTIONS[8][0], SOLUTIONS[8][1], size)
        pool.add(SOLUTIONS[9][0], SOLUTIONS[9][1], size)
        self.assertLessEqual(pool.length(), 2, msg="Pool must keep one solution (length=%s)" % pool.length())
        size = 3
        pool = BestSolutionArchiver()
        pool.add(SOLUTIONS[0][0], SOLUTIONS[0][1], size)
        pool.add(SOLUTIONS[1][0], SOLUTIONS[1][1], size)
        pool.add(SOLUTIONS[2][0], SOLUTIONS[2][1], size)
        pool.add(SOLUTIONS[3][0], SOLUTIONS[3][1], size)
        pool.add(SOLUTIONS[4][0], SOLUTIONS[4][1], size)
        pool.add(SOLUTIONS[5][0], SOLUTIONS[5][1], size)
        pool.add(SOLUTIONS[6][0], SOLUTIONS[6][1], size)
        pool.add(SOLUTIONS[7][0], SOLUTIONS[7][1], size)
        pool.add(SOLUTIONS[8][0], SOLUTIONS[8][1], size)
        pool.add(SOLUTIONS[9][0], SOLUTIONS[9][1], size)
        self.assertLessEqual(pool.length(), 3, msg="Pool must keep one solution (length=%s)" % pool.length())
        size = 4
        pool = BestSolutionArchiver()
        pool.add(SOLUTIONS[0][0], SOLUTIONS[0][1], size)
        pool.add(SOLUTIONS[1][0], SOLUTIONS[1][1], size)
        pool.add(SOLUTIONS[2][0], SOLUTIONS[2][1], size)
        pool.add(SOLUTIONS[3][0], SOLUTIONS[3][1], size)
        pool.add(SOLUTIONS[4][0], SOLUTIONS[4][1], size)
        pool.add(SOLUTIONS[5][0], SOLUTIONS[5][1], size)
        pool.add(SOLUTIONS[6][0], SOLUTIONS[6][1], size)
        pool.add(SOLUTIONS[7][0], SOLUTIONS[7][1], size)
        pool.add(SOLUTIONS[8][0], SOLUTIONS[8][1], size)
        pool.add(SOLUTIONS[9][0], SOLUTIONS[9][1], size)
        self.assertLessEqual(pool.length(), 4, msg="Pool must keep one solution (length=%s)" % pool.length())

    def test_callable_pool(self):
        pool = BestSolutionArchiver()
        size = 3
        args = {}
        args.setdefault('max_archive_size', size)
        population = [SolutionTuple(SOLUTIONS[0][0], SOLUTIONS[0][1]),
                      SolutionTuple(SOLUTIONS[1][0], SOLUTIONS[1][1]),
                      SolutionTuple(SOLUTIONS[2][0], SOLUTIONS[2][1]),
                      SolutionTuple(SOLUTIONS[3][0], SOLUTIONS[3][1]),
                      SolutionTuple(SOLUTIONS[4][0], SOLUTIONS[4][1]),
                      SolutionTuple(SOLUTIONS[5][0], SOLUTIONS[5][1]),
                      SolutionTuple(SOLUTIONS[6][0], SOLUTIONS[6][1])]
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

        fitness = of(None, solution, None)
        self.assertAlmostEqual((0.6 * 2)/10, fitness)

        solution.set_primal('substrate', 0)

        fitness = of(None, solution, None)
        self.assertEquals(0, fitness)

    def test_yield(self):
        solution = self._MockupSolution()
        solution.set_primal('biomass', 0.6)
        solution.set_primal('product', 2)
        solution.set_primal('substrate', -10)

        of = product_yield("product", "substrate")
        fitness = of(None, solution, None)
        self.assertAlmostEqual(2.0/10.0, fitness)

        solution.set_primal('substrate', 0)
        fitness = of(None, solution, None)
        self.assertEquals(0, fitness)

    def test_number_of_knockouts(self):

        of_max = number_of_knockouts(sense='max')
        of_min = number_of_knockouts(sense='min')

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
        self.assertItemsEqual(reactions1, reactions2)

    def test_gene_knockout_decoder(self):
        decoder = GeneKnockoutDecoder([g.id for g in self.model.genes], self.model)
        reactions1, genes = decoder([1, 2, 3, 4])
        reactions2 = find_gene_knockout_reactions(self.model, genes)
        self.assertItemsEqual(reactions1, reactions2)


class TestGeneratos(unittest.TestCase):
    def setUp(self):
        self.model = TEST_MODEL
        self.args = {}
        self.args.setdefault('representation', [r.id for r in self.model.reactions])
        self.random = Random()

    def test_fixed_size_generator(self):
        self.args.setdefault('variable_candidate_size', False)

        self.args['candidate_size'] = 10
        for _ in xrange(10000):
            candidate = set_generator(self.random, self.args)
            self.assertEqual(len(candidate), 10)
            candidate = unique_set_generator(self.random, self.args)
            self.assertEqual(len(candidate), 10)

        self.args['candidate_size'] = 20
        for _ in xrange(10000):
            candidate = set_generator(self.random, self.args)
            self.assertEqual(len(candidate), 20)
            candidate = unique_set_generator(self.random, self.args)
            self.assertEqual(len(candidate), 20)

    def test_variable_size_generator(self):
        self.args.setdefault('variable_candidate_size', True)

        self.args['candidate_size'] = 10
        for _ in xrange(10000):
            candidate = set_generator(self.random, self.args)
            self.assertLessEqual(len(candidate), 10)
            candidate = unique_set_generator(self.random, self.args)
            self.assertLessEqual(len(candidate), 10)

        self.args['candidate_size'] = 20
        for _ in xrange(10000):
            candidate = set_generator(self.random, self.args)
            self.assertLessEqual(len(candidate), 20)
            candidate = unique_set_generator(self.random, self.args)
            self.assertLessEqual(len(candidate), 20)


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


class TestReactionKnockoutOptimization(unittest.TestCase):
    def setUp(self):
        self.model = TEST_MODEL
        self.essential_reactions = set([r.id for r in self.model.essential_reactions()])

    def test_initialize(self):
        rko = ReactionKnockoutOptimization(model=self.model,
                                           simulation_method=fba,
                                           seed=SEED)

        self.assertItemsEqual(self.essential_reactions, rko.essential_reactions)
        self.assertEqual(rko._ko_type, "reaction")
        self.assertTrue(isinstance(rko._decoder, ReactionKnockoutDecoder))

    @unittest.skip('Not deterministic when seeded')
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

        self.assertEqual(rko.random.random(), 0.1915194503788923)

        results = rko.run(max_evaluations=3000, pop_size=10, view=SequentialView())

        self.assertEqual(rko.random.random(), 0.04225378600400298)

        # with open(result_file, 'w') as f:
        #     pickle.dump(results, f)

        with open(result_file, 'r') as f:
            expected_results = pickle.load(f)

        assert_frame_equal(results.solutions, expected_results.solutions)

    @unittest.skip('Not deterministic when seeded')
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

        results = rko.run(max_evaluations=3000, pop_size=10, view=SequentialView())

        with open(result_file, 'w') as file:
            pickle.dump(results, file)

        with open(result_file, 'r') as file:
            expected_results = pickle.load(file)


        assert_frame_equal(results.solutions, expected_results.solutions)

    def test_evaluator(self):
        pass
