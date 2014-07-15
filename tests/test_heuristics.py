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
import os

import unittest
from cameo import load_model
from cameo.strain_design.heuristic.archivers import SolutionTuple, BestSolutionArchiver
from cameo.strain_design.heuristic.objective_functions import biomass_product_coupled_yield, product_yield, \
    number_of_knockouts
from cobra.manipulation.delete import find_gene_knockout_reactions
from cameo.util import TimeMachine

iJO1366_PATH = os.path.join(os.path.dirname(__file__), "data/iJO1366.xml")

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
    def test_solution_comparison(self):
        sol1 = SolutionTuple(SOLUTIONS[0][0], SOLUTIONS[0][1])
        sol2 = SolutionTuple(SOLUTIONS[1][0], SOLUTIONS[1][1])
        sol3 = SolutionTuple(SOLUTIONS[2][0], SOLUTIONS[2][1])

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


class TestObjectiveFunctions(unittest.TestCase):
    lvaline_sol = ["b0115", "b3236", "b3916"] #aceF, mdh, pfkA from Park et al. 2007 (PNAS)

    def test_biomass_product_coupled_yield(self):
        model = load_model(iJO1366_PATH)
        of = biomass_product_coupled_yield("Ec_biomass_iJO1366_core_53p95M", "EX_succ_LPAREN_e_RPAREN_",
                                           "EX_glc_LPAREN_e_RPAREN_")
        reactions = find_gene_knockout_reactions(model, [model.genes.get_by_id(gene) for gene in self.lvaline_sol])
        tm = TimeMachine()
        for reaction in reactions:
            tm(do=partial(setattr, reaction, 'lower_bound', 0),
               undo=partial(setattr, reaction, 'lower_bound', reaction.lower_bound))
            tm(do=partial(setattr, reaction, 'upper_bound', 0),
               undo=partial(setattr, reaction, 'upper_bound', reaction.upper_bound))

        solution = model.optimize()
        tm.reset()

    def test_yield(self):
        model = load_model(iJO1366_PATH)
        of = product_yield("EX_succ_LPAREN_e_RPAREN_", "EX_glc_LPAREN_e_RPAREN_")
        reactions = find_gene_knockout_reactions(model, [model.genes.get_by_id(gene) for gene in self.lvaline_sol])
        tm = TimeMachine()
        for reaction in reactions:
            tm(do=partial(setattr, reaction, 'lower_bound', 0),
               undo=partial(setattr, reaction, 'lower_bound', reaction.lower_bound))
            tm(do=partial(setattr, reaction, 'upper_bound', 0),
               undo=partial(setattr, reaction, 'upper_bound', reaction.upper_bound))

        solution = model.optimize()
        tm.reset()

    def test_number_of_knockouts(self):
        model = load_model(iJO1366_PATH)
        of = number_of_knockouts(sense='min')
        reactions = find_gene_knockout_reactions(model, [model.genes.get_by_id(gene) for gene in self.lvaline_sol])
        tm = TimeMachine()
        for reaction in reactions:
            tm(do=partial(setattr, reaction, 'lower_bound', 0),
               undo=partial(setattr, reaction, 'lower_bound', reaction.lower_bound))
            tm(do=partial(setattr, reaction, 'upper_bound', 0),
               undo=partial(setattr, reaction, 'upper_bound', reaction.upper_bound))

        solution = model.optimize()
        tm.reset()
