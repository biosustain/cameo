# @formatter:off
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
# @formatter:on

from __future__ import absolute_import, print_function

import os
import copy
import unittest

from cobra import Metabolite
import numpy
from optlang import Objective
from cobra.io import read_sbml_model
import optlang
import pandas

from cameo import load_model, Reaction, Model
from cameo.config import solvers
from cameo.exceptions import UndefinedSolution
from cameo.core.solver_based_model import Reaction
from cameo.util import TimeMachine
import six

TRAVIS = os.getenv('TRAVIS', False)
TESTDIR = os.path.dirname(__file__)
REFERENCE_FVA_SOLUTION_ECOLI_CORE = pandas.read_csv(os.path.join(TESTDIR, 'data/REFERENCE_flux_ranges_EcoliCore.csv'),
                                                    index_col=0)
TESTMODEL = load_model(os.path.join(TESTDIR, 'data/EcoliCore.xml'), sanitize=False)
COBRAPYTESTMODEL = read_sbml_model(os.path.join(TESTDIR, 'data/EcoliCore.xml'))
ESSENTIAL_GENES = ['b2779', 'b1779', 'b0720', 'b0451', 'b2416', 'b2926', 'b1136', 'b2415']
ESSENTIAL_REACTIONS = ['GLNS', 'Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2', 'PIt2r', 'GAPD', 'ACONTb',
                       'EX_nh4_LPAREN_e_RPAREN_', 'ENO', 'EX_h_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_', 'ICDHyr',
                       'CS', 'NH4t', 'GLCpts', 'PGM', 'EX_pi_LPAREN_e_RPAREN_', 'PGK', 'RPI', 'ACONTa']


class CommonGround(object):
    def setUp(self):
        self.model = TESTMODEL.copy()
        self.model.optimize()


class AbstractTestLazySolution(CommonGround):
    def setUp(self):
        super(AbstractTestLazySolution, self).setUp()
        self.solution = self.model.optimize()

    def test_self_invalidation(self):
        solution = self.model.solve()
        self.assertAlmostEqual(solution.f, 0.873921506968431, delta=0.000001)
        self.model.optimize()
        self.assertRaises(UndefinedSolution, getattr, solution, 'f')

    def test_solution_contains_only_reaction_specific_values(self):
        reaction_IDs = set([reaction.id for reaction in self.model.reactions])
        self.assertEqual(set(self.solution.x_dict.keys()).difference(reaction_IDs), set())
        self.assertEqual(set(self.solution.y_dict.keys()).difference(reaction_IDs), set())
        self.assertEqual(set(self.solution.reduced_costs.keys()).difference(reaction_IDs), set())
        metabolite_IDs = set([metabolite.id for metabolite in self.model.metabolites])
        self.assertEqual(set(self.solution.shadow_prices.keys()).difference(metabolite_IDs), set())


class TestLazySolutionGLPK(AbstractTestLazySolution, unittest.TestCase):
    def setUp(self):
        super(TestLazySolutionGLPK, self).setUp()
        self.model.solver = 'glpk'
        self.solution = self.model.optimize()


@unittest.skipIf(TRAVIS, 'CPLEX not available on Travis.')
class TestLazySolutionCPLEX(AbstractTestLazySolution, unittest.TestCase):
    def setUp(self):
        super(TestLazySolutionCPLEX, self).setUp()
        self.model.solver = 'glpk'
        self.solution = self.model.optimize()


class AbstractTestReaction(object):
    def test_clone_cobrapy_reaction(self):
        for reaction in self.cobrapy_model.reactions:
            cloned_reaction = Reaction.clone(reaction)
            self.assertEqual(cloned_reaction.objective_coefficient, reaction.objective_coefficient)
            self.assertEqual(cloned_reaction.gene_reaction_rule, reaction.gene_reaction_rule)
            self.assertEqual(cloned_reaction.genes, reaction.genes)
            self.assertEqual(cloned_reaction.metabolites, reaction.metabolites)
            self.assertEqual(cloned_reaction.products, reaction.products)
            self.assertEqual(cloned_reaction.reactants, reaction.reactants)

    def test_str(self):
        self.assertTrue(self.model.reactions[0].__str__().startswith('ACALD'))

    def test_add_metabolite(self):
        model = self.model
        pgi_reaction = model.reactions.PGI
        test_met = model.metabolites[0]
        pgi_reaction.add_metabolites({test_met: 42}, combine=False)
        self.assertEqual(pgi_reaction.metabolites[test_met], 42)
        self.assertEqual(
            model.solver.constraints[test_met.id].expression.as_coefficients_dict()[pgi_reaction.forward_variable], 42)
        self.assertEqual(
            model.solver.constraints[test_met.id].expression.as_coefficients_dict()[pgi_reaction.reverse_variable], -42)

        pgi_reaction.add_metabolites({test_met: -10}, combine=True)
        self.assertEqual(pgi_reaction.metabolites[test_met], 32)
        self.assertEqual(
            model.solver.constraints[test_met.id].expression.as_coefficients_dict()[pgi_reaction.forward_variable], 32)
        self.assertEqual(
            model.solver.constraints[test_met.id].expression.as_coefficients_dict()[pgi_reaction.reverse_variable], -32)

        pgi_reaction.add_metabolites({test_met: 0}, combine=False)
        with self.assertRaises(KeyError):
            pgi_reaction.metabolites[test_met]
        self.assertEqual(
            model.solver.constraints[test_met.id].expression.as_coefficients_dict()[pgi_reaction.forward_variable], 0)
        self.assertEqual(
            model.solver.constraints[test_met.id].expression.as_coefficients_dict()[pgi_reaction.reverse_variable], 0)

        # test_met_2 = Metabolite("Test2", compartment="c")
        # pgi_reaction.add_metabolites({test_met_2: 43}, combine=False)
        # self.assertEqual(pgi_reaction.metabolites[test_met], 43)
        # self.assertEqual(
        #    model.solver.constraints[test_met.id].expression.as_coefficients_dict()[pgi_reaction.forward_variable], 43)
        # self.assertEqual(
        #    model.solver.constraints[test_met.id].expression.as_coefficients_dict()[pgi_reaction.reverse_variable], -43)

    def test_removal_from_model_retains_bounds(self):
        model_cp = self.model.copy()
        reaction = model_cp.reactions.ACALD
        self.assertEqual(reaction.model, model_cp)
        self.assertEqual(reaction.lower_bound, -999999.0)
        self.assertEqual(reaction.upper_bound, 999999.0)
        self.assertEqual(reaction._lower_bound, -999999.0)
        self.assertEqual(reaction._upper_bound, 999999.0)
        model_cp.remove_reactions([reaction])
        self.assertEqual(reaction.model, None)
        self.assertEqual(reaction.lower_bound, -999999.0)
        self.assertEqual(reaction.upper_bound, 999999.0)
        self.assertEqual(reaction._lower_bound, -999999.0)
        self.assertEqual(reaction._upper_bound, 999999.0)

    def test_set_bounds_scenario_1(self):
        model = self.model
        acald_reaction = model.reactions.ACALD
        self.assertEqual(acald_reaction.lower_bound, -999999.)
        self.assertEqual(acald_reaction.upper_bound, 999999.)
        self.assertEqual(acald_reaction.variable.lb, 0.)
        self.assertEqual(acald_reaction.variable.ub, 999999.)
        self.assertEqual(acald_reaction.reverse_variable.lb, 0)
        self.assertEqual(acald_reaction.reverse_variable.ub, 999999.)
        acald_reaction.upper_bound = acald_reaction.lower_bound - 100
        self.assertEqual(acald_reaction.lower_bound, -1000099.0)
        self.assertEqual(acald_reaction.upper_bound, -1000099.0)
        self.assertEqual(acald_reaction.variable.lb, 0)
        self.assertEqual(acald_reaction.variable.ub, 0)
        self.assertEqual(acald_reaction.reverse_variable.lb, 1000099.)
        self.assertEqual(acald_reaction.reverse_variable.ub, 1000099.)
        acald_reaction.upper_bound = 100
        self.assertEqual(acald_reaction.lower_bound, -1000099.0)
        self.assertEqual(acald_reaction.upper_bound, 100)
        self.assertEqual(acald_reaction.variable.lb, 0)
        self.assertEqual(acald_reaction.variable.ub, 100)
        self.assertEqual(acald_reaction.reverse_variable.lb, 0)
        self.assertEqual(acald_reaction.reverse_variable.ub, 1000099.0)

    def test_set_upper_before_lower_bound_to_0(self):
        model = self.model
        model.reactions.GAPD.upper_bound = 0
        model.reactions.GAPD.lower_bound = 0
        self.assertEqual(model.reactions.GAPD.lower_bound, 0)
        self.assertEqual(model.reactions.GAPD.upper_bound, 0)
        self.assertEqual(model.reactions.GAPD.variable.lb, 0)
        self.assertEqual(model.reactions.GAPD.variable.ub, 0)
        self.assertEqual(model.reactions.GAPD.reverse_variable.lb, 0)
        self.assertEqual(model.reactions.GAPD.reverse_variable.ub, 0)

    def test_set_bounds_scenario_2(self):
        model = self.model
        acald_reaction = model.reactions.ACALD
        self.assertEqual(acald_reaction.lower_bound, -999999.)
        self.assertEqual(acald_reaction.upper_bound, 999999.)
        self.assertEqual(acald_reaction.variable.lb, 0.)
        self.assertEqual(acald_reaction.variable.ub, 999999.)
        self.assertEqual(acald_reaction.reverse_variable.lb, 0)
        self.assertEqual(acald_reaction.reverse_variable.ub, 999999.)
        acald_reaction.lower_bound = acald_reaction.upper_bound + 100
        self.assertEqual(acald_reaction.lower_bound, 1000099.0)
        self.assertEqual(acald_reaction.upper_bound, 1000099.0)
        self.assertEqual(acald_reaction.variable.lb, 1000099.0)
        self.assertEqual(acald_reaction.variable.ub, 1000099.0)
        self.assertEqual(acald_reaction.reverse_variable.lb, 0)
        self.assertEqual(acald_reaction.reverse_variable.ub, 0)
        acald_reaction.lower_bound = -100
        self.assertEqual(acald_reaction.lower_bound, -100.)
        self.assertEqual(acald_reaction.upper_bound, 1000099.)
        self.assertEqual(acald_reaction.variable.lb, 0)
        self.assertEqual(acald_reaction.variable.ub, 1000099.)
        self.assertEqual(acald_reaction.reverse_variable.lb, 0)
        self.assertEqual(acald_reaction.reverse_variable.ub, 100)

    def test_make_irreversible(self):
        model = self.model
        acald_reaction = model.reactions.ACALD
        self.assertEqual(acald_reaction.lower_bound, -999999.)
        self.assertEqual(acald_reaction.upper_bound, 999999.)
        self.assertEqual(acald_reaction.variable.lb, 0.)
        self.assertEqual(acald_reaction.variable.ub, 999999.)
        self.assertEqual(acald_reaction.reverse_variable.lb, 0)
        self.assertEqual(acald_reaction.reverse_variable.ub, 999999.)
        acald_reaction.lower_bound = 0
        self.assertEqual(acald_reaction.lower_bound, 0)
        self.assertEqual(acald_reaction.upper_bound, 999999.)
        self.assertEqual(acald_reaction.variable.lb, 0)
        self.assertEqual(acald_reaction.variable.ub, 999999.0)
        self.assertEqual(acald_reaction.reverse_variable.lb, 0)
        self.assertEqual(acald_reaction.reverse_variable.ub, 0)
        acald_reaction.lower_bound = -100
        self.assertEqual(acald_reaction.lower_bound, -100.)
        self.assertEqual(acald_reaction.upper_bound, 999999.)
        self.assertEqual(acald_reaction.variable.lb, 0)
        self.assertEqual(acald_reaction.variable.ub, 999999.)
        self.assertEqual(acald_reaction.reverse_variable.lb, 0)
        self.assertEqual(acald_reaction.reverse_variable.ub, 100)

    def test_make_reversible(self):
        model = self.model
        pfk_reaction = model.reactions.PFK
        self.assertEqual(pfk_reaction.lower_bound, 0.)
        self.assertEqual(pfk_reaction.upper_bound, 999999.)
        self.assertEqual(pfk_reaction.variable.lb, 0.)
        self.assertEqual(pfk_reaction.variable.ub, 999999.)
        self.assertEqual(pfk_reaction.reverse_variable.lb, 0)
        self.assertEqual(pfk_reaction.reverse_variable.ub, 0)
        pfk_reaction.lower_bound = -100.
        self.assertEqual(pfk_reaction.lower_bound, -100.)
        self.assertEqual(pfk_reaction.upper_bound, 999999.)
        self.assertEqual(pfk_reaction.variable.lb, 0)
        self.assertEqual(pfk_reaction.variable.ub, 999999.0)
        self.assertEqual(pfk_reaction.reverse_variable.lb, 0)
        self.assertEqual(pfk_reaction.reverse_variable.ub, 100.)
        pfk_reaction.lower_bound = 0
        self.assertEqual(pfk_reaction.lower_bound, 0)
        self.assertEqual(pfk_reaction.upper_bound, 999999.)
        self.assertEqual(pfk_reaction.variable.lb, 0)
        self.assertEqual(pfk_reaction.variable.ub, 999999.)
        self.assertEqual(pfk_reaction.reverse_variable.lb, 0)
        self.assertEqual(pfk_reaction.reverse_variable.ub, 0)

    def test_make_irreversible_irreversible_to_the_other_side(self):
        model = self.model
        pfk_reaction = model.reactions.PFK
        self.assertEqual(pfk_reaction.lower_bound, 0.)
        self.assertEqual(pfk_reaction.upper_bound, 999999.)
        self.assertEqual(pfk_reaction.variable.lb, 0.)
        self.assertEqual(pfk_reaction.variable.ub, 999999.)
        self.assertEqual(pfk_reaction.reverse_variable.lb, 0)
        self.assertEqual(pfk_reaction.reverse_variable.ub, 0)
        pfk_reaction.upper_bound = -100.
        self.assertEqual(pfk_reaction.forward_variable.lb, 0)
        self.assertEqual(pfk_reaction.forward_variable.ub, 0)
        self.assertEqual(pfk_reaction.reverse_variable.lb, 100)
        self.assertEqual(pfk_reaction.reverse_variable.ub, 100)
        pfk_reaction.lower_bound = -999999.
        self.assertEqual(pfk_reaction.lower_bound, -999999.)
        self.assertEqual(pfk_reaction.upper_bound, -100.)
        self.assertEqual(pfk_reaction.variable.lb, 0)
        self.assertEqual(pfk_reaction.variable.ub, 0)
        self.assertEqual(pfk_reaction.reverse_variable.lb, 100)
        self.assertEqual(pfk_reaction.reverse_variable.ub, 999999.)
        # self.assertEqual(pfk_reaction.reverse_variable.lb, 0.)
        # self.assertEqual(pfk_reaction.reverse_variable.ub, 0.)

    def test_make_lhs_irreversible_reversible(self):
        model = self.model
        rxn = Reaction('test')
        rxn.add_metabolites({model.metabolites[0]: -1., model.metabolites[1]: 1.})
        rxn.lower_bound = -999999.
        rxn.upper_bound = -100
        model.add_reaction(rxn)
        self.assertEqual(rxn.lower_bound, -999999.)
        self.assertEqual(rxn.upper_bound, -100.)
        self.assertEqual(rxn.variable.lb, 0.)
        self.assertEqual(rxn.variable.ub, 0.)
        self.assertEqual(rxn.reverse_variable.lb, 100.)
        self.assertEqual(rxn.reverse_variable.ub, 999999.)
        rxn.upper_bound = 666.
        self.assertEqual(rxn.lower_bound, -999999.)
        self.assertEqual(rxn.upper_bound, 666.)
        self.assertEqual(rxn.variable.lb, 0.)
        self.assertEqual(rxn.variable.ub, 666)
        self.assertEqual(rxn.reverse_variable.lb, 0.)
        self.assertEqual(rxn.reverse_variable.ub, 999999.)

    def test_model_less_reaction(self):
        # self.model.solver.configuration.verbosity = 3
        self.model.solve()
        print(self.model.reactions.ACALD.flux)
        for reaction in self.model.reactions:
            self.assertTrue(isinstance(reaction.flux, float))
            self.assertTrue(isinstance(reaction.reduced_cost, float))
        for reaction in self.model.reactions:
            self.model.remove_reactions([reaction])
            self.assertEqual(reaction.flux, None)
            self.assertEqual(reaction.reduced_cost, None)

    def test_knockout(self):
        original_bounds = dict()
        for reaction in self.model.reactions:
            original_bounds[reaction.id] = (reaction.lower_bound, reaction.upper_bound)
            reaction.knock_out()
            self.assertEqual(reaction.lower_bound, 0)
            self.assertEqual(reaction.upper_bound, 0)
        for k, (lb, ub) in six.iteritems(original_bounds):
            self.model.reactions.get_by_id(k).lower_bound = lb
            self.model.reactions.get_by_id(k).upper_bound = ub
        for reaction in self.model.reactions:
            self.assertEqual(reaction.lower_bound, original_bounds[reaction.id][0])
            self.assertEqual(reaction.upper_bound, original_bounds[reaction.id][1])
        with TimeMachine() as tm:
            for reaction in self.model.reactions:
                original_bounds[reaction.id] = (reaction.lower_bound, reaction.upper_bound)
                reaction.knock_out(time_machine=tm)
                self.assertEqual(reaction.lower_bound, 0)
                self.assertEqual(reaction.upper_bound, 0)
        for reaction in self.model.reactions:
            self.assertEqual(reaction.lower_bound, original_bounds[reaction.id][0])
            self.assertEqual(reaction.upper_bound, original_bounds[reaction.id][1])

    def test_weird_left_to_right_reaction_issue(self):

        model = Model("Toy Model")

        m1 = Metabolite("M1")
        d1 = Reaction("ex1")
        d1.add_metabolites({m1: -1})
        d1.upper_bound = 0
        d1.lower_bound = -1000
        # print d1.reaction, d1.lower_bound, d1.upper_bound
        model.add_reactions([d1])
        self.assertFalse(d1.reversibility)
        self.assertEqual(d1.lower_bound, -1000)
        self.assertEqual(d1._lower_bound, -1000)
        self.assertEqual(d1.upper_bound, 0)
        self.assertEqual(d1._upper_bound, 0)
        with TimeMachine() as tm:
            d1.knock_out(time_machine=tm)
            self.assertEqual(d1.lower_bound, 0)
            self.assertEqual(d1._lower_bound, 0)
            self.assertEqual(d1.upper_bound, 0)
            self.assertEqual(d1._upper_bound, 0)
        self.assertEqual(d1.lower_bound, -1000)
        self.assertEqual(d1._lower_bound, -1000)
        self.assertEqual(d1.upper_bound, 0)
        self.assertEqual(d1._upper_bound, 0)

    def test_one_left_to_right_reaction_set_positive_ub(self):

        model = Model("Toy Model")

        m1 = Metabolite("M1")
        d1 = Reaction("ex1")
        d1.add_metabolites({m1: -1})
        d1.upper_bound = 0
        d1.lower_bound = -1000
        model.add_reactions([d1])
        self.assertEqual(d1.reverse_variable.lb, 0)
        self.assertEqual(d1.reverse_variable.ub, 1000)
        self.assertEqual(d1._lower_bound, -1000)
        self.assertEqual(d1.lower_bound, -1000)
        self.assertEqual(d1._upper_bound, 0)
        self.assertEqual(d1.upper_bound, 0)
        self.assertEqual(d1.variable.lb, 0)
        self.assertEqual(d1.variable.ub, 0)
        d1.upper_bound = .1
        self.assertEqual(d1.variable.lb, 0)
        self.assertEqual(d1.variable.ub, .1)
        self.assertEqual(d1.reverse_variable.lb, 0)
        self.assertEqual(d1.reverse_variable.ub, 1000)
        self.assertEqual(d1._lower_bound, -1000)
        self.assertEqual(d1.upper_bound, .1)
        self.assertEqual(d1._lower_bound, -1000)
        self.assertEqual(d1.upper_bound, .1)

    def test_irrev_reaction_set_negative_lb(self):
        self.assertFalse(self.model.reactions.PFK.reversibility)
        self.assertEqual(self.model.reactions.PFK.lower_bound, 0)
        self.assertEqual(self.model.reactions.PFK.upper_bound, 999999.0)
        self.assertEqual(self.model.reactions.PFK.variable.lb, 0)
        self.assertEqual(self.model.reactions.PFK.variable.ub, 999999.0)
        self.assertEqual(self.model.reactions.PFK.reverse_variable.lb, 0)
        self.assertEqual(self.model.reactions.PFK.reverse_variable.ub, 0)
        self.model.reactions.PFK.lower_bound = -1000000000
        self.assertEqual(self.model.reactions.PFK.lower_bound, -1000000000)
        self.assertEqual(self.model.reactions.PFK.upper_bound, 999999.0)
        self.assertEqual(self.model.reactions.PFK.variable.lb, 0)
        self.assertEqual(self.model.reactions.PFK.variable.ub, 999999.0)
        self.assertEqual(self.model.reactions.PFK.reverse_variable.lb, 0)
        self.assertEqual(self.model.reactions.PFK.reverse_variable.ub, 1000000000)

    def test_twist_irrev_right_to_left_reaction_to_left_to_right(self):
        self.assertFalse(self.model.reactions.PFK.reversibility)
        self.assertEqual(self.model.reactions.PFK.lower_bound, 0)
        self.assertEqual(self.model.reactions.PFK.upper_bound, 999999.0)
        self.assertEqual(self.model.reactions.PFK.variable.lb, 0)
        self.assertEqual(self.model.reactions.PFK.variable.ub, 999999.0)
        self.assertEqual(self.model.reactions.PFK.reverse_variable.lb, 0)
        self.assertEqual(self.model.reactions.PFK.reverse_variable.ub, 0)
        self.model.reactions.PFK.lower_bound = -1000000000
        self.model.reactions.PFK.upper_bound = 0
        self.assertEqual(self.model.reactions.PFK.lower_bound, -1000000000)
        self.assertEqual(self.model.reactions.PFK.upper_bound, 0)
        self.assertEqual(self.model.reactions.PFK.variable.lb, 0)
        self.assertEqual(self.model.reactions.PFK.variable.ub, 0)
        self.assertEqual(self.model.reactions.PFK.reverse_variable.lb, 0)
        self.assertEqual(self.model.reactions.PFK.reverse_variable.ub, 1000000000)

    def test_iMM904_4HGLSDm_problem(self):
        model = load_model(os.path.join(TESTDIR, 'data/iMM904.xml'))
        # set upper bound before lower bound after knockout
        cp = model.copy()
        rxn = cp.reactions.get_by_id('4HGLSDm')
        prev_lb, prev_ub = rxn.lower_bound, rxn.upper_bound
        rxn.lower_bound = 0
        rxn.upper_bound = 0
        rxn.upper_bound = prev_ub
        rxn.lower_bound = prev_lb
        self.assertEquals(rxn.lower_bound, prev_lb)
        self.assertEquals(rxn.upper_bound, prev_ub)
        # set lower bound before upper bound after knockout
        cp = model.copy()
        rxn = cp.reactions.get_by_id('4HGLSDm')
        prev_lb, prev_ub = rxn.lower_bound, rxn.upper_bound
        rxn.lower_bound = 0
        rxn.upper_bound = 0
        rxn.lower_bound = prev_lb
        rxn.upper_bound = prev_ub
        self.assertEquals(rxn.lower_bound, prev_lb)
        self.assertEquals(rxn.upper_bound, prev_ub)

    def test_setting_lower_bound_higher_than_higher_bound_sets_higher_bound_to_new_lower_bound(self):
        for reaction in self.model.reactions:
            self.assertTrue(reaction.lower_bound <= reaction.upper_bound)
            reaction.lower_bound = reaction.upper_bound + 100
            self.assertEqual(reaction.lower_bound, reaction.upper_bound)

    def test_setting_higher_bound_lower_than_lower_bound_sets_lower_bound_to_new_higher_bound(self):
        for reaction in self.model.reactions:
            self.assertTrue(reaction.lower_bound <= reaction.upper_bound)
            reaction.upper_bound = reaction.lower_bound - 100
            self.assertEqual(reaction.lower_bound, reaction.upper_bound)

    def test_add_metabolites_combine_true(self):
        test_metabolite = Metabolite('test')
        for reaction in self.model.reactions:
            reaction.add_metabolites({test_metabolite: -66}, combine=True)
            self.assertEqual(reaction.metabolites[test_metabolite], -66)
            self.assertIn(-66. * reaction.forward_variable, self.model.solver.constraints['test'].expression)
            self.assertIn(66. * reaction.reverse_variable, self.model.solver.constraints['test'].expression)
            already_included_metabolite = list(reaction.metabolites.keys())[0]
            previous_coefficient = reaction.get_coefficient(already_included_metabolite.id)
            reaction.add_metabolites({already_included_metabolite: 10}, combine=True)
            new_coefficient = previous_coefficient + 10
            self.assertEqual(reaction.metabolites[already_included_metabolite], new_coefficient)
            self.assertIn(new_coefficient * reaction.forward_variable,
                          self.model.solver.constraints[already_included_metabolite.id].expression)
            self.assertIn(-1 * new_coefficient * reaction.reverse_variable,
                          self.model.solver.constraints[already_included_metabolite.id].expression)

    @unittest.skipIf(TRAVIS, 'This test behaves non-deterministic on travis-ci')
    def test_add_metabolites_combine_false(self):
        test_metabolite = Metabolite('test')
        for reaction in self.model.reactions:
            reaction.add_metabolites({test_metabolite: -66}, combine=False)
            self.assertEqual(reaction.metabolites[test_metabolite], -66)
            self.assertIn(-66. * reaction.forward_variable, self.model.solver.constraints['test'].expression)
            self.assertIn(66. * reaction.reverse_variable, self.model.solver.constraints['test'].expression)
            already_included_metabolite = list(reaction.metabolites.keys())[0]
            reaction.add_metabolites({already_included_metabolite: 10}, combine=False)
            self.assertEqual(reaction.metabolites[already_included_metabolite], 10)
            self.assertIn(10 * reaction.forward_variable,
                          self.model.solver.constraints[already_included_metabolite.id].expression)
            self.assertIn(-10 * reaction.reverse_variable,
                          self.model.solver.constraints[already_included_metabolite.id].expression)

    @unittest.skip('Not implemented yet.')
    def test_change_id_is_reflected_in_solver(self):
        for i, reaction in enumerate(self.model.reactions):
            old_reaction_id = reaction.id
            self.assertTrue(self.model.solver.variables[old_reaction_id].name, old_reaction_id)
            self.assertIn(old_reaction_id, self.model.solver.variables)
            self.assertTrue(old_reaction_id in self.model.solver)
            new_reaction_id = reaction.id + '_' + str(i)
            reaction.id = new_reaction_id
            self.assertEqual(reaction.id, new_reaction_id)
            self.assertFalse(old_reaction_id in self.model.solver)
            self.assertTrue(new_reaction_id in self.model.solver)
            self.assertTrue(self.model.solver.variables[new_reaction_id].name, new_reaction_id)


class TestReactionGLPK(AbstractTestReaction, unittest.TestCase):
    def setUp(self):
        self.cobrapy_model = COBRAPYTESTMODEL.copy()
        self.model = TESTMODEL.copy()
        self.model.solver = 'glpk'


@unittest.skipIf(TRAVIS, 'CPLEX not available on Travis.')
class TestReactionCPLEX(AbstractTestReaction, unittest.TestCase):
    def setUp(self):
        self.cobrapy_model = COBRAPYTESTMODEL.copy()
        self.model = TESTMODEL.copy()
        self.model.solver = 'cplex'


class AbstractTestSolverBasedModel(object):
    def setUp(self):
        self.model = TESTMODEL.copy()
        self.model.solve()

    def test_model_from_other_cameo_model(self):
        model = Model(description=self.model)
        for reaction in model.reactions:
            self.assertEqual(reaction, self.model.reactions.get_by_id(reaction.id))

    # def test_reactions_and_variables_match(self):
    #     self.model.reversible_encoding = 'unsplit'
    #     reactions = self.model.reactions
    #     for reaction in reactions:
    #         self.assertIn(reaction.id, list(self.model.solver.variables.keys()))
    #         self.assertEqual(reaction.lower_bound, self.model.solver.variables[reaction.id].lb)
    #         self.assertEqual(reaction.upper_bound, self.model.solver.variables[reaction.id].ub)

    def test_add_reactions(self):
        r1 = Reaction('r1')
        r1.add_metabolites({Metabolite('A'): -1, Metabolite('B'): 1})
        r1.lower_bound, r1.upper_bound = -999999., 999999.
        r2 = Reaction('r2')
        r2.add_metabolites({Metabolite('A'): -1, Metabolite('C'): 1, Metabolite('D'): 1})
        r2.lower_bound, r2.upper_bound = 0., 999999.
        r2.objective_coefficient = 3.
        self.model.add_reactions([r1, r2])
        self.assertEqual(self.model.reactions[-2], r1)
        self.assertEqual(self.model.reactions[-1], r2)
        self.assertTrue(isinstance(self.model.reactions[-2].reverse_variable, self.model.solver.interface.Variable))
        self.assertEqual(self.model.objective.expression.coeff(
            self.model.reactions.Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2.forward_variable), 1.)
        self.assertEqual(self.model.objective.expression.coeff(
            self.model.reactions.Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2.reverse_variable), -1.)
        self.assertEqual(self.model.objective.expression.coeff(self.model.reactions.r2.forward_variable), 3.)
        self.assertEqual(self.model.objective.expression.coeff(self.model.reactions.r2.reverse_variable), -3.)

    def test_all_objects_point_to_all_other_correct_objects(self):
        model = load_model(os.path.join(TESTDIR, 'data/EcoliCore.xml'))
        for reaction in model.reactions:
            self.assertEqual(reaction.model, model)
            for gene in reaction.genes:
                self.assertEqual(gene, model.genes.get_by_id(gene.id))
                self.assertEqual(gene.model, model)
                for reaction2 in gene.reactions:
                    self.assertEqual(reaction2.model, model)
                    self.assertEqual(reaction2, model.reactions.get_by_id(reaction2.id))

            for metabolite in reaction.metabolites:
                self.assertEqual(metabolite.model, model)
                self.assertEqual(metabolite, model.metabolites.get_by_id(metabolite.id))
                for reaction2 in metabolite.reactions:
                    self.assertEqual(reaction2.model, model)
                    self.assertEqual(reaction2, model.reactions.get_by_id(reaction2.id))

    def test_all_objects_point_to_all_other_correct_objects_after_copy(self):
        model = load_model(os.path.join(TESTDIR, 'data/EcoliCore.xml'))
        model = model.copy()
        for reaction in model.reactions:
            self.assertEqual(reaction.model, model)
            for gene in reaction.genes:
                self.assertEqual(gene, model.genes.get_by_id(gene.id))
                self.assertEqual(gene.model, model)
                for reaction2 in gene.reactions:
                    self.assertEqual(reaction2.model, model)
                    self.assertEqual(reaction2, model.reactions.get_by_id(reaction2.id))

            for metabolite in reaction.metabolites:
                self.assertEqual(metabolite.model, model)
                self.assertEqual(metabolite, model.metabolites.get_by_id(metabolite.id))
                for reaction2 in metabolite.reactions:
                    self.assertEqual(reaction2.model, model)
                    self.assertEqual(reaction2, model.reactions.get_by_id(reaction2.id))

    def test_remove_reactions(self):
        reactions_to_remove = self.model.reactions[10:30]
        self.assertTrue(all([reaction.model is self.model for reaction in reactions_to_remove]))
        self.assertTrue(
            all([self.model.reactions.get_by_id(reaction.id) == reaction for reaction in reactions_to_remove]))

        self.model.remove_reactions(reactions_to_remove)
        self.assertTrue(all([reaction.model is None for reaction in reactions_to_remove]))
        for reaction in reactions_to_remove:
            self.assertNotIn(reaction.id, list(self.model.solver.variables.keys()))

        self.model.add_reactions(reactions_to_remove)
        for reaction in reactions_to_remove:
            self.assertIn(reaction, self.model.reactions)

    def test_add_demand(self):
        for metabolite in self.model.metabolites:
            demand_reaction = self.model.add_demand(metabolite, prefix="DemandReaction_")
            self.assertEqual(self.model.reactions.get_by_id(demand_reaction.id), demand_reaction)
            self.assertEqual(demand_reaction.reactants, [metabolite])
            self.assertTrue(
                self.model.solver.variables["DemandReaction_" + metabolite.id] in self.model.solver.constraints[
                    metabolite.id].expression)

    def test_objective(self):
        obj = self.model.objective
        self.assertEqual(
            obj.__str__(),
            'Maximize\n-1.0*Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2_reverse_9ebcd + 1.0*Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2')

    def test_change_objective(self):
        expression = 1.0 * self.model.solver.variables['ENO'] + 1.0 * self.model.solver.variables['PFK']
        self.model.objective = Objective(expression)
        self.assertEqual(str(self.model.objective.expression), str(expression))

    def test_set_reaction_objective(self):
        self.model.objective = self.model.reactions.ACALD
        self.assertEqual(str(self.model.objective.expression), str(
            1.0 * self.model.reactions.ACALD.variable - 1.0 * self.model.reactions.ACALD.reverse_variable))

    def test_set_reaction_objective_str(self):
        self.model.objective = self.model.reactions.ACALD.id
        self.assertEqual(str(self.model.objective.expression), str(
            1.0 * self.model.reactions.ACALD.variable - 1.0 * self.model.reactions.ACALD.reverse_variable))

    def test_invalid_objective_raises(self):
        self.assertRaises(Exception, setattr, self.model, 'objective', 'This is not a valid objective!')
        self.assertRaises(Exception, setattr, self.model, 'objective', 3.)

    def test_solver_change(self):
        solver_id = id(self.model.solver)
        problem_id = id(self.model.solver.problem)
        solution = self.model.solve().x_dict
        self.model.solver = 'glpk'
        self.assertNotEqual(id(self.model.solver), solver_id)
        self.assertNotEqual(id(self.model.solver.problem), problem_id)
        new_solution = self.model.solve()
        for key in list(solution.keys()):
            self.assertAlmostEqual(new_solution.x_dict[key], solution[key])

    def test_solver_change_with_optlang_interface(self):
        solver_id = id(self.model.solver)
        problem_id = id(self.model.solver.problem)
        solution = self.model.solve().x_dict
        self.model.solver = optlang.glpk_interface
        self.assertNotEqual(id(self.model.solver), solver_id)
        self.assertNotEqual(id(self.model.solver.problem), problem_id)
        new_solution = self.model.solve()
        for key in list(solution.keys()):
            self.assertAlmostEqual(new_solution.x_dict[key], solution[key])

    def test_invalid_solver_change_raises(self):
        self.assertRaises(ValueError, setattr, self.model, 'solver', [1, 2, 3])
        self.assertRaises(ValueError, setattr, self.model, 'solver', 'ThisIsDefinitelyNotAvalidSolver')
        self.assertRaises(ValueError, setattr, self.model, 'solver', os)

    @unittest.skipIf('cplex' not in solvers, "No cplex interface available")
    def test_change_solver_to_cplex_and_check_copy_works(self):
        # First, load model from scratch
        model = load_model(os.path.join(TESTDIR, 'data/EcoliCore.xml'), solver_interface='cplex')
        self.assertAlmostEqual(model.optimize().f, 0.8739215069684306)
        model_copy = model.copy()
        self.assertAlmostEqual(model_copy.optimize().f, 0.8739215069684306)
        # Second, change existing glpk based model to cplex
        self.model.solver = 'cplex'
        self.assertAlmostEqual(self.model.optimize().f, 0.8739215069684306)
        model_copy = copy.copy(self.model)
        self.assertAlmostEqual(model_copy.optimize().f, 0.8739215069684306)

    def test_copy_preserves_existing_solution(self):
        self.model.solve()  # TODO: not sure why the model has to be solved here because it is already in setUp
        model_cp = copy.copy(self.model)
        primals_original = [variable.primal for variable in self.model.solver.variables]
        primals_copy = [variable.primal for variable in model_cp.solver.variables]
        abs_diff = abs(numpy.array(primals_copy) - numpy.array(primals_original))
        self.assertFalse(any(abs_diff > 1e-6))

    def test_essential_genes(self):
        essential_genes = [g.id for g in self.model.essential_genes()]
        self.assertTrue(sorted(essential_genes) == sorted(ESSENTIAL_GENES))

    def test_essential_reactions(self):
        essential_reactions = [r.id for r in self.model.essential_reactions()]
        self.assertTrue(sorted(essential_reactions) == sorted(ESSENTIAL_REACTIONS))

    def test_effective_bounds(self):
        self.model.reactions.Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2.lower_bound = 0.873921
        for reaction in self.model.reactions:
            self.assertAlmostEqual(reaction.effective_lower_bound,
                                   REFERENCE_FVA_SOLUTION_ECOLI_CORE['lower_bound'][reaction.id], delta=0.000001)
            self.assertAlmostEqual(reaction.effective_upper_bound,
                                   REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound'][reaction.id], delta=0.000001)

    def test_add_demand_for_non_existing_metabolite(self):
        metabolite = Metabolite(id="a_metabolite")
        self.model.add_demand(metabolite)
        self.assertTrue(self.model.solver.variables["DM_" + metabolite.id]
                        in self.model.solver.constraints[metabolite.id].expression)

    def test_add_ratio_constraint(self):
        solution = self.model.solve()
        self.assertAlmostEqual(solution.f, 0.873921506968)
        self.assertNotEqual(2 * solution.x_dict['PGI'], solution.x_dict['G6PDH2r'])
        cp = self.model.copy()
        ratio_constr = cp.add_ratio_constraint(cp.reactions.PGI, cp.reactions.G6PDH2r, 0.5)
        self.assertEqual(ratio_constr.name, 'ratio_constraint_PGI_G6PDH2r')
        solution = cp.solve()
        self.assertAlmostEqual(solution.f, 0.870407873712)
        self.assertAlmostEqual(2 * solution.x_dict['PGI'], solution.x_dict['G6PDH2r'])
        cp = self.model.copy()

        ratio_constr = cp.add_ratio_constraint(cp.reactions.PGI, cp.reactions.G6PDH2r, 0.5)
        self.assertEqual(ratio_constr.name, 'ratio_constraint_PGI_G6PDH2r')
        solution = cp.solve()
        self.assertAlmostEqual(solution.f, 0.870407873712)
        self.assertAlmostEqual(2 * solution.x_dict['PGI'], solution.x_dict['G6PDH2r'])


class TestSolverBasedModelGLPK(AbstractTestSolverBasedModel, unittest.TestCase):
    def setUp(self):
        super(TestSolverBasedModelGLPK, self).setUp()
        self.model.solver = 'glpk'


@unittest.skipIf(TRAVIS, 'CPLEX not available on Travis.')
class TestSolverBasedModelCPLEX(AbstractTestSolverBasedModel, unittest.TestCase):
    def setUp(self):
        super(TestSolverBasedModelCPLEX, self).setUp()
        self.model.solver = 'glpk'


if __name__ == '__main__':
    import nose

    nose.runmodule()
