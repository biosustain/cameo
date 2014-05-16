# Copyright (c) 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.
# See LICENSE for details.
import copy

import unittest

import os
from optlang import Objective

from cameo import load_model


TESTDIR = os.path.dirname(__file__)
TESTMODEL = load_model(os.path.join(TESTDIR, 'data/EcoliCore.xml'))

class CommonGround(unittest.TestCase):

    def setUp(self):
        self.model = TESTMODEL.copy()
        self.model.optimize()

class TestLazySolution(CommonGround):

    def test_self_invalidation(self):
        solution = self.model.optimize()
        self.assertAlmostEqual(solution.f, 0.873921506968431)
        self.model.optimize()
        self.assertRaises(Exception, solution, f)


class TestOptlangBasedModel(unittest.TestCase):
    def setUp(self):
        self.model = TESTMODEL.copy()
        self.model.optimize()

    def test_reactions_and_variables_match(self):
        reactions = self.model.reactions
        for reaction in reactions:
            self.assertIn(reaction.id, self.model.solver.variables.keys())
            self.assertEqual(reaction.lower_bound, self.model.solver.variables[reaction.id].lb)
            self.assertEqual(reaction.upper_bound, self.model.solver.variables[reaction.id].ub)

    def test_remove_reactions(self):
        reactions_to_remove = self.model.reactions[10:30]
        self.assertTrue(all([reaction.get_model() is self.model for reaction in reactions_to_remove]))
        self.assertTrue(
            all([self.model.reactions.get_by_id(reaction.id) == reaction for reaction in reactions_to_remove]))

        self.model.remove_reactions(reactions_to_remove)
        self.assertTrue(all([reaction.get_model() is None for reaction in reactions_to_remove]))
        for reaction in reactions_to_remove:
            self.assertNotIn(reaction.id, self.model.solver.variables.keys())

    def test_add_demand(self):
        for metabolite in self.model.metabolites:
            demand_reaction = self.model.add_demand(metabolite, prefix="DemandReaction_")
            self.assertEqual(self.model.reactions.get_by_id(demand_reaction.id), demand_reaction)
            self.assertEqual(demand_reaction.reactants, [metabolite])

    def test_objective(self):
        obj = self.model.objective
        print obj
        self.assertEqual(
            obj.__str__(), 'Maximize\n1.0*Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2')

    def test_change_objective(self):
        self.model.objective = Objective(
            self.model.solver.variables['ENO'] + self.model.solver.variables['PFK'])
        self.assertEqual(self.model.objective.__str__(),
                         'Maximize\n1.0*ENO + 1.0*PFK')

    def test_change_solver(self):
        pass

    def test_copy_preserves_existing_solution(self):
        model_cp = copy.copy(self.model)
        primals_original = [variable.primal for variable in self.model.solver.variables.values()]
        primals_copy = [variable.primal for variable in model_cp.solver.variables.values()]
        self.assertEqual(primals_copy, primals_original)



if __name__ == '__main__':
    import nose

    nose.runmodule()
    
