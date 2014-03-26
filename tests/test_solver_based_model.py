# Copyright (c) 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.
# See LICENSE for details.

import os
import unittest
from cameo.io import load_model
from optlang import Objective

TESTDIR = os.path.dirname(__file__)

class TestOptlangBasedModel(unittest.TestCase):

    def setUp(self):
        self.model = load_model(os.path.join(TESTDIR, 'data/EcoliCore.xml'))

    def test_reactions_and_variables_match(self):
        reactions = self.model.reactions
        for reaction in reactions:
            self.assertIn(reaction.id, self.model.solver.variables.keys())
            self.assertEqual(reaction.lower_bound, self.model.solver.variables[reaction.id].lb)
            self.assertEqual(reaction.upper_bound, self.model.solver.variables[reaction.id].ub)

    def test_remove_reactions(self):
        reactions_to_remove = self.model.reactions[10:30]
        self.assertTrue(all([reaction.get_model() is self.model for reaction in reactions_to_remove]))
        self.assertTrue(all([self.model.reactions.get_by_id(reaction.id) == reaction for reaction in reactions_to_remove]))
        
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
        self.assertEqual(
            obj.__str__(), 'Maximize\n1.0*Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2')

    def test_change_objective(self):
        self.model.objective = Objective(
            self.model.solver.variables['ENO'] + self.model.solver.variables['PFK'])
        self.assertEqual(self.model.objective.__str__(),
                         'Maximize\n1.0*ENO + 1.0*PFK')

if __name__ == '__main__':
    import nose
    nose.runmodule()
    
