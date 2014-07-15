# Copyright (c) 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.
# See LICENSE for details.
import copy

import unittest

import os
from optlang import Objective
from cameo import load_model
from cameo.exceptions import UndefinedSolution
from cameo.solver_based_model import Reaction
from cobra.io import read_sbml_model


TESTDIR = os.path.dirname(__file__)
TESTMODEL = load_model(os.path.join(TESTDIR, 'data/EcoliCore.xml'))
ESSENTIAL_GENES = ['b2779', 'b1779', 'b0720', 'b0451', 'b2416', 'b2926', 'b1136', 'b2415']
ESSENTIAL_REACTIONS = ['GLNS', 'Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2', 'PIt2r', 'GAPD', 'ACONTb',
                       'EX_nh4_LPAREN_e_RPAREN_', 'ENO', 'EX_h_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_', 'ICDHyr',
                       'CS', 'NH4t', 'GLCpts', 'PGM', 'EX_pi_LPAREN_e_RPAREN_', 'PGK', 'RPI', 'ACONTa']


class CommonGround(unittest.TestCase):
    def setUp(self):
        self.model = TESTMODEL.copy()
        # self.model = TESTMODEL
        self.model.optimize()


class TestLazySolution(CommonGround):
    def test_self_invalidation(self):
        solution = self.model.optimize()
        self.assertAlmostEqual(solution.f, 0.873921506968431, delta=0.000001)
        self.model.optimize()
        self.assertRaises(UndefinedSolution, getattr, solution, 'f')


class TestReaction(unittest.TestCase):
    def setUp(self):
        self.cobrapy_model = read_sbml_model(os.path.join(TESTDIR, 'data/EcoliCore.xml'))

    def test_clone_cobrapy_reaction(self):
        for reaction in self.cobrapy_model.reactions:
            cloned_reaction = Reaction.clone(reaction)
            self.assertEqual(cloned_reaction.objective_coefficient, reaction.objective_coefficient)
            self.assertEqual(cloned_reaction.gene_reaction_rule, reaction.gene_reaction_rule)
            self.assertEqual(cloned_reaction.genes, reaction.genes)
            self.assertEqual(cloned_reaction.metabolites, reaction.metabolites)
            self.assertEqual(cloned_reaction.get_products(), reaction.get_products())
            self.assertEqual(cloned_reaction.get_reactants(), reaction.get_reactants())
            self.assertEqual(cloned_reaction.get_model(), None)
            self.assertEqual(cloned_reaction.variable, None)


class TestSolverBasedModel(CommonGround):
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
        self.assertEqual(
            obj.__str__(), 'Maximize\n1.0*Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2')

    def test_change_objective(self):
        self.model.objective = Objective(
            self.model.solver.variables['ENO'] + self.model.solver.variables['PFK'])
        self.assertEqual(self.model.objective.__str__(),
                         'Maximize\n1.0*ENO + 1.0*PFK')

    def test_change_solver_change(self):
        self.model.solver = 'glpk'

    def test_copy_preserves_existing_solution(self):
        model_cp = copy.copy(self.model)
        primals_original = [variable.primal for variable in self.model.solver.variables.values()]
        primals_copy = [variable.primal for variable in model_cp.solver.variables.values()]
        self.assertEqual(primals_copy, primals_original)

    def test_essential_genes(self):
        essential_genes = [g.id for g in self.model.essential_genes()]
        self.assertItemsEqual(essential_genes, ESSENTIAL_GENES)

    def test_essential_reactions(self):
        essential_reactions = [r.id for r in self.model.essential_reactions()]
        self.assertItemsEqual(essential_reactions, ESSENTIAL_REACTIONS)


if __name__ == '__main__':
    import nose
    nose.runmodule()