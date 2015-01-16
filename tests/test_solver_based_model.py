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
import copy
import unittest

from cobra import Metabolite
from optlang import Objective
from cameo import load_model, solvers
from cameo.exceptions import UndefinedSolution
from cameo.solver_based_model import Reaction
from cobra.io import read_sbml_model
import pandas

TESTDIR = os.path.dirname(__file__)
REFERENCE_FVA_SOLUTION_ECOLI_CORE = pandas.read_csv(os.path.join(TESTDIR, 'data/REFERENCE_flux_ranges_EcoliCore.csv'),
                                                    index_col=0)
TESTMODEL = load_model(os.path.join(TESTDIR, 'data/EcoliCore.xml'), sanitize=False)
COBRAPYTESTMODEL = read_sbml_model(os.path.join(TESTDIR, 'data/EcoliCore.xml'))
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
    def setUp(self):
        super(TestLazySolution, self).setUp()
        self.solution = self.model.optimize()

    def test_self_invalidation(self):
        solution = self.model.solve()
        self.assertAlmostEqual(solution.f, 0.873921506968431, delta=0.000001)
        self.model.optimize()
        self.assertRaises(UndefinedSolution, getattr, solution, 'f')

    def test_solution_contains_only_reaction_specific_values(self):
        reaction_IDs = set([reaction.id for reaction in self.model.reactions])
        for attr in ('x_dict', 'y_dict'):
            self.assertEqual(set(getattr(self.solution, attr).keys()).intersection(reaction_IDs), reaction_IDs)
            self.assertEqual(set(getattr(self.solution, attr).keys()).difference(reaction_IDs), set())

class TestReaction(unittest.TestCase):
    def setUp(self):
        self.cobrapy_model = COBRAPYTESTMODEL.copy()
        self.model = TESTMODEL.copy()

    def test_clone_cobrapy_reaction(self):
        for reaction in self.cobrapy_model.reactions:
            cloned_reaction = Reaction.clone(reaction)
            self.assertEqual(cloned_reaction.objective_coefficient, reaction.objective_coefficient)
            self.assertEqual(cloned_reaction.gene_reaction_rule, reaction.gene_reaction_rule)
            self.assertEqual(cloned_reaction.genes, reaction.genes)
            self.assertEqual(cloned_reaction.metabolites, reaction.metabolites)
            self.assertEqual(cloned_reaction.products, reaction.products)
            self.assertEqual(cloned_reaction.reactants, reaction.reactants)
            self.assertEqual(cloned_reaction.get_model(), None)
            self.assertEqual(cloned_reaction.variable, None)

    def test_knockout(self):
        for reaction in self.model.reactions:
            reaction.knock_out()
            self.assertEqual(reaction.lower_bound, 0)
            self.assertEqual(reaction.upper_bound, 0)
            self.assertEqual(self.model.solver.variables[reaction.id].lb, 0)
            self.assertEqual(self.model.solver.variables[reaction.id].ub, 0)

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
        self.assertEqual(acald_reaction.variable.lb, -1000099.0)
        self.assertEqual(acald_reaction.variable.ub, -1000099.0)
        self.assertEqual(acald_reaction.reverse_variable.lb, 0)
        self.assertEqual(acald_reaction.reverse_variable.ub, 0)
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
            print reaction.id
            self.assertTrue(reaction.lower_bound <= reaction.upper_bound)
            reaction.lower_bound = reaction.upper_bound + 100
            self.assertEqual(reaction.lower_bound, reaction.upper_bound)

    def test_setting_higher_bound_lower_than_lower_bound_sets_lower_bound_to_new_higher_bound(self):
        for reaction in self.model.reactions:
            self.assertTrue(reaction.lower_bound <= reaction.upper_bound)
            print reaction.id, reaction.lower_bound, reaction.upper_bound
            reaction.upper_bound = reaction.lower_bound - 100
            print reaction.id, reaction.lower_bound, reaction.upper_bound
            self.assertEqual(reaction.lower_bound, reaction.upper_bound)

    def test_add_metabolites(self):
        for reaction in self.model.reactions:
            reaction.add_metabolites({Metabolite('test'):-66})
            self.assertIn("66 test", str(reaction))
            self.assertIn(-66.*reaction.variable, self.model.solver.constraints['test'].expression)
            already_included_metabolite = reaction.metabolites.keys()[0]
            previous_coefficient = reaction.get_coefficient(already_included_metabolite.id)
            reaction.add_metabolites({already_included_metabolite: 10})
            new_coefficient = previous_coefficient + 10
            new_coefficient2 = new_coefficient
            if new_coefficient < 0:
                new_coefficient *= -1
            if new_coefficient % 1 == 0:
                new_coefficient = new_coefficient
            self.assertIn(str(new_coefficient)+" "+already_included_metabolite.id, str(reaction))
            self.assertIn(new_coefficient2*reaction.variable, self.model.solver.constraints[already_included_metabolite.id].expression)

class TestSolverBasedModel(CommonGround):
    def test_reactions_and_variables_match(self):
        self.model.reversible_encoding = 'unsplit'
        reactions = self.model.reactions
        for reaction in reactions:
            self.assertIn(reaction.id, self.model.solver.variables.keys())
            self.assertEqual(reaction.lower_bound, self.model.solver.variables[reaction.id].lb)
            self.assertEqual(reaction.upper_bound, self.model.solver.variables[reaction.id].ub)

    def test_all_objects_point_to_all_other_correct_objects(self):
        model = load_model(os.path.join(TESTDIR, 'data/EcoliCore.xml'))
        for reaction in model.reactions:
            self.assertEqual(reaction.get_model(), model)
            for gene in reaction.genes:
                self.assertEqual(gene, model.genes.get_by_id(gene.id))
                self.assertEqual(gene.model, model)
                for reaction2 in gene.reactions:
                    self.assertEqual(reaction2.get_model(), model)
                    self.assertEqual(reaction2, model.reactions.get_by_id(reaction2.id))

            for metabolite in reaction.metabolites:
                self.assertEqual(metabolite.model, model)
                self.assertEqual(metabolite, model.metabolites.get_by_id(metabolite.id))
                for reaction2 in metabolite.reactions:
                    self.assertEqual(reaction2.get_model(), model)
                    self.assertEqual(reaction2, model.reactions.get_by_id(reaction2.id))

    def test_all_objects_point_to_all_other_correct_objects_after_copy(self):
        model = load_model(os.path.join(TESTDIR, 'data/EcoliCore.xml'))
        model = model.copy()
        for reaction in model.reactions:
            self.assertEqual(reaction.get_model(), model)
            for gene in reaction.genes:
                self.assertEqual(gene, model.genes.get_by_id(gene.id))
                self.assertEqual(gene.model, model)
                for reaction2 in gene.reactions:
                    self.assertEqual(reaction2.get_model(), model)
                    self.assertEqual(reaction2, model.reactions.get_by_id(reaction2.id))

            for metabolite in reaction.metabolites:
                self.assertEqual(metabolite.model, model)
                self.assertEqual(metabolite, model.metabolites.get_by_id(metabolite.id))
                for reaction2 in metabolite.reactions:
                    self.assertEqual(reaction2.get_model(), model)
                    self.assertEqual(reaction2, model.reactions.get_by_id(reaction2.id))

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
            self.assertTrue(
                self.model.solver.variables["DemandReaction_" + metabolite.id] in self.model.solver.constraints[
                    metabolite.id].expression)

    def test_add_demand_for_non_existing_metabolite(self):
        metabolite = Metabolite(id="a_metabolite")
        self.model.add_demand(metabolite)
        self.assertTrue(self.model.solver.variables["DM_" + metabolite.id]
                        in self.model.solver.constraints[metabolite.id].expression)

    def test_objective(self):
        obj = self.model.objective
        self.assertEqual(
            obj.__str__(), 'Maximize\n1.0*Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2')

    def test_change_objective(self):
        self.model.objective = Objective(
            self.model.solver.variables['ENO'] + self.model.solver.variables['PFK'])

    def test_solver_change(self):
        solver_id = id(self.model.solver)
        problem_id = id(self.model.solver.problem)
        solution = self.model.solve().x_dict
        self.model.solver = 'glpk'
        self.assertNotEqual(id(self.model.solver), solver_id)
        self.assertNotEqual(id(self.model.solver.problem), problem_id)
        new_solution = self.model.solve()
        for key in solution.keys():
            self.assertAlmostEqual(new_solution.x_dict[key], solution[key])

    @unittest.skipIf(not solvers.has_key('cplex'), "No cplex interface available")
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
        model_cp = copy.copy(self.model)
        primals_original = [variable.primal for variable in self.model.solver.variables]
        primals_copy = [variable.primal for variable in model_cp.solver.variables]
        self.assertEqual(primals_copy, primals_original)

    def test_essential_genes(self):
        self.model.reversible_encoding = 'split'
        essential_genes = [g.id for g in self.model.essential_genes()]
        self.assertItemsEqual(essential_genes, ESSENTIAL_GENES)

    def test_essential_reactions(self):
        essential_reactions = [r.id for r in self.model.essential_reactions()]
        self.assertItemsEqual(essential_reactions, ESSENTIAL_REACTIONS)

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
        self.assertNotEqual(2*solution.x_dict['PGI'], solution.x_dict['G6PDH2r'])
        cp = self.model.copy()
        ratio_constr = cp.add_ratio_constraint(cp.reactions.PGI, cp.reactions.G6PDH2r, 0.5)
        self.assertEqual(ratio_constr.name, 'ratio_constraint_PGI_G6PDH2r')
        solution = cp.solve()
        self.assertAlmostEqual(solution.f, 0.870407873712)
        self.assertAlmostEqual(2*solution.x_dict['PGI'], solution.x_dict['G6PDH2r'])
        cp = self.model.copy()
        cp.reversible_encoding = 'unsplit'
        ratio_constr = cp.add_ratio_constraint(cp.reactions.PGI, cp.reactions.G6PDH2r, 0.5)
        self.assertEqual(ratio_constr.name, 'ratio_constraint_PGI_G6PDH2r')
        solution = cp.solve()
        self.assertAlmostEqual(solution.f, 0.870407873712)
        self.assertAlmostEqual(2*solution.x_dict['PGI'], solution.x_dict['G6PDH2r'])


if __name__ == '__main__':
    import nose

    nose.runmodule()