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

import copy
import os
import pickle

import cobra.test
from cobra.util import SolverNotFound

import numpy
import pandas
import pytest

from cobra.util import fix_objective_as_constraint
from cobra import Model, Reaction, Metabolite
from cobra.exceptions import OptimizationError

from cameo import load_model
from cameo.config import solvers
from cameo.core.utils import get_reaction_for, load_medium, medium
from cameo.flux_analysis.structural import create_stoichiometric_array
from cameo.flux_analysis.analysis import find_essential_metabolites
from cobra.flux_analysis import find_essential_genes, find_essential_reactions


TRAVIS = bool(os.getenv('TRAVIS', False))
TESTDIR = os.path.dirname(__file__)
REFERENCE_FVA_SOLUTION_ECOLI_CORE = pandas.read_csv(os.path.join(TESTDIR, 'data/REFERENCE_flux_ranges_EcoliCore.csv'),
                                                    index_col=0)
ESSENTIAL_GENES = ['b2779', 'b1779', 'b0720', 'b0451', 'b2416', 'b2926', 'b1136', 'b2415']
ESSENTIAL_METABOLITES = ['13dpg_c', '2pg_c', '3pg_c', 'accoa_c', 'acon_DASH_C_c', 'adp_c', 'akg_c', 'atp_c', 'cit_c',
                         'coa_c', 'e4p_c', 'f6p_c', 'g3p_c', 'g6p_c', 'glc_DASH_D_e', 'gln_DASH_L_c', 'glu_DASH_L_c',
                         'h2o_c', 'h_c', 'h_e', 'icit_c', 'nad_c', 'nadh_c', 'nadp_c', 'nadph_c', 'nh4_c', 'nh4_e',
                         'oaa_c', 'pep_c', 'pi_c', 'pi_e', 'pyr_c', 'r5p_c']
ESSENTIAL_REACTIONS = ['GLNS', 'Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2', 'PIt2r', 'GAPD', 'ACONTb',
                       'EX_nh4_LPAREN_e_RPAREN_', 'ENO', 'EX_h_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_', 'ICDHyr',
                       'CS', 'NH4t', 'GLCpts', 'PGM', 'EX_pi_LPAREN_e_RPAREN_', 'PGK', 'RPI', 'ACONTa']


@pytest.fixture(scope="function", params=list(solvers))
def solved_model(request, data_directory):
    core_model = load_model(os.path.join(data_directory, 'EcoliCore.xml'), sanitize=False)
    core_model.solver = request.param
    solution = core_model.optimize()
    return solution, core_model


@pytest.fixture(scope="module", params=list(solvers))
def tiny_toy_model(request):
    tiny = Model("Toy Model")
    m1 = Metabolite("M1")
    d1 = Reaction("ex1")
    d1.add_metabolites({m1: -1})
    d1.bounds = -1000, 0
    tiny.add_reactions([d1])
    tiny.solver = request.param
    return tiny


# class TestLazySolution:
#     def test_self_invalidation(self, solved_model):
#         solution, model = solved_model
#         assert abs(solution.objective_value - 0.873921506968431) < 0.000001
#         model.optimize()
#         with pytest.raises(UndefinedSolution):
#             getattr(solution, 'f')
#
#     def test_solution_contains_only_reaction_specific_values(self, solved_model):
#         solution, model = solved_model
#         reaction_ids = set([reaction.id for reaction in model.reactions])
#         assert set(solution.fluxes.keys()).difference(reaction_ids) == set()
#         assert set(solution.reduced_costs.keys()).difference(reaction_ids) == set()
#         assert set(solution.reduced_costs.keys()).difference(reaction_ids) == set()
#         metabolite_ids = set([metabolite.id for metabolite in model.metabolites])
#         assert set(solution.shadow_prices.keys()).difference(metabolite_ids) == set()


class TestReaction:
    # def test_clone_cobrapy_reaction(self):
    #     model = cobra.test.create_test_model('textbook')
    #     for reaction in model.reactions:
    #         cloned_reaction = Reaction.clone(reaction)
    #         assert cloned_reaction.gene_reaction_rule == reaction.gene_reaction_rule
    #         assert set([gene.id for gene in cloned_reaction.genes]) == set([gene.id for gene in reaction.genes])
    #         assert all(isinstance(gene, cobra.Gene) for gene in list(cloned_reaction.genes))
    #         assert {metabolite.id for metabolite in cloned_reaction.metabolites} == {metabolite.id for metabolite in
    #                                                                                  reaction.metabolites}
    #         assert all(isinstance(metabolite, cameo.core.Metabolite) for metabolite in cloned_reaction.metabolites)
    #         assert {metabolite.id for metabolite in cloned_reaction.products} == {metabolite.id for metabolite in
    #                                                                               reaction.products}
    #         assert {metabolite.id for metabolite in cloned_reaction.reactants} == {metabolite.id for metabolite in
    #                                                                                reaction.reactants}
    #         assert reaction.id == cloned_reaction.id
    #         assert reaction.name == cloned_reaction.name
    #         assert reaction.upper_bound == cloned_reaction.upper_bound
    #         assert reaction.lower_bound == cloned_reaction.lower_bound

    # test moved to cobra
    # def test_gene_reaction_rule_setter(self, core_model):
    #     rxn = Reaction('rxn')
    #     rxn.add_metabolites({Metabolite('A'): -1, Metabolite('B'): 1})
    #     rxn.gene_reaction_rule = 'A2B1 or A2B2 and A2B3'
    #     assert hasattr(list(rxn.genes)[0], 'knock_out')
    #     core_model.add_reaction(rxn)
    #     with cameo.util.TimeMachine() as tm:
    #         core_model.genes.A2B1.knock_out(time_machine=tm)
    #         assert not core_model.genes.A2B1.functional
    #         core_model.genes.A2B3.knock_out(time_machine=tm)
    #         assert not rxn.functional
    #     assert core_model.genes.A2B3.functional
    #     assert rxn.functional
    #     core_model.genes.A2B1.knock_out()
    #     assert not core_model.genes.A2B1.functional
    #     assert core_model.reactions.rxn.functional
    #     core_model.genes.A2B3.knock_out()
    #     assert not core_model.reactions.rxn.functional
    #     non_functional = [gene.id for gene in core_model.non_functional_genes]
    #     assert all(gene in non_functional for gene in ['A2B3', 'A2B1'])

    def test_gene_reaction_rule_setter_reaction_already_added_to_model(self, core_model):
        rxn = Reaction('rxn')
        rxn.add_metabolites({Metabolite('A'): -1, Metabolite('B'): 1})
        core_model.add_reaction(rxn)
        rxn.gene_reaction_rule = 'A2B'
        assert hasattr(list(rxn.genes)[0], 'knock_out')

    def test_str(self, core_model):
        assert core_model.reactions[0].__str__().startswith('ACALD')

    def test_add_metabolite(self, solved_model):
        solution, model = solved_model
        pgi_reaction = model.reactions.PGI
        test_met = model.metabolites[0]
        pgi_reaction.add_metabolites({test_met: 42}, combine=False)
        constraints = model.solver.constraints
        assert pgi_reaction.metabolites[test_met] == 42
        assert constraints[test_met.id].expression.as_coefficients_dict()[pgi_reaction.forward_variable] == 42
        assert constraints[test_met.id].expression.as_coefficients_dict()[pgi_reaction.reverse_variable] == -42

        pgi_reaction.add_metabolites({test_met: -10}, combine=True)
        assert pgi_reaction.metabolites[test_met] == 32
        assert constraints[test_met.id].expression.as_coefficients_dict()[pgi_reaction.forward_variable] == 32
        assert constraints[test_met.id].expression.as_coefficients_dict()[pgi_reaction.reverse_variable] == -32

        pgi_reaction.add_metabolites({test_met: 0}, combine=False)
        with pytest.raises(KeyError):
            assert pgi_reaction.metabolites[test_met]
        assert constraints[test_met.id].expression.as_coefficients_dict()[pgi_reaction.forward_variable] == 0
        assert constraints[test_met.id].expression.as_coefficients_dict()[pgi_reaction.reverse_variable] == 0

    def test_removal_from_model_retains_bounds(self, core_model):
        core_model_cp = core_model.copy()
        reaction = core_model_cp.reactions.ACALD
        assert reaction.model == core_model_cp
        assert reaction.lower_bound == -1000.0
        assert reaction.upper_bound == 1000.0
        assert reaction._lower_bound == -1000.0
        assert reaction._upper_bound == 1000.0
        core_model_cp.remove_reactions([reaction])
        assert reaction.model is None
        assert reaction.lower_bound == -1000.0
        assert reaction.upper_bound == 1000.0
        assert reaction._lower_bound == -1000.0
        assert reaction._upper_bound == 1000.0

    def test_set_bounds_scenario_1(self, core_model):
        acald_reaction = core_model.reactions.ACALD
        assert acald_reaction.lower_bound == -1000.
        assert acald_reaction.upper_bound == 1000.
        assert acald_reaction.forward_variable.lb == 0.
        assert acald_reaction.forward_variable.ub == 1000.
        assert acald_reaction.reverse_variable.lb == 0
        assert acald_reaction.reverse_variable.ub == 1000.
        acald_reaction.upper_bound = acald_reaction.lower_bound - 100
        assert acald_reaction.lower_bound == -1100.0
        assert acald_reaction.upper_bound == -1100.0
        assert acald_reaction.forward_variable.lb == 0
        assert acald_reaction.forward_variable.ub == 0
        assert acald_reaction.reverse_variable.lb == 1100.
        assert acald_reaction.reverse_variable.ub == 1100.
        acald_reaction.upper_bound = 100
        assert acald_reaction.lower_bound == -1100.0
        assert acald_reaction.upper_bound == 100
        assert acald_reaction.forward_variable.lb == 0
        assert acald_reaction.forward_variable.ub == 100
        assert acald_reaction.reverse_variable.lb == 0
        assert acald_reaction.reverse_variable.ub == 1100.0

    def test_set_bounds_scenario_3(self, core_model):
        reac = core_model.reactions.ACALD
        reac.bounds = -10, -10
        assert reac.lower_bound == -10
        assert reac.upper_bound == -10
        reac.lower_bound = -9
        assert reac.lower_bound == -9
        assert reac.upper_bound == -9
        reac.lower_bound = 2
        assert reac.lower_bound == 2
        assert reac.upper_bound == 2
        reac.upper_bound = -10
        assert reac.lower_bound == -10
        assert reac.upper_bound == -10
        reac.upper_bound = -11
        assert reac.lower_bound == -11
        assert reac.upper_bound == -11
        reac.upper_bound = 2
        assert reac.lower_bound == -11
        assert reac.upper_bound == 2

    def test_set_bounds_scenario_4(self, core_model):
        reac = core_model.reactions.ACALD
        reac.lower_bound = reac.upper_bound = 0
        reac.lower_bound = 2
        assert reac.lower_bound == 2
        assert reac.upper_bound == 2
        assert reac.forward_variable.lb == 2
        assert reac.forward_variable.ub == 2
        reac.knock_out()
        reac.upper_bound = -2
        assert reac.lower_bound == -2
        assert reac.upper_bound == -2
        assert reac.reverse_variable.lb == 2
        assert reac.reverse_variable.ub == 2

    def test_set_upper_before_lower_bound_to_0(self, core_model):
        core_model.reactions.GAPD.bounds = 0, 0
        core_model.reactions.GAPD.lower_bound = 0
        assert core_model.reactions.GAPD.lower_bound == 0
        assert core_model.reactions.GAPD.upper_bound == 0
        assert core_model.reactions.GAPD.forward_variable.lb == 0
        assert core_model.reactions.GAPD.forward_variable.ub == 0
        assert core_model.reactions.GAPD.reverse_variable.lb == 0
        assert core_model.reactions.GAPD.reverse_variable.ub == 0

    def test_set_bounds_scenario_2(self, core_model):
        acald_reaction = core_model.reactions.ACALD
        assert acald_reaction.lower_bound == -1000.
        assert acald_reaction.upper_bound == 1000.
        assert acald_reaction.forward_variable.lb == 0.
        assert acald_reaction.forward_variable.ub == 1000.
        assert acald_reaction.reverse_variable.lb == 0
        assert acald_reaction.reverse_variable.ub == 1000.
        acald_reaction.lower_bound = acald_reaction.upper_bound + 100
        assert acald_reaction.lower_bound == 1100.0
        assert acald_reaction.upper_bound == 1100.0
        assert acald_reaction.forward_variable.lb == 1100.0
        assert acald_reaction.forward_variable.ub == 1100.0
        assert acald_reaction.reverse_variable.lb == 0
        assert acald_reaction.reverse_variable.ub == 0
        acald_reaction.lower_bound = -100
        assert acald_reaction.lower_bound == -100.
        assert acald_reaction.upper_bound == 1100.
        assert acald_reaction.forward_variable.lb == 0
        assert acald_reaction.forward_variable.ub == 1100.
        assert acald_reaction.reverse_variable.lb == 0
        assert acald_reaction.reverse_variable.ub == 100

    def test_change_bounds(self, core_model):
        reac = core_model.reactions.ACALD
        reac.bounds = (2, 2)
        assert reac.lower_bound == 2
        assert reac.upper_bound == 2
        with core_model:
            reac.lower_bound = 5
            assert reac.lower_bound == 5
            assert reac.upper_bound == 5
        assert reac.lower_bound == 2
        assert reac.upper_bound == 2

    def test_make_irreversible(self, core_model):
        acald_reaction = core_model.reactions.ACALD
        assert acald_reaction.lower_bound == -1000.
        assert acald_reaction.upper_bound == 1000.
        assert acald_reaction.forward_variable.lb == 0.
        assert acald_reaction.forward_variable.ub == 1000.
        assert acald_reaction.reverse_variable.lb == 0
        assert acald_reaction.reverse_variable.ub == 1000.
        acald_reaction.lower_bound = 0
        assert acald_reaction.lower_bound == 0
        assert acald_reaction.upper_bound == 1000.
        assert acald_reaction.forward_variable.lb == 0
        assert acald_reaction.forward_variable.ub == 1000.0
        assert acald_reaction.reverse_variable.lb == 0
        assert acald_reaction.reverse_variable.ub == 0
        acald_reaction.lower_bound = -100
        assert acald_reaction.lower_bound == -100.
        assert acald_reaction.upper_bound == 1000.
        assert acald_reaction.forward_variable.lb == 0
        assert acald_reaction.forward_variable.ub == 1000.
        assert acald_reaction.reverse_variable.lb == 0
        assert acald_reaction.reverse_variable.ub == 100

    def test_make_reversible(self, core_model):
        pfk_reaction = core_model.reactions.PFK
        assert pfk_reaction.lower_bound == 0.
        assert pfk_reaction.upper_bound == 1000.
        assert pfk_reaction.forward_variable.lb == 0.
        assert pfk_reaction.forward_variable.ub == 1000.
        assert pfk_reaction.reverse_variable.lb == 0
        assert pfk_reaction.reverse_variable.ub == 0
        pfk_reaction.lower_bound = -100.
        assert pfk_reaction.lower_bound == -100.
        assert pfk_reaction.upper_bound == 1000.
        assert pfk_reaction.forward_variable.lb == 0
        assert pfk_reaction.forward_variable.ub == 1000.0
        assert pfk_reaction.reverse_variable.lb == 0
        assert pfk_reaction.reverse_variable.ub == 100.
        pfk_reaction.lower_bound = 0
        assert pfk_reaction.lower_bound == 0
        assert pfk_reaction.upper_bound == 1000.
        assert pfk_reaction.forward_variable.lb == 0
        assert pfk_reaction.forward_variable.ub == 1000.
        assert pfk_reaction.reverse_variable.lb == 0
        assert pfk_reaction.reverse_variable.ub == 0

    def test_make_irreversible_irreversible_to_the_other_side(self, core_model):
        pfk_reaction = core_model.reactions.PFK
        assert pfk_reaction.lower_bound == 0.
        assert pfk_reaction.upper_bound == 1000.
        assert pfk_reaction.forward_variable.lb == 0.
        assert pfk_reaction.forward_variable.ub == 1000.
        assert pfk_reaction.reverse_variable.lb == 0
        assert pfk_reaction.reverse_variable.ub == 0
        pfk_reaction.upper_bound = -100.
        assert pfk_reaction.forward_variable.lb == 0
        assert pfk_reaction.forward_variable.ub == 0
        assert pfk_reaction.reverse_variable.lb == 100
        assert pfk_reaction.reverse_variable.ub == 100
        pfk_reaction.lower_bound = -1000.
        assert pfk_reaction.lower_bound == -1000.
        assert pfk_reaction.upper_bound == -100.
        assert pfk_reaction.forward_variable.lb == 0
        assert pfk_reaction.forward_variable.ub == 0
        assert pfk_reaction.reverse_variable.lb == 100
        assert pfk_reaction.reverse_variable.ub == 1000.

    def test_make_lhs_irreversible_reversible(self, core_model):
        rxn = Reaction('test')
        rxn.add_metabolites(
            {core_model.metabolites[0]: -1., core_model.metabolites[1]: 1.})
        rxn.lower_bound = -1000.
        rxn.upper_bound = -100
        core_model.add_reaction(rxn)
        assert rxn.lower_bound == -1000.
        assert rxn.upper_bound == -100.
        assert rxn.forward_variable.lb == 0.
        assert rxn.forward_variable.ub == 0.
        assert rxn.reverse_variable.lb == 100.
        assert rxn.reverse_variable.ub == 1000.
        rxn.upper_bound = 666.
        assert rxn.lower_bound == -1000.
        assert rxn.upper_bound == 666.
        assert rxn.forward_variable.lb == 0.
        assert rxn.forward_variable.ub == 666
        assert rxn.reverse_variable.lb == 0.
        assert rxn.reverse_variable.ub == 1000.

    def test_model_less_reaction(self, solved_model):
        solution, model = solved_model
        for reaction in model.reactions:
            assert isinstance(reaction.flux, float)
            assert isinstance(reaction.reduced_cost, float)
        for reaction in model.reactions:
            model.remove_reactions([reaction])
            with pytest.raises(RuntimeError):
                assert reaction.flux
            with pytest.raises(RuntimeError):
                assert reaction.reduced_cost

    def test_knockout(self, core_model):
        original_bounds = dict()
        for reaction in core_model.reactions:
            original_bounds[reaction.id] = (
                reaction.lower_bound, reaction.upper_bound)
            reaction.knock_out()
            assert reaction.lower_bound == 0
            assert reaction.upper_bound == 0
        for k, (lb, ub) in original_bounds.items():
            core_model.reactions.get_by_id(k).lower_bound = lb
            core_model.reactions.get_by_id(k).upper_bound = ub
        for reaction in core_model.reactions:
            assert reaction.lower_bound == original_bounds[reaction.id][0]
            assert reaction.upper_bound == original_bounds[reaction.id][1]
        with core_model:
            for reaction in core_model.reactions:
                original_bounds[reaction.id] = (
                    reaction.lower_bound, reaction.upper_bound)
                reaction.knock_out()
                assert reaction.lower_bound == 0
                assert reaction.upper_bound == 0
        for reaction in core_model.reactions:
            assert reaction.lower_bound == original_bounds[reaction.id][0]
            assert reaction.upper_bound == original_bounds[reaction.id][1]

    @pytest.mark.xfail(reason="to be implemented in cobra")
    def test_repr_html_(self, core_model):
        assert '<table>' in core_model.reactions[0]._repr_html_()

    def test_reaction_without_model(self):
        r = Reaction('blub')
        assert r.flux_expression is None
        assert r.forward_variable is None
        assert r.reverse_variable is None

    def test_weird_left_to_right_reaction_issue(self, tiny_toy_model):
        d1 = tiny_toy_model.reactions.get_by_id('ex1')
        assert not d1.reversibility
        assert d1.lower_bound == -1000
        assert d1._lower_bound == -1000
        assert d1.upper_bound == 0
        assert d1._upper_bound == 0
        with tiny_toy_model:
            d1.knock_out()
            assert d1.lower_bound == 0
            assert d1._lower_bound == 0
            assert d1.upper_bound == 0
            assert d1._upper_bound == 0
        assert d1.lower_bound == -1000
        assert d1._lower_bound == -1000
        assert d1.upper_bound == 0
        assert d1._upper_bound == 0

    def test_one_left_to_right_reaction_set_positive_ub(self, tiny_toy_model):
        d1 = tiny_toy_model.reactions.get_by_id('ex1')
        assert d1.reverse_variable.lb == 0
        assert d1.reverse_variable.ub == 1000
        assert d1._lower_bound == -1000
        assert d1.lower_bound == -1000
        assert d1._upper_bound == 0
        assert d1.upper_bound == 0
        assert d1.forward_variable.lb == 0
        assert d1.forward_variable.ub == 0
        d1.upper_bound = .1
        assert d1.forward_variable.lb == 0
        assert d1.forward_variable.ub == .1
        assert d1.reverse_variable.lb == 0
        assert d1.reverse_variable.ub == 1000
        assert d1._lower_bound == -1000
        assert d1.upper_bound == .1
        assert d1._lower_bound == -1000
        assert d1.upper_bound == .1

    def test_irrev_reaction_set_negative_lb(self, core_model):
        assert not core_model.reactions.PFK.reversibility
        assert core_model.reactions.PFK.lower_bound == 0
        assert core_model.reactions.PFK.upper_bound == 1000.0
        assert core_model.reactions.PFK.forward_variable.lb == 0
        assert core_model.reactions.PFK.forward_variable.ub == 1000.0
        assert core_model.reactions.PFK.reverse_variable.lb == 0
        assert core_model.reactions.PFK.reverse_variable.ub == 0
        core_model.reactions.PFK.lower_bound = -1000
        assert core_model.reactions.PFK.lower_bound == -1000
        assert core_model.reactions.PFK.upper_bound == 1000.0
        assert core_model.reactions.PFK.forward_variable.lb == 0
        assert core_model.reactions.PFK.forward_variable.ub == 1000.0
        assert core_model.reactions.PFK.reverse_variable.lb == 0
        assert core_model.reactions.PFK.reverse_variable.ub == 1000

    def test_twist_irrev_right_to_left_reaction_to_left_to_right(self, core_model):
        assert not core_model.reactions.PFK.reversibility
        assert core_model.reactions.PFK.lower_bound == 0
        assert core_model.reactions.PFK.upper_bound == 1000.0
        assert core_model.reactions.PFK.forward_variable.lb == 0
        assert core_model.reactions.PFK.forward_variable.ub == 1000.0
        assert core_model.reactions.PFK.reverse_variable.lb == 0
        assert core_model.reactions.PFK.reverse_variable.ub == 0
        core_model.reactions.PFK.lower_bound = -1000
        core_model.reactions.PFK.upper_bound = 0
        assert core_model.reactions.PFK.lower_bound == -1000
        assert core_model.reactions.PFK.upper_bound == 0
        assert core_model.reactions.PFK.forward_variable.lb == 0
        assert core_model.reactions.PFK.forward_variable.ub == 0
        assert core_model.reactions.PFK.reverse_variable.lb == 0
        assert core_model.reactions.PFK.reverse_variable.ub == 1000

    @pytest.mark.skipif(TRAVIS, reason='too slow for ci')
    def test_imm904_4hglsdm_problem(self, imm904):
        # set upper bound before lower bound after knockout
        cp = imm904.copy()
        rxn = cp.reactions.get_by_id('4HGLSDm')
        prev_lb, prev_ub = rxn.lower_bound, rxn.upper_bound
        rxn.lower_bound = 0
        rxn.upper_bound = 0
        rxn.upper_bound = prev_ub
        rxn.lower_bound = prev_lb
        assert rxn.lower_bound == prev_lb
        assert rxn.upper_bound == prev_ub
        # set lower bound before upper bound after knockout
        cp = imm904.copy()
        rxn = cp.reactions.get_by_id('4HGLSDm')
        prev_lb, prev_ub = rxn.lower_bound, rxn.upper_bound
        rxn.lower_bound = 0
        rxn.upper_bound = 0
        rxn.lower_bound = prev_lb
        rxn.upper_bound = prev_ub
        assert rxn.lower_bound == prev_lb
        assert rxn.upper_bound == prev_ub

    def test_set_lb_higher_than_ub_sets_ub_to_new_lb(self, core_model):
        for reaction in core_model.reactions:
            assert reaction.lower_bound <= reaction.upper_bound
            reaction.lower_bound = reaction.upper_bound + 100
            assert reaction.lower_bound == reaction.upper_bound

    def test_set_ub_lower_than_lb_sets_lb_to_new_ub(self, core_model):
        for reaction in core_model.reactions:
            assert reaction.lower_bound <= reaction.upper_bound
            reaction.upper_bound = reaction.lower_bound - 100
            assert reaction.lower_bound == reaction.upper_bound

    def test_add_metabolites_combine_true(self, core_model):
        test_metabolite = Metabolite('test')
        for reaction in core_model.reactions:
            reaction.add_metabolites({test_metabolite: -66}, combine=True)
            assert reaction.metabolites[test_metabolite] == -66
            assert core_model.solver.constraints['test'].expression.has(-66. * reaction.forward_variable)
            assert core_model.solver.constraints['test'].expression.has(66. * reaction.reverse_variable)
            already_included_metabolite = list(reaction.metabolites.keys())[0]
            previous_coefficient = reaction.get_coefficient(
                already_included_metabolite.id)
            reaction.add_metabolites({already_included_metabolite: 10},
                                     combine=True)
            new_coefficient = previous_coefficient + 10
            assert reaction.metabolites[already_included_metabolite] == new_coefficient
            assert core_model.solver.constraints[already_included_metabolite.id].expression.has(
                new_coefficient * reaction.forward_variable)
            assert core_model.solver.constraints[already_included_metabolite.id].expression.has(
                -1 * new_coefficient * reaction.reverse_variable)

    @pytest.mark.skipif(TRAVIS, reason='non-deterministic')
    def test_add_metabolites_combine_false(self, core_model):
        test_metabolite = Metabolite('test')
        for reaction in core_model.reactions:
            reaction.add_metabolites({test_metabolite: -66}, combine=False)
            assert reaction.metabolites[test_metabolite] == -66
            assert core_model.solver.constraints['test'].expression.has(-66. * reaction.forward_variable)
            assert core_model.solver.constraints['test'].expression.has(66. * reaction.reverse_variable)
            already_included_metabolite = list(reaction.metabolites.keys())[0]
            reaction.add_metabolites({already_included_metabolite: 10}, combine=False)
            assert reaction.metabolites[already_included_metabolite] == 10
            assert core_model.solver.constraints[already_included_metabolite.id].expression.has(
                10 * reaction.forward_variable)
            assert core_model.solver.constraints[already_included_metabolite.id].expression.has(
                -10 * reaction.reverse_variable)

    # def test_pop(self, core_model):
    #     pgi = core_model.reactions.PGI
    #     g6p = core_model.metabolites.get_by_id("g6p_c")
    #     f6p = core_model.metabolites.get_by_id("f6p_c")
    #     g6p_expr = core_model.solver.constraints["g6p_c"].expression
    #     g6p_coef = pgi.pop("g6p_c")
    #     assert g6p not in pgi.metabolites
    #     actual = core_model.solver.constraints["g6p_c"].expression.as_coefficients_dict()
    #     expected = (g6p_expr - g6p_coef * pgi.flux_expression).as_coefficients_dict()
    #     assert actual == expected
    #     assert pgi.metabolites[f6p] == 1
    #
    #     f6p_expr = core_model.solver.constraints["f6p_c"].expression
    #     f6p_coef = pgi.pop(f6p)
    #     assert f6p not in pgi.metabolites
    #     assert core_model.solver.constraints["f6p_c"].expression.as_coefficients_dict() == (
    #         f6p_expr - f6p_coef * pgi.flux_expression
    #     ).as_coefficients_dict()

    def test_remove_from_model(self, core_model):
        pgi = core_model.reactions.PGI
        pgi.remove_from_model()
        assert pgi.model is None
        assert not ("PGI" in core_model.reactions)
        assert not (pgi.id in core_model.solver.variables)
        assert not (pgi.reverse_id in core_model.solver.variables)

    def test_delete(self, core_model):
        pgi = core_model.reactions.PGI
        pgi.delete()
        assert pgi.model is None
        assert not ("PGI" in core_model.reactions)
        assert not (pgi.id in core_model.solver.variables)
        assert not (pgi.reverse_id in core_model.solver.variables)

    def test_change_id_is_reflected_in_solver(self, core_model):
        for i, reaction in enumerate(core_model.reactions):
            old_reaction_id = reaction.id
            assert core_model.solver.variables[old_reaction_id].name == old_reaction_id
            assert old_reaction_id in core_model.solver.variables
            new_reaction_id = reaction.id + '_' + str(i)
            reaction.id = new_reaction_id
            assert reaction.id == new_reaction_id
            assert not (old_reaction_id in core_model.solver.variables)
            assert reaction.id in core_model.solver.variables
            assert reaction.reverse_id in core_model.solver.variables
            name = core_model.solver.variables[reaction.id].name
            assert name == reaction.id


class TestModel:
    # def test_model_is_subclassed(self, core_model):
    #     assert isinstance(core_model, cobra.Model)
    #     for reac in core_model.reactions:
    #         assert isinstance(reac, Reaction)
    #         for met in reac.metabolites:
    #             assert isinstance(met, Metabolite)
    #             assert met in core_model.metabolites
    #             assert met is core_model.metabolites.get_by_id(met.id)
    #         for gene in reac.genes:
    #             assert isinstance(gene, Gene)
    #             assert gene in core_model.genes
    #             assert gene is core_model.genes.get_by_id(gene.id)
    #
    #     for gene in core_model.genes:
    #         assert isinstance(gene, Gene)
    #         for reac in gene.reactions:
    #             assert isinstance(reac, Reaction)
    #             assert reac in core_model.reactions
    #             assert reac is core_model.reactions.get_by_id(reac.id)
    #
    #     for met in core_model.metabolites:
    #         assert isinstance(met, Metabolite)
    #         for reac in met.reactions:
    #             assert isinstance(reac, Reaction)
    #             assert reac in core_model.reactions
    #             assert reac is core_model.reactions.get_by_id(reac.id)

    def test_objective_coefficient_reflects_changed_objective(self, core_model):
        biomass_r = core_model.reactions.get_by_id('Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2')
        assert biomass_r.objective_coefficient == 1
        core_model.objective = "PGI"
        assert biomass_r.objective_coefficient == 0
        assert core_model.reactions.PGI.objective_coefficient == 1

    def test_change_objective_through_objective_coefficient(self, core_model):
        biomass_r = core_model.reactions.get_by_id('Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2')
        pgi = core_model.reactions.PGI
        pgi.objective_coefficient = 2
        coef_dict = core_model.solver.objective.expression.as_coefficients_dict()
        # Check that objective has been updated
        assert coef_dict[pgi.forward_variable] == 2
        assert coef_dict[pgi.reverse_variable] == -2
        # Check that original objective is still in there
        assert coef_dict[biomass_r.forward_variable] == 1
        assert coef_dict[biomass_r.reverse_variable] == -1

    # def test_model_from_other_model(self, core_model):
    #     core_model = Model(id_or_model=core_model)
    #     for reaction in core_model.reactions:
    #         assert reaction == core_model.reactions.get_by_id(reaction.id)

    def test_add_reactions(self, core_model):
        r1 = Reaction('r1')
        r1.add_metabolites({Metabolite('A'): -1, Metabolite('B'): 1})
        r1.lower_bound, r1.upper_bound = -999999., 999999.
        r2 = Reaction('r2')
        r2.add_metabolites(
            {Metabolite('A'): -1, Metabolite('C'): 1, Metabolite('D'): 1})
        r2.lower_bound, r2.upper_bound = 0., 999999.
        core_model.add_reactions([r1, r2])
        r2.objective_coefficient = 3.
        assert r2.objective_coefficient == 3.
        assert core_model.reactions[-2] == r1
        assert core_model.reactions[-1] == r2
        assert isinstance(core_model.reactions[-2].reverse_variable, core_model.solver.interface.Variable)
        coefficients_dict = core_model.solver.objective.expression.as_coefficients_dict()
        biomass_r = core_model.reactions.get_by_id('Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2')
        assert coefficients_dict[biomass_r.forward_variable] == 1.
        assert coefficients_dict[biomass_r.reverse_variable] == -1.
        assert coefficients_dict[core_model.reactions.r2.forward_variable] == 3.
        assert coefficients_dict[core_model.reactions.r2.reverse_variable] == -3.

    def test_remove_reactions_1(self, core_model):
        core_model.remove_reactions([core_model.reactions.PGI, core_model.reactions.PGK])
        assert "PGI" not in core_model.reactions
        assert "PGK" not in core_model.reactions
        assert "PGI" not in core_model.reactions
        assert "PGK" not in core_model.reactions

    def test_remove_reactions_2(self, core_model):
        reactions_to_remove = core_model.reactions[10:30]
        assert all([reaction.model is core_model for reaction in reactions_to_remove])
        assert all([core_model.reactions.get_by_id(reaction.id) == reaction for reaction in reactions_to_remove])

        core_model.remove_reactions(reactions_to_remove)
        assert all([reaction.model is None for reaction in reactions_to_remove])
        for reaction in reactions_to_remove:
            assert reaction.id not in list(core_model.solver.variables.keys())

        core_model.add_reactions(reactions_to_remove)
        for reaction in reactions_to_remove:
            assert reaction in core_model.reactions

    def test_remove_and_add_reactions(self, core_model):
        model_copy = core_model.copy()
        pgi, pgk = model_copy.reactions.PGI, model_copy.reactions.PGK
        model_copy.remove_reactions([pgi, pgk])
        assert "PGI" not in model_copy.reactions
        assert "PGK" not in model_copy.reactions
        assert "PGI" in core_model.reactions
        assert "PGK" in core_model.reactions
        model_copy.add_reactions([pgi, pgk])
        assert "PGI" in core_model.reactions
        assert "PGK" in core_model.reactions
        assert "PGI" in model_copy.reactions
        assert "PGK" in model_copy.reactions

    # def test_add_cobra_reaction(self, core_model):
    #     r = cobra.Reaction(id="c1")
    #     core_model.add_reaction(r)
    #     assert isinstance(core_model.reactions.c1, Reaction)

    def test_all_objects_point_to_all_other_correct_objects(self, core_model):
        for reaction in core_model.reactions:
            assert reaction.model == core_model
            for gene in reaction.genes:
                assert gene == core_model.genes.get_by_id(gene.id)
                assert gene.model == core_model
                for reaction2 in gene.reactions:
                    assert reaction2.model == core_model
                    assert reaction2 == core_model.reactions.get_by_id(reaction2.id)

            for metabolite in reaction.metabolites:
                assert metabolite.model == core_model
                assert metabolite == core_model.metabolites.get_by_id(metabolite.id)
                for reaction2 in metabolite.reactions:
                    assert reaction2.model == core_model
                    assert reaction2 == core_model.reactions.get_by_id(reaction2.id)

    def test_objects_point_to_correct_other_after_copy(self, core_model):
        for reaction in core_model.reactions:
            assert reaction.model == core_model
            for gene in reaction.genes:
                assert gene == core_model.genes.get_by_id(gene.id)
                assert gene.model == core_model
                for reaction2 in gene.reactions:
                    assert reaction2.model == core_model
                    assert reaction2 == core_model.reactions.get_by_id(reaction2.id)

            for metabolite in reaction.metabolites:
                assert metabolite.model == core_model
                assert metabolite == core_model.metabolites.get_by_id(metabolite.id)
                for reaction2 in metabolite.reactions:
                    assert reaction2.model == core_model
                    assert reaction2 == core_model.reactions.get_by_id(reaction2.id)

    def test_objective(self, core_model):
        obj = core_model.objective
        assert {var.name: coef for var, coef in obj.expression.as_coefficients_dict().items()} == \
               {'Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2_reverse_9ebcd': -1,
                'Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2': 1}
        assert obj.direction == "max"

    def test_change_objective(self, core_model):
        expression = 1.0 * core_model.solver.variables['ENO'] + 1.0 * core_model.solver.variables['PFK']
        core_model.solver.objective = core_model.solver.interface.Objective(expression)
        assert core_model.solver.objective.expression == expression

        core_model.objective = "ENO"
        eno_obj = core_model.solver.interface.Objective(
            core_model.reactions.ENO.flux_expression, direction="max")
        pfk_obj = core_model.solver.interface.Objective(
            core_model.reactions.PFK.flux_expression, direction="max")
        assert core_model.solver.objective == eno_obj

        with core_model:
            core_model.objective = "PFK"
            assert core_model.solver.objective == pfk_obj
        assert core_model.solver.objective == eno_obj

    def test_set_reaction_objective(self, core_model):
        core_model.objective = core_model.reactions.ACALD
        assert str(core_model.solver.objective.expression) == str(
            1.0 * core_model.reactions.ACALD.forward_variable -
            1.0 * core_model.reactions.ACALD.reverse_variable)

    def test_set_reaction_objective_str(self, core_model):
        core_model.objective = core_model.reactions.ACALD.id
        assert str(core_model.solver.objective.expression) == str(
            1.0 * core_model.reactions.ACALD.forward_variable -
            1.0 * core_model.reactions.ACALD.reverse_variable)

    def test_invalid_objective_raises(self, core_model):
        with pytest.raises(ValueError):
            core_model.objective = 'This is not a valid objective!'
        with pytest.raises(TypeError):
            setattr(core_model, 'objective', 3.)
    #
    # def test_solver_change(self, core_model):
    #     solver_id = id(core_model.solver)
    #     problem_id = id(core_model.solver.problem)
    #     solution = core_model.optimize().fluxes
    #     core_model.solver = 'glpk'
    #     assert id(core_model.solver) != solver_id
    #     assert id(core_model.solver.problem) != problem_id
    #     new_solution = core_model.optimize()
    #     for key in list(solution.keys()):
    #         assert round(abs(new_solution.fluxes[key] - solution[key]), 7) == 0
    #
    # def test_solver_change_with_optlang_interface(self, core_model):
    #     solver_id = id(core_model.solver)
    #     problem_id = id(core_model.solver.problem)
    #     solution = core_model.optimize().fluxes
    #     core_model.solver = optlang.glpk_interface
    #     assert id(core_model.solver) != solver_id
    #     assert id(core_model.solver.problem) != problem_id
    #     new_solution = core_model.optimize()
    #     for key in list(solution.keys()):
    #         assert round(abs(new_solution.fluxes[key] - solution[key]), 7) == 0

    def test_invalid_solver_change_raises(self, core_model):
        with pytest.raises(SolverNotFound):
            setattr(core_model, 'solver', [1, 2, 3])
        with pytest.raises(SolverNotFound):
            setattr(core_model, 'solver', 'ThisIsDefinitelyNotAvalidSolver')
        with pytest.raises(SolverNotFound):
            setattr(core_model, 'solver', os)

    @pytest.mark.skipif('cplex' not in solvers, reason='no cplex')
    def test_change_solver_to_cplex_and_check_copy_works(self, core_model):
        assert round(abs(core_model.slim_optimize() - 0.8739215069684306), 7) == 0
        core_model_copy = core_model.copy()
        assert round(abs(core_model_copy.slim_optimize() - 0.8739215069684306), 7) == 0
        # Second, change existing glpk based model to cplex
        core_model.solver = 'cplex'
        assert round(abs(core_model.slim_optimize() - 0.8739215069684306), 7) == 0
        core_model_copy = copy.copy(core_model)
        assert round(abs(core_model_copy.slim_optimize() - 0.8739215069684306), 7) == 0

    def test_copy_preserves_existing_solution(self, solved_model):
        solution, model = solved_model
        model_cp = copy.copy(model)
        primals_original = [variable.primal for variable in model.solver.variables]
        primals_copy = [variable.primal for variable in model_cp.solver.variables]
        abs_diff = abs(numpy.array(primals_copy) - numpy.array(primals_original))
        assert not any(abs_diff > 1e-6)

    def test_essential_genes(self, core_model):
        observed_essential_genes = [g.id for g in find_essential_genes(core_model)]
        assert sorted(observed_essential_genes) == sorted(ESSENTIAL_GENES)
        with pytest.raises(OptimizationError):
            core_model.reactions.Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2.lower_bound = 999999.
            find_essential_genes(core_model)

    def test_essential_reactions(self, core_model):
        observed_essential_reactions = [r.id for r in find_essential_reactions(core_model)]
        assert sorted(observed_essential_reactions) == sorted(ESSENTIAL_REACTIONS)
        with pytest.raises(OptimizationError):
            core_model.reactions.Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2.lower_bound = 999999.
            find_essential_reactions(core_model)

    def test_essential_metabolites_steady_state(self, core_model):
        essential_metabolites_balanced = [m.id for m in find_essential_metabolites(core_model,
                                                                                   force_steady_state=True)]
        assert sorted(essential_metabolites_balanced) == sorted(ESSENTIAL_METABOLITES)

        with pytest.raises(OptimizationError):
            core_model.reactions.Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2.lower_bound = 999999.
            find_essential_metabolites(core_model, force_steady_state=True)

    @pytest.mark.xfail(reason='needs some refactoring, uses missing bounds, not allowed by cplex')
    def test_essential_metabolites(self, core_model):
        essential_metabolites_unbalanced = [m.id for m in find_essential_metabolites(core_model,
                                                                                     force_steady_state=False)]
        assert sorted(essential_metabolites_unbalanced) == sorted(ESSENTIAL_METABOLITES)

        with pytest.raises(OptimizationError):
            core_model.reactions.Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2.lower_bound = 999999.
            find_essential_metabolites(core_model, force_steady_state=False)

    # def test_effective_bounds(self, core_model):
    #     core_model.reactions.Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2.lower_bound = 0.873921
    #     for reaction in core_model.reactions:
    #         assert abs(reaction.effective_lower_bound - REFERENCE_FVA_SOLUTION_ECOLI_CORE['lower_bound'][
    #             reaction.id]) < 0.000001
    #         assert abs(reaction.effective_upper_bound - REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound'][
    #             reaction.id]) < 0.000001

    # def test_add_ratio_constraint(self, solved_model):
    #     solution, model = solved_model
    #     assert round(abs(solution.objective_value - 0.873921506968), 7) == 0
    #     assert 2 * solution.fluxes['PGI'] != solution.fluxes['G6PDH2r']
    #     cp = model.copy()
    #     ratio_constr = cp.add_ratio_constraint(cp.reactions.PGI, cp.reactions.G6PDH2r, 0.5)
    #     assert ratio_constr.name == 'ratio_constraint_PGI_G6PDH2r'
    #     solution = cp.optimize()
    #     assert round(abs(solution.objective_value - 0.870407873712), 7) == 0
    #     assert round(abs(2 * solution.fluxes['PGI'] - solution.fluxes['G6PDH2r']), 7) == 0
    #     cp = model.copy()
    #
    #     ratio_constr = cp.add_ratio_constraint(cp.reactions.PGI, cp.reactions.G6PDH2r, 0.5)
    #     assert ratio_constr.name == 'ratio_constraint_PGI_G6PDH2r'
    #     solution = cp.optimize()
    #     assert round(abs(solution.objective_value - 0.870407873712), 7) == 0
    #     assert round(abs(2 * solution.fluxes['PGI'] - solution.fluxes['G6PDH2r']), 7) == 0
    #
    #     cp = model.copy()
    #     ratio_constr = cp.add_ratio_constraint('PGI', 'G6PDH2r', 0.5)
    #     assert ratio_constr.name == 'ratio_constraint_PGI_G6PDH2r'
    #     solution = cp.optimize()
    #     assert abs(solution.objective_value - 0.870407) < 1e-6
    #     assert abs(2 * solution.fluxes['PGI'] - solution.fluxes['G6PDH2r']) < 1e-6
    #
    #     cp = model.copy()
    #     ratio_constr = cp.add_ratio_constraint([cp.reactions.PGI, cp.reactions.ACALD],
    #                                            [cp.reactions.G6PDH2r, cp.reactions.ACONTa], 0.5)
    #     assert ratio_constr.name == 'ratio_constraint_PGI+ACALD_G6PDH2r+ACONTa'
    #     solution = cp.optimize()
    #     assert abs(solution.objective_value - 0.872959) < 1e-6
    #     assert abs((solution.fluxes['PGI'] + solution.fluxes['ACALD']) -
    #                0.5 * (solution.fluxes['G6PDH2r'] + solution.fluxes['ACONTa'])) < 1e-5

    def test_fix_objective_as_constraint(self, core_model):
        # with TimeMachine
        with core_model:
            fix_objective_as_constraint(core_model)
            constraint_name = core_model.solver.constraints[-1]
            assert core_model.solver.constraints[-1].expression - core_model.objective.expression == 0
        assert constraint_name not in core_model.solver.constraints
        # without TimeMachine
        fix_objective_as_constraint(core_model)
        constraint_name = core_model.solver.constraints[-1]
        assert core_model.solver.constraints[-1].expression - core_model.objective.expression == 0
        assert constraint_name in core_model.solver.constraints

    def test_get_reaction_for(self, core_model):
        with core_model:
            for r in core_model.reactions:
                assert isinstance(get_reaction_for(core_model, r.id), cobra.Reaction)
                assert isinstance(get_reaction_for(core_model, r), cobra.Reaction)
            for m in core_model.metabolites:
                assert isinstance(get_reaction_for(core_model, m.id), cobra.Reaction)
                assert isinstance(get_reaction_for(core_model, m), cobra.Reaction)

        with pytest.raises(TypeError):
            get_reaction_for(core_model, None)
        with pytest.raises(KeyError):
            get_reaction_for(core_model, "blablabla")
        with pytest.raises(KeyError):
            get_reaction_for(core_model, "accoa_lp_c_lp_", add=False)

    def test_stoichiometric_matrix(self, core_model):
        stoichiometric_matrix = create_stoichiometric_array(core_model)
        assert len(core_model.reactions) == stoichiometric_matrix.shape[1]
        assert len(core_model.metabolites) == stoichiometric_matrix.shape[0]

        for i, reaction in enumerate(core_model.reactions):
            for j, metabolite in enumerate(core_model.metabolites):
                if metabolite in reaction.metabolites:
                    coefficient = reaction.metabolites[metabolite]
                else:
                    coefficient = 0
                assert stoichiometric_matrix[j, i] == coefficient

    def test_set_medium(self, core_model):
        this_medium = medium(core_model)
        for reaction in core_model.exchanges:
            if reaction.lower_bound == 0:
                assert reaction.id not in this_medium.reaction_id.values
            if reaction.lower_bound < 0:
                assert reaction.id in this_medium.reaction_id.values
        load_medium(core_model, this_medium)
        for rid in medium(core_model).reaction_id:
            assert len(this_medium[this_medium.reaction_id == rid]) == 1

    def test_solver_change_preserves_non_metabolic_constraints(self, core_model):
        with core_model:
            constraint = core_model.problem.Constraint(core_model.reactions.PGK.flux_expression -
                                                       0.5 * core_model.reactions.PFK.flux_expression,
                                                       lb=0, ub=0)
            core_model.add_cons_vars(constraint)
            all_constraint_ids = core_model.solver.constraints.keys()
            assert all_constraint_ids[-1], 'ratio_constraint_PGK_PFK'
            resurrected = pickle.loads(pickle.dumps(core_model))
            assert resurrected.solver.constraints.keys() == all_constraint_ids


class TestMetabolite:
    def test_set_id(self, core_model):
        met = Metabolite("test")
        with pytest.raises(TypeError):
            setattr(met, 'id', 1)
        core_model.add_metabolites([met])
        with pytest.raises(ValueError):
            setattr(met, "id", 'g6p_c')
        met.id = "test2"
        assert "test2" in core_model.metabolites
        assert "test" not in core_model.metabolites

    def test_remove_from_model(self, core_model):
        met = core_model.metabolites.get_by_id("g6p_c")
        met.remove_from_model()
        assert not (met.id in core_model.metabolites)
        assert not (met.id in core_model.solver.constraints)

    @pytest.mark.xfail(reason='to be implemented in cobra')
    def test_notebook_repr(self):
        met = Metabolite(id="test", name="test metabolites", formula="CH4")
        expected = """
        <table>
            <tr>
                <td><strong>Id</strong></td><td>test</td>
            </tr>
            <tr>
                <td><strong>Name</strong></td><td>test metabolites</td>
            </tr>
            <tr>
                <td><strong>Formula</strong></td><td>CH4</td>
             </tr>
        </table>""".replace(' ', '')
        assert met._repr_html_().replace(' ', '') == expected
