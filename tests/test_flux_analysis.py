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

from __future__ import absolute_import

import copy
import os
import re

import numpy as np
import pandas
import pytest
from cobra import Metabolite, Reaction
from cobra.exceptions import OptimizationError
from cobra.flux_analysis import find_essential_reactions
from cobra.flux_analysis.parsimonious import add_pfba
from cobra.util import create_stoichiometric_matrix, fix_objective_as_constraint
from sympy import Add

from cameo.flux_analysis import remove_infeasible_cycles, structural
from cameo.flux_analysis.analysis import (find_blocked_reactions,
                                          flux_variability_analysis,
                                          phenotypic_phase_plane,
                                          fix_pfba_as_constraint)
from cameo.flux_analysis.simulation import fba, lmoma, moma, pfba, room
from cameo.flux_analysis.structural import nullspace
from cameo.parallel import MultiprocessingView, SequentialView
from cameo.util import current_solver_name, pick_one, ProblemCache

TRAVIS = 'TRAVIS' in os.environ
TEST_DIR = os.path.dirname(__file__)

REFERENCE_FVA_SOLUTION_ECOLI_CORE = pandas.read_csv(os.path.join(TEST_DIR, 'data/REFERENCE_flux_ranges_EcoliCore.csv'),
                                                    index_col=0)
REFERENCE_PPP_o2_EcoliCore = pandas.read_csv(os.path.join(TEST_DIR, 'data/REFERENCE_PPP_o2_EcoliCore.csv'))
REFERENCE_PPP_o2_EcoliCore_ac = pandas.read_csv(os.path.join(TEST_DIR, 'data/REFERENCE_PPP_o2_EcoliCore_ac.csv'))
REFERENCE_PPP_o2_glc_EcoliCore = pandas.read_csv(os.path.join(TEST_DIR, 'data/REFERENCE_PPP_o2_glc_EcoliCore.csv'))


def assert_data_frames_equal(obj, expected, delta=0.001, sort_by=None):
    df = obj.data_frame
    expected_names = [name for name in expected.columns.values if not re.match(r'^Unnamed.*', name)]
    df = df.fillna(value=0)
    expected = expected.fillna(value=0)
    if sort_by:
        df = df.sort_values(sort_by).reset_index(drop=True)
        expected = expected.sort_values(sort_by).reset_index(drop=True)
    for column in expected_names:
        for key in df.index:
            assert abs(df[column][key] - expected[column][key]) < delta


def test_find_blocked_reactions(core_model):
    core_model.reactions.PGK.knock_out()
    blocked_reactions = find_blocked_reactions(core_model)
    assert blocked_reactions == {core_model.reactions.GAPD, core_model.reactions.PGK}


class TestFluxVariabilityAnalysis:
    def test_flux_variability_parallel(self, core_model):
        original_objective = core_model.objective
        mp_view = MultiprocessingView(2)
        fva_solution = flux_variability_analysis(core_model, fraction_of_optimum=0.999999419892,
                                                 remove_cycles=False, view=mp_view)
        pfba_fva = flux_variability_analysis(core_model, fraction_of_optimum=1, pfba_factor=1,
                                             view=mp_view).data_frame
        mp_view.shutdown()
        assert_data_frames_equal(fva_solution, REFERENCE_FVA_SOLUTION_ECOLI_CORE)
        assert original_objective == core_model.objective
        assert sum(abs(pfba_fva.lower_bound)) - 518.422 < .001
        assert sum(abs(pfba_fva.upper_bound)) - 518.422 < .001

    def test_add_remove_pfba(self, core_model):
        with core_model:
            add_pfba(core_model)
            assert '_pfba_objective' == core_model.objective.name
        assert '_pfba_objective' != core_model.solver.constraints
        with core_model:
            fix_pfba_as_constraint(core_model)
            assert '_fixed_pfba_constraint' in core_model.solver.constraints
        assert '_fixed_pfba_constraint' not in core_model.solver.constraints

    @pytest.mark.skipif(TRAVIS, reason='Skip multiprocessing on Travis')
    def test_flux_variability_parallel_remove_cycles(self, core_model):
        original_objective = core_model.objective
        fva_solution = flux_variability_analysis(core_model, fraction_of_optimum=0.999999419892,
                                                 remove_cycles=True, view=MultiprocessingView())
        assert REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound']['FRD7'] > 666.
        assert round(abs(fva_solution['upper_bound']['FRD7'] - 0.), 7) == 0
        for key in fva_solution.data_frame.index:
            if REFERENCE_FVA_SOLUTION_ECOLI_CORE['lower_bound'][key] > -666:
                assert abs(
                    fva_solution['lower_bound'][key] - REFERENCE_FVA_SOLUTION_ECOLI_CORE['lower_bound'][key]) < 0.0001
            if REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound'][key] < 666:
                assert abs(
                    fva_solution['upper_bound'][key] - REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound'][key]) < 0.0001
        assert original_objective == core_model.objective

    def test_flux_variability_sequential(self, core_model):
        original_objective = core_model.objective
        fva_solution = flux_variability_analysis(core_model, fraction_of_optimum=0.999999419892,
                                                 remove_cycles=False, view=SequentialView())
        pfba_fva = flux_variability_analysis(core_model, fraction_of_optimum=1, pfba_factor=1).data_frame
        assert_data_frames_equal(fva_solution, REFERENCE_FVA_SOLUTION_ECOLI_CORE)
        assert_data_frames_equal(fva_solution, REFERENCE_FVA_SOLUTION_ECOLI_CORE)
        assert original_objective == core_model.objective
        assert sum(abs(pfba_fva.lower_bound)) - 518.422 < 0.001
        assert sum(abs(pfba_fva.upper_bound)) - 518.422 < 0.001

    def test_flux_variability_sequential_remove_cycles(self, core_model):
        original_objective = core_model.objective
        fva_solution = flux_variability_analysis(core_model, fraction_of_optimum=0.999999419892,
                                                 remove_cycles=True,
                                                 view=SequentialView())
        assert REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound']['FRD7'] > 666.
        assert round(abs(fva_solution['upper_bound']['FRD7'] - 0.), 7) == 0
        for key in fva_solution.data_frame.index:
            if REFERENCE_FVA_SOLUTION_ECOLI_CORE['lower_bound'][key] > -666:
                assert abs(
                    fva_solution['lower_bound'][key] - REFERENCE_FVA_SOLUTION_ECOLI_CORE['lower_bound'][key]) < 0.0001
            if REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound'][key] < 666:
                assert abs(
                    fva_solution['upper_bound'][key] - REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound'][key]) < 0.0001

        cycle_reac = Reaction("minus_PGI")  # Create fake cycle
        cycle_reac.lower_bound = -1000
        core_model.add_reaction(cycle_reac)
        cycle_reac.add_metabolites({met: -c for met, c in core_model.reactions.PGI.metabolites.items()})
        fva_solution = flux_variability_analysis(core_model, remove_cycles=False, reactions=["PGI"])
        assert fva_solution.data_frame.loc["PGI", "upper_bound"] == 1000
        fva_solution = flux_variability_analysis(core_model, remove_cycles=True, reactions=["PGI"])
        assert fva_solution.data_frame.loc["PGI", "upper_bound"] < 666
        assert original_objective == core_model.objective


class TestPhenotypicPhasePlane:

    @pytest.mark.skipif(TRAVIS, reason='Running in Travis')
    def test_one_variable_parallel(self, core_model):
        ppp = phenotypic_phase_plane(core_model, ['EX_o2_LPAREN_e_RPAREN_'], view=MultiprocessingView())
        assert_data_frames_equal(ppp, REFERENCE_PPP_o2_EcoliCore, sort_by=['EX_o2_LPAREN_e_RPAREN_'])
        ppp = phenotypic_phase_plane(core_model, 'EX_o2_LPAREN_e_RPAREN_', view=MultiprocessingView())
        assert_data_frames_equal(ppp, REFERENCE_PPP_o2_EcoliCore, sort_by=['EX_o2_LPAREN_e_RPAREN_'])

    def test_one_variable_sequential(self, core_model):
        ppp = phenotypic_phase_plane(core_model, ['EX_o2_LPAREN_e_RPAREN_'], view=SequentialView())
        assert_data_frames_equal(ppp, REFERENCE_PPP_o2_EcoliCore, sort_by=['EX_o2_LPAREN_e_RPAREN_'])
        ppp = phenotypic_phase_plane(core_model, 'EX_o2_LPAREN_e_RPAREN_', view=SequentialView())
        assert_data_frames_equal(ppp, REFERENCE_PPP_o2_EcoliCore, sort_by=['EX_o2_LPAREN_e_RPAREN_'])

    def test_one_variable_sequential_yield(self, core_model):
        ppp = phenotypic_phase_plane(core_model, ['EX_o2_LPAREN_e_RPAREN_'], view=SequentialView())
        assert_data_frames_equal(ppp, REFERENCE_PPP_o2_EcoliCore, sort_by=['EX_o2_LPAREN_e_RPAREN_'])
        ppp = phenotypic_phase_plane(core_model, 'EX_o2_LPAREN_e_RPAREN_', view=SequentialView())
        assert_data_frames_equal(ppp, REFERENCE_PPP_o2_EcoliCore, sort_by=['EX_o2_LPAREN_e_RPAREN_'])

    @pytest.mark.skipif(TRAVIS, reason='Running in Travis')
    def test_two_variables_parallel(self, core_model):
        ppp2d = phenotypic_phase_plane(core_model, ['EX_o2_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_'],
                                       view=MultiprocessingView())
        assert_data_frames_equal(ppp2d, REFERENCE_PPP_o2_glc_EcoliCore,
                                 sort_by=['EX_o2_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_'])

    def test_two_variables_sequential(self, core_model):
        ppp2d = phenotypic_phase_plane(core_model, ['EX_o2_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_'],
                                       view=SequentialView())
        assert_data_frames_equal(ppp2d, REFERENCE_PPP_o2_glc_EcoliCore,
                                 sort_by=['EX_o2_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_'])

    def test_one_variable_sequential_metabolite(self, core_model):
        ppp = phenotypic_phase_plane(core_model, ['EX_o2_LPAREN_e_RPAREN_'], core_model.metabolites.ac_c,
                                     view=SequentialView())
        assert_data_frames_equal(ppp, REFERENCE_PPP_o2_EcoliCore_ac, sort_by=['EX_o2_LPAREN_e_RPAREN_'])


class TestSimulationMethods:
    def test_fba(self, core_model):
        solution = fba(core_model)
        original_objective = core_model.objective
        assert abs(solution.objective_value - 0.873921) < 0.000001
        assert len(solution.fluxes) == len(core_model.reactions)
        assert core_model.objective.expression == original_objective.expression

    def test_fba_with_reaction_filter(self, core_model):
        original_objective = core_model.objective
        solution = fba(core_model, reactions=['EX_o2_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_'])
        assert abs(solution.objective_value - 0.873921) < 0.000001
        assert len(solution.fluxes) == 2
        assert core_model.objective.expression == original_objective.expression

    def test_pfba(self, core_model):
        original_objective = core_model.objective
        fba_solution = fba(core_model)
        fba_flux_sum = sum((abs(val) for val in list(fba_solution.fluxes.values)))
        pfba_solution = pfba(core_model)
        pfba_flux_sum = sum((abs(val) for val in list(pfba_solution.fluxes.values)))
        # looks like GLPK finds a parsimonious solution without the flux minimization objective
        assert (pfba_flux_sum - fba_flux_sum) < 1e-6, \
            "FBA sum is suppose to be lower than PFBA (was %f)" % (pfba_flux_sum - fba_flux_sum)
        assert core_model.objective.expression == original_objective.expression

    def test_pfba_with_reaction_filter(self, core_model):
        original_objective = core_model.objective
        pfba_solution = pfba(core_model, reactions=['EX_o2_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_'])
        assert len(pfba_solution.fluxes) == 2
        assert core_model.objective.expression == original_objective.expression

    def test_pfba_ijo1366(self, ijo1366):
        original_objective = ijo1366.objective
        fba_solution = fba(ijo1366)
        fba_flux_sum = sum((abs(val) for val in fba_solution.fluxes.values))
        pfba_solution = pfba(ijo1366)
        pfba_flux_sum = sum((abs(val) for val in pfba_solution.fluxes.values))
        assert (pfba_flux_sum - fba_flux_sum) < 1e-6, \
            "FBA sum is suppose to be lower than PFBA (was %f)" % (pfba_flux_sum - fba_flux_sum)
        assert ijo1366.objective.expression == original_objective.expression

    def test_lmoma(self, core_model):
        original_objective = core_model.objective
        pfba_solution = pfba(core_model)
        solution = lmoma(core_model, reference=pfba_solution)
        distance = sum((abs(solution[v] - pfba_solution[v]) for v in pfba_solution.keys()))
        assert abs(0 - distance) < 1e-6, "lmoma distance without knockouts must be 0 (was %f)" % distance
        assert core_model.objective.expression == original_objective.expression

    def test_lmoma_change_ref(self, core_model):
        original_objective = core_model.objective
        pfba_solution = pfba(core_model)
        fluxes = {rid: 10 * flux for rid, flux in pfba_solution.items()}
        solution = lmoma(core_model, reference=fluxes)
        distance = sum((abs(solution[v] - pfba_solution[v]) for v in pfba_solution.keys()))
        assert abs(0 - distance) > 1e-6, "lmoma distance without knockouts must be 0 (was %f)" % distance
        assert core_model.objective.expression == original_objective.expression
        assert not any(v.name.startswith("u_") for v in core_model.solver.variables)
        assert not any(c.name.startswith("lmoma_const_") for c in core_model.solver.constraints)

    def test_lmoma_with_reaction_filter(self, core_model):
        original_objective = core_model.objective
        pfba_solution = pfba(core_model)
        solution = lmoma(core_model, reference=pfba_solution,
                         reactions=['EX_o2_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_'])
        assert len(solution.fluxes) == 2
        assert core_model.objective.expression == original_objective.expression
        assert not any(v.name.startswith("u_") for v in core_model.solver.variables)
        assert not any(c.name.startswith("lmoma_const_") for c in core_model.solver.constraints)

    def test_moma(self, core_model):
        if current_solver_name(core_model) == 'glpk':
            pytest.skip('glpk does not support qp')
        original_objective = core_model.objective
        pfba_solution = pfba(core_model)
        solution = moma(core_model, reference=pfba_solution)
        distance = sum((abs(solution[v] - pfba_solution[v]) for v in pfba_solution.keys()))
        assert abs(0 - distance) < 1e-6, "moma distance without knockouts must be 0 (was %f)" % distance
        assert core_model.objective.expression == original_objective.expression
        assert not any(v.name.startswith("moma_aux_") for v in core_model.solver.variables)
        assert not any(c.name.startswith("moma_const_") for c in core_model.solver.constraints)

    def test_room(self, core_model):
        original_objective = core_model.objective
        pfba_solution = pfba(core_model)
        solution = room(core_model, reference=pfba_solution)
        assert abs(0 - solution.objective_value) < 1e-6, \
            "room objective without knockouts must be 0 (was %f)" % solution.objective_value
        assert core_model.objective.expression == original_objective.expression
        assert not any(v.name.startswith("y_") for v in core_model.solver.variables)
        assert not any(c.name.startswith("moma_const_") for c in core_model.solver.constraints)

    def test_room_with_reaction_filter(self, core_model):
        original_objective = core_model.objective
        pfba_solution = pfba(core_model)
        solution = room(core_model, reference=pfba_solution,
                        reactions=['EX_o2_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_'])
        assert len(solution.fluxes) == 2
        assert core_model.objective.expression == original_objective.expression
        assert not any(v.name.startswith("y_") for v in core_model.solver.variables)

    def test_moma_with_cache(self, core_model):
        if current_solver_name(core_model) == 'glpk':
            pytest.skip('glpk does not support qp')
        original_objective = core_model.objective
        pfba_solution = pfba(core_model)
        essential_reactions = find_essential_reactions(core_model)
        cache = ProblemCache(core_model)
        for r in core_model.reactions:
            if r not in essential_reactions:
                with core_model:
                    r.knock_out()
                    moma(core_model, reference=pfba_solution, cache=cache)
                    assert any(v.name.startswith("moma_aux_") for v in core_model.solver.variables)
                    assert any(c.name.startswith("moma_const_") for c in core_model.solver.constraints)
        cache.reset()
        assert core_model.objective.expression == original_objective.expression
        assert not any(v.name.startswith("moma_aux_") for v in core_model.solver.variables)
        assert not any(c.name.startswith("moma_const_") for c in core_model.solver.constraints)

    def test_lmoma_with_cache(self, core_model):
        original_objective = core_model.objective
        pfba_solution = pfba(core_model)
        essential_reactions = find_essential_reactions(core_model)
        cache = ProblemCache(core_model)
        for r in core_model.reactions:
            if r not in essential_reactions:
                with core_model:
                    r.knock_out()
                    lmoma(core_model, reference=pfba_solution, cache=cache)
                    assert any(v.name.startswith("u_") for v in core_model.solver.variables)
                    assert any(c.name.startswith("lmoma_const_") for c in core_model.solver.constraints)
        cache.reset()
        assert core_model.objective.expression == original_objective.expression
        assert not any(v.name.startswith("u_") for v in core_model.solver.variables)
        assert not any(c.name.startswith("lmoma_const_") for c in core_model.solver.constraints)

    def test_room_with_cache(self, core_model):
        original_objective = core_model.objective
        pfba_solution = pfba(core_model)
        essential_reactions = find_essential_reactions(core_model)
        cache = ProblemCache(core_model)
        infeasible = 0
        for r in core_model.reactions:
            if r not in essential_reactions:
                with core_model:
                    r.knock_out()
                    try:
                        room(core_model, reference=pfba_solution, cache=cache)
                        assert any(v.name.startswith("y_") for v in core_model.solver.variables)
                        assert any(c.name.startswith("room_const_") for c in core_model.solver.constraints)
                    except OptimizationError:  # TODO: room shouldn't return infeasible for non-essential reactions
                        infeasible += 1
                        continue
        assert infeasible < len(core_model.reactions)
        cache.reset()
        assert core_model.objective.expression == original_objective.expression
        assert not any(v.name.startswith("y_") for v in core_model.solver.variables)
        assert not any(c.name.startswith("room_const_") for c in core_model.solver.constraints)

    def test_room_shlomi_2005(self, toy_model):
        if current_solver_name(toy_model) == "glpk":
            pytest.xfail("this test doesn't work with glpk")
        original_objective = toy_model.objective
        reference = {"b1": 10, "v1": 10, "v2": 5, "v3": 0, "v4": 0, "v5": 0, "v6": 5, "b2": 5, "b3": 5}
        expected = {'b1': 10.0, 'b2': 5.0, 'b3': 5.0, 'v1': 10.0,
                    'v2': 5.0, 'v3': 0.0, 'v4': 5.0, 'v5': 5.0, 'v6': 0.0}
        assert not any(v.name.startswith("y_") for v in toy_model.solver.variables)

        with toy_model:
            toy_model.reactions.v6.knock_out()
            result = room(toy_model, reference=reference, delta=0, epsilon=0)

        for k in reference.keys():
            assert abs(expected[k] - result.fluxes[k]) < 0.1, "%s: %f | %f"
        assert toy_model.objective.expression == original_objective.expression
        assert not any(v.name.startswith("y_") for v in toy_model.variables)

    def test_moma_shlomi_2005(self, toy_model):
        if current_solver_name(toy_model) == 'glpk':
            pytest.skip('glpk does not support qp')

        original_objective = toy_model.objective
        reference = {"b1": 10, "v1": 10, "v2": 5, "v3": 0, "v4": 0, "v5": 0, "v6": 5, "b2": 5, "b3": 5}
        expected = {'b1': 8.8, 'b2': 4.4, 'b3': 4.4, 'v1': 8.8,
                    'v2': 3.1, 'v3': 1.3, 'v4': 4.4, 'v5': 3.1, 'v6': 0.0}

        with toy_model:
            toy_model.reactions.v6.knock_out()
            result = moma(toy_model, reference=reference)

        for k in reference.keys():
            assert abs(expected[k] - result.fluxes[k]) < 0.1, "%s: %f | %f"
        assert toy_model.objective.expression == original_objective.expression
        assert not any(v.name.startswith("u_") for v in toy_model.solver.variables)

    def test_moma_shlomi_2005_change_ref(self, toy_model):
        if current_solver_name(toy_model) == 'glpk':
            pytest.skip('glpk does not support qp')

        original_objective = toy_model.objective
        reference = {"b1": 10, "v1": 10, "v2": 5, "v3": 0, "v4": 0, "v5": 0, "v6": 5, "b2": 5, "b3": 5}
        expected = {'b1': 8.8, 'b2': 4.4, 'b3': 4.4, 'v1': 8.8,
                    'v2': 3.1, 'v3': 1.3, 'v4': 4.4, 'v5': 3.1, 'v6': 0.0}

        with toy_model:
            toy_model.reactions.v6.knock_out()
            result = moma(toy_model, reference=reference)

        for k in reference.keys():
            assert abs(expected[k] - result.fluxes[k]) < 0.1, "%s: %f | %f"
        assert toy_model.objective.expression == original_objective.expression
        assert not any(v.name.startswith("u_") for v in toy_model.solver.variables)

    # TODO: this test should be merged with the one above but problem cache is not resetting the model properly anymore.
    def test_moma_shlomi_2005_change_ref_1(self, toy_model):
        if current_solver_name(toy_model) == 'glpk':
            pytest.skip('glpk does not support qp')
        expected = {'b1': 8.8, 'b2': 4.4, 'b3': 4.4, 'v1': 8.8,
                    'v2': 3.1, 'v3': 1.3, 'v4': 4.4, 'v5': 3.1, 'v6': 0.0}

        reference_changed = {"b1": 5, "v1": 5, "v2": 5, "v3": 0, "v4": 0, "v5": 0, "v6": 5, "b2": 5, "b3": 5}
        with toy_model:
            toy_model.reactions.v6.knock_out()
            result_changed = moma(toy_model, reference=reference_changed)
        assert np.all([expected != result_changed.fluxes])
        assert not any(v.name.startswith("u_") for v in toy_model.solver.variables)


class TestRemoveCycles:
    def test_remove_cycles(self, core_model):
        with core_model:
            fix_objective_as_constraint(core_model)
            original_objective = copy.copy(core_model.objective)
            core_model.objective = core_model.solver.interface.Objective(
                Add(*core_model.solver.variables.values()), name='Max_all_fluxes')
            solution = core_model.optimize()
            assert abs(solution.to_frame().fluxes.abs().sum() - 2508.293334) < 1e-6
            fluxes = solution.fluxes
            core_model.objective = original_objective
        clean_fluxes = remove_infeasible_cycles(core_model, fluxes)
        assert abs(sum(abs(pandas.Series(clean_fluxes))) - 518.42208550050827) < 1e-6


class TestStructural:
    def test_find_blocked_reactions(self, core_model):
        assert "PGK" in core_model.reactions
        core_model.reactions.PGK.knock_out()  # there are no blocked reactions in EcoliCore
        blocked_reactions = structural.find_blocked_reactions_nullspace(core_model)
        assert len(blocked_reactions) == 0

    def test_find_dead_end_reactions(self, core_model):
        assert len(structural.find_dead_end_reactions(core_model)) == 0
        met1 = Metabolite("fake_metabolite_1")
        met2 = Metabolite("fake_metabolite_2")
        reac = Reaction("fake_reac")
        reac.add_metabolites({met1: -1, met2: 1})
        core_model.add_reaction(reac)
        assert structural.find_dead_end_reactions(core_model) == {reac}

    def test_find_coupled_reactions(self, core_model):
        couples = structural.find_coupled_reactions(core_model)
        fluxes = core_model.optimize().fluxes
        for coupled_set in couples:
            coupled_set = list(coupled_set)
            assert round(abs(fluxes[coupled_set[0].id] - fluxes[coupled_set[1].id]), 7) == 0

        couples, blocked = structural.find_coupled_reactions(core_model, return_dead_ends=True)
        assert blocked == structural.find_dead_end_reactions(core_model)

    @pytest.mark.skipif(TRAVIS, reason="ShortestElementaryFluxModes needs refactor")
    def test_shortest_elementary_flux_modes(self, core_model):
        if current_solver_name(core_model) == 'glpk':
            pytest.skip('sefm not supported for glpk')
        sefm = structural.ShortestElementaryFluxModes(core_model)
        ems = []
        for i, em in enumerate(sefm):
            if i > 10:
                break
            ems.append(em)
        assert list(map(len, ems)) == sorted(map(len, ems))

    def test_dead_end_metabolites_are_in_dead_end_reactions(self, core_model):
        dead_end_reactions = structural.find_dead_end_reactions(core_model)
        dead_end_metabolites = {m for m in core_model.metabolites if len(m.reactions) == 1}
        for dead_end_metabolite in dead_end_metabolites:
            assert any(dead_end_metabolite in r.metabolites for r in dead_end_reactions)

    def test_coupled_reactions(self, core_model):
        # If a reaction is essential, all coupled reactions are essential
        essential_reactions = find_essential_reactions(core_model)
        coupled_reactions = structural.find_coupled_reactions_nullspace(core_model)
        for essential_reaction in essential_reactions:
            for group in coupled_reactions:
                assert isinstance(group, dict)
                if essential_reaction in group:
                    assert all(group_reaction in essential_reactions for group_reaction in group)

    # # FIXME: this test has everything to run, but sometimes removing the reactions doesn't seem to work.
    # @pytest.mark.skipif(TRAVIS, reason="Inconsistent behaviour (bug)")
    def test_reactions_in_group_become_blocked_if_one_is_removed(self, core_model):
        essential_reactions = find_essential_reactions(core_model)
        coupled_reactions = structural.find_coupled_reactions_nullspace(core_model)
        for group in coupled_reactions:
            representative = pick_one(group)
            if representative not in essential_reactions:
                with core_model:
                    assert core_model == representative.model
                    core_model.remove_reactions([representative])
                    # # FIXME: Hack because of optlang queue issues with GLPK
                    # core_model.solver.update()
                    assert representative not in core_model.reactions
                    assert representative.forward_variable not in core_model.solver.variables
                    assert representative.reverse_variable not in core_model.solver.variables
                    assert representative not in core_model.reactions
                    assert representative.model is None
                    blocked_reactions = find_blocked_reactions(core_model)
                    assert all(r in blocked_reactions for r in group if r != representative)
                assert representative in core_model.reactions

        coupled_reactions = structural.find_coupled_reactions(core_model)
        for group in coupled_reactions:
            representative = pick_one(group)
            if representative not in essential_reactions:
                with core_model:
                    fwd_var_name = representative.forward_variable.name
                    rev_var_name = representative.reverse_variable.name
                    assert core_model == representative.model
                    core_model.remove_reactions([representative])
                    # # FIXME: Hack because of optlang queue issues with GLPK
                    # core_model.solver.update()
                    assert representative not in core_model.reactions
                    assert fwd_var_name not in core_model.solver.variables
                    assert rev_var_name not in core_model.solver.variables
                    assert representative not in core_model.reactions
                    assert representative.model is None
                    blocked_reactions = find_blocked_reactions(core_model)
                    assert representative not in core_model.reactions
                    assert all(r in blocked_reactions for r in group if r != representative)
                assert representative in core_model.reactions


class TestNullSpace:
    def test_wikipedia_toy(self):
        a = np.array([[2, 3, 5], [-4, 2, 3]])
        ns = nullspace(a)
        assert round(abs(np.dot(ns.T, a[0])[0] - 0), 10) == 0
        assert round(abs(np.dot(ns.T, a[1])[0] - 0), 10) == 0

    def test_with_core_model(self, core_model):
        s = create_stoichiometric_matrix(core_model)
        ns = nullspace(s)
        assert round(abs(np.abs(s.dot(ns)).max() - 0), 10) == 0
