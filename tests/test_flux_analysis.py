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
import unittest
from functools import partial

import numpy as np
import pandas
from nose.tools import assert_almost_equal
from optlang.exceptions import IndicatorConstraintsNotSupported
from sympy import Add

import cameo
from cameo.config import solvers
from cameo.flux_analysis import remove_infeasible_cycles, fix_pfba_as_constraint
from cameo.flux_analysis import structural
from cameo.flux_analysis.analysis import flux_variability_analysis, phenotypic_phase_plane, find_blocked_reactions
from cameo.flux_analysis.simulation import fba, pfba, lmoma, room, moma, add_pfba
from cameo.flux_analysis.structural import nullspace
from cameo.io import load_model
from cameo.parallel import SequentialView, MultiprocessingView
from cameo.util import TimeMachine, pick_one


def assert_data_frames_equal(obj, expected, delta=0.001, sort_by=None, nan=0):
    df = obj.data_frame
    expected_names = [name for name in expected.columns.values if not re.match(r'^Unnamed.*', name)]
    df = df.fillna(value=0)
    expected = expected.fillna(value=0)
    if sort_by:
        df = df.sort_values(sort_by).reset_index(drop=True)
        expected = expected.sort_values(sort_by).reset_index(drop=True)
    for column in expected_names:
        for key in df.index:
            assert_almost_equal(df[column][key], expected[column][key], delta=delta)

TRAVIS = os.getenv('TRAVIS', False)
TEST_DIR = os.path.dirname(__file__)

REFERENCE_FVA_SOLUTION_ECOLI_CORE = pandas.read_csv(os.path.join(TEST_DIR, 'data/REFERENCE_flux_ranges_EcoliCore.csv'),
                                                    index_col=0)
REFERENCE_PPP_o2_EcoliCore = pandas.read_csv(os.path.join(TEST_DIR, 'data/REFERENCE_PPP_o2_EcoliCore.csv'))
REFERENCE_PPP_o2_EcoliCore_ac = pandas.read_csv(os.path.join(TEST_DIR, 'data/REFERENCE_PPP_o2_EcoliCore_ac.csv'))
REFERENCE_PPP_o2_glc_EcoliCore = pandas.read_csv(os.path.join(TEST_DIR, 'data/REFERENCE_PPP_o2_glc_EcoliCore.csv'))


class Wrapper:
    class CommonGround(unittest.TestCase):
        @classmethod
        def setUpClass(cls):
            cls.ecoli_core = load_model(os.path.join(TEST_DIR, 'data/EcoliCore.xml'), sanitize=False)
            cls.toy_model = load_model(os.path.join(TEST_DIR, "data/toy_model_Papin_2003.xml"))

        def __init__(self, *args, **kwargs):
            super(Wrapper.CommonGround, self).__init__(*args, **kwargs)
            self._iJO1366 = None

        @property
        def iJO1366(self):
            if self._iJO1366 is None:
                self._iJO1366 = load_model(os.path.join(TEST_DIR, 'data/iJO1366.xml'), sanitize=False)
            if self._iJO1366.solver.interface.__name__ != self.ecoli_core.solver.interface.__name__:
                self._iJO1366.solver = self.ecoli_core.solver.interface
            return self._iJO1366

    class AbstractTestFindBlockedReactions(CommonGround):
        def test_find_blocked_reactions(self):
            self.ecoli_core.reactions.PGK.knock_out()  # there are no blocked reactions in EcoliCore
            blocked_reactions = find_blocked_reactions(self.ecoli_core)
            self.assertEqual(blocked_reactions, {self.ecoli_core.reactions.GAPD, self.ecoli_core.reactions.PGK})

    class AbstractTestFluxVariabilityAnalysis(CommonGround):

        def test_flux_variability_parallel(self):
            original_objective = self.ecoli_core.objective
            mp_view = MultiprocessingView(2)
            fva_solution = flux_variability_analysis(self.ecoli_core, fraction_of_optimum=0.999999419892,
                                                     remove_cycles=False, view=mp_view)
            pfba_fva = flux_variability_analysis(self.ecoli_core, fraction_of_optimum=1, pfba_factor=1,
                                                 view=mp_view).data_frame
            mp_view.shutdown()
            assert_data_frames_equal(fva_solution, REFERENCE_FVA_SOLUTION_ECOLI_CORE)
            self.assertEqual(original_objective, self.ecoli_core.objective)
            self.assertAlmostEqual(sum(abs(pfba_fva.lower_bound)), 518.422, delta=0.001)
            self.assertAlmostEqual(sum(abs(pfba_fva.upper_bound)), 518.422, delta=0.001)

        def test_add_remove_pfb(self):
            with TimeMachine() as tm:
                add_pfba(self.ecoli_core, time_machine=tm)
                self.assertEquals('_pfba_objective', self.ecoli_core.objective.name)
            self.assertNotEqual('_pfba_objective', self.ecoli_core.solver.constraints)
            with TimeMachine() as tm:
                fix_pfba_as_constraint(self.ecoli_core, time_machine=tm)
                self.assertTrue('_fixed_pfba_constraint' in self.ecoli_core.solver.constraints)
            self.assertTrue('_fixed_pfba_constraint' not in self.ecoli_core.solver.constraints)


        @unittest.skipIf(TRAVIS, 'Skip multiprocessing in Travis')
        def test_flux_variability_parallel_remove_cycles(self):
            original_objective = self.ecoli_core.objective
            fva_solution = flux_variability_analysis(self.ecoli_core, fraction_of_optimum=0.999999419892,
                                                     remove_cycles=True, view=MultiprocessingView())
            self.assertGreater(REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound']['FRD7'], 666.)
            self.assertAlmostEqual(fva_solution['upper_bound']['FRD7'], 0.)
            for key in fva_solution.data_frame.index:
                if REFERENCE_FVA_SOLUTION_ECOLI_CORE['lower_bound'][key] > -666:
                    self.assertAlmostEqual(fva_solution['lower_bound'][key],
                                           REFERENCE_FVA_SOLUTION_ECOLI_CORE['lower_bound'][key], delta=0.0001)
                if REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound'][key] < 666:
                    self.assertAlmostEqual(fva_solution['upper_bound'][key],
                                           REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound'][key], delta=0.0001)
            self.assertEqual(original_objective, self.ecoli_core.objective)

        def test_flux_variability_sequential(self):
            original_objective = self.ecoli_core.objective
            fva_solution = flux_variability_analysis(self.ecoli_core, fraction_of_optimum=0.999999419892,
                                                     remove_cycles=False, view=SequentialView())
            pfba_fva = flux_variability_analysis(self.ecoli_core, fraction_of_optimum=1, pfba_factor=1).data_frame
            assert_data_frames_equal(fva_solution, REFERENCE_FVA_SOLUTION_ECOLI_CORE)
            self.assertEqual(original_objective, self.ecoli_core.objective)
            self.assertAlmostEqual(sum(abs(pfba_fva.lower_bound)), 518.422, delta=0.001)
            self.assertAlmostEqual(sum(abs(pfba_fva.upper_bound)), 518.422, delta=0.001)

        def test_flux_variability_sequential_remove_cycles(self):
            original_objective = self.ecoli_core.objective
            fva_solution = flux_variability_analysis(self.ecoli_core, fraction_of_optimum=0.999999419892, remove_cycles=True,
                                                     view=SequentialView())
            self.assertGreater(REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound']['FRD7'], 666.)
            self.assertAlmostEqual(fva_solution['upper_bound']['FRD7'], 0.)
            for key in fva_solution.data_frame.index:
                if REFERENCE_FVA_SOLUTION_ECOLI_CORE['lower_bound'][key] > -666:
                    self.assertAlmostEqual(fva_solution['lower_bound'][key],
                                           REFERENCE_FVA_SOLUTION_ECOLI_CORE['lower_bound'][key], delta=0.0001)
                if REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound'][key] < 666:
                    self.assertAlmostEqual(fva_solution['upper_bound'][key],
                                           REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound'][key], delta=0.0001)

            cycle_reac = cameo.Reaction("minus_PGI")  # Create fake cycle
            cycle_reac.lower_bound = -1000
            self.ecoli_core.add_reaction(cycle_reac)
            cycle_reac.add_metabolites({met: -c for met, c in self.ecoli_core.reactions.PGI.metabolites.items()})
            fva_solution = flux_variability_analysis(self.ecoli_core, remove_cycles=False, reactions=["PGI"])
            self.assertEqual(fva_solution.data_frame.loc["PGI", "upper_bound"], 1000)
            fva_solution = flux_variability_analysis(self.ecoli_core, remove_cycles=True, reactions=["PGI"])
            self.assertTrue(fva_solution.data_frame.loc["PGI", "upper_bound"] < 666)
            self.assertEqual(original_objective, self.ecoli_core.objective)

    class AbstractTestPhenotypicPhasePlane(CommonGround):
        @unittest.skipIf(TRAVIS, 'Running in Travis')
        def test_one_variable_parallel(self):
            ppp = phenotypic_phase_plane(self.ecoli_core, ['EX_o2_LPAREN_e_RPAREN_'], view=MultiprocessingView())
            assert_data_frames_equal(ppp, REFERENCE_PPP_o2_EcoliCore, sort_by=['EX_o2_LPAREN_e_RPAREN_'])
            ppp = phenotypic_phase_plane(self.ecoli_core, 'EX_o2_LPAREN_e_RPAREN_', view=MultiprocessingView())
            assert_data_frames_equal(ppp, REFERENCE_PPP_o2_EcoliCore, sort_by=['EX_o2_LPAREN_e_RPAREN_'])

        def test_one_variable_sequential(self):
            ppp = phenotypic_phase_plane(self.ecoli_core, ['EX_o2_LPAREN_e_RPAREN_'], view=SequentialView())
            assert_data_frames_equal(ppp, REFERENCE_PPP_o2_EcoliCore, sort_by=['EX_o2_LPAREN_e_RPAREN_'])
            ppp = phenotypic_phase_plane(self.ecoli_core, 'EX_o2_LPAREN_e_RPAREN_', view=SequentialView())
            assert_data_frames_equal(ppp, REFERENCE_PPP_o2_EcoliCore, sort_by=['EX_o2_LPAREN_e_RPAREN_'])

        def test_one_variable_sequential_yield(self):
            ppp = phenotypic_phase_plane(self.ecoli_core, ['EX_o2_LPAREN_e_RPAREN_'], view=SequentialView())
            assert_data_frames_equal(ppp, REFERENCE_PPP_o2_EcoliCore, sort_by=['EX_o2_LPAREN_e_RPAREN_'])
            ppp = phenotypic_phase_plane(self.ecoli_core, 'EX_o2_LPAREN_e_RPAREN_', view=SequentialView())
            assert_data_frames_equal(ppp, REFERENCE_PPP_o2_EcoliCore, sort_by=['EX_o2_LPAREN_e_RPAREN_'])

        @unittest.skipIf(TRAVIS, 'Running in Travis')
        def test_two_variables_parallel(self):
            ppp2d = phenotypic_phase_plane(self.ecoli_core, ['EX_o2_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_'],
                                           view=MultiprocessingView())
            assert_data_frames_equal(ppp2d, REFERENCE_PPP_o2_glc_EcoliCore,
                                     sort_by=['EX_o2_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_'])

        def test_two_variables_sequential(self):
            ppp2d = phenotypic_phase_plane(self.ecoli_core, ['EX_o2_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_'],
                                           view=SequentialView())
            assert_data_frames_equal(ppp2d, REFERENCE_PPP_o2_glc_EcoliCore,
                                     sort_by=['EX_o2_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_'])

    class AbstractTestSimulationMethods(CommonGround):
        def test_fba(self):
            solution = fba(self.ecoli_core)
            original_objective = self.ecoli_core.objective
            self.assertAlmostEqual(solution.objective_value, 0.873921, delta=0.000001)
            self.assertEqual(len(solution.fluxes), len(self.ecoli_core.reactions))
            self.assertIs(self.ecoli_core.objective, original_objective)

        def test_fba_with_reaction_filter(self):
            original_objective = self.ecoli_core.objective
            solution = fba(self.ecoli_core, reactions=['EX_o2_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_'])
            self.assertAlmostEqual(solution.objective_value, 0.873921, delta=0.000001)
            self.assertEqual(len(solution.fluxes), 2)
            self.assertIs(self.ecoli_core.objective, original_objective)

        def test_pfba(self):
            original_objective = self.ecoli_core.objective
            fba_solution = fba(self.ecoli_core)
            fba_flux_sum = sum((abs(val) for val in list(fba_solution.fluxes.values())))
            pfba_solution = pfba(self.ecoli_core)
            pfba_flux_sum = sum((abs(val) for val in list(pfba_solution.fluxes.values())))
            # looks like GLPK finds a parsimonious solution without the flux minimization objective
            self.assertTrue((pfba_flux_sum - fba_flux_sum) < 1e-6,
                            msg="FBA sum is suppose to be lower than PFBA (was %f)" % (pfba_flux_sum - fba_flux_sum))
            self.assertIs(self.ecoli_core.objective, original_objective)

        def test_pfba_with_reaction_filter(self):
            original_objective = self.ecoli_core.objective
            pfba_solution = pfba(self.ecoli_core, reactions=['EX_o2_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_'])
            self.assertEqual(len(pfba_solution.fluxes), 2)
            self.assertIs(self.ecoli_core.objective, original_objective)

        def test_pfba_iJO(self):
            original_objective = self.iJO1366.objective
            fba_solution = fba(self.iJO1366)
            fba_flux_sum = sum((abs(val) for val in fba_solution.fluxes.values()))
            pfba_solution = pfba(self.iJO1366)
            pfba_flux_sum = sum((abs(val) for val in pfba_solution.fluxes.values()))
            self.assertTrue((pfba_flux_sum - fba_flux_sum) < 1e-6,
                            msg="FBA sum is suppose to be lower than PFBA (was %f)" % (pfba_flux_sum - fba_flux_sum))
            self.assertIs(self.iJO1366.objective, original_objective)

        def test_lmoma(self):
            original_objective = self.ecoli_core.objective
            pfba_solution = pfba(self.ecoli_core)
            solution = lmoma(self.ecoli_core, reference=pfba_solution)
            distance = sum((abs(solution[v] - pfba_solution[v]) for v in pfba_solution.keys()))
            self.assertAlmostEqual(0, distance,
                                   delta=1e-6,
                                   msg="lmoma distance without knockouts must be 0 (was %f)" % distance)
            self.assertIs(self.ecoli_core.objective, original_objective)

        def test_lmoma_change_ref(self):
            original_objective = self.ecoli_core.objective
            pfba_solution = pfba(self.ecoli_core)
            fluxes = {rid: 10*flux for rid, flux in pfba_solution.items()}
            solution = lmoma(self.ecoli_core, reference=fluxes)
            distance = sum((abs(solution[v] - pfba_solution[v]) for v in pfba_solution.keys()))
            self.assertNotAlmostEqual(0, distance,
                                   delta=1e-6,
                                   msg="lmoma distance without knockouts must be 0 (was %f)" % distance)
            self.assertIs(self.ecoli_core.objective, original_objective)

        def test_lmoma_with_reaction_filter(self):
            original_objective = self.ecoli_core.objective
            pfba_solution = pfba(self.ecoli_core)
            solution = lmoma(self.ecoli_core, reference=pfba_solution,
                             reactions=['EX_o2_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_'])
            self.assertEqual(len(solution.fluxes), 2)
            self.assertIs(self.ecoli_core.objective, original_objective)

        def test_moma(self):
            original_objective = self.ecoli_core.objective
            pfba_solution = pfba(self.ecoli_core)
            solution = moma(self.ecoli_core, reference=pfba_solution)
            distance = sum((abs(solution[v] - pfba_solution[v]) for v in pfba_solution.keys()))
            self.assertAlmostEqual(0, distance,
                                   delta=1e-6,
                                   msg="moma distance without knockouts must be 0 (was %f)" % distance)
            self.assertIs(self.ecoli_core.objective, original_objective)

        def test_room(self):
            original_objective = self.ecoli_core.objective
            pfba_solution = pfba(self.ecoli_core)
            solution = room(self.ecoli_core, reference=pfba_solution)
            self.assertAlmostEqual(0, solution.objective_value,
                                   delta=1e-6,
                                   msg="room objective without knockouts must be 0 (was %f)" % solution.objective_value)
            self.assertIs(self.ecoli_core.objective, original_objective)

        def test_room_with_reaction_filter(self):
            original_objective = self.ecoli_core.objective
            pfba_solution = pfba(self.ecoli_core)
            solution = room(self.ecoli_core, reference=pfba_solution,
                            reactions=['EX_o2_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_'])
            self.assertEqual(len(solution.fluxes), 2)
            self.assertIs(self.ecoli_core.objective, original_objective)

        def test_room_shlomi_2005(self):
            original_objective = self.toy_model.objective
            reference = {"b1": 10, "v1": 10, "v2": 5, "v3": 0, "v4": 0, "v5": 0, "v6": 5, "b2": 5, "b3": 5}
            expected = {'b1': 10.0, 'b2': 5.0, 'b3': 5.0, 'v1': 10.0,
                        'v2': 5.0, 'v3': 0.0, 'v4': 5.0, 'v5': 5.0, 'v6': 0.0}
            with TimeMachine() as tm:
                self.toy_model.reactions.v6.knock_out(tm)
                result = room(self.toy_model, reference=reference, delta=0, epsilon=0)

            for k in reference.keys():
                self.assertAlmostEqual(expected[k], result.fluxes[k], delta=0.1, msg="%s: %f | %f")
            self.assertIs(self.toy_model.objective, original_objective)

        def test_moma_shlomi_2005(self):
            original_objective = self.toy_model.objective
            reference = {"b1": 10, "v1": 10, "v2": 5, "v3": 0, "v4": 0, "v5": 0, "v6": 5, "b2": 5, "b3": 5}
            expected = {'b1': 8.8, 'b2': 4.4, 'b3': 4.4, 'v1': 8.8,
                        'v2': 3.1, 'v3': 1.3, 'v4': 4.4, 'v5': 3.1, 'v6': 0.0}

            with TimeMachine() as tm:
                self.toy_model.reactions.v6.knock_out(tm)
                result = moma(self.toy_model, reference=reference)

            for k in reference.keys():
                self.assertAlmostEqual(expected[k], result.fluxes[k], delta=0.1, msg="%s: %f | %f")
            self.assertIs(self.toy_model.objective, original_objective)

        def test_moma_shlomi_2005_change_ref(self):
            original_objective = self.toy_model.objective
            reference = {"b1": 10, "v1": 10, "v2": 5, "v3": 0, "v4": 0, "v5": 0, "v6": 5, "b2": 5, "b3": 5}
            expected = {'b1': 8.8, 'b2': 4.4, 'b3': 4.4, 'v1': 8.8,
                        'v2': 3.1, 'v3': 1.3, 'v4': 4.4, 'v5': 3.1, 'v6': 0.0}

            with TimeMachine() as tm:
                self.toy_model.reactions.v6.knock_out(tm)
                result = moma(self.toy_model, reference=reference)

            for k in reference.keys():
                self.assertAlmostEqual(expected[k], result.fluxes[k], delta=0.1, msg="%s: %f | %f")
            self.assertIs(self.toy_model.objective, original_objective)

            reference_changed = {"b1": 5, "v1": 5, "v2": 5, "v3": 0, "v4": 0, "v5": 0, "v6": 5, "b2": 5, "b3": 5}
            with TimeMachine() as tm:
                self.toy_model.reactions.v6.knock_out(tm)
                result_changed = moma(self.toy_model, reference=reference_changed)

            self.assertNotEqual(expected, result_changed.fluxes)

    class AbstractTestRemoveCycles(CommonGround):
        def test_remove_cyles(self):
            model = self.ecoli_core
            with TimeMachine() as tm:
                model.fix_objective_as_constraint(time_machine=tm)
                original_objective = copy.copy(self.ecoli_core.objective)
                model.objective = self.ecoli_core.solver.interface.Objective(Add(*model.solver.variables.values()),
                                                                             name='Max all fluxes')
                solution = model.solve()
                self.assertAlmostEqual(solution.data_frame.fluxes.abs().sum(), 2508.293334, delta=1e-6)
                fluxes = solution.fluxes
                model.objective = original_objective
            clean_fluxes = remove_infeasible_cycles(model, fluxes)
            self.assertAlmostEqual(pandas.Series(clean_fluxes).abs().sum(), 518.42208550050827, delta=1e-6)

    class AbstractTestStructural(CommonGround):
        def test_find_blocked_reactions(self):
            self.assertIn("PGK", self.ecoli_core.reactions)
            self.ecoli_core.reactions.PGK.knock_out()  # there are no blocked reactions in EcoliCore
            blocked_reactions = structural.find_blocked_reactions_nullspace(self.ecoli_core)
            self.assertEqual(len(blocked_reactions), 0)

        def test_find_dead_end_reactions(self):
            self.assertEqual(len(structural.find_dead_end_reactions(self.ecoli_core)), 0)
            met1 = cameo.Metabolite("fake_metabolite_1")
            met2 = cameo.Metabolite("fake_metabolite_2")
            reac = cameo.Reaction("fake_reac")
            reac.add_metabolites({met1: -1, met2: 1})
            self.ecoli_core.add_reaction(reac)
            self.assertEqual(structural.find_dead_end_reactions(self.ecoli_core), {reac})

        def test_find_coupled_reactions(self):
            couples = structural.find_coupled_reactions(self.ecoli_core)
            fluxes = self.ecoli_core.solve().fluxes
            for coupled_set in couples:
                coupled_set = list(coupled_set)
                self.assertAlmostEqual(fluxes[coupled_set[0].id], fluxes[coupled_set[1].id])

            couples, blocked = structural.find_coupled_reactions(self.ecoli_core, return_dead_ends=True)
            self.assertEqual(blocked, structural.find_dead_end_reactions(self.ecoli_core))

        @unittest.skipIf(TRAVIS, "ShortestElementaryFluxModes needs refactor")
        def test_shortest_elementary_flux_modes(self):
            sefm = structural.ShortestElementaryFluxModes(self.ecoli_core)
            ems = []
            for i, em in enumerate(sefm):
                if i > 10:
                    break
                ems.append(em)
            self.assertEqual(list(map(len, ems)), sorted(map(len, ems)))

        def test_dead_end_metabolites_are_in_dead_end_reactions(self):
            dead_end_reactions = structural.find_dead_end_reactions(self.ecoli_core)
            dead_end_metabolites = {m for m in self.ecoli_core.metabolites if len(m.reactions) == 1}
            for dead_end_metabolite in dead_end_metabolites:
                self.assertTrue(any(dead_end_metabolite in r.metabolites for r in dead_end_reactions))

        def test_coupled_reactions(self):
            # If a reaction is essential, all coupled reactions are essential
            essential_reactions = self.ecoli_core.essential_reactions()
            coupled_reactions = structural.find_coupled_reactions_nullspace(self.ecoli_core)
            for essential_reaction in essential_reactions:
                for group in coupled_reactions:
                    self.assertIsInstance(group, frozenset)
                    if essential_reaction in group:
                        self.assertTrue(all(group_reaction in essential_reactions for group_reaction in group))

        # FIXME: this test has everything to run, but sometimes removing the reactions doesn't seem to work.
        @unittest.skipIf(TRAVIS, "Inconsistent behaviour")
        def test_reactions_in_group_become_blocked_if_one_is_removed(self):
            essential_reactions = self.ecoli_core.essential_reactions()
            coupled_reactions = structural.find_coupled_reactions_nullspace(self.ecoli_core)
            for group in coupled_reactions:
                representative = pick_one(group)
                if representative not in essential_reactions:
                    with TimeMachine() as tm:
                        self.assertEqual(self.ecoli_core, representative.model)
                        tm(do=partial(self.ecoli_core.remove_reactions, [representative], delete=False),
                           undo=partial(self.ecoli_core.add_reactions, [representative]))
                        # FIXME: Hack because of optlang queue issues with GLPK
                        self.ecoli_core.solver.update()
                        self.assertNotIn(representative, self.ecoli_core.reactions)
                        self.assertNotIn(representative.forward_variable, self.ecoli_core.solver.variables)
                        self.assertNotIn(representative.reverse_variable, self.ecoli_core.solver.variables)
                        self.assertNotIn(representative, self.ecoli_core.reactions)
                        self.assertEqual(representative.model, None)
                        blocked_reactions = find_blocked_reactions(self.ecoli_core)
                        self.assertTrue(all(r in blocked_reactions for r in group if r != representative))
                    self.assertIn(representative, self.ecoli_core.reactions)

            coupled_reactions = structural.find_coupled_reactions(self.ecoli_core)
            for group in coupled_reactions:
                representative = pick_one(group)
                if representative not in essential_reactions:
                    with TimeMachine() as tm:
                        fwd_var_name = representative.forward_variable.name
                        rev_var_name = representative.reverse_variable.name
                        self.assertEqual(self.ecoli_core, representative.model)
                        tm(do=partial(self.ecoli_core.remove_reactions, [representative], delete=False),
                           undo=partial(self.ecoli_core.add_reactions, [representative]))
                        # FIXME: Hack because of optlang queue issues with GLPK
                        self.ecoli_core.solver.update()
                        self.assertNotIn(representative, self.ecoli_core.reactions)
                        self.assertNotIn(fwd_var_name, self.ecoli_core.solver.variables)
                        self.assertNotIn(rev_var_name, self.ecoli_core.solver.variables)
                        self.assertNotIn(representative, self.ecoli_core.reactions)
                        self.assertEqual(representative.model, None)
                        blocked_reactions = find_blocked_reactions(self.ecoli_core)
                        self.assertNotIn(representative, self.ecoli_core.reactions)
                        self.assertTrue(all(r in blocked_reactions for r in group if r != representative))
                    self.assertIn(representative, self.ecoli_core.reactions)


class TestFindBlockedReactionsGLPK(Wrapper.AbstractTestFindBlockedReactions):
    def setUp(self):
        self.ecoli_core.solver = 'glpk'
        self.toy_model.solver = 'glpk'


@unittest.skipIf('cplex' not in solvers, "No cplex interface available")
class TestFindBlockedReactionsCPLEX(Wrapper.AbstractTestFindBlockedReactions):
    def setUp(self):
        self.ecoli_core.solver = 'cplex'
        self.toy_model.solver = 'cplex'


class TestFluxVariabilityAnalysisGLPK(Wrapper.AbstractTestFluxVariabilityAnalysis):
    def setUp(self):
        self.ecoli_core.solver = 'glpk'
        self.toy_model.solver = 'glpk'


@unittest.skipIf('cplex' not in solvers, "No cplex interface available")
class TestFluxVariabilityAnalysisCPLEX(Wrapper.AbstractTestFluxVariabilityAnalysis):
    def setUp(self):
        self.ecoli_core.solver = 'cplex'
        self.toy_model.solver = 'cplex'


class TestPhenotypicPhasePlaneGLPK(Wrapper.AbstractTestPhenotypicPhasePlane):
    def setUp(self):
        self.ecoli_core.solver = 'glpk'
        self.toy_model.solver = 'glpk'

    def test_one_variable_sequential_metabolite(self):
        ppp = phenotypic_phase_plane(self.ecoli_core, ['EX_o2_LPAREN_e_RPAREN_'], self.ecoli_core.metabolites.ac_c,
                                     view=SequentialView())
        assert_data_frames_equal(ppp, REFERENCE_PPP_o2_EcoliCore_ac, sort_by=['EX_o2_LPAREN_e_RPAREN_'])


@unittest.skipIf('cplex' not in solvers, "No cplex interface available")
class TestPhenotypicPhasePlaneCPLEX(Wrapper.AbstractTestPhenotypicPhasePlane):
    def setUp(self):
        self.ecoli_core.solver = 'cplex'
        self.toy_model.solver = 'cplex'


class TestSimulationMethodsGLPK(Wrapper.AbstractTestSimulationMethods):
    def setUp(self):
        self.ecoli_core.solver = 'glpk'
        self.toy_model.solver = 'glpk'

    def test_moma(self):
        self.assertRaises(ValueError, super(TestSimulationMethodsGLPK, self).test_moma)  # GLPK has no QP support

    def test_moma_shlomi_2005(self):
        self.assertRaises(ValueError, super(TestSimulationMethodsGLPK, self).test_moma)  # GLPK has no QP support

    def test_moma_shlomi_2005_change_ref(self):
        self.assertRaises(ValueError, super(TestSimulationMethodsGLPK, self).test_moma)  # GLPK has no QP support


@unittest.skipIf('cplex' not in solvers, "No cplex interface available")
class TestSimulationMethodsCPLEX(Wrapper.AbstractTestSimulationMethods):
    def setUp(self):
        self.ecoli_core.solver = 'cplex'
        self.toy_model.solver = 'cplex'


class TestStructuralMethodsGLPK(Wrapper.AbstractTestStructural):
    def setUp(self):
        self.ecoli_core.solver = 'glpk'
        self.toy_model.solver = 'glpk'

    def test_shortest_elementary_flux_modes(self):
        self.assertRaises(IndicatorConstraintsNotSupported, structural.ShortestElementaryFluxModes, self.ecoli_core)


@unittest.skipIf('cplex' not in solvers, "No cplex interface available")
class TestStructuralMethodsCPLEX(Wrapper.AbstractTestStructural):
    def setUp(self):
        self.ecoli_core.solver = 'cplex'
        self.toy_model.solver = 'cplex'


class TestRemoveCyclesGLPK(Wrapper.AbstractTestRemoveCycles):
    def setUp(self):
        self.ecoli_core.solver = 'glpk'
        self.toy_model.solver = 'glpk'


@unittest.skipIf('cplex' not in solvers, "No cplex interface available")
class TestRemoveCyclesCPLEX(Wrapper.AbstractTestRemoveCycles):
    def setUp(self):
        self.ecoli_core.solver = 'cplex'
        self.toy_model.solver = 'cplex'


class NullSpaceTestCase(Wrapper.CommonGround):
    # toy from https://en.wikipedia.org/wiki/Kernel_(linear_algebra)#Illustration
    def test_wikipedia_toy(self):
        A = np.array([[2, 3, 5], [-4, 2, 3]])
        NS = nullspace(A)
        self.assertAlmostEqual(np.dot(NS.T, A[0])[0], 0, places=10)
        self.assertAlmostEqual(np.dot(NS.T, A[1])[0], 0, places=10)

    def test_with_core_model(self):
        S = self.ecoli_core.S
        NS = nullspace(S)
        self.assertAlmostEqual(np.abs(S.dot(NS)).max(), 0, places=10)


if __name__ == '__main__':
    import nose

    nose.runmodule()
