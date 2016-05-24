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
import unittest

import pandas
from pandas.util.testing import assert_frame_equal
from sympy import Add

from cameo.flux_analysis import remove_infeasible_cycles
from cameo.flux_analysis.analysis import flux_variability_analysis, phenotypic_phase_plane, find_blocked_reactions
from cameo.flux_analysis.simulation import fba, pfba, lmoma, room, moma
from cameo.io import load_model
from cameo.parallel import SequentialView, MultiprocessingView
from cameo.util import TimeMachine

TRAVIS = os.getenv('TRAVIS', False)


def assert_dataframes_equal(df, expected):
    try:
        assert_frame_equal(df, expected, check_names=False)
        return True
    except AssertionError:
        return False


TESTDIR = os.path.dirname(__file__)
REFERENCE_FVA_SOLUTION_ECOLI_CORE = pandas.read_csv(os.path.join(TESTDIR, 'data/REFERENCE_flux_ranges_EcoliCore.csv'),
                                                    index_col=0)
REFERENCE_PPP_o2_EcoliCore = pandas.read_csv(os.path.join(TESTDIR, 'data/REFERENCE_PPP_o2_EcoliCore.csv'))
REFERENCE_PPP_o2_glc_EcoliCore = pandas.read_csv(os.path.join(TESTDIR, 'data/REFERENCE_PPP_o2_glc_EcoliCore.csv'))

CORE_MODEL = load_model(os.path.join(TESTDIR, 'data/EcoliCore.xml'), sanitize=False)
iJO_MODEL = load_model(os.path.join(TESTDIR, 'data/iJO1366.xml'), sanitize=False)

iJO_MODEL_COBRAPY = load_model(os.path.join(TESTDIR, 'data/iJO1366.xml'), solver_interface=None, sanitize=False)

TOY_MODEL_PAPIN_2003 = load_model(os.path.join(TESTDIR, "data/toy_model_Papin_2003.xml"))


class Wrapper:
    class AbstractTestFindBlockedReactions(unittest.TestCase):
        def test_find_blocked_reactions(self):
            self.model.reactions.PGK.knock_out()  # there are no blocked reactions in EcoliCore
            blocked_reactions = find_blocked_reactions(self.model)
            self.assertEqual(blocked_reactions, [self.model.reactions.GAPD, self.model.reactions.PGK])

    class AbstractTestFluxVariabilityAnalysis(unittest.TestCase):

        def test_flux_variability_parallel(self):
            mp_view = MultiprocessingView(2)
            fva_solution = flux_variability_analysis(self.model, remove_cycles=False, view=mp_view)
            mp_view.shutdown()
            assert_dataframes_equal(fva_solution, REFERENCE_FVA_SOLUTION_ECOLI_CORE)

        @unittest.skipIf(TRAVIS, 'Skip multiprocessing in Travis')
        def test_flux_variability_parallel_remove_cycles(self):
            fva_solution = flux_variability_analysis(self.model, fraction_of_optimum=0.999999419892, remove_cycles=True,
                                                     view=MultiprocessingView())
            self.assertGreater(REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound']['FRD7'], 666.)
            self.assertAlmostEqual(fva_solution['upper_bound']['FRD7'], 0.)
            for key in fva_solution.data_frame.index:
                if abs(REFERENCE_FVA_SOLUTION_ECOLI_CORE['lower_bound'][key]) > -666:
                    self.assertAlmostEqual(fva_solution['lower_bound'][key],
                                           REFERENCE_FVA_SOLUTION_ECOLI_CORE['lower_bound'][key], delta=0.0001)
                if abs(REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound'][key]) < 666:
                    self.assertAlmostEqual(fva_solution['upper_bound'][key],
                                           REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound'][key], delta=0.0001)

        def test_flux_variability_sequential(self):
            fva_solution = flux_variability_analysis(self.model, fraction_of_optimum=0.999999419892,
                                                     remove_cycles=False, view=SequentialView())
            assert_dataframes_equal(fva_solution, REFERENCE_FVA_SOLUTION_ECOLI_CORE)
            for key in fva_solution.data_frame.index:
                self.assertAlmostEqual(fva_solution['lower_bound'][key],
                                       REFERENCE_FVA_SOLUTION_ECOLI_CORE['lower_bound'][key], delta=0.00001)
                self.assertAlmostEqual(fva_solution['upper_bound'][key],
                                       REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound'][key], delta=0.00001)
            assert_dataframes_equal(fva_solution, REFERENCE_FVA_SOLUTION_ECOLI_CORE)

        def test_flux_variability_sequential_remove_cycles(self):
            fva_solution = flux_variability_analysis(self.model, fraction_of_optimum=0.999999419892, remove_cycles=True,
                                                     view=SequentialView())
            self.assertGreater(REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound']['FRD7'], 666.)
            self.assertAlmostEqual(fva_solution['upper_bound']['FRD7'], 0.)
            for key in fva_solution.data_frame.index:
                if abs(REFERENCE_FVA_SOLUTION_ECOLI_CORE['lower_bound'][key]) > -666:
                    self.assertAlmostEqual(fva_solution['lower_bound'][key],
                                           REFERENCE_FVA_SOLUTION_ECOLI_CORE['lower_bound'][key], delta=0.0001)
                if abs(REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound'][key]) < 666:
                    self.assertAlmostEqual(fva_solution['upper_bound'][key],
                                           REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound'][key], delta=0.0001)

    class AbstractTestPhenotypicPhasePlane(unittest.TestCase):
        @unittest.skipIf(TRAVIS, 'Running in Travis')
        def test_one_variable_parallel(self):
            ppp = phenotypic_phase_plane(self.model, ['EX_o2_LPAREN_e_RPAREN_'], view=MultiprocessingView())
            assert_dataframes_equal(ppp, REFERENCE_PPP_o2_EcoliCore)
            ppp = phenotypic_phase_plane(self.model, 'EX_o2_LPAREN_e_RPAREN_', view=MultiprocessingView())
            assert_dataframes_equal(ppp, REFERENCE_PPP_o2_EcoliCore)

        def test_one_variable_sequential(self):
            ppp = phenotypic_phase_plane(self.model, ['EX_o2_LPAREN_e_RPAREN_'], view=SequentialView())
            assert_dataframes_equal(ppp, REFERENCE_PPP_o2_EcoliCore)
            ppp = phenotypic_phase_plane(self.model, 'EX_o2_LPAREN_e_RPAREN_', view=SequentialView())
            assert_dataframes_equal(ppp, REFERENCE_PPP_o2_EcoliCore)

        @unittest.skipIf(TRAVIS, 'Running in Travis')
        def test_two_variables_parallel(self):
            ppp2d = phenotypic_phase_plane(self.model, ['EX_o2_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_'],
                                           view=MultiprocessingView())
            assert_dataframes_equal(ppp2d, REFERENCE_PPP_o2_glc_EcoliCore)

        def test_two_variables_sequential(self):
            ppp2d = phenotypic_phase_plane(self.model, ['EX_o2_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_'],
                                           view=SequentialView())
            assert_dataframes_equal(ppp2d, REFERENCE_PPP_o2_glc_EcoliCore)

    class AbstractTestSimulationMethods(unittest.TestCase):
        def setUp(self):
            self.model = CORE_MODEL

        def test_fba(self):
            solution = fba(self.model)
            original_objective = self.model.objective
            self.assertAlmostEqual(solution.objective_value, 0.873921, delta=0.000001)
            self.assertEqual(len(solution.fluxes), len(self.model.reactions))
            self.assertIs(self.model.objective, original_objective)

        def test_fba_with_reaction_filter(self):
            original_objective = self.model.objective
            solution = fba(self.model, reactions=['EX_o2_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_'])
            self.assertAlmostEqual(solution.objective_value, 0.873921, delta=0.000001)
            self.assertEqual(len(solution.fluxes), 2)
            self.assertIs(self.model.objective, original_objective)

        def test_pfba(self):
            original_objective = self.model.objective
            fba_solution = fba(self.model)
            fba_flux_sum = sum((abs(val) for val in list(fba_solution.fluxes.values())))
            pfba_solution = pfba(self.model)
            pfba_flux_sum = sum((abs(val) for val in list(pfba_solution.fluxes.values())))
            # looks like GLPK finds a parsimonious solution without the flux minimization objective
            self.assertTrue((pfba_flux_sum - fba_flux_sum) < 1e-6,
                            msg="FBA sum is suppose to be lower than PFBA (was %f)" % (pfba_flux_sum - fba_flux_sum))
            self.assertIs(self.model.objective, original_objective)

        def test_pfba_with_reaction_filter(self):
            original_objective = self.model.objective
            pfba_solution = pfba(self.model, reactions=['EX_o2_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_'])
            self.assertEqual(len(pfba_solution.fluxes), 2)
            self.assertIs(self.model.objective, original_objective)

        def test_pfba_iJO(self):
            original_objective = self.model.objective
            fba_solution = fba(iJO_MODEL)
            fba_flux_sum = sum((abs(val) for val in fba_solution.fluxes.values()))
            pfba_solution = pfba(iJO_MODEL)
            pfba_flux_sum = sum((abs(val) for val in pfba_solution.fluxes.values()))
            self.assertTrue((pfba_flux_sum - fba_flux_sum) < 1e-6,
                            msg="FBA sum is suppose to be lower than PFBA (was %f)" % (pfba_flux_sum - fba_flux_sum))
            self.assertIs(self.model.objective, original_objective)

        def test_lmoma(self):
            original_objective = self.model.objective
            pfba_solution = pfba(self.model)
            solution = lmoma(self.model, reference=pfba_solution)
            distance = sum((abs(solution[v] - pfba_solution[v]) for v in pfba_solution.keys()))
            self.assertAlmostEqual(0, distance,
                                   delta=1e-6,
                                   msg="lmoma distance without knockouts must be 0 (was %f)" % distance)
            self.assertIs(self.model.objective, original_objective)

        def test_lmoma_with_reaction_filter(self):
            original_objective = self.model.objective
            pfba_solution = pfba(self.model)
            solution = lmoma(self.model, reference=pfba_solution,
                             reactions=['EX_o2_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_'])
            self.assertEqual(len(solution.fluxes), 2)
            self.assertIs(self.model.objective, original_objective)

        def test_moma(self):
            original_objective = self.model.objective
            pfba_solution = pfba(self.model)
            solution = moma(self.model, reference=pfba_solution)
            distance = sum((abs(solution[v] - pfba_solution[v]) for v in pfba_solution.keys()))
            self.assertAlmostEqual(0, distance,
                                   delta=1e-6,
                                   msg="moma distance without knockouts must be 0 (was %f)" % distance)
            self.assertIs(self.model.objective, original_objective)

        def test_room(self):
            original_objective = self.model.objective
            pfba_solution = pfba(self.model)
            solution = room(self.model, reference=pfba_solution)
            self.assertAlmostEqual(0, solution.objective_value,
                                   delta=1e-6,
                                   msg="room objective without knockouts must be 0 (was %f)" % solution.objective_value)
            self.assertIs(self.model.objective, original_objective)

        def test_room_with_reaction_filter(self):
            original_objective = self.model.objective
            pfba_solution = pfba(self.model)
            solution = room(self.model, reference=pfba_solution,
                            reactions=['EX_o2_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_'])
            self.assertEqual(len(solution.fluxes), 2)
            self.assertIs(self.model.objective, original_objective)

        def test_room_shlomi_2005(self):
            original_objective = self.model.objective
            reference = {"b1": 10, "v1": 10, "v2": 5, "v3": 0, "v4": 0, "v5": 0, "v6": 5, "b2": 5, "b3": 5}
            expected = {'b1': 10.0, 'b2': 5.0, 'b3': 5.0, 'v1': 10.0,
                        'v2': 5.0, 'v3': 0.0, 'v4': 5.0, 'v5': 5.0, 'v6': 0.0}
            TOY_MODEL_PAPIN_2003.solver = self.model.solver.interface
            with TimeMachine() as tm:
                TOY_MODEL_PAPIN_2003.reactions.v6.knock_out(tm)
                result = room(TOY_MODEL_PAPIN_2003, reference=reference, delta=0, epsilon=0)

            for k in reference.keys():
                self.assertAlmostEqual(expected[k], result.fluxes[k], delta=0.1, msg="%s: %f | %f")
            self.assertIs(self.model.objective, original_objective)

        def test_moma_shlomi_2005(self):
            original_objective = self.model.objective
            reference = {"b1": 10, "v1": 10, "v2": 5, "v3": 0, "v4": 0, "v5": 0, "v6": 5, "b2": 5, "b3": 5}
            expected = {'b1': 8.8, 'b2': 4.4, 'b3': 4.4, 'v1': 8.8,
                        'v2': 3.1, 'v3': 1.3, 'v4': 4.4, 'v5': 3.1, 'v6': 0.0}

            TOY_MODEL_PAPIN_2003.solver = self.model.solver.interface
            with TimeMachine() as tm:
                TOY_MODEL_PAPIN_2003.reactions.v6.knock_out(tm)
                result = moma(TOY_MODEL_PAPIN_2003, reference=reference)

            for k in reference.keys():
                self.assertAlmostEqual(expected[k], result.fluxes[k], delta=0.1, msg="%s: %f | %f")
            self.assertIs(self.model.objective, original_objective)


class TestFindBlockedReactionsGLPK(Wrapper.AbstractTestFindBlockedReactions):
    def setUp(self):
        self.model = CORE_MODEL.copy()
        self.model.solver = 'glpk'


class TestFindBlockedReactionsCPLEX(Wrapper.AbstractTestFindBlockedReactions):
    def setUp(self):
        self.model = CORE_MODEL.copy()
        self.model.solver = 'cplex'


class TestFluxVariabilityAnalysisGLPK(Wrapper.AbstractTestFluxVariabilityAnalysis):
    def setUp(self):
        self.model = CORE_MODEL.copy()
        self.model.solver = 'glpk'


class TestFluxVariabilityAnalysisCPLEX(Wrapper.AbstractTestFluxVariabilityAnalysis):
    def setUp(self):
        self.model = CORE_MODEL.copy()
        self.model.solver = 'cplex'



class TestRemoveCycles(unittest.TestCase):
    def setUp(self):
        self.model = CORE_MODEL.copy()

    def test_remove_cyles(self):
        with TimeMachine() as tm:
            self.model.fix_objective_as_constraint(time_machine=tm)
            original_objective = copy.copy(self.model.objective)
            self.model.objective = self.model.solver.interface.Objective(Add(*self.model.solver.variables.values()),
                                                                         name='Max all fluxes')
            solution = self.model.solve()
            self.assertAlmostEqual(solution.data_frame.fluxes.abs().sum(), 2508.293334, delta=1e-6)
            self.model.objective = original_objective
            fluxes = solution.fluxes
        clean_fluxes = remove_infeasible_cycles(self.model, fluxes)
        self.assertAlmostEqual(pandas.Series(clean_fluxes).abs().sum(), 518.42208550050827, delta=1e-6)


class TestPhenotypicPhasePlaneGLPK(Wrapper.AbstractTestPhenotypicPhasePlane):
    def setUp(self):
        self.model = CORE_MODEL.copy()
        self.model.solver = 'glpk'

    def test_one_variable_sequential_metabolite(self):
        ppp = phenotypic_phase_plane(self.model, ['EX_o2_LPAREN_e_RPAREN_'], self.model.metabolites.o2_c,
                                     view=SequentialView())
        assert_dataframes_equal(ppp, REFERENCE_PPP_o2_EcoliCore)


class TestPhenotypicPhasePlaneCPLEX(Wrapper.AbstractTestPhenotypicPhasePlane):
    def setUp(self):
        self.model = CORE_MODEL.copy()
        self.model.solver = 'cplex'


class TestSimulationMethodsGLPK(Wrapper.AbstractTestSimulationMethods):
    def setUp(self):
        self.model = CORE_MODEL
        self.model.solver = 'glpk'

    def test_moma(self):
        self.assertRaises(ValueError, super(TestSimulationMethodsGLPK, self).test_moma)  # GLPK has no QP support

    def test_moma_shlomi_2005(self):
        self.assertRaises(ValueError, super(TestSimulationMethodsGLPK, self).test_moma)  # GLPK has no QP support


class TestSimulationMethodsCPLEX(Wrapper.AbstractTestSimulationMethods):
    def setUp(self):
        self.model = CORE_MODEL
        self.model.solver = 'cplex'


if __name__ == '__main__':
    import nose

    nose.runmodule()
