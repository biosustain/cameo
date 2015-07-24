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
import unittest

from cameo.flux_analysis.simulation import fba, pfba, lmoma, room
from cameo.parallel import SequentialView, MultiprocessingView
from cameo.io import load_model
from cameo.flux_analysis.analysis import flux_variability_analysis, phenotypic_phase_plane, _cycle_free_fva, \
    find_blocked_reactions

import pandas
from pandas.util.testing import assert_frame_equal
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

TOY_MODEL_PAPIN_2004 = load_model(os.path.join(TESTDIR, "data/toy_model_Papin_2003.xml"))


class AbstractTestFindBlockedReactions(object):
    def test_find_blocked_reactions(self):
        self.model.reactions.PGK.knock_out()  # there are no blocked reactions in EcoliCore
        blocked_reactions = find_blocked_reactions(self.model)
        self.assertEqual(blocked_reactions, [self.model.reactions.GAPD, self.model.reactions.PGK])


class TestFindBlockedReactionsGLPK(AbstractTestFindBlockedReactions, unittest.TestCase):
    def setUp(self):
        self.model = CORE_MODEL.copy()
        self.model.solver = 'glpk'


@unittest.skipIf(TRAVIS, 'CPLEX not available on Travis.')
class TestFindBlockedReactionsCPLEX(AbstractTestFindBlockedReactions, unittest.TestCase):
    def setUp(self):
        self.model = CORE_MODEL.copy()
        self.model.solver = 'cplex'
        self.model.solver.problem.parameters.lpmethod = self.model.solver.problem.parameters.lpmethod.values.dual


class AbstractTestFluxVariabilityAnalysis(object):
    def test_flux_variability_sequential(self):
        fva_solution = flux_variability_analysis(self.model, remove_cycles=False, view=SequentialView())
        assert_dataframes_equal(fva_solution, REFERENCE_FVA_SOLUTION_ECOLI_CORE)
        print(REFERENCE_FVA_SOLUTION_ECOLI_CORE)
        for key in fva_solution.data_frame.index:
            self.assertAlmostEqual(fva_solution['lower_bound'][key],
                                   REFERENCE_FVA_SOLUTION_ECOLI_CORE['lower_bound'][key], delta=0.000001)
            self.assertAlmostEqual(fva_solution['upper_bound'][key],
                                   REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound'][key], delta=0.000001)
        assert_dataframes_equal(fva_solution, REFERENCE_FVA_SOLUTION_ECOLI_CORE)

    @unittest.skip('Removing cycles is still not robust.')
    def test_flux_variability_sequential_remove_cycles(self):
        fva_solution = flux_variability_analysis(self.model, remove_cycles=True, view=SequentialView())
        print(fva_solution.min())
        print(fva_solution.max())
        for key in fva_solution.data_frame.index:
            if abs(REFERENCE_FVA_SOLUTION_ECOLI_CORE['lower_bound'][key]) < 999993:
                self.assertAlmostEqual(fva_solution['lower_bound'][key],
                                       REFERENCE_FVA_SOLUTION_ECOLI_CORE['lower_bound'][key], delta=0.000001)
            if abs(REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound'][key]) < 999993:
                print(key)
                self.assertAlmostEqual(fva_solution['upper_bound'][key],
                                       REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound'][key], delta=0.000001)

    @unittest.skipIf(TRAVIS, 'Running multiprocess in Travis breaks')
    def test_flux_variability_parallel(self):
        mp_view = MultiprocessingView()
        fva_solution = flux_variability_analysis(self.model, remove_cycles=False, view=mp_view)
        mp_view.shutdown()
        assert_dataframes_equal(fva_solution, REFERENCE_FVA_SOLUTION_ECOLI_CORE)

    @unittest.skip("multiprocessing doesn't work with cycle_free_fva yet")
    def test_flux_variability_parallel_remove_cycles(self):
        fva_solution = flux_variability_analysis(self.model, remove_cycles=True, view=MultiprocessingView())
        for key in fva_solution.index:
            if abs(REFERENCE_FVA_SOLUTION_ECOLI_CORE['lower_bound'][key]) < 999993:
                self.assertAlmostEqual(fva_solution['lower_bound'][key],
                                       REFERENCE_FVA_SOLUTION_ECOLI_CORE['lower_bound'][key], delta=0.000001)
            if abs(REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound'][key]) < 999993:
                self.assertAlmostEqual(fva_solution['upper_bound'][key],
                                       REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound'][key], delta=0.000001)


class TestFluxVariabilityAnalysisGLPK(AbstractTestFluxVariabilityAnalysis, unittest.TestCase):
    def setUp(self):
        self.model = CORE_MODEL.copy()
        self.model.solver = 'glpk'
        self.biomass_flux = 0.873921
        self.model.reactions.Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2.lower_bound = self.biomass_flux


@unittest.skipIf(TRAVIS, 'CPLEX not available on Travis.')
class TestFluxVariabilityAnalysisCPLEX(AbstractTestFluxVariabilityAnalysis, unittest.TestCase):
    def setUp(self):
        self.model = CORE_MODEL.copy()
        self.model.solver = 'cplex'
        self.model.solver.problem.parameters.lpmethod = self.model.solver.problem.parameters.lpmethod.values.dual
        self.biomass_flux = 0.873921
        self.model.reactions.Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2.lower_bound = self.biomass_flux


class AbstractTestPhenotypicPhasePlane(object):
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


class TestPhenotypicPhasePlaneGLPK(AbstractTestPhenotypicPhasePlane, unittest.TestCase):
    def setUp(self):
        self.model = CORE_MODEL.copy()
        self.model.solver = 'glpk'


@unittest.skipIf(TRAVIS, 'CPLEX not available on Travis.')
class TestPhenotypicPhasePlaneCPLEX(AbstractTestPhenotypicPhasePlane, unittest.TestCase):
    def setUp(self):
        self.model = CORE_MODEL.copy()
        self.model.solver = 'cplex'


class AbstractTestSimulationMethods(object):
    def setUp(self):
        self.model = CORE_MODEL

    def test_fba(self):
        solution = fba(self.model)
        self.assertAlmostEqual(solution.objective_value, 0.873921, delta=0.000001)

    def test_pfba(self):
        fba_solution = fba(self.model)
        fba_flux_sum = sum((abs(val) for val in list(fba_solution.fluxes.values())))
        pfba_solution = pfba(self.model)
        pfba_flux_sum = sum((abs(val) for val in list(pfba_solution.fluxes.values())))
        # looks like GLPK finds a parsimonious solution without the flux minimization objective
        self.assertTrue((pfba_flux_sum - fba_flux_sum) < 1e-6,
                        msg="FBA sum is suppose to be lower than PFBA (was %f)" % (pfba_flux_sum - fba_flux_sum))

    def test_pfba_iJO(self):
        fba_solution = fba(iJO_MODEL)
        fba_flux_sum = sum((abs(val) for val in list(fba_solution.fluxes.values())))
        pfba_solution = pfba(iJO_MODEL)
        pfba_flux_sum = sum((abs(val) for val in list(pfba_solution.fluxes.values())))
        print(pfba_flux_sum)
        self.assertTrue((pfba_flux_sum - fba_flux_sum) < 1e-6,
                        msg="FBA sum is suppose to be lower than PFBA (was %f)" % (pfba_flux_sum - fba_flux_sum))

    @unittest.skip('quadratic moma not implemented yet.')
    def test_moma(self):
        pass

    def test_lmoma(self):
        pfba_solution = pfba(self.model)
        solution = lmoma(self.model, reference=pfba_solution)
        distance = sum([abs(solution[v] - pfba_solution[v]) for v in list(pfba_solution.keys())])
        self.assertAlmostEqual(0, distance,
                               delta=1e-6,
                               msg="lmoma distance without knockouts must be 0 (was %f)" % distance)

    def test_room(self):
        pfba_solution = pfba(self.model)
        solution = room(self.model, reference=pfba_solution)
        self.assertAlmostEqual(0, solution.objective_value,
                               delta=1e-6,
                               msg="room objective without knockouts must be 0 (was %f)" % solution.objective_value)

    def test_room_shlomi_2005(self):
        reference = {"b1": -10, "v1": 10, "v2": 5, "v3": 0, "v4": 0, "v5": 0, "v6": 5, "b2": 5, "b3": 5}
        TOY_MODEL_PAPIN_2004.solver = self.model.solver.interface
        with TimeMachine() as tm:
            TOY_MODEL_PAPIN_2004.reactions.v6.knock_out(tm)
            result = room(TOY_MODEL_PAPIN_2004, reference=reference, delta=0, epsilon=0)

        self.assertEquals(
            result.fluxes,
            {'b1': 10.0, 'b2': 5.0, 'b3': 5.0, 'v1': 5.0, 'v2': 5.0, 'v3': 0.0, 'v4': 5.0, 'v5': 5.0, 'v6': 0.0})


class TestSimulationMethodsGLPK(AbstractTestSimulationMethods, unittest.TestCase):
    def setUp(self):
        self.model = CORE_MODEL
        self.model.solver = 'glpk'


@unittest.skipIf(TRAVIS, 'CPLEX not available on Travis.')
class TestSimulationMethodsCPLEX(AbstractTestSimulationMethods, unittest.TestCase):
    def setUp(self):
        self.model = CORE_MODEL
        self.model.solver = 'cplex'


if __name__ == '__main__':
    import nose

    nose.runmodule()
