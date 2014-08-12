# Copyright (c) 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.
# See LICENSE for details.

import unittest

import os

from cameo.flux_analysis.simulation import fba, pfba, lmoma
from cameo.parallel import SequentialView, MultiprocessingView
from cameo.io import load_model
from cameo.flux_analysis.analysis import flux_variability_analysis, phenotypic_phase_plane, _cycle_free_fva

import pandas
from pandas.util.testing import assert_frame_equal
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

CORE_MODEL = load_model(os.path.join(TESTDIR, 'data/EcoliCore.xml'))
iJO_MODEL = load_model(os.path.join(TESTDIR, 'data/iJO1366.xml'))

iJO_MODEL_COBRAPY = load_model(os.path.join(TESTDIR, 'data/iJO1366.xml'), solver_interface=None)


class TestFluxVariabilityAnalysis(unittest.TestCase):
    def setUp(self):
        self.model = CORE_MODEL
        self.biomass_flux = 0.873921
        self.model.reactions.Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2.lower_bound = self.biomass_flux

    def test_flux_variability_sequential(self):
        fva_solution = flux_variability_analysis(self.model, remove_cycles=False, view=SequentialView())
        assert_dataframes_equal(fva_solution, REFERENCE_FVA_SOLUTION_ECOLI_CORE)
        print REFERENCE_FVA_SOLUTION_ECOLI_CORE
        for key in fva_solution.index:
            self.assertAlmostEqual(fva_solution['lower_bound'][key],
                                   REFERENCE_FVA_SOLUTION_ECOLI_CORE['lower_bound'][key], delta=0.000001)
            self.assertAlmostEqual(fva_solution['upper_bound'][key],
                                   REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound'][key], delta=0.000001)
        assert_dataframes_equal(fva_solution, REFERENCE_FVA_SOLUTION_ECOLI_CORE)

    def test_flux_variability_sequential_remove_cycles(self):
        fva_solution = flux_variability_analysis(self.model, remove_cycles=True, view=SequentialView())
        for key in fva_solution.index:
            if abs(REFERENCE_FVA_SOLUTION_ECOLI_CORE['lower_bound'][key]) < 999993:
                self.assertAlmostEqual(fva_solution['lower_bound'][key],
                                       REFERENCE_FVA_SOLUTION_ECOLI_CORE['lower_bound'][key], delta=0.000001)
            if abs(REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound'][key]) < 999993:
                self.assertAlmostEqual(fva_solution['upper_bound'][key],
                                       REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound'][key], delta=0.000001)

    @unittest.skipIf(TRAVIS, 'Running multiprocess in Travis breaks')
    def test_flux_variability_parallel(self):
        fva_solution = flux_variability_analysis(self.model, remove_cycles=False, view=MultiprocessingView())
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


class TestPhenotypicPhasePlane(unittest.TestCase):
    def setUp(self):
        self.model = CORE_MODEL.copy()

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


class TestSimulationMethods(unittest.TestCase):
    def setUp(self):
        self.model = CORE_MODEL.copy()

    def test_fba(self):
        solution = fba(self.model)
        print self.model.objective
        print self.model.solver
        self.assertAlmostEqual(solution.f, 0.873921, delta=0.000001)

    def test_pfba(self):
        fba_solution = fba(iJO_MODEL)
        fba_flux_sum = sum((abs(val) for val in fba_solution.x_dict.values()))
        pfba_solution = pfba(iJO_MODEL)
        pfba_flux_sum = sum((abs(val) for val in pfba_solution.x_dict.values()))
        print pfba_flux_sum
        self.assertTrue(pfba_flux_sum < fba_flux_sum)

    def test_moma(self):
        pass

    def test_lmoma(self):
        pfba_solution = pfba(self.model)
        ref = pfba_solution.x_dict
        lmoma_solution = lmoma(self.model, reference=ref)
        res = lmoma_solution.x_dict
        distance = sum([abs(res[v] - ref[v]) for v in res.keys()])
        self.assertAlmostEqual(0, distance, delta=0.000001, msg="moma distance without knockouts must be 0")

if __name__ == '__main__':
    import nose

    nose.runmodule()
