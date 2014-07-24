# Copyright (c) 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.
# See LICENSE for details.

import unittest

import os

from cameo.flux_analysis.simulation import fba, pfba
from cameo.parallel import SequentialView, MultiprocessingView
from cameo.io import load_model
from cameo.flux_analysis.analysis import flux_variability_analysis, phenotypic_phase_plane, _cycle_free_fva

import pandas
from pandas.util.testing import assert_frame_equal


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


class TestFluxVariabilityAnalysis(unittest.TestCase):
    def setUp(self):
        self.model = CORE_MODEL
        self.biomass_flux = 0.873921
        self.model.reactions.Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2.lower_bound = self.biomass_flux

    def test_cycle_free_fva(self):
        fva_solution = _cycle_free_fva(self.model)
        assert_dataframes_equal(fva_solution, REFERENCE_FVA_SOLUTION_ECOLI_CORE)

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

    def test_one_variable(self):
        ppp = phenotypic_phase_plane(self.model, ['EX_o2_LPAREN_e_RPAREN_'])
        assert_dataframes_equal(ppp, REFERENCE_PPP_o2_EcoliCore)
        ppp = phenotypic_phase_plane(self.model, 'EX_o2_LPAREN_e_RPAREN_')
        assert_dataframes_equal(ppp, REFERENCE_PPP_o2_EcoliCore)

    def test_two_variables(self):
        ppp2d = phenotypic_phase_plane(self.model, ['EX_o2_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_'])
        assert_dataframes_equal(ppp2d, REFERENCE_PPP_o2_glc_EcoliCore)


class TestSimulationMethods(unittest.TestCase):
    def setUp(self):
        self.model = CORE_MODEL.copy()

    def test_fba(self):
        solution = fba(self.model)
        self.assertAlmostEqual(0.873921, solution.f, delta=0.000001)

    def test_pfba(self):
        # fba_solution = fba(self.model)
        # fba_flux_sum = sum((abs(val) for val in fba_solution.x_dict.values()))
        # solution = pfba(self.model)
        # pfba_flux_sum = sum((abs(val) for val in solution['fluxes'].values()))
        # self.assertTrue(pfba_flux_sum < fba_flux_sum)
        pass

    def test_moma(self):
        pass

    def test_lmoma(self):
        pass


if __name__ == '__main__':
    import nose

    nose.runmodule()
