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
REFERENCE_FVA_SOLUTION_ECOLI_CORE = pandas.DataFrame.from_dict(
    {'G6PDH2r': {'lower_bound': 4.959736180498324, 'upper_bound': 4.960214031421614},
     'AKGDH': {'lower_bound': 5.063908797401745, 'upper_bound': 5.064466288000558},
     'EX_succ_LPAREN_e_RPAREN_': {'lower_bound': 0.0, 'upper_bound': 9.712413884699345e-06},
     'GLNS': {'lower_bound': 0.22346159969999999, 'upper_bound': 0.22356115230783047},
     'ADK1': {'lower_bound': 0.0, 'upper_bound': 9.955227569169445e-05},
     'Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2': {'lower_bound': 0.873921,
                                                                'upper_bound': 0.8739215069685606},
     'PYRt2r': {'lower_bound': -1.474847971749682e-05, 'upper_bound': -0.0},
     'ATPM': {'lower_bound': 8.39, 'upper_bound': 8.39},
     'SUCCt2_2': {'lower_bound': 0.0, 'upper_bound': 0.000132736367599288},
     'PIt2r': {'lower_bound': 3.214893182699999, 'upper_bound': 3.2148950476852747},
     'EX_glu_L_LPAREN_e_RPAREN_': {'lower_bound': 0.0, 'upper_bound': 7.3742397468379295e-06},
     'EX_lac_D_LPAREN_e_RPAREN_': {'lower_bound': 0.0, 'upper_bound': 1.2444034462433251e-05},
     'GAPD': {'lower_bound': 16.02345112836896, 'upper_bound': 16.023610412200544},
     'ACALDt': {'lower_bound': -1.4748485284599155e-05, 'upper_bound': -0.0},
     'TPI': {'lower_bound': 7.4773061039689654, 'upper_bound': 7.477465387800544},
     'CO2t': {'lower_bound': -22.809854884864762, 'upper_bound': -22.80978851636553},
     'SUCOAS': {'lower_bound': -5.064466287964024, 'upper_bound': -5.063908791674233},
     'ACONTb': {'lower_bound': 6.007180371259459, 'upper_bound': 6.007339654900543},
     'TKT2': {'lower_bound': 1.1814154455994548, 'upper_bound': 1.1815747293876484},
     'TKT1': {'lower_bound': 1.4969009265994102, 'upper_bound': 1.4970602102405381},
     'O2t': {'lower_bound': 21.799470570869744, 'upper_bound': 21.79951481637545},
     'EX_etoh_LPAREN_e_RPAREN_': {'lower_bound': 0.0, 'upper_bound': 1.2845454928189781e-05},
     'PYK': {'lower_bound': 1.7580783563546856, 'upper_bound': 1.7585031154729684},
     'EX_nh4_LPAREN_e_RPAREN_': {'lower_bound': -4.765323803061619,
                                 'upper_bound': -4.765316428849474},
     'FORt2': {'lower_bound': 0.0, 'upper_bound': 0.0003982088343690293},
     'SUCCt3': {'lower_bound': 0.0, 'upper_bound': 0.000132736367599288},
     'GLUN': {'lower_bound': 0.0, 'upper_bound': 9.95522393998638e-05},
     'FUMt2_2': {'lower_bound': 0.0, 'upper_bound': -0.0},
     'FUM': {'lower_bound': 5.064307003747672, 'upper_bound': 5.064466288000558},
     'ACALD': {'lower_bound': -1.4748540706932545e-05, 'upper_bound': -0.0},
     'MALS': {'lower_bound': 0.0, 'upper_bound': 0.00039821327663957964},
     'CYTBD': {'lower_bound': 43.5989411418559, 'upper_bound': 43.5990296327509},
     'ENO': {'lower_bound': 14.716065311744945, 'upper_bound': 14.716224596200558},
     'EX_h_LPAREN_e_RPAREN_': {'lower_bound': 17.53085525990589,
                               'upper_bound': 17.530921628128894},
     'FORti': {'lower_bound': 0.0, 'upper_bound': 0.00039820883439745103},
     'EX_glc_LPAREN_e_RPAREN_': {'lower_bound': -10.0, 'upper_bound': -9.999994469318016},
     'H2Ot': {'lower_bound': -29.17584501698604, 'upper_bound': -29.175778648663876},
     'LDH_D': {'lower_bound': -1.2444034462433251e-05, 'upper_bound': -0.0},
     'EX_pyr_LPAREN_e_RPAREN_': {'lower_bound': 0.0, 'upper_bound': 1.474847971749682e-05},
     'FRD7': {'lower_bound': 0.0, 'upper_bound': 999993.9356929963},
     'ME1': {'lower_bound': 0.0, 'upper_bound': 9.955227566883915e-05},
     'ME2': {'lower_bound': 0.0, 'upper_bound': 0.00013273636762990902},
     'FRUpts2': {'lower_bound': 0.0, 'upper_bound': -0.0},
     'MDH': {'lower_bound': 5.064280456879112, 'upper_bound': 5.0647052159765735},
     'GLUDy': {'lower_bound': -4.5418622033397495, 'upper_bound': -4.541755276848562},
     'ETOHt2r': {'lower_bound': -1.2845454924637068e-05, 'upper_bound': -0.0},
     'MALt2_2': {'lower_bound': 0.0, 'upper_bound': -0.0},
     'PTAr': {'lower_bound': 0.0, 'upper_bound': 2.2122727937556874e-05},
     'EX_fum_LPAREN_e_RPAREN_': {'lower_bound': 0.0, 'upper_bound': -0.0},
     'TALA': {'lower_bound': 1.4969009265994413, 'upper_bound': 1.4970602102405381},
     'THD2': {'lower_bound': 0.0, 'upper_bound': 0.00039820910280496946},
     'ICDHyr': {'lower_bound': 6.006782158574233, 'upper_bound': 6.007339654900558},
     'ALCD2x': {'lower_bound': -1.2845454928189781e-05, 'upper_bound': -0.0},
     'PDH': {'lower_bound': 9.282461495041089, 'upper_bound': 9.28285970667297},
     'CS': {'lower_bound': 6.007180371241089, 'upper_bound': 6.007339654900543},
     'GLUt2r': {'lower_bound': -7.3742397468379295e-06, 'upper_bound': -0.0},
     'GLNabc': {'lower_bound': 0.0, 'upper_bound': -0.0},
     'EX_for_LPAREN_e_RPAREN_': {'lower_bound': 0.0, 'upper_bound': 6.636836750217867e-05},
     'FBP': {'lower_bound': 0.0, 'upper_bound': 9.95522756710443e-05},
     'ATPS4r': {'lower_bound': 45.513801109744236, 'upper_bound': 45.51411967708108},
     'AKGt2r': {'lower_bound': -8.296018624065482e-06, 'upper_bound': -0.0},
     'FBA': {'lower_bound': 7.47730610356237, 'upper_bound': 7.477465387800558},
     'EX_h2o_LPAREN_e_RPAREN_': {'lower_bound': 29.17577864858322,
                                 'upper_bound': 29.175845016850896},
     'EX_fru_LPAREN_e_RPAREN_': {'lower_bound': 0.0, 'upper_bound': -0.0},
     'NH4t': {'lower_bound': 4.7653164288, 'upper_bound': 4.765323803039749},
     'D_LACt2': {'lower_bound': -1.2444034462433251e-05, 'upper_bound': -0.0},
     'EX_mal_L_LPAREN_e_RPAREN_': {'lower_bound': 0.0, 'upper_bound': -0.0},
     'GLCpts': {'lower_bound': 9.999994469318016, 'upper_bound': 10.0},
     'GND': {'lower_bound': 4.959736180498324, 'upper_bound': 4.960214031993163},
     'ACKr': {'lower_bound': -2.2122727937556874e-05, 'upper_bound': -0.0},
     'NADH16': {'lower_bound': 38.53456334479944, 'upper_bound': 38.53472262844052},
     'SUCDi': {'lower_bound': 5.064307004379993, 'upper_bound': 999999.0},
     'PPC': {'lower_bound': 2.5039098042219976, 'upper_bound': 2.5044407542688543},
     'GLUSy': {'lower_bound': 0.0, 'upper_bound': 9.955224896884829e-05},
     'PGL': {'lower_bound': 4.959736180498364, 'upper_bound': 4.9602140303491264},
     'PGM': {'lower_bound': -14.716224596200558, 'upper_bound': -14.716065312559458},
     'EX_pi_LPAREN_e_RPAREN_': {'lower_bound': -3.2148950476852747,
                                'upper_bound': -3.2148931827396154},
     'PGK': {'lower_bound': -16.02361041220059, 'upper_bound': -16.023451128559458},
     'PGI': {'lower_bound': 4.860632163523265, 'upper_bound': 4.8611100145016355},
     'PPS': {'lower_bound': 0.0, 'upper_bound': 9.955227569169445e-05},
     'RPE': {'lower_bound': 2.678316372198883, 'upper_bound': 2.6786349394810767},
     'RPI': {'lower_bound': -2.2815790921310337, 'upper_bound': -2.2814198082994546},
     'NADTRHD': {'lower_bound': 0.0, 'upper_bound': 0.0003982091027996404},
     'EX_gln_L_LPAREN_e_RPAREN_': {'lower_bound': 0.0, 'upper_bound': -0.0},
     'ICL': {'lower_bound': 0.0, 'upper_bound': 0.00039821327663957964},
     'EX_acald_LPAREN_e_RPAREN_': {'lower_bound': 0.0, 'upper_bound': 1.4748485284599155e-05},
     'EX_o2_LPAREN_e_RPAREN_': {'lower_bound': -21.799514816754773,
                                'upper_bound': -21.79947057092795},
     'ACONTa': {'lower_bound': 6.007180370546308, 'upper_bound': 6.007339654900558},
     'PPCK': {'lower_bound': 0.0, 'upper_bound': 9.955227568836378e-05},
     'EX_ac_LPAREN_e_RPAREN_': {'lower_bound': 0.0, 'upper_bound': 2.2122727937556874e-05},
     'PFK': {'lower_bound': 7.477306103344945, 'upper_bound': 7.4774852982556705},
     'PFL': {'lower_bound': 0.0, 'upper_bound': 6.636836750217867e-05},
     'EX_co2_LPAREN_e_RPAREN_': {'lower_bound': 22.80978851628419,
                                 'upper_bound': 22.809854884744922},
     'EX_akg_LPAREN_e_RPAREN_': {'lower_bound': 0.0, 'upper_bound': 8.296018624065482e-06},
     'ACt2r': {'lower_bound': -2.2122727921569663e-05, 'upper_bound': -0.0}}, orient='index')

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
        self.model = CORE_MODEL


class TestSimulationMethods(unittest.TestCase):
    def setUp(self):
        self.model = CORE_MODEL

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
