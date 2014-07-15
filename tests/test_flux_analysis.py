# Copyright (c) 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.
# See LICENSE for details.

import unittest

import os

import cameo
from cameo.flux_analysis.simulation import fba, pfba
from cameo.parallel import SequentialView, MultiprocessingView


cameo.config.default_view = SequentialView()
from cameo.io import load_model
from cameo.flux_analysis.analysis import flux_variability_analysis, _flux_variability_analysis
from cobra.io import read_sbml_model


TESTDIR = os.path.dirname(__file__)
REFERENCE_FVA_SOLUTION_ECOLI_CORE = {'G6PDH2r': {'minimum': 4.959736180498324, 'maximum': 4.960214031421614},
                                     'AKGDH': {'minimum': 5.063908797401745, 'maximum': 5.064466288000558},
                                     'EX_succ_LPAREN_e_RPAREN_': {'minimum': 0.0, 'maximum': 9.712413884699345e-06},
                                     'GLNS': {'minimum': 0.22346159969999999, 'maximum': 0.22356115230783047},
                                     'ADK1': {'minimum': 0.0, 'maximum': 9.955227569169445e-05},
                                     'Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2': {'minimum': 0.873921,
                                                                                                'maximum': 0.8739215069685606},
                                     'PYRt2r': {'minimum': -1.474847971749682e-05, 'maximum': -0.0},
                                     'ATPM': {'minimum': 8.39, 'maximum': 8.39},
                                     'SUCCt2_2': {'minimum': 0.0, 'maximum': 0.000132736367599288},
                                     'PIt2r': {'minimum': 3.214893182699999, 'maximum': 3.2148950476852747},
                                     'EX_glu_L_LPAREN_e_RPAREN_': {'minimum': 0.0, 'maximum': 7.3742397468379295e-06},
                                     'EX_lac_D_LPAREN_e_RPAREN_': {'minimum': 0.0, 'maximum': 1.2444034462433251e-05},
                                     'GAPD': {'minimum': 16.02345112836896, 'maximum': 16.023610412200544},
                                     'ACALDt': {'minimum': -1.4748485284599155e-05, 'maximum': -0.0},
                                     'TPI': {'minimum': 7.4773061039689654, 'maximum': 7.477465387800544},
                                     'CO2t': {'minimum': -22.809854884864762, 'maximum': -22.80978851636553},
                                     'SUCOAS': {'minimum': -5.064466287964024, 'maximum': -5.063908791674233},
                                     'ACONTb': {'minimum': 6.007180371259459, 'maximum': 6.007339654900543},
                                     'TKT2': {'minimum': 1.1814154455994548, 'maximum': 1.1815747293876484},
                                     'TKT1': {'minimum': 1.4969009265994102, 'maximum': 1.4970602102405381},
                                     'O2t': {'minimum': 21.799470570869744, 'maximum': 21.79951481637545},
                                     'EX_etoh_LPAREN_e_RPAREN_': {'minimum': 0.0, 'maximum': 1.2845454928189781e-05},
                                     'PYK': {'minimum': 1.7580783563546856, 'maximum': 1.7585031154729684},
                                     'EX_nh4_LPAREN_e_RPAREN_': {'minimum': -4.765323803061619,
                                                                 'maximum': -4.765316428849474},
                                     'FORt2': {'minimum': 0.0, 'maximum': 0.0003982088343690293},
                                     'SUCCt3': {'minimum': 0.0, 'maximum': 0.000132736367599288},
                                     'GLUN': {'minimum': 0.0, 'maximum': 9.95522393998638e-05},
                                     'FUMt2_2': {'minimum': 0.0, 'maximum': -0.0},
                                     'FUM': {'minimum': 5.064307003747672, 'maximum': 5.064466288000558},
                                     'ACALD': {'minimum': -1.4748540706932545e-05, 'maximum': -0.0},
                                     'MALS': {'minimum': 0.0, 'maximum': 0.00039821327663957964},
                                     'CYTBD': {'minimum': 43.5989411418559, 'maximum': 43.5990296327509},
                                     'ENO': {'minimum': 14.716065311744945, 'maximum': 14.716224596200558},
                                     'EX_h_LPAREN_e_RPAREN_': {'minimum': 17.53085525990589,
                                                               'maximum': 17.530921628128894},
                                     'FORti': {'minimum': 0.0, 'maximum': 0.00039820883439745103},
                                     'EX_glc_LPAREN_e_RPAREN_': {'minimum': -10.0, 'maximum': -9.999994469318016},
                                     'H2Ot': {'minimum': -29.17584501698604, 'maximum': -29.175778648663876},
                                     'LDH_D': {'minimum': -1.2444034462433251e-05, 'maximum': -0.0},
                                     'EX_pyr_LPAREN_e_RPAREN_': {'minimum': 0.0, 'maximum': 1.474847971749682e-05},
                                     'FRD7': {'minimum': 0.0, 'maximum': 999993.9356929963},
                                     'ME1': {'minimum': 0.0, 'maximum': 9.955227566883915e-05},
                                     'ME2': {'minimum': 0.0, 'maximum': 0.00013273636762990902},
                                     'FRUpts2': {'minimum': 0.0, 'maximum': -0.0},
                                     'MDH': {'minimum': 5.064280456879112, 'maximum': 5.0647052159765735},
                                     'GLUDy': {'minimum': -4.5418622033397495, 'maximum': -4.541755276848562},
                                     'ETOHt2r': {'minimum': -1.2845454924637068e-05, 'maximum': -0.0},
                                     'MALt2_2': {'minimum': 0.0, 'maximum': -0.0},
                                     'PTAr': {'minimum': 0.0, 'maximum': 2.2122727937556874e-05},
                                     'EX_fum_LPAREN_e_RPAREN_': {'minimum': 0.0, 'maximum': -0.0},
                                     'TALA': {'minimum': 1.4969009265994413, 'maximum': 1.4970602102405381},
                                     'THD2': {'minimum': 0.0, 'maximum': 0.00039820910280496946},
                                     'ICDHyr': {'minimum': 6.006782158574233, 'maximum': 6.007339654900558},
                                     'ALCD2x': {'minimum': -1.2845454928189781e-05, 'maximum': -0.0},
                                     'PDH': {'minimum': 9.282461495041089, 'maximum': 9.28285970667297},
                                     'CS': {'minimum': 6.007180371241089, 'maximum': 6.007339654900543},
                                     'GLUt2r': {'minimum': -7.3742397468379295e-06, 'maximum': -0.0},
                                     'GLNabc': {'minimum': 0.0, 'maximum': -0.0},
                                     'EX_for_LPAREN_e_RPAREN_': {'minimum': 0.0, 'maximum': 6.636836750217867e-05},
                                     'FBP': {'minimum': 0.0, 'maximum': 9.95522756710443e-05},
                                     'ATPS4r': {'minimum': 45.513801109744236, 'maximum': 45.51411967708108},
                                     'AKGt2r': {'minimum': -8.296018624065482e-06, 'maximum': -0.0},
                                     'FBA': {'minimum': 7.47730610356237, 'maximum': 7.477465387800558},
                                     'EX_h2o_LPAREN_e_RPAREN_': {'minimum': 29.17577864858322,
                                                                 'maximum': 29.175845016850896},
                                     'EX_fru_LPAREN_e_RPAREN_': {'minimum': 0.0, 'maximum': -0.0},
                                     'NH4t': {'minimum': 4.7653164288, 'maximum': 4.765323803039749},
                                     'D_LACt2': {'minimum': -1.2444034462433251e-05, 'maximum': -0.0},
                                     'EX_mal_L_LPAREN_e_RPAREN_': {'minimum': 0.0, 'maximum': -0.0},
                                     'GLCpts': {'minimum': 9.999994469318016, 'maximum': 10.0},
                                     'GND': {'minimum': 4.959736180498324, 'maximum': 4.960214031993163},
                                     'ACKr': {'minimum': -2.2122727937556874e-05, 'maximum': -0.0},
                                     'NADH16': {'minimum': 38.53456334479944, 'maximum': 38.53472262844052},
                                     'SUCDi': {'minimum': 5.064307004379993, 'maximum': 999999.0},
                                     'PPC': {'minimum': 2.5039098042219976, 'maximum': 2.5044407542688543},
                                     'GLUSy': {'minimum': 0.0, 'maximum': 9.955224896884829e-05},
                                     'PGL': {'minimum': 4.959736180498364, 'maximum': 4.9602140303491264},
                                     'PGM': {'minimum': -14.716224596200558, 'maximum': -14.716065312559458},
                                     'EX_pi_LPAREN_e_RPAREN_': {'minimum': -3.2148950476852747,
                                                                'maximum': -3.2148931827396154},
                                     'PGK': {'minimum': -16.02361041220059, 'maximum': -16.023451128559458},
                                     'PGI': {'minimum': 4.860632163523265, 'maximum': 4.8611100145016355},
                                     'PPS': {'minimum': 0.0, 'maximum': 9.955227569169445e-05},
                                     'RPE': {'minimum': 2.678316372198883, 'maximum': 2.6786349394810767},
                                     'RPI': {'minimum': -2.2815790921310337, 'maximum': -2.2814198082994546},
                                     'NADTRHD': {'minimum': 0.0, 'maximum': 0.0003982091027996404},
                                     'EX_gln_L_LPAREN_e_RPAREN_': {'minimum': 0.0, 'maximum': -0.0},
                                     'ICL': {'minimum': 0.0, 'maximum': 0.00039821327663957964},
                                     'EX_acald_LPAREN_e_RPAREN_': {'minimum': 0.0, 'maximum': 1.4748485284599155e-05},
                                     'EX_o2_LPAREN_e_RPAREN_': {'minimum': -21.799514816754773,
                                                                'maximum': -21.79947057092795},
                                     'ACONTa': {'minimum': 6.007180370546308, 'maximum': 6.007339654900558},
                                     'PPCK': {'minimum': 0.0, 'maximum': 9.955227568836378e-05},
                                     'EX_ac_LPAREN_e_RPAREN_': {'minimum': 0.0, 'maximum': 2.2122727937556874e-05},
                                     'PFK': {'minimum': 7.477306103344945, 'maximum': 7.4774852982556705},
                                     'PFL': {'minimum': 0.0, 'maximum': 6.636836750217867e-05},
                                     'EX_co2_LPAREN_e_RPAREN_': {'minimum': 22.80978851628419,
                                                                 'maximum': 22.809854884744922},
                                     'EX_akg_LPAREN_e_RPAREN_': {'minimum': 0.0, 'maximum': 8.296018624065482e-06},
                                     'ACt2r': {'minimum': -2.2122727921569663e-05, 'maximum': -0.0}}

CORE_MODEL = load_model(os.path.join(TESTDIR, 'data/EcoliCore.xml'))


class TestFluxVariabilityAnalysis(unittest.TestCase):
    def setUp(self):
        self.model = CORE_MODEL
        self.biomass_flux = 0.873921
        self.model.reactions.Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2.lower_bound = self.biomass_flux

    def test_flux_variability_sequential(self):
        fva_solution = flux_variability_analysis(self.model, view=SequentialView())
        for key, val in fva_solution.iteritems():
            self.assertAlmostEqual(val['maximum'], REFERENCE_FVA_SOLUTION_ECOLI_CORE[key]['maximum'], delta=0.000001)
            self.assertAlmostEqual(val['minimum'], REFERENCE_FVA_SOLUTION_ECOLI_CORE[key]['minimum'], delta=0.000001)

    def test_flux_variability_parallel(self):
        fva_solution = flux_variability_analysis(self.model, view=MultiprocessingView())
        for key, val in fva_solution.iteritems():
            self.assertAlmostEqual(val['maximum'], REFERENCE_FVA_SOLUTION_ECOLI_CORE[key]['maximum'], delta=0.000001)
            self.assertAlmostEqual(val['minimum'], REFERENCE_FVA_SOLUTION_ECOLI_CORE[key]['minimum'], delta=0.000001)


class TestSimulationMethods(unittest.TestCase):
    def setUp(self):
        self.model = CORE_MODEL
        self.biomass_flux = 0.873921
        self.model.reactions.Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2.lower_bound = self.biomass_flux

    def test_fba(self):
        model = self.model.copy()
        solution = fba(model)
        self.assertAlmostEqual(self.biomass_flux, solution.f, delta=0.000001)

    def test_pfba(self):
        pass
        #model = self.model.copy()
        #solution = pfba(model)
        #self.assertAlmostEqual(self.biomass_flux, solution.f, delta=0.000001)
        #assert net conversion

    def test_moma(self):
        pass

    def test_lmoma(self):
        pass


if __name__ == '__main__':
    import nose
    nose.runmodule()
