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

import os
import unittest
from cameo.bioreactor import Organism, BioReactor
from cameo.flux_analysis.dfba import combinatorial_dfba, dfba

from cameo.flux_analysis.simulation import fba, pfba, lmoma
from cameo.parallel import SequentialView, MultiprocessingView
from cameo.io import load_model
from cameo.flux_analysis.analysis import flux_variability_analysis, phenotypic_phase_plane

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

CORE_MODEL = load_model(os.path.join(TESTDIR, 'data/EcoliCore.xml'), sanitize=False)
CORE_MODEL_SANE_NAMES = load_model(os.path.join(TESTDIR, 'data/ecoli_core_model.xml'))
iJO_MODEL = load_model(os.path.join(TESTDIR, 'data/iJO1366.xml'), sanitize=False)

iJO_MODEL_COBRAPY = load_model(os.path.join(TESTDIR, 'data/iJO1366.xml'), solver_interface=None, sanitize=False)


class AbstractTestFluxVariabilityAnalysis(object):
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

    @unittest.skip('Removing cycles is still not robust.')
    def test_flux_variability_sequential_remove_cycles(self):
        fva_solution = flux_variability_analysis(self.model, remove_cycles=True, view=SequentialView())
        print fva_solution.min()
        print fva_solution.max()
        for key in fva_solution.index:
            if abs(REFERENCE_FVA_SOLUTION_ECOLI_CORE['lower_bound'][key]) < 999993:
                self.assertAlmostEqual(fva_solution['lower_bound'][key],
                                       REFERENCE_FVA_SOLUTION_ECOLI_CORE['lower_bound'][key], delta=0.000001)
            if abs(REFERENCE_FVA_SOLUTION_ECOLI_CORE['upper_bound'][key]) < 999993:
                print key
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
        self.model = CORE_MODEL
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
        self.assertAlmostEqual(solution.f, 0.873921, delta=0.000001)

    def test_pfba(self):
        fba_solution = fba(self.model)
        fba_flux_sum = sum((abs(val) for val in fba_solution.x_dict.values()))
        pfba_solution = pfba(self.model)
        pfba_flux_sum = sum((abs(val) for val in pfba_solution.x_dict.values()))
        self.assertTrue((pfba_flux_sum - fba_flux_sum) < 1e-6)  # looks like GLPK finds a parsimonious solution without the flux minimization objective

    def test_pfba_iJO(self):
        fba_solution = fba(iJO_MODEL)
        fba_flux_sum = sum((abs(val) for val in fba_solution.x_dict.values()))
        pfba_solution = pfba(iJO_MODEL)
        pfba_flux_sum = sum((abs(val) for val in pfba_solution.x_dict.values()))
        print pfba_flux_sum
        self.assertTrue(pfba_flux_sum < fba_flux_sum)

    @unittest.skip('quadratic moma not implemented yet.')
    def test_moma(self):
        pass

    def test_lmoma(self):
        pfba_solution = pfba(self.model)
        ref = pfba_solution.x_dict
        lmoma_solution = lmoma(self.model, reference=ref)
        res = lmoma_solution.x_dict
        distance = sum([abs(res[v] - ref[v]) for v in res.keys()])
        print distance
        self.assertAlmostEqual(0, distance, delta=0.000001, msg="moma distance without knockouts must be 0")


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


class DynamicFBATestCase(unittest.TestCase):
    def setUp(self):
        assert "EX_o2_e" in [r.id for r in CORE_MODEL_SANE_NAMES.reactions]
        self.ec = Ecoli(CORE_MODEL_SANE_NAMES, 'Ecoli')
        self.ec_glc = GlucoseUser(CORE_MODEL_SANE_NAMES, 'Ecoli_glucose_user')
        self.ec_ac = AcetateUser(CORE_MODEL_SANE_NAMES, 'Ecoli_acetate_user')

    def test_dynamic_fba_1_organism(self):
        br = BioReactor([self.ec], ['EX_glc_e', 'EX_ac_e', 'EX_o2_e'])
        init = [1, 0.01, 10, 0, 100]
        br.initial_conditions = init

        t0 = 0
        tf = 20
        dt = 1

        result = dfba(br, t0, tf, dt)
        self.assertEqual(result.keys(),
                         ['time', 'volume', 'Ecoli', 'EX_glc_e', 'EX_ac_e', 'EX_o2_e'])

    def test_dynamic_2_organisms(self):
        br = BioReactor()
        br.organisms = [self.ec_glc, self.ec_ac]
        br.metabolites = ['EX_glc_e', 'EX_ac_e']
        y0 = [1, 0.01, 0.01, 10, 0]
        t0 = 0
        tf = 20
        dt = 1

        result = dfba(br, t0, tf, dt, y0)
        self.assertEqual(result.keys(),
                         ['time', 'volume', 'Ecoli_glucose_user', 'Ecoli_acetate_user', 'EX_glc_e', 'EX_ac_e'])

    def test_dynamic_fba_combinations(self):
        br1 = BioReactor(metabolites=['EX_glc_e', 'EX_ac_e'], id='Br1')
        br2 = BioReactor(metabolites=['EX_glc_e', 'EX_ac_e', 'EX_o2_e'], id='Br2')

        br1.initial_conditions = [1, 0.01, 10, 0]
        br2.initial_conditions = [1, 0.01, 10, 0, 10]

        t0 = 0
        tf = 20
        dt = 1

        result = combinatorial_dfba([self.ec_glc, self.ec_ac], [br1, br2], t0, tf, dt)
        self.assertEqual(result.keys(), [('Ecoli_glucose_user', 'Br1'), ('Ecoli_glucose_user', 'Br2'),
                                         ('Ecoli_acetate_user', 'Br1'), ('Ecoli_acetate_user', 'Br2')])


class Ecoli(Organism):
    """
    Fixture class for testing dynamic simulations
    """

    def update(self, volume, growth_rate, substrates):
        env = self.environment

        index = env.metabolites.index('EX_glc_e')
        vlb_glc = float(-10 * substrates[index] / (substrates[index] + 1))
        self.constraints['EX_glc_e'] = (vlb_glc, 0)

        # calculating and updating the acetate uptake constraint
        index = env.metabolites.index('EX_ac_e')
        vlb_ac = float(-10 * substrates[index] / (substrates[index] + 1))
        self.constraints['EX_ac_e'] = (vlb_ac, None)
        self.constraints['EX_o2_e'] = (-10, None)

        index = env.metabolites.index('EX_o2_e')
        vlb_o2 = float(-15 * substrates[index] / (substrates[index] + 0.001))
        self.constraints['EX_o2_e'] = (vlb_o2, None)


class GlucoseUser(Organism):
    """
    Fixture class for testing dynamic simulations
    """

    def update(self, volume, growth_rate, substrates):
        env = self.environment

        index = env.metabolites.index('EX_glc_e')
        vlb_glc = float(-10 * substrates[index] / (substrates[index] + 1))
        self.constraints['EX_glc_e'] = (vlb_glc, 0)

        self.constraints['EX_ac_e'] = (0, None)
        self.constraints['EX_o2_e'] = (-10, None)


class AcetateUser(Organism):
    """
    Fixture class for testing dynamic simulations
    """

    def update(self, volume, growth_rate, substrates):
        env = self.environment

        self.constraints['EX_glc_e'] = (0, 0)

        index = env.metabolites.index('EX_ac_e')
        vlb_ac = float(-10 * substrates[index] / (substrates[index] + 1))
        self.constraints['EX_ac_e'] = (vlb_ac, None)
        self.constraints['EX_o2_e'] = (-10, None)

if __name__ == '__main__':
    import nose

    nose.runmodule()
