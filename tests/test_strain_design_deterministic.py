# Copyright 2015 Novo Nordisk Foundation Center for Biosustainability, DTU.
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

from __future__ import absolute_import, print_function

import os
import unittest

import pandas
import six
from pandas import DataFrame
from pandas.util.testing import assert_frame_equal

import cameo
from cameo import load_model
from cameo.strain_design.deterministic.flux_variability_based import FSEOF, FSEOFResult, DifferentialFVA
from cameo.strain_design.deterministic.linear_programming import OptKnock

TRAVIS = os.getenv('TRAVIS', False)
TESTDIR = os.path.dirname(__file__)
ECOLICORE = load_model(os.path.join(TESTDIR, 'data/EcoliCore.xml'))


def assert_dataframes_equal(df, expected):
    try:
        assert_frame_equal(df, expected, check_names=False)
        return True
    except AssertionError:
        return False


class TestFSEOF(unittest.TestCase):
    def setUp(self):
        self.model = ECOLICORE.copy()
        self.model.solver = 'glpk'

    def test_fseof(self):
        objective = self.model.objective
        fseof = FSEOF(self.model)
        fseof_result = fseof.run(target="EX_succ_lp_e_rp_")
        self.assertIsInstance(fseof_result, FSEOFResult)
        self.assertIs(objective, self.model.objective)

    def test_fseof_result(self):
        fseof = FSEOF(self.model)
        fseof_result = fseof.run(target=self.model.reactions.EX_ac_lp_e_rp_)
        self.assertIsInstance(fseof_result.data_frame, DataFrame)
        self.assertIs(fseof_result.target, self.model.reactions.EX_ac_lp_e_rp_)
        self.assertIs(fseof_result.model, self.model)


# if six.PY2:  # Make these test cases work with PY3 as well
class TestDifferentialFVA(unittest.TestCase):
    def setUp(self):
        self.model = ECOLICORE.copy()

    def test_minimal_input(self):
        result = DifferentialFVA(self.model, self.model.reactions.EX_succ_lp_e_rp_, points=5).run()

        # result.data_frame.iloc[0].to_csv(os.path.join(TESTDIR, 'data/REFERENCE_DiffFVA1.csv'))
        ref_df = pandas.read_csv(os.path.join(TESTDIR, 'data/REFERENCE_DiffFVA1.csv'), index_col=0).astype(
            "O").sort_index(axis=1)
        pandas.util.testing.assert_frame_equal(result.data_frame.iloc[0].sort_index(axis=1), ref_df)

    def test_with_reference_model(self):
        reference_model = self.model.copy()
        biomass_rxn = reference_model.reactions.Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2
        biomass_rxn.lower_bound = 0.3
        target = reference_model.reactions.EX_succ_lp_e_rp_
        target.lower_bound = 2
        result = DifferentialFVA(self.model, target, reference_model=reference_model, points=5).run()
        # result.data_frame.iloc[0].to_csv(os.path.join(TESTDIR, 'data/REFERENCE_DiffFVA2.csv'))
        ref_df = pandas.read_csv(os.path.join(TESTDIR, 'data/REFERENCE_DiffFVA2.csv'), index_col=0).astype(
            "O").sort_index(axis=1)
        pandas.util.testing.assert_frame_equal(result.data_frame.iloc[0].sort_index(axis=1), ref_df)


class TestOptKnock(unittest.TestCase):
    def setUp(self):
        self.model = ECOLICORE.copy()
        self.model.reactions.Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2.lower_bound = 0.1
        self.model.solver = "cplex"
        self.optknock = OptKnock(self.model)

    def test_optknock_runs(self):
        result = self.optknock.run(max_knockouts=0, target="EX_ac_lp_e_rp_",
                                   biomass="Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2", max_results=1)
        self.assertEqual(len(result), 1)
        self.assertEqual(len(result.knockouts[0]), 0)
        self.assertEqual(len(list(result)), 1)
        self.assertIsInstance(result.data_frame, DataFrame)

    def test_result_is_correct(self):
        result = self.optknock.run(max_knockouts=1, target="EX_ac_lp_e_rp_",
                                   biomass="Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2", max_results=1)
        production = result.production[0]
        knockouts = result.knockouts[0]
        for knockout in knockouts:
            self.model.reactions.get_by_id(knockout).knock_out()
        fva = cameo.flux_variability_analysis(self.model, fraction_of_optimum=1, remove_cycles=False,
                                              reactions=["EX_ac_lp_e_rp_"])
        self.assertAlmostEqual(fva["upper_bound"][0], production)

    def test_invalid_input(self):
        self.assertRaises(KeyError, self.optknock.run, target="EX_ac_lp_e_rp_")
        self.assertRaises(KeyError, self.optknock.run, biomass="Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2")


if __name__ == "__main__":
    import nose

    nose.runmodule()
