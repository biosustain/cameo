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

import six

import os
import unittest

from cameo import load_model
from cameo.strain_design.deterministic.flux_variability_based import Fseof, FseofResult, DifferentialFVA

from pandas import DataFrame, pandas
from pandas.util.testing import assert_frame_equal

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
        fseof = Fseof(self.model, enforced_reaction="EX_succ_lp_e_rp_")
        fseof_result = fseof.run()
        self.assertIsInstance(fseof_result, FseofResult)
        self.assertIs(objective, self.model.objective)

    def test_fseof_result(self):
        fseof = Fseof(self.model, self.model.reactions.EX_ac_lp_e_rp_)
        fseof_result = fseof.run()
        self.assertIsInstance(fseof_result.data_frame, DataFrame)
        self.assertIs(fseof_result.enforced_reaction, self.model.reactions.EX_ac_lp_e_rp_)
        self.assertIs(fseof_result.model, self.model)
        self.assertEqual(list(fseof_result), list(fseof_result.reactions))

if six.PY2:  # Make these test cases work with PY3 as well
    class TestDifferentialFVA(unittest.TestCase):
        def setUp(self):
            self.model = ECOLICORE

        def test_minimal_input(self):
            result = DifferentialFVA(self.model, self.model.reactions.EX_succ_lp_e_rp_, points=5).run()
            # result.data_frame.iloc[0].to_pickle(os.path.join(TESTDIR, 'data/REFERENCE_DiffFVA1.pickle'))
            pandas.util.testing.assert_frame_equal(result.data_frame.iloc[0], pandas.read_pickle(os.path.join(TESTDIR, 'data/REFERENCE_DiffFVA1.pickle')))

        def test_with_reference_model(self):
            reference_model = self.model.copy()
            biomass_rxn = reference_model.reactions.Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2
            biomass_rxn.lower_bound = 0.3
            target = reference_model.reactions.EX_succ_lp_e_rp_
            target.lower_bound = 2
            result = DifferentialFVA(self.model, target, reference_model=reference_model, points=5).run()
            # result.data_frame.iloc[0].to_pickle(os.path.join(TESTDIR, 'data/REFERENCE_DiffFVA2.pickle'))
            pandas.util.testing.assert_frame_equal(result.data_frame.iloc[0], pandas.read_pickle(os.path.join(TESTDIR, 'data/REFERENCE_DiffFVA2.pickle')))

if __name__ == "__main__":
    import nose

    nose.runmodule()
