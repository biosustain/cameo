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

from cameo import load_model
from cameo.strain_design.deterministic.flux_variability_based import fseof, FseofResult
from six.moves import range
from pandas import DataFrame

TESTDIR = os.path.dirname(__file__)
iJO_MODEL = load_model(os.path.join(TESTDIR, 'data/iJO1366.xml'), sanitize=False)


class TestFSEOF(unittest.TestCase):
    def setUp(self):
        self.model = iJO_MODEL.copy()
        self.model.solver = 'glpk'

    def test_fseof(self):
        objective = self.model.objective
        fseof_result = fseof(self.model, enforced_reaction="EX_succ_LPAREN_e_RPAREN_")
        self.assertIsInstance(fseof_result, FseofResult)
        self.assertIs(objective, self.model.objective)

    def test_fseof_result(self):
        fseof_result = fseof(self.model, self.model.reactions.EX_ac_LPAREN_e_RPAREN_, 0.8, exclude=["PGI"])
        self.assertIsInstance(fseof_result.data_frame, DataFrame)
        self.assertIs(fseof_result.objective, self.model.reactions.EX_ac_LPAREN_e_RPAREN_)
        self.assertIs(fseof_result.model, self.model)
        self.assertEqual(list(fseof_result), list(fseof_result.reactions))


if __name__ == "__main__":
    import nose

    nose.runmodule()
