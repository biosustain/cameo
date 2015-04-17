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
from cameo.util import RandomGenerator as Random
from cameo.strain_design.deterministic import fseof
from cameo.parallel import SequentialView, MultiprocessingView
from six.moves import range

TESTDIR = os.path.dirname(__file__)
iJO_MODEL = load_model(os.path.join(TESTDIR, 'data/iJO1366.xml'), sanitize=False)


class TestFSEOF(unittest.TestCase):
    def setUp(self):
        self.model = iJO_MODEL.copy()
        self.model.solver = 'glpk'

    def test_fseof(self):
        fseof_result = fseof(self.model, enforced_reaction="EX_succ_LPAREN_e_RPAREN_")
        self.assertIsInstance(fseof_result, list)


if __name__ == "__main__":
    import nose
    nose.runmodule()