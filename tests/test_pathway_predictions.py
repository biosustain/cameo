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
from cameo.strain_design.pathway_prediction import PathwayPredictor

TESTDIR = os.path.dirname(__file__)
TESTMODEL = load_model(os.path.join(TESTDIR, 'data/EcoliCore.xml'))
UNIVERSALMODEL = load_model(os.path.join(TESTDIR, 'data/iJO1366.xml'))

TRAVIS = os.getenv('TRAVIS', False)

PATHWAYPREDICTOR = PathwayPredictor(TESTMODEL, universal_model=UNIVERSALMODEL)


class TestPathwayPredictor(unittest.TestCase):
    def setUp(self):
        self.pathway_predictor = PATHWAYPREDICTOR

    def test_setting_incorrect_universal_model_raises(self):
        with self.assertRaisesRegexp(ValueError, 'Provided universal_model.*'):
            PathwayPredictor(TESTMODEL, universal_model='Mickey_Mouse')

    # def test_predict_native_compound_returns_shorter_alternatives(self):
    #     result = self.pathway_predictor.run(product='Phosphoenolpyruvate', max_predictions=1)
    #     self.assertTrue(len(result.pathways) == 1)
    #     self.assertTrue(len(result.pathways[0].pathway) == 3)
    #     self.assertTrue(len(result.pathways[0].adapters) == 0)

    def test_predict_non_native_compound(self):
        result = self.pathway_predictor.run(product='L-Serine', max_predictions=1)
        self.assertTrue(len(result.pathways) == 1)
        self.assertTrue(len(result.pathways[0].reactions) == 3)
        self.assertTrue(len(result.pathways[0].adapters) == 0)
