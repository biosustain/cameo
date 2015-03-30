# Copyright 2015 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
from unittest import TestCase
from cameo import load_model
from cameo.dynamic import batch_dfba
import os

MODEL = load_model(os.path.join(os.path.dirname(__file__), "data/iJO1366.xml"))


class DFBACase(TestCase):
    def assertFloatListAlmostEqual(self, list1, list2, places=None, delta=None):
        self.assertEqual(len(list1), len(list2), "Lists must have the same size")
        for i, v in enumerate(list1):
            self.assertAlmostEqual(v, list2[i], places=places, delta=delta)

    def setUp(self):
        self.times = [0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.]
        self.batch_volumes = [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.]
        self.batch_metabolites = {
            "EX_glc_lp_e_rp_": [10.,
                                8.84881022,
                                7.05239141,
                                4.34479529,
                                0.95071594,
                                0.01358383,
                                0.01358383,
                                0.01358383,
                                0.01358383,
                                0.01358383,
                                0.01358383],
            "EX_ac_lp_e_rp_": [0.,
                               1.34850845,
                               3.46513311,
                               6.68721992,
                               10.31430285,
                               10.66605171,
                               10.66605171,
                               10.66605171,
                               10.66605171,
                               10.66605171,
                               10.66605171]
        }
        self.batch_biomass = {
            "iJO1366": [0.1,
                        0.15920909,
                        0.25237768,
                        0.39578862,
                        0.59210645,
                        0.66802286,
                        0.66802286,
                        0.66802286,
                        0.66802286,
                        0.66802286,
                        0.66802286]
        }

    def test_batch(self):
        def update(volume, growth_rate, metabolites, substrates):
            index = metabolites.index('EX_glc_lp_e_rp_')

            vlb_glc = float(-10 * substrates[index] / (substrates[index] + 1))

            return {'EX_glc_lp_e_rp_': (vlb_glc, 0), 'EX_o2_lp_e_rp_': (-5, 0)}, None

        # organism = Organism(MODEL, constraints={'EX_o2_lp_e_rp_': (-5, 0), 'EX_ac_lp_e_rp_': (19, 20)})
        # organism.update = update.__get__(organism, Organism)
        # reactor = IdealBatch(organisms=[organism],
        #                      metabolites=['EX_glc_lp_e_rp_', 'EX_ac_lp_e_rp_'],
        #                      initial_conditions=[1, 0.1, 10, 0])

        result = batch_dfba(metabolites=['EX_glc_lp_e_rp_', 'EX_ac_lp_e_rp_'], initial_concentrations=[10, 0],
                            models=[MODEL], initial_biomass=[0.1], models_dynamics=[update],
                            initial_volume=1., t0=0, tf=10, dt=1)

        self.assertFloatListAlmostEqual(result["EX_glc_lp_e_rp_"], self.batch_metabolites["EX_glc_lp_e_rp_"], delta=0.0001)
        self.assertFloatListAlmostEqual(result["EX_ac_lp_e_rp_"], self.batch_metabolites["EX_ac_lp_e_rp_"], delta=0.0001)
        self.assertFloatListAlmostEqual(result["iJO1366"], self.batch_biomass["iJO1366"], delta=0.001)