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

# import types
#
from cobra.test import create_test_model
from cameo.core.solver_based_model import to_solver_based_model
import unittest
import six
from cobra.flux_analysis import *


class TestCobrapyCompatibility(unittest.TestCase):
    def setUp(self):
        # Make Model pickable and then load a solver based version of test_pickle
        self.model = to_solver_based_model(create_test_model("textbook"))

    def test_cobra_phenotypic_phase_plane(self):
        data = calculate_phenotype_phase_plane(self.model, "EX_glc__D_e", "EX_o2_e",
                                               reaction1_npoints=20, reaction2_npoints=20)
        self.assertEquals(data.growth_rates.shape, (20, 20))
        self.assertTrue(abs(data.growth_rates.max()) - 1.20898 < 0.001)
        self.assertTrue(abs(data.growth_rates[0, :].max()) < 0.0001)

    def test_single_gene_deletion_fba(self):
        growth_dict = {"b0008": 0.87, "b0114": 0.80, "b0116": 0.78,
                       "b2276": 0.21, "b1779": 0.00}
        rates, statuses = single_gene_deletion(self.model,
                                               gene_list=growth_dict.keys(),
                                               method="fba")
        for gene, expected_value in six.iteritems(growth_dict):
            self.assertTrue(statuses[gene] == 'optimal')
            self.assertTrue(abs(rates[gene] - expected_value) < 0.01)
