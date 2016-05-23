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
from tempfile import _TemporaryFileWrapper
from unittest import TestCase

import requests
from pandas import DataFrame

from cameo.models.webmodels import index_models_minho, get_sbml_file, NotFoundException


class WebmodelsTestCase(TestCase):
    def test_invalid_host(self):
        self.assertRaises(requests.ConnectionError, index_models_minho, host="http://blabla")
        self.assertRaises(requests.ConnectionError, get_sbml_file, 1, host="http://blabla")

    def test_index(self):
        index = index_models_minho()
        self.assertIsInstance(index, DataFrame)
        self.assertListEqual(list(index.columns),
                             ["id", "name", "doi", "author", "year", "formats", "organism", "taxonomy", "validated"])

    def test_get_sbml(self):
        tmp = get_sbml_file(1)
        self.assertIsInstance(tmp, _TemporaryFileWrapper)
        self.assertRaises(NotFoundException, get_sbml_file, -1)
