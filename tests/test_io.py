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

from __future__ import absolute_import, print_function

import os

import cobra
import pytest

import cameo
from cameo import load_model
from cameo.config import solvers

try:
    import libsbml
except ImportError:
    libsbml = None

TESTDIR = os.path.dirname(__file__)


@pytest.fixture(scope="module", params=list(solvers))
def solver_interface(request):
    return solvers[request.param]


class TestModelLoading(object):
    def test_load_model_pickle_path(self, solver_interface):
        model = load_model(os.path.join(TESTDIR, 'data/iJO1366.pickle'), solver_interface=solver_interface)
        assert abs(model.optimize().objective_value - 0.9823718127269768) < 10e-6

    def test_load_model_pickle_handle(self, solver_interface):
        with open(os.path.join(TESTDIR, 'data/iJO1366.pickle'), 'rb') as handle:
            model = load_model(handle, solver_interface=solver_interface)
        assert abs(model.slim_optimize() - 0.9823718127269768) < 10e-6

    def test_load_model_sbml_path(self, solver_interface):
        model = load_model(os.path.join(TESTDIR, 'data/iJO1366.xml'), solver_interface=solver_interface)
        assert abs(model.slim_optimize() - 0.9823718127269768) < 10e-6

    def test_load_model_sbml_handle(self, solver_interface):
        with open(os.path.join(TESTDIR, 'data/iJO1366.xml')) as handle:
            model = load_model(handle, solver_interface=solver_interface)
        assert abs(model.slim_optimize() - 0.9823718127269768) < 10e-6

    def test_load_model_sbml_path_set_none_interface(self):
        model = load_model(os.path.join(TESTDIR, 'data/EcoliCore.xml'), solver_interface=None)
        assert abs(model.slim_optimize() - 0.8739215069684306) < 10e-6
        assert isinstance(model, cobra.Model)

    def test_import_model_bigg(self):
        model = cameo.models.bigg.e_coli_core
        assert model.id == 'e_coli_core'

    @pytest.mark.skipif(libsbml is None, reason="minho has fbc < 2, requiring missing lisbml")
    def test_import_model_minho(self):
        model = cameo.models.minho
        if model.status != 'indexed':
            pytest.skip('failed to index minho db')
        assert model.__getattr__('Ecoli core Model').id == 'Ecoli_core_model'

    def test_invalid_path(self):
        with pytest.raises(Exception):
            load_model("blablabla_model")
