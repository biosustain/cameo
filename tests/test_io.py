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

import six

import os
import unittest

import cobra
import optlang

from cameo import load_model
from cameo.config import solvers
from cameo.core.solver_based_model import SolverBasedModel

TESTDIR = os.path.dirname(__file__)


class AbstractTestModelLoading(object):
    def test_load_model_pickle_path(self):
        model = load_model(os.path.join(TESTDIR, 'data/iJO1366.pickle'), solver_interface=self.interface)
        self.assertAlmostEqual(model.optimize().f, 0.9823718127269768)

    def test_load_model_pickle_handle(self):
        with open(os.path.join(TESTDIR, 'data/iJO1366.pickle'), 'rb') as handle:
            model = load_model(handle, solver_interface=self.interface)
        self.assertAlmostEqual(model.optimize().f, 0.9823718127269768)

    @unittest.skipIf(six.PY3, 'cobra.io.read_sbml_model broken in py3.')
    def test_load_model_sbml_path(self):
        model = load_model(os.path.join(TESTDIR, 'data/iJO1366.xml'), solver_interface=self.interface)
        self.assertAlmostEqual(model.optimize().f, 0.9823718127269768)

    @unittest.skipIf(six.PY3, 'cobra.io.read_sbml_model broken in py3.')
    def test_load_model_sbml_handle(self):
        with open(os.path.join(TESTDIR, 'data/iJO1366.xml')) as handle:
            model = load_model(handle, solver_interface=self.interface)
        self.assertAlmostEqual(model.optimize().f, 0.9823718127269768)

    @unittest.skipIf(six.PY3, 'cobra.io.read_sbml_model broken in py3.')
    def test_load_model_sbml_path_set_none_interface(self):
        model = load_model(os.path.join(TESTDIR, 'data/EcoliCore.xml'), solver_interface=None)
        self.assertAlmostEqual(model.optimize().f, 0.8739215069684306)
        self.assertTrue(isinstance(model, cobra.core.Model))
        self.assertFalse(hasattr(model, 'solver'))


class TestModelLoadingGLPK(AbstractTestModelLoading, unittest.TestCase):
    def setUp(self):
        self.interface = optlang.glpk_interface


@unittest.skipIf('cplex' not in solvers, "No cplex interface available")
class TestModelLoadingCPLEX(AbstractTestModelLoading, unittest.TestCase):
    def setUp(self):
        self.interface = optlang.glpk_interface


if __name__ == '__main__':
    import nose

    nose.runmodule()
