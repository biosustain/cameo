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
import unittest

import cobra
import optlang

import cameo
from cameo import load_model

TESTDIR = os.path.dirname(__file__)
TRAVIS = os.getenv('TRAVIS', False)


class AbstractTestModelLoading(object):
    def test_load_model_pickle_path(self):
        model = load_model(os.path.join(TESTDIR, 'data/iJO1366.pickle'), solver_interface=self.interface)
        self.assertAlmostEqual(model.optimize().f, 0.9823718127269768, places=6)

    def test_load_model_pickle_handle(self):
        with open(os.path.join(TESTDIR, 'data/iJO1366.pickle'), 'rb') as handle:
            model = load_model(handle, solver_interface=self.interface)
        self.assertAlmostEqual(model.optimize().f, 0.9823718127269768, places=6)

    # @unittest.skipIf(six.PY3, 'cobra.io.read_sbml_model broken in py3.')
    def test_load_model_sbml_path(self):
        model = load_model(os.path.join(TESTDIR, 'data/iJO1366.xml'), solver_interface=self.interface)
        self.assertAlmostEqual(model.optimize().f, 0.9823718127269768, places=6)

    # @unittest.skipIf(six.PY3, 'cobra.io.read_sbml_model broken in py3.')
    def test_load_model_sbml_handle(self):
        with open(os.path.join(TESTDIR, 'data/iJO1366.xml')) as handle:
            model = load_model(handle, solver_interface=self.interface)
        self.assertAlmostEqual(model.optimize().f, 0.9823718127269768, places=6)

    # @unittest.skipIf(six.PY3, 'cobra.io.read_sbml_model broken in py3.')
    def test_load_model_sbml_path_set_none_interface(self):
        model = load_model(os.path.join(TESTDIR, 'data/EcoliCore.xml'), solver_interface=None)
        self.assertAlmostEqual(model.optimize().f, 0.8739215069684306, places=6)
        self.assertTrue(isinstance(model, cobra.core.Model))
        self.assertFalse(hasattr(model, 'solver'))

    def test_import_model_bigg(self):
        model = cameo.models.bigg.e_coli_core
        self.assertEqual(model.id, 'e_coli_core')

    def test_import_model_minho(self):
        model = cameo.models.bigg.e_coli_core
        self.assertEqual(model.id, 'e_coli_core')

    def test_invalid_path(self):
        self.assertRaises(Exception, load_model, "blablabla_model")


class TestModelLoadingGLPK(AbstractTestModelLoading, unittest.TestCase):
    def setUp(self):
        self.interface = optlang.glpk_interface


# @unittest.skipIf(six.PY2, 'Build stalling in python 2.7.')
class TestModelLoadingCPLEX(AbstractTestModelLoading, unittest.TestCase):
    def setUp(self):
        self.interface = optlang.cplex_interface


if __name__ == '__main__':
    import nose

    nose.runmodule()
