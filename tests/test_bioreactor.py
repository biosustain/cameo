"""
Unit testing module for dFBA and MdFBA.
"""
import os

from cameo import load_model
from cameo.dynamic.bioreactor import Organism, BioReactor
from cameo.dynamic.bioreactor.base import Environment


__author__ = 'kaizhuang'

import unittest

TESTDIR = os.path.dirname(__file__)
TESTMODEL = load_model(os.path.join(TESTDIR, 'data/EcoliCore.xml'), sanitize=False)


class OrganismTest(unittest.TestCase):
    def setUp(self):
        self.ec1 = Organism(TESTMODEL)

    def test_initialization(self):
        self.assertEqual(TESTMODEL, self.ec1.model)
        self.assertEqual(TESTMODEL.id, self.ec1.model.id)

    def test_update(self):
        self.assertRaises(NotImplementedError, self.ec1.update)
        self.ec1.update = update_organism
        self.assertTrue(self.ec1.update(self.ec1) == self.ec1)

    def tearDown(self):
        del self.ec1


class EnvironmentTest(unittest.TestCase):
    def setUp(self):
        self.o1 = Organism(TESTMODEL)
        self.o2 = Organism(TESTMODEL)
        self.env = Environment()

    def test_initialization(self):
        assert type(self.env) == Environment

    def test_add_organisms(self):
        self.env._add_organism(self.o1)
        self.assertEqual(self.env.organisms, [self.o1])
        self.env._add_organism(self.o2)
        self.assertEqual(self.env.organisms, [self.o1, self.o2])
        self.env._add_organism(self.o1)
        self.env._add_organism(self.o2)
        self.assertEqual(self.env.organisms, [self.o1, self.o2, self.o1, self.o2])

    def test_add_metabolites(self):
        self.env.metabolites.append('EX_glc(e)')
        self.assertEqual(self.env.metabolites, ['EX_glc(e)'])
        self.env.metabolites += ['EX_ac(e)', 'EX_o2(e)']
        self.assertEqual(self.env.metabolites, ['EX_glc(e)', 'EX_ac(e)', 'EX_o2(e)'])

    def tearDown(self):
        del self.o1
        del self.o2
        del self.env


class BioreactorTest(unittest.TestCase):
    def setUp(self):
        self.o1 = Organism(TESTMODEL)
        self.o2 = Organism(TESTMODEL)
        self.br = BioReactor([self.o1, self.o2], ['R_EX_glc_e', 'R_EX_ac_e', 'R_EX_o2_e'])
        self.br2 = BioReactor()

    def test_initialization(self):
        self.assertEqual(self.br.organisms,[self.o1, self.o2])
        self.assertEqual(self.br.metabolites, ['R_EX_glc_e', 'R_EX_ac_e', 'R_EX_o2_e'])
        self.assertEqual(len(self.br.organisms), 2)

    def test_set_organisms(self):
        self.br2.organisms = [self.o1, self.o2]
        assert self.br2.organisms == self.br.organisms

    def test_set_metabolites(self):
        self.br2.metabolites = ['R_EX_glc_e', 'R_EX_ac_e', 'R_EX_o2_e']
        self.assertEqual(self.br2.metabolites, self.br.metabolites)

    def test_set_initial_conditions(self):
        self.br.initial_conditions = [1, 0.1, 0.1, 10, 1, 0]
        self.assertEqual(self.br.initial_conditions, [1, 0.1, 0.1, 10, 1, 0])

    def tearDown(self):
        del self.br
        del self.o1
        del self.o2


def update_organism(self):
    """
    fixture function for testing Organism.update()
    """
    return self


def suite():
    tests = [OrganismTest, EnvironmentTest, BioreactorTest]

    test_suite = unittest.TestSuite()
    for test in tests:
        test_suite.addTest(unittest.makeSuite(test))
    return test_suite


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())
