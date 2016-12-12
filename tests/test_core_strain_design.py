# Copyright 2016 The Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import os
import unittest

import six

from cameo import load_model, Reaction, Metabolite
from cameo.core.strain_design import StrainDesign
from cameo.core.target import ReactionKnockoutTarget, ReactionKnockinTarget, ReactionModulationTarget
from cameo.exceptions import IncompatibleTargets

TESTDIR = os.path.dirname(__file__)


class StrainDesignAbstractionTestCase(unittest.TestCase):
    model = None
    cad_reaction = None

    @classmethod
    def setUpClass(cls):
        cls.model = load_model(os.path.join(TESTDIR, 'data', 'EcoliCore.xml'))
        cls.cad_reaction = Reaction(id="CAD", name="Cis-Aconitate Decarboxylase")
        acon_C_c = cls.model.metabolites.acon_dsh_C_c
        co2_c = cls.model.metabolites.co2_c
        ita_c = Metabolite(id="ita_c", name="Itaconate", compartment="c")
        cls.cad_reaction.add_metabolites({acon_C_c: -1, co2_c: 1, ita_c: 1})

    def test_create_strain_design(self):
        t1 = ReactionKnockoutTarget('PGI')
        t2 = ReactionKnockoutTarget('GAPD')
        t3 = ReactionKnockinTarget("CAD", self.cad_reaction)

        strain_design = StrainDesign([t1, t2, t3])

        self.assertEqual(len(strain_design), 3)
        strain_design2 = StrainDesign([t1, t2, t3])
        strain_design3 = StrainDesign([t2, t1, t3])

        self.assertEqual(strain_design, strain_design2)
        self.assertEqual(strain_design, strain_design3)
        self.assertEqual(strain_design3, strain_design2)

        self.assertIn(t1, strain_design)
        self.assertIn(t2, strain_design)
        self.assertIn(t3, strain_design)

    def test_add_strain_design(self):
        t1 = ReactionKnockoutTarget('PGI')
        t2 = ReactionKnockoutTarget('GAPD')
        t3 = ReactionKnockinTarget("CAD", self.cad_reaction)

        strain_design1 = StrainDesign([t1, t2, t3])

        t4 = ReactionModulationTarget("PGI", 5, 1)

        strain_design2 = StrainDesign([t4])

        self.assertRaises(IncompatibleTargets, strain_design1.__add__, strain_design2)
        self.assertRaises(IncompatibleTargets, strain_design2.__add__, strain_design1)

        self.assertRaises(IncompatibleTargets, strain_design1.__iadd__, strain_design2)
        self.assertRaises(IncompatibleTargets, strain_design2.__iadd__, strain_design1)

        t5 = ReactionModulationTarget("RPI", 2, 0)
        strain_design3 = StrainDesign([t5])

        strain_design4 = strain_design3 + strain_design1
        self.assertIn(t1, strain_design4)
        self.assertIn(t2, strain_design4)
        self.assertIn(t3, strain_design4)
        self.assertNotIn(t4, strain_design4)
        self.assertIn(t5, strain_design4)

        strain_design3 += strain_design1

        self.assertIn(t1, strain_design3)
        self.assertIn(t2, strain_design3)
        self.assertIn(t3, strain_design3)
        self.assertNotIn(t4, strain_design3)
        self.assertIn(t5, strain_design3)

    @unittest.skipIf(six.PY2, "Gnomic is not compatible with python 2")
    def test_design_to_gnomic(self):
        from gnomic import Genotype
        t1 = ReactionKnockoutTarget('PGI')
        t2 = ReactionKnockoutTarget('GAPD')
        t3 = ReactionKnockinTarget("CAD", self.cad_reaction)

        strain_design1 = StrainDesign([t1, t2, t3])

        sd_gnomic = strain_design1.to_gnomic()

        self.assertIsInstance(sd_gnomic, Genotype)
        self.assertEqual(len(sd_gnomic.added_features), 1)
        self.assertEqual(len(sd_gnomic.removed_features), 2)
