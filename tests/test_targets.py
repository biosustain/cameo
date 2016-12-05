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

from cameo import load_model, Reaction, Metabolite
from cameo.core.target import *
from cameo.util import TimeMachine

TESTDIR = os.path.dirname(__file__)


class TargetsTestCase(unittest.TestCase):
    model = None

    @classmethod
    def setUpClass(cls):
        cls.model = load_model(os.path.join(TESTDIR, 'data', 'EcoliCore.xml'))

    def test_reaction_knockout_target(self):
        knockout_target = ReactionKnockoutTarget("ACALD")
        with TimeMachine() as tm:
            knockout_target.apply(self.model, time_machine=tm)
            self.assertEqual(self.model.reactions.ACALD.lower_bound, 0)
            self.assertEqual(self.model.reactions.ACALD.upper_bound, 0)

        self.assertEqual(self.model.reactions.ACALD.lower_bound, -1000)
        self.assertEqual(self.model.reactions.ACALD.upper_bound, 1000)

    def test_reaction_down_regulation_target(self):
        reaction_id = "PGI"
        ref_val = 4.86
        value = 3.4

        # (B - A) / A
        fold_change = 0.300411

        down_reg_target = ReactionModulationTarget(reaction_id, value, ref_val)
        self.assertAlmostEqual(down_reg_target.fold_change, fold_change, places=5)
        with TimeMachine() as tm:
            down_reg_target.apply(self.model, time_machine=tm)
            self.assertEqual(self.model.reactions.PGI.upper_bound, 3.4)
            self.assertEqual(self.model.reactions.PGI.lower_bound, -1000)
            self.assertAlmostEqual(self.model.solve().f, 0.8706, delta=0.0001)

        self.assertEqual(self.model.reactions.PGI.upper_bound, 1000)
        self.assertEqual(self.model.reactions.PGI.lower_bound, -1000)

        reaction_id = "RPI"
        ref_val = -2.28150
        value = -1.5

        fold_change = 0.342537

        down_reg_target = ReactionModulationTarget(reaction_id, value, ref_val)
        self.assertAlmostEqual(down_reg_target.fold_change, fold_change, places=5)
        with TimeMachine() as tm:
            down_reg_target.apply(self.model, time_machine=tm)
            self.assertEqual(self.model.reactions.RPI.lower_bound, -1.5)
            self.assertEqual(self.model.reactions.RPI.upper_bound, 1000)
            self.assertAlmostEqual(self.model.solve().f, 0.8691, delta=0.0001)

        self.assertEqual(self.model.reactions.RPI.lower_bound, -1000)
        self.assertEqual(self.model.reactions.RPI.upper_bound, 1000)

    def test_reaction_knock_in_target(self):
        reaction = Reaction(id="atpzase", name="Cosmic ATP generator")
        atp_z = Metabolite(id="atp_z", name="Cosmic ATP", compartment="c")

        reaction.add_metabolites({self.model.metabolites.atp_c: 1, atp_z: -1})
        knockin_target = ReactionKnockinTarget("atpzase", reaction)
        with TimeMachine() as tm:
            knockin_target.apply(self.model, time_machine=tm)
            self.assertIn(atp_z, self.model.metabolites)
            self.assertIn(reaction, self.model.reactions)

        self.assertNotIn(atp_z, self.model.metabolites)
        self.assertNotIn(reaction, self.model.reactions)

    def test_reaction_cofactor_swap_target(self):
        cofactor_id_swaps = [("nad_c", "nadh_c"), ("nadp_c", "nadph_c")]

        swap_pairs = ([self.model.metabolites.get_by_id(m) for m in cofactor_id_swaps[0]],
                      [self.model.metabolites.get_by_id(m) for m in cofactor_id_swaps[1]])

        swap_target = ReactionCofactorSwapTarget("GAPD", swap_pairs)
        with TimeMachine() as tm:
            swap_target.apply(self.model, time_machine=tm)
            self.assertNotIn(self.model.metabolites.nad_c, self.model.reactions.GAPD.metabolites)
            self.assertNotIn(self.model.metabolites.nadh_c, self.model.reactions.GAPD.metabolites)
            self.assertIn(self.model.metabolites.nadp_c, self.model.reactions.GAPD.metabolites)
            self.assertIn(self.model.metabolites.nadph_c, self.model.reactions.GAPD.metabolites)

        self.assertNotIn(self.model.metabolites.nadp_c, self.model.reactions.GAPD.metabolites)
        self.assertNotIn(self.model.metabolites.nadph_c, self.model.reactions.GAPD.metabolites)
        self.assertIn(self.model.metabolites.nad_c, self.model.reactions.GAPD.metabolites)
        self.assertIn(self.model.metabolites.nadh_c, self.model.reactions.GAPD.metabolites)

        swap_target = ReactionCofactorSwapTarget("GND", swap_pairs)
        with TimeMachine() as tm:
            swap_target.apply(self.model, time_machine=tm)
            self.assertIn(self.model.metabolites.nad_c, self.model.reactions.GND.metabolites)
            self.assertIn(self.model.metabolites.nadh_c, self.model.reactions.GND.metabolites)
            self.assertNotIn(self.model.metabolites.nadp_c, self.model.reactions.GND.metabolites)
            self.assertNotIn(self.model.metabolites.nadph_c, self.model.reactions.GND.metabolites)

        self.assertIn(self.model.metabolites.nadp_c, self.model.reactions.GND.metabolites)
        self.assertIn(self.model.metabolites.nadph_c, self.model.reactions.GND.metabolites)
        self.assertNotIn(self.model.metabolites.nad_c, self.model.reactions.GND.metabolites)
        self.assertNotIn(self.model.metabolites.nadh_c, self.model.reactions.GND.metabolites)

    def test_invalid_reaction_knockout_target(self):
        knockout_target = ReactionKnockoutTarget("ACALDXYZ")

        self.assertRaises(KeyError, knockout_target.apply, self.model)

    def test_invalid_reaction_modulation_target(self):

        reaction_id = "PGI_XY"
        ref_val = 4.86
        value = 4

        down_reg_target = ReactionModulationTarget(reaction_id, value, ref_val)

        self.assertRaises(KeyError, down_reg_target.apply, self.model)

        reaction_id = "RPI_Z"
        ref_val = -2.28150
        value = -2.0

        down_reg_target = ReactionModulationTarget(reaction_id, value, ref_val)

        self.assertRaises(KeyError, down_reg_target.apply, self.model)

    def test_invalid_reaction_cofactor_swap_target(self):
        cofactor_id_swaps = [("nad_c", "nadh_c"), ("nadp_c", "nadph_c")]

        swap_pairs = ([self.model.metabolites.get_by_id(m) for m in cofactor_id_swaps[0]],
                      [self.model.metabolites.get_by_id(m) for m in cofactor_id_swaps[1]])

        swap_target = ReactionCofactorSwapTarget("GAPD_124", swap_pairs)

        self.assertRaises(KeyError, swap_target.apply, self.model)

        swap_target = ReactionCofactorSwapTarget("ACKr", swap_pairs)

        self.assertRaises(ValueError, swap_target.apply, self.model)

    def test_gene_knockout_target(self):
        gene = "b4025"

        knockout_target = GeneKnockoutTarget(gene)
        with TimeMachine() as tm:
            knockout_target.apply(self.model, time_machine=tm)
            self.assertEqual(self.model.reactions.PGI.lower_bound, 0)
            self.assertEqual(self.model.reactions.PGI.upper_bound, 0)
            self.assertAlmostEqual(self.model.solve().f, 0.8631, delta=0.0001)

        self.assertEqual(self.model.reactions.PGI.lower_bound, -1000)
        self.assertEqual(self.model.reactions.PGI.upper_bound, 1000)

    @unittest.skip("Gene Overexpression not implemente yet")
    def test_gene_over_express_target(self):
        raise NotImplementedError

    @unittest.skip("Gene Downregulation not implemente yet")
    def test_gene_down_regulation_target(self):
        raise NotImplementedError