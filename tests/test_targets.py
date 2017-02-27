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
from cameo.core.target import *
from cameo.core.target import Target, FluxModulationTarget, EnsembleTarget, ReactionInversionTarget
from cameo.exceptions import IncompatibleTargets
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
        fold_change = -0.30041

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

        fold_change = -0.342537

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

    def test_reaction_inversion_target(self):
        inversion_target = ReactionInversionTarget("GND", value=-10, reference_value=10)
        self.assertEqual(inversion_target.fold_change, 0)
        lower_bound = self.model.reactions.GND.lower_bound
        with TimeMachine() as tm:
            inversion_target.apply(self.model, time_machine=tm)
            self.assertEqual(self.model.reactions.GND.lower_bound, -10)
        self.assertEqual(self.model.reactions.GND.lower_bound, lower_bound)

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

    @unittest.skip("Gene Overexpression not implemented yet")
    def test_gene_over_express_target(self):
        raise NotImplementedError

    @unittest.skip("Gene Downregulation not implemented yet")
    def test_gene_down_regulation_target(self):
        raise NotImplementedError

    @unittest.skipIf(six.PY2, 'gnomic is not compatible with python 2')
    def test_gnomic_integration(self):
        from gnomic.models import Accession, Feature, Mutation, FeatureTree
        abstract_target = Target("test")
        abstract_target_gnomic = abstract_target.to_gnomic()
        self.assertIsInstance(abstract_target_gnomic, Accession)
        self.assertEqual(abstract_target_gnomic.identifier, abstract_target.id)

        flux_modulation_target = FluxModulationTarget("test", 1, 0)
        flux_modulation_target_gnomic = flux_modulation_target.to_gnomic()
        self.assertIsInstance(flux_modulation_target_gnomic, Mutation)
        self.assertIsInstance(flux_modulation_target_gnomic.old, FeatureTree)
        self.assertIsInstance(flux_modulation_target_gnomic.old[0], Feature)
        self.assertEqual(flux_modulation_target_gnomic.old[0].accession.identifier, flux_modulation_target.id)
        self.assertEqual(flux_modulation_target_gnomic.old[0].variant, None)
        self.assertEqual(flux_modulation_target_gnomic.old[0].type, 'flux')
        self.assertIsInstance(flux_modulation_target_gnomic.new, FeatureTree)
        self.assertIsInstance(flux_modulation_target_gnomic.new[0], Feature)
        self.assertEqual(flux_modulation_target_gnomic.new[0].accession.identifier, flux_modulation_target.id)
        self.assertEqual(flux_modulation_target_gnomic.new[0].type, 'flux')
        self.assertEqual(flux_modulation_target_gnomic.new[0].variant, "over-expression(%.3f)" % flux_modulation_target.fold_change)

        flux_modulation_target = FluxModulationTarget("test", 0.5, 1)
        flux_modulation_target_gnomic = flux_modulation_target.to_gnomic()
        self.assertIsInstance(flux_modulation_target_gnomic, Mutation)
        self.assertIsInstance(flux_modulation_target_gnomic.old, FeatureTree)
        self.assertIsInstance(flux_modulation_target_gnomic.old[0], Feature)
        self.assertEqual(flux_modulation_target_gnomic.old[0].accession.identifier, flux_modulation_target.id)
        self.assertEqual(flux_modulation_target_gnomic.old[0].variant, None)
        self.assertEqual(flux_modulation_target_gnomic.old[0].type, 'flux')
        self.assertIsInstance(flux_modulation_target_gnomic.new, FeatureTree)
        self.assertIsInstance(flux_modulation_target_gnomic.new[0], Feature)
        self.assertEqual(flux_modulation_target_gnomic.new[0].accession.identifier, flux_modulation_target.id)
        self.assertEqual(flux_modulation_target_gnomic.new[0].type, 'flux')
        self.assertEqual(flux_modulation_target_gnomic.new[0].variant, "down-regulation(%.3f)" % flux_modulation_target.fold_change)

        flux_modulation_target = FluxModulationTarget("test", 0, 1)
        flux_modulation_target_gnomic = flux_modulation_target.to_gnomic()
        self.assertIsInstance(flux_modulation_target_gnomic, Mutation)
        self.assertIsInstance(flux_modulation_target_gnomic.old, FeatureTree)
        self.assertIsInstance(flux_modulation_target_gnomic.old[0], Feature)
        self.assertEqual(flux_modulation_target_gnomic.old[0].accession.identifier, flux_modulation_target.id)
        self.assertEqual(flux_modulation_target_gnomic.old[0].variant, None)
        self.assertEqual(flux_modulation_target_gnomic.old[0].type, 'flux')
        self.assertIs(flux_modulation_target_gnomic.new, None)

        reaction = Reaction(id="atpzase", name="Cosmic ATP generator")
        atp_z = Metabolite(id="atp_z", name="Cosmic ATP", compartment="c")

        reaction.add_metabolites({self.model.metabolites.atp_c: 1, atp_z: -1})
        knockin_target = ReactionKnockinTarget("atpzase", reaction)
        knockin_target_gnomic = knockin_target.to_gnomic()
        self.assertIsInstance(knockin_target_gnomic, Mutation)
        self.assertIsInstance(knockin_target_gnomic.new, FeatureTree)
        self.assertIsInstance(knockin_target_gnomic.new[0], Feature)
        self.assertEqual(knockin_target_gnomic.new[0].accession.identifier, knockin_target.id)
        self.assertEqual(knockin_target_gnomic.new[0].variant, None)
        self.assertEqual(knockin_target_gnomic.new[0].type, 'reaction')
        self.assertIs(knockin_target_gnomic.old, None)

        cofactor_id_swaps = [("nad_c", "nadh_c"), ("nadp_c", "nadph_c")]

        swap_pairs = ([self.model.metabolites.get_by_id(m) for m in cofactor_id_swaps[0]],
                      [self.model.metabolites.get_by_id(m) for m in cofactor_id_swaps[1]])

        swap_target = ReactionCofactorSwapTarget("GAPD", swap_pairs)
        swap_target_gnomic = swap_target.to_gnomic()

        self.assertIsInstance(swap_target_gnomic, Mutation)
        self.assertIsInstance(swap_target_gnomic.old, FeatureTree)
        self.assertIsInstance(swap_target_gnomic.old[0], Feature)
        self.assertEqual(swap_target_gnomic.old[0].accession.identifier, swap_target.id)
        self.assertEqual(swap_target_gnomic.old[0].variant, None)
        self.assertEqual(swap_target_gnomic.old[0].type, 'reaction')
        self.assertIsInstance(swap_target_gnomic.new, FeatureTree)
        self.assertIsInstance(swap_target_gnomic.new[0], Feature)
        self.assertEqual(swap_target_gnomic.new[0].accession.identifier, swap_target.id + swap_target.swap_str)
        self.assertEqual(swap_target_gnomic.new[0].variant, None)
        self.assertEqual(swap_target_gnomic.new[0].type, 'reaction')


class EnsembleTargetsTestCase(unittest.TestCase):
    def test_compatible_targets(self):
        modulation_target = ReactionModulationTarget("a", 1, 0)
        ki_target = ReactionKnockinTarget("a", None)

        ensemble_1 = EnsembleTarget("a", [modulation_target, ki_target])
        ensemble_2 = EnsembleTarget("a", [ki_target, modulation_target])

        self.assertEqual(ensemble_1.targets, ensemble_2.targets)

        modulation_target = ReactionModulationTarget("b", 1, 0)
        swap_target = ReactionCofactorSwapTarget("b", [("nad_c", "nadh_c"), ("nadp_c", "nadph_c")])

        ensemble_1 = EnsembleTarget("b", [modulation_target, swap_target])
        ensemble_2 = EnsembleTarget("b", [swap_target, modulation_target])

        self.assertEqual(ensemble_1.targets, ensemble_2.targets)

        ki_target = ReactionKnockinTarget("c", None)
        modulation_target = ReactionModulationTarget("c", 1, 0)
        swap_target = ReactionCofactorSwapTarget("c", [("nad_c", "nadh_c"), ("nadp_c", "nadph_c")])

        ensemble = EnsembleTarget("c", [modulation_target, swap_target, ki_target])
        self.assertEqual(ensemble.targets[0], ki_target)
        self.assertEqual(ensemble.targets[1], swap_target)
        self.assertEqual(ensemble.targets[2], modulation_target)

    def test_incompatible_targets(self):
        ko_target = ReactionKnockoutTarget("a")
        ki_target = ReactionKnockinTarget("a", None)

        self.assertRaises(IncompatibleTargets, EnsembleTarget, "a", [ko_target, ki_target])
        self.assertRaises(IncompatibleTargets, EnsembleTarget, "a", [ki_target, ko_target])

        ko_target = ReactionKnockoutTarget("b")
        swap_target = ReactionCofactorSwapTarget("b", [("nad_c", "nadh_c"), ("nadp_c", "nadph_c")])

        self.assertRaises(IncompatibleTargets, EnsembleTarget, "b", [ko_target, swap_target])
        self.assertRaises(IncompatibleTargets, EnsembleTarget, "b", [swap_target, ko_target])

        modulation_target = ReactionModulationTarget("c", 0, 0)
        ki_target = ReactionKnockinTarget("c", None)

        self.assertRaises(IncompatibleTargets, EnsembleTarget, "c", [modulation_target, ki_target])
        self.assertRaises(IncompatibleTargets, EnsembleTarget, "c", [ki_target, modulation_target])

