# Copyright 2016 The Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software


# See the License for the specific language governing permissions and
# limitations under the License.
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
import pytest


from gnomic.genotype import Genotype

from gnomic.types import Accession, Feature, Change
from gnomic.utils import genotype_to_string, genotype_to_text
from cobra import Metabolite, Reaction

from cameo.core.target import (ReactionKnockoutTarget, ReactionModulationTarget, GeneKnockoutTarget,
                               ReactionCofactorSwapTarget, ReactionKnockinTarget)
from cameo.core.target import (EnsembleTarget, FluxModulationTarget,
                               ReactionInversionTarget, Target)
from cameo.exceptions import IncompatibleTargets


class TestTargets:
    def test_hashable(self):
        knockout_target1 = ReactionKnockoutTarget("ACALD")
        knockout_target2 = ReactionKnockoutTarget("ACALD")

        assert knockout_target1 == knockout_target2
        assert hash(knockout_target1) == hash(knockout_target2)

        assert len({knockout_target1, knockout_target2}) == 1

    def test_reaction_knockout_target(self, model):
        knockout_target = ReactionKnockoutTarget("ACALD")
        with model:
            knockout_target.apply(model)
            assert model.reactions.ACALD.lower_bound == 0
            assert model.reactions.ACALD.upper_bound == 0

        assert model.reactions.ACALD.lower_bound == -1000
        assert model.reactions.ACALD.upper_bound == 1000

    def test_reaction_down_regulation_target(self, model):
        reaction_id = "PGI"
        ref_val = 4.86
        value = 3.4

        # (B - A) / A
        fold_change = -0.30041

        down_reg_target = ReactionModulationTarget(reaction_id, value, ref_val)
        assert round(abs(down_reg_target.fold_change - fold_change), 5) == 0
        with model:
            down_reg_target.apply(model)
            assert model.reactions.PGI.upper_bound == 3.4
            assert model.reactions.PGI.lower_bound == -1000
            assert abs(model.slim_optimize() - 0.8706) < 0.0001

        assert model.reactions.PGI.upper_bound == 1000
        assert model.reactions.PGI.lower_bound == -1000

        reaction_id = "RPI"
        ref_val = -2.28150
        value = -1.5

        fold_change = -0.342537

        down_reg_target = ReactionModulationTarget(reaction_id, value, ref_val)
        assert round(abs(down_reg_target.fold_change - fold_change), 5) == 0
        with model:
            down_reg_target.apply(model)
            assert model.reactions.RPI.lower_bound == -1.5
            assert model.reactions.RPI.upper_bound == 1000
            assert abs(model.slim_optimize() - 0.8691) < 0.0001

        assert model.reactions.RPI.lower_bound == -1000
        assert model.reactions.RPI.upper_bound == 1000

    def test_reaction_knock_in_target(self, model):
        reaction = Reaction(id="atpzase", name="Cosmic ATP generator")
        atp_z = Metabolite(id="atp_z", name="Cosmic ATP", compartment="c")

        reaction.add_metabolites({model.metabolites.atp_c: 1, atp_z: -1})
        knockin_target = ReactionKnockinTarget("atpzase", reaction)
        with model:
            knockin_target.apply(model)
            assert atp_z in model.metabolites
            assert reaction in model.reactions

        assert atp_z not in model.metabolites
        assert reaction not in model.reactions

    def test_reaction_cofactor_swap_target(self, model):
        cofactor_id_swaps = [("nad_c", "nadh_c"), ("nadp_c", "nadph_c")]

        swap_pairs = ([model.metabolites.get_by_id(m) for m in cofactor_id_swaps[0]],
                      [model.metabolites.get_by_id(m) for m in cofactor_id_swaps[1]])

        swap_target = ReactionCofactorSwapTarget("GAPD", swap_pairs)
        with model:
            swap_target.apply(model)
            assert model.metabolites.nad_c not in model.reactions.GAPD.metabolites
            assert model.metabolites.nadh_c not in model.reactions.GAPD.metabolites
            assert model.metabolites.nadp_c in model.reactions.GAPD.metabolites
            assert model.metabolites.nadph_c in model.reactions.GAPD.metabolites

        assert model.metabolites.nadp_c not in model.reactions.GAPD.metabolites
        assert model.metabolites.nadph_c not in model.reactions.GAPD.metabolites
        assert model.metabolites.nad_c in model.reactions.GAPD.metabolites
        assert model.metabolites.nadh_c in model.reactions.GAPD.metabolites

        swap_target = ReactionCofactorSwapTarget("GND", swap_pairs)
        with model:
            swap_target.apply(model)
            assert model.metabolites.nad_c in model.reactions.GND.metabolites
            assert model.metabolites.nadh_c in model.reactions.GND.metabolites
            assert model.metabolites.nadp_c not in model.reactions.GND.metabolites
            assert model.metabolites.nadph_c not in model.reactions.GND.metabolites

        assert model.metabolites.nadp_c in model.reactions.GND.metabolites
        assert model.metabolites.nadph_c in model.reactions.GND.metabolites
        assert model.metabolites.nad_c not in model.reactions.GND.metabolites
        assert model.metabolites.nadh_c not in model.reactions.GND.metabolites

    def test_reaction_inversion_target(self, model):
        inversion_target = ReactionInversionTarget("GND", value=-10, reference_value=10)
        assert inversion_target.fold_change == 0
        lower_bound = model.reactions.GND.lower_bound
        with model:
            inversion_target.apply(model)
            assert model.reactions.GND.lower_bound == -10
        assert model.reactions.GND.lower_bound == lower_bound

    def test_invalid_reaction_knockout_target(self, model):
        knockout_target = ReactionKnockoutTarget("ACALDXYZ")

        with pytest.raises(KeyError):
            knockout_target.apply(model)

    def test_invalid_reaction_modulation_target(self, model):
        reaction_id = "PGI_XY"
        ref_val = 4.86
        value = 4

        down_reg_target = ReactionModulationTarget(reaction_id, value, ref_val)

        with pytest.raises(KeyError):
            down_reg_target.apply(model)

        reaction_id = "RPI_Z"
        ref_val = -2.28150
        value = -2.0

        down_reg_target = ReactionModulationTarget(reaction_id, value, ref_val)

        with pytest.raises(KeyError):
            down_reg_target.apply(model)

    def test_invalid_reaction_cofactor_swap_target(self, model):
        cofactor_id_swaps = [("nad_c", "nadh_c"), ("nadp_c", "nadph_c")]

        swap_pairs = ([model.metabolites.get_by_id(m) for m in cofactor_id_swaps[0]],
                      [model.metabolites.get_by_id(m) for m in cofactor_id_swaps[1]])

        swap_target = ReactionCofactorSwapTarget("GAPD_124", swap_pairs)

        with pytest.raises(KeyError):
            swap_target.apply(model)

        swap_target = ReactionCofactorSwapTarget("ACKr", swap_pairs)

        with pytest.raises(ValueError):
            swap_target.apply(model)

    def test_gene_knockout_target(self, model):
        gene = "b4025"
        knockout_target = GeneKnockoutTarget(gene, accession_id=gene, accession_db='bigg')
        knockout_target_gnomic = knockout_target.to_gnomic()
        assert genotype_to_string(Genotype([knockout_target_gnomic])) == "-b4025#bigg:b4025"
        with model:
            knockout_target.apply(model)
            assert model.reactions.PGI.lower_bound == 0
            assert model.reactions.PGI.upper_bound == 0
            assert abs(model.slim_optimize() - 0.8631) < 0.0001

        assert model.reactions.PGI.lower_bound == -1000
        assert model.reactions.PGI.upper_bound == 1000

    @pytest.mark.skip("Gene Overexpression not implemented yet")
    def test_gene_over_express_target(self):
        raise NotImplementedError

    @pytest.mark.skip("Gene Downregulation not implemented yet")
    def test_gene_down_regulation_target(self):
        raise NotImplementedError

    def test_gnomic_integration_Target(self, model):
        abstract_target = Target("test")
        abstract_target_gnomic = abstract_target.to_gnomic()
        assert abstract_target_gnomic is None
        abstract_target = Target("test", accession_id='test', accession_db='test')
        abstract_target_gnomic = abstract_target.to_gnomic()
        assert isinstance(abstract_target_gnomic, Accession)
        assert abstract_target_gnomic.identifier == abstract_target.id

    def test_gnomic_integration_FluxModulationTarget(self, model):
        # with pytest.raises(ValueError):
        #     FluxModulationTarget("test", 0, 0)  # TODO: this should really not be possible
        flux_modulation_target = FluxModulationTarget("test", 1, 0)
        flux_modulation_target_gnomic = flux_modulation_target.to_gnomic()
        expected = "flux.test(value=1)"
        assert genotype_to_string(Genotype([flux_modulation_target_gnomic])) == expected

        flux_modulation_target = FluxModulationTarget("PGK", 0.5, 1, accession_id="PGK", accession_db="bigg")
        flux_modulation_target_gnomic = flux_modulation_target.to_gnomic()
        expected = "flux.PGK#bigg:PGK(value=0.5)"
        assert genotype_to_string(Genotype([flux_modulation_target_gnomic])) == expected

    def test_gnomic_integration_ReactionKnockinTarget(self, model):
        reaction = Reaction(id="atpzase", name="Cosmic ATP generator")
        atp_z = Metabolite(id="atp_z", name="Cosmic ATP", compartment="c")

        reaction.add_metabolites({model.metabolites.atp_c: 1, atp_z: -1})
        knockin_target = ReactionKnockinTarget("atpzase", reaction)
        knockin_target_gnomic = knockin_target.to_gnomic()
        assert genotype_to_string(Genotype([knockin_target_gnomic])) == "+reaction.atpzase"

        reaction.add_metabolites({model.metabolites.atp_c: 1, atp_z: -1})
        knockin_target = ReactionKnockinTarget("atpzase", reaction, accession_id='atpzase', accession_db='unicorn')
        knockin_target_gnomic = knockin_target.to_gnomic()
        assert genotype_to_string(Genotype([knockin_target_gnomic])) == "+reaction.atpzase#unicorn:atpzase"

    def test_gnomic_integration_ReactionCofactorSwapTarget(self, model):

        cofactor_id_swaps = [("nad_c", "nadh_c"), ("nadp_c", "nadph_c")]

        swap_pairs = ([model.metabolites.get_by_id(m) for m in cofactor_id_swaps[0]],
                      [model.metabolites.get_by_id(m) for m in cofactor_id_swaps[1]])

        swap_target = ReactionCofactorSwapTarget("GAPD", swap_pairs)
        swap_target_gnomic = swap_target.to_gnomic()
        expected = "reaction.GAPD(cofactors=nadp_c,nadph_c)"
        assert genotype_to_string(Genotype([swap_target_gnomic])) == expected


class TestEnsembleTargets:
    def test_compatible_targets(self):
        modulation_target = ReactionModulationTarget("a", 1, 0)
        ki_target = ReactionKnockinTarget("a", None)

        ensemble_1 = EnsembleTarget("a", [modulation_target, ki_target])
        ensemble_2 = EnsembleTarget("a", [ki_target, modulation_target])

        assert ensemble_1.targets == ensemble_2.targets

        modulation_target = ReactionModulationTarget("b", 1, 0)
        swap_target = ReactionCofactorSwapTarget("b", [("nad_c", "nadh_c"), ("nadp_c", "nadph_c")])

        ensemble_1 = EnsembleTarget("b", [modulation_target, swap_target])
        ensemble_2 = EnsembleTarget("b", [swap_target, modulation_target])

        assert ensemble_1.targets == ensemble_2.targets

        ki_target = ReactionKnockinTarget("c", None)
        modulation_target = ReactionModulationTarget("c", 1, 0)
        swap_target = ReactionCofactorSwapTarget("c", [("nad_c", "nadh_c"), ("nadp_c", "nadph_c")])

        ensemble = EnsembleTarget("c", [modulation_target, swap_target, ki_target])
        assert ensemble.targets[0] == ki_target
        assert ensemble.targets[1] == swap_target
        assert ensemble.targets[2] == modulation_target

    def test_incompatible_targets(self):
        ko_target = ReactionKnockoutTarget("a")
        ki_target = ReactionKnockinTarget("a", None)

        with pytest.raises(IncompatibleTargets):
            EnsembleTarget("a", [ko_target, ki_target])
        with pytest.raises(IncompatibleTargets):
            EnsembleTarget("a", [ki_target, ko_target])

        ko_target = ReactionKnockoutTarget("b")
        swap_target = ReactionCofactorSwapTarget("b", [("nad_c", "nadh_c"), ("nadp_c", "nadph_c")])

        with pytest.raises(IncompatibleTargets):
            EnsembleTarget("b", [ko_target, swap_target])
        with pytest.raises(IncompatibleTargets):
            EnsembleTarget("b", [swap_target, ko_target])

        modulation_target = ReactionModulationTarget("c", 0, 0)
        ki_target = ReactionKnockinTarget("c", None)

        with pytest.raises(IncompatibleTargets):
            EnsembleTarget("c", [modulation_target, ki_target])
        with pytest.raises(IncompatibleTargets):
            EnsembleTarget("c", [ki_target, modulation_target])
