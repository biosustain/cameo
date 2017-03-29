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
import six

from cobra import Metabolite, Reaction

from cameo.core.strain_design import StrainDesign
from cameo.core.target import (ReactionKnockinTarget, ReactionKnockoutTarget,
                               ReactionModulationTarget)
from cameo.exceptions import IncompatibleTargets


@pytest.fixture(scope="function")
def cad_reaction(core_model):
    reaction = Reaction(id="CAD", name="Cis-Aconitate Decarboxylase")
    acon = core_model.metabolites.acon_DASH_C_c
    co2_c = core_model.metabolites.co2_c
    ita_c = Metabolite(id="ita_c", name="Itaconate", compartment="c")
    reaction.add_metabolites({acon: -1, co2_c: 1, ita_c: 1})
    return reaction


class TestStrainDesign:

    def test_create_strain_design(self, cad_reaction):
        t1 = ReactionKnockoutTarget('PGI')
        t2 = ReactionKnockoutTarget('GAPD')
        t3 = ReactionKnockinTarget("CAD", cad_reaction)

        strain_design = StrainDesign([t1, t2, t3])

        assert len(strain_design) == 3
        strain_design2 = StrainDesign([t1, t2, t3])
        strain_design3 = StrainDesign([t2, t1, t3])

        assert strain_design == strain_design2
        assert strain_design == strain_design3
        assert strain_design3 == strain_design2

        assert t1 in strain_design
        assert t2 in strain_design
        assert t3 in strain_design

    def test_add_strain_design(self, cad_reaction):
        t1 = ReactionKnockoutTarget('PGI')
        t2 = ReactionKnockoutTarget('GAPD')
        t3 = ReactionKnockinTarget("CAD", cad_reaction)

        strain_design1 = StrainDesign([t1, t2, t3])

        t4 = ReactionModulationTarget("PGI", 5, 1)

        strain_design2 = StrainDesign([t4])

        with pytest.raises(IncompatibleTargets):
            strain_design1.__add__(strain_design2)
        with pytest.raises(IncompatibleTargets):
            strain_design2.__add__(strain_design1)

        with pytest.raises(IncompatibleTargets):
            strain_design1.__iadd__(strain_design2)
        with pytest.raises(IncompatibleTargets):
            strain_design2.__iadd__(strain_design1)

        t5 = ReactionModulationTarget("RPI", 2, 0)
        strain_design3 = StrainDesign([t5])

        strain_design4 = strain_design3 + strain_design1
        assert t1 in strain_design4
        assert t2 in strain_design4
        assert t3 in strain_design4
        assert t4 not in strain_design4
        assert t5 in strain_design4

        strain_design3 += strain_design1

        assert t1 in strain_design3
        assert t2 in strain_design3
        assert t3 in strain_design3
        assert t4 not in strain_design3
        assert t5 in strain_design3

    @pytest.mark.skipif(six.PY2, reason="Gnomic is not compatible with python 2")
    def test_design_to_gnomic(self, cad_reaction):
        from gnomic import Genotype
        t1 = ReactionKnockoutTarget('PGI')
        t2 = ReactionKnockoutTarget('GAPD')
        t3 = ReactionKnockinTarget("CAD", cad_reaction)

        strain_design1 = StrainDesign([t1, t2, t3])

        sd_gnomic = strain_design1.to_gnomic()

        assert isinstance(sd_gnomic, Genotype)
        assert len(sd_gnomic.added_features) == 1
        assert len(sd_gnomic.removed_features) == 2
