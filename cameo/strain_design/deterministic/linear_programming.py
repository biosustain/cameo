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
from cobra import Metabolite
from cameo import Reaction
from cameo.core.solver_based_model_dual import SolverBasedModel, to_dual_model
from cameo.strain_design import StrainDesignMethod


class OptKnock(StrainDesignMethod):
    def __init__(self, model, essential_reactions=None, biomass_reaction=None,
                 biomass_exchange=None, *args, **kwargs):
        super(OptKnock, self).__init__(*args, **kwargs)
        assert isinstance(model, SolverBasedModel)
        self._model = to_dual_model(model.copy())
        if biomass_exchange is None:
            self.biomass_exchange = biomass_exchange
            self.biomass_metabolite = biomass_exchange.metabolites.keys()[0]
        else:
            self.biomass_exchange, self.biomass_metabolite = self._add_biomass_metabolite(biomass_reaction)
        self._build_problem(essential_reactions)

    def _add_biomass_metabolite(self, biomass_reaction):
        assert isinstance(biomass_reaction, Reaction)
        biomass = Metabolite("biomass_c", name="Biomass")
        biomass_reaction.add_metabolites({biomass: 1})
        ex_biomass = Reaction("EX_biomass_c")
        ex_biomass.name = "Exchange Biomass"
        ex_biomass.add_metabolites({biomass: -1})
        ex_biomass.upper_bound = 1000
        self._model.add_reactions([biomass_reaction])
        return biomass_reaction, biomass

    def _build_problem(self, essential_reactions):
        self.essential_reactions = self._model.essential_reactions()
        self.essential_reactions += essential_reactions or []

    def run(self, carbon_source, target, *args, **kwargs):
        pass


class RobustKnock(StrainDesignMethod):
    pass