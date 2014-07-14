# Copyright 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import inspyred

from cameo.strain_design.heuristic.algorithms import ReactionKnockoutOptimization
from cameo import load_model
from cameo.strain_design.heuristic.objective_functions import bpcy
from cameo.flux_analysis.simulation import fba


#cobrapy_model =  read_sbml_model("/Users/joao/Downloads/iJO1366.xml")
model = load_model("../tests/data/iJO1366.xml")
of = bpcy("Ec_biomass_iJO1366_core_53p95M", "EX_succ_LPAREN_e_RPAREN_", "EX_glc_LPAREN_e_RPAREN_")

ko = ReactionKnockoutOptimization(model=model, objective_function=of, simulation_method=fba)
ko.heuristic_method = inspyred.ec.SA
ko.run(max_generations=300, max_evaluations=30000)