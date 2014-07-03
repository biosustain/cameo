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

from cameo.strain_design.heuristics import ReactionKnockoutOptimization
from cameo.strain_design.heuristics.multiprocess import MultiprocessHeuristicOptimization
from cameo.strain_design.heuristics.objective_functions import bpcy
from cameo.flux_analysis.simulation import fba
from cameo.solver_based_model import to_solver_based_model
from cobra.io import read_sbml_model
from optlang import glpk_interface
import inspyred

model = read_sbml_model("/Users/joao/Documents/repos/cameo/tests/data/iJO1366.xml")
model = to_solver_based_model(model, solver_interface=glpk_interface)

of = bpcy("Ec_biomass_iJO1366_core_53p95M", "EX_ac_LPAREN_e_RPAREN_", "EX_glc_LPAREN_e_RPAREN_")

ko1 = ReactionKnockoutOptimization(model=model, objective_function=of,
                                   simulation_method=fba, heuristic_method=inspyred.ec.GA)
ko2 = ReactionKnockoutOptimization(model=model, objective_function=of,
                                   simulation_method=fba, heuristic_method=inspyred.ec.GA)
ko3 = ReactionKnockoutOptimization(model=model, objective_function=of,
                                   simulation_method=fba, heuristic_method=inspyred.ec.GA)
ko4 = ReactionKnockoutOptimization(model=model, objective_function=of,
                                   simulation_method=fba, heuristic_method=inspyred.ec.GA)

mp = MultiprocessHeuristicOptimization(islands=[ko1, ko2, ko3, ko4])
mp.run(max_evaluations=30000, n=2)