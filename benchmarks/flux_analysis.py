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


#############################
# Flux Variability Analysis #
#############################
from vbench.benchmark import Benchmark
from common import *

fva_setup = """
from cameo.flux_analysis.analysis import flux_variability_analysis
from cameo import config
from cameo.parallel import SequentialView

config.default_view = SequentialView()
"""

fva_statement = """
flux_variability_analysis(model)
"""

#GLPK
glpk_fva_benchmark = Benchmark(fva_statement, common_setup + fva_setup + glpk_model_setup)

#CPLEX
cplex_fva_benchmark = Benchmark(fva_statement, common_setup + fva_setup + cplex_model_setup)

######################
# Simulation methods #
######################

simulation_setup = """
from cameo.flux_analysis.simulation import *
"""

#FBA
fba_statement = """
fba(model)
"""

#GLPK
glpk_fba_benchmark = Benchmark(fba_statement, common_setup + simulation_setup + glpk_model_setup)

#CPLEX
cplex_fba_benchmark = Benchmark(fba_statement, common_setup + simulation_setup + cplex_model_setup)