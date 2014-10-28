# Copyright 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


#############################
# Flux Variability Analysis #
#############################
from datetime import datetime
from vbench.benchmark import Benchmark
from common import *

fva_setup = """
from cameo.flux_analysis.analysis import flux_variability_analysis
from cameo import config
from cameo.parallel import SequentialView, MultiprocessingView
"""

sequential_fva_statement = """
flux_variability_analysis(model, view=SequentialView())
"""

#GLPK
glpk_sequential_fva_benchmark = Benchmark(sequential_fva_statement,
                               common_setup + fva_setup + glpk_model_setup,
                               start_date=datetime(2014, 7, 15))

#CPLEX
cplex_sequential_fva_benchmark = Benchmark(sequential_fva_statement,
                                common_setup + fva_setup + cplex_model_setup,
                                start_date=datetime(2014, 7, 15))


multiprocessing_fva_statement = """
flux_variability_analysis(model, view=MultiprocessingView())
"""

#GLPK
glpk_multiprocessing_fva_benchmark = Benchmark(multiprocessing_fva_statement,
                               common_setup + fva_setup + glpk_model_setup,
                               start_date=datetime(2014, 7, 15))

#CPLEX
cplex_multiprocessing_fva_benchmark = Benchmark(multiprocessing_fva_statement,
                                common_setup + fva_setup + cplex_model_setup,
                                start_date=datetime(2014, 7, 15))


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
glpk_fba_benchmark = Benchmark(fba_statement,
                               common_setup + simulation_setup + glpk_model_setup,
                               start_date=datetime(2014, 7, 15))

#CPLEX
cplex_fba_benchmark = Benchmark(fba_statement,
                                common_setup + simulation_setup + cplex_model_setup,
                                start_date=datetime(2014, 7, 15))