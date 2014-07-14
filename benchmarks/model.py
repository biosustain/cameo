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

import os
from vbench.benchmark import Benchmark
from common import *

BENCHMARKS_DIR = os.path.dirname(__file__)

MODEL_DIR = os.path.join(BENCHMARKS_DIR, "../examples/data/iJO1366.xml")

#########################
# Load model benchmarks #
#########################


#GLPK
glpk_load_statemment = """
model = to_solver_based_model(cobra_model, solve='glpk')
"""

glpk_load_benchmark = Benchmark(glpk_load_statemment, common_setup + read_sbml_model)


#CPLEX
cplex_load_statemment = """
model = to_solver_based_model(cobra_model, solve='cplex')
"""

cplex_load_benchmark = Benchmark(cplex_load_statemment, common_setup + read_sbml_model)


##################
# Solve model #
##################

simulate_statement = """
model.solve()
"""
#GLPK
glpk_simulate_benchmark = Benchmark(simulate_statement, common_setup + glpk_model_setup)

#CPLEX
cplex_simulate_benchmark = Benchmark(simulate_statement, common_setup + cplex_model_setup)

######################
# Critical reactions #
######################

critical_reactions_statement = """
model.critical_reactions()
"""

#GLPK
glpk_critical_reactions_benchmark = Benchmark(critical_reactions_statement, common_setup + glpk_model_setup)

#CPLEX
cplex_critical_reactions_benchmark = Benchmark(critical_reactions_statement, common_setup + cplex_model_setup)

##################
# Critical genes #
##################

critical_genes_statement = """
model.critical_genes()
"""

#GLPK
glpk_critical_genes_benchmark = Benchmark(critical_genes_statement, common_setup + glpk_model_setup)

#CPLEX
cplex_critical_reactions_benchmark = Benchmark(critical_genes_statement, common_setup + cplex_model_setup)
