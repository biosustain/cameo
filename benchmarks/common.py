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

import os

BENCHMARKS_DIR = os.path.dirname(__file__)

MODEL_DIR = os.path.join(BENCHMARKS_DIR, "../tests/data/iJO1366.xml")

common_setup = "MODEL_DIR = '%s'" % MODEL_DIR
common_setup = common_setup + """
from cameo.solver_based_model import to_solver_based_model
from cobra.io import read_sbml_model
""".format(BENCHMARKS_DIR)

#cobrapy model
read_sbml_model = """
cobra_model = read_sbml_model(MODEL_DIR)
"""

#GLPK model
glpk_model_setup = read_sbml_model + """
model = to_solver_based_model(cobra_model, solver_interface='glpk')
"""

#CPLEX model
cplex_model_setup = read_sbml_model + """
model = to_solver_based_model(cobra_model, solver_interface='cplex')
"""

