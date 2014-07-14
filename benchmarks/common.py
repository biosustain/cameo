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

common_setup = """
from cameo.solver_base_model import to_solver_based_model
from cobra.io import read_sbml_model
"""

#cobrapy model
read_sbml_model = """
cobra_model = read_sbml_model(MODEL_DIR)
"""

#GLPK model
glpk_model_setup = read_sbml_model + """
model = to_solver_based_model(cobra_model, solve='glpk')
"""

#CPLEX model
cplex_model_setup = read_sbml_model + """
model = to_solver_based_model(cobra_model, solve='cplex')
"""

