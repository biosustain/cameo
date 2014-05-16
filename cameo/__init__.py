# Copyright 2013 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


__version__ = 'v0.0.0'

from .io import load_model
import config

# from cameo.io import load_model

# from .io import load_model
# from .solver_based_model import Reaction
# from .solver_based_model import OptlangBasedModel as Model


# import os
# from cameo.io import load_model
# debug_model = load_model(os.path.join(os.path.dirname(__file__), '../tests/data/EcoliCore.xml'))