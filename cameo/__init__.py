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


"""
CAMEO: Computer Assisted Metabolic Engineering & Optimization

from cameo import load_model

model = load_model('../tests/data/EcoliCore.xml')
solution = model.solve()
print 'Objective value:', solution.f
"""

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

try:
    from .io import load_model
except ImportError:
    pass

from .flux_analysis.analysis import flux_variability_analysis, phenotypic_phase_plane
from .flux_analysis.simulation import fba, pfba

import config