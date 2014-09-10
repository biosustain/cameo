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

ID_SANITIZE_RULES = [('_DASH_', '__'), ('_FSLASH_', '/'),('_BSLASH_', "\\"), ('_LPAREN_', '('), ('_LSQBKT_', '['),
                     ('_RSQBKT_', ']'), ('_RPAREN_', ')'), ('_COMMA_', ','), ('_PERIOD_', '.'), ('_APOS_', "'"),
                     ('&amp;', '&'), ('&lt;', '<'), ('&gt;', '>'), ('&quot;', '"'), ('__', '-')]


#TODO: describe cameo with some examples
"""
CAMEO: Computer Assisted Metabolic Engineering & Optimization


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


def sanitize_ids(model):
    for metabolite in model.metabolites:
        metabolite.id = _apply_sanitize_rules(metabolite.id)
        metabolite.name = _apply_sanitize_rules(metabolite.name)

    for reaction in model.reactions:
        reaction.id = _apply_sanitize_rules(reaction.id)
        reaction.name = _apply_sanitize_rules(reaction.name)

    for variable in model.solver.variables:
        variable.name = _apply_sanitize_rules(variable.name)

    for constraint in model.solver.constraints:
        constraint.name = _apply_sanitize_rules(constraint.name)


def _apply_sanitize_rules(id):
    for rule in ID_SANITIZE_RULES:
        id = id.replace(*rule)

    return id