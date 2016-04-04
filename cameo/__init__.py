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


"""
CAMEO: Computer Aided Metabolic Engineering & Optimization

Cameo is a high-level python library developed to aid the in silico
strain design process in metabolic engineering projects. The library
provides a modular architecture that enables the efficient construction
of custom analysis workflows.

Example
-------

from cameo import load_model

# load a model from SBML format (can be found under cameo/tests/data)
model = load_model('EcoliCore.xml')

# solve the model and print the objective value
solution = model.solve()
print 'Objective value:', solution.f

# Determine a set of gene deletions that will optimize the production
# of a desired compound
from cameo.strain_design.heuristic import GeneKnockoutOptimization
from cameo.strain_design.heuristic.objective_functions import biomass_product_coupled_yield
from cameo.flux_analysis.simulation import fba

objective = biomass_product_coupled_yield("Ec_biomass_iJO1366_core_53p95M",
                                          "EX_succ_lp_e_rp_", "EX_glc_lp_e_rp_")
optimization = GeneKnockoutOptimization(model=model, objective_function=of,
                              simulation_method=fba, heuristic_method=inspyred.ec.GA)
optimization.run(max_evaluations=2000, n=1,
       mutation_rate=0.3, view=cameo.parallel.SequentialView(),
       product="EX_succ_lp_e_rp_", num_elites=1)
"""
from __future__ import absolute_import, print_function

import os
import sys
from cameo import config
from cameo.util import get_system_info, in_ipnb

if sys.version_info[0] == 2:
    import imp


    def find_module(name):
        try:
            imp.find_module(name)
            return True
        except ImportError:
            return False

elif sys.version_info[0] == 3:
    if sys.version_info[1] <= 3:
        from importlib import find_loader as find
    else:
        from importlib.util import find_spec as find


    def find_module(name):
        return find(name) is not None

_cameo_path = __path__[0]
_cameo_data_path = os.path.join(_cameo_path, 'data')

# fix - if matplotlib is installed it is not possible to import cameo without importing matplotlib on jupyter notebook.
if find_module("matplotlib") and in_ipnb():
    from IPython import get_ipython

    ipython = get_ipython()
    ipython.magic("matplotlib inline")

system_info = get_system_info()

from ._version import get_versions

__version__ = get_versions()['version']
del get_versions

from cameo.io import load_model
#except ImportError:
#    raise
#    pass

from cameo import models

from cameo.core.solver_based_model import SolverBasedModel as Model
from cameo.core.reaction import Reaction
from cobra.core import Metabolite

from .flux_analysis.analysis import flux_variability_analysis, phenotypic_phase_plane
from .flux_analysis.simulation import fba, pfba
