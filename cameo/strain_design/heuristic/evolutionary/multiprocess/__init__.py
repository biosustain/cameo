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
Multiprocess heuristic optimization.

It implements an islands model. Multiple processes run in parallel and each single process is an heuristic optimization.
All of those processes are connected to a Queue and every iteration an individual from a process migrates through the
Queue to another process enriching the results of the receiving by introducing variability.

For information on how to run Gene Knockout or Reaction Knockout optimizations refer to
~cameo.strain_design.heuristic.optimization

The result is the same as in the single objective. The knockouts solutions resulting from all processes are merged.

Examples
--------
>>> from cameo import models
>>> from cameo.strain_design.heuristic.evolutionary.multiprocess import MultiprocessGeneKnockoutOptimization
>>> from cameo.strain_design.heuristic.evolutionary.objective_functions import biomass_product_coupled_yield
>>> import inspyred
>>> model = models.bigg.iJO1366
>>> objective_function = biomass_product_coupled_yield('Ec_biomass_iJO1366_core_53p95M', 'EX_succ_e', 'EX_glc__D_e')
>>> opt = MultiprocessGeneKnockoutOptimization(model=model, objective_function=objective_function,
>>>                                            heuristic_method=inspyred.ec.GA, max_migrants=1)
>>> result = opt.run()


"""
from .optimization import MultiprocessGeneKnockoutOptimization, MultiprocessReactionKnockoutOptimization
