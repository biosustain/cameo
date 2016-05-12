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

Implementation of the islands model. The islands model consists of having multiple processes running in parallel where
individual process is running an isolated heuristic optimization, therefor called island.
At certain time points, individuals in an island can migrate through a Queue. That individual will arrive in a another
island and introduce variability in population of that island.

For information on how to run Gene Knockout or Reaction Knockout optimizations refer to
cameo.strain_design.heuristic.optimization

The result object is the same as in the single objective optimization. The knockouts solutions resulting from all
processes are merged.

Examples
--------
>>> from cameo import models
>>> from cameo.strain_design.heuristic.evolutionary.multiprocess import MultiprocessGeneKnockoutOptimization
>>> from cameo.strain_design.heuristic.evolutionary.objective_functions import biomass_product_coupled_yield
>>> import inspyred
>>>
>>> model = models.bigg.iJO1366
>>> objective_function = biomass_product_coupled_yield('BIOMASS_iJO1366_core_53p95M', 'EX_succ_e', 'EX_glc__D_e')
>>> opt = MultiprocessGeneKnockoutOptimization(model=model,
>>>                                            objective_function=objective_function,
>>>                                            heuristic_method=inspyred.ec.GA,
>>>                                            max_migrants=1)
>>> result = opt.run()


"""
from .optimization import MultiprocessGeneKnockoutOptimization, MultiprocessReactionKnockoutOptimization
