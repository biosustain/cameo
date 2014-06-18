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

from cameo.flux_analysis.simulation import pfba

class HeuristicOptimization(object):
    def __init__(self, model=model, heuristic_method=heuristic_method, objective_function=objective_function, *args,
                 **kwargs):
        super(HeuristicOptimization, self).__init__(*args, **kwargs)
        self.model = model
        self.heuristic_method = heuristic_method
        self.objective_function = objective_function


class KnockoutOptimization(HeuristicOptimization):
    def __init__(self, simulation_method=pfba, *args, **kwargs):
        super(HeuristicOptimization, self).__init__(*args, **kwargs)
        #check reactions


class GeneKnockoutOptimization(KnockoutOptimization):
    pass


