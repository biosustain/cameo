# Copyright 2015 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
from time import sleep
from IPython.core.display import display
from cameo.reporters import AbstractReporter
from cameo.visualization.escher_ext import NotebookBuilder


class FluxDistributionReporter(AbstractReporter):
    def close(self):
        raise NotImplementedError

    def __call__(self, model, flux_distribution_result):
        raise NotImplementedError


class EscherFluxReporter(FluxDistributionReporter):
    def __init__(self, pathway):
        self.pathway = pathway
        self.builders = {}

    def __call__(self, model, flux_distribution_result):
        if model.id not in self.builders.keys():
            self.builders[model.id] = NotebookBuilder(model=model,
                                                      map_name=self.pathway,
                                                      reaction_styles=['color'],
                                                      reaction_scale=[{'type': 'min', 'color': '#cc0033'},
                                                                      {'type': 'mean', 'color': '#66ccff'},
                                                                      {'type': 'max', 'color': '#66ff00'}])
            builder = self.builders[model.id]
            display(builder.display_in_notebook(height=650))
        else:
            builder = self.builders[model.id]

        if flux_distribution_result is not None:
            builder.update(reaction_data=flux_distribution_result.fluxes)

    def close(self):
        self.builders = {}