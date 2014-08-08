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
from bokeh.objects import Range1d, Plot
from inspyred.ec.emo import Pareto
import numpy as np
from cameo import config

if config.use_bokeh:
    from bokeh.plotting import *


class GenericStatsData(object):
    def __init__(self, solution, *args, **kwargs):
        super(GenericStatsData, self).__init__(*args, **kwargs)
        self.solution = solution
        self.knockouts_hist, self.knockouts_edges = np.histogram(solution.solutions['Size'])

    def display(self):
        raise NotImplementedError


class CLIStatsData(GenericStatsData):
    def __init__(self, *args, **kwargs):
        super(CLIStatsData, self).__init__(*args, **kwargs)

    def display(self):
        pass


class BokehStatsData(GenericStatsData):
    def __init__(self, *args, **kwargs):
        super(BokehStatsData, self).__init__(*args, **kwargs)
        self.xdr = Range1d(start=-0.5, end=20.5)
        self.ydr = Range1d(start=-0.5, end=20.5)

    def display(self):
        output_notebook()
        hold()
        figure(title="Knockout size distribution")
        quad(top=self.knockouts_hist, bottom=np.zeros(len(self.knockouts_hist)),
             left=self.knockouts_edges[:-1], right=self.knockouts_edges[1:],
             title="Knockout size distribution")
        xaxis()[0].axis_label = "Number of knockouts"
        yaxis()[0].axis_label = "Number of solutions"

        figure(title="Correlation between number of knockouts and fitness")
        fitness = self.solution.solutions['Fitness']
        if isinstance(fitness[0], Pareto):
            pass
        else:
            scatter(self.solution.solutions['Size'], self.solution.solutions['Fitness'])
            xaxis()[0].axis_label = "Number of knockouts"
            yaxis()[0].axis_label = "Fitness"
        show()