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

from __future__ import absolute_import, print_function

from six.moves import zip

from inspyred.ec.emo import Pareto
import numpy as np
from cameo import config

if config.use_bokeh:
    from bokeh.plotting import output_notebook, figure, show
    from bokeh.models import Range1d


class GenericStatsData(object):
    def __init__(self, solution, *args, **kwargs):
        super(GenericStatsData, self).__init__(*args, **kwargs)
        self.solution = solution
        self.knockouts_hist, self.knockouts_edges = np.histogram(solution.solutions['Size'])

    def display(self):
        raise NotImplementedError


try:
    from bashplotlib.scatterplot import plot_scatter
    from bashplotlib.histogram import plot_hist
except ImportError:
    pass
else:
    class CLIStatsData(GenericStatsData):
        def __init__(self, *args, **kwargs):
            super(CLIStatsData, self).__init__(*args, **kwargs)

        def display(self):
            plot_hist(list(self.knockouts_hist), title="Knockout size distribution", colour="blue")
            lines = ["%s, %s" % (x, y) for x, y in
                     zip(self.solution.solutions['Size'], self.solution.solutions['Fitness'])]
            plot_scatter(lines, None, None, 20, "*", "blue", "Correlation between number of knockouts and fitness")


class BokehStatsData(GenericStatsData):
    def __init__(self, *args, **kwargs):
        super(BokehStatsData, self).__init__(*args, **kwargs)
        self.xdr = Range1d(start=-0.5, end=20.5)
        self.ydr = Range1d(start=-0.5, end=20.5)

    def display(self):
        output_notebook(hide_banner=True)
        plot = figure(title="Knockout size distribution")
        plot.quad(top=self.knockouts_hist, bottom=np.zeros(len(self.knockouts_hist)),
                  left=self.knockouts_edges[:-1], right=self.knockouts_edges[1:],
                  title="Knockout size distribution")
        plot.xaxis.axis_label = "Number of knockouts"
        plot.yaxis.axis_label = "Number of solutions"
        show(plot)

        plot = figure(title="Correlation between number of knockouts and fitness")
        fitness = self.solution.solutions['Fitness']
        if isinstance(fitness[0], Pareto):
            pass
        else:
            plot.scatter(self.solution.solutions['Size'], self.solution.solutions['Fitness'])
            plot.xaxis.axis_label = "Number of knockouts"
            plot.yaxis.axis_label = "Fitness"
        show(plot)
