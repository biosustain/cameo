# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import absolute_import, print_function

import scipy
import numpy as np

from cameo import config

if config.use_bokeh:
    from bokeh.plotting import show, figure
    from bokeh.models import ColumnDataSource
    from bokeh.io import push_notebook


class IPythonBokehFitnessPlotter(object):
    __name__ = "IPython Bokeh Fitness Plot"

    def __init__(self, window_size=1000):
        self.iteration = 0
        self.window_size = window_size
        self.iterations = []
        self.fitness = []
        self.uuid = None
        self.plotted = False

    def _set_plot(self):
        self.plot = figure(title="Fitness plot", tools='', plot_height=400, plot_width=650)
        self.plot.xaxis.axis_label = "Iteration"
        self.plot.yaxis.axis_label = "Fitness"
        self.ds = ColumnDataSource(data=dict(x=[], y=[]))
        self.plot.circle('x', 'y', source=self.ds)
        show(self.plot)
        self.plotted = True

    def __call__(self, population, num_generations, num_evaluations, args):
        if not self.plotted:
            self._set_plot()

        self.iteration += 1
        self.iterations.append(self.iteration)
        if len(population) > 0:
            best = max(population)
            self.fitness.append(best.fitness)
        else:
            self.fitness.append(None)

        if self.iteration % args.get('n', 1) == 0:
            self._update_plot()

    def _update_plot(self):
        self.ds.data['x'] = self.iterations[-self.window_size:]
        self.ds.data['y'] = self.fitness[-self.window_size:]
        push_notebook()

    def reset(self):
        self.iteration = 0
        self.iterations = []
        self.fitness = []
        self.plotted = False

    def end(self):
        # end is not necessary
        pass


class IPythonBokehParetoPlotter(object):
    __name__ = "IPython Bokeh Pareto Plotter"

    def __init__(self, ofs=None, x=0, y=1):
        self.x = x
        self.y = y
        self.ofs = ofs
        self.fitness = []
        self.uuid = None
        self.plotted = False

    def _set_plot(self):
        self.plot = figure(title="Multi-objective Fitness Plot", tools='', plot_height=400, plot_width=650)
        self.plot.xaxis.axis_label = self.ofs[self.x].name
        self.plot.yaxis.axis_label = self.ofs[self.y].name
        self.ds = ColumnDataSource(data=dict(x=[], y=[]))
        self.plot.circle('x', 'y', source=self.ds)

        show(self.plot)
        self.plotted = True

    def __call__(self, population, num_generations, num_evaluations, args):
        if not self.plotted:
            self._set_plot()

        self.fitness = [i.fitness for i in population]
        if num_evaluations % args.get('n', 1) == 0:
            self._update()

    def _update(self):
        self.ds.data['x'] = [e[self.x] for e in self.fitness]
        self.ds.data['y'] = [e[self.y] for e in self.fitness]
        self.ds.push_notebook()

    def reset(self):
        self.fitness = []
        self.plotted = False

    def end(self):
        # does not need end
        pass


class GeneFrequencyPlotter(object):
    def __init__(self, solutions, url='default'):
        self.solutions = solutions
        self.url = url
        self.freqs = self.frequencies()

    def frequencies(self):
        kos = []
        for solution in self.solutions:
            for ko in solution.knockout_list:
                kos.append(ko.id)

        return scipy.stats.itemfreq(kos)

    def plot(self):
        plot = figure()

        plot.quad(top=self.freqs[:, 1], left=self.freqs[:, 1], bottom=np.zeros(len(self.freqs[:, 1])),
                  right=self.freqs[:, 1], x_range=list(self.freqs[:, 0]))
        plot.xaxis().major_label_orientation = np.pi / 3
        show(plot)
