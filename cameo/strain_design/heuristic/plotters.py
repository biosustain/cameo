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


from uuid import uuid1

import scipy
import numpy as np

from cameo import config

if config.use_bokeh:
    from bokeh.plotting import *
    from bokeh.models import GlyphRenderer

class IPythonBokehFitnessPlotter(object):
    __name__ = "IPython Bokeh Fitness Plot"

    def __init__(self, window_size=1000):
        self.iteration = 0
        self.window_size = window_size
        self.iterations = []
        self.fitness = []
        self.uuid = None
        self.plotted = False
        self.can_plot = True

    def _set_plot(self):
        self.uuid = uuid1()
        try:
            output_notebook(url=config.bokeh_url, docname=str(self.uuid))
            p = figure(title="Best solution fitness plot", tools='')
            p.scatter([], [])
            p.xaxis.axis_label = "Iteration"
            p.yaxis.axis_label = "Fitness"
            self.plot = p
            renderer = p.select(dict(type=GlyphRenderer))
            self.ds = renderer[0].data_source
            show(p)
            self.plotted = True
        except SystemExit as e:
            logger.info("Bokeh-server is not running. Skipping plotting.")
            logger.error(e)
            self.can_plot = False

    def __call__(self, population, num_generations, num_evaluations, args):
        if not self.plotted and self.can_plot:
            self._set_plot()

        self.iteration += 1
        self.iterations.append(self.iteration)
        if len(population) > 0:
            best = max(population)
            self.fitness.append(best.fitness)
        else:
            self.fitness.append(None)

        if self.iteration % args.get('n', 20) == 0:
            self._update_plot()

    def _update_plot(self):
        if self.can_plot:
            self.ds.data['x'] = self.iterations[-self.window_size:]
            self.ds.data['y'] = self.fitness[-self.window_size:]
            cursession().store_objects(self.ds)

    def reset(self):
        self.iteration = 0
        self.iterations = []
        self.fitness = []
        self.plotted = False

    def end(self):
        pass


class IPythonBokehParetoPlotter(object):
    __name__ = "IPython Bokeh Pareto Plotter"

    def __init__(self, ofs=None, x=0, y=1, url='default'):
        self.url = url
        self.x = x
        self.y = y
        self.ofs = ofs
        self.fitness = []
        self.uuid = None
        self.plotted = False
        self.can_plot = True

    def _set_plot(self):
        try:
            self.uuid = uuid1()
            output_notebook(url=config.bokeh_url, docname=str(self.uuid))
            p = figure(tools='', title="Multi-objective Pareto Fitness Plot")
            p.scatter([], [])
            p.xaxis.axis_label = self.ofs[self.x].name
            p.yaxis.axis_label = self.ofs[self.y].name
            self.plot = p
            renderer = p.select(dict(type=GlyphRenderer))
            self.ds = renderer[0].data_source
            show(p)
            self.plotted = True
        except SystemExit as e:
            logger.info("Bokeh-server is not running. Skipping plotting.")
            logger.error(e)
            self.can_plot = False

    def __call__(self, population, num_generations, num_evaluations, args):
        if not self.plotted and self.can_plot:
            self._set_plot()

        self.fitness = [i.fitness for i in population]
        if num_evaluations % args.get('n', 20) == 0:
            self._update()

    def _update(self):
        if self.can_plot:
            self.ds.data['x'] = [e[self.x] for e in self.fitness]
            self.ds.data['y'] = [e[self.y] for e in self.fitness]
            cursession().store_objects(self.ds)

    def reset(self):
        self.fitness = []
        self.plotted = False

    def end(self):
        pass


class GeneFrequencyPlotter():
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
        self.uuid = uuid1()
        output_notebook(url=config.bokeh_url, docname=str(self.uuid))
        p = figure()

        p.quad(top=self.freqs[:, 1], left=self.freqs[:, 1], bottom=np.zeros(len(self.freqs[:, 1])),
             right=self.freqs[:, 1], x_range=list(self.freqs[:, 0]))
        p.xaxis.major_label_orientation = np.pi / 3
        show(p)
