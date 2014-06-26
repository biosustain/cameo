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

from uuid import uuid1
from bokeh.plotting import *
import pandas


class PlotObserver(object):
    def __init__(self, url='default'):
        self.i = 0
        self.url = url
        self.iterations = []
        self.fitness = []
        self.uuid = None
        self.in_ipnb = pandas.core.common.in_ipnb()
        self.plotted = False

    def _set_plot(self):
        if self.in_ipnb:
            self.uuid = uuid1()
            output_notebook(url=self.url, docname=str(self.uuid))
            scatter([], [], tools='', title="Best solution convergence plot")
            xaxis().axis_label = "Iteration"
            yaxis().axis_label = "Fitness"
            self.plot = curplot()
            renderer = [r for r in self.plot.renderers if isinstance(r, Glyph)][0]
            self.ds = renderer.data_source
            show()
            self.plotted = True

    def __call__(self, population, num_generations, num_evaluations, args):
        if not self.plotted:
            self._set_plot()

        self.i += 1
        self.iterations.append(self.i)
        if len(population) > 0:
            best = max(population)
            self.fitness.append(best.fitness)
        else:
            self.fitness.append(None)

        if self.i % args.get('n', 20) == 0:
            self._ipnb_plot()

    def __name__(self):
        return "Fitness Plot"

    def _ipnb_plot(self):
        if self.in_ipnb:
            self.ds.data['x'] = self.iterations
            self.ds.data['y'] = self.fitness

            session().store_obj(self.ds)

    def reset(self):
        self.i = 0
        self.iterations = []
        self.fitness = []
        self.plotted = False


class ParetoPlotObserver(object):
    def __init__(self, ofs=None, x=0, y=1, url='default'):
        self.url = url
        self.x = x
        self.y = y
        self.ofs = ofs
        self.fitness = []
        self.uuid = None
        self.in_ipnb = pandas.core.common.in_ipnb()
        self.plotted = False

    def _set_plot(self):
        if self.in_ipnb:
            self.uuid = uuid1()
            output_notebook(url=self.url, docname=str(self.uuid))
            scatter([], [], tools='', title="Multi-objective Pareto Fitness Plot")
            xaxis().axis_label = self.ofs[self.x].__name__
            yaxis().axis_label = self.ofs[self.y].__name__
            self.plot = curplot()
            renderer = [r for r in self.plot.renderers if isinstance(r, Glyph)][0]
            self.ds = renderer.data_source
            show()
            self.plotted = True

    def __call__(self, population, num_generations, num_evaluations, args):
        if not self.plotted:
            self._set_plot()

        self.fitness = [i.fitness for i in population]
        if num_evaluations % args.get('n', 20) == 0:
            self._ipnb_plot()

    def __name__(self):
        return "Pareto Plot"

    def _ipnb_plot(self):
        if self.in_ipnb:
            self.ds.data['x'] = [e[self.x] for e in self.fitness]
            self.ds.data['y'] = [e[self.y] for e in self.fitness]
            session().store_obj(self.ds)

    def reset(self):
        self.fitness = []
        self.plotted = False