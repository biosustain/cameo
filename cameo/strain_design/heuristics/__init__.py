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
from bisect import bisect, insort
from uuid import uuid1
from bokeh.plotting import *
import pandas.core.common


class PlotObserver(object):
    def __init__(self, url='default'):
        self.i = 0
        self.iterations = []
        self.fitness = []
        self.uuid = uuid1()
        self.in_ipnb = pandas.core.common.in_ipnb()
        if self.in_ipnb:
            output_notebook(url=url, docname=str(self.uuid))
            scatter([], [], tools='', title="Convergence")
            self.plot = curplot()
            renderer = [r for r in self.plot.renderers if isinstance(r, Glyph)][0]
            self.ds = renderer.data_source
            show()

    def __call__(self, population, num_generations, num_evaluations, args):
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
        self._ipnb_plot()



class BestSolutionPool(object):
    class SolutionTuple(object):
        def __init__(self, solution, fitness):
            self.solution = set(solution)
            self.fitness = fitness

        def __eq__(self, other):
            return self.solution == other.solution and self.fitness == other.fitness

        def __cmp__(self, other):
            if self.fitness > other.fitness:
                return -1
            elif self.fitness == other.fitness:
                if self.improves(other):
                    return -1
                elif self == other:
                    return 0
                else:
                    return 1
            else:
                return 1

        def __str__(self):
            return "%s - %s" % (list(self.solution), self.fitness)

        def __repr__(self):
            return "SolutionTuple #%s: %s" % (id(self), self.__str__())

        def issubset(self, other):
            return self.solution.issubset(other.solution)

        def symmetric_difference(self, other):
            return self.solution.symmetric_difference(other.solution)

        def improves(self, other):
            return self.issubset(other) and len(self.symmetric_difference(other)) > 0 and self.fitness >= other.fitness

    def __init__(self, size):
        self.pool = []
        bisect(self.pool, 0)
        self.size = size
        self.worst_fitness = 0

    def add(self, solution, fitness):
        if fitness >= self.worst_fitness:

            solution = self.SolutionTuple(solution, fitness)

            for sol in self.pool:
                if solution.improves(sol) or solution == sol:
                    self.pool.remove(sol)

            insort(self.pool, solution)

            while len(self.pool) > self.size:
                self.pool.pop()

            self.worst_fitness = self.pool[len(self.pool) - 1].fitness

    def length(self):
        return len(self.pool)

    def get(self, index):
        return self.pool[index]

    def __iter__(self):
        for solution in self.pool:
            yield solution
