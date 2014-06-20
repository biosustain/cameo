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
import matplotlib.pyplot as plt

try:
    from IPython.display import display, clear_output
    USE_IPYTHON = True
except:
    USE_IPYTHON = False


class PlotObserver(object):
    def __init__(self):
        self.i = 0
        self.iterations = []
        self.fitness = []
        self.f, self.ax = plt.subplots()

    def __call__(self, population, num_generations, num_evaluations, args):
        best = max(population)
        self.i += 1
        self.iterations.append(self.i)
        self.fitness.append(best.fitness)
        if self.i % 20 == 0:
            self.ax.plot(self.iterations, self.fitness, 'ro')
            self.ax.axis([0, self.i+1, 0, max(self.fitness)+0.5])
            if USE_IPYTHON:
                clear_output()
                display(self.f)

    def __name__(self):
        return "Fitness Plot"


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


