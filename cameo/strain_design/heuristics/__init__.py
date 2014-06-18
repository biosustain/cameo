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
from bisect import bisect


class BestSolutionPool(object):

    class SolutionTuple(object):
        def __init__(self, solution, fitness):
            self.solution = set(solution)
            self.fitness = fitness

        def __cmp__(self, other):
            if self.fitness > other.fitness:
                return 1
            elif self.fitness == other.fitness:
                return 0
            else:
                return -1

        def issubset(self, other):
            return self.solution.issubset(other.solution)

        def difference(self, other):
            return self.solution.difference(other.solution)

    def __init__(self, size):
        self.solutions = bisect([])
        self.size = size

    def add(self, solution, fitness):
        if fitness > self.worst_fitness:

            solution = self.SolutionTuple(solution, fitness)
            self.solutions.insort(solution)
            worst_fitness = fitness
            worst_solution = None

            for sol, fit in self.solutions:

                if solution.issubset(sol) and len(solution.difference(fit)) == 0 and fit <= fitness:
                    del self.solutions[sol]

                if fit < worst_fitness:
                    worst_fitness = fit
                    worst_solution = sol

            if len(self.solutions) > self.size:
                del self.solutions[worst_solution]
