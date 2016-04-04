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

import time

from bisect import insort

from inspyred.ec import Individual as OriginalIndividual

__all__ = ['BestSolutionArchive', 'ProductionStrainArchive']


class BestSolutionArchive(object):
    def __init__(self):
        self.__name__ = self.__class__.__name__
        self.worst_fitness = None
        self.archive = []

    def __call__(self, random, population, archive, args):
        self.archive = archive
        maximize = args.get("maximize", True)
        max_size = args.get('max_archive_size', 100)
        [self.add(i.candidate, i.fitness, i.birthdate, maximize, max_size) for i in population]
        return self.archive

    def add(self, candidate, fitness, birthdate, maximize, max_size):
        if self.worst_fitness is None:
            self.worst_fitness = fitness

        if (maximize and fitness >= self.worst_fitness) or (not maximize and fitness <= self.worst_fitness):

            candidate = Individual(candidate, fitness, maximize, birthdate)
            add = True
            for c in self.archive:
                if c == candidate:
                    add = False
                elif c.improves(candidate) and candidate.fitness == c.fitness:
                    add = False
                elif candidate.improves(c) and candidate.fitness == c.fitness:
                    self.archive.remove(c)

            if add:
                insort(self.archive, candidate)

            while self.length() > max_size:
                self.archive.pop()

            self.worst_fitness = self.archive[len(self.archive) - 1].fitness

    def length(self):
        return len(self.archive)

    def get(self, index):
        return self.archive[index]

    def __iter__(self):
        for solution in self.archive:
            yield solution

    def __len__(self):
        return self.length()


class ProductionStrainArchive(BestSolutionArchive):
    def __call__(self, random, population, archive, args):
        self.archive = archive
        max_size = args.get('max_archive_size', 100)
        [self.add(i.candidate, i.fitness, i.birthdate, True, max_size) for i in population if i.fitness > 0]
        return self.archive


class Individual(OriginalIndividual):
    def __init__(self, candidate, fitness, maximize=True, birthdate=time.time()):
        super(Individual, self).__init__(set(candidate), maximize)
        self.fitness = fitness
        self.birthdate = birthdate

    def __iter__(self):
        for v in self.candidate:
            yield v

    def __eq__(self, other):
        return self.candidate == other.candidate and self.fitness == other.fitness

    def __cmp__(self, other):
        if self.fitness > other.fitness:
            return -1 if self.maximize else 1
        elif self.fitness == other.fitness:
            if self.improves(other):
                return -1
            elif self == other:
                return 0
            else:
                return 1
        else:
            return 1 if self.maximize else -1

    def __lt__(self, other):
        if self.fitness > other.fitness:
            return self.maximize
        elif self.fitness == other.fitness:
            if self.improves(other):
                return True
            elif self == other:
                return False
            else:
                return False
        else:
            return not self.maximize

    def __ge__(self, other):
        return self.__eq__(other) or self.__gt__(other)

    def __le__(self, other):
        return self.__eq__(other) or self.__lt__(other)

    def __gt__(self, other):
        if self.fitness > other.fitness:
            return not self.maximize
        elif self.fitness == other.fitness:
            if self.improves(other):
                return False
            elif self == other:
                return False
            else:
                return True
        else:
            return self.maximize

    def __str__(self):
        sense = "max" if self.maximize else "min"
        return "%s - %s sense: %s" % (list(self.candidate), self.fitness, sense)

    def __repr__(self):
        return "SolutionTuple #%s: %s" % (id(self), self.__str__())

    def issubset(self, other):
        return self.candidate.issubset(other.candidate)

    def symmetric_difference(self, other):
        return self.candidate.symmetric_difference(other.candidate)

    def improves(self, other):
        assert isinstance(other, Individual)
        if self.maximize:
            return self.issubset(other) and len(self.symmetric_difference(other)) > 0 and self.fitness >= other.fitness
        else:
            return self.issubset(other) and len(self.symmetric_difference(other)) > 0 and self.fitness <= other.fitness
