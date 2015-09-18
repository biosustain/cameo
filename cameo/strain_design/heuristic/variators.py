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

__all__ = ['set_mutation', 'set_indel']

from six.moves import range

from inspyred.ec.variators import mutator, crossover
from ordered_set import OrderedSet
from cameo.strain_design.heuristic.genomes import MultipleChromosomeGenome
from numpy import float32 as float


def _do_set_n_point_crossover(representation, mom, dad, points, random, candidate_size):
    chunks = []
    i = 0
    for point in points:
        chunks.append(representation[i:point])
        i = point
    chunks.append(representation[i:])

    bro = OrderedSet()
    sis = OrderedSet()

    cross = True
    for variables in chunks:
        for v in variables:
            if v in mom:
                bro.append(v) if cross else sis.append(v)
            if v in dad:
                sis.append(v) if cross else bro.append(v)
        cross = not cross

    if len(bro) > candidate_size:
        bro = random.sample(bro, candidate_size)

    if len(sis) > candidate_size:
        sis = random.sample(sis, candidate_size)
    return bro, sis


@crossover
def set_n_point_crossover(random, mom, dad, args):
    representation = list(set(mom).union(set(dad)))
    crossover_rate = args.setdefault('crossover_rate', 1.0)
    num_crossover_points = args.setdefault('num_crossover_points', 1)
    candidate_size = args.setdefault('candidate_size', 9)
    children = []
    if random.random() <= crossover_rate:
        points = random.sample(representation, num_crossover_points)
        bro, sis = _do_set_n_point_crossover(representation, mom, dad, points, random, candidate_size)

        # ensure number of knockouts > 0 or do not add individual
        if len(bro) > 0:
            children.append(bro)
        if len(sis) > 0:
            children.append(sis)
    else:
        children.append(mom)
        children.append(dad)
    return children


@mutator
def set_mutation(random, individual, args):
    """
    Mutates a given set based on the entries available on the representation.

    Parameters
    ----------

    random: Random
    individual: list
        with unique integers
    args: dict
        must contain the representation

    Returns
    -------
    list
        created based on an ordered set

    """
    representation = args.get('representation')
    new_individual = []
    mutation_rate = float(args.get('mutation_rate', .1))
    for index in individual:
        if random.random() < mutation_rate:
            new_individual.append(random.randint(0, len(representation) - 1))
        else:
            new_individual.append(index)

    return list(OrderedSet(new_individual))


@mutator
def set_indel(random, individual, args):
    """
    Creates a random insertion or deletion in the individual.

    Parameters
    ----------

    random: Random
    individual: list
        with unique integers
    args: dict
        must contain the representation

    Returns
    -------
    list
        created based on an ordered set

    """
    representation = args.get('representation')
    indel_rate = float(args.get('indel_rate', .1))
    new_individual = list(individual)
    if random.random() < indel_rate:
        if random.random() > 0.5:
            if len(individual) > 1:
                new_individual = random.sample(new_individual, len(new_individual) - 1)
        else:
            new_individual.append(random.sample(range(len(representation)), 1)[0])

    return list(OrderedSet(new_individual))


@mutator
def multiple_chromosome_set_mutation(random, individual, args):
    """
    Mutates a given set based on the entries available on the representation.

    Parameters
    ----------

    random: Random
    individual: MultipleChromosomeGenome
        with unique integers in each chromosome
    args: dict
        must contain the representation of each chromosome

    Returns
    -------
    MultipleChromosomeGenome
        A mutated individual

    """
    new_individual = MultipleChromosomeGenome(keys=individual.keys)

    for key in individual.keys:
        representation = args.get('%s_representation' % key)
        mutation_rate = args.get('%s_mutation_rate' % key, .1)
        for index in individual[key]:
            if random.random() < mutation_rate:
                new_individual[key].append(random.randint(0, len(representation) - 1))
            else:
                new_individual[key].append(index)

    return new_individual


@mutator
def multiple_chromosome_set_indel(random, individual, args):
    """
    Creates a random insertion or deletion in the individual.

    Parameters
    ----------

    random: Random
    individual: MultipleChromosomeGenome
        with unique integers in each chromosome
    args: dict
        must contain the representation of each chromosome

    Returns
    -------
    MultipleChromosomeGenome
        A mutated individual
    """
    new_individual = individual.copy()

    for key in individual.keys:
        representation = args.get('%s_representation' % key)
        indel_rate = args.get('%s_indel_rate' % key, .1)
        if random.random() < indel_rate:
            if random.random() > 0.5:
                if len(individual[key]) > 1:
                    new_individual[key] = random.sample(new_individual[key], len(new_individual[key]) - 1)
            else:
                new_individual[key].append(random.sample(range(len(representation)), 1)[0])

    return new_individual


@crossover
def multiple_chromosome_n_point_crossover(random, mom, dad, args):
    children = MultipleChromosomeGenome(keys=mom.keys)
    for key in children.keys:
        children[key] = set_n_point_crossover(random, [mom[key], dad[key]], args)

    return children
