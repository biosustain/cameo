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
from inspyred.ec.variators import mutator, crossover, n_point_crossover
from ordered_set import OrderedSet
from cameo.strain_design.heuristic.genomes import MultipleChromosomeGenome
from numpy import float32 as float

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
                new_individual.pop(random.randint(0, len(new_individual) - 1))
        else:
            new_individual.append(random.sample(xrange(len(representation)), 1)[0])

    return list(OrderedSet(new_individual))


@mutator
def multiple_chromosome_set_mutation(random, individual, args):
    new_individual = MultipleChromosomeGenome(keys=individual.keys)

    for key in individual.keys:
        representation = args.get('%s_representation' % key)
        for index in individual:
            mutation_rate = float(args.get('%s_mutation_rate' % key, .1))
            if random.random() < mutation_rate:
                new_individual[key].append(random.randint(0, len(representation) - 1))
            else:
                new_individual[key].append(index)

    return new_individual


@mutator
def multiple_chromosome_set_indel(random, individual, args):
    new_individual = individual.copy()

    for key in individual.keys:
        representation = args.get('%s_representation' % key)
        indel_rate = float(args.get('%s_indel_rate' % key, .1))
        if random.random() < indel_rate:
            if random.random() > float(0.5):
                if len(individual) > 1:
                    new_individual[key].pop(random.randint(0, len(new_individual) - 1))
            else:
                new_individual[key].append(random.sample(xrange(len(representation)), 1)[0])

    return individual


@crossover
def multiple_chromosome_n_point_crossover(random, mom, dad, args):
    children = MultipleChromosomeGenome(keys=mom.keys)
    for key in children.keys:
        children[key] = n_point_crossover(random, mom[key], dad[key], args)

    return children