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

from collections import OrderedDict

from inspyred.ec.generators import diversify
from cameo.strain_design.heuristic.evolutionary.genomes import MultipleChromosomeGenome
from six.moves import range
from six.moves import zip

__all__ = ['set_generator', 'unique_set_generator']


def set_generator(random, args):
    """
    Generates a list containing non-repeated elements of a discrete or
    continuous representation.

    Parameters
    ----------

    random : Random
    args : dict
        representation: set containing the possible values
        max_candidate_size: int, default: 9
        variable_candidate_size: bool, default: True


    Returns
    -------
    list
        A list containing a sample of the elements. If variable_candidate_size is
        True the list size is up to max_candidate_size, otherwise the candidate
        size equals candidate_size
    """
    representation = args.get('representation')
    indices = list(range(len(representation)))
    max_size = args.get('max_size', 9)
    variable_size = args.get('variable_size', True)
    if variable_size:
        size = random.randint(1, max_size)
    else:
        size = max_size
    candidate = random.sample(indices, size)
    return sorted(candidate)


@diversify
def unique_set_generator(random, args):
    """
    Generates a list containing non-repeated elements of a discrete or
    continuous representation. When this generator is used, the population
    will contain unique individuals.

    See Also
    --------
    inspyred.ec.generators.diversify

    Parameters
    ----------
    random : Random
    args : dict
        representation: set containing the possible values
        max_candidate_size: int, default: 9
        variable_candidate_size: bool, default: True

    Returns
    -------
    list
        A list containing a sample of the elements. If variable_candidate_size is
        True the list size is up to max_candidate_size, otherwise the candidate
        size equals candidate_size
    """
    representation = args.get('representation')
    indices = list(range(len(representation)))
    max_size = args.get('max_size', 9)
    variable_size = args.get('variable_size', True)
    if variable_size:
        size = random.randint(1, max_size)
    else:
        size = max_size
    candidate = random.sample(indices, size)
    return sorted(candidate)


def multiple_chromosome_set_generator(random, args):
    """
    Generates a candidate in with a genome containing multiple chromosomes.

    Parameters
    ----------
    random : Random
    args : dict
        this dictionary contains an extra key: "keys", which will be the names of the chromosomes.
        for each argument on <unique_set>, this dictionary must contain the same fields, starting with
        the values in "keys".

    Returns
    -------
    MultipleChromosomeGenome

    """
    keys = args.get('keys')
    candidate = MultipleChromosomeGenome(keys=keys)
    for key in keys:
        key_args = {
            'representation': args.get("%s_representation" % key),
            'max_size': args.get("%s_max_size" % key),
            'variable_size': args.get('variable_size')
        }
        candidate[key] = unique_set_generator(random, key_args)

    return candidate


def linear_set_generator(random, args):
    """
    Generates a list continuous values of the size of a representation.
    This function requires that a bounder is defined on the EvolutionaryAlgorithm.

    See Also
    --------
    inspyred.ec


    Parameters
    ----------
    random : Random
    args : dict
        representation: set containing the possible values
        max_candidate_size: int, default: 9
        variable_candidate_size: bool, default: True

    Returns
    -------
    list
        A list containing tuples - sample of the elements and linear value.
        If variable_candidate_size is True the list size is up to max_candidate_size,
        otherwise the candidate size equals candidate_size
    """
    bounder = args.get("_ec").bounder
    representation = args.get('representation')
    max_size = args.get('max_size', 9)
    variable_size = args.get('variable__size', True)
    if variable_size:
        size = random.randint(1, max_size)
    else:
        size = max_size

    indices = random.sample(range(len(representation)), size)
    values = random.uniform(next(bounder.lower_bound), next(bounder.upper_bound), len(indices))
    return OrderedDict({i: v for i, v in zip(indices, values)})
