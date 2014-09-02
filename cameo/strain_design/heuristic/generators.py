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
from inspyred.ec.generators import diversify
from cameo.strain_design.heuristic.genomes import MultipleChromosomeGenome


def set_generator(random, args):
    """
    Generates a list containing non-repeated elements of a discrete or
    continuous representation.

    Parameters
    ----------

    random: Random
    args: dict
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
    max_size = args.get('candidate_size', 9)
    variable_size = args.get('variable_candidate_size', True)
    if variable_size:
        size = random.randint(1, max_size)
    else:
        size = max_size
    candidate = random.sample(xrange(len(representation)), size)
    return candidate


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
    random: Random
    args: dict
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
    max_size = args.get('candidate_size', 9)
    variable_size = args.get('variable_candidate_size', True)
    if variable_size:
        size = random.randint(1, max_size)
    else:
        size = max_size
    candidate = random.sample(xrange(len(representation)), size)
    return candidate


def multiple_chromosome_set_generator(random, args):
    """
    Generates a candidate in with a genome containing multiple chromosomes.

    Parameters
    ----------
    random: Random
    args: dict
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
            'representation': args.get("%s_representation"),
            'candidate_size': args.get("%s_candidate_size"),
            'variable_candidate_size': args.get('variable_candidate_size')
        }
        candidate[key] = unique_set_generator(random, key_args)

    return candidate