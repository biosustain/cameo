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
from itertools import combinations, chain
from ordered_set import OrderedSet
import pprint

N = "-"


def identify_currency_metabolites_by_pattern(model, top=20, min_combination=3, max_combination=5):
    """
    Identify currency metabolites and their patterns in reactions.

    Parameters
    ---------
    model : cobra.Model
    top : int
        Top max connected metabolites in the network, default: 20
    min_combination : int
        Size of the smallest pattern, default: 3
    max_combination : int
        Size of the longest pattern, default: 5

    Returns
    -------
    possible_motifs : dict
        The identified motifs and their rank
    currency_metabolites : list
        The ids of the most connected metabolites
    """


    metabolites_degree = dict([[m.id, len(m.reactions)] for m in model.metabolites])
    metabolite_ids = metabolites_degree.keys()
    metabolite_ids = sorted(metabolite_ids, key=metabolites_degree.get)
    metabolite_ids.reverse()
    currency_metabolites = metabolite_ids[:top]

    possible_motifs = {}

    for reaction in model.reactions:
        metabolites = [m.id for m in reaction.metabolites.keys() if m.id in currency_metabolites]
        r = xrange(min_combination, min(len(metabolites), max_combination)+1)
        motifs = list(set(chain(*[[frozenset(OrderedSet(c)) for c in combinations(metabolites, n)] for n in r])))
        if len(motifs) > 0:
            while len(motifs) > 0:
                motif = motifs.pop()
                add = True
                for other_motif in motifs:
                    if motif.issubset(other_motif):
                        add = False
                    if other_motif.issubset(motif):
                        motifs.remove(other_motif)

                if add:
                    if motif in possible_motifs:
                        possible_motifs[motif] += 1
                    else:
                        possible_motifs[motif] = 1

    motifs = possible_motifs.keys()
    motifs = sorted(motifs, key=possible_motifs.get)
    motifs.reverse()

    for m in motifs[:10]:
        print m, ": ", possible_motifs[m]

    return possible_motifs, currency_metabolites