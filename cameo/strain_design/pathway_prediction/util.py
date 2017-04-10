# Copyright 2015 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import absolute_import, print_function

__all__ = ['create_adapter_reactions', 'display_pathway']

import re

from cameo import Reaction
from cameo.ui import notice
from cameo.util import in_ipnb

try:
    from IPython.display import display
except:
    pass


def create_adapter_reactions(original_metabolites, universal_model, mapping, compartment_regexp):
    """Create adapter reactions that connect host and universal model.

    Arguments
    ---------
    original_metabolites : list
        List of host metababolites.
    universal_model : cobra.Model
        The universal model.
    mapping : dict
        A mapping between between host and universal model metabolite IDs.
    compartment_regexp : regex
        A compile regex that matches metabolites that should be connected to the universal model.

    Returns
    -------
    reactions : list
        The list of adapter reactions.
    """
    adapter_reactions = []
    for metabolite in original_metabolites:  # model is the original host model
        if compartment_regexp.match(metabolite.compartment):
            name = re.sub('_{}$'.format(metabolite.compartment), '', metabolite.id)  # TODO: still a hack
            mapped_name = None
            for prefix in ['bigg:', 'kegg:', 'rhea:', 'brenda:', '']:  # try no prefix at last
                if prefix + name in mapping:
                    mapped_name = mapping[prefix + name]
                    break
            if mapped_name is not None:
                adapter_reaction = Reaction(str('adapter_' + metabolite.id + '_' + mapped_name))
                adapter_reaction.lower_bound = -1000
                try:
                    adapter_reaction.add_metabolites({metabolite: -1,
                                                      universal_model.metabolites.get_by_id(mapped_name): 1})
                except KeyError:
                    pass
                else:
                    adapter_reactions.append(adapter_reaction)

    return adapter_reactions


def display_pathway(pathway, i):
    notice("Pathway %i" % i)
    if in_ipnb():
        display(pathway.data_frame)
    else:
        print(pathway.data_frame)
