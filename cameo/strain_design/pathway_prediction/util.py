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
from cameo import Reaction


def create_adaptor_reactions(original_metabolites, expanded_model, mapping, compartment_regexp):
    adapter_reactions = []
    for metabolite in original_metabolites:  # model is the original host model
        if not compartment_regexp.match(metabolite.id):
            continue

        name = metabolite.id[0:-2]
        try:
            mnx_name = mapping[name]
        except KeyError:
            continue
            # print name, 'N/A'
        adapter_reaction = Reaction('adapter_' + metabolite.id + '_' + mnx_name)
        adapter_reaction.lower_bound = -1000
        try:
            adapter_reaction.add_metabolites({expanded_model.metabolites.get_by_id(metabolite.id): -1,
                                              expanded_model.metabolites.get_by_id(mnx_name): 1})
        except KeyError:
            pass
        else:
            adapter_reactions.append(adapter_reaction)

    return adapter_reactions