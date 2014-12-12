# Copyright 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import difflib
from cameo.api import _METANETX as METANETX

class Compound(object):

    def __init__(self, metabolite):
        self.metabolite = metabolite
        self.__dict__.update(METANETX['chem_prop'].loc[metabolite.id].to_dict())


class Products(object):

    def __init__(self):
        # for metabolite in METANETX['universal_model'].metabolites:
        #     setattr(self, metabolite.id, Compound(metabolite))
        metabolite_ids = [metabolite.id for metabolite in METANETX['universal_model'].metabolites]
        self.data_frame = METANETX['chem_prop'].loc[metabolite_ids]


    def search(self, name):
        matches = difflib.get_close_matches(name, self.data_frame.name, n=20, cutoff=.6)
        ranks = dict([(match, i) for i, match in enumerate(matches)])
        selection = self.data_frame[self.data_frame.name.isin(matches)]
        selection['search_rank'] = selection.name.map(ranks)
        return selection.sort('search_rank')


products = Products()