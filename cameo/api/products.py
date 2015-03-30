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
from cameo.data import metanetx


def inchi_to_svg(inchi, file=None):
    try:
        import openbabel
    except ImportError, e:
        print(e)
        raise ImportError("OpenBabel seems to be not installed.")
    convert = openbabel.OBConversion()
    convert.SetInFormat("inchi")
    convert.SetOutFormat("svg")
    mol = openbabel.OBMol()
    if not convert.ReadString(mol, inchi):
        raise Exception("%s could not be parsed as an inchi string.")
    return convert.WriteString(mol)


class Compound(object):

    def __init__(self, inchi):
        self.inchi

    def _repr_svg_(self):
        try:
            return inchi_to_svg(self.InChI)
        except ImportError:
            return self.__repr__()

    def _repr_html_(self):
        return self._repr_svg_()


class Products(object):

    def __init__(self):
        self.data_frame = metanetx.chem_prop

    def search(self, query):
        matches = self._search_by_source(query)
        if len(matches) > 0:
            return matches
        matches = self._search_by_inchi(query)
        if len(matches) > 0:
            return matches
        matches = self._search_by_name_fuzzy(query)
        if len(matches) > 0:
            return matches
        matches = self._search_by_inchi_fuzzy(query)
        if len(matches) > 0:
            return matches
        else:
            raise Exception("No compound matches found for query %s" % query)

    def _search_by_name_fuzzy(self, name):
        matches = difflib.get_close_matches(name, self.data_frame.name, n=5, cutoff=.6)
        ranks = dict([(match, i) for i, match in enumerate(matches)])
        selection = self.data_frame[self.data_frame.name.isin(matches)]
        selection['search_rank'] = selection.name.map(ranks)
        return selection.sort('search_rank')

    def _search_by_source(self, source_id):
        return self.data_frame[self.data_frame.source == source_id.lower()]

    def _search_by_inchi(self, inchi):
        return self.data_frame[self.data_frame.InChI == inchi]

    def _search_by_inchi_fuzzy(self, inchi):
        # TODO: use openbabel if available
        matches = difflib.get_close_matches(inchi, self.data_frame.InChI, n=5, cutoff=.6)
        ranks = dict([(match, i) for i, match in enumerate(matches)])
        selection = self.data_frame[self.data_frame.InChI.isin(matches)]
        selection['search_rank'] = selection.name.map(ranks)
        return selection.sort('search_rank')
        return self.data_frame[self.data_frame.InChI == inchi]


products = Products()