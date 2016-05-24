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

"""This module provides access to supported target products."""

from __future__ import absolute_import, print_function

__all__ = ['products']

import difflib
from pandas import DataFrame
from cameo.data import metanetx
from cameo.visualization import inchi_to_svg


class Compound(object):
    def __init__(self, inchi):
        self.InChI = inchi

    def _repr_svg_(self):
        try:
            return inchi_to_svg(self.InChI)
        except ImportError:
            return self.__repr__()

    def _repr_html_(self):
        return self._repr_svg_()


class Products(object):
    """Supported target products."""

    def __init__(self):
        self.data_frame = metanetx.chem_prop

    def search(self, query):
        """Fuzzy search of available target products.

        Parameters
        ----------
        query : str
            Compound ID, name or InChI string.

        Returns
        -------
        pandas.DataFrame
            A dataframe containing the scored search results.
        """
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
        matches = difflib.get_close_matches(name, self.data_frame.name.dropna(), n=5, cutoff=.8)
        ranks = dict([(match, i) for i, match in enumerate(matches)])
        selection = DataFrame(self.data_frame[self.data_frame.name.isin(matches)])
        selection['search_rank'] = selection.name.map(ranks)
        return selection.sort_values('search_rank')

    def _search_by_source(self, source_id):
        return self.data_frame[self.data_frame.source == source_id.lower()]

    def _search_by_inchi(self, inchi):
        return self.data_frame[self.data_frame.InChI == inchi]

    def _search_by_inchi_fuzzy(self, inchi):
        # TODO: use openbabel if available
        matches = difflib.get_close_matches(inchi, self.data_frame.InChI.dropna(), n=5, cutoff=.8)
        ranks = dict([(match, i) for i, match in enumerate(matches)])
        selection = DataFrame(self.data_frame[self.data_frame.InChI.isin(matches)])
        selection['search_rank'] = selection.name.map(ranks)
        return selection.sort_values('search_rank')


products = Products()
