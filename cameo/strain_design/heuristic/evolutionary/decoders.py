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

from __future__ import absolute_import

__all__ = ['ReactionKnockoutDecoder', 'GeneKnockoutDecoder']

from cobra.manipulation.delete import find_gene_knockout_reactions


class KnockoutDecoder(object):
    def __init__(self, representation, model, *args, **kwargs):
        super(KnockoutDecoder, self).__init__(*args, **kwargs)
        self.representation = representation
        self.model = model

    def __call__(self, individual, flat=False):
        raise NotImplementedError


class ReactionKnockoutDecoder(KnockoutDecoder):
    """
    Decoder for set representation. Converts an integer set into reactions available for knockout

    Parameters
    ----------

    representation : list
        reactions to knockout
    model : SolverBasedModel

    """

    def __init__(self, representation, model, *args, **kwargs):
        super(ReactionKnockoutDecoder, self).__init__(representation, model, *args, **kwargs)

    def __call__(self, individual, flat=False):
        """
        Parameters
        ----------

        individual: list
            a list of integers
        flat: bool
            if True, returns strings. Otherwise returns Reaction

        Returns
        -------
        list
            [knockouts, decoded representation]

        """
        reactions = [self.model.reactions.get_by_id(self.representation[index]) for index in individual]
        if flat:
            return [tuple(r.id for r in reactions), tuple(r.id for r in reactions)]
        return [tuple(reactions), tuple(reactions)]


class GeneKnockoutDecoder(KnockoutDecoder):
    """
    Decoder for set representation. Converts an integer set into genes available for knockout
    Parameters
    ----------

    representation : list
        genes to knockout
    model : SolverBasedModel
    """

    def __init__(self, representation, model, *args, **kwargs):
        super(GeneKnockoutDecoder, self).__init__(representation, model, *args, **kwargs)

    def __call__(self, individual, flat=False):
        """
        Parameters
        ----------

        individual: list
            a list of integers
        flat: bool
            if True, returns strings. Otherwise returns Gene

        Returns
        -------
        list
            [knockouts, decoded representation]

        """
        genes = [self.model.genes.get_by_id(self.representation[index]) for index in individual]
        reactions = find_gene_knockout_reactions(self.model, genes)

        if flat:
            return [tuple(r.id for r in reactions), tuple(g.id for g in genes)]
        return [tuple(reactions), tuple(genes)]
