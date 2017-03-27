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

from cameo.util import decompose_reaction_groups

__all__ = ['ReactionSetDecoder', 'GeneSetDecoder']


class SetDecoder(object):
    """
    Decoder for set representation. Decodes integer into string.
    """

    def __init__(self, representation, model, *args, **kwargs):
        super(SetDecoder, self).__init__(*args, **kwargs)
        self.representation = representation
        self.model = model

    def __call__(self, individual, flat=False, decompose=False):
        return [[self.representation[index] for index in individual]]


class ReactionSetDecoder(SetDecoder):
    """
    Decoder for set representation. Converts an integer set into reactions.

    Attributes
    ----------

    representation : list
        Reactions.
    model : cobra.Model
    groups : list
        A list of dict(reaction: relative_coefficient) where the dict contains coupled reactions.

    """

    def __init__(self, representation, model, groups=None, *args, **kwargs):
        super(ReactionSetDecoder, self).__init__(representation, model, *args, **kwargs)
        self.groups = groups

    def __call__(self, individual, flat=False, decompose=False):
        """
        Parameters
        ----------

        individual : list
            a list of integers.
        flat : bool
            if True, returns strings. Otherwise returns Reaction.
        decompose : bool
            If groups are available returns all possible substitutions.

        Returns
        -------
        list
            list of decoded representation (more then one combination is possible with decompositions)

        """
        reactions = [self.model.reactions.get_by_id(self.representation[index]) for index in individual]

        if decompose and self.groups:
            combinations = decompose_reaction_groups(self.groups, reactions)
        else:
            combinations = [reactions]

        if flat:
            return [tuple(r.id for r in reactions) for reactions in combinations]
        return [tuple(reactions) for reactions in combinations]


class GeneSetDecoder(SetDecoder):
    """
    Decoder for set representation. Converts an integer set into genes.

    Attributes
    ----------

    representation : list
        Genes to knockout.
    model : cobra.Model
    groups : list
        A list of dict(gene: relative_coefficient) where the dict contains coupled genes.
    """

    def __init__(self, representation, model, groups=None, *args, **kwargs):
        super(GeneSetDecoder, self).__init__(representation, model, *args, **kwargs)
        self.groups = groups

    def __call__(self, individual, flat=False, decompose=False):
        """
        Parameters
        ----------

        individual : list
            a list of integers.
        flat : bool
            if True, returns strings. Otherwise returns Gene.
        decompose : bool
            If groups are available returns all possible substitutions.


        Returns
        -------
        list
            list of decoded representation (more then one combination is possible with decompositions).

        """
        genes = [self.model.genes.get_by_id(self.representation[index]) for index in individual]

        if flat:
            return [tuple(g.id for g in genes)]
        return [tuple(genes)]
