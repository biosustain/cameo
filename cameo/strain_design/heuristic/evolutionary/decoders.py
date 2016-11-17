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

__all__ = ['ReactionSetDecoder', 'GeneSetDecoder']


class SetDecoder(object):
    """
    Decoder for set representation. Decodes integer into string.
    """

    def __init__(self, representation, model, *args, **kwargs):
        super(SetDecoder, self).__init__(*args, **kwargs)
        self.representation = representation
        self.model = model

    def __call__(self, individual, flat=False):
        return [self.representation[index] for index in individual]


class ReactionSetDecoder(SetDecoder):
    """
    Decoder for set representation. Converts an integer set into reactions.

    Parameters
    ----------

    representation : list
        Reactions.
    model : SolverBasedModel

    """

    def __init__(self, representation, model, *args, **kwargs):
        super(ReactionSetDecoder, self).__init__(representation, model, *args, **kwargs)

    def __call__(self, individual, flat=False):
        """
        Parameters
        ----------

        individual: list
            a list of integers
        flat: bool
            if True, returns strings. Otherwise returns Reaction.

        Returns
        -------
        list
            Decoded representation

        """
        reactions = [self.model.reactions.get_by_id(self.representation[index]) for index in individual]
        if flat:
            return tuple(r.id for r in reactions)
        return tuple(reactions)


class GeneSetDecoder(SetDecoder):
    """
    Decoder for set representation. Converts an integer set into genes.
    Parameters
    ----------

    representation : list
        Genes to knockout.
    model : SolverBasedModel
    """

    def __init__(self, representation, model, *args, **kwargs):
        super(GeneSetDecoder, self).__init__(representation, model, *args, **kwargs)

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
            Decoded representation.

        """
        genes = [self.model.genes.get_by_id(self.representation[index]) for index in individual]

        if flat:
            return tuple(g.id for g in genes)
        return tuple(genes)
