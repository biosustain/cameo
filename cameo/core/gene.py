# Copyright 2016 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from functools import partial

import cobra
import logging
import six
from cameo.util import inheritdocstring

logger = logging.getLogger(__name__)


@six.add_metaclass(inheritdocstring)
class Gene(cobra.core.Gene):

    @classmethod
    def clone(cls, gene, model=None):
        new_gene = cls(id=gene.id)
        for attribute, value in gene.__dict__.items():
            try:
                setattr(new_gene, attribute, value)
            except AttributeError:
                logger.info(
                    "Can't set attribute %s for gene %s (while cloning it to a cameo style gene). Skipping it ..." %
                    (attribute, gene)
                )
        if model is not None:
            new_gene._model = model
        return new_gene

    @property
    def id(self):
        return getattr(self, "_id", None)  # Returns None if _id is not set

    @id.setter
    def id(self, value):
        if value == self.id:
            pass
        elif not isinstance(value, six.string_types):
            raise TypeError("ID must be a string")
        elif getattr(self, "_model", None) is not None:  # (= if hasattr(self, "_model") and self._model is not None)
            if value in self.model.genes:
                raise ValueError("The model already contains a gene with the id:", value)

            self._id = value
            self.model.genes._generate_index()
        else:
            self._id = value

    def knock_out(self, time_machine=None):
        """Knockout gene by setting all its affected reactions' bounds to zero.

        Parameters
        ----------
        time_machine = TimeMachine
            A time TimeMachine instance can be provided to undo the knockout eventually.

        Returns
        -------
        None
        """

        from cobra.manipulation.delete import find_gene_knockout_reactions
        for reaction in find_gene_knockout_reactions(self.model, [self]):
            def _(reaction, lb, ub):
                reaction.upper_bound = ub
                reaction.lower_bound = lb

            if time_machine is not None:
                time_machine(do=reaction.knock_out, undo=partial(_, reaction, reaction.lower_bound, reaction.upper_bound))
            else:
                reaction.knock_out()
