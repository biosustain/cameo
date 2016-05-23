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

import cobra
import logging
import six
from cameo.util import inheritdocstring

logger = logging.getLogger(__name__)


@six.add_metaclass(inheritdocstring)
class Metabolite(cobra.core.Metabolite):

    @classmethod
    def clone(cls, metabolite, model=None):
        new_metabolite = cls(id=metabolite.id)
        for attribute, value in metabolite.__dict__.items():
            try:
                setattr(new_metabolite, attribute, value)
            except AttributeError:
                logger.info(
                    "Can't set attribute %s for metabolite %s (while cloning it to a cameo style metabolite). Skipping it ..." %
                    (attribute, metabolite)
                )
        if model is not None:
            new_metabolite._model = model
        return new_metabolite

    def remove_from_model(self, method="subtractive", **kwargs):
        model = self.model
        super(Metabolite, self).remove_from_model(method, **kwargs)
        model.solver.remove(model.solver.constraints[self.id])

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
            if value in self.model.metabolites:
                raise ValueError("The model already contains a metabolite with the id:", value)
            self.model.solver.constraints[self.id].name = value

            self._id = value
            self.model.metabolites._generate_index()
        else:
            self._id = value
