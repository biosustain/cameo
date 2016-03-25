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

logger = logging.getLogger(__name__)


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
