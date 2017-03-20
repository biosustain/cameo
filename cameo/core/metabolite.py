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
import logging
from functools import partial

import cobra
import six
from cameo.util import inheritdocstring

logger = logging.getLogger(__name__)

__all__ = ['Metabolite']


@six.add_metaclass(inheritdocstring)
class Metabolite(cobra.core.Metabolite):
    # TODO: figure out how to handle the _reaction attribute
    @classmethod
    def clone(cls, metabolite, model=None):
        new_metabolite = cls(id=metabolite.id)
        for attribute, value in metabolite.__dict__.items():
            try:
                setattr(new_metabolite, attribute, value)
            except AttributeError:
                logger.info(
                    "Can't set attribute %s for metabolite %s (while cloning it to a "
                    "cameo style metabolite). Skipping it ..." %
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
    def n_carbon(self):
        """number of carbon atoms

        Returns
        -------
        int
            number of carbons in this metabolite
        """
        return self.elements.get('C', 0)

    @property
    def constraint(self):
        if self.model is not None:
            return self.model.solver.constraints[self.id]
        else:
            return None

    def _relax_mass_balance_constrain(self, time_machine):
        if time_machine:
            time_machine(do=partial(setattr, self.constraint, "lb", -1000 * len(self.reactions)),
                         undo=partial(setattr, self.constraint, "lb", self.constraint.lb))
            time_machine(do=partial(setattr, self.constraint, "ub", 1000 * len(self.reactions)),
                         undo=partial(setattr, self.constraint, "ub", self.constraint.ub))
        else:
            self.constraint.lb = None
            self.constraint.ub = None

    def knock_out(self, time_machine=None, force_steady_state=False):
        """'Knockout' a metabolite. This can be done in 2 ways:

        1. Implementation follows the description in [1]
            "All fluxes around the metabolite M should be restricted to only produce the metabolite,
             for which balancing constraint of mass conservation is relaxed to allow nonzero values
             of the incoming fluxes whereas all outgoing fluxes are limited to zero."

        2. Force steady state
            All reactions consuming the metabolite are restricted to only produce the metabolite. A demand
            reaction is added to sink the metabolite produced to keep the problem feasible under
            the S.v = 0 constraint.


        Knocking out a metabolite overrules the constraints set on the reactions producing the metabolite.

        Parameters
        ----------
        time_machine : TimeMachine
            An action stack to reverse actions
        force_steady_state: bool
            If True, uses approach 2.

        References
        ----------
        .. [1] Kim, P.-J., Lee, D.-Y., Kim, T. Y., Lee, K. H., Jeong, H., Lee, S. Y., & Park, S. (2007).
          Metabolite essentiality elucidates robustness of Escherichia coli metabolism. PNAS, 104(34), 13638-13642
        """
        # restrict reactions to produce metabolite
        for reaction in self.reactions:
            if reaction.metabolites[self] > 0:  # for positive stoichiometric coefficient set lb to 0
                if reaction.upper_bound < 0:
                    reaction.change_bounds(lb=0, ub=0, time_machine=time_machine)
                else:
                    reaction.change_bounds(lb=0, time_machine=time_machine)
            elif reaction.metabolites[self] < 0:  # for negative stoichiometric coefficient set ub to 0
                if reaction.lower_bound > 0:
                    reaction.change_bounds(lb=0, ub=0, time_machine=time_machine)
                else:
                    reaction.change_bounds(ub=0, time_machine=time_machine)
        if force_steady_state:
            self.model.add_exchange(self, prefix="KO_", time_machine=time_machine)
        else:
            self._relax_mass_balance_constrain(time_machine)

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

    def _repr_html_(self):
        return """
        <table>
            <tr>
                <td><strong>Id</strong></td><td>%s</td>
            </tr>
            <tr>
                <td><strong>Name</strong></td><td>%s</td>
            </tr>
            <tr>
                <td><strong>Formula</strong></td><td>%s</td>
            </tr>
        </table>""" % (self.id, self.name, self.formula)
