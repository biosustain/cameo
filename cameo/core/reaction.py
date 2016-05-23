# -*- coding: utf-8 -*-
# Copyright 2013 Novo Nordisk Foundation Center for Biosustainability, DTU.
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

from __future__ import absolute_import, print_function

from functools import partial
import hashlib
import cobra as _cobrapy

import cameo
from cameo import flux_analysis
from cameo.parallel import SequentialView

import logging
import six
from cameo.util import inheritdocstring

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


@six.add_metaclass(inheritdocstring)
class Reaction(_cobrapy.core.Reaction):
    """This class extends the cobrapy Reaction class to work with SolverBasedModel.

    Notes
    -----
    See also documentation for cobra.core.Reaction.Reaction for an extensive list of inherited attributes.
    """

    @classmethod
    def clone(cls, reaction, model=None):
        """Clone a reaction.

        Parameters
        ----------
        reaction : Reaction, cobra.core.Reaction.Reaction
        model : model, optional

        Returns
        -------
        Reaction

        """
        new_reaction = cls(reaction.id, name=reaction.name)
        for attribute, value in six.iteritems(reaction.__dict__):
            if attribute == '_model':  # _model needs to be set after all other attributes have been set.
                continue
            try:
                setattr(new_reaction, attribute, value)
            except AttributeError:
                logger.info(
                    "Can't set attribute %s for reaction %s (while cloning it to a cameo style reaction). Skipping it ..." % (
                        attribute, reaction))
        if not isinstance(reaction.model, cameo.core.solver_based_model.SolverBasedModel):
            new_reaction._model = None
        if model is not None:
            new_reaction._model = model
        # for gene in new_reaction.genes:
        #     gene._reaction.remove(reaction)
        #     gene._reaction.add(new_reaction)
        # for metabolite in new_reaction.metabolites:
        #     metabolite._reaction.remove(reaction)
        #     metabolite._reaction.add(new_reaction)
        return new_reaction

    def __init__(self, id=None, name='', subsystem="", lower_bound=0, upper_bound=1000):
        """
        Parameters
        ----------
        name : str, optional
            The name of the reaction.
        """
        super(Reaction, self).__init__(id=id, name=name, subsystem=subsystem)
        self._lower_bound = lower_bound
        self._upper_bound = upper_bound
        self._model = None
        self._reverse_variable = None
        self._forward_variable = None

    def __str__(self):
        return ''.join((self.id, ": ", self.build_reaction_string()))

    @property
    def reversibility(self):
        return self._lower_bound < 0 < self._upper_bound

    def _get_reverse_id(self):
        """Generate the id of reverse_variable from the reaction's id."""
        return '_'.join((self.id, 'reverse', hashlib.md5(self.id.encode('utf-8')).hexdigest()[0:5]))

    def _get_forward_id(self):
        """Generate the id of forward_variable from the reaction's id."""
        return self.id
        # return '_'.join((self.id, 'forward', hashlib.md5(self.id.encode('utf-8')).hexdigest()[0:5]))

    @property
    def flux_expression(self):
        """An optlang variable representing the forward flux (if associated with model), otherwise None.
        Representing the net flux if model.reversible_encoding == 'unsplit'"""
        model = self.model
        if model is not None:
            return 1. * self.forward_variable - 1. * self.reverse_variable
        else:
            return None

    @property
    def forward_variable(self):
        """An optlang variable representing the forward flux (if associated with model), otherwise None."""
        model = self.model
        if model is not None:
            return model.solver.variables[self._get_forward_id()]
        else:
            return None

    @property
    def reverse_variable(self):
        """An optlang variable representing the reverse flux (if associated with model), otherwise None."""
        model = self.model
        if model is not None:
            return model.solver.variables[self._get_reverse_id()]
        else:
            return None

    @property
    def lower_bound(self):
        return self._lower_bound

    @lower_bound.setter
    def lower_bound(self, value):
        model = self.model

        if model is not None:

            forward_variable, reverse_variable = self.forward_variable, self.reverse_variable
            if self._lower_bound < 0 < self._upper_bound:  # reversible
                if value < 0:
                    reverse_variable.ub = -1 * value
                elif value >= 0:
                    reverse_variable.ub = 0
                    try:
                        forward_variable.lb = value
                    except ValueError:
                        forward_variable.ub = value
                        self._upper_bound = value
                        forward_variable.lb = value
            elif self._lower_bound == 0 and self._upper_bound == 0:  # knockout
                if value < 0:
                    reverse_variable.ub = -1 * value
                elif value >= 0:
                    forward_variable.ub = value
                    forward_variable.lb = value
            elif self._lower_bound >= 0:  # forward irreversible
                if value < 0:
                    reverse_variable.ub = -1 * value
                    forward_variable.lb = 0
                else:
                    try:
                        forward_variable.lb = value
                    except ValueError:
                        forward_variable.ub = value
                        self._upper_bound = value
                        forward_variable.lb = value

            elif self._upper_bound <= 0:  # reverse irreversible
                if value > 0:
                    reverse_variable.lb = 0
                    reverse_variable.ub = 0
                    forward_variable.ub = value
                    self._upper_bound = value
                    forward_variable.lb = value
                else:
                    try:
                        reverse_variable.ub = -1 * value
                    except ValueError:
                        reverse_variable.lb = -1 * value
                        self._upper_bound = value
                        reverse_variable.ub = -1 * value
            else:
                print({'value': value, 'self._lower_bound': self._lower_bound, 'self._upper_bound': self._upper_bound})
                raise ValueError('lower_bound issue')

        self._lower_bound = value

    @property
    def upper_bound(self):
        return self._upper_bound

    @upper_bound.setter
    def upper_bound(self, value):
        model = self.model
        if model is not None:

            forward_variable, reverse_variable = self.forward_variable, self.reverse_variable
            if self._lower_bound < 0 < self._upper_bound:  # reversible
                if value > 0:
                    forward_variable.ub = value
                elif value <= 0:
                    forward_variable.ub = 0
                    try:
                        reverse_variable.lb = -1 * value
                    except ValueError:
                        reverse_variable.ub = -1 * value
                        self._lower_bound = value
                        reverse_variable.lb = -1 * value
            elif self._lower_bound == 0 and self._upper_bound == 0:  # knockout
                if value > 0:
                    forward_variable.ub = value
                elif value <= 0:
                    reverse_variable.ub = -1 * value
            elif self._lower_bound >= 0:  # forward irreversible
                if value > 0:
                    try:
                        forward_variable.ub = value
                    except ValueError:
                        forward_variable.lb = value
                        self._lower_bound = value
                        forward_variable.ub = value
                else:
                    forward_variable.lb = 0
                    forward_variable.ub = 0
                    reverse_variable.ub = -1 * value
                    self._lower_bound = value
                    reverse_variable.lb = -1 * value

            elif self._upper_bound <= 0:  # reverse irreversible
                if value < 0:
                    try:
                        reverse_variable.lb = -1 * value
                    except ValueError:
                        reverse_variable.ub = -1 * value
                        self._lower_bound = value
                        reverse_variable.lb = -1 * value
                else:
                    forward_variable.ub = value
                    reverse_variable.lb = 0
            else:
                print({'value': value, 'self._lower_bound': self._lower_bound, 'self._upper_bound': self._upper_bound})
                raise ValueError('upper_bound issue')

        self._upper_bound = value

    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, value):
        if value is None:
            self.objective_coefficient  # Get objective coefficient from model
        elif not isinstance(value, cameo.core.SolverBasedModel):
            raise ValueError("Must be an instance of cameo.core.SolverBasedModel, not %s" % type(value))
        self._model = value

    @property
    def objective_coefficient(self):
        if self.model is not None and isinstance(self.model, cameo.core.SolverBasedModel) and self.model.objective is not None:
            coefficients_dict = self.model.objective.expression.as_coefficients_dict()
            forw_coef = coefficients_dict.get(self.forward_variable, 0)
            rev_coef = coefficients_dict.get(self.reverse_variable, 0)
            if forw_coef == -rev_coef:
                self._objective_coefficient = float(forw_coef)
            else:
                self._objective_coefficient = 0
        return self._objective_coefficient

    @objective_coefficient.setter
    def objective_coefficient(self, value):
        model = self.model
        if model is not None:
            coef_difference = value - self.objective_coefficient
            model.objective += coef_difference * self.flux_expression
        self._objective_coefficient = value

    @property
    def effective_lower_bound(self):
        model = self.model
        fva_result = flux_analysis.flux_variability_analysis(model, reactions=[self], view=SequentialView(),
                                                             remove_cycles=False)
        return fva_result['lower_bound'][self.id]

    @property
    def effective_upper_bound(self):
        model = self.model
        fva_result = flux_analysis.flux_variability_analysis(model, reactions=[self], view=SequentialView(),
                                                             remove_cycles=False)
        return fva_result['upper_bound'][self.id]

    @property
    def flux(self):
        if self.model is not None:
            return self.forward_variable.primal - self.reverse_variable.primal
        else:
            return None

    @property
    def reduced_cost(self):
        if self.model is not None:
            return self.forward_variable.dual - self.reverse_variable.dual
        else:
            return None

    def add_metabolites(self, metabolites, combine=True, **kwargs):
        if not combine:
            old_coefficients = self.metabolites
        super(Reaction, self).add_metabolites(metabolites, combine=combine, **kwargs)
        model = self.model
        if model is not None:
            for metabolite, coefficient in six.iteritems(metabolites):
                if isinstance(metabolite, six.string_types):  # support metabolites added as strings.
                    metabolite = model.metabolites.get_by_id(metabolite)
                if not combine:
                    try:
                        old_coefficient = old_coefficients[metabolite]
                    except KeyError:
                        pass
                    else:
                        coefficient = coefficient - old_coefficient
                model.solver.constraints[metabolite.id] += coefficient * self.flux_expression

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
            if value in self.model.reactions:
                raise ValueError("The model already contains a reaction with the id:", value)
            forward_variable = self.forward_variable
            reverse_variable = self.reverse_variable

            self._id = value
            self.model.reactions._generate_index()

            forward_variable.name = self._get_forward_id()
            reverse_variable.name = self._get_reverse_id()
        else:
            self._id = value

    def knock_out(self, time_machine=None):
        """Knockout reaction by setting its bounds to zero.

        Parameters
        ----------
        time_machine = TimeMachine
            A time TimeMachine instance can be provided to undo the knockout eventually.

        Returns
        -------
        None
        """

        def _(reaction, lb, ub):
            reaction.upper_bound = ub
            reaction.lower_bound = lb

        if time_machine is not None:
            time_machine(do=super(Reaction, self).knock_out, undo=partial(_, self, self.lower_bound, self.upper_bound))
        else:
            super(Reaction, self).knock_out()

    def pop(self, metabolite_id):
        """Removes a given metabolite from the reaction stoichiometry, and returns the coefficient.
        """
        if self._model is None:
            return super(Reaction, self).pop(metabolite_id)
        else:
            if isinstance(metabolite_id, six.string_types):
                met = self.model.metabolites.get_by_id(metabolite_id)
            else:
                met = metabolite_id
            coef = self.metabolites[met]
            self.add_metabolites({met: -coef}, combine=True)
            return coef

    def remove_from_model(self, model=None, remove_orphans=False):
        reaction_model = self.model
        forward = self.forward_variable
        reverse = self.reverse_variable
        super(Reaction, self).remove_from_model(model, remove_orphans)
        reaction_model.solver.remove([forward, reverse])

    def delete(self, remove_orphans=False):
        model = self.model
        forward = self.forward_variable
        reverse = self.reverse_variable
        super(Reaction, self).delete(remove_orphans)
        model.solver.remove([forward, reverse])
        # if remove_orphans:
        #     model.solver.remove([metabolite.model.solver for metabolite in self.metabolites.keys()])

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
                <td><strong>Stoichiometry</strong></td><td>%s</td>
            </tr>
            <tr>
                <td><strong>Lower bound</strong></td><td>%f</td>
            </tr>
            <tr>
                <td><strong>Upper bound</strong></td><td>%f</td>
            </tr>
        </table>
        """ % (self.id, self.name, self.reaction, self.lower_bound, self.upper_bound)
