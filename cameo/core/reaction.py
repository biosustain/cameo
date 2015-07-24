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

import warnings

from functools import partial
import hashlib
import cobra as _cobrapy

import cameo
from cameo import flux_analysis
from cameo.parallel import SequentialView

import logging
import six

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


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
        new_reaction = cls(name=reaction.name)
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
        for gene in new_reaction.genes:
            gene._reaction.remove(reaction)
            gene._reaction.add(new_reaction)
        for metabolite in new_reaction.metabolites:
            metabolite._reaction.remove(reaction)
            metabolite._reaction.add(new_reaction)
        return new_reaction

    def __init__(self, name=None):
        """
        Parameters
        ----------
        name : str, optional
            The name of the reaction.
        """
        super(Reaction, self).__init__(name=name)
        self._lower_bound = 0
        self._upper_bound = 1000.
        self._objective_coefficient = 0.

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
    def variable(self):
        warnings.warn('reaction.variable is deprecated. Please use reaction.forward_variable.', DeprecationWarning)
        return self.forward_variable

    @property
    def forward_variable(self):
        """An optlang variable representing the forward flux (if associated with model), otherwise None.
        Representing the net flux if model.reversible_encoding == 'unsplit'"""
        model = self.model
        if model is not None:
            aux_id = self._get_forward_id()
            try:
                return model.solver.variables[aux_id]
            except KeyError:
                return None
        else:
            return None

    @property
    def reverse_variable(self):
        """An optlang variable representing the reverse flux (if associated with model), otherwise None."""
        model = self.model
        if model is not None:
            aux_id = self._get_reverse_id()
            try:
                return model.solver.variables[aux_id]
            except KeyError:
                return None
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
    def objective_coefficient(self):
        return self._objective_coefficient

    @objective_coefficient.setter
    def objective_coefficient(self, value):
        model = self.model
        if model is not None:
            model.solver._set_linear_objective_term(self.forward_variable, value)
            model.solver._set_linear_objective_term(self.reverse_variable, -1 * value)
        self._objective_coefficient = value

    @property
    def effective_lower_bound(self):
        model = self.model
        return \
            flux_analysis.flux_variability_analysis(model, reactions=[self], view=SequentialView(),
                                                    remove_cycles=False)[
                'lower_bound'][
                self.id]

    @property
    def effective_upper_bound(self):
        model = self.model
        return \
            flux_analysis.flux_variability_analysis(model, reactions=[self], view=SequentialView(),
                                                    remove_cycles=False)[
                'upper_bound'][
                self.id]

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
                if not combine:
                    try:
                        old_coefficient = old_coefficients[metabolite]
                    except KeyError:
                        pass
                    else:
                        coefficient = coefficient - old_coefficient
                model.solver.constraints[metabolite.id] += coefficient * self.flux_expression

    def knock_out(self, time_machine=None):
        def _(reaction, lb, ub):
            reaction.upper_bound = ub
            reaction.lower_bound = lb

        if time_machine is not None:
            time_machine(do=super(Reaction, self).knock_out, undo=partial(_, self, self.lower_bound, self.upper_bound))
        else:
            super(Reaction, self).knock_out()

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
