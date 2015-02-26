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

from functools import partial
import hashlib
import cobra as _cobrapy
import sympy

import cameo
from cameo import flux_analysis
from cameo.parallel import SequentialView

import logging
logger = logging.getLogger(__name__)


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
        for attribute, value in reaction.__dict__.iteritems():
            try:
                setattr(new_reaction, attribute, value)
            except AttributeError:
                logger.debug("Can't set attribute %s for reaction %s (while cloning it to a cameo style reaction). Skipping it ..." % (attribute, reaction))
        if not isinstance(reaction.model, cameo.core.solver_based_model.SolverBasedModel):
            new_reaction._model = None
        if model is not None:
            new_reaction._model = model
        for gene in new_reaction.genes:
            # print gene._reaction
            # for r in list(gene._reaction):
            #     print r, type(r)
            gene._reaction.remove(reaction)
            # print gene._reaction
            gene._reaction.add(new_reaction)
            # print gene._reaction
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
    def variable(self):
        """An optlang variable representing the forward flux (if associated with model), otherwise None.
        Representing the net flux if model.reversible_encoding == 'unsplit'"""
        model = self.model
        if model is not None:
            return model.solver.variables[self.id]
        else:
            return None

    @property
    def reversibility(self):
        return self._lower_bound < 0 and self._upper_bound > 0

    def _get_reverse_id(self):
        """Generate the id of revers_variable from the reaction's id."""
        return '_'.join((self.id, 'reverse', hashlib.md5(self.id).hexdigest()[0:5]))

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
        # model = self.model
        # if model is not None:
        #     if model.reversible_encoding == 'split' and self.reversibility:
        #         return -1 * self.reverse_variable.ub
        #     else:
        #         return self.variable.lb
        # else:
        #     return self._lower_bound

    @lower_bound.setter
    def lower_bound(self, value):
        model = self.model

        if model is not None:

            if value >= 0 and self._lower_bound < 0 and self._upper_bound > 0:
                reverse_variable = self.reverse_variable
                reverse_variable.lb, reverse_variable.ub = 0, 0
            elif value < 0 and self._lower_bound >= 0 and self.reverse_variable is None:  # self._lower_bound >= 0 implies self._upper_bound >= 0
                reverse_variable = model.solver._add_variable(
                    model.solver.interface.Variable(self._get_reverse_id(), lb=0, ub=0))
                for met, coeff in self._metabolites.iteritems():
                    model.solver.constraints[met.id] += sympy.Mul._from_args((-1 * sympy.RealNumber(coeff), reverse_variable))

            variable = self.variable
            reverse_variable = self.reverse_variable

            if model.reversible_encoding == 'split' and value < 0 and self._upper_bound > 0:
                if self._lower_bound > 0:
                    variable.lb = 0
                try:
                    reverse_variable.ub = -1 * value
                except ValueError:
                    reverse_variable.lb = -1 * value
                    reverse_variable.ub = -1 * value
            else:
                try:
                    variable.lb = value
                except ValueError:
                    variable.ub = value
                    self._upper_bound = value
                    variable.lb = value

        self._lower_bound = value

    @property
    def upper_bound(self):
        return self._upper_bound
        # if self.model is not None:
        #     return self.variable.ub
        # else:
        #     return self._upper_bound

    @upper_bound.setter
    def upper_bound(self, value):
        model = self.model
        if model is not None:
            # Remove auxiliary variable if not needed anymore
            reverse_variable = self.reverse_variable
            variable = self.variable
            if value <= 0 and self._upper_bound > 0 and self._lower_bound < 0:
                reverse_variable.lb, reverse_variable.ub = 0, 0

            # Add auxiliary variable if needed
            elif value > 0 and self._upper_bound <= 0 and self.reverse_variable is None:  # self._upper_bound < 0 implies self._lower_bound < 0
                reverse_variable = model.solver._add_variable(
                    model.solver.interface.Variable(self._get_reverse_id(), lb=0, ub=0))
                for met, coeff in self._metabolites.iteritems():
                    model.solver.constraints[met.id] += sympy.Mul._from_args((-1 * sympy.RealNumber(coeff), reverse_variable))

            if model.reversible_encoding == 'split' and value > 0 and self._lower_bound < 0:
                variable.ub = value
                if self._upper_bound <= 0:
                    reverse_variable.ub = -1 * variable.lb
                    variable.lb = 0
            else:
                try:
                    variable.ub = value
                except ValueError:
                    variable.lb = value
                    self._lower_bound = value
                    variable.ub = value

        self._upper_bound = value

    @property
    def objective_coefficient(self):
        return self._objective_coefficient

    @objective_coefficient.setter
    def objective_coefficient(self, value):
        model = self.model
        if model is not None:
            model.solver._set_linear_objective_term(self.variable, value)

        self._objective_coefficient = value

    @property
    def effective_lower_bound(self):
        model = self.model
        return \
            flux_analysis.flux_variability_analysis(model, reactions=[self], view=SequentialView(), remove_cycles=False)[
                'lower_bound'][
                self.id]

    @property
    def effective_upper_bound(self):
        model = self.model
        return \
            flux_analysis.flux_variability_analysis(model, reactions=[self], view=SequentialView(), remove_cycles=False)[
                'upper_bound'][
                self.id]

    @property
    def flux(self):
        if self.variable is not None:
            primal = self.variable.primal
            if self.reversibility:
                primal -= self.reverse_variable.primal
            return primal
        else:
            return None

    @property
    def reduced_cost(self):
        if self.variable is not None:
            dual = self.variable.dual
            if dual is None:  # cplex cannot determine reduced costs for MILP problems
                return None
            if self.reversibility:
                dual -= self.reverse_variable.dual
            return dual
        else:
            return None

    def add_metabolites(self, metabolites, **kwargs):
        super(Reaction, self).add_metabolites(metabolites, **kwargs)
        model = self.model
        if model is not None:
            for metabolite, coefficient in metabolites.iteritems():
                model.solver.constraints[metabolite.id] += coefficient*self.variable

    def knock_out(self, time_machine=None):
        def _(reaction, lb, ub):
            # print 'reaction, lb, ub', reaction, lb, ub
            # print 'reaction1, reaction.lower_bound, reaction.upper_bound', reaction, reaction.lower_bound, reaction.upper_bound
            reaction.upper_bound = ub
            reaction.lower_bound = lb
            # print 'reaction2, reaction.lower_bound, reaction.upper_bound', reaction, reaction.lower_bound, reaction.upper_bound
        if time_machine is not None:
            time_machine(do=super(Reaction, self).knock_out, undo=partial(_, self, self.lower_bound, self.upper_bound))
        else:
            super(Reaction, self).knock_out()