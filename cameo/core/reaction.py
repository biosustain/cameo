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

import logging

import cobra
import six

import cameo
from cameo import flux_analysis
from cameo.parallel import SequentialView
from cameo.util import inheritdocstring

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


@six.add_metaclass(inheritdocstring)
class Reaction(cobra.core.Reaction):
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
                    "Can't set attribute %s for reaction %s (while cloning it to a "
                    "cameo style reaction). Skipping it ..." % (
                        attribute, reaction))
        if not isinstance(reaction.model, cameo.core.solver_based_model.SolverBasedModel):
            new_reaction._model = None
        if model is not None:
            new_reaction._model = model
        new_reaction._clone_genes(model)
        new_reaction._clone_metabolites(model)
        return new_reaction

    def _clone_genes(self, model):
        cloned_genes = []
        for gene in self._genes:
            if isinstance(gene, cameo.core.Gene):
                cloned_genes.append(gene)
            else:
                cloned_gene = cameo.core.Gene.clone(gene)
                cloned_genes.append(cloned_gene)
                if model is not None:
                    model.genes._replace_on_id(cloned_gene)
        self._genes = set(cloned_genes)

    def _clone_metabolites(self, model):
        cloned_metabolites = {}
        for metabolite, coeff in self.metabolites.items():
            if isinstance(metabolite, cameo.core.Metabolite):
                cloned_metabolites[metabolite] = coeff
            else:
                cloned_metabolite = cameo.core.Metabolite.clone(metabolite)
                cloned_metabolites[cloned_metabolite] = coeff
                if model is not None:
                    model.metabolites._replace_on_id(cloned_metabolite)
        self._metabolites = cloned_metabolites

    def __init__(self, id=None, name='', subsystem="", lower_bound=0, upper_bound=1000):
        """
        Parameters
        ----------
        name : str, optional
            The name of the reaction.
        """
        super(Reaction, self).__init__(id=id, name=name, subsystem=subsystem,
                                       lower_bound=lower_bound, upper_bound=upper_bound)
        # self._lower_bound = lower_bound
        # self._upper_bound = upper_bound
        # self._model = None
        # self._reverse_variable = None
        # self._forward_variable = None

    def __str__(self):
        return ''.join((self.id, ": ", self.build_reaction_string()))

    # @_cobrapy.core.Reaction.gene_reaction_rule.setter
    # def gene_reaction_rule(self, rule):
    #     _cobrapy.core.Reaction.gene_reaction_rule.fset(self, rule)
    #     self._clone_genes(self.model)

    # @property
    # def reversibility(self):
    #     return self._lower_bound < 0 < self._upper_bound

    # @property
    # def flux_expression(self):
    #     """An optlang variable representing the forward flux (if associated with model), otherwise None.
    #     Representing the net flux if model.reversible_encoding == 'unsplit'"""
    #     model = self.model
    #     if model is not None:
    #         return 1. * self.forward_variable - 1. * self.reverse_variable
    #     else:
    #         return None
    #
    # @property
    # def forward_variable(self):
    #     """An optlang variable representing the forward flux (if associated with model), otherwise None."""
    #     model = self.model
    #     if model is not None:
    #         if self._forward_variable is None:
    #             self._forward_variable = model.solver.variables[self.id]
    #         assert self._forward_variable.problem is self.model.solver
    #
    #         return self._forward_variable
    #     else:
    #         return None
    #
    # @property
    # def reverse_variable(self):
    #     """An optlang variable representing the reverse flux (if associated with model), otherwise None."""
    #     model = self.model
    #     if model is not None:
    #         if self._reverse_variable is None:
    #             self._reverse_variable = model.solver.variables[self.reverse_id]
    #         assert self._reverse_variable.problem is self.model.solver
    #         return self._reverse_variable
    #     else:
    #         return None

    # def __copy__(self):
    #     cop = copy(super(Reaction, self))
    #     cop._reset_var_cache()
    #     return cop
    #
    # def __deepcopy__(self, memo):
    #     cop = deepcopy(super(Reaction, self), memo)
    #     cop._reset_var_cache()
    #     return cop

    # @property
    # def lower_bound(self):
    #     return self._lower_bound

    # @property
    # def functional(self):
    #     """ reaction is functional
    #
    #     Returns
    #     -------
    #     bool
    #         True if the gene-protein-reaction (GPR) rule is fulfilled for this reaction, or if reaction is not
    #         associated to a model, otherwise False.
    #     """
    #     if self._model:
    #         tree, _ = parse_gpr(self.gene_reaction_rule)
    #         return eval_gpr(tree, {gene.id for gene in self.genes if not gene.functional})
    #     return True

    # @lower_bound.setter
    # def lower_bound(self, value):
    #     model = self.model
    #
    #     if model is not None:
    #
    #         forward_variable, reverse_variable = self.forward_variable, self.reverse_variable
    #         if self._lower_bound < 0 < self._upper_bound:  # reversible
    #             if value < 0:
    #                 reverse_variable.ub = -1 * value
    #             elif value >= 0:
    #                 reverse_variable.ub = 0
    #                 try:
    #                     forward_variable.lb = value
    #                 except ValueError:
    #                     forward_variable.ub = value
    #                     self._upper_bound = value
    #                     forward_variable.lb = value
    #         elif self._lower_bound == 0 and self._upper_bound == 0:  # knockout
    #             if value < 0:
    #                 reverse_variable.ub = -1 * value
    #             elif value >= 0:
    #                 forward_variable.ub = value
    #                 self._upper_bound = value
    #                 forward_variable.lb = value
    #         elif self._lower_bound >= 0:  # forward irreversible
    #             if value < 0:
    #                 reverse_variable.ub = -1 * value
    #                 forward_variable.lb = 0
    #             else:
    #                 try:
    #                     forward_variable.lb = value
    #                 except ValueError:
    #                     forward_variable.ub = value
    #                     self._upper_bound = value
    #                     forward_variable.lb = value
    #
    #         elif self._upper_bound <= 0:  # reverse irreversible
    #             if value > 0:
    #                 reverse_variable.lb = 0
    #                 reverse_variable.ub = 0
    #                 forward_variable.ub = value
    #                 self._upper_bound = value
    #                 forward_variable.lb = value
    #             else:
    #                 try:
    #                     reverse_variable.ub = -1 * value
    #                 except ValueError:
    #                     reverse_variable.lb = -1 * value
    #                     self._upper_bound = value
    #                     reverse_variable.ub = -1 * value
    #         else:
    #             raise ValueError('lower_bound issue')
    #
    #     self._lower_bound = value
    #
    # @property
    # def upper_bound(self):
    #     return self._upper_bound
    #
    # @upper_bound.setter
    # def upper_bound(self, value):
    #     model = self.model
    #     if model is not None:
    #
    #         forward_variable, reverse_variable = self.forward_variable, self.reverse_variable
    #         if self._lower_bound < 0 < self._upper_bound:  # reversible
    #             if value > 0:
    #                 forward_variable.ub = value
    #             elif value <= 0:
    #                 forward_variable.ub = 0
    #                 try:
    #                     reverse_variable.lb = -1 * value
    #                 except ValueError:
    #                     reverse_variable.ub = -1 * value
    #                     self._lower_bound = value
    #                     reverse_variable.lb = -1 * value
    #         elif self._lower_bound == 0 and self._upper_bound == 0:  # knockout
    #             if value > 0:
    #                 forward_variable.ub = value
    #             elif value <= 0:
    #                 reverse_variable.ub = -1 * value
    #                 self._lower_bound = value
    #                 reverse_variable.lb = -1 * value
    #         elif self._lower_bound >= 0:  # forward irreversible
    #             if value > 0:
    #                 try:
    #                     forward_variable.ub = value
    #                 except ValueError:
    #                     forward_variable.lb = value
    #                     self._lower_bound = value
    #                     forward_variable.ub = value
    #             else:
    #                 forward_variable.lb = 0
    #                 forward_variable.ub = 0
    #                 reverse_variable.ub = -1 * value
    #                 self._lower_bound = value
    #                 reverse_variable.lb = -1 * value
    #
    #         elif self._upper_bound <= 0:  # reverse irreversible
    #             if value < 0:
    #                 try:
    #                     reverse_variable.lb = -1 * value
    #                 except ValueError:
    #                     reverse_variable.ub = -1 * value
    #                     self._lower_bound = value
    #                     reverse_variable.lb = -1 * value
    #             else:
    #                 forward_variable.ub = value
    #                 reverse_variable.lb = 0
    #         else:
    #             raise ValueError('upper_bound issue')
    #
    #     self._upper_bound = value

    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, value):
        # if value is None:
        #     self.objective_coefficient  # Get objective coefficient from model
        # elif not isinstance(value, cameo.core.SolverBasedModel):
        #     raise ValueError("Must be an instance of cameo.core.SolverBasedModel, not %s" % type(value))
        self._reset_var_cache()
        self._model = value

    # @property
    # def objective_coefficient(self):
    #     if self.model is not None and isinstance(self.model,
    #                                              cameo.core.SolverBasedModel) and self.model.objective is not None:
    #         coefficients_dict = self.model.objective.expression.as_coefficients_dict()
    #         forw_coef = coefficients_dict.get(self.forward_variable, 0)
    #         rev_coef = coefficients_dict.get(self.reverse_variable, 0)
    #         if forw_coef == -rev_coef:
    #             self._objective_coefficient = float(forw_coef)
    #         else:
    #             self._objective_coefficient = 0
    #     return self._objective_coefficient
    #
    # @objective_coefficient.setter
    # def objective_coefficient(self, value):
    #     model = self.model
    #     if model is not None:
    #         coef_difference = value - self.objective_coefficient
    #         model.objective += coef_difference * self.flux_expression
    #     self._objective_coefficient = value

    # @property
    # def effective_lower_bound(self):
    #     model = self.model
    #     fva_result = flux_analysis.flux_variability_analysis(model, reactions=[self], view=SequentialView(),
    #                                                          remove_cycles=False)
    #     return fva_result['lower_bound'][self.id]
    #
    # @property
    # def effective_upper_bound(self):
    #     model = self.model
    #     fva_result = flux_analysis.flux_variability_analysis(model, reactions=[self], view=SequentialView(),
    #                                                          remove_cycles=False)
    #     return fva_result['upper_bound'][self.id]

    # @property
    # def flux(self):
    #     if self.model is not None:
    #         return self.forward_variable.primal - self.reverse_variable.primal
    #     else:
    #         return None
    #
    # @property
    # def reduced_cost(self):
    #     if self.model is not None and self.forward_variable.dual is not None:
    #         return self.forward_variable.dual - self.reverse_variable.dual
    #     else:
    #         return None

    @property
    def is_exchange(self):
        return (len(self.reactants) == 0 or len(self.products) == 0) and len(self.metabolites) == 1

    # def add_metabolites(self, metabolites, combine=True, **kwargs):
    #     if combine:
    #         old_coefficients = self.metabolites
    #     super(Reaction, self).add_metabolites(metabolites, combine=combine, **kwargs)
    #     model = self.model
    #     if model is not None:
    #         for metabolite, coefficient in six.iteritems(metabolites):
    #
    #             if isinstance(metabolite, six.string_types):  # support metabolites added as strings.
    #                 metabolite = model.metabolites.get_by_id(metabolite)
    #             if combine:
    #                 try:
    #                     old_coefficient = old_coefficients[metabolite]
    #                 except KeyError:
    #                     pass
    #                 else:
    #                     coefficient = coefficient + old_coefficient
    #
    #             model.solver.constraints[metabolite.id].set_linear_coefficients({
    #                 self.forward_variable: coefficient,
    #                 self.reverse_variable: -coefficient
    #             })

    # def knock_out(self, time_machine=None):
    #     """Knockout reaction by setting its bounds to zero.
    #
    #     Parameters
    #     ----------
    #     time_machine = TimeMachine
    #         A time TimeMachine instance can be provided to undo the knockout eventually.
    #
    #     Returns
    #     -------
    #     None
    #     """
    #
    #     def _revert(reaction, lb, ub):
    #         """
    #         Set reaction bounds.
    #         """
    #         reaction.upper_bound = ub
    #         reaction.lower_bound = lb
    #
    #     if time_machine is not None:
    #         time_machine(do=super(Reaction, self).knock_out,
    #                      undo=partial(_revert, self, self.lower_bound, self.upper_bound))
    #     else:
    #         super(Reaction, self).knock_out()

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

    # def remove_from_model(self, model=None, remove_orphans=False):
    #     reaction_model = self.model
    #     forward = self.forward_variable
    #     reverse = self.reverse_variable
    #     super(Reaction, self).remove_from_model(model, remove_orphans)
    #     reaction_model.solver.remove([forward, reverse])
    #     self.model = None  # Trigger model setter, since cobrapy only sets _model

    # def delete(self, remove_orphans=False):
    #     model = self.model
    #     forward = self.forward_variable
    #     reverse = self.reverse_variable
    #     super(Reaction, self).delete(remove_orphans)
    #     model.solver.remove([forward, reverse])
    #     self.model = None  # Trigger model setter, since cobrapy only sets _model
    #     # if remove_orphans:
        #     model.solver.remove([metabolite.model.solver for metabolite in self.metabolites.keys()])

    # def change_bounds(self, lb=None, ub=None, time_machine=None):
    #     """Changes one or both of the reaction bounds and allows the changes to be reversed with a TimeMachine"""
    #     if time_machine is not None:
    #         time_machine(do=int,
    #                      undo=partial(setattr, self, "lower_bound", self.lower_bound))
    #         time_machine(do=int,
    #                      undo=partial(setattr, self, "upper_bound", self.upper_bound))
    #     if lb is not None:
    #         self.lower_bound = lb
    #     if ub is not None:
    #         self.upper_bound = ub

    @property
    def n_carbon(self):
        """number of carbon atoms

        Returns
        -------
        int
           number of carbons for all metabolites involved in a reaction
        """
        return sum(metabolite.n_carbon for metabolite in self.metabolites)

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
                <td><strong>GPR</strong></td><td>%s</td>
            </tr>
            <tr>
                <td><strong>Lower bound</strong></td><td>%f</td>
            </tr>
            <tr>
                <td><strong>Upper bound</strong></td><td>%f</td>
            </tr>
        </table>
        """ % (self.id, self.name, self.reaction, self.gene_reaction_rule, self.lower_bound, self.upper_bound)
