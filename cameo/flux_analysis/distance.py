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

from __future__ import absolute_import, print_function

from types import DictType
from sympy import Add
from cameo.core.result import FluxDistributionResult
from cameo.core.solution import SolutionBase, Solution
import six

add = Add._from_args


class ManhattanDistance(object):

    """Compute steady-state fluxes that maximize/minimize the Manhattan distance (L1 norm)
    to a reference flux distribution.

    Parameters
    ----------
    model : Model
    reference : dict or solution
        A

    Attributes
    ----------
    model : Model
    reference : dict
    """

    @staticmethod
    def __check_valid_reference(reference):
        if not isinstance(reference, DictType) and not isinstance(reference, SolutionBase):
            raise ValueError('%s is not a valid reference flux distribution. Needs to be either a dict or Solution object.')

    def __init__(self, model, reference=None, *args, **kwargs):
        super(ManhattanDistance, self).__init__(*args, **kwargs)
        self._aux_variables = dict()
        self._deviation_constraints = dict()
        self.__check_valid_reference(reference)
        self._reference = reference
        self._model = model.copy()
        self.__prep_model()

    @property  # read-only
    def model(self):
        return self._model

    @property
    def reference(self):
        return self._reference

    @reference.setter
    def reference(self, value):
        self.__check_valid_reference(value)
        self.__set_new_reference(value)

    def __prep_model(self):
        for rid, flux_value in six.iteritems(self.reference):
            self.__add_deviavtion_constraint(rid, flux_value)
        objective = self.model.solver.interface.Objective(add(list(self._aux_variables.values())), name='deviations')
        self.model.objective = objective

    def __set_new_reference(self, reference):
        # remove unnecessary constraints
        constraints_to_remove = list()
        aux_vars_to_remove = list()
        for key in list(self._deviation_constraints.keys()):
            if key not in reference:
                constraints_to_remove.extend(self._deviation_constraints.pop(key))
                aux_vars_to_remove.append(self._aux_variables[key])
        self.model.solver._remove_constraints(constraints_to_remove)
        self.model.solver._remove_variables(aux_vars_to_remove)
        # Add new or adapt existing constraints
        for key, value in six.iteritems(reference):
            try:
                (lb_constraint, ub_constraint) = self._deviation_constraints[key]
                lb_constraint.lb = value
                ub_constraint.ub = value
            except KeyError:
                self.__add_deviavtion_constraint(key, value)

    def __add_deviavtion_constraint(self, reaction_id, flux_value):
        reaction = self.model.reactions.get_by_id(reaction_id)
        aux_var = self.model.solver.interface.Variable('aux_'+reaction_id, lb=0)
        self._aux_variables[reaction_id] = aux_var
        self.model.solver._add_variable(aux_var)
        if reaction.reverse_variable is None:
            expression = reaction.variable - aux_var
        else:
            expression = reaction.variable - reaction.reverse_variable - aux_var
        constraint_lb = self.model.solver.interface.Constraint(expression, ub=flux_value, name='deviation_lb_'+reaction_id)
        self.model.solver._add_constraint(constraint_lb, sloppy=True)
        if reaction.reverse_variable is None:
            expression = reaction.variable + aux_var
        else:
            expression = reaction.variable - reaction.reverse_variable + aux_var
        constraint_ub = self.model.solver.interface.Constraint(expression, lb=flux_value, name='deviation_ub_'+reaction_id)
        self.model.solver._add_constraint(constraint_ub, sloppy=True)
        self._deviation_constraints[reaction_id] = (constraint_lb, constraint_ub)

    def minimize_L1(self, *args, **kwargs):
        self.model.objective.direction = 'min'
        solution = self.model.solve()
        result = FluxDistributionResult(solution)
        return result

    # def maximize_L1(self, *args, **kwargs):
    #     self.model.objective.direction = 'max'
    #     return self.model.solve(solution_type=Solution)

