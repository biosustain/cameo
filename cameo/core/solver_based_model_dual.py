# Copyright 2015 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import six
from cameo.core.solver_based_model import SolverBasedModel


def to_dual_model(solver_based_model, solver_interface=None):
    assert isinstance(solver_based_model, solver_based_model)
    solver_interface = solver_interface if solver_interface is not None else solver_based_model.solver.interface
    SolverBasedModelDual(solver_based_model, solver_interface=solver_interface)


class SolverBasedModelDual(SolverBasedModel):
    def __init__(self, *args, **kwargs):
        super(SolverBasedModelDual, self).__init__(*args, **kwargs)

    def _populate_solver(self, reactions):
        metabolites = reduce(set.union, {r.metabolites for r in reactions})
        self._populate_metabolites(metabolites)
        objective_coefficients = self.objective.expression.as_coefficients_dict()
        maximization = self.objective.direction == "max"
        for reaction in reactions:
            forward_coeff = objective_coefficients.get(reaction.forward_variable, 0)
            reverse_coeff = objective_coefficients.get(reaction.reverse_variable, 0)
            self._add_reaction_dual_constraint(reaction, forward_coeff, maximization, "fwd")
            self._add_reaction_dual_constraint(reaction, reverse_coeff, maximization, "rvs")

    def _add_reaction_dual_constraint(self, reaction, coefficient, maximization, prefix):
        stoichiometry = {self.solver.variables[m.id + "_u"]: c for m, c in six.iteritems(reaction.metabolites)}
        if maximization:
            constraint = self._dual_solver.interface.Constraint(sum([c * v for v, c in six.iteritems(stoichiometry)]),
                                                                name="r_%s_%s" % (reaction.id, prefix),
                                                                lb=coefficient)
        else:
            constraint = self._dual_solver.interface.Constraint(sum([c * v for v, c in six.iteritems(stoichiometry)]),
                                                                name="r_%s_%s" % (reaction.id, prefix),
                                                                ub=coefficient)
        self._dual_solver._add_constraint(constraint)

    @property
    def objective(self):
        return self._solver.objective

    @objective.setter
    def objective(self, objective):
        objective = self.solver.interface.Objective(objective)
        objective_coefficients = objective.expression.as_coefficients_dict()
        maximization = objective.direction == "max"
        for reaction in self.reactions:
            forward_coeff = objective_coefficients.get(reaction.forward_variable, 0)
            reverse_coeff = objective_coefficients.get(reaction.reverse_variable, 0)
            self._update_dual_reaction_constraint(reaction, forward_coeff, maximization, "fwd")
            self._update_dual_reaction_constraint(reaction, reverse_coeff, maximization, "rvs")

        self.solver.objective = sum([self.solver.variables[m.id + "_u"] for m in self.metabolites])
        self.objective.direction = "min" if maximization else "max"

    def _update_dual_reaction_constraint(self, reaction, coefficient, maximization, prefix):
        constraint = self.solver.constraints["r_%s_%s" % (reaction.id, prefix)]
        if coefficient == 0:
            constraint.lb = None
            constraint.ub = None
        else:
            if maximization:
                constraint.lb = coefficient
                constraint.ub = None
            else:
                constraint.lb = None
                constraint.ub = coefficient

    def _populate_metabolites(self, metabolites):
        for met in metabolites:
            self._add_dual_variable(met.id, "_u")

    def _add_dual_variable(self, identifier, prefix):
        dual_variable_id = prefix + "_" + identifier
        dual_variable = self.solver.interface.Variable(dual_variable_id)
        self._dual_variables[dual_variable_id] = dual_variable
        self.solver._add_variable(dual_variable)

    @property
    def objective(self):
        return self.solver.objective

    @objective.setter
    def objective(self, objective):
        self.solver.objective = objective
        self._update_constraints()