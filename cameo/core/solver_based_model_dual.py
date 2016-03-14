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
from sympy import Add

from cameo.core.solver_based_model import SolverBasedModel


def convert_to_dual(model):
    dual_model = model.interface.Model()
    maximization = model.objective.direction == "max"

    if maximization:
        sign = 1
    else:
        sign = -1

    coefficients = {}
    dual_objective = {}

    # Add dual variables from primal constraints:
    for constraint in model.constraints:
        if constraint.expression == 0:
            continue
        if not constraint.is_Linear:
            raise NotImplementedError("Non-linear problems are currently not supported: " + str(constraint))
        if constraint.lb is None and constraint.ub is None:
            continue
        if constraint.lb == constraint.ub:
            const_var = model.interface.Variable("dual_" + constraint.name + "_constraint", lb=-1000, ub=1000)
            dual_model._add_variable(const_var)
            if constraint.lb != 0:
                dual_objective[const_var] = sign * constraint.lb
            for variable, coef in constraint.expression.as_coefficients_dict().items():
                coefficients.setdefault(variable.name, {})[const_var] = sign * coef
        else:
            if constraint.lb is not None:
                lb_var = model.interface.Variable("dual_" + constraint.name + "_constraint_lb", lb=0, ub=1000)
                dual_model._add_variable(lb_var)
                if constraint.lb != 0:
                    dual_objective[lb_var] = -sign * constraint.lb
            if constraint.ub is not None:
                ub_var = model.interface.Variable("dual_" + constraint.name + "_constraint_ub", lb=0, ub=1000)
                dual_model._add_variable(ub_var)
                if constraint.ub != 0:
                    dual_objective[ub_var] = sign * constraint.ub

            for variable, coef in constraint.expression.as_coefficients_dict().items():
                if constraint.lb is not None:
                    coefficients.setdefault(variable.name, {})[lb_var] = -sign * coef
                if constraint.ub is not None:
                    coefficients.setdefault(variable.name, {})[ub_var] = sign * coef

    # Add dual variables from primal bounds
    for variable in model.variables:
        if variable.type != "continuous":
            raise NotImplementedError("Integer variables are currently not supported: " + str(variable))
        if variable.lb is None or variable.lb < 0:
            raise ValueError("Problem is not in standard form (" + variable.name + " can be negative)")
        if variable.lb > 0:
            bound_var = model.interface.Variable("dual_" + variable.name + "_lb", lb=0, ub=1000)
            dual_model._add_variable(bound_var)
            coefficients.setdefault(variable.name, {})[bound_var] = -sign * 1
            dual_objective[bound_var] = -sign * variable.lb
        if variable.ub is not None:
            bound_var = model.interface.Variable("dual_" + variable.name + "_ub", lb=0, ub=1000)
            dual_model._add_variable(bound_var)
            coefficients.setdefault(variable.name, {})[bound_var] = sign * 1
            if variable.ub != 0:
                dual_objective[bound_var] = sign * variable.ub

    # Add dual constraints from primal objective
    primal_objective_dict = model.objective.expression.as_coefficients_dict()
    for variable in model.variables:
        expr = Add(*((coef * dual_var) for dual_var, coef in coefficients[variable.name].items()))
        obj_coef = primal_objective_dict[variable]
        if maximization:
            const = model.interface.Constraint(expr, lb=obj_coef, name="dual_" + variable.name)
        else:
            const = model.interface.Constraint(expr, ub=obj_coef)
        dual_model._add_constraint(const)

    # Make dual objective
    expr = Add(*((coef * dual_var) for dual_var, coef in dual_objective.items() if coef != 0))
    if maximization:
        objective = model.interface.Objective(expr, direction="min")
    else:
        objective = model.interface.Objective(expr, direction="max")
    dual_model.objective = objective

    return dual_model


def to_dual_model(solver_based_model, solver_interface=None):
    if not isinstance(solver_based_model, SolverBasedModel):
        raise TypeError("Input model must be of type 'SolverBasedModel'")
    if solver_interface is None:
        solver_interface = solver_based_model.solver.interface
    return SolverBasedModelDual(solver_based_model, solver_interface=solver_interface)


class SolverBasedModelDual(SolverBasedModel):
    """
    A SolverBasedModel that also contains the dual variables and constraints, allowing primal and dual
    problems to be combined.

    Dual variables corresponding to stoichiometric constraints are prefixed by lambda
    Dual variables corresponding to flux bounds are prefixed by mu
    Other constraints are not supported at the moment!

    Dual constraints will be set according to the original primal objective.
    The objective can be changed subsequently to optimize an outer problem.
    """

    def __init__(self, *args, **kwargs):
        self._dual_variables = {}
        super(SolverBasedModelDual, self).__init__(*args, **kwargs)

    def _populate_solver(self, reactions):
        super(SolverBasedModelDual, self)._populate_solver(reactions)
        metabolites = set.union(*(set(r.metabolites) for r in reactions))
        self._populate_metabolites(metabolites)
        objective_coefficients = self.objective.expression.as_coefficients_dict()
        maximization = self.objective.direction == "max"
        for reaction in reactions:
            forward_coeff = objective_coefficients.get(reaction.forward_variable, 0)
            reverse_coeff = objective_coefficients.get(reaction.reverse_variable, 0)
            self._add_reaction_dual_constraint(reaction, forward_coeff, maximization, "fwd")
            self._add_reaction_dual_constraint(reaction, reverse_coeff, maximization, "rvs")
            self._add_flux_bound_dual_variable(reaction, forward_coeff, maximization, )

    def _add_reaction_dual_constraint(self, reaction, coefficient, maximization, prefix):
        """Add a dual constraint corresponding to the reaction's objective coefficient"""
        stoichiometry = {self.solver.variables["lambda_" + m.id]: c for m, c in six.iteritems(reaction.metabolites)}
        if maximization:
            constraint = self.solver.interface.Constraint(
                Add._from_args(tuple(c * v for v, c in six.iteritems(stoichiometry))),
                name="r_%s_%s" % (reaction.id, prefix),
                lb=coefficient)
        else:
            constraint = self._dual_solver.interface.Constraint(
                Add._from_args(tuple(c * v for v, c in six.iteritems(stoichiometry))),
                name="r_%s_%s" % (reaction.id, prefix),
                ub=coefficient)
        self.solver._add_constraint(constraint)

    @property
    def objective(self):
        return self.solver.objective

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

        self.solver.objective = sum(self.solver.variables["lambda_" + m.id] for m in self.metabolites)
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
            self._add_dual_variable(met.id, "lambda")

    def _add_flux_bound_dual_variable(self, reaction, coefficient, maximization, prefix):
        pass

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

    @property
    def dual_objective(self):
        raise NotImplementedError(
            "This is not yet implemented, but should return an expression (?) that describes the dual objective"
        )

    def primal_objective(self):
        raise NotImplementedError(
            "This is not yet implemented, but should return an expression that describes the primal objective"
        )
