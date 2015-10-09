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

from __future__ import print_function

from cameo.util import TimeMachine
from cobra import Metabolite
from cameo import Reaction
from cameo.core.solver_based_model_dual import SolverBasedModel, to_dual_model, convert_to_dual
from cameo.strain_design import StrainDesignMethod, StrainDesignResult
from sympy import Add, Mul
import pandas as pd
import logging

logger = logging.getLogger(__name__)


# TODO: Implement integer cuts to return multiple knockout strategies
class OptKnock(StrainDesignMethod):
    def __init__(self, model, exclude_reactions=None, *args, **kwargs):
        super(OptKnock, self).__init__(*args, **kwargs)
        self._model = model.copy()

        self._build_problem(exclude_reactions)

    def _make_dual(self):
        dual_problem = convert_to_dual(self._model.solver)
        self._dual_problem = dual_problem
        logger.debug("Dual problem successfully created")

    def _build_problem(self, essential_reactions):
        logger.debug("Starting to formulate OptKnock problem")

        self.essential_reactions = self._model.essential_reactions() + self._model.exchanges
        if essential_reactions:
            self.essential_reactions += essential_reactions

        self._make_dual()

        self._combine_primal_and_dual()
        logger.debug("Primal and dual successfully combined")

        y_vars = {}
        constrained_dual_vars = set()
        for reaction in self._model.reactions:
            if reaction not in self.essential_reactions and reaction.lower_bound <= 0 <= reaction.upper_bound:
                y_var, constrained_vars = self._add_knockout_constraints(reaction)
                y_vars[y_var] = reaction
                constrained_dual_vars.update(constrained_vars)
        self._y_vars = y_vars

        primal_objective = self._model.solver.objective
        dual_objective = self._model.solver.interface.Objective.clone(
            self._dual_problem.objective, model=self._model.solver)

        reduced_expression = Add(*((c * v) for v, c in dual_objective.expression.as_coefficients_dict().items()
                                   if v not in constrained_dual_vars))
        dual_objective = self._model.solver.interface.Objective(reduced_expression, direction=dual_objective.direction)

        optimality_constraint = self._model.solver.interface.Constraint(
            primal_objective.expression - dual_objective.expression,
            lb=0, ub=0, name="inner_optimality")
        self._model.solver.add(optimality_constraint)
        logger.debug("Inner optimality constrained")

        logger.debug("Adding constraint for number of knockouts")
        knockout_number_constraint = self._model.solver.interface.Constraint(
            Add(*y_vars), lb=len(y_vars), ub=len(y_vars)
        )
        self._model.solver.add(knockout_number_constraint)
        self._number_of_knockouts_constraint = knockout_number_constraint

    def _combine_primal_and_dual(self):
        primal_problem = self._model.solver
        dual_problem = self._dual_problem

        for var in dual_problem.variables:
            var = primal_problem.interface.Variable.clone(var)
            primal_problem.add(var)
        for const in dual_problem.constraints:
            const = primal_problem.interface.Constraint.clone(const, model=primal_problem)
            primal_problem.add(const)

    def _add_knockout_constraints(self, reaction):
        interface = self._model.solver.interface
        y_var = interface.Variable("y_"+reaction.id, type="binary")

        self._model.solver.add(interface.Constraint(reaction.flux_expression-1000*y_var, ub=0))
        self._model.solver.add(interface.Constraint(reaction.flux_expression+1000*y_var, lb=0))

        constrained_vars = []

        if reaction.upper_bound != 0:
            dual_forward_ub = self._model.solver.variables["dual_"+reaction.forward_variable.name+"_ub"]
            self._model.solver.add(interface.Constraint(dual_forward_ub-1000*(1-y_var), ub=0))
            constrained_vars.append(dual_forward_ub)
        if reaction.lower_bound != 0:
            dual_reverse_ub = self._model.solver.variables["dual_"+reaction.reverse_variable.name+"_ub"]
            self._model.solver.add(interface.Constraint(dual_reverse_ub - 1000*(1-y_var), ub=0))
            constrained_vars.append(dual_reverse_ub)

        return y_var, constrained_vars

    def run(self, k, target, *args, **kwargs):
        """
        Perform the OptKnock simulation
        :param k: The maximal allowed number of knockouts
        :param target: The reaction to be optimized
        :return: KnockoutResult
        """
        self._model.objective = target
        self._number_of_knockouts_constraint.lb = self._number_of_knockouts_constraint.ub - k
        solution = self._model.solve()

        knockouts = set(reac for y, reac in self._y_vars.items() if round(y.primal, 3) == 0)
        fluxes = solution.fluxes
        production = solution.f

        return KnockoutResult(knockouts, fluxes, production, target)


class RobustKnock(StrainDesignMethod):
    pass


class KnockoutResult(StrainDesignResult):
    def __init__(self, knockouts, fluxes, production, target, *args, **kwargs):
        super(KnockoutResult, self).__init__(*args, **kwargs)
        self._knockouts = knockouts
        self._fluxes = fluxes
        self._production = production
        self._target = target

    @property
    def knockouts(self):
        return self._knockouts

    @property
    def fluxes(self):
        return self._fluxes

    @property
    def production(self):
        return self._production

    @property
    def target(self):
        return self._target

    @property
    def data_frame(self):
        return pd.DataFrame((self. knockouts, self.production), names=["Knockouts", "Production"])

    def _repr_html_(self):
        html_string = """
<b>Target:</b> %s</br>
%s""" % (self._production, self.data_frame._repr_html_())
        return html_string
