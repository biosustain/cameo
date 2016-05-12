# Copyright 2015 Novo Nordisk Foundation Center for Biosustainability, DTU.
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

"""Methods for stoichiometric analaysis"""

from __future__ import print_function

from copy import copy

import sympy
from cameo.core import SolverBasedModel
from cameo import Reaction
from cobra import Metabolite
import optlang
import six
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

from cameo.exceptions import SolveError, Infeasible
from cameo.util import TimeMachine

__all__ = ['find_dead_end_reactions', 'find_coupled_reactions', 'ShortestElementaryFluxModes']


def find_dead_end_reactions(model):
    """
    Identify reactions that are structurally prevented from carrying flux (dead ends).
    """
    stoichiometries = {}
    for reaction in model.reactions:
        for met, coef in reaction.metabolites.items():
            stoichiometries.setdefault(met.id, {})[reaction.id] = coef

    blocked_reactions = set()
    while True:
        new_blocked = set()
        for met_id, stoichiometry in stoichiometries.items():
            if len(stoichiometry) == 1:  # Metabolite is only associated with 1 reaction, which can thus not be active
                new_blocked.add(list(stoichiometry)[0])
        if len(new_blocked) == 0:
            break  # No more blocked reactions

        # Remove blocked reactions from stoichiometries
        stoichiometries = {
            met_id: {reac_id: coef for reac_id, coef in stoichiometry.items() if reac_id not in new_blocked}
            for met_id, stoichiometry in stoichiometries.items()
        }
        blocked_reactions.update(new_blocked)

    return blocked_reactions


def find_coupled_reactions(model, return_dead_ends=False):
    """Find reaction sets that are structurally forced to carry equal flux"""
    blocked = find_dead_end_reactions(model)
    stoichiometries = {}
    for reaction in model.reactions:
        if reaction.id in blocked:
            continue
        for met, coef in reaction.metabolites.items():
            stoichiometries.setdefault(met.id, {})[reaction.id] = coef

    # Find reaction pairs that are constrained to carry equal flux
    couples = []
    for met_id, stoichiometry in stoichiometries.items():
        if len(stoichiometry) == 2 and set(stoichiometry.values()) == {1, -1}:
            couples.append(set(stoichiometry.keys()))

    # Aggregate the pairs into groups
    coupled_groups = []
    for couple in couples:
        for group in coupled_groups:
            if len(couple & group) != 0:
                group.update(couple)
                break
        else:
            coupled_groups.append(couple)

    if return_dead_ends:
        return coupled_groups, blocked
    else:
        return coupled_groups


class ShortestElementaryFluxModes(six.Iterator):
    def __init__(self, model, reactions=None, c=1e-5, copy=True, change_bounds=True):
        self._indicator_variables = None
        if copy:
            self._model = model.copy()
        else:
            self._model = model
        if reactions is None:
            self._reactions = self.model.reactions
        else:
            self._reactions = []
            for reaction in reactions:
                if isinstance(reaction, six.string_types):
                    self._reactions.append(self.model.reactions.get_by_id(reaction))
                else:
                    self._reactions.append(reaction)
        if change_bounds:
            self.__set_exchange_bounds()
        self.__set_up_constraints_and_objective(c)
        if type(self.model.solver) == optlang.cplex_interface.Model:
            self._elementary_mode_generator = self.__generate_elementary_modes_via_fixed_size_constraint()
        else:
            self._elementary_mode_generator = self.__generate_elementary_modes()

    def __set_exchange_bounds(self):
        exchanges = self.model.exchanges
        min_bound = min(exchange.lower_bound for exchange in exchanges)
        max_bound = max(exchange.upper_bound for exchange in exchanges)
        for exchange in exchanges:
            if 0 < exchange.upper_bound < max_bound:
                exchange.upper_bound = max_bound
            if min_bound < exchange.lower_bound < 0:
                exchange.lower_bound = min_bound

    def __set_up_constraints_and_objective(self, c=1):
        indicator_variables = list()
        for reaction in self._reactions:
            y_fwd = self.model.solver.interface.Variable('y_fwd_' + reaction.id, type='binary')
            reaction._indicator_variable_fwd = y_fwd
            indicator_constraint_fwd_1 = self.model.solver.interface.Constraint(
                reaction.forward_variable,
                indicator_variable=y_fwd, active_when=0, lb=0,
                ub=0,
                name='indicator_constraint_fwd_1_{}'.format(reaction.id)
            )

            self.model.solver._add_variable(y_fwd)
            self.model.solver._add_constraint(indicator_constraint_fwd_1, sloppy=True)
            indicator_constraint_fwd_2 = self.model.solver.interface.Constraint(
                reaction.forward_variable,
                indicator_variable=y_fwd, active_when=1, lb=c,
                name='indicator_constraint_fwd_2_{}'.format(reaction.id)
            )
            self.model.solver._add_constraint(indicator_constraint_fwd_2, sloppy=True)

            y_rev = self.model.solver.interface.Variable('y_rev_' + reaction.id, type='binary')
            reaction._indicator_variable_rev = y_rev
            indicator_constraint_rev = self.model.solver.interface.Constraint(
                reaction.reverse_variable,
                indicator_variable=y_rev, active_when=0, lb=0,
                ub=0,
                name='indicator_constraint_rev_1_{}'.format(reaction.id)
            )

            self.model.solver._add_variable(y_rev)
            self.model.solver._add_constraint(indicator_constraint_rev, sloppy=True)
            indicator_constraint_rev_2 = self.model.solver.interface.Constraint(
                reaction.reverse_variable,
                indicator_variable=y_rev, active_when=1, lb=c,
                name='indicator_constraint_rev_2_{}'.format(reaction.id)
            )
            self.model.solver._add_constraint(indicator_constraint_rev_2, sloppy=True)

            one_direction_constraint = self.model.solver.interface.Constraint(
                y_fwd + y_rev, lb=0, ub=1, name='one_direction_constraint_{}'.format(reaction.id)
            )
            self.model.solver._add_constraint(one_direction_constraint, sloppy=True)
            indicator_variables.append(y_fwd)
            indicator_variables.append(y_rev)
        at_least_one_active = self.model.solver.interface.Constraint(
            sympy.Add(*indicator_variables), lb=1, name='an_EM_must_constain_at_least_one_active_reaction'
        )
        self.model.solver._add_constraint(at_least_one_active, sloppy=True)
        self.model.objective = self.model.solver.interface.Objective(sympy.Add(*indicator_variables), direction='min')
        self._indicator_variables = indicator_variables

    def __generate_elementary_modes(self):
        while True:
            try:
                self.model.solve()
            except SolveError:
                raise StopIteration
            elementary_flux_mode = list()
            exclusion_list = list()
            for reaction in self._reactions:
                if reaction._indicator_variable_fwd.primal == 1.:
                    reaction_copy = copy(reaction)
                    reaction_copy.lower_bound = 0
                    elementary_flux_mode.append(reaction_copy)
                    exclusion_list.append(reaction._indicator_variable_fwd)
                elif reaction._indicator_variable_rev.primal == 1.:
                    reaction_copy = copy(reaction)
                    reaction_copy.upper_bound = 0
                    elementary_flux_mode.append(reaction_copy)
                    exclusion_list.append(reaction._indicator_variable_rev)
            exclusion_constraint = self.model.solver.interface.Constraint(sympy.Add(*exclusion_list), ub=len(exclusion_list) - 1)
            self.model.solver._add_constraint(exclusion_constraint, sloppy=True)
            yield elementary_flux_mode

    def __generate_elementary_modes_via_fixed_size_constraint(self):
        fixed_size_constraint = self.model.solver.interface.Constraint(
            sympy.Add(*self.indicator_variables), name='fixed_size_constraint', lb=1, ub=1
        )
        self.model.solver._add_constraint(fixed_size_constraint, sloppy=True)

        while True:
            logger.debug("Looking for solutions with cardinality " + str(fixed_size_constraint.lb))
            try:
                self.model.solver.problem.populate_solution_pool()
            except Exception as e:
                raise e
            solution_num = self.model.solver.problem.solution.pool.get_num()
            if solution_num > 0:
                exclusion_lists = list()
                for i in range(solution_num):
                    elementary_flux_mode = list()
                    exclusion_list = list()
                    for reaction in self._reactions:
                        if self.model.solver.problem.solution.pool.get_values(i, reaction._indicator_variable_fwd.name) == 1.:
                            reaction_copy = copy(reaction)
                            reaction_copy.lower_bound = 0
                            elementary_flux_mode.append(reaction_copy)
                            exclusion_list.append(reaction._indicator_variable_fwd)
                        elif self.model.solver.problem.solution.pool.get_values(i, reaction._indicator_variable_rev.name) == 1.:
                            reaction_copy = copy(reaction)
                            reaction_copy.upper_bound = 0
                            elementary_flux_mode.append(reaction_copy)
                            exclusion_list.append(reaction._indicator_variable_rev)
                    exclusion_lists.append(exclusion_list)
                    yield elementary_flux_mode
                for exclusion_list in exclusion_lists:
                    exclusion_constraint = self.model.solver.interface.Constraint(sympy.Add(*exclusion_list), ub=len(exclusion_list) - 1)
                    self.model.solver._add_constraint(exclusion_constraint, sloppy=True)
            new_fixed_size = fixed_size_constraint.ub + 1
            if new_fixed_size > len(self._reactions):
                break
            fixed_size_constraint.ub, fixed_size_constraint.lb = new_fixed_size, new_fixed_size


    @property
    def model(self):
        return self._model

    @property
    def indicator_variables(self):
        return self._indicator_variables

    def __iter__(self):
        return self

    def __next__(self):
        return six.next(self._elementary_mode_generator)


class MinimalCutSetsEnumerator(ShortestElementaryFluxModes):
    def __init__(self, model, targets, constraints=None, exclude=None, c=1e-5):
        if exclude is not None:
            exclude_copy = []
            for re in exclude:
                if isinstance(re, six.string_types):
                    exclude_copy.append(re)
                else:
                    exclude_copy.append(re.id)
        else:
            exclude_copy = ()
        exclude = exclude_copy

        self._primal_model = model.copy()
        self._constraints = self._construct_constraints(constraints)
        self._target_constraints = self._convert_target_to_constraints(targets)
        self._dual_model = self._make_dual_model(model)
        self._add_target_constraints(c)

        super(MinimalCutSetsEnumerator, self).__init__(
            self._dual_model, copy=False, change_bounds=False, reactions=self._v_reactions, c=c)

        def iterator_wrapper(generator, func):
            """Convert MCS's to primal and filter according to constraints"""
            for mcs in generator:
                allowed_mcs = func(mcs)
                if allowed_mcs is not None:
                    yield allowed_mcs

        self._elementary_mode_generator = iterator_wrapper(self._elementary_mode_generator, self._allowed_mcs)

    def _allowed_mcs(self, dual_em):
        mcs = self._convert_mcs_to_primal(dual_em)
        if len(self._constraints) > 0:
            with TimeMachine() as tm:
                for reac_id in mcs:
                    self._primal_model.reactions.get_by_id(reac_id).knock_out(tm)
                try:
                    self._primal_model.solve()
                except Infeasible:
                    return None
                else:
                    return mcs
        else:
            return mcs

    def _convert_mcs_to_primal(self, dual_em):
        primal_mcs = []
        for reac in dual_em:
            name = "_".join(reac.id.split("_")[1:])
            if name in self._primal_model.reactions:
                primal_mcs.append(name)
            elif "_reverse_" in name and "_".join(name.split("_")[:-2]) in self._primal_model.reactions:
                primal_mcs.append("_".join(name.split("_")[:-2]))
            else:
                raise RuntimeError("Primal reaction could not be found for " + reac.id)
        return set(primal_mcs)

    def _make_dual_model(self, model):
        dual_model = SolverBasedModel(solver_interface=model.solver.interface)

        # Add dual metabolites
        dual_metabolite_names = []
        irreversibles = []  # Metabolites for which a z-reaction must be added
        self._dual_to_primal_mapping = {}
        for re in model.reactions:
            forward_id = re._get_forward_id()
            reverse_id = re._get_reverse_id()
            if forward_id in self._split_vars or reverse_id in self._split_vars:
                if re.upper_bound > 0 or forward_id in self._split_vars:
                    dual_metabolite_names.append(forward_id)
                    irreversibles.append(forward_id)
                    self._dual_to_primal_mapping[forward_id] = re.id
                if re.lower_bound < 0 or reverse_id in self._split_vars:
                    dual_metabolite_names.append(reverse_id)
                    irreversibles.append(reverse_id)
                    self._dual_to_primal_mapping[reverse_id] = re.id
            else:
                dual_metabolite_names.append(re.id)
                if re.lower_bound >= 0:
                    irreversibles.append(re.id)
                self._dual_to_primal_mapping[re.id] = re.id
        dual_model.add_metabolites([Metabolite(name) for name in dual_metabolite_names])

        # Add dual "u-reactions"
        transposed_stoichiometry = {}
        for reaction in model.reactions:
            for met, coef in reaction.metabolites.items():
                if reaction._get_reverse_id() in dual_metabolite_names:
                    transposed_stoichiometry.setdefault(met, {})[
                        dual_model.metabolites.get_by_id(reaction._get_reverse_id())] = -coef
                if reaction._get_forward_id() in dual_metabolite_names:
                    transposed_stoichiometry.setdefault(met, {})[
                        dual_model.metabolites.get_by_id(reaction._get_forward_id())] = coef
                elif reaction.id in dual_metabolite_names:  # This should be the same as forward_var.name but in general it might not be
                    transposed_stoichiometry.setdefault(met, {})[
                        dual_model.metabolites.get_by_id(reaction.id)] = coef

        u_reactions = []
        for met, stoichiometry in transposed_stoichiometry.items():
            met_id = met.id
            reac = Reaction("u_" + met_id)
            reac.lower_bound = -reac.upper_bound  # Make reversible
            reac.add_metabolites(stoichiometry)
            u_reactions.append(reac)

        dual_model.add_reactions(u_reactions)

        # Add dual "v-reactions"
        v_reactions = []
        for dual_met in dual_model.metabolites:
            reac = Reaction("v_" + dual_met.id)
            reac.lower_bound = -reac.upper_bound  # Make reversible
            reac.add_metabolites({dual_met: 1})
            v_reactions.append(reac)
        dual_model.add_reactions(v_reactions)
        self._v_reactions = v_reactions

        # Add dual "z-reactions"
        z_reactions = []
        for dual_met in dual_model.metabolites:
            if dual_met.id in irreversibles:
                reac = Reaction("z_" + dual_met.id)
                reac.lower_bound = 0
                reac.add_metabolites({dual_met: -1})
                z_reactions.append(reac)
        dual_model.add_reactions(z_reactions)

        return dual_model

    def _convert_target_to_constraints(self, targets):
        """Make a list of constraints that describe the target polytope."""

        def _convert_string_to_target_constraints(target_string):
            """Convert a string to constraints, e.g. a reaction name to constraints that target that reaction"""
            if target_string in self._primal_model.reactions:
                raise NotImplementedError  # TODO
            else:
                raise NotImplementedError
            return NotImplemented  # return constraints

        if isinstance(targets, optlang.interface.Constraint):
            targets = [targets]
        elif isinstance(targets, six.string_types):
            targets = [targets]
        if isinstance(targets[0], six.string_types):
            new_targets = []
            for target in targets:
                new_targets.extend(_convert_string_to_target_constraints(target))
            targets = new_targets

        # Find the variables that are used to describe the target polytope (these have to be split).
        split_vars = set()
        for target in targets:
            split_vars.update(set(var.name for var in target.variables))
        self._split_vars = split_vars

        return targets

    def _add_target_constraints(self, c):
        """Add the target constraints to the dual model"""
        targets = self._target_constraints

        w_reactions = []
        for i, target in enumerate(targets):
            assert isinstance(target, optlang.interface.Constraint)
            if not target.is_Linear:
                raise ValueError("Target constraints must be linear.")
            if (target.lb is None and target.ub is None) or (target.lb is not None and target.ub is not None):
                raise ValueError("Target constraints must be one-sided inequalities.")
            coefficients_dict = target.expression.as_coefficients_dict()
            w_reac = Reaction("w_" + str(i))
            w_reac.lower_bound = c
            if target.ub is not None:
                coefficients = {
                    self._dual_model.metabolites.get_by_id(var.name): coef for var, coef in coefficients_dict.items()
                    if var.name in self._dual_model.metabolites
                }
            elif target.lb is not None:
                coefficients = {
                    self._dual_model.metabolites.get_by_id(var.name): -coef for var, coef in coefficients_dict.items()
                    if var.name in self._dual_model.metabolites
                }
            w_reac.add_metabolites(coefficients)
            w_reactions.append(w_reac)
        self._dual_model.add_reactions(w_reactions)

        return None

    def _construct_constraints(self, constraints):
        if constraints is None:
            self._illegal_knockouts = set()
            return ()
        else:
            cloned_constraints = [
                self._primal_model.solver.interface.Constraint.clone(constraint, model=self._primal_model.solver)
                for constraint in constraints
            ]
            self._primal_model.solver.add(cloned_constraints)
            illegal_knockouts = []
            for reaction in self._primal_model.reactions:
                # If single knockout causes the constrained model to become infeasible, then no superset
                # of knockouts can be feasible either.
                with TimeMachine() as tm:
                    reaction.knock_out(tm)
                    try:
                        self._primal_model.solve()
                    except Infeasible:
                        illegal_knockouts.append(reaction.id)
            self._illegal_knockouts = illegal_knockouts
            return cloned_constraints


if __name__ == '__main__':
    from cameo import load_model

    model = load_model("e_coli_core")
    # model = load_model('../../tests/data/EcoliCore.xml')
    model.reactions.ATPM.lower_bound, model.reactions.ATPM.upper_bound = 0, 1000.
    model.reactions.BIOMASS_Ecoli_core_w_GAM.lower_bound = 1
    model.solver = 'cplex'
    shortest_emo = ShortestElementaryFluxModes(model)
    # s.model.solver.configuration.verbosity = 3
    count = 0
    for emo in shortest_emo:
        if count == 1000:
            break
        count += 1
        print(str(count) + " " + 80 * "#")
        print(len(emo))
        for reaction in emo:
            print(reaction, reaction.lower_bound, reaction.upper_bound)
