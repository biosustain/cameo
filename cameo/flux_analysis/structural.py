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

from future.utils import implements_iterator

from copy import copy

import sympy
import optlang

from cameo.exceptions import SolveError


@implements_iterator
class ShortestElementaryFluxModes(object):
    def __init__(self, model):
        self._indicator_variables = None
        self._model = model.copy()
        self.__set_exchange_bounds()
        self.__set_up_constraints_and_objective()
        if self.model.solver == optlang.cplex_interface.Model:
            fixed_size_constraint = self.model.solver.interface.Constraint(sympy.Add(*self.indicator_variables), name='fixed_size_constraint', lb=1, ub=1)
            self.model.solver._add_constraint(fixed_size_constraint, sloppy=True)
            self._elementary_mode_generator = self.__generate_elementary_modes_via_fixed_size_constraint()
        else:
            self._elementary_mode_generator = self.__generate_elementary_modes()

    def __set_exchange_bounds(self):
        exchanges = self.model.exchanges
        min_bound = min(exchange.lower_bound for exchange in exchanges)
        max_bound = max(exchange.lower_bound for exchange in exchanges)
        for exchange in exchanges:
            if 0 < exchange.upper_bound < max_bound:
                exchange.upper_bound = max_bound
            if min_bound < exchange.lower_bound < 0:
                exchange.lower_bound = min_bound

    def __set_up_constraints_and_objective(self, c=1):
        indicator_variables = list()
        for reaction in self.model.reactions:
            y_fwd = self.model.solver.interface.Variable('y_fwd_' + reaction.id, type='binary')
            reaction._indicator_variable_fwd = y_fwd
            indicator_constraint_fwd_1 = self.model.solver.interface.Constraint(reaction.forward_variable,
                                                                          indicator_variable=y_fwd, active_when=0, lb=0,
                                                                          ub=0,
                                                                          name='indicator_constraint_fwd_1_{}'.format(reaction.id))

            self.model.solver._add_variable(y_fwd)
            self.model.solver._add_constraint(indicator_constraint_fwd_1, sloppy=True)
            indicator_constraint_fwd_2 = self.model.solver.interface.Constraint(reaction.forward_variable,
                                                                          indicator_variable=y_fwd, active_when=1, lb=c,
                                                                          name='indicator_constraint_fwd_2_{}'.format(reaction.id))
            self.model.solver._add_constraint(indicator_constraint_fwd_2, sloppy=True)


            y_rev = self.model.solver.interface.Variable('y_ref_' + reaction.id, type='binary')
            reaction._indicator_variable_rev = y_rev
            indicator_constraint_rev = self.model.solver.interface.Constraint(reaction.reverse_variable,
                                                                          indicator_variable=y_rev, active_when=0, lb=0,
                                                                          ub=0,
                                                                          name='indicator_constraint_rev_1_{}'.format(reaction.id))

            self.model.solver._add_variable(y_rev)
            self.model.solver._add_constraint(indicator_constraint_rev, sloppy=True)
            indicator_constraint_rev_2 = self.model.solver.interface.Constraint(reaction.reverse_variable,
                                                                          indicator_variable=y_rev, active_when=1, lb=c,
                                                                          name='indicator_constraint_rev_2_{}'.format(reaction.id))
            self.model.solver._add_constraint(indicator_constraint_rev_2, sloppy=True)


            one_direction_constraint = self.model.solver.interface.Constraint(y_fwd + y_rev, lb=0, ub=1, name='one_direction_constraint_{}'.format(reaction.id))
            self.model.solver._add_constraint(one_direction_constraint, sloppy=True)
            indicator_variables.append(y_fwd)
            indicator_variables.append(y_rev)
        at_least_one_active = self.model.solver.interface.Constraint(sympy.Add(*indicator_variables), lb=1, name='an_EM_must_constain_at_least_one_active_reaction')
        self.model.solver._add_constraint(at_least_one_active, sloppy=True)
        self.model.objective = self.model.solver.interface.Objective(sympy.Add(*indicator_variables), direction='min')
        self._indicator_variables = indicator_variables

    def __generate_elementary_modes(self):
        while True:
            try:
                self.model.solve()
            except SolveError as e:
                raise StopIteration
            elementary_flux_mode = list()
            exclusion_list = list()
            for reaction in self.model.reactions:
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
        fixed_size_constraint = self.model.solver.constraints['fixed_size_constraint']
        while True:
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
                    for reaction in self.model.reactions:
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
            if new_fixed_size > len(self.model.reactions):
                break
            fixed_size_constraint.ub, fixed_size_constraint.lb = new_fixed_size, new_fixed_size
            new_fixed_size = fixed_size_constraint.ub + 1
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
        return next(self._elementary_mode_generator)


class MetabolicCutSetsEnumerator(object):
    def __init__(self, model):
        raise NotImplementedError


if __name__ == '__main__':
    from cameo import load_model

    model = load_model('../../tests/data/EcoliCore.xml')
    model.reactions.ATPM.lower_bound, model.reactions.ATPM.upper_bound = 0, 9999999.
    model.solver = 'cplex'
    shortest_emo = ShortestElementaryFluxModes(model)
    # s.model.solver.configuration.verbosity = 3
    count = 0
    for emo in shortest_emo:
        # if count == 1000:
        #     break
        count +=1
        print(str(count)+" "+80*"#")
        print(len(emo))
        for reaction in emo:
            print(reaction, reaction.lower_bound, reaction.upper_bound)