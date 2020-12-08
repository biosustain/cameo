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

import logging
from copy import copy
from itertools import product

import numpy as np
import optlang
from optlang.interface import OPTIMAL
import pandas
import sympy
from cobra import Metabolite, Reaction, Model
from numpy.linalg import svd
from scipy.sparse import dok_matrix, lil_matrix

from cobra.exceptions import OptimizationError

__all__ = ['find_dead_end_reactions', 'find_coupled_reactions', 'ShortestElementaryFluxModes']

logger = logging.getLogger(__name__)


def create_stoichiometric_array(model, array_type='dense', dtype=None):
    """Return a stoichiometric array representation of the given model.
    The the columns represent the reactions and rows represent
    metabolites. S[i,j] therefore contains the quantity of metabolite `i`
    produced (negative for consumed) by reaction `j`.
    Parameters
    ----------
    model : cobra.Model
        The cobra model to construct the matrix for.
    array_type : string
        The type of array to construct. if 'dense', return a standard
        numpy.array, 'dok', or 'lil' will construct a sparse array using
        scipy of the corresponding type and 'data_frame' will give a
        pandas `DataFrame` with metabolite and reaction identifiers as indices.
    dtype : data-type
        The desired data-type for the array. If not given, defaults to float.
    Returns
    -------
    matrix of class `dtype`
        The stoichiometric matrix for the given model.
    """

    if dtype is None:
        dtype = np.float64

    def data_frame(_, dtype):
        metabolite_ids = [met.id for met in model.metabolites]
        reaction_ids = [rxn.id for rxn in model.reactions]
        index = pandas.MultiIndex.from_tuples(
            list(product(metabolite_ids, reaction_ids)))
        return pandas.DataFrame(data=0, index=index, columns=['stoichiometry'], dtype=dtype)

    array_constructor = {
        'dense': np.zeros, 'dok': dok_matrix, 'lil': lil_matrix,
        'data_frame': data_frame
    }

    n_metabolites = len(model.metabolites)
    n_reactions = len(model.reactions)
    array = array_constructor[array_type]((n_metabolites, n_reactions), dtype=dtype)

    m_ind = model.metabolites.index
    r_ind = model.reactions.index

    for reaction in model.reactions:
        for metabolite, stoich in reaction.metabolites.items():
            if array_type == 'data_frame':
                array.set_value((metabolite.id, reaction.id), 'stoichiometry', stoich)
            else:
                array[m_ind(metabolite), r_ind(reaction)] = stoich

    return array


# Taken from http://wiki.scipy.org/Cookbook/RankNullspace
def nullspace(matrix, atol=1e-13, rtol=0):
    """Compute an approximate basis for the nullspace of A.

    The algorithm used by this function is based on the singular value
    decomposition of `A`.

    Parameters
    ----------
    matrix : ndarray
        A should be at most 2-D.  A 1-D array with length k will be treated
        as a 2-D with shape (1, k)
    atol : float
        The absolute tolerance for a zero singular value.  Singular values
        smaller than `atol` are considered to be zero.
    rtol : float
        The relative tolerance.  Singular values less than rtol*smax are
        considered to be zero, where smax is the largest singular value.

    If both `atol` and `rtol` are positive, the combined tolerance is the
    maximum of the two; that is::
        tol = max(atol, rtol * smax)
    Singular values smaller than `tol` are considered to be zero.

    Return value
    ------------
    ns : ndarray
        If `A` is an array with shape (m, k), then `ns` will be an array
        with shape (k, n), where n is the estimated dimension of the
        nullspace of `A`.  The columns of `ns` are a basis for the
        nullspace; each element in numpy.dot(A, ns) will be approximately
        zero.
    """
    matrix = np.atleast_2d(matrix)
    u, s, vh = svd(matrix)
    tol = max(atol, rtol * s[0])
    nnz = (s >= tol).sum()
    ns = vh[nnz:].conj().T
    return ns


def find_blocked_reactions_nullspace(model, ns=None, tol=1e-10):
    """Identify reactions that can't carry flux based on the nullspace of the stoichiometric matrix.
    A blocked reaction will correspond to an all-zero row in N(S)."""
    if ns is None:
        ns = nullspace(create_stoichiometric_array(model))
    mask = (np.abs(ns) <= tol).all(1)
    blocked = frozenset(reac for reac, b in zip(model.reactions, mask) if b is True)
    return blocked


def find_coupled_reactions_nullspace(model, ns=None, tol=1e-10):
    """
    Find groups of reactions whose fluxes are forced to be multiples of each other.

    Parameters
    ----------
    model: SolverBasedModel
        A constraint-based model.
    ns: numpy.array
        The nullspace if already computed somewhere else.
    tol: float
        Value after which a value is treated as zero.

    Returns
    -------
    list:
        A list with groups: dictionaries {reaction: relative_coefficient}

    """
    if ns is None:
        ns = nullspace(create_stoichiometric_array(model))
    mask = (np.abs(ns) <= tol).all(axis=1)  # Mask for blocked reactions
    non_blocked_ns = ns[~mask]
    non_blocked_reactions = np.array(list(model.reactions))[~mask]

    corr_mat = np.corrcoef(non_blocked_ns)
    dist_mat = 1 - np.abs(corr_mat)
    groups = []

    reaction_index = {r: i for i, r in enumerate(non_blocked_reactions)}

    for i, reaction_i in enumerate(non_blocked_reactions):
        left = non_blocked_ns[i]
        group = next((g for g in groups if reaction_i in g), None)
        if group:
            reaction_i = next(reac for reac, c in group.items() if c == 1)
            left = non_blocked_ns[reaction_index[reaction_i]]
        else:
            group = {reaction_i: 1}
            groups.append(group)

        good_corr_index = np.argwhere(dist_mat[i] < tol).reshape(-1)
        for j in good_corr_index:
            if j == i:
                continue
            reaction_j = non_blocked_reactions[j]
            right = non_blocked_ns[j]
            ratio = np.apply_along_axis(lambda x: x[0] / x[1] if abs(x[1]) > 0. else np.inf, 0, np.array([left, right]))

            # special case:
            # if ratio is 1 (a/b == 1) then a == b.
            # but if a = 0 and b = 0, then a/b = np.inf.
            # solution:
            # mask inf from ratio
            # check if non-inf elements ratio is ~1
            # check if left and right values are 0 for indices with inf ratio
            # if yes, replace with 1
            inf_mask = np.isinf(ratio)
            non_inf = ratio[~inf_mask]

            if (abs((non_inf - 1)) < tol * 100).all():
                right_is_zero = (abs(right[inf_mask]) < tol).all()
                left_is_zero = (abs(left[inf_mask]) < tol).all()

                if right_is_zero and left_is_zero:
                    ratio[inf_mask] = 1

            if abs(max(ratio) - min(ratio)) < tol * 100:
                group[reaction_j] = round(ratio.mean(), 10)

    groups = [g for g in groups if len(g) > 1]

    return groups


def find_dead_end_reactions(model):
    """
    Identify reactions that are structurally prevented from carrying flux (dead ends).
    """
    stoichiometries = {}
    for reaction in model.reactions:
        for met, coef in reaction.metabolites.items():
            stoichiometries.setdefault(met.id, {})[reaction] = coef

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
            met_id: {reac: coef for reac, coef in stoichiometry.items() if reac not in new_blocked}
            for met_id, stoichiometry in stoichiometries.items()}
        blocked_reactions.update(new_blocked)

    return frozenset(blocked_reactions)


def find_coupled_reactions(model, return_dead_ends=False):
    """Find reaction sets that are structurally forced to carry equal flux"""
    blocked = find_dead_end_reactions(model)
    stoichiometries = {}
    for reaction in model.reactions:
        if reaction in blocked:
            continue
        for met, coef in reaction.metabolites.items():
            stoichiometries.setdefault(met.id, {})[reaction] = coef

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


class ShortestElementaryFluxModes():
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
                if isinstance(reaction, str):
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
        exchanges = self.model.boundary
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
                self.model.slim_optimize(error_value=None)
            except OptimizationError:
                raise StopIteration
            elementary_flux_mode = list()
            exclusion_list = list()
            for reaction in self._reactions:
                if reaction._indicator_variable_fwd.primal >= 0.9:
                    reaction_copy = copy(reaction)
                    reaction_copy.lower_bound = 0
                    elementary_flux_mode.append(reaction_copy)
                    exclusion_list.append(reaction._indicator_variable_fwd)
                elif reaction._indicator_variable_rev.primal >= 0.9:
                    reaction_copy = copy(reaction)
                    reaction_copy.upper_bound = 0
                    elementary_flux_mode.append(reaction_copy)
                    exclusion_list.append(reaction._indicator_variable_rev)
            exclusion_constraint = self.model.solver.interface.Constraint(sympy.Add(*exclusion_list),
                                                                          ub=len(exclusion_list) - 1)
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
                    problem_solution = self.model.solver.problem.solution
                    for reaction in self._reactions:
                        if problem_solution.pool.get_values(i, reaction._indicator_variable_fwd.name) >= 0.9:
                            reaction_copy = copy(reaction)
                            reaction_copy.lower_bound = 0
                            elementary_flux_mode.append(reaction_copy)
                            exclusion_list.append(reaction._indicator_variable_fwd)
                        elif problem_solution.pool.get_values(i, reaction._indicator_variable_rev.name) >= 0.9:
                            reaction_copy = copy(reaction)
                            reaction_copy.upper_bound = 0
                            elementary_flux_mode.append(reaction_copy)
                            exclusion_list.append(reaction._indicator_variable_rev)
                    exclusion_lists.append(exclusion_list)
                    yield elementary_flux_mode
                for exclusion_list in exclusion_lists:
                    exclusion_constraint = self.model.solver.interface.Constraint(sympy.Add(*exclusion_list),
                                                                                  ub=len(exclusion_list) - 1)
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
        return next(self._elementary_mode_generator)


class MinimalCutSetsEnumerator(ShortestElementaryFluxModes):  # pragma: no cover # don't test until it works
    def __init__(self, model, targets, constraints=None, exclude=None, c=1e-5):
        if exclude is not None:
            exclude_copy = []
            for re in exclude:
                if isinstance(re, str):
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
            with self._primal_model:
                for reac_id in mcs:
                    self._primal_model.reactions.get_by_id(reac_id).knock_out()
                self._primal_model.solver.optimize()
                if self._primal_model.solver.status == OPTIMAL:
                    return mcs
                else:
                    return None
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
        dual_model = Model(solver_interface=model.solver.interface)

        # Add dual metabolites
        dual_metabolite_names = []
        irreversibles = []  # Metabolites for which a z-reaction must be added
        self._dual_to_primal_mapping = {}
        for re in model.reactions:
            forward_id = re.id()
            reverse_id = re.reverse_id()
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
                if reaction.reverse_id() in dual_metabolite_names:
                    transposed_stoichiometry.setdefault(met, {})[
                        dual_model.metabolites.get_by_id(reaction.reverse_id())] = -coef
                if reaction.id() in dual_metabolite_names:
                    transposed_stoichiometry.setdefault(met, {})[
                        dual_model.metabolites.get_by_id(reaction.id())] = coef
                # This should be the same as forward_var.name but in general it might not be
                elif reaction.id in dual_metabolite_names:
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
        elif isinstance(targets, str):
            targets = [targets]
        if isinstance(targets[0], str):
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
                    if var.name in self._dual_model.metabolites}
            elif target.lb is not None:
                coefficients = {
                    self._dual_model.metabolites.get_by_id(var.name): -coef for var, coef in coefficients_dict.items()
                    if var.name in self._dual_model.metabolites}
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
                for constraint in constraints]
            self._primal_model.solver.add(cloned_constraints)
            illegal_knockouts = []
            for reaction in self._primal_model.reactions:
                # If single knockout causes the constrained model to become infeasible, then no superset
                # of knockouts can be feasible either.
                with self._primal_model:
                    reaction.knock_out()
                    self._primal_model.solver.optimize()
                    if self._primal_model.solver.status != OPTIMAL:
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
