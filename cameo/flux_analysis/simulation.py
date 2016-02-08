# -*- coding: utf-8 -*-
#  Copyright 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.
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
"""
Methods for simulating flux distributions.
Currently implements:
* fba - Flux Balance Analysis
* pfba - Parsimonious Flux Balance Analysis
* lmoma - (Linear) Minimization of Metabolic Adjustment
* room - Regulatory On/Off Minimization

"""

from __future__ import absolute_import, print_function
from optlang.interface import OptimizationExpression
import pandas
from sympy.parsing.sympy_parser import parse_expr
import cameo
from cameo.core.result import Result
from cameo.core.solution import SolutionBase

__all__ = ['fba', 'pfba', 'moma', 'lmoma', 'room']

import six

from functools import partial

import sympy
from sympy import Add
from sympy import Mul

from cameo.util import TimeMachine, ProblemCache
from cameo.exceptions import SolveError

import logging

logger = logging.getLogger(__name__)

add = Add._from_args
mul = Mul._from_args
NegativeOne = sympy.singleton.S.NegativeOne
One = sympy.singleton.S.One
FloatOne = sympy.Float(1)
RealNumber = sympy.RealNumber


def fba(model, objective=None, reactions=None, *args, **kwargs):
    """Flux Balance Analysis.

    Parameters
    ----------
    model: SolverBasedModel
    objective: a valid objective - see SolverBaseModel.objective (optional)

    Returns
    -------
    FluxDistributionResult
        Contains the result of the linear solver.

    """
    with TimeMachine() as tm:
        if objective is not None:
            tm(do=partial(setattr, model, 'objective', objective),
               undo=partial(setattr, model, 'objective', model.objective))
        solution = model.solve()
        if reactions is not None:
            result = FluxDistributionResult({r: solution.get_primal_by_id(r) for r in reactions}, solution.f)
        else:
            result = FluxDistributionResult.from_solution(solution)
        return result


def pfba(model, objective=None, reactions=None, *args, **kwargs):
    """Parsimonious Enzyme Usage Flux Balance Analysis [1].

    Parameters
    ----------
    model: SolverBasedModel
    objective: str or reaction or optlang.Objective
        An objective to be minimized/maximized for

    Returns
    -------
    FluxDistributionResult
        Contains the result of the linear solver.

    References
    ----------
    .. [1] Lewis, N. E., Hixson, K. K., Conrad, T. M., Lerman, J. A., Charusanti, P., Polpitiya, A. D., …
     Palsson, B. Ø. (2010). Omic data from evolved E. coli are consistent with computed optimal growth from
     genome-scale models. Molecular Systems Biology, 6, 390. doi:10.1038/msb.2010.47

    """
    with TimeMachine() as tm:
        original_objective = model.objective
        if objective is not None:
            tm(do=partial(setattr, model, 'objective', objective),
               undo=partial(setattr, model, 'objective', original_objective))
        try:
            obj_val = model.solve().f
        except SolveError as e:
            logger.debug("pfba could not determine maximum objective value for\n%s." % model.objective)
            raise e
        if model.objective.direction == 'max':
            fix_obj_constraint = model.solver.interface.Constraint(model.objective.expression, lb=obj_val)
        else:
            fix_obj_constraint = model.solver.interface.Constraint(model.objective.expression, ub=obj_val)
        tm(do=partial(model.solver._add_constraint, fix_obj_constraint),
           undo=partial(model.solver._remove_constraint, fix_obj_constraint))
        pfba_obj = model.solver.interface.Objective(add(
            [mul((sympy.singleton.S.One, variable)) for variable in list(model.solver.variables.values())]),
            direction='min', sloppy=True)
        tm(do=partial(setattr, model, 'objective', pfba_obj),
           undo=partial(setattr, model, 'objective', original_objective))
        try:
            solution = model.solve()
            if reactions is not None:
                result = FluxDistributionResult({r: solution.get_primal_by_id(r) for r in reactions}, solution.f)
            else:
                result = FluxDistributionResult.from_solution(solution)
            tm.reset()
            return result
        except SolveError as e:
            logger.error("pfba could not determine an optimal solution for objective %s" % model.objective)
            raise e


def moma(model, reference=None, *args, **kwargs):
    """
    Minimization of Metabolic Adjustment
    """
    with TimeMachine() as tm:
        aux_vars = {}
        for reac_id in reference.keys():
            var = model.solver.interface.Variable("moma_aux_"+reac_id)
            aux_vars[reac_id] = var
        tm(do=partial(model.solver.add, aux_vars.values()),
           undo=partial(model.solver.remove, aux_vars.values()))

        constraints = {}
        for reac_id, aux_var in aux_vars.items():
            reac = model.reactions.get_by_id(reac_id)
            const = model.solver.interface.Constraint(
                add([mul([One, reac.forward_variable]),
                     mul([NegativeOne, reac.reverse_variable]),
                     mul([NegativeOne, aux_var])
                ]),
                lb=reference[reac_id],
                ub=reference[reac_id],
                sloppy=True
            )
            constraints[reac_id] = const
        tm(do=partial(model.solver.add, constraints.values()),
           undo=partial(model.solver.remove, constraints.values()))

        obj = model.solver.interface.Objective(Add(*(FloatOne*var**2 for var in aux_vars.values())), direction="min", sloppy=True)
        tm(do=partial(setattr, model, "objective", obj),
           undo=partial(setattr, model, "objective", model.objective))

        sol = model.solve()
        result = FluxDistributionResult.from_solution(sol)

    return result


def lmoma(model, reference=None, cache=None, reactions=None, *args, **kwargs):
    """Linear Minimization Of Metabolic Adjustment [1].

    Parameters
    ----------
    model: SolverBasedModel
    reference: dict
    objective: str or reaction or optlang.Objective
        An objective to be minimized/maximized for
    volatile: boolean
    cache: dict

    Returns
    -------
    FluxDistributionResult
        Contains the result of the linear solver.

    References
    ----------
    .. [1] Becker, S. A., Feist, A. M., Mo, M. L., Hannum, G., Palsson, B. Ø., & Herrgard, M. J. (2007).
    Quantitative prediction of cellular metabolism with constraint-based models: the COBRA Toolbox.
    Nature Protocols, 2(3), 727–38. doi:10.1038/nprot.2007.99

    """
    volatile = False
    if cache is None:
        volatile = True
        cache = ProblemCache(model)

    cache.begin_transaction()

    if not isinstance(reference, (dict, FluxDistributionResult)):
        raise TypeError("reference must be a flux distribution (dict or FluxDistributionResult")

    try:
        for rid, flux_value in six.iteritems(reference):
            reaction = model.reactions.get_by_id(rid)

            def create_variable(model, var_id, lb):
                var = model.solver.interface.Variable(var_id, lb=lb)
                return var

            pos_var_id = "u_%s_pos" % rid
            cache.add_variable(pos_var_id, create_variable, None, 0)

            neg_var_id = "u_%s_neg" % rid
            cache.add_variable(neg_var_id, create_variable, None, 0)

            # ui = vi - wt
            def update_upper_constraint(model, constraint, var, reaction, flux_value):
                constraint.lb = flux_value

            def create_upper_constraint(model, constraint_id, var, reaction, flux_value):
                constraint = model.solver.interface.Constraint(reaction.flux_expression + var,
                                                               lb=flux_value,
                                                               sloppy=True,
                                                               name=constraint_id)
                return constraint

            cache.add_constraint("c_%s_ub" % rid, create_upper_constraint, update_upper_constraint,
                                 cache.variables[pos_var_id], reaction, flux_value)

            def update_lower_constraint(model, constraint, var, reaction, flux_value):
                constraint.ub = flux_value

            def create_lower_constraint(model, constraint_id, var, reaction, flux_value):
                constraint = model.solver.interface.Constraint(reaction.flux_expression - var,
                                                               ub=flux_value,
                                                               sloppy=True,
                                                               name=constraint_id)
                return constraint

            cache.add_constraint("c_%s_lb" % rid, create_lower_constraint, update_lower_constraint,
                                 cache.variables[neg_var_id], reaction, flux_value)

        model.objective = model.solver.interface.Objective(add([mul([One, var]) for var in cache.variables.values()]),
                                                           direction="min")

        try:

            solution = model.solve()
            if reactions is not None:
                result = FluxDistributionResult({r: solution.get_primal_by_id(r) for r in reactions}, solution.f)
            else:
                result = FluxDistributionResult.from_solution(solution)
            return result
        except SolveError as e:
            raise e
    except Exception as e:
        cache.rollback()
        raise e

    finally:
        if volatile:
            cache.reset()


def room(model, reference=None, cache=None, delta=0.03, epsilon=0.001, reactions=None, *args, **kwargs):
    """Regulatory On/Off Minimization [1].

    Parameters
    ----------
    model: SolverBasedModel
    reference: dict
    cache: ProblemCache

    Returns
    -------
    FluxDistributionResult
        Contains the result of the linear solver.

    References
    ----------
    .. [1] Tomer Shlomi, Omer Berkman and Eytan Ruppin, "Regulatory on/off minimization of metabolic
    flux changes after genetic perturbations", PNAS 2005 102 (21) 7695-7700; doi:10.1073/pnas.0406346102
    """
    volatile = False
    if cache is None:
        volatile = True
        cache = ProblemCache(model)
    elif not isinstance(cache, ProblemCache):
        raise TypeError("Invalid cache object (must be a cameo.util.ProblemCache)")

    cache.begin_transaction()

    if not isinstance(reference, (dict, FluxDistributionResult)):
        raise TypeError("reference must be a flux distribution (dict or FluxDistributionResult")

    try:
        for rid, flux_value in six.iteritems(reference):
            reaction = model.reactions.get_by_id(rid)

            def create_variable(model, var_id):
                return model.solver.interface.Variable(var_id, type="binary")

            cache.add_variable("y_%s" % rid, create_variable, None)

            def create_upper_constraint(model, constraint_id, reaction, variable, flux_value, epsilon):
                w_u = flux_value + delta * abs(flux_value) + epsilon
                return model.solver.interface.Constraint(
                    reaction.flux_expression - variable * (reaction.upper_bound - w_u),
                    ub=w_u,
                    sloppy=True,
                    name=constraint_id)

            def update_upper_constraint(model, constraint, reaction, variable, flux_value, epsilon):
                w_u = flux_value + delta * abs(flux_value) + epsilon
                constraint._set_coefficients_low_level({variable: reaction.upper_bound - w_u})
                constraint.ub = w_u

            cache.add_constraint("c_%s_upper" % rid, create_upper_constraint, update_upper_constraint,
                                 reaction, cache.variables["y_%s" % rid], flux_value, epsilon)

            def create_lower_constraint(model, constraint_id, reaction, variable, flux_value, epsilon):
                w_l = flux_value - delta * abs(flux_value) - epsilon
                return model.solver.interface.Constraint(
                    reaction.flux_expression - variable * (reaction.lower_bound - w_l),
                    lb=w_l,
                    sloppy=True,
                    name=constraint_id)

            def update_lower_constraint(model, constraint, reaction, variable, flux_value, epsilon):
                w_l = flux_value - delta * abs(flux_value) - epsilon
                constraint._set_coefficients_low_level({variable: reaction.lower_bound - w_l})
                constraint.lb = w_l

            cache.add_constraint("c_%s_lower" % rid, create_lower_constraint, update_lower_constraint,
                                 reaction, cache.variables["y_%s" % rid], flux_value, epsilon)

        model.objective = model.solver.interface.Objective(add([mul([One, var]) for var in cache.variables.values()]),
                                                           direction='min')
        try:
            solution = model.solve()
            if reactions is not None:
                result = FluxDistributionResult({r: solution.get_primal_by_id(r) for r in reactions}, solution.f)
            else:
                result = FluxDistributionResult.from_solution(solution)
            return result
        except SolveError as e:
            logger.error("room could not determine an optimal solution for objective %s" % model.objective)
            raise e

    except Exception as e:
        cache.rollback()
        raise e

    finally:
        if volatile:
            cache.reset()


class FluxDistributionResult(Result):
    @classmethod
    def from_solution(cls, solution, *args, **kwargs):
        return  cls(solution.fluxes, solution.f, *args, **kwargs)

    def __init__(self, fluxes, objective_value, *args, **kwargs):
        super(FluxDistributionResult, self).__init__(*args, **kwargs)
        self._fluxes = fluxes
        self._objective_value = objective_value

    def __getitem__(self, item):
        if isinstance(item, cameo.Reaction):
            return self.fluxes[item.id]
        elif isinstance(item, str):
            try:
                return self.fluxes[item]
            except KeyError:
                exp = parse_expr(item)
        elif isinstance(item, OptimizationExpression):
            exp = item.expression
        elif isinstance(item, sympy.Expr):
            exp = item
        else:
            raise KeyError(item)

        return exp.evalf(subs={v: self.fluxes[v.name] for v in exp.atoms(sympy.Symbol)})

    @property
    def data_frame(self):
        return pandas.DataFrame(list(self._fluxes.values()), index=list(self._fluxes.keys()), columns=['flux'])

    @property
    def fluxes(self):
        return self._fluxes

    @property
    def objective_value(self):
        return self._objective_value

    def plot(self, grid=None, width=None, height=None, title=None):
        # TODO: Add barchart or something similar.
        pass

    def iteritems(self):
        return six.iteritems(self.fluxes)

    def items(self):
        return six.iteritems(self.fluxes)

    def keys(self):
        return self.fluxes.keys()

    def values(self):
        return self.fluxes.values()

    def _repr_html_(self):
        return "<strong>objective value: %s</strong>" % self.objective_value


if __name__ == '__main__':
    import time
    from cobra.io import read_sbml_model
    from cobra.flux_analysis.parsimonious import optimize_minimal_flux
    from cameo import load_model

    # sbml_path = '../../tests/data/EcoliCore.xml'
    sbml_path = '../../tests/data/iJO1366.xml'

    cb_model = read_sbml_model(sbml_path)
    model = load_model(sbml_path)

    # model.solver = 'glpk'

    # print("cobra fba")
    # tic = time.time()
    # cb_model.optimize(solver='cglpk')
    # print("flux sum:", sum([abs(val) for val in list(cb_model.solution.x_dict.values())]))
    # print("cobra fba runtime:", time.time() - tic)

    # print("cobra pfba")
    # tic = time.time()
    # optimize_minimal_flux(cb_model, solver='cglpk')
    # print("flux sum:", sum([abs(val) for val in list(cb_model.solution.x_dict.values())]))
    # print("cobra pfba runtime:", time.time() - tic)

    print("pfba")
    tic = time.time()
    ref = pfba(model)
    print("flux sum:", end=' ')
    print(sum([abs(val) for val in list(ref.values())]))
    print("cameo pfba runtime:", time.time() - tic)

    print("lmoma")
    cache = ProblemCache(model)
    tic = time.time()
    res = lmoma(model, reference=ref, cache=cache)
    print("flux distance:", end=' ')
    print(sum([abs(res[v] - ref[v]) for v in list(res.keys())]))
    print("cameo lmoma runtime:", time.time() - tic)

    print("room")
    tic = time.time()
    res = room(model, reference=ref)
    print(sum([abs(res[v] - ref[v]) for v in list(res.keys())]))
    print("cameo room runtime:", time.time() - tic)

    print("flux distance:", end=' ')
    print(sum([abs(res[v] - ref[v]) for v in list(res.keys())]))
    print("sum yi:", res.objective_value)
    print("cameo room runtime:", time.time() - tic)

    print("lmoma w/ ko")
    tic = time.time()
    model.reactions.PGI.lower_bound = 0
    model.reactions.PGI.upper_bound = 0
    res = lmoma(model, reference=ref, cache=cache)
    print("flux distance:", end=' ')
    print(sum([abs(res[v] - ref[v]) for v in list(res.keys())]))
    print("cameo lmoma runtime:", time.time() - tic)

    # print model.solver
