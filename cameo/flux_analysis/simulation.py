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

import logging
import os

import numpy
import pandas
import six
import sympy
from cobra import Reaction
from cobra.flux_analysis import pfba as cobrapy_pfba
from optlang.interface import OptimizationExpression
from sympy import Add
from sympy import Mul
from sympy.parsing.sympy_parser import parse_expr

from cameo.config import ndecimals
from cameo.core.result import Result
from cameo.util import ProblemCache, in_ipnb
from cameo.visualization.palette import mapper, Palette

__all__ = ['fba', 'pfba', 'moma', 'lmoma', 'room']

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
    model: cobra.Model
    objective: a valid objective - see SolverBaseModel.objective (optional)

    Returns
    -------
    FluxDistributionResult
        Contains the result of the linear solver.

    """
    with model:
        if objective is not None:
            model.objective = objective
        solution = model.optimize(raise_error=True)
        if reactions is not None:
            result = FluxDistributionResult({r: solution[r] for r in reactions}, solution.objective_value)
        else:
            result = FluxDistributionResult.from_solution(solution)
        return result


def pfba(model, objective=None, reactions=None, fraction_of_optimum=1, *args, **kwargs):
    """Parsimonious Enzyme Usage Flux Balance Analysis [1].

    Parameters
    ----------
    model : cobra.Model
        The model to perform pFBA with
    objective: str or reaction or optlang.Objective
        An objective to be minimized/maximized for
    reactions : list
        list of reactions to get results for. Getting fluxes from solution can be time consuming so if not all are
        needed it may be faster to request specific reactions.
    fraction_of_optimum : float
        Fix the value of the current objective to a fraction of is maximum.

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
    solution = cobrapy_pfba(model, objective=objective, fraction_of_optimum=fraction_of_optimum, reactions=reactions)
    return FluxDistributionResult.from_solution(solution)


def moma(model, reference=None, cache=None, reactions=None, *args, **kwargs):
    """
    Minimization of Metabolic Adjustment[1]

    Parameters
    ----------
    model: cobra.Model
    reference: FluxDistributionResult, dict
    cache: ProblemCache
    reactions: list

    Returns
    -------
    FluxDistributionResult
        Contains the result of the solver.

    References
    ----------
    .. [1] Segrè, D., Vitkup, D., & Church, G. M. (2002). Analysis of optimality in natural and perturbed metabolic
     networks. Proceedings of the National Academy of Sciences of the United States of America, 99(23), 15112–7.
     doi:10.1073/pnas.232349399

    """
    volatile = False
    if cache is None:
        volatile = True
        cache = ProblemCache(model)

    cache.begin_transaction()
    try:
        for rid, flux_value in six.iteritems(reference):

            def create_variable(model, variable_id):
                var = model.solver.interface.Variable(variable_id)
                return var

            var_id = "moma_aux_%s" % rid

            cache.add_variable(var_id, create_variable, None)

            def create_constraint(model, constraint_id, var, reaction, flux_value):
                constraint = model.solver.interface.Constraint(reaction.flux_expression - var,
                                                               lb=flux_value,
                                                               ub=flux_value,
                                                               name=constraint_id)
                return constraint

            def update_constraint(model, constraint, var, reaction, flux_value):
                if constraint.lb != flux_value:
                    constraint.lb = flux_value
                    constraint.ub = flux_value

            constraint_id = "moma_const_%s" % rid
            reaction = model.reactions.get_by_id(rid)
            cache.add_constraint(constraint_id, create_constraint, update_constraint,
                                 cache.variables[var_id], reaction, flux_value)

        def create_objective(model, variables):
            return model.solver.interface.Objective(Add(*[FloatOne * var ** 2 for var in variables]),
                                                    direction="min",
                                                    sloppy=True)

        cache.add_objective(create_objective, None, cache.variables.values())

        solution = model.optimize(raise_error=True)

        if reactions is not None:
            result = FluxDistributionResult(
                {r: solution.get_primal_by_id(r) for r in reactions}, solution.objective_value)
        else:
            result = FluxDistributionResult.from_solution(solution)
        return result
    except Exception as e:
        cache.rollback()
        raise e
    finally:
        if volatile:
            cache.reset()


# TODO: limiting step seems to be the update function. If the reference and cache can be linked
# TODO: there is no need for update
def lmoma(model, reference=None, cache=None, reactions=None, *args, **kwargs):
    """Linear Minimization Of Metabolic Adjustment [1].

    Parameters
    ----------
    model: cobra.Model
    reference: FluxDistributionResult, dict
    cache: ProblemCache
    reactions: list

    Returns
    -------
    FluxDistributionResult
        Contains the result of the solver.

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

    if not isinstance(reference, (dict, pandas.Series, FluxDistributionResult)):
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
                if constraint.lb != flux_value:
                    constraint.lb = flux_value

            def create_upper_constraint(model, constraint_id, var, reaction, flux_value):
                constraint = model.solver.interface.Constraint(reaction.flux_expression + var,
                                                               lb=flux_value,
                                                               sloppy=True,
                                                               name=constraint_id)
                return constraint

            cache.add_constraint("lmoma_const_%s_ub" % rid, create_upper_constraint, update_upper_constraint,
                                 cache.variables[pos_var_id], reaction, flux_value)

            def update_lower_constraint(model, constraint, var, reaction, flux_value):
                if constraint.ub != flux_value:
                    constraint.ub = flux_value

            def create_lower_constraint(model, constraint_id, var, reaction, flux_value):
                constraint = model.solver.interface.Constraint(reaction.flux_expression - var,
                                                               ub=flux_value,
                                                               sloppy=True,
                                                               name=constraint_id)
                return constraint

            cache.add_constraint("lmoma_const_%s_lb" % rid, create_lower_constraint, update_lower_constraint,
                                 cache.variables[neg_var_id], reaction, flux_value)

        def create_objective(model, variables):
            return model.solver.interface.Objective(add([mul((One, var)) for var in variables]),
                                                    direction="min",
                                                    sloppy=False)
        cache.add_objective(create_objective, None, cache.variables.values())

        solution = model.optimize(raise_error=True)
        if reactions is not None:
            result = FluxDistributionResult(
                {r: solution.get_primal_by_id(r) for r in reactions}, solution.objective_value)
        else:
            result = FluxDistributionResult.from_solution(solution)
        return result

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
    model: cobra.Model
    reference: FluxDistributionResult, dict
    delta: float
    epsilon: float
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

    if not isinstance(reference, (dict, pandas.Series, FluxDistributionResult)):
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
                constraint.set_linear_coefficients({variable: reaction.upper_bound - w_u})
                constraint.ub = w_u

            cache.add_constraint("room_const_%s_upper" % rid, create_upper_constraint, update_upper_constraint,
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
                constraint.set_linear_coefficients({variable: reaction.lower_bound - w_l})
                constraint.lb = w_l

            cache.add_constraint("room_const_%s_lower" % rid, create_lower_constraint, update_lower_constraint,
                                 reaction, cache.variables["y_%s" % rid], flux_value, epsilon)

        model.objective = model.solver.interface.Objective(add([mul([One, var]) for var in cache.variables.values()]),
                                                           direction='min')

        solution = model.optimize(raise_error=True)
        if reactions is not None:
            result = FluxDistributionResult(
                {r: solution.get_primal_by_id(r) for r in reactions}, solution.objective_value)
        else:
            result = FluxDistributionResult.from_solution(solution)
        return result

    except Exception as e:
        cache.rollback()
        raise e

    finally:
        if volatile:
            cache.reset()


class FluxDistributionResult(Result):
    """
    Contains a flux distribution of a simulation method.


    """

    @classmethod
    def from_solution(cls, solution, *args, **kwargs):
        return cls(solution.fluxes, solution.objective_value, *args, **kwargs)

    def __init__(self, fluxes, objective_value, *args, **kwargs):
        super(FluxDistributionResult, self).__init__(*args, **kwargs)
        self._fluxes = fluxes
        self._objective_value = objective_value

    def __getitem__(self, item):
        if isinstance(item, Reaction):
            return self.fluxes[item.id]
        elif isinstance(item, six.string_types):
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
        return pandas.DataFrame(list(self._fluxes.values), index=list(self._fluxes.keys()), columns=['flux'])

    @property
    def fluxes(self):
        return self._fluxes

    @property
    def objective_value(self):
        return self._objective_value

    def plot(self, grid=None, width=None, height=None, title=None):
        # TODO: Add barchart or something similar.
        raise NotImplementedError

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

    def plot_scale(self, values, palette="YlGnBu"):
        """
        Generates a color scale based on the flux distribution.
        It makes an array containing the absolute values and minus absolute values.

        The colors set as follows (p standsfor palette colors array):
        min   -2*std  -std      0      std      2*std   max
        |-------|-------|-------|-------|-------|-------|
        p[2]  p[2] .. p[1] .. p[0] .. p[1] ..  p[2]   p[2]


        Arguments
        ---------
        palette: Palette, list, str
            A Palette from palettable of equivalent, a list of colors (size 3) or a palette name

        Returns
        -------
        tuple
            ((-2*std, color), (-std, color) (0 color) (std, color) (2*std, color))
        """
        if isinstance(palette, six.string_types):
            palette = mapper.map_palette(palette, 3)
            palette = palette.hex_colors

        elif isinstance(palette, Palette):
            palette = palette.hex_colors

        std = numpy.std(values)

        return (-2 * std, palette[2]), (-std, palette[1]), (0, palette[0]), (std, palette[1]), (2 * std, palette[2])

    def display_on_map(self, map_name=None, palette="YlGnBu"):
        try:
            import escher
            if os.path.exists(map_name):
                map_json = map_name
                map_name = None
            else:
                map_json = None

            active_fluxes = {rid: flux for rid, flux in six.iteritems(self.fluxes) if abs(flux) > 10 ** -ndecimals}

            values = [abs(v) for v in active_fluxes.values()]
            values += [-v for v in values]

            scale = self.plot_scale(values, palette)

            reaction_scale = [dict(type='min', color=scale[0][1], size=24),
                              dict(type='value', value=scale[0][0], color=scale[0][1], size=21),
                              dict(type='value', value=scale[1][0], color=scale[1][1], size=16),
                              dict(type='value', value=scale[2][0], color=scale[2][1], size=8),
                              dict(type='value', value=scale[3][0], color=scale[3][1], size=16),
                              dict(type='value', value=scale[4][0], color=scale[4][1], size=21),
                              dict(type='max', color=scale[4][1], size=24)]

            active_fluxes = {rid: round(flux, ndecimals) for rid, flux in six.iteritems(self.fluxes)
                             if abs(flux) > 10 ** -ndecimals}

            active_fluxes['min'] = min(values)
            active_fluxes['max'] = max(values)

            builder = escher.Builder(map_name=map_name, map_json=map_json, reaction_data=active_fluxes,
                                     reaction_scale=reaction_scale)

            if in_ipnb():
                from IPython.display import display
                display(builder.display_in_notebook())
            else:
                builder.display_in_browser()

        except ImportError:
            print("Escher must be installed in order to visualize maps")


if __name__ == '__main__':
    import time
    from cobra.io import read_sbml_model
    from cameo import load_model

    # sbml_path = '../../tests/data/EcoliCore.xml'
    sbml_path = '../../tests/data/iJO1366.xml'

    cb_model = read_sbml_model(sbml_path)
    model = load_model(sbml_path)

    # model.solver = 'glpk'

    # print("cobra fba")
    # tic = time.time()
    # cb_model.optimize(solver='cglpk')
    # print("flux sum:", sum([abs(val) for val in list(cb_model.solution.fluxes.values())]))
    # print("cobra fba runtime:", time.time() - tic)

    # print("cobra pfba")
    # tic = time.time()
    # optimize_minimal_flux(cb_model, solver='cglpk')
    # print("flux sum:", sum([abs(val) for val in list(cb_model.solution.fluxes.values())]))
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
