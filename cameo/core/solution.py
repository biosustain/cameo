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

from collections import OrderedDict

import time
import datetime
from pandas import DataFrame, Series

import cobra

import cameo
from cameo.exceptions import UndefinedSolution

import logging

logger = logging.getLogger(__name__)


class SolutionBase(object):
    def __new__(cls, *args, **kwargs):
        # this is a cobrapy compatibility hack
        if len(args) == 1 and not isinstance(args[0], cameo.core.solver_based_model.SolverBasedModel):
            cobrapy_solution = super(SolutionBase, cls).__new__(cobra.core.Solution)
            cobrapy_solution.__init__(*args, **kwargs)
            return cobrapy_solution
        else:
            return super(SolutionBase, cls).__new__(cls)

    def __init__(self, model):
        self.model = model
        self._x = None
        self._y = None
        self._x_dict = None
        self._y_dict = None

    @property
    def data_frame(self):
        return DataFrame({'fluxes': Series(self.fluxes), 'reduced_costs': Series(self.reduced_costs)})

    def __str__(self):
        """A pandas DataFrame representation of the solution.

        Returns
        -------
        pandas.DataFrame
        """
        return str(self.data_frame)

    def _repr_html_(self):
        return self.data_frame._repr_html_()

    def as_cobrapy_solution(self):
        """Convert into a cobrapy Solution.

        Returns
        -------
        cobra.core.Solution.Solution
        """
        return Solution(self.f, x=self.x,
                        x_dict=self.x_dict, y=self.y, y_dict=self.y_dict,
                        the_solver=None, the_time=0, status=self.status)

    def get_primal_by_id(self, reaction_id):
        """Return a flux/primal value for a reaction.

        Parameters
        ----------
        reaction_id : str
            A reaction ID.
        """
        return self.x_dict[reaction_id]

    @property
    def x_dict(self):
        if self._x_dict is None:
            return self.fluxes
        else:
            return self._x_dict

    @x_dict.setter
    def x_dict(self, value):
        self._x_dict = value

    @property
    def x(self):
        if self._x is None:
            return self.fluxes.values()
        else:
            return self._x

    @x.setter
    def x(self, value):
        self._x = value

    @property
    def y_dict(self):
        if self._y_dict is None:
            return self.reduced_costs
        else:
            return self._y_dict

    @y_dict.setter
    def y_dict(self, value):
        self._y_dict = value

    @property
    def y(self):
        if self._y is None:
            return self.reduced_costs.values()
        else:
            return self._y

    @y.setter
    def y(self, value):
        self._y = value

    @property
    def objective_value(self):
        return self.f


class Solution(SolutionBase):
    """This class mimicks the cobrapy Solution class.

    Attributes
    ----------
    fluxes : OrderedDict
        A dictionary of flux values.
    reduced_costs : OrderedDict
        A dictionary of reduced costs.

    Notes
    -----
    See also documentation for cobra.core.Solution.Solution for an extensive list of inherited attributes.
    """

    def __init__(self, model, *args, **kwargs):
        """
        Parameters
        ----------
        model : SolverBasedModel
        """
        super(Solution, self).__init__(model, *args, **kwargs)
        self.f = model.solver.objective.value
        self.fluxes = OrderedDict()
        self.shadow_prices = OrderedDict()
        self.reduced_costs = OrderedDict()
        for reaction in model.reactions:
            self.fluxes[reaction.id] = reaction.flux
            self.reduced_costs[reaction.id] = reaction.reduced_cost
        for metabolite in model.metabolites:
            self.shadow_prices[metabolite.id] = self.model.solver.constraints[metabolite.id].dual
        self.status = model.solver.status
        self._reaction_ids = [r.id for r in self.model.reactions]
        self._metabolite_ids = [m.id for m in self.model.metabolites]

    def __dir__(self):
        # Hide 'cobrapy' attributes and methods from user.
        fields = sorted(dir(type(self)) + list(self.__dict__.keys()))
        fields.remove('x')
        fields.remove('y')
        fields.remove('x_dict')
        fields.remove('y_dict')
        return fields

    def _repr_html_(self):
        return "%s: %f" % (self.model.objective.expression, self.f)


class LazySolution(SolutionBase):
    """This class implements a lazy evaluating version of the cobrapy Solution class.

    Attributes
    ----------
    model : SolverBasedModel
    fluxes : OrderedDict
        A dictionary of flux values.
    reduced_costs : OrderedDict
        A dictionary of reduced costs.

    Notes
    -----
    See also documentation for cobra.core.Solution.Solution for an extensive list of inherited attributes.
    """

    def __init__(self, model, *args, **kwargs):
        """
        Parameters
        ----------
        model : SolverBasedModel
        """
        super(LazySolution, self).__init__(model, *args, **kwargs)
        if self.model._timestamp_last_optimization is not None:
            self._time_stamp = self.model._timestamp_last_optimization
        else:
            self._time_stamp = time.time()
        self._f = None

    def _check_freshness(self):
        """Raises an exceptions if the solution might have become invalid due to re-optimization of the attached model.

        Raises
        ------
        UndefinedSolution
            If solution has become invalid.
        """
        if self._time_stamp != self.model._timestamp_last_optimization:
            def timestamp_formatter(timestamp):
                datetime.datetime.fromtimestamp(timestamp).strftime(
                    "%Y-%m-%d %H:%M:%S:%f")

            raise UndefinedSolution(
                'The solution (captured around %s) has become invalid as the model has been re-optimized recently (%s).' % (
                    timestamp_formatter(self._time_stamp),
                    timestamp_formatter(self.model._timestamp_last_optimization))
            )

    @property
    def status(self):
        self._check_freshness()
        return self.model.solver.status

    @property
    def f(self):
        self._check_freshness()
        if self._f is None:
            return self.model.solver.objective.value
        else:
            return self._f

    @f.setter
    def f(self, value):
        self._f = value

    @property
    def fluxes(self):
        self._check_freshness()
        primals = OrderedDict()
        for reaction in self.model.reactions:
            primals[reaction.id] = reaction.flux
        return primals

    @property
    def reduced_costs(self):
        self._check_freshness()
        duals = OrderedDict()
        for reaction in self.model.reactions:
            duals[reaction.id] = reaction.reduced_cost
        return duals

    @property
    def shadow_prices(self):
        self._check_freshness()
        duals = OrderedDict()
        for metabolite in self.model.metabolites:
            duals[metabolite.id] = self.model.solver.constraints[metabolite.id].dual
        return duals

    def get_primal_by_id(self, reaction_id):
        """Return a flux/primal value for a reaction.

        Parameters
        ----------
        reaction_id : str
            A reaction ID.
        """
        self._check_freshness()
        return self.model.reactions.get_by_id(reaction_id).flux

    def __dir__(self):
        # Hide 'cobrapy' attributes and methods from user.
        fields = sorted(dir(type(self)) + list(self.__dict__.keys()))
        fields.remove('x')
        fields.remove('y')
        fields.remove('x_dict')
        fields.remove('y_dict')
        return fields
