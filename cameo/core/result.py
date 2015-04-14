# -*- coding: utf-8 -*-
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

from __future__ import absolute_import, print_function

import time
from datetime import datetime
from optlang.interface import OptimizationExpression
import sympy
from sympy.parsing.sympy_parser import parse_expr
from cameo import system_info, config
import getpass

import pandas
from cameo.visualization import plotting


class MetaInformation(object):

    def __init__(self, *args, **kwargs):
        super(MetaInformation, self).__init__(*args, **kwargs)
        self._system_info = system_info
        self._responsible = getpass.getuser()
        self._timestamp = time.time()

    @property
    def system_info(self):
        return self._system_info

    @property
    def responsible(self):
        return self._responsible

    @property
    def timestamp(self):
        """doc string"""
        return self._timestamp

    @property
    def human_readable_timestamp(self):
        dt = datetime.fromtimestamp(self.timestamp)
        return dt.strftime('%Y-%m-%d %H:%M:%S.%f')


class Result(object):

    def __init__(self, *args, **kwargs):
        super(Result, self).__init__(*args, **kwargs)
        self._meta_information = MetaInformation()

    @property
    def meta_information(self):
        return self._meta_information

    @property
    def data_frame(self):
        raise NotImplementedError

    @property
    def plot(self):
        raise NotImplementedError


class FluxDistributionResult(Result):
    def __init__(self, solution, *args, **kwargs):
        super(FluxDistributionResult, self).__init__(*args, **kwargs)
        self._fluxes = solution.fluxes
        self._objective_value = solution.f

    def __getitem__(self, item):
        if isinstance(item, str):
            exp = parse_expr(item)
        elif isinstance(item, OptimizationExpression):
            exp = item.expression
        elif isinstance(item, sympy.Expr):
            exp = item
        else:
            raise KeyError(item)

        return exp.evalf(n=config.ndecimals,
                         subs={v: self.fluxes[v.name] for v in exp.atoms() if isinstance(v, sympy.Symbol)})

    @property
    def data_frame(self):
        return pandas.DataFrame(list(self._fluxes.values()), index=list(self._fluxes.keys()), columns=['flux'])

    @property
    def fluxes(self):
        return self._fluxes

    @property
    def objective_value(self):
        return self._objective_value

    @property
    def plot(self):
        pass


class PhenotypicPhasePlaneResult(Result):
    def __init__(self, phase_plane, reaction_ids, *args, **kwargs):
        super(PhenotypicPhasePlaneResult, self).__init__(*args, **kwargs)
        self._phase_plane = phase_plane
        self.reaction_ids = reaction_ids

    @property
    def data_frame(self):
        return pandas.DataFrame(self._phase_plane)

    @property
    def plot(self):
        for r_id in self.reaction_ids:
            plotting.plot_production_envelope(self._phase_plane, key=r_id)

    def __getitem__(self, item):
        return self._phase_plane[item]

    def iterrows(self):
        return self._phase_plane.iterrows()


if __name__ == '__main__':

    from cameo import load_model
    model = load_model('iJO1366')
    solution = model.solve()
    from cameo.util import Timer
    with Timer():
        result = FluxDistributionResult(solution)
        print(result.meta_information.system_info)
        print(result.meta_information.responsible)
        print(result.meta_information.timestamp)
        print(result.meta_information.human_readable_timestamp)
        print(result.data_frame)
