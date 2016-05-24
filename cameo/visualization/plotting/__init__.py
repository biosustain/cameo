# Copyright 2016 Novo Nordisk Foundation Center for Biosustainability, DTU.
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
from __future__ import absolute_import

from warnings import warn

import six

from cameo.visualization.plotting.abstract import AbstractPlotter

__all__ = ["plotter"]

_engine = None
_engines = {}

try:
    from cameo.visualization.plotting.with_bokeh import BokehPlotter

    _engine = BokehPlotter() if _engine is None else _engine
    _engines["bokeh"] = BokehPlotter
except ImportError:
    pass

try:
    from cameo.visualization.plotting.with_plotly import PlotlyPlotter

    _engine = PlotlyPlotter(mode='offline') if _engine is None else _engine
    _engines["plotly"] = PlotlyPlotter
except ImportError:
    pass

try:
    from cameo.visualization.plotting.with_plotly import GGPlotPlotter

    _engine = GGPlotPlotter() if _engine is None else _engine
    _engines["ggplot"] = GGPlotPlotter
except ImportError:
    pass

if _engine is None:
    warn("Cannot import any plotting library. Please install one of 'plotly', 'bokeh' or 'ggplot' if you want to"
         " use any plotting function.")
    _engine = AbstractPlotter()


class _plotting:
    def __init__(self, engine):
        self.__dict__['_engine'] = engine

    def __getattr__(self, item):
        if item not in ["_engine", "engine"]:
            return getattr(self.__dict__['_engine'], item)
        else:
            return self.__dict__['_engine']

    def __setattr__(self, key, item):
        if key not in ["_engine", "engine"]:
            raise KeyError(key)
        else:
            if isinstance(item, six.string_types):
                item = _engines[item]()

            if not isinstance(item, AbstractPlotter):
                raise AssertionError("Invalid engine %s" % item)

            self.__dict__['_engine'] = item


plotter = _plotting(_engine)
