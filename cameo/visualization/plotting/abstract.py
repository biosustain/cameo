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

import six
import collections

from cameo.visualization.palette import mapper, Palette

GOLDEN_RATIO = 1.618033988


class Grid(object):
    def __init__(self, engine, n_rows=1, width=None, height=None, title=None):
        self.engine = engine
        self.plots = []
        if n_rows <= 1:
            n_rows = 1
        self.n_rows = n_rows
        self.title = title

        if width is None or height is None:
            width, height = AbstractPlotter.golden_ratio(width, height)

        self.width = width
        self.height = height

    def add_plot(self, plot):
        self.plots.append(plot)

    def __enter__(self):
        return self

    def append(self, plot):
        self.plots.append(plot)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.plot()

    def plot(self):
        self.engine.display(self)


class AbstractPlotter(object):
    __default_options__ = {
        'color': 'blue',
        'palette': 'Pastel2',
        'width': 700,
        'alpha': 0.3
    }

    def __init__(self, **options):
        self.__options__ = {}
        self.__options__.update(options)

    def set_option(self, key, value):
        self.__options__[key] = value

    def get_option(self, key):
        return self.__options__.get(key, self.__default_options__.get(key, None))

    def production_envelope(self, dataframe, grid=None, width=None, height=None, title=None,
                            points=None, points_colors=None, palette=None, x_axis_label=None, y_axis_label=None):
        """
        Plots production envelopes from a pandas.DataFrame.

        The DataFrame format is:
            ub     lb     strain   value
            10     0      WT         0.4
            5      0      WT         0.5
            10     0      MT         0.4
            2      0      MT         0.5

        Arguments
        ---------
        dataframe: pandas.DataFrame
            The data to plot.

        grid: AbstractGrid
            An instance of AbstractGrid compatible with the plotter implementation. see AbstractPlotter.grid
        width: int
            Plot width
        height: int
            Plot height
        title:
            Plot title
        points: list
            list of (x, y) points to add to the plot
        points_colors: list
            list of colors for each point
        palette: list
            see color brewer
        x_axis_label: str
            X-axis label
        y_axis_label: str
            Y-axis label

        Returns
        -------
        a plottable object
            see AbstractPlotter.display
        """
        raise NotImplementedError

    def flux_variability_analysis(self, dataframe, grid=None, width=None, height=None, title=None,
                                  palette=None, x_axis_label=None, y_axis_label=None):
        """
        Plots flux variability analysis bars from a pandas.DataFrame.

        The DataFrame format is:
            ub     lb     variable   reaction
            10     0      WT         PFK1
            5      0      WT         ADT
            10     -10    MT         PFK1
            2      0      MT         ADT

        Arguments
        ---------
        dataframe: pandas.DataFrame
            The data to plot.
        grid: AbstractGrid
            An instance of AbstractGrid compatible with the plotter implementation. see AbstractPlotter.grid
        width: int
            Plot width
        height: int
            Plot height
        title:
            Plot title
        palette: list
            see color brewer
        x_axis_label: str
            X-axis label
        y_axis_label: str
            Y-axis label

        Returns
        -------
        a displayable object.
            If a grid is given, returns the grid with the plot in it's specific position.

        See Also
        --------
        AbstractPlotter.display
        """
        raise NotImplementedError

    @property
    def _display(self):
        raise NotImplementedError

    @staticmethod
    def _make_grid(grid):
        raise NotImplementedError

    def display(self, plot):
        """
        Displays an object using the implemented library

        """

        if isinstance(plot, Grid):
            plot = self._make_grid(plot)

        self._display(plot)

    @staticmethod
    def golden_ratio(width, height):
        if width is None and height is None:
            return width, height
        if width is None:
            width = int(height + height / GOLDEN_RATIO)

        elif height is None:
            height = int(width / GOLDEN_RATIO)

        return width, height

    @staticmethod
    def _palette(palette, number):
        if isinstance(palette, six.string_types):
            palette = mapper.map_palette(palette, number)

        if isinstance(palette, collections.Iterable):
            return palette
        elif isinstance(palette, Palette):
            return palette.hex_colors
        else:
            raise ValueError("Invalid palette %s" % palette)

    def grid(self, n_rows=1, width=None, height=None, title=None):
        return Grid(self, n_rows=n_rows, width=width, height=height, title=title)
