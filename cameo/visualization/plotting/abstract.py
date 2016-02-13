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


GOLDEN_RATIO = 1.618033988


class AbstractPlotter(object):

    __default_options__ = {
        'color': 'blue',
        'palette': 'Pastel2'
    }

    def __init__(self, **options):
        self.__options__ = {}
        self.__options__.update(options)

    def set_option(self, key, value):
        self.__options__[key] = value

    def get_option(self, key):
        if key in self.__options__:
            return self.__options__[key]
        elif key in self.__default_options__[key]:
            return self.__default_options__[key]
        else:
            return None

    def production_envelope(self, dataframe, values=None, grid=None, width=None, height=None, title=None,
                             points=None, points_colors=None, palette=None, xaxis_label=None, yaxis_label=None):
        """
        Plots production envelopes from a pandas.DataFrame.

        The DataFrame format is:
            ub     lb     variable   value
            10     0      WT         0.4
            5      0      WT         0.5
            10     0      MT         0.4
            2      0      MT         0.5

        Arguments
        ---------
        dataframe: pandas.DataFrame
            The data to plot.

        values: str
            The column of the DataFrame that identifies different envelopes
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

    def flux_variability_analysis(self, dataframe, values=None, grid=None, width=None, height=None, title=None,
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
        values: str
            The column of the DataFrame that identifies different FVA results
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

    @classmethod
    def display(cls, plot):
        """
        Displays an object using the implemented library

        """
        raise NotImplementedError

    @staticmethod
    def golden_ratio(width, height):
        if width is None:
            width = int(height + height/GOLDEN_RATIO)

        elif height is None:
            height = int(width/GOLDEN_RATIO)

        return width, height


class AbstractGrid(object):
    def __init__(self, n_rows=1, width=None, height=None, title=None):
        self.plots = []
        if n_rows <= 1:
            n_rows = 1
        self.n_rows = n_rows
        self.title = title
        self.width = width
        self.height = height

    def add_plot(self, plot):
        self.plots.append(plot)

    @property
    def plot(self):
        raise NotImplementedError
