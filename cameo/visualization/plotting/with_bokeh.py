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

from bokeh.models import GridPlot, FactorRange
from bokeh.plotting import figure, show
from bokeh.palettes import brewer

from cameo.util import partition
from cameo.visualization.plotting.abstract import AbstractPlotter, AbstractGrid


class BokehPlotter(AbstractPlotter):

    def __init__(self, **options):
        if 'alpha' not in options:
            options['alpha'] = 0.3
        if 'palette' not in options:
            options['palette'] = 'RdYlBu'
        super(BokehPlotter, self).__init__(**options)

    def _add_production_envelope(self, plot, dataframe, values, variable, patch_color=None, line_color=None):
        patch_alpha = self.get_option('alpha')
        ub = dataframe["ub"].values.tolist()
        lb = dataframe["lb"].values.tolist()
        var = dataframe[values].values

        x = [v for v in var] + [v for v in reversed(var)]
        y = [v for v in lb] + [v for v in reversed(ub)]

        plot.patch(x=x, y=y, color=patch_color, alpha=patch_alpha, legend=variable)

        if "label" in dataframe.columns:
            plot.text(dataframe["label"].values, var, ub)

        plot.line(var, ub, color=line_color)
        plot.line(var, lb, color=line_color)
        if ub[-1] != lb[-1]:
            plot.line((var[-1], var[-1]), (ub[-1], lb[-1]), color=line_color)

    def _add_fva_bars(self, plot, factors, dataframe, h_bottom, h_top, color, variable):
        alpha = self.get_option('alpha')
        left = list(dataframe.ub)
        right = list(dataframe.lb)
        bottom = [i+h_bottom for i in range(0, len(factors))]
        top = [i+h_top for i in range(0, len(factors))]

        plot.quad(top=top, bottom=bottom, left=left, right=right, color=color, alpha=alpha, legend=variable)

    def flux_variability_analysis(self, dataframe, values=None, grid=None, width=None, height=None, title=None,
                                  palette=None, x_axis_label=None, y_axis_label=None):
        palette = self.get_option('palette') if palette is None else palette
        width = self.get_option('width') if width is None else width
        variables = dataframe[values].unique()
        factors = dataframe['reaction'].unique().tolist()
        x_range = [min(dataframe.lb), max(dataframe.ub)]

        plot = figure(title=title, plot_width=width, plot_height=height,
                      x_range=x_range, y_range=FactorRange(factors=factors))

        n = len(variables)
        h = 1.0/float(len(variables))
        for variable, i, color in zip(variables, range(1, n+1), self._palette(palette, n)):
            _dataframe = dataframe[dataframe[values] == variable]
            self._add_fva_bars(plot, factors, _dataframe, i*h, (i+1)*h, color, variable)

        if x_axis_label:
            plot.xaxis.axis_label = x_axis_label
        if y_axis_label:
            plot.yaxis.axis_label = y_axis_label

        return plot

    def production_envelope(self, dataframe, values=None, grid=None, width=None, height=None, title=None, points=None,
                            points_colors=None, palette='RdYlBu', x_axis_label=None, y_axis_label=None):
        variables = dataframe[values].unique()
        palette = self.get_option('palette') if palette is None else palette
        width = self.get_option('width') if width is None else width
        if not height:
            width, height = self.golden_ratio(width, height)

        plot = figure(title=title, plot_width=width, plot_height=height)
        for variable, color in zip(variables, self._palette(palette, len(variables))):
            _dataframe = dataframe[dataframe[values] == variable]
            self._add_production_envelope(plot, _dataframe, values, variable, patch_color=color, line_color=color)

        if points is not None:
            plot.scatter(*zip(*points), color="green" if points_colors is None else points_colors)

        if x_axis_label:
            plot.xaxis.axis_label = x_axis_label
        if y_axis_label:
            plot.yaxis.axis_label = y_axis_label

        if grid is not None:
            grid.append(plot)
            return grid

        return plot

    @staticmethod
    def _palette(palette, number, **kwargs):
        if isinstance(palette, six.string_types):
            n = 3 if number < 3 else number
            return brewer[palette][n]
        elif isinstance(palette, collections.Iterable):
            return palette
        else:
            raise ValueError("Invalid palette %s" % palette)

    @classmethod
    def display(cls, plot):
        if isinstance(plot, BokehGrid):
            plot = plot.plot
        show(plot)


class BokehGrid(AbstractGrid):
    @property
    def plot(self):
        width, height = AbstractPlotter.golden_ratio(self.width, self.height)

        return GridPlot(children=partition(self.plots, self.n_rows),
                        title=self.title,
                        plot_width=width,
                        plot_heigt=height)
