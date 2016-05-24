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
from bokeh.models import GridPlot, FactorRange
from bokeh.plotting import figure, show

from cameo.util import partition, inheritdocstring
from cameo.visualization.plotting.abstract import AbstractPlotter


@six.add_metaclass(inheritdocstring)
class BokehPlotter(AbstractPlotter):
    def __init__(self, **options):
        super(BokehPlotter, self).__init__(**options)

    def _add_fva_bars(self, plot, factors, dataframe, height, step, color, strain):
        alpha = self.get_option('alpha')

        left = dataframe.ub.tolist()
        right = dataframe.lb.tolist()
        bottom = [(i + 1.5) + height for i in range(len(factors))]
        top = [(i + 1.5) + height + step for i in range(len(factors))]

        plot.quad(top=top, bottom=bottom, left=left, right=right, color=color, legend=strain, fill_alpha=alpha)

    def flux_variability_analysis(self, dataframe, grid=None, width=None, height=None, title=None,
                                  palette=None, x_axis_label=None, y_axis_label=None):
        palette = self.get_option('palette') if palette is None else palette
        width = self.get_option('width') if width is None else width
        strains = dataframe["strain"].unique()
        factors = dataframe['reaction'].unique().tolist()
        x_range = [min(dataframe.lb) - 5, max(dataframe.ub) + 5]

        width, height = self.golden_ratio(width, height)

        plot = figure(title=title, plot_width=width, plot_height=height,
                      x_range=x_range, y_range=FactorRange(factors=[""] + factors + [""]),
                      min_border=5)

        plot.line([0, 0], [0, len(factors) + 3], line_color="black", line_width=1.5, line_alpha=1)

        n = len(strains)
        step = 1.0 / float(len(strains))
        for strain, i, color in zip(strains, range(n), self._palette(palette, n)):
            _dataframe = dataframe[dataframe["strain"] == strain]
            self._add_fva_bars(plot, factors, _dataframe, i * step, step, color, strain)

        if x_axis_label:
            plot.xaxis.axis_label = x_axis_label
        if y_axis_label:
            plot.yaxis.axis_label = y_axis_label

        return plot

    def _add_production_envelope(self, plot, dataframe, strain, color=None):
        patch_alpha = self.get_option('alpha')
        ub = dataframe["ub"].values.tolist()
        lb = dataframe["lb"].values.tolist()
        var = dataframe["value"].values

        x = [v for v in var] + [v for v in reversed(var)]
        y = [v for v in lb] + [v for v in reversed(ub)]

        plot.patch(x=x, y=y, fill_color=color, fill_alpha=patch_alpha, legend=strain)

        if "label" in dataframe.columns:
            plot.text(dataframe["label"].values, var, ub)

        plot.line(var, ub, color=color)
        plot.line(var, lb, color=color)
        if ub[-1] != lb[-1]:
            plot.line((var[-1], var[-1]), (ub[-1], lb[-1]), color=color)

    def production_envelope(self, dataframe, grid=None, width=None, height=None, title=None, points=None,
                            points_colors=None, palette='RdYlBu', x_axis_label=None, y_axis_label=None):
        strains = dataframe["strain"].unique()
        palette = self.get_option('palette') if palette is None else palette
        width = self.get_option('width') if width is None else width
        if not height:
            width, height = self.golden_ratio(width, height)

        plot = figure(title=title, plot_width=width, plot_height=height)
        for strain, color in zip(strains, self._palette(palette, len(strains))):
            _dataframe = dataframe[dataframe["strain"] == strain]
            self._add_production_envelope(plot, _dataframe, strain, color=color)

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

    @property
    def _display(self):
        return show

    @staticmethod
    def _make_grid(grid):
        return GridPlot(children=partition(grid.plots, grid.n_rows),
                        name=grid.title)
