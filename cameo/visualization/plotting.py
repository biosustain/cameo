# Copyright 2015 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


__all__ = ['plot_production_envelope']
from cameo.util import partition
from cameo import config, util

GOLDEN_RATIO = 1.618033988


try:
    # import matplotlib.pyplot as plt

    def plot_production_envelope_ipython_matplotlib(envelope, objective, key, grid=None, width=None, height=None,
                                                    title=None, points=None, points_colors=None, axis_font_size=None):
        pass
        # plt.plot(envelope["objective_upper_bound"], envelope[key], title="Production envelop")
        # plt.xlabel("growth")
        # plt.ylabel(key)

    def plot_flux_variability_analysis_ipython_matplotlib(fva_result, grid=None, width=None, height=None, title=None,
                                                          axis_font_size=None):
        pass

except ImportError:

    def plot_production_envelope_ipython_matplotlib(envelope, objective, key, grid=None, width=None, height=None,
                                                    title=None, points=None, points_colors=None, axis_font_size=None):
        pass

    def plot_flux_variability_analysis_ipython_matplotlib(fva_result, grid=None, width=None, height=None, title=None,
                                                          axis_font_size=None):
        pass

try:
    from bokeh import plotting
    from bokeh.models import GridPlot

    def plot_production_envelope_ipython_bokeh(envelope, objective, key, grid=None, width=None, height=None,
                                               title=None, points=None, points_colors=None, axis_font_size=None):

        p = plotting.figure(title=title if title is not None else "Production envelope",
                            tools="save",
                            plot_width=width if width is not None else 700,
                            plot_height=height if height is not None else 700)
        p.xaxis.axis_label = key
        p.yaxis.axis_label = objective

        ub = envelope["objective_upper_bound"].values
        lb = envelope["objective_lower_bound"].values
        var = envelope[key].values

        x = [0] + [v for v in var] + list(reversed([v for v in var]))
        y = [0] + [v for v in lb] + list(reversed([v for v in ub]))

        p.patch(x=x, y=y, color="#99d8c9", alpha=0.3)

        if "label" in envelope.columns:
            p.text(envelope["label"].values, var, ub)
        p.line(var, ub, color="blue")
        p.line(var, lb, color="blue")
        if ub[-1] != lb[-1]:
            p.line((var[-1], var[-1]), (ub[-1], lb[-1]), color="blue")

        if axis_font_size is not None:
            p.xaxis.axis_label_text_font_size = axis_font_size
            p.yaxis.axis_label_text_font_size = axis_font_size

        if points is not None:
            p.scatter(*zip(*points), color="green" if points_colors is None else points_colors)

        if grid is not None:
            grid.append(p)
        else:
            plotting.show(p)

    def plot_flux_variability_analysis_ipython_bokeh(fva_result, grid=None, width=None, height=None, title=None,
                                                     axis_font_size=None):

        factors = list(fva_result.index)
        x0 = fva_result['lower_bound'].values
        x1 = fva_result['upper_bound'].values

        x_range = [min([min(x0), min(x1)]) -5, max([max(x0), max(x1)]) + 5]

        p = plotting.figure(title=title if title is not None else "Production envelope",
                            tools="save",
                            plot_width=width if width is not None else 700,
                            plot_height=height if height is not None else 700,
                            y_range=factors, x_range=x_range)

        p.segment(x0, factors, x1, factors, line_width=10, line_color="#99d8c9")
        p.line([0, 0], [-1, len(factors)+1], line_width=2, color="black")

        if grid is not None:
            grid.append(p)
        else:
            plotting.show(p)

except ImportError:

    def plot_production_envelope_ipython_bokeh(envelope, objective, key, grid=None, width=None, height=None,
                                               title=None, points=None, points_colors=None, axis_font_size=None):
        pass

    def plot_flux_variability_analysis_ipython_bokeh(fva_result, grid=None, width=None, height=None, title=None,
                                                     axis_font_size=None):
        pass

try:
    from bashplotlib import scatterplot

    def plot_production_envelope_cli(envelope, objective, key, grid=None, width=None, height=None,
                                     title=None, points=None, points_colors=None, axis_font_size=None):
        scatterplot.plot_scatter(None, envelope[key], envelope["objective_upper_bound"], "*")

    def plot_flux_variability_analysis_cli(fva_result, grid=None, width=None, height=None, title=None,
                                           axis_font_size=None):
        pass


except ImportError:
    def plot_production_envelope_cli(envelope, objective, key, grid=None, width=None, height=None,
                                     title=None, points=None, points_colors=None, axis_font_size=None):
        pass

    def plot_flux_variability_analysis_cli(fva_result, grid=None, width=None, height=None, title=None,
                                           axis_font_size=None):
        pass


def plot_production_envelope(envelope, objective, key, grid=None, width=None, height=None, title=None,
                             points=None, points_colors=None, axis_font_size=None):

    if width is None and height is None:
        width = 700
    if width is None or height is None:
        width, height = _golden_ratio(width, height)
    if util.in_ipnb():
        if config.use_bokeh:
            plot_production_envelope_ipython_bokeh(envelope, objective, key, grid=grid, width=width, height=height,
                                                   title=title, points=points, points_colors=points_colors,
                                                   axis_font_size=axis_font_size)
        elif config.use_matplotlib:
            plot_production_envelope_ipython_matplotlib(envelope, objective, key, grid=grid, width=width, height=height,
                                                        title=title, points=points, points_colors=points_colors,
                                                        axis_font_size=axis_font_size)
    else:
        plot_production_envelope_cli(envelope, objective, key, width=width, height=height, title=title, points=points,
                                     points_colors=points_colors, axis_font_size=axis_font_size)


def plot_flux_variability_analysis(fva_result, grid=None, width=None, height=None, title=None, axis_font_size=None):
    if width is None and height is None:
        width = 700
    if width is None or height is None:
        width, height = _golden_ratio(width, height)
    if util.in_ipnb():
        if config.use_bokeh:
            plot_flux_variability_analysis_ipython_bokeh(fva_result, grid=grid, width=width, height=height,
                                                         title=title, axis_font_size=axis_font_size)
        elif config.use_matplotlib:
            plot_flux_variability_analysis_matplotlib(fva_result, grid=grid, width=width, height=height,
                                                      title=title, axis_font_size=axis_font_size)
    else:
        plot_flux_variability_analysis_cli(fva_result, grid=grid, width=width, height=height,
                                           title=title, axis_font_size=axis_font_size)


class Grid(object):
    def __init__(self, nrows=1, title=None):
        self.plots = []
        if nrows <= 1:
            nrows = 1
        self.nrows = nrows
        self.title = title

    def __enter__(self):
        return self

    def append(self, plot):
        self.plots.append(plot)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._plot_grid()

    def _plot_grid(self):
        if util.in_ipnb():
            if config.use_bokeh:
                self._plot_bokeh_grid()
            elif config.use_matplotlib:
                self._plot_matplotlib_grid()
        else:
            self._plot_cli_grid()

    def _plot_bokeh_grid(self):
        if len(self.plots) > 0:
            grid = GridPlot(children=partition(self.plots, self.nrows), title=self.title)
            plotting.show(grid)


def _golden_ratio(width, height):
    if width is None:
        width = int(height + height/GOLDEN_RATIO)

    elif height is None:
        height = int(width/GOLDEN_RATIO)

    return width, height