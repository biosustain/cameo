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

from cameo.util import partition
from cameo import config, util

import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)

__all__ = ['plot_production_envelope', 'plot_2_production_envelopes', 'plot_flux_variability_analysis']

MISSING_PLOTTING_BACKEND_MESSAGE = "No supported plotting backend could be found. Please install bokeh if you'd like to generate plots (other backends will be supported in the future)."
GOLDEN_RATIO = 1.618033988

try:
    from bokeh import plotting
    from bokeh.models import GridPlot
except ImportError:
    pass
else:
    def _figure(title, width, height, **kwargs):
        return plotting.figure(title=title,
                               tools="save",
                               plot_width=width if width is not None else 700,
                               plot_height=height if height is not None else 700,
                               **kwargs)


    def _add_production_envelope(plot, envelope, key, patch_color="#99d8c9", patch_alpha=0.3, line_color="blue"):
        ub = envelope["objective_upper_bound"].values
        lb = envelope["objective_lower_bound"].values
        var = envelope[key].values

        x = [v for v in var] + [v for v in reversed(var)]
        y = [v for v in lb] + [v for v in reversed(ub)]

        plot.patch(x=x, y=y, color=patch_color, alpha=patch_alpha)

        if "label" in envelope.columns:
            plot.text(envelope["label"].values, var, ub)

        plot.line(var, ub, color=line_color)
        plot.line(var, lb, color=line_color)
        if ub[-1] != lb[-1]:
            plot.line((var[-1], var[-1]), (ub[-1], lb[-1]), color=line_color)


    def plot_flux_variability_analysis_bokeh(fva_result, grid=None, width=None, height=None, title=None,
                                             axis_font_size=None, color="blue"):

        title = "Flux Variability Analysis" if title is None else title

        factors = list(fva_result.index)
        x0 = list(fva_result.lower_bound.values)
        x1 = list(fva_result.upper_bound.values)
        min_value = min([min(x0), min(x1)])
        max_value = max([max(x0), max(x1)])

        p = _figure(title, width, height, y_range=factors, x_range=[min_value - 5, max_value + 5])

        line_width = (width - 30) / len(factors) / 2

        p.segment(x0, factors, x1, factors, line_width=line_width, line_color=color)
        for x_0, x_1, f in zip(x0, x1, factors):
            if x_0 == x_1:
                p.segment([x_0 - 0.01], [f], [x_1 + 0.01], [f], line_width=line_width, color=color)
        p.line([0, 0], [0, len(factors) + 1], line_color="black", line_width=1, line_alpha=0.3)

        if axis_font_size is not None:
            p.xaxis.axis_label_text_font_size = axis_font_size
            p.yaxis.axis_label_text_font_size = axis_font_size

        if grid is not None:
            grid.append(p)
        else:
            plotting.show(p)


    def plot_2_flux_variability_analysis_bokeh(fva_result1, fva_result2, grid=None, width=None, height=None, title=None,
                                               axis_font_size=None, color1="blue", color2="orange"):

        factors = list(fva_result1.index)

        left1 = list(fva_result1.upper_bound)
        right1 = list(fva_result1.lower_bound)
        top1 = [i for i in range(0, len(factors))]
        bottom1 = [i + 0.5 for i in range(0, len(factors))]

        left2 = list(fva_result2.upper_bound)
        right2 = list(fva_result2.lower_bound)
        top2 = [i for i in range(0, len(factors))]
        bottom2 = [i - 0.5 for i in range(0, len(factors))]

        x_range = [min([min(bottom1), min(bottom2)]) - 5, max([max(top1), max(top2)]) + 5]

        title = "Comparing Flux Variability Results" if title is None else title

        p = _figure(title, width, height, x_range=x_range, y_range=factors)

        p.quad(top=top1, bottom=bottom1, left=left1, right=right1, color=color1)
        p.quad(top=top2, bottom=bottom2, left=left2, right=right2, color=color2)

        if axis_font_size is not None:
            p.xaxis.axis_label_text_font_size = axis_font_size
            p.yaxis.axis_label_text_font_size = axis_font_size

        if grid is not None:
            grid.append(p)
        else:
            plotting.show(p)


    def plot_production_envelope_bokeh(envelope, objective, key, grid=None, width=None, height=None, title=None,
                                       points=None, points_colors=None, axis_font_size=None, color="blue"):

        title = "Production Envelope" if title is None  else title
        p = _figure(title, width, height)

        p.xaxis.axis_label = key
        p.yaxis.axis_label = objective

        _add_production_envelope(p, envelope, key, patch_color=color, patch_alpha=0.3, line_color=color)

        if axis_font_size is not None:
            p.xaxis.axis_label_text_font_size = axis_font_size
            p.yaxis.axis_label_text_font_size = axis_font_size

        if points is not None:
            p.scatter(*zip(*points), color="green" if points_colors is None else points_colors)

        if grid is not None:
            grid.append(p)
        else:
            plotting.show(p)


    def plot_2_production_envelopes_bokeh(envelope1, envelope2, objective, key, grid=None, width=None,
                                          height=None, title=None, points=None, points_colors=None,
                                          axis_font_size=None, color1="blue", color2="orange"):

        title = "2 Production Envelopes" if title is None else title
        p = _figure(title, width, height)

        p.xaxis.axis_label = key
        p.yaxis.axis_label = objective

        _add_production_envelope(p, envelope1, key, patch_color=color1, patch_alpha=0.3, line_color=color1)
        _add_production_envelope(p, envelope2, key, patch_color=color2, patch_alpha=0.3, line_color=color2)

        if axis_font_size is not None:
            p.xaxis.axis_label_text_font_size = axis_font_size
            p.yaxis.axis_label_text_font_size = axis_font_size

        if points is not None:
            p.scatter(*zip(*points), color="green" if points_colors is None else points_colors)

        if grid is not None:
            grid.append(p)
        else:
            plotting.show(p)

try:
    from bashplotlib import scatterplot
except ImportError:
    pass
else:
    def plot_production_envelope_cli(envelope, objective, key, grid=None, width=None, height=None, title=None,
                                     points=None, points_colors=None, axis_font_size=None, color="blue"):
        scatterplot.plot_scatter(None, envelope[key], envelope["objective_upper_bound"], "*")


def plot_production_envelope(envelope, objective, key, grid=None, width=None, height=None, title=None,
                             points=None, points_colors=None, axis_font_size=None, color="blue"):
    if width is None and height is None:
        width = 700
    if width is None or height is None:
        width, height = _golden_ratio(width, height)
    try:
        if config.use_bokeh:
            plot_production_envelope_bokeh(envelope, objective, key, grid=grid, width=width, height=height,
                                           title=title, points=points, points_colors=points_colors,
                                           axis_font_size=axis_font_size, color=color)
        else:
            plot_production_envelope_cli(envelope, objective, key, width=width, height=height, title=title,
                                         points=points,
                                         points_colors=points_colors, axis_font_size=axis_font_size, color=color)
    except NameError:
        logger.logger.warn(MISSING_PLOTTING_BACKEND_MESSAGE)


def plot_2_production_envelopes(envelope1, envelope2, objective, key, grid=None, width=None, height=None, title=None,
                                points=None, points_colors=None, axis_font_size=None, color1="blue", color2="orange"):
    if width is None and height is None:
        width = 700
    if width is None or height is None:
        width, height = _golden_ratio(width, height)
    try:
        if config.use_bokeh:
            plot_2_production_envelopes_bokeh(envelope1, envelope2, objective, key, grid=grid, width=width,
                                              height=height,
                                              title=title, points=points, points_colors=points_colors,
                                              axis_font_size=axis_font_size, color1=color1, color2=color2)
        else:
            logger.warn(MISSING_PLOTTING_BACKEND_MESSAGE)
    except NameError:
        logger.warn(MISSING_PLOTTING_BACKEND_MESSAGE)


def plot_flux_variability_analysis(fva_result, grid=None, width=None, height=None, title=None, axis_font_size=None,
                                   color="blue"):
    if width is None and height is None:
        width = 700
    if width is None or height is None:
        width, height = _golden_ratio(width, height)
    try:
        if config.use_bokeh:
            plot_flux_variability_analysis_bokeh(fva_result, grid=grid, width=width, height=height, title=title,
                                                 axis_font_size=axis_font_size, color=color)
        else:
            logger.warn(MISSING_PLOTTING_BACKEND_MESSAGE)
    except NameError:
        logger.warn(MISSING_PLOTTING_BACKEND_MESSAGE)


def plot_2_flux_variability_analysis(fva_result1, fva_result2, grid=None, width=None, height=None, title=None,
                                     axis_font_size=None, color1="blue", color2="orange"):
    if width is None and height is None:
        width = 700
    if width is None or height is None:
        width, height = _golden_ratio(width, height)
    try:
        if config.use_bokeh:
            plot_2_flux_variability_analysis_bokeh(fva_result1, fva_result2, grid=grid, width=width, height=height,
                                                   title=title, axis_font_size=axis_font_size, color1=color1,
                                                   color2=color2)
        else:
            logger.warn(MISSING_PLOTTING_BACKEND_MESSAGE)
    except NameError:
        logger.warn(MISSING_PLOTTING_BACKEND_MESSAGE)


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
        width = int(height + height / GOLDEN_RATIO)

    elif height is None:
        height = int(width / GOLDEN_RATIO)

    return width, height
