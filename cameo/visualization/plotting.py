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

__all__ = ['plot_production_envelope']

from cameo import config, util

try:
    # import matplotlib.pyplot as plt

    def plot_production_envelope_ipython_matplotlib(envelope, objective, key, grid=None, width=None, height=None,
                                                    title=None, points=None, points_colors=None, axis_font_size=None):
        pass
        # plt.plot(envelope["objective_upper_bound"], envelope[key], title="Production envelop")
        # plt.xlabel("growth")
        # plt.ylabel(key)

except:

    def plot_production_envelope_ipython_matplotlib(envelope, objective, key, grid=None, width=None, height=None,
                                                    title=None, points=None, points_colors=None, axis_font_size=None):
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
            p.scatter(*points, color="green" if points_colors is None else points_colors)

        if grid is not None:
            grid.append(p)
        else:
            plotting.show(p)

except:

    def plot_production_envelope_ipython_bokeh(envelope, objective, key, grid=None, width=None, height=None,
                                               title=None, points=None, points_colors=None, axis_font_size=None):
        pass

try:
    from bashplotlib import scatterplot

    def plot_production_envelope_cli(envelope, objective, key, grid=None, width=None, height=None,
                                     title=None, points=None, points_colors=None, axis_font_size=None):
        scatterplot.plot_scatter(None, envelope[key], envelope["objective_upper_bound"], "*")

except:
    def plot_production_envelope_cli(envelope, objective, key, grid=None, width=None, height=None,
                                     title=None, points=None, points_colors=None, axis_font_size=None):
        pass


def plot_production_envelope(envelope, objective, key, grid=None, width=None, height=None, title=None,
                             points=None, points_colors=None, axis_font_size=None):
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
