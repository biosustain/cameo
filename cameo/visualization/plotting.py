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

    def plot_production_envelope_ipython_matplotlib(envelope, key, highligt=[], grid=None,
                                                    width=None, height=None, title=None):
        pass
        # plt.plot(envelope["objective_upper_bound"], envelope[key], title="Production envelop")
        # plt.xlabel("growth")
        # plt.ylabel(key)

except ImportError:

    def plot_production_envelope_ipython_matplotlib(envelope, key, highligt=[], grid=None,
                                                    width=None, height=None, title=None):
        pass

try:
    from bokeh import plotting
    from bokeh.models import GridPlot

    def plot_production_envelope_ipython_bokeh(envelope, key, highligt=[], grid=None,
                                               width=None, height=None, title=None):
        p = plotting.figure(title=title if title is not None else "Production envelope",
                            tools="",
                            plot_width=width if width is not None else 700,
                            plot_height=height if height is not None else 700)
        p.xaxis.axis_label = 'growth'
        p.yaxis.axis_label = key
        colors = ["blue" if i in highligt else "grey" for i in envelope.index]
        if "label" in envelope.columns:
            p.text(envelope["label"], envelope["objective_upper_bound"], envelope[key])
        p.circle(envelope["objective_upper_bound"], envelope[key], color=colors, size=5)
        p.line(envelope["objective_upper_bound"], envelope[key], color="green")
        if grid is not None:
            grid.append(p)
        else:
            plotting.show(p)

except ImportError:

    def plot_production_envelope_ipython_bokeh(envelope, key, highligt=[], grid=None,
                                               width=None, height=None, title=None):
        pass

try:
    from bashplotlib import scatterplot

    def plot_production_envelope_cli(envelope, key, highligt=[], grid=None, width=None, height=None, title=None):
        scatterplot.plot_scatter(None, envelope["objective_upper_bound"], envelope[key], "*")

except ImportError:
    def plot_production_envelope_cli(envelope, key, highligt=[], grid=None, width=None, height=None, title=None):
        pass


def plot_production_envelope(envelope, key, highligt=[], grid=None, width=None, height=None, title=None):
    if util.in_ipnb():
        if config.use_bokeh:
            plot_production_envelope_ipython_bokeh(envelope, key, highligt=highligt, grid=grid,
                                                   width=width, height=height, title=title)
        elif config.use_matplotlib:
            plot_production_envelope_ipython_matplotlib(envelope, key, highligt=highligt, grid=grid,
                                                        width=width, height=height, title=title)
    else:
        plot_production_envelope_cli(envelope, key, highligt=highligt, width=width, height=height, title=title)


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

