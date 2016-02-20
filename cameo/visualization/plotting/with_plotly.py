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

import math

from plotly import tools
import plotly.graph_objs as go

from cameo.util import zip_repeat, in_ipnb
from cameo.visualization.plotting.abstract import AbstractPlotter, AbstractGrid


class PlotlyGrid(AbstractGrid):
    def plot(self):
        rows = self.n_rows
        columns = math.ceil(len(self.plots)/rows)

        plot = tools.make_subplots(rows=rows, cols=columns, subplot_titles=[p.layout['title'] for p in self.plots])
        plot['layout']['width'] = self.width
        plot['layout']['heigth'] = self.height
        for i, subplot in enumerate(self.plots):
            plot.append_trace(go.Data(data=subplot.data))
            plot['layout']['xaxis%i' % (i+1)].update(subplot.layout['xaxis'])
            plot['layout']['yaxis%i' % (i+1)].update(subplot.layout['yaxis'])
            plot['layout']['shapes%i' % (i+1)].update(subplot.layout['shapes'])

        return plot


class PlotlyPlotter(AbstractPlotter):

    __grid__class__ = PlotlyGrid

    class Figure(object):
        def __init__(self, data=None, layout=None):
            self.data = data
            self.layout = layout

    class Path(object):
        def __init__(self, points, close=True, line_color=None, line_width=1, fill_color=None, fill_alpha=None):
            self.points = points
            self.close = close
            self.line_color = line_color
            self.line_width = line_width
            self.fill_color = fill_color
            self.fill_alpha = fill_alpha

        def to_dict(self):
            res = {
                'type': 'path',
                'path': self.path,
                'line': {
                    'width': self.line_width
                }
            }
            if self.line_color:
                res['line']['color'] = self.line_color
            if self.fill_alpha:
                res['opacity'] = self.fill_alpha
            if self.fill_color:
                res['fillcolor'] = self.fill_color
            return res

        @property
        def path(self):
            first = next(self.points)
            return "M %.3f, %.3f" % first + " ".join("L%.3f,%.3f" % p for p in self.points) + "Z" if self.close else "z"

    class Rectangle(object):
        def __init__(self, x0, x1, y0, y1, line_color=None, line_width=1, fill_color=None, fill_alpha=None):
            self.x0 = x0
            self.x1 = x1
            self.y0 = y0
            self.y1 = y1
            self.line_color = line_color
            self.line_width = line_width
            self.fill_color = fill_color
            self.fill_alpha = fill_alpha

        def to_dict(self):
            res = {
                'type': 'rect',
                'x0': self.x0,
                'y0': self.y0,
                'x1': self.x1,
                'y1': self.y1,
                'line': {
                    'width': self.line_width,
                }
            }
            if self.line_color:
                res['line']['color'] = self.line_color
            if self.fill_alpha:
                res['opacity'] = self.fill_alpha
            if self.fill_color:
                res['fillcolor'] = self.fill_color
            return res

    def __init__(self, **options):
        if 'mode' not in options:
            options['mode'] = 'offline'
        super(PlotlyPlotter, self).__init__(**options)
        if self.get_option('mode') is 'offline' and in_ipnb():
            from plotly.offline import init_notebook_mode
            init_notebook_mode()

    def _make_production_envelope(self, dataframe, variable, color=None):
        alpha = self.get_option('alpha')
        ub = dataframe["ub"].values.tolist()
        lb = dataframe["lb"].values.tolist()
        var = dataframe["value"].values

        x = [v for v in var] + [v for v in reversed(var)]
        y = [v for v in lb] + [v for v in reversed(ub)]

        scatter = go.Scatter(x=x, y=y,
                             mode="line",
                             name=variable,
                             hoverinfo='none',
                             fillcolor=color,
                             opacity=alpha,
                             marker=dict(line=dict(color=color), opacity=alpha))
        path = self.Path(zip(scatter.x, scatter.y),
                         close=True,
                         fill_alpha=alpha,
                         line_color=color,
                         line_width=0,
                         fill_color=color)

        return scatter, path

    def production_envelope(self, dataframe, grid=None, width=None, height=None, title=None, points=None,
                            points_colors=None, palette='RdYlBu', x_axis_label=None, y_axis_label=None):
        variables = dataframe["strain"].unique()
        palette = self.get_option('palette') if palette is None else palette
        width = self.get_option('width') if width is None else width

        width, height = self.golden_ratio(width, height)

        data = []
        shapes = []
        palette = self._palette(palette, len(variables))
        for variable, color in zip_repeat(variables, palette):
            _dataframe = dataframe[dataframe["strain"] == variable]
            d, s = self._make_production_envelope(_dataframe, variable, color=color)
            data.append(d)
            shapes.append(s)

        if points is not None:
            x, y = zip(*points)
            scatter = go.Scatter(x=x, y=y, mode="marker", name="Data Points",
                                 marker=dict(color="green" if points_colors is None else points_colors))
            data.append(scatter)

        layout = go.Layout(
            title=title,
            xaxis=dict(title=x_axis_label),
            yaxis=dict(title=y_axis_label),
            shapes=[s.to_dict() for s in shapes],
            width=width,
            height=height
        )

        if grid is not None:
            plot = self.Figure(data=data, layout=layout)
            grid.append(plot)
            return grid
        else:
            plot = go.Figure(data=data, layout=layout)
        return plot

    def _make_fva_bars(self, factores, dataframe, height, step, color, variable):
        alpha = self.get_option('alpha')
        ub = dataframe['ub'].values
        lb = dataframe['lb'].values

        rectangles = []
        scatter = go.Scatter(
            x=[(x0+x1)/2.0 for x0, x1 in zip(ub, lb)],
            y=[(y+1+height) for y in range(len(factores))],
            name=variable,
            fillcolor=color,
            mode="markers",
            opacity=alpha,
            hoverinfo='none'
        )
        for x0, x1, y in zip(ub, lb, range(len(factores))):
            y_ = y+1
            rect = self.Rectangle(x0, x1, y_+height-step/2, y_+height+step/2,
                                  fill_color=color,
                                  line_color=color,
                                  line_width=0,
                                  fill_alpha=alpha)
            rectangles.append(rect)

        return scatter, rectangles

    def flux_variability_analysis(self, dataframe, grid=None, width=None, height=None, title=None,
                                  palette=None, x_axis_label=None, y_axis_label=None):
        palette = self.get_option('palette') if palette is None else palette
        width = self.get_option('width') if width is None else width

        width, height = self.golden_ratio(width, height)

        variables = dataframe["strain"].unique()
        factors = dataframe['reaction'].unique().tolist()
        data = []
        shapes = []
        n = len(variables)
        step = 1.0/float(len(variables))
        for variable, i, color in zip(variables, range(n), self._palette(palette, n)):
            _dataframe = dataframe[dataframe["strain"] == variable]
            scatter_, shapes_ = self._make_fva_bars(factors, _dataframe, i*step, step, color, variable)
            data.append(scatter_)
            shapes += [s.to_dict() for s in shapes_]

        layout = go.Layout(
            title=title,
            xaxis=dict(title=x_axis_label),
            yaxis=dict(title=y_axis_label, ticktext=[""]+factors, tickvals=[i for i in range(len(factors)+1)]),
            width=width,
            height=height,
            shapes=shapes
        )

        if grid is not None:
            plot = self.Figure(data=data, layout=layout)
            grid.append(plot)
            return grid
        else:
            plot = go.Figure(data=data, layout=layout)
        return plot

    def display(self, plot):
        if self.get_option('mode') is "offline":
            from plotly.offline import iplot
        else:
            from plotly.plotly import iplot

        iplot(plot)


