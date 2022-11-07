from __future__ import absolute_import

import math

import matplotlib.pyplot as plt
import seaborn as sns
from cameo.util import zip_repeat, inheritdocstring, partition
from cameo.visualization.plotting.abstract import AbstractPlotter


class SeabornPlotter(AbstractPlotter, metaclass=inheritdocstring):
    class Figure(object):
        def __init__(self, data=None, layout=None):
            self.data = data
            self.layout = layout

    class Layout(object):
        def __init__(self, title=None, xaxis=None, yaxis=None, width=None, height=None):
            self.title = title
            self.xaxis = xaxis
            self.yaxis = yaxis
            self.width = width
            self.height = height

    class Rectangle(object):
        def __init__(
            self,
            x0,
            x1,
            y0,
            y1,
            line_color=None,
            line_width=1,
            fill_color=None,
            fill_alpha=None,
        ):
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
                "type": "rect",
                "x0": self.x0,
                "y0": self.y0,
                "x1": self.x1,
                "y1": self.y1,
                "line": {
                    "width": self.line_width,
                },
            }
            if self.line_color:
                res["line"]["color"] = self.line_color
            if self.fill_alpha:
                res["opacity"] = self.fill_alpha
            if self.fill_color:
                res["fillcolor"] = self.fill_color
            return res

    def __init__(self, **options):
        super(SeabornPlotter, self).__init__(**options)

    def _make_production_envelope(self, dataframe, variable, color=None, ax=None):
        alpha = self.get_option("alpha")
        ub = dataframe["ub"].values.tolist()
        lb = dataframe["lb"].values.tolist()
        var = dataframe["value"].values.tolist()

        x1 = [v for v in reversed(var)]
        y1 = [v for v in reversed(ub)]
        x = [v for v in var] + x1
        y = [v for v in lb] + y1

        if lb[0] != ub[0]:
            x.extend([var[0], var[0]])
            y.extend([lb[0], ub[0]])

        ax = sns.lineplot(
            x=x1,
            y=y1,
            errorbar=None,
            ax=ax,
            alpha=alpha,
        )
        ax.fill_between(
            x=x,
            y1=y,
            color=color,
            alpha=alpha,
        )
        return ax

    def production_envelope(
        self,
        dataframe,
        grid=None,
        width=None,
        height=None,
        title=None,
        points=None,
        points_colors=None,
        palette="RdYlBu",
        x_axis_label=None,
        y_axis_label=None,
    ):
        variables = dataframe["strain"].unique()
        palette = self.get_option("palette") if palette is None else palette
        width = self.get_option("width") if width is None else width

        width, height = self.golden_ratio(width, height)

        sns.set_theme()
        ax = None
        palette = self._palette(palette, len(variables))
        for variable, color in zip_repeat(variables, palette):
            _dataframe = dataframe[dataframe["strain"] == variable]
            ax = self._make_production_envelope(
                _dataframe, variable, color=color, ax=ax
            )

        if points is not None:
            x, y = zip(*points)
            color = "green" if points_colors is None else points_colors
            ax = sns.Scatter(
                x=x,
                y=y,
                color=color,
                ax=ax,
            )

        ax.set_title(title)
        ax.set_xlabel(xlabel=x_axis_label)
        ax.set_ylabel(ylabel=y_axis_label)
        # ax.set_size_inches(width, height)

        if grid is not None:
            grid.append(ax)
            return grid
        return ax

    """
    def _make_production_envelope_3d(self, dataframe, variable, color=None):
        ub_data = dataframe.pivot("value1", "value2", "ub")

        surface = go.Surface(
            x=ub_data.index.tolist(),
            y=ub_data.columns.tolist(),
            z=ub_data.as_matrix(),
            name=variable,
            hoverinfo="none",
            surfacecolor=color,
        )

        return surface

    def production_envelope_3d(
        self,
        dataframe,
        grid=None,
        width=None,
        height=None,
        title=None,
        points=None,
        points_colors=None,
        palette=None,
        x_axis_label=None,
        y_axis_label=None,
        z_axis_label=None,
    ):

        variables = dataframe["strain"].unique()
        palette = self.get_option("palette") if palette is None else palette
        width = self.get_option("width") if width is None else width

        width, height = self.golden_ratio(width, height)
        data = []
        palette = self._palette(palette, len(variables))
        for variable, color in zip_repeat(variables, palette):
            _dataframe = dataframe[dataframe["strain"] == variable]
            surface = self._make_production_envelope_3d(
                _dataframe, variable, color=color
            )
            data.append(surface)

        if points is not None:
            x, y, z = zip(*points)
            scatter = go.Scatter3d(
                x=x,
                y=y,
                z=z,
                mode="markers",
                name="Data Points",
                marker=dict(color="green" if points_colors is None else points_colors),
            )
            data.append(scatter)

        layout = go.Layout(
            title=title,
            scene=go.Scene(
                xaxis=dict(title=x_axis_label),
                yaxis=dict(title=y_axis_label),
                zaxis=dict(title=z_axis_label),
            ),
            width=width,
            height=height,
        )

        if grid is not None:
            plot = self.Figure(data=data, layout=layout)
            grid.append(plot)
            return grid
        else:
            plot = go.Figure(data=data, layout=layout)
        return plot

    def _make_fva_bars(self, factores, dataframe, height, step, color, variable):
        alpha = self.get_option("alpha")
        ub = dataframe["ub"].values
        lb = dataframe["lb"].values

        rectangles = []
        scatter = go.Scatter(
            x=[(x0 + x1) / 2.0 for x0, x1 in zip(ub, lb)],
            y=[(y + 1 + height) for y in range(len(factores))],
            name=variable,
            fillcolor=color,
            mode="markers",
            opacity=alpha,
            hoverinfo="none",
        )
        for x0, x1, y in zip(ub, lb, range(len(factores))):
            y_ = y + 1
            rect = self.Rectangle(
                x0,
                x1,
                y_ + height - step / 2,
                y_ + height + step / 2,
                fill_color=color,
                line_color=color,
                line_width=0,
                fill_alpha=alpha,
            )
            rectangles.append(rect)

        return scatter, rectangles

    def flux_variability_analysis(
        self,
        dataframe,
        grid=None,
        width=None,
        height=None,
        title=None,
        palette=None,
        x_axis_label=None,
        y_axis_label=None,
    ):
        palette = self.get_option("palette") if palette is None else palette
        width = self.get_option("width") if width is None else width

        width, height = self.golden_ratio(width, height)

        variables = dataframe["strain"].unique()
        factors = dataframe["reaction"].unique().tolist()
        data = []
        shapes = []
        n = len(variables)
        step = 1.0 / float(len(variables))
        for variable, i, color in zip(variables, range(n), self._palette(palette, n)):
            _dataframe = dataframe[dataframe["strain"] == variable]
            scatter_, shapes_ = self._make_fva_bars(
                factors, _dataframe, i * step, step, color, variable
            )
            data.append(scatter_)
            shapes += [s.to_dict() for s in shapes_]

        layout = go.Layout(
            title=title,
            xaxis=dict(title=x_axis_label),
            yaxis=dict(
                title=y_axis_label,
                ticktext=[""] + factors,
                tickvals=[i for i in range(len(factors) + 1)],
            ),
            width=width,
            height=height,
            shapes=shapes,
        )

        if grid is not None:
            plot = self.Figure(data=data, layout=layout)
            grid.append(plot)
            return grid
        else:
            plot = go.Figure(data=data, layout=layout)
        return plot

    def line(
        self,
        dataframe,
        width=None,
        height=None,
        palette=None,
        title="Line",
        x_axis_label=None,
        y_axis_label=None,
        grid=None,
    ):

        palette = self.get_option("palette") if palette is None else palette
        width = self.get_option("width") if width is None else width

        width, height = self.golden_ratio(width, height)

        traces = []

        for i, color in zip(
            dataframe.index, self._palette(palette, len(dataframe.index))
        ):
            y = dataframe.loc[i]
            x = list(range(len(y)))
            traces.append(go.Scatter(x=x, y=y, mode="lines", name=i))

        layout = go.Layout(
            title=title,
            xaxis=dict(title=x_axis_label),
            yaxis=dict(title=y_axis_label),
            width=width,
            height=height,
        )

        if grid is not None:
            plot = self.Figure(data=traces, layout=layout)
            grid.append(plot)
            return grid
        else:
            plot = go.Figure(data=traces, layout=layout)
        return plot

    def frequency(
        self,
        dataframe,
        width=None,
        height=None,
        palette=None,
        title="Frequency plot",
        x_axis_label=None,
        y_axis_label="Frequency",
        grid=None,
    ):

        palette = self.get_option("palette") if palette is None else palette
        width = self.get_option("width") if width is None else width

        width, height = self.golden_ratio(width, height)

        bar = go.Bar(
            x=dataframe.index.tolist(),
            y=dataframe["frequency"].tolist(),
            marker=dict(self._palette(palette, len(dataframe.index))),
        )

        layout = go.Layout(
            title=title,
            xaxis=dict(title=x_axis_label),
            yaxis=dict(title=y_axis_label),
            width=width,
            height=height,
        )

        if grid is not None:
            plot = self.Figure(data=[bar], layout=layout)
            grid.append(plot)
            return grid
        else:
            plot = go.Figure(data=[bar], layout=layout)
        return plot
    """

    @property
    def _display(self):
        plt.show()
        return lambda x: x

    @staticmethod
    def _make_grid(grid):
        rows = grid.n_rows
        columns = math.ceil(len(grid.plots) / rows)

        fig = plt.figure(figsize=(grid.width, grid.height))
        gs = fig.add_gridspec(rows, columns)

        actual = 0
        for i in range(rows):
            for j in range(columns):
                ax = fig.add_subplot(gs[i, j])
                ax = grid.plots[actual]
                actual += 1
        return fig
