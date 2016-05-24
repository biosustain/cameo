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
from __future__ import absolute_import

from math import ceil
import six

from warnings import warn
from ggplot import scale_colour_manual, geom_area, geom_tile, scale_x_continuous, scale_y_continuous, aes, facet_grid

from cameo.util import in_ipnb, inheritdocstring
from cameo.visualization.plotting import AbstractPlotter


@six.add_metaclass(inheritdocstring)
class GGPlotPlotter(AbstractPlotter):
    def __init__(self, **options):
        warn("ggplot interface is under construction...")

        super(GGPlotPlotter, self).__init__(**options)

    def production_envelope(self, dataframe, grid=None, width=None, height=None, title=None, points=None,
                            points_colors=None, palette=None, x_axis_label=None, y_axis_label=None):

        palette = self.get_option('palette') if palette is None else palette
        width = self.get_option('width') if width is None else width
        colors = self._palette(palette, len(dataframe.strain.unique()))

        plot = aes(data=dataframe, ymin="lb", ymax="ub", x="value", color=scale_colour_manual(colors)) + geom_area()
        if title:
            plot += geom_tile(title)
        if x_axis_label:
            plot += scale_x_continuous(name=x_axis_label)
        if y_axis_label:
            plot += scale_y_continuous(name=y_axis_label)

        return plot

    def flux_variability_analysis(self, dataframe, grid=None, width=None, height=None, title=None, palette=None,
                                  x_axis_label=None, y_axis_label=None):

        return aes(data=dataframe, )

    @property
    def _display(self):
        if in_ipnb():
            from IPython.display import display
            return display

    @staticmethod
    def _make_grid(grid):
        columns = ceil(grid.n_rows / len(grid.plots()))
        return grid.plot[0] + facet_grid(grid.n_rows, columns, scales="fixed")
