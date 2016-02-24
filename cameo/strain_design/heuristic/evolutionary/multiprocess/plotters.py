# Copyright 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.
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

from __future__ import absolute_import, print_function

__all__ = ['IPythonNotebookBokehMultiprocessPlotObserver']

import six.moves.queue
from pandas import DataFrame
from cameo import config
from cameo.strain_design.heuristic.evolutionary.multiprocess.observers import AbstractParallelObserver, \
    AbstractParallelObserverClient

if config.use_bokeh:
    from bokeh.plotting import figure, show
    from bokeh.models import ColumnDataSource
    from bokeh.io import push_notebook


class IPythonNotebookBokehMultiprocessPlotObserver(AbstractParallelObserver):
    __name__ = "IPython Notebook Bokeh Multiprocess Plot Observer"

    def __init__(self, color_map={}, *args, **kwargs):
        super(IPythonNotebookBokehMultiprocessPlotObserver, self).__init__(*args, **kwargs)
        self.connections = {}
        self.color_map = color_map
        self.data_frame = DataFrame(columns=['iteration', 'island', 'color', 'fitness'])
        self.plotted = False

    def _create_client(self, i):
        self.clients[i] = IPythonNotebookBokehMultiprocessPlotObserverClient(queue=self.queue, index=i)

    def start(self):
        self._plot()
        AbstractParallelObserver.start(self)

    def _plot(self):
        self.plot = figure(title="Fitness plot", tools='', plot_height=400, plot_width=650)
        self.plot.xaxis.axis_label = "Iteration"
        self.plot.yaxis.axis_label = "Fitness"
        self.ds = ColumnDataSource(data=dict(x=[], y=[], island=[]))
        self.plot.circle('x', 'y', source=self.ds)

        show(self.plot)
        self.plotted = True

    def _process_message(self, message):
        if not self.plotted:
            self._plot()

        index = message['index']
        df = DataFrame({
            'iteration': [message['iteration']],
            'fitness': [message['fitness']],
            'color': [self.color_map[index]],
            'island': [index]
        })
        self.data_frame = self.data_frame.append(df, ignore_index=True)
        if message['iteration'] % message['n'] == 0:
            self._update_plot()

    def _update_plot(self):
        self.ds.data['x'] = self.data_frame['iteration']
        self.ds.data['y'] = self.data_frame['fitness']
        self.ds.data['fill_color'] = self.data_frame['color']
        self.ds.data['line_color'] = self.data_frame['color']
        self.ds._dirty = True
        push_notebook()

    def stop(self):
        self.data_frame = DataFrame(columns=['iteration', 'island', 'color', 'fitness'])
        self.plotted = False


class IPythonNotebookBokehMultiprocessPlotObserverClient(AbstractParallelObserverClient):
    __name__ = "IPython Notebook Bokeh Multiprocess Plot Observer"

    def __init__(self, *args, **kwargs):
        super(IPythonNotebookBokehMultiprocessPlotObserverClient, self).__init__(*args, **kwargs)
        self.iteration = 0

    def __call__(self, population, num_generations, num_evaluations, args):
        self.iteration += 1
        best = max(population)
        try:
            self._queue.put_nowait({
                'fitness': best.fitness,
                'iteration': self.iteration,
                'index': self.index,
                'n': args.get('n', 1)})
        except six.moves.queue.Full:
            pass

    def reset(self):
        self.iteration = 0

    def close(self):
        self.connection.close()
