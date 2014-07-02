# Copyright 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import colorsys
from multiprocessing import Pipe
from threading import Thread, current_thread

from uuid import uuid1
from IPython.core.display import display
from bokeh.plotting import *
from pandas import DataFrame
from pandas.core.common import in_ipnb


class IslandsPlotObserver(object):
    def __init__(self, url='default', number_of_islands=None, color_map={}):
        self.url = url
        self.iterations = []
        self.fitness = []
        self.uuid = None
        self.in_ipnb = in_ipnb()
        self.plotted = False
        self.color_map = {}
        self.connections = {}
        self.clients = {}
        self.color_map = color_map
        for i in xrange(number_of_islands):
            parent_conn, child_conn = Pipe()
            self.connections[i] = parent_conn
            self.clients[i] = self.Client(child_conn)

        self.data_frame = DataFrame(columns=['iteration', 'island', 'color', 'fitness'])

    def _set_plot(self):
        if self.in_ipnb:
            self.uuid = uuid1()
            output_notebook(url=self.url, docname=str(self.uuid))
            figure()
            scatter([], [], title="Best solution convergence plot", tools='', x_axis_label="Iteration",
                    y_axis_label="Fitness", color=self.color_map, fill_alpha=0.2, size=7)

            self.plot = curplot()
            renderer = [r for r in self.plot.renderers if isinstance(r, Glyph)][0]
            self.ds = renderer.data_source
            show()

    def _listen(self, connection, index, n):
        while True:
            try:
                message = connection.recv()
                df = DataFrame({'iteration': [message['iteration']], 'fitness': [message['fitness']],
                                'color': [self.color_map[index]], 'island': index})
                self.data_frame = self.data_frame.append(df, ignore_index=True)
                if message['iteration'] % n == 0:
                    self._ipnb_plot()
            except EOFError:
                print "Broken pipe"
                break
            except Exception:
                pass
        current_thread().exit()

    def _ipnb_plot(self):
        if self.in_ipnb:
            self.ds.data['x'] = self.data_frame['iteration']
            self.ds.data['y'] = self.data_frame['fitness']
            self.ds.data['fill_color'] = self.data_frame['color']
            self.ds.data['line_color'] = self.data_frame['color']
            self.ds.data['colors'] = self.data_frame['color']
            self.ds._dirty = True
            session().store_obj(self.ds)

    def start(self, n):
        self._set_plot()
        for i, connection in self.connections.iteritems():
            t = Thread(target=self._listen, args=(connection, i, n))
            t.start()

    def reset(self):
        self.data_frame = DataFrame(columns=['iteration', 'island', 'color', 'fitness'])
        self.plotted = False

    def close(self):
        for i, connection in self.connections.iteritems():
            client = self.clients[i]
            try:
                client.close()
            except Exception as e:
                print e
            try:
                connection.close()
            except Exception as e:
                print e


    class Client():
        def __init__(self, connection):
            self.connection = connection
            self.iteration = 0

        def __call__(self, population, num_generations, num_evaluations, args):
            self.iteration += 1
            best = max(population)
            try:
                self.connection.send({'fitness': best.fitness, 'iteration': self.iteration})
            except Exception:
                pass

        def __name__(self):
            return "Fitness Plot Client"

        def reset(self):
            self.iteration = 0

        def close(self):
            self.connection.close()