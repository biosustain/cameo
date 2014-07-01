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
import multiprocessing
from threading import Thread, current_thread
import threading
from ipython_notebook_utils import ProgressBar


class IPythonMultiprocessObserver():
    def __init__(self, number_of_islands, color_map={}):
        self.progress = {}
        self.connections = {}
        self.clients = {}
        self.color_map = color_map
        self.thread_pool = []
        for i in xrange(number_of_islands):
            label = "<span style='color:%s;'>Island %i </span>" % (self.color_map[i], i+1)
            self.progress[i] = ProgressBar(label=label)
            self.progress[i].start()
            parent_conn, child_conn = multiprocessing.Pipe()
            self.connections[i] = parent_conn
            self.clients[i] = self.Client(child_conn)

    def start(self):
        for i, connection in self.connections.iteritems():
            t = Thread(target=self._listen, args=(connection, self.progress[i]))
            t.start()
            self.thread_pool.append(t)

    def _listen(self, connection, progress):
        while True:
            try:
                res = connection.recv()
                progress.set(res)
            except EOFError as e:
                print "Broken pipe"
                break
            except Exception:
                pass
        threading.exit()

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

        def __call__(self, population, num_generations, num_evaluations, args):
            p = (float(num_evaluations) / float(args.get('max_evaluations', 50000))) * 100.0
            try:
                self.connection.send(p)
            except Exception:
                pass

        def reset(self):
            pass

        def __name__(self):
            return "Multiprocess IPython Notebook progress"

        def close(self):
            self.connection.close()
