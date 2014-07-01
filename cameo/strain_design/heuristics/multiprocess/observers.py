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
import threading
from ipython_notebook_utils import ProgressBar


class IPythonMultiprocessObserver():
    def __init__(self, number_of_islands):
        self.progress = {}
        self.connections = {}
        self.clients = {}
        self.thread_pool = []
        for i in xrange(number_of_islands):
            self.progress[i] = ProgressBar(label="Island %i " % (i+1))
            self.progress[i].start()
            parent_conn, child_conn = multiprocessing.Pipe()
            self.connections[i] = parent_conn
            self.clients[i] = self.Client(child_conn)

    def start(self):
        for i, connection in self.connections.iteritems():
            t = threading.Thread(target=self._listen, args=(connection, self.progress[i]))
            t.start()
            self.thread_pool.append(t)

    def _listen(self, connection, progress):
        while True:
            try:
                res = connection.recv()
                progress.set(res)
            except EOFError, e:
                print "Broken pipe: %s" % e.message
                break
            except Exception, e:
                print "Something went wrong %s" % e.message
        connection.close()

    class Client():
        def __init__(self, connection):
            self.connection = connection

        def __call__(self, population, num_generations, num_evaluations, args):
            if num_evaluations % args.get('n', 1) == 0:
                p = (float(num_evaluations) / float(args.get('max_evaluations', 50000))) * 100.0
                self.connection.send(p)

        def reset(self):
            pass

        def __name__(self):
            return "Multiprocess IPython Notebook progress"

        def close(self):
            self.connection.close()
