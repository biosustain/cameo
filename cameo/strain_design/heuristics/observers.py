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

from ipython_notebook_utils import ProgressBar


class IPythonNotebookObserver():
    def __init__(self):
        self.progress = ProgressBar()

    def __call__(self, population, num_generations, num_evaluations, args):
        if num_evaluations % args.get('n', 1) == 0:
            p = (float(num_evaluations) / float(args.get('max_evaluations', 50000))) * 100.0
            self.progress.set(p)

    def __name__(self):
        return "IPython Notebook progress"

    def reset(self, args):
        self.progress.start()

class IPythonMultiprocessObserver():
    def __init__(self, number_of_islands):
        self.progress = dict([(i, ProgressBar(label="Island %i" % i)) for i in xrange(number_of_islands)])
        self._lock = multiprocessing.Lock()

    def __call__(self, population, num_generations, num_evaluations, args):
        with self._lock:
            i = args.get('i')
            if num_evaluations % args.get('n', 1) == 0:
                p = (float(num_evaluations) / float(args.get('max_evaluations', 50000))) * 100.0
                self.progress[i].set(p)

    def __name__(self):
        return "Multiprocess IPython Notebook progress"

    def reset(self, args):
        i = args.get('i')
        with self._lock:
            self.progress[i].start()
