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
from multiprocessing import Process

import inspyred
from pandas.core.common import in_ipnb

from cameo import config
from cameo.strain_design.heuristics.multiprocess.observers import IPythonNotebookMultiprocessProgressObserver, \
    CliMultiprocessProgressObserver
from cameo.strain_design.heuristics.multiprocess.plotters import IPythonNotebookBokehMultiprocessPlotObserver


class MultiprocessHeuristicOptimization(object):
    def __init__(self, view=config.default_view, islands=[], number_of_migrators=1):
        self.view = view
        self.islands = islands
        self.number_of_migrators = number_of_migrators
        self.color_map = {}
        self._generate_colors(len(islands))

    def _generate_colors(self, n):
        hsv_tuples = [(v*1.0/n, 0.5, 0.5) for v in xrange(n)]
        for i in xrange(n):
            color = colorsys.hsv_to_rgb(*hsv_tuples[i])
            color = tuple(map(lambda v: v*256, color))
            self.color_map[i] = '#%02x%02x%02x' % color

    def run(self, *args, **kwargs):
        migrator = inspyred.ec.migrators.MultiprocessingMigrator(self.number_of_migrators)
        jobs = []
        i = 0
        if in_ipnb():
            progress_observer = IPythonNotebookMultiprocessProgressObserver(number_of_islands=len(self.islands),
                                                                            color_map=self.color_map)
            plotting_observer = IPythonNotebookBokehMultiprocessPlotObserver(number_of_islands=len(self.islands),
                                                                             color_map=self.color_map,
                                                                             n=kwargs.get('n', 20))
            progress_observer.start()
            plotting_observer.start()
        else:
            progress_observer = CliMultiprocessProgressObserver(number_of_islands=len(self.islands))
            plotting_observer = None
            progress_observer.start()
            plotting_observer.start()

        for island in self.islands:
            island.observer = []
            if not plotting_observer is None:
                island.observer.append(plotting_observer.clients[i])
            if not progress_observer is None:
                island.observer.append(progress_observer.clients[i])
            island.heuristic_method.migrator = migrator
            p = Process(target=self._run_island, args=(island, config.SequentialView(), i, kwargs))
            p.start()
            jobs.append(p)
            i += 1
        try:
            for job in jobs:
                job.join()
        except KeyboardInterrupt:
            for job in jobs:
                job.terminate()

        finally:
            plotting_observer.close()
            progress_observer.stop()

    def _run_island(self, island, view, i, kwargs):
        island.run(view=view, i=i, evaluate_migrant=False, **kwargs)