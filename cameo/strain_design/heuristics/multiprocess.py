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
from multiprocessing import Process, Pool
from cameo import config
import inspyred
from cameo.strain_design.heuristics.observers import IPythonMultiprocessObserver


class MultiprocessHeuristicOptimization(object):
    def __init__(self, view=config.default_view, islands=[], number_of_migrators=1):
        self.view = view
        self.islands = islands
        self.number_of_migrators = number_of_migrators

    def run(self, *args, **kwargs):
        migrator = inspyred.ec.migrators.MultiprocessingMigrator(self.number_of_migrators)
        jobs = []
        i = 1
        observer = IPythonMultiprocessObserver(len(self.islands))
        for island in self.islands:
            island.observer = [observer]
            island.heuristic_method.migrator = migrator
            p = Process(target=self._run_island, args=(island, config.SequentialView(), i, kwargs))
            p.start()
            jobs.append(p)
            i += 1

        for job in jobs:
            print job.is_alive()
            job.join()

        for island in self.islands:
            island.heuristic_method.migrator = None

    def _run_island(self, island, view, i, kwargs):
        print "Launching island %i with **%s" % (i, kwargs)
        try:
            island.run(view=view, i=i, evaluate_migrant=False, **kwargs)
        except Exception, e:
            print e

        print "Run finished"
