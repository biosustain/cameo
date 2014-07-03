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
from cameo import parallel

import inspyred
from pandas.core.common import in_ipnb
from cameo import util
from cameo import config
from cameo.flux_analysis.simulation import pfba
from cameo.strain_design.heuristics import HeuristicOptimization, ReactionKnockoutOptimization, GeneKnockoutOptimization
from cameo.strain_design.heuristics.multiprocess.observers import IPythonNotebookMultiprocessProgressObserver, \
    CliMultiprocessProgressObserver
from cameo.strain_design.heuristics.multiprocess.plotters import IPythonNotebookBokehMultiprocessPlotObserver


class MultiprocessRunner(object):
    def __init__(self, kwargs):
        self._kwargs = kwargs

    def __call__(self, island):
        return island.evolve(**self._kwargs)


class MultiprocessHeuristicOptimization(HeuristicOptimization):
    def __init__(self, view=config.default_view, number_of_islands=4, number_of_migrators=1, *args, **kwargs):
        self.islands = []
        self.number_of_islands = number_of_islands
        super(MultiprocessHeuristicOptimization, self).__init__(*args, **kwargs)
        self.view = view
        self.color_map = util.generate_colors(number_of_islands)
        self._migrator = inspyred.ec.migrators.MultiprocessingMigrator(number_of_migrators)
        self._build_islands(*args, **kwargs)

    def _build_islands(self, *args, **kwargs):
        raise NotImplementedError

    def run(self,  **kwargs):
        kwargs['view'] = parallel.SequentialView()
        runner = MultiprocessRunner(kwargs)
        try:
            results = self.view.map(runner, self.islands)
        except KeyboardInterrupt as e:
            self.view.shutdown()
            raise e

        return results

    @HeuristicOptimization.heuristic_method.setter
    def heuristic_method(self, heuristic_method):
        for island in self.islands:
            island.heuristic_method = heuristic_method

    @HeuristicOptimization.objective_function.setter
    def objective_function(self, objective_function):
        for island in self.islands:
            island.objective_function = objective_function


class MultiprocessKnockoutOptimization(MultiprocessHeuristicOptimization):
    def __init__(self, simulation_method=pfba, *args, **kwargs):
        self.simulation_method = simulation_method
        super(MultiprocessKnockoutOptimization, self).__init__(*args, **kwargs)
        self.observers = self._set_observers(**kwargs)

    def _set_observers(self, **kwargs):
        observers = []
        progress_observer = None
        plotting_observer = None
        if in_ipnb():
            progress_observer = IPythonNotebookMultiprocessProgressObserver(number_of_islands=self.number_of_islands,
                                                                            color_map=self.color_map)
            if config.use_bokeh:
                plotting_observer = IPythonNotebookBokehMultiprocessPlotObserver(number_of_islands=self.number_of_islands,
                                                                                 color_map=self.color_map,
                                                                                 n=kwargs.get('n', 20))
            elif config.use_matplotlib:
                pass
        else:
            progress_observer = CliMultiprocessProgressObserver(number_of_islands=self.number_of_islands)

        if not progress_observer is None:
            observers.append(progress_observer)
        if not plotting_observer is None:
            observers.append(plotting_observer)

        return observers

    def run(self,  **kwargs):
        for observer in self.observers:
            observer.start()

        for i in xrange(self.number_of_islands):
            for observer in self.observers:
                self.islands[i].observers = observer.clients[i]

        results = MultiprocessHeuristicOptimization.run(self, **kwargs)

        for observer in self.observers:
            observer.finish()

        return results


class MultiprocessReactionKnockoutOptimization(MultiprocessKnockoutOptimization):
    def __init__(self, reactions=None, essential_reactions=None, *args, **kwargs):
        self.essential_reactions = essential_reactions
        self.reactions = reactions
        super(MultiprocessReactionKnockoutOptimization, self).__init__(*args, **kwargs)

    def _build_islands(self, *args, **kwargs):
        for i in xrange(self.number_of_islands):
            self.islands.append(ReactionKnockoutOptimization(simulation_method=self.simulation_method,
                                                             reactions=self.reactions,
                                                             essential_reactions=self.essential_reactions,
                                                             *args, **kwargs))
            self.islands[i].migrator = self._migrator
            if self.reactions is None:
                self.reactions = self.islands[i].reactions
                self.essential_reactions = self.islands[i].essential_reactions


class MultiprocessGeneKnockoutOptimization(MultiprocessKnockoutOptimization):
    def __init__(self, genes=None, essential_genes=None, *args, **kwargs):
        self.essential_genes = essential_genes
        self.genes = genes
        super(MultiprocessGeneKnockoutOptimization, self).__init__(*args, **kwargs)

    def _build_islands(self, *args, **kwargs):
        for i in xrange(self.number_of_islands):
            self.islands.append(GeneKnockoutOptimization(simulation_method=self.simulation_method,
                                                         genes=self.genes,
                                                         essential_genes=self.essential_genes,
                                                         *args, **kwargs))
            self.islands[i].migrator = self._migrator
            if self.genes is None:
                self.genes = self.islands[i].genes
                self.essential_genes = self.islands[i].essential_genes
