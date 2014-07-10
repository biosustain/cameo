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
import inspyred

from pandas.core.common import in_ipnb
from cameo import util
from cameo import config
from cameo import parallel
from cameo.flux_analysis.simulation import pfba
from cameo.strain_design.heuristic import ReactionKnockoutOptimization, GeneKnockoutOptimization, \
    KnockoutOptimizationResult
from cameo.strain_design.heuristic.multiprocess.observers import IPythonNotebookMultiprocessProgressObserver, \
    CliMultiprocessProgressObserver
from cameo.strain_design.heuristic.multiprocess.plotters import IPythonNotebookBokehMultiprocessPlotObserver
from cameo.strain_design.heuristic.multiprocess.migrators import MultiprocessingMigrator


def run_island(island_class, init_kwargs, clients, migrator, run_kwargs):
    island = island_class(**init_kwargs)
    island.migrator = migrator
    island.observer = clients
    return island.run(**run_kwargs)


class MultiprocessHeuristicOptimization(object):

    _island_class = None

    def __init__(self, model=None, objective_function=None, heuristic_method=inspyred.ec.GA,
                 view=config.default_view, number_of_islands=4, max_migrants=1, *args, **kwargs):
        super(MultiprocessHeuristicOptimization, self).__init__(*args, **kwargs)
        self.model = model
        self.objective_function = objective_function
        self.heuristic_method = heuristic_method
        self.number_of_islands = number_of_islands
        self.view = view
        self.color_map = util.generate_colors(number_of_islands)
        self.migrator = MultiprocessingMigrator(max_migrants)
        self.observers = []

    def _init_kwargs(self):
        return {
            'model': self.model,
            'objective_function': self.objective_function,
            'heuristic_method': self.heuristic_method
        }

    def run(self, **run_kwargs):
        run_kwargs['view'] = parallel.SequentialView()

        results = []
        try:
            async_res = []
            for i in xrange(self.number_of_islands):
                clients = [o.clients[i] for o in self.observers]
                res = self.view.apply_async(run_island,
                                            self._island_class,
                                            self._init_kwargs(),
                                            clients,
                                            self.migrator,
                                            run_kwargs)
                async_res.append(res)

            results = [res.get() for res in async_res]

        except KeyboardInterrupt as e:
            self.view.shutdown()
            raise e

        return results


class MultiprocessKnockoutOptimization(MultiprocessHeuristicOptimization):
    def __init__(self, simulation_method=pfba, *args, **kwargs):
        super(MultiprocessKnockoutOptimization, self).__init__(*args, **kwargs)
        self.simulation_method = simulation_method
        self.observers = self._set_observers()

    def _init_kwargs(self):
        init_kwargs = MultiprocessHeuristicOptimization._init_kwargs(self)
        init_kwargs['simulation_method'] = self.simulation_method
        return init_kwargs

    def _set_observers(self):
        observers = []
        progress_observer = None
        plotting_observer = None
        if in_ipnb():
            progress_observer = IPythonNotebookMultiprocessProgressObserver(number_of_islands=self.number_of_islands,
                                                                            color_map=self.color_map)
            if config.use_bokeh:
                plotting_observer = IPythonNotebookBokehMultiprocessPlotObserver(number_of_islands=self.number_of_islands,
                                                                                 color_map=self.color_map)
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

        results = MultiprocessHeuristicOptimization.run(self, **kwargs)

        for observer in self.observers:
            observer.finish()

        return map(KnockoutOptimizationResult.merge, results)


class MultiprocessReactionKnockoutOptimization(MultiprocessKnockoutOptimization):

    _island_class = ReactionKnockoutOptimization

    def __init__(self, reactions=None, essential_reactions=None, *args, **kwargs):
        super(MultiprocessReactionKnockoutOptimization, self).__init__(*args, **kwargs)
        if reactions is None:
            self.reactions = set([r.id for r in self.model.reactions])
        else:
            self.reactions = reactions

        if essential_reactions is None:
            self.essential_reactions = set([r.id for r in self.model.essential_reactions()])
        else:
            self.essential_reactions = essential_reactions

    def _init_kwargs(self):
        init_kwargs = MultiprocessKnockoutOptimization._init_kwargs(self)
        init_kwargs['essential_reactions'] = self.essential_reactions
        init_kwargs['reactions'] = self.reactions
        return init_kwargs


class MultiprocessGeneKnockoutOptimization(MultiprocessKnockoutOptimization):

    _island_class = GeneKnockoutOptimization

    def __init__(self, genes=None, essential_genes=None, *args, **kwargs):
        super(MultiprocessGeneKnockoutOptimization, self).__init__(*args, **kwargs)
        if genes is None:
            self.genes = set([g.id for g in self.model.genes])
        else:
            self.genes = genes

        if essential_genes is None:
            self.essential_genes = set([g.id for g in self.model.essential_genes()])
        else:
            self.essential_genes = essential_genes

    def _init_kwargs(self):
        init_kwargs = MultiprocessKnockoutOptimization._init_kwargs(self)
        init_kwargs['essential_genes'] = self.essential_genes
        init_kwargs['genes'] = self.genes
        return init_kwargs