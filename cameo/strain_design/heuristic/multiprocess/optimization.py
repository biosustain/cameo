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

from __future__ import absolute_import, print_function

__all__ = ['MultiprocessReactionKnockoutOptimization', 'MultiprocessGeneKnockoutOptimization']

from six.moves import range
from functools import reduce

import inspyred

from cameo import util
from cameo import config
from cameo import parallel
from cameo.flux_analysis.simulation import pfba
from cameo.strain_design import StrainDesignMethod
from cameo.strain_design.heuristic import ReactionKnockoutOptimization, GeneKnockoutOptimization
from cameo.strain_design.heuristic.optimization import KnockoutOptimizationResult
from cameo.strain_design.heuristic.multiprocess.observers import IPythonNotebookMultiprocessProgressObserver, \
    CliMultiprocessProgressObserver
from cameo.strain_design.heuristic.multiprocess.plotters import IPythonNotebookBokehMultiprocessPlotObserver
from cameo.strain_design.heuristic.multiprocess.migrators import MultiprocessingMigrator


class MultiprocessRunner(object):
    """
    Runner for multiprocessing model. It generates the non-pickable
    objects on the beginning of the process.

    """

    def __init__(self, island_class, init_kwargs, migrator, run_kwargs):
        self.island_class = island_class
        self.init_kwargs = init_kwargs
        self.migrator = migrator
        self.run_kwargs = run_kwargs

    def __call__(self, clients):
        island = self.island_class(**self.init_kwargs)
        island.migrator = self.migrator
        island.observer = clients
        return island.run(**self.run_kwargs)


class MultiprocessHeuristicOptimization(StrainDesignMethod):
    """
    Heuristic Optimization abstract implementation.

    Arguments
    =========

    model: SolverBasedModel
        A model to simulate.
    objective_function: a list of or one objective_function
        The objective for the algorithm to optimize.
    heuristic_method: inspyred.ec instance
        The method using for search (default: inspyred.ec.GA).
    max_migrants: int
        The number of individuals travelling between islands (different processes) at the same time (default: 1).

    """
    _island_class = None

    def __init__(self, model=None, objective_function=None, heuristic_method=inspyred.ec.GA, max_migrants=1, *args,
                 **kwargs):
        super(MultiprocessHeuristicOptimization, self).__init__(*args, **kwargs)
        self.model = model
        self.objective_function = objective_function
        self.heuristic_method = heuristic_method
        self.migrator = MultiprocessingMigrator(max_migrants)
        self.observers = []

    def _init_kwargs(self):
        return {
            'model': self.model,
            'objective_function': self.objective_function,
            'heuristic_method': self.heuristic_method
        }

    def run(self, view=config.default_view, number_of_islands=None, **run_kwargs):
        if number_of_islands is None:
            number_of_islands = len(view)
        run_kwargs['view'] = parallel.SequentialView()
        runner = MultiprocessRunner(self._island_class, self._init_kwargs(), self.migrator, run_kwargs)
        clients = [[o.clients[i] for o in self.observers] for i in range(number_of_islands)]
        try:
            results = view.map(runner, clients)
        except KeyboardInterrupt as e:
            view.shutdown()
            raise e
        return results


class MultiprocessKnockoutOptimization(MultiprocessHeuristicOptimization):
    """
    Heuristic Knockout Optimization Abstract implementation.

    Arguments
    =========

    model: SolverBasedModel
        A model to simulate.
    objective_function: a list of or one objective_function
        The objective for the algorithm to optimize.
    heuristic_method: inspyred.ec instance
        The method using for search (default: inspyred.ec.GA).
    max_migrants: int
        The number of individuals travelling between islands (different processes) at the same time (default: 1).
    simulation_method: a function from flux_analysis.simulation
        The method to simulate the model (default: pfba).
    """

    def __init__(self, simulation_method=pfba, *args, **kwargs):
        super(MultiprocessKnockoutOptimization, self).__init__(*args, **kwargs)
        self.simulation_method = simulation_method

    def _init_kwargs(self):
        init_kwargs = MultiprocessHeuristicOptimization._init_kwargs(self)
        init_kwargs['simulation_method'] = self.simulation_method
        return init_kwargs

    def _set_observers(self, number_of_islands):
        observers = []
        progress_observer = None
        plotting_observer = None
        if util.in_ipnb():
            color_map = util.generate_colors(number_of_islands)
            progress_observer = IPythonNotebookMultiprocessProgressObserver(number_of_islands=number_of_islands,
                                                                            color_map=color_map)
            if config.use_bokeh:
                plotting_observer = IPythonNotebookBokehMultiprocessPlotObserver(number_of_islands=number_of_islands,
                                                                                 color_map=color_map)
            elif config.use_matplotlib:
                pass
        else:
            progress_observer = CliMultiprocessProgressObserver(number_of_islands=number_of_islands)

        if not progress_observer is None:
            observers.append(progress_observer)
        if not plotting_observer is None:
            observers.append(plotting_observer)

        return observers

    def run(self, view=config.default_view, number_of_islands=None, **kwargs):
        if number_of_islands is None:
            number_of_islands = len(view)
        self.observers = self._set_observers(number_of_islands)
        for observer in self.observers:
            observer.start()

        results = MultiprocessHeuristicOptimization.run(self, view=view, number_of_islands=number_of_islands, **kwargs)

        for observer in self.observers:
            observer.finish()

        return reduce(KnockoutOptimizationResult.merge, results)


class MultiprocessReactionKnockoutOptimization(MultiprocessKnockoutOptimization):
    """
    Heuristic Knockout Optimization Abstract implementation.

    Arguments
    =========

    model: SolverBasedModel
        A model to simulate.
    objective_function: a list of or one objective_function
        The objective for the algorithm to optimize.
    heuristic_method: inspyred.ec instance
        The method using for search (default: inspyred.ec.GA).
    max_migrants: int
        The number of individuals travelling between islands (different processes) at the same time (default: 1).
    simulation_method: a function from flux_analysis.simulation
        The method to simulate the model (default: pfba).
    """
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
