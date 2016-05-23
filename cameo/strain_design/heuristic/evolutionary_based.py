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

import logging
import inspyred
import numpy

from pandas import DataFrame

from cameo.exceptions import SolveError
from cameo.visualization.plotting import plotter
from cameo.util import ProblemCache, TimeMachine

from cameo.flux_analysis.simulation import fba
from cameo.flux_analysis.analysis import phenotypic_phase_plane

from cameo.core.solver_based_model import SolverBasedModel

from cameo.strain_design.strain_design import StrainDesignMethod, StrainDesignResult, StrainDesign

from cameo.strain_design.heuristic.evolutionary.archives import ProductionStrainArchive
from cameo.strain_design.heuristic.evolutionary.objective_functions import biomass_product_coupled_min_yield, \
    biomass_product_coupled_yield
from cameo.strain_design.heuristic.evolutionary.optimization import GeneKnockoutOptimization, \
    ReactionKnockoutOptimization
from cameo.strain_design.heuristic.evolutionary.processing import process_knockout_solution

from IProgress.progressbar import ProgressBar
from IProgress.widgets import Bar, Percentage

__all__ = ["OptGene"]


logger = logging.getLogger(__name__)


class OptGene(StrainDesignMethod):
    def __init__(self, model, evolutionary_algorithm=inspyred.ec.GA, manipulation_type="genes", essential_genes=None,
                 essential_reactions=None, plot=True, exclude_non_gene_reactions=True, seed=None, *args, **kwargs):
        if not isinstance(model, SolverBasedModel):
            raise TypeError("Argument 'model' should be of type 'cameo.core.SolverBasedModel'.")

        super(OptGene, self).__init__(*args, **kwargs)

        if exclude_non_gene_reactions:
            essential_reactions = essential_reactions or []
            essential_reactions += [r for r in model.reactions if not r.genes]

        self._model = model
        self._algorithm = evolutionary_algorithm
        self._optimization_algorithm = None
        self._manipulation_type = None
        self._essential_genes = essential_genes
        self._essential_reactions = essential_reactions
        self._plot = plot
        self._seed = None
        self.manipulation_type = manipulation_type

    @property
    def seed(self):
        return self._seed

    @seed.setter
    def seed(self, seed):
        self._seed = seed
        self._optimization_algorithm.seed = seed

    @property
    def manipulation_type(self):
        return self._manipulation_type

    @property
    def plot(self):
        return self._plot

    @plot.setter
    def plot(self, plot):
        self._plot = plot
        if self._optimization_algorithm is not None:
            self._optimization_algorithm.plot = plot

    @manipulation_type.setter
    def manipulation_type(self, manipulation_type):
        self._manipulation_type = manipulation_type
        if manipulation_type is "genes":
            self._optimization_algorithm = GeneKnockoutOptimization(
                model=self._model,
                heuristic_method=self._algorithm,
                essential_genes=self._essential_genes,
                plot=self.plot,
                seed=self._seed)
        elif manipulation_type is "reactions":
            self._optimization_algorithm = ReactionKnockoutOptimization(
                model=self._model,
                heuristic_method=self._algorithm,
                essential_reactions=self._essential_reactions,
                plot=self.plot,
                seed=self._seed)
        else:
            raise ValueError("Invalid manipulation type %s" % manipulation_type)

    def run(self, target=None, biomass=None, substrate=None, max_knockouts=5, simulation_method=fba,
            growth_coupled=False, max_evaluations=20000, population_size=200, time_machine=None,
            max_results=50, **kwargs):
        """
        Parameters
        ----------
        target: str, Metabolite or Reaction
            The design target
        biomass: str, Metabolite or Reaction
            The biomass definition in the model
        substrate: str, Metabolite or Reaction
            The main carbon source
        max_knockouts: int
            Max number of knockouts allowed
        simulation_method: function
            Any method from cameo.flux_analysis.simulation or equivalent
        robust: bool
            If true will use the minimum flux rate to compute the fitness
        max_evaluations: int
            Number of evaluations before stop
        population_size: int
            Number of individuals in each generation
        time_machine: TimeMachine
            See TimeMachine
        max_results: int
            Max number of different designs to return if found
        kwargs: dict
            Arguments for the simulation method.


        Returns
        -------
        OptGeneResult
        """

        target = self._model.reaction_for(target, time_machine=time_machine)
        biomass = self._model.reaction_for(biomass, time_machine=time_machine)
        substrate = self._model.reaction_for(substrate, time_machine=time_machine)

        if growth_coupled:
            objective_function = biomass_product_coupled_min_yield(biomass, target, substrate)
        else:
            objective_function = biomass_product_coupled_yield(biomass, target, substrate)
        self._optimization_algorithm.objective_function = objective_function
        self._optimization_algorithm.simulation_kwargs = kwargs
        self._optimization_algorithm.simulation_method = simulation_method
        self._optimization_algorithm.archiver = ProductionStrainArchive()
        result = self._optimization_algorithm.run(max_evaluations=max_evaluations,
                                                  pop=population_size,
                                                  max_size=max_knockouts,
                                                  maximize=True,
                                                  max_archive_size=max_results,
                                                  **kwargs)

        kwargs.update(self._optimization_algorithm.simulation_kwargs)

        return OptGeneResult(self._model, result, objective_function, simulation_method, self._manipulation_type,
                             biomass, target, substrate, kwargs)


class OptGeneResult(StrainDesignResult):
    __method_name__ = "OptGene"

    __aggregation_function = {
        "genes": lambda x: tuple(tuple(e for e in elements) for elements in x.values)
    }

    def __init__(self, model, designs, objective_function, simulation_method, manipulation_type,
                 biomass, target, substrate, simulation_kwargs, *args, **kwargs):
        super(OptGeneResult, self).__init__(*args, **kwargs)
        assert isinstance(model, SolverBasedModel)

        self._model = model
        self._designs = designs
        self._objective_function = objective_function
        self._simulation_method = simulation_method
        self._manipulation_type = manipulation_type
        self._biomass = biomass
        self._target = target
        self._substrate = substrate
        self._processed_solutions = None
        self._simulation_kwargs = simulation_kwargs

    def __iter__(self):
        for solution in self._designs:
            yield StrainDesign(knockouts=solution, manipulation_type=self._manipulation_type)

    def __len__(self):
        if self._processed_solutions is None:
            self._process_solutions()
        return len(self._processed_solutions)

    def _repr_html_(self):
        return """
        <h3>OptGene Result</h3>
        <ul>
            <li>Simulation: %s<br/></li>
            <li>Objective Function: %s<br/></li>
        </ul>
        %s
        """ % (self._simulation_method.__name__,
               self._objective_function._repr_latex_(),
               self.data_frame._repr_html_())

    @property
    def data_frame(self):
        if self._processed_solutions is None:
            self._process_solutions()

        if self._manipulation_type == "reactions":
            data_frame = DataFrame(self._processed_solutions)
        else:
            columns = self._processed_solutions.columns.difference(["reactions", "size"])
            aggregation_functions = {k: self.__aggregation_function.get(k, lambda x: x.values[0]) for k in columns}
            data_frame = self._processed_solutions.groupby(["reactions", "size"], as_index=False) \
                .aggregate(aggregation_functions)
            data_frame = data_frame[self._processed_solutions.columns]

        data_frame.sort_values("size", inplace=True)
        data_frame.index = [i for i in range(len(data_frame))]
        return data_frame

    def _process_solutions(self):
        processed_solutions = DataFrame(columns=["reactions", "genes", "size", "fva_min", "fva_max",
                                                 "target_flux", "biomass_flux", "yield", "fitness"])

        cache = ProblemCache(self._model)
        progress = ProgressBar(maxval=len(self._designs), widgets=["Processing solutions: ", Bar(), Percentage()])
        for i, solution in progress(enumerate(self._designs)):
            try:
                processed_solutions.loc[i] = process_knockout_solution(
                    self._model, solution, self._simulation_method, self._simulation_kwargs, self._biomass,
                    self._target, self._substrate, [self._objective_function], cache=cache)
            except SolveError as e:
                logger.error(e)
                processed_solutions.loc[i] = [numpy.nan for _ in processed_solutions.columns]

        if self._manipulation_type == "reactions":
            processed_solutions.drop('genes', axis=1, inplace=True)

        self._processed_solutions = processed_solutions

    def display_on_map(self, index=0, map_name=None, palette="YlGnBu"):
        with TimeMachine() as tm:
            for ko in self.data_frame.loc[index, "reactions"]:
                self._model.reactions.get_by_id(ko).knock_out(tm)
            fluxes = self._simulation_method(self._model, **self._simulation_kwargs)
            fluxes.display_on_map(map_name=map_name, palette=palette)

    def plot(self, index=0, grid=None, width=None, height=None, title=None, palette=None, **kwargs):
        wt_production = phenotypic_phase_plane(self._model, objective=self._target, variables=[self._biomass])
        with TimeMachine() as tm:
            for ko in self.data_frame.loc[index, "reactions"]:
                self._model.reactions.get_by_id(ko).knock_out(tm)
            mt_production = phenotypic_phase_plane(self._model, objective=self._target, variables=[self._biomass])

        if title is None:
            title = "Production Envelope"

        dataframe = DataFrame(columns=["ub", "lb", "value", "strain"])
        for _, row in wt_production.iterrows():
            _df = DataFrame([[row['objective_upper_bound'], row['objective_lower_bound'], row[self._biomass.id], "WT"]],
                            columns=dataframe.columns)
            dataframe = dataframe.append(_df)
        for _, row in mt_production.iterrows():
            _df = DataFrame([[row['objective_upper_bound'], row['objective_lower_bound'], row[self._biomass.id], "MT"]],
                            columns=dataframe.columns)
            dataframe = dataframe.append(_df)

        plot = plotter.production_envelope(dataframe, grid=grid, width=width, height=height, title=title,
                                           x_axis_label=self._biomass.id, y_axis_label=self._target.id, palette=palette)
        plotter.display(plot)
