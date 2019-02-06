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
from IProgress.progressbar import ProgressBar
from IProgress.widgets import Bar, Percentage
from pandas import DataFrame

from cobra import Model
from cameo.core.strain_design import StrainDesignMethod, StrainDesignMethodResult, StrainDesign
from cameo.core.target import ReactionKnockoutTarget, GeneKnockoutTarget, ReactionCofactorSwapTarget
from cameo.core.manipulation import swap_cofactors
from cobra.exceptions import OptimizationError
from cameo.flux_analysis.analysis import phenotypic_phase_plane
from cameo.flux_analysis.simulation import fba
from cameo.strain_design.heuristic.evolutionary.archives import ProductionStrainArchive
from cameo.strain_design.heuristic.evolutionary.objective_functions import biomass_product_coupled_min_yield, \
    biomass_product_coupled_yield
from cameo.strain_design.heuristic.evolutionary.optimization import GeneKnockoutOptimization, \
    ReactionKnockoutOptimization, CofactorSwapOptimization, NADH_NADPH
from cameo.strain_design.heuristic.evolutionary.processing import process_reaction_knockout_solution, \
    process_gene_knockout_solution, process_reaction_swap_solution
from cameo.util import TimeMachine
from cameo.visualization.plotting import plotter
from cameo.core.utils import get_reaction_for

__all__ = ["OptGene"]


logger = logging.getLogger(__name__)


class OptGene(StrainDesignMethod):
    def __init__(self, model, evolutionary_algorithm=inspyred.ec.GA, manipulation_type="genes", essential_genes=None,
                 essential_reactions=None, plot=True, exclude_non_gene_reactions=True, *args, **kwargs):
        if not isinstance(model, Model):
            raise TypeError("Argument 'model' should be of type 'cobra.Model'.")

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
        self._manipulation_type = manipulation_type

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

    def run(self, target=None, biomass=None, substrate=None, max_knockouts=5, variable_size=True,
            simulation_method=fba, growth_coupled=False, max_evaluations=20000, population_size=200,
            max_results=50, use_nullspace_simplification=True, seed=None, **kwargs):
        """
        Parameters
        ----------
        target : str, Metabolite or Reaction
            The design target
        biomass : str, Metabolite or Reaction
            The biomass definition in the model
        substrate : str, Metabolite or Reaction
            The main carbon source
        max_knockouts : int
            Max number of knockouts allowed
        variable_size : bool
            If true, all candidates have the same size. Otherwise the candidate size can be from 1 to max_knockouts.
        simulation_method : function
            Any method from cameo.flux_analysis.simulation or equivalent
        growth_coupled : bool
            If true will use the minimum flux rate to compute the fitness
        max_evaluations : int
            Number of evaluations before stop
        population_size : int
            Number of individuals in each generation
        max_results : int
            Max number of different designs to return if found.
        kwargs : dict
            Arguments for the simulation method.
        seed : int
            A seed for random.
        use_nullspace_simplification : Boolean (default True)
            Use a basis for the nullspace to find groups of reactions whose fluxes are multiples of each other and dead
            end reactions. From each of these groups only 1 reaction will be included as a possible knockout.


        Returns
        -------
        OptGeneResult
        """

        target = get_reaction_for(self._model, target)
        biomass = get_reaction_for(self._model, biomass)
        substrate = get_reaction_for(self._model, substrate)

        if growth_coupled:
            objective_function = biomass_product_coupled_min_yield(biomass, target, substrate)
        else:
            objective_function = biomass_product_coupled_yield(biomass, target, substrate)
        if self.manipulation_type == "genes":
            optimization_algorithm = GeneKnockoutOptimization(
                model=self._model,
                heuristic_method=self._algorithm,
                essential_genes=self._essential_genes,
                plot=self.plot,
                objective_function=objective_function,
                use_nullspace_simplification=use_nullspace_simplification)
        elif self.manipulation_type == "reactions":
            optimization_algorithm = ReactionKnockoutOptimization(
                model=self._model,
                heuristic_method=self._algorithm,
                essential_reactions=self._essential_reactions,
                plot=self.plot,
                objective_function=objective_function,
                use_nullspace_simplification=use_nullspace_simplification)
        else:
            raise ValueError("Invalid manipulation type %s" % self.manipulation_type)
        optimization_algorithm.simulation_kwargs = kwargs
        optimization_algorithm.simulation_method = simulation_method
        optimization_algorithm.archiver = ProductionStrainArchive()

        result = optimization_algorithm.run(max_evaluations=max_evaluations,
                                            pop_size=population_size,
                                            max_size=max_knockouts,
                                            variable_size=variable_size,
                                            maximize=True,
                                            max_archive_size=max_results,
                                            seed=seed,
                                            **kwargs)

        kwargs.update(optimization_algorithm.simulation_kwargs)

        return OptGeneResult(self._model, result, objective_function, simulation_method, self.manipulation_type,
                             biomass, target, substrate, kwargs)


class OptGeneResult(StrainDesignMethodResult):
    __method_name__ = "OptGene"

    __aggregation_function = {
        "genes": lambda x: tuple(tuple(e for e in elements) for elements in x.values)
    }

    def __init__(self, model, knockouts, objective_function, simulation_method, manipulation_type,
                 biomass, target, substrate, simulation_kwargs, *args, **kwargs):
        super(OptGeneResult, self).__init__(self._generate_designs(knockouts, manipulation_type), *args, **kwargs)
        assert isinstance(model, Model)

        self._model = model
        self._knockouts = knockouts
        self._objective_function = objective_function
        self._simulation_method = simulation_method
        self._manipulation_type = manipulation_type
        self._biomass = biomass
        self._target = target
        self._substrate = substrate
        self._processed_solutions = None
        self._simulation_kwargs = simulation_kwargs

    @staticmethod
    def _generate_designs(knockouts, manipulation_type):
        designs = []
        if manipulation_type == "reactions":
            target_class = ReactionKnockoutTarget
        elif manipulation_type == "genes":
            target_class = GeneKnockoutTarget
        else:
            raise ValueError("Invalid 'manipulation_type' %s" % manipulation_type)

        for knockout_design, _ in knockouts:
            designs.append(StrainDesign([target_class(ko) for ko in knockout_design]))

        return designs

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
        if self._manipulation_type == "reactions":
            self._process_reaction_knockout_solutions()
        elif self._manipulation_type == "genes":
            self._process_gene_knockout_solutions()

    def _process_gene_knockout_solutions(self):
        processed_solutions = DataFrame(columns=["reactions", "genes", "size", "fva_min", "fva_max",
                                                 "target_flux", "biomass_flux", "yield", "fitness"])

        if len(self._knockouts) == 0:
            logger.warn("No solutions found")
            self._processed_solutions = processed_solutions

        else:
            progress = ProgressBar(maxval=len(self._knockouts), widgets=["Processing solutions: ", Bar(), Percentage()])
            for i, solution in progress(enumerate(self._knockouts)):
                try:
                    processed_solutions.loc[i] = process_gene_knockout_solution(
                        self._model, solution[0], self._simulation_method, self._simulation_kwargs, self._biomass,
                        self._target, self._substrate, self._objective_function)
                except OptimizationError as e:
                    logger.error(e)
                    processed_solutions.loc[i] = [numpy.nan for _ in processed_solutions.columns]

            self._processed_solutions = processed_solutions

    def _process_reaction_knockout_solutions(self):
        processed_solutions = DataFrame(columns=["reactions", "size", "fva_min", "fva_max",
                                                 "target_flux", "biomass_flux", "yield", "fitness"])

        if len(self._knockouts) == 0:
            logger.warn("No solutions found")
            self._processed_solutions = processed_solutions

        else:
            progress = ProgressBar(maxval=len(self._knockouts), widgets=["Processing solutions: ", Bar(), Percentage()])
            for i, solution in progress(enumerate(self._knockouts)):
                try:
                    processed_solutions.loc[i] = process_reaction_knockout_solution(
                        self._model, solution[0], self._simulation_method, self._simulation_kwargs, self._biomass,
                        self._target, self._substrate, self._objective_function)
                except OptimizationError as e:
                    logger.error(e)
                    processed_solutions.loc[i] = [numpy.nan for _ in processed_solutions.columns]

            self._processed_solutions = processed_solutions

    def display_on_map(self, index=0, map_name=None, palette="YlGnBu"):
        with self._model:
            for ko in self.data_frame.loc[index, "reactions"]:
                self._model.reactions.get_by_id(ko).knock_out()
            fluxes = self._simulation_method(self._model, **self._simulation_kwargs)
            fluxes.display_on_map(map_name=map_name, palette=palette)

    def plot(self, index=0, grid=None, width=None, height=None, title=None, palette=None, **kwargs):
        wt_production = phenotypic_phase_plane(self._model, objective=self._target, variables=[self._biomass])
        with self._model:
            for ko in self.data_frame.loc[index, "reactions"]:
                self._model.reactions.get_by_id(ko).knock_out()
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


class HeuristicOptSwap(StrainDesignMethod):
    def __init__(self, model, evolutionary_algorithm=inspyred.ec.GA, plot=True, cofactor_id_swaps=NADH_NADPH,
                 exclude_non_gene_reactions=True, *args, **kwargs):
        super(HeuristicOptSwap, self).__init__(*args, **kwargs)
        self._skip_reactions = []
        if exclude_non_gene_reactions:
            self._skip_reactions += [r for r in model.reactions if not r.genes]

        self._algorithm = evolutionary_algorithm
        self._swap_pairs = cofactor_id_swaps
        self._optimization_algorithm = None
        self._model = self._optimization_algorithm.model
        self._plot = plot

    def run(self, target=None, biomass=None, substrate=None, max_swaps=5, variable_size=True,
            simulation_method=fba, growth_coupled=False, max_evaluations=20000, population_size=200,
            time_machine=None, max_results=50, seed=None, **kwargs):
        """
        Parameters
        ----------
        target : str, Metabolite or Reaction
            The design target.
        biomass : str, Metabolite or Reaction
            The biomass definition in the model.
        substrate : str, Metabolite or Reaction
            The main carbon source.
        max_swaps : int
            Max number of swaps allowed.
        variable_size : bool
            If true, all candidates have the same size. Otherwise the candidate size can be from 1 to max_knockouts.
        simulation_method : function
            Any method from cameo.flux_analysis.simulation or equivalent.
        growth_coupled : bool
            If true will use the minimum flux rate to compute the fitness.
        max_evaluations : int
            Number of evaluations before stop.
        population_size : int
            Number of individuals in each generation.
        time_machine : TimeMachine
            See TimeMachine.
        max_results : int
            Max number of different designs to return if found.
        kwargs : dict
            Arguments for the simulation method.
        seed : int
            A seed for random.


        Returns
        -------
        HeuristicOptSwapResult
        """

        target = get_reaction_for(self._model, target)
        biomass = get_reaction_for(self._model, biomass)
        substrate = get_reaction_for(self._model, substrate)

        if growth_coupled:
            objective_function = biomass_product_coupled_min_yield(biomass, target, substrate)
        else:
            objective_function = biomass_product_coupled_yield(biomass, target, substrate)

        optimization_algorithm = CofactorSwapOptimization(model=self._model,
                                                          cofactor_id_swaps=self._cofactor_id_swaps,
                                                          skip_reactions=self._skip_reactions,
                                                          objective_function=objective_function)

        optimization_algorithm.simulation_kwargs = kwargs
        optimization_algorithm.simulation_method = simulation_method
        optimization_algorithm.archiver = ProductionStrainArchive()

        result = optimization_algorithm.run(max_evaluations=max_evaluations,
                                            pop_size=population_size,
                                            max_size=max_swaps,
                                            variable_size=variable_size,
                                            maximize=True,
                                            max_archive_size=max_results,
                                            seed=seed,
                                            **kwargs)

        kwargs.update(optimization_algorithm.simulation_kwargs)

        return HeuristicOptSwapResult(self._model, result, self._swap_pairs, objective_function,
                                      simulation_method, biomass, target, substrate, kwargs)


class HeuristicOptSwapResult(StrainDesignMethodResult):
    __method_name__ = "HeuristicOptSwap"

    def __init__(self, model, swaps, swap_pairs, objective_function, simulation_method, biomass, target,
                 substrate, simulation_kwargs, *args, **kwargs):
        super(HeuristicOptSwapResult, self).__init__(self._generate_designs(swaps, swap_pairs), *args, **kwargs)
        assert isinstance(model, Model)

        self._model = model
        self._swaps = swaps
        self._swap_pairs = swap_pairs
        self._objective_function = objective_function
        self._simulation_method = simulation_method
        self._biomass = biomass
        self._target = target
        self._substrate = substrate
        self._processed_solutions = None
        self._simulation_kwargs = simulation_kwargs

    @staticmethod
    def _generate_designs(swaps, swap_pair):
        designs = []
        for swap_design, _ in swaps:
            designs.append(StrainDesign([ReactionCofactorSwapTarget(swap, swap_pair) for swap in swap_design]))

        return designs

    def _repr_html_(self):
        return """
        <h3>OptSwap Result</h3>
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

        data_frame = DataFrame(self._processed_solutions)

        data_frame.sort_values("size", inplace=True)
        data_frame.index = [i for i in range(len(data_frame))]
        return data_frame

    def _process_solutions(self):
        processed_solutions = DataFrame(columns=["reactions", "size", "fva_min", "fva_max",
                                                 "target_flux", "biomass_flux", "yield", "fitness"])

        if len(self._swaps) == 0:
            logger.warn("No solutions found")
            self._processed_solutions = processed_solutions

        else:
            progress = ProgressBar(maxval=len(self._swaps), widgets=["Processing solutions: ", Bar(), Percentage()])
            for i, solution in progress(enumerate(self._swaps)):
                try:
                    processed_solutions.loc[i] = process_reaction_swap_solution(
                        self._model, solution[0], self._simulation_method, self._simulation_kwargs, self._biomass,
                        self._target, self._substrate, self._objective_function, self._swap_pairs)
                except OptimizationError as e:
                    logger.error(e)
                    processed_solutions.loc[i] = [numpy.nan for _ in processed_solutions.columns]

            self._processed_solutions = processed_solutions

    def display_on_map(self, index=0, map_name=None, palette="YlGnBu"):
        with self._model:
            for ko in self.data_frame.loc[index, "reactions"]:
                swap_cofactors(self._model.reactions.get_by_id(ko), self._model, self._swap_pairs)
            fluxes = self._simulation_method(self._model, **self._simulation_kwargs)
            fluxes.display_on_map(map_name=map_name, palette=palette)

    def plot(self, index=0, grid=None, width=None, height=None, title=None, palette=None, **kwargs):
        wt_production = phenotypic_phase_plane(self._model, objective=self._target, variables=[self._biomass])
        with self._model:
            for ko in self.data_frame.loc[index, "reactions"]:
                swap_cofactors(self._model.reactions.get_by_id(ko), self._model, self._swap_pairs)
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
