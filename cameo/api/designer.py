# Copyright 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""This module implements the high-level interface function `design`.
"""

from __future__ import absolute_import, print_function

import numpy as np

from cobra import Metabolite, Model

from cameo.core.strain_design import StrainDesign
from cameo.strain_design.pathway_prediction.pathway_predictor import PathwayResult

try:
    from IPython.core.display import display
    from IPython.core.display import HTML
except ImportError:
    def display(*args, **kwargs):
        print(*args, **kwargs)

    def HTML(*args, **kwargs):
        print(*args, **kwargs)

from pandas import DataFrame


from cameo import fba
from cameo import config, util
from cameo.api.hosts import hosts as HOSTS, Host
from cameo.api.products import products
from cobra.exceptions import OptimizationError
from cameo.strain_design import OptGene, DifferentialFVA
from cameo.ui import notice, searching, stop_loader
from cameo.strain_design import pathway_prediction
from cameo.models import universal
from cameo.strain_design.heuristic.evolutionary.objective_functions import biomass_product_coupled_min_yield
from cameo.strain_design.heuristic.evolutionary.objective_functions import product_yield

from cameo.visualization import visualization

import logging

__all__ = ['design']

logger = logging.getLogger(__name__)


# TODO: implement cplex preference (if available)


class _OptimizationRunner(object):

    def __init__(self, options):
        self.options = options

    def __call__(self, strategy):
        raise NotImplementedError


class _OptGeneRunner(_OptimizationRunner):

    def __call__(self, strategy):
        max_evaluations = 15000
        max_time = (self.options.heuristic_optimization_timeout, 0)

        model, pathway, aerobic = strategy['model'], strategy['pathway'], strategy['aerobic']
        model = model.copy()
        assert isinstance(model, Model)
        assert isinstance(pathway, PathwayResult)
        assert isinstance(aerobic, bool)

        with model:
            if not aerobic and 'EX_o2_e' in model.reactions:
                model.reactions.EX_o2_e.lower_bound = 0
            pathway.apply(model)
            model.objective = model.biomass
            opt_gene = OptGene(model=model, plot=False)
            logger.info('Setting time limit for heuristic optimization to {}'.format(max_time))
            designs = opt_gene.run(target=pathway.product.id, biomass=model.biomass, substrate=model.carbon_source,
                                   max_evaluations=max_evaluations, max_knockouts=15, max_time=max_time)
            logger.debug('Termination cause: {}'.format(designs._knockouts.heuristic_method.termination_cause))

            return designs


class _DifferentialFVARunner(_OptimizationRunner):
    def __call__(self, strategy):

        model, pathway, aerobic = strategy['model'], strategy['pathway'], strategy['aerobic']
        model = model.copy()
        assert isinstance(model, Model)
        assert isinstance(pathway, PathwayResult)
        assert isinstance(aerobic, bool)

        with model:
            if not aerobic and 'EX_o2_e' in model.reactions:
                model.reactions.EX_o2_e.lower_bound = 0

            pathway.apply(model)
            model.objective = model.biomass
            diff_fva = DifferentialFVA(design_space_model=model,
                                       objective=pathway.product.id,
                                       variables=[model.biomass],
                                       points=self.options.differential_fva_points)
            designs = diff_fva.run(improvements_only=True, surface_only=True, progress=False)

            return designs


class DesignerOptions(object):

    def __init__(self):
        self.differential_fva = True
        self.heuristic_optimization = True
        self.max_pathway_predictions = 1
        self.differential_fva_points = 10
        self.pathway_prediction_timeout = 10      # min
        self.heuristic_optimization_timeout = 45  # min
        self.verbose = True


class Designer(object):
    """High-level strain design functionality.

    Example
    -------
    design = Designer()
    designs = design(product='L-glutamate')
    """

    def __init__(self):
        self.options = DesignerOptions()

    def __call__(self, product='L-glutamate',
                 hosts=HOSTS,
                 database=None,
                 aerobic=True,
                 view=config.default_view):
        """The works.

        The following workflow will be followed to determine suitable
        metabolic engineering strategies for a desired product:

        - Determine production pathways for desired product and host organisms.
          Try a list of default hosts if no hosts are specified.
        - Determine maximum theoretical yields and production envelopes for
          all routes.
        - Determine if production routes can be coupled to growth.
        - Determine over-expression, down-regulation, and KO targets.

        Parameters
        ----------
        product : str or Metabolite
            The desired product.
        hosts : list or Model or Host
            A list of hosts (e.g. cameo.api.hosts), models, mixture thereof, or a single model or host.
        aerobic: bool
            If false, sets `model.reaction.EX_o2_e.lower_bound = 0`

        Returns
        -------
        Designs
        """
        if database is None:
            database = universal.metanetx_universal_model_bigg

        notice("Starting searching for compound %s" % product)
        try:
            product = self.translate_product_to_universal_reactions_model_metabolite(product, database)
        except KeyError:
            raise KeyError("Product %s is not in the %s database" % (product, database.id))
        pathways = self.predict_pathways(product, hosts=hosts, database=database, aerobic=aerobic)
        optimization_reports = self.optimize_strains(pathways, view, aerobic=aerobic)
        return optimization_reports

    def optimize_strains(self, pathways, view, aerobic=True):
        """
        Optimize targets for the identified pathways. The optimization will only run if the pathway can be optimized.

        Parameters
        ----------
        pathways : dict
            A dictionary with information of pathways to optimize ([Host, Model] -> PredictedPathways).
        view : object
            A view for multi, single os distributed processing.
        aerobic : bool
            If True, it will set `model.reactions.EX_o2_e.lower_bound` to 0.

        Returns
        -------
        pandas.DataFrame
            A data frame with strain designs processed and ranked.

        """
        strategies = [{'host': host.name, 'model': model, 'pathway': pathway, 'aerobic': aerobic} for (host, model) in
                      pathways for pathway in pathways[host, model] if
                      pathway.needs_optimization(model, objective=model.biomass)]

        print("Optimizing %i pathways" % len(strategies))

        results = DataFrame(columns=["host", "model", "manipulations", "heterologous_pathway",
                                     "fitness", "yield", "product", "biomass", "method"])

        if self.options.heuristic_optimization:
            logger.info('Running heuristic optimizations on predicted pathways.')
            opt_gene_runner = _OptGeneRunner(options=self.options)
            designs = view.map(opt_gene_runner, strategies)
            logger.info('Processing heuristic optimization results.')
            results = self.build_results_data(designs, strategies, results, 'PathwayPredictor+OptGene')
        if self.options.differential_fva:
            logger.info('Running differential flux variability analysis on predicted pathways.')
            differential_fva_runner = _DifferentialFVARunner(options=self.options)
            designs = view.map(differential_fva_runner, strategies)
            logger.info('Processing differential flux variability analysis results.')
            results = self.build_results_data(designs, strategies, results, 'PathwayPredictor+DifferentialFVA')

        return results

    def build_results_data(self, strategy_designs, strategies, results, method_name):
        """
        Process the designs and add them to the `results` DataFrame.

        Parameters
        ----------
        strategy_designs : list
            A list of list[StrainDesign]. Each list corresponds to a strategy.
        strategies : list
            List of [(Host, Model, PathwayResult, Boolean)]. Note: variables: host, model, pathway, anaerobic
        results : pandas.DataFrame
            An existing DataFrame to be extended.

        Returns
        -------
        pandas.DataFrame
            A data frame with the processed designs appended

        """
        for i, strain_designs in enumerate(strategy_designs):
            strategy = strategies[i]
            _results = DataFrame(columns=results.columns, index=[j for j in range(len(strain_designs))])
            manipulations, fitness, yields, target_flux, biomass = \
                self.process_strain_designs(strain_designs, **strategy)
            for j, strain_design in enumerate(manipulations):
                _results.loc[j, 'manipulations'] = strain_design
                _results.loc[j, 'heterologous_pathway'] = strategy['pathway']
            _results['host'] = strategy['host']
            _results['model'] = strategy['model'].id
            _results['fitness'] = fitness
            _results['yield'] = yields
            _results['biomass'] = biomass
            _results['product'] = target_flux
            _results['method'] = method_name  # TODO: determine method name from targets

            results = results.append(_results, ignore_index=True)
        return results

    def process_strain_designs(self, strain_designs, model=None, pathway=None, aerobic=None, **kwargs):
        model = model.copy()
        assert isinstance(pathway, StrainDesign)
        assert isinstance(pathway, PathwayResult)
        final_strain_designs = []
        fitness = []
        yields = []
        biomass = []
        target_flux = []
        pyield = product_yield(pathway.product, model.carbon_source)
        bpcy = biomass_product_coupled_min_yield(model.biomass, pathway.product, model.carbon_source)
        for strain_design in strain_designs:
            assert isinstance(strain_design, StrainDesign)
            _fitness, _yield, _target_flux, _biomass = self.evaluate_design(model, strain_design,
                                                                            pathway, aerobic, bpcy, pyield)
            fitness.append(_fitness)
            yields.append(_yield)
            final_strain_designs.append(strain_design)
            biomass.append(_biomass)
            target_flux.append(_target_flux)

        return final_strain_designs, fitness, yields, target_flux, biomass

    @staticmethod
    def evaluate_design(model, strain_design, pathway, aerobic, bpcy, pyield):
        with model:
            if not aerobic and 'EX_o2_e' in model.reactions:
                model.reactions.EX_o2_e.lower_bound = 0
            pathway.apply(model)
            strain_design.apply(model)
            try:
                solution = fba(model, objective=model.biomass)
                _bpcy = bpcy(model, solution, strain_design.targets)
                _pyield = pyield(model, solution, strain_design.targets)
                target_flux = solution.fluxes[pyield.product]
                biomass = solution.fluxes[bpcy.biomass]
            except OptimizationError:
                _bpcy, _pyield, target_flux, biomass = np.nan, np.nan, np.nan, np.nan
            return _bpcy, _pyield, target_flux, biomass

    def predict_pathways(self, product, hosts, database, aerobic=True):
        """Predict production routes for a desired product and host spectrum.
        Parameters
        ----------
        product : str or Metabolite
            The desired product.
        hosts : list
            A list of hosts (e.g. cameo.api.hosts), models, mixture thereof, or a single model or host.
        database: Model
            A model to use as database. See also: cameo.models.universal
        aerobic: bool
            If True, it will set `model.reactions.EX_o2_e.lower_bound` to 0.

        Returns
        -------
        dict
            ([Host, Model] -> PredictedPathways)
        """
        pathways = dict()
        product = self.translate_product_to_universal_reactions_model_metabolite(product, database)
        for host in hosts:
            logging.debug('Processing host {}'.format(host.name))
            if isinstance(host, Model):
                host = Host(name='UNKNOWN_HOST', models=[host])
            for model in list(host.models):
                with model:
                    if not aerobic and "EX_o2_e" in model.reactions:
                        model.reactions.EX_o2_e.lower_bound = 0
                    identifier = searching()
                    logging.debug('Processing model {} for host {}'.format(model.id, host.name))
                    notice('Predicting pathways for product %s in %s (using model %s).'
                           % (product.name, host, model.id))
                    logging.debug('Predicting pathways for model {}'.format(model.id))
                    pathway_predictor = pathway_prediction.PathwayPredictor(model, universal_model=database)
                    predicted_pathways = pathway_predictor.run(product,
                                                               max_predictions=self.options.max_pathway_predictions,
                                                               timeout=self.options.pathway_prediction_timeout,
                                                               silent=not self.options.verbose)
                    stop_loader(identifier)
                    pathways[(host, model)] = predicted_pathways
                    self.__display_pathways_information(predicted_pathways, host, model)
        return pathways

    def translate_product_to_universal_reactions_model_metabolite(self, product, database):
        if isinstance(product, Metabolite):
            return product
        elif isinstance(product, str):
            search_result = products.search(product)
            search_result = search_result.loc[[i for i in search_result.index if i in database.metabolites]]
            if len(search_result) == 0:
                raise KeyError("No compound matches found for query %s" % product)
            notice("Found %d compounds that match query '%s'" % (len(search_result), product))
            self.__display_product_search_result(search_result)
            notice("Choosing best match (%s) ... please interrupt if this is not the desired compound."
                   % search_result.name.values[0])
            self.__display_compound(search_result.InChI.values[0])
            return database.metabolites.get_by_id(search_result.index[0])

    @staticmethod
    def __display_product_search_result(search_result):
        if util.in_ipnb():
            Designer.__display_product_search_results_html(search_result)
        else:
            Designer.__display_product_search_results_cli(search_result)

    @staticmethod
    def __display_compound(inchi):
        if util.in_ipnb():
            Designer.__display_compound_html(inchi)

    @staticmethod
    def __display_compound_html(inchi):
        svg = Designer.__generate_svg(inchi)
        display(HTML("""
        <p>
            %s
        </p>
        """ % svg))

    @staticmethod
    def __display_product_search_results_html(search_result):
        rows = []
        for index, row in search_result.iterrows():
            name = row["name"]
            formula = row["formula"]
            rows.append("<tr><td>%s</td><td>%s</td><td>%s</td></tr>" % (index, name, formula))

        display(HTML(
            """
            <table>
                <thead>
                    <th>Id</th>
                    <th>Name</th>
                    <th>Formula</th>
                </thead>
                <tbody>
                    %s
                </tbody>
            </table>
            """ % "\n".join(rows)
        ))

    @staticmethod
    def __display_product_search_results_cli(search_result):
        rows = np.ndarray((len(search_result), 3), dtype=object)
        for i, index in enumerate(search_result.index):
            row = search_result.loc[index]
            name = row["name"]
            formula = row["formula"]
            rows[i, ] = [index, name, formula]
            i += 1

        display(DataFrame(rows, columns=["Id", "Name", "Formula"]))

    @staticmethod
    def __generate_svg(inchi):
        if isinstance(inchi, float) or inchi is None:
            return ""
        else:
            return visualization.inchi_to_svg(inchi, three_d=False)

    @staticmethod
    def __generate_ascii(inchi):
        if isinstance(inchi, float) or inchi is None:
            return ""
        else:
            return visualization.inchi_to_ascii(inchi)

    @staticmethod
    def __display_pathways_information(predicted_pathways, host, original_model):
        if util.in_ipnb():
            title = "Production envelopes for %s (%s)" % (host.name, original_model.id)
            predicted_pathways.plot_production_envelopes(original_model, title=title, objective=original_model.biomass)

    @staticmethod
    def calculate_yield(model, source, product):
        try:
            flux_dist = fba(model, objective=product)
            return flux_dist[product.id] / abs(flux_dist[source.id])
        except OptimizationError:
            return 0.0


design = Designer()
