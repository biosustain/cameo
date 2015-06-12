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

from __future__ import absolute_import, print_function
from IPython.core.display import display
from IPython.core.display import HTML
from pandas import DataFrame
import re

import numpy as np

from cameo import Metabolite, Model, phenotypic_phase_plane, load_model
from cameo import config, util
from cameo.core.result import Result
from cameo.api.hosts import hosts, Host
from cameo.api.products import products
from cameo.strain_design.heuristic import GeneKnockoutOptimization
from cameo.strain_design.heuristic.objective_functions import biomass_product_coupled_yield
from cameo.ui import notice, bold, loading, stop_loader
from cameo.strain_design import pathway_prediction
from cameo.util import TimeMachine
from cameo.data import universal_models

import logging
from cameo.visualization import visualization
from cameo.visualization.plotting import Grid

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# TODO: implement cplex preference (if available)


class _OptimizationRunner(object):

    def __call__(self, strategy, *args, **kwargs):
        (host, model, pathway) = strategy
        with TimeMachine() as tm:
            pathway.plug_model(model, tm)
            objective = biomass_product_coupled_yield(model.biomass,
                                                      pathway.exchange,
                                                      model.carbon_source)
            opt = GeneKnockoutOptimization(model, objective_function=objective, progress=True, plot=False)
            return opt.run(product=pathway.exchange.id, max_evaluations=10000)


class StrainDesigns(Result):

    def __init__(self, *args, **kwargs):
        super(StrainDesigns, self).__init__(*args, **kwargs)


class Designer(object):
    """High-level strain design functionality.

    Example
    -------
    design = Designer()
    designs = design(product='L-glutamate')
    """

    def __init__(self):
        """"""
        pass

    def __call__(self, product='L-glutamate', hosts=hosts, database=None, view=config.default_view):
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

        Returns
        -------
        Designs
        """
        if database is None:
            database = universal_models.metanetx_universal_model_bigg_rhea

        notice("Starting searching for compound %s" % product)
        product = self.__translate_product_to_universal_reactions_model_metabolite(product, database)
        pathways = self.predict_pathways(product, hosts=hosts, database=database)
        optimization_reports = self.optimize_strains(pathways, view)
        return pathways

    def optimize_strains(self, pathways, view):
        runner = _OptimizationRunner()
        return view.apply(runner, [(host, model, pathway) for (host, model) in pathways
                                   for pathway in pathways[host, model]])


    def predict_pathways(self, product, hosts=None, database=None):  #TODO: make this work with a single host or model
        """Predict production routes for a desired product and host spectrum.
        Parameters
        ----------
        product : str or Metabolite
            The desired product.
        hosts : list or Model or Host
            A list of hosts (e.g. cameo.api.hosts), models, mixture thereof, or a single model or host.

        Returns
        -------
        dict
            ...
        """
        pathways = dict()

        product = self.__translate_product_to_universal_reactions_model_metabolite(product, database)
        for host in hosts:
            if isinstance(host, Model):
                host = Host(name='UNKNOWN_HOST', models=[host])
            for model in list(host.models):
                notice('Predicting pathways for product %s in %s (using model %s).'
                       % (product.name, host, model.id))
                identifier = loading()
                try:
                    logger.debug('Trying to set solver to cplex for pathway predictions.')
                    model.solver = 'cplex'  # CPLEX is better predicting pathways
                except ValueError:
                    logger.debug('Could not set solver to cplex for pathway predictions.')
                    pass
                pathway_predictor = pathway_prediction.PathwayPredictor(model,
                                                                        universal_model=database,
                                                                        compartment_regexp=re.compile(".*_c$"))
                # TODO adjust these numbers to something reasonable
                predicted_pathways = pathway_predictor.run(product, max_predictions=4, timeout=3*60)
                pathways[(host, model)] = predicted_pathways
                stop_loader(identifier)
                self.__display_pathways_information(predicted_pathways, host, model)
        return pathways

    def __translate_product_to_universal_reactions_model_metabolite(self, product, database):
        if isinstance(product, Metabolite):
            return product
        elif isinstance(product, str):
            search_result = products.search(product)
            notice("Found %d compounds that match query '%s'" % (len(search_result), product))
            self.__display_product_search_result(search_result)
            notice("Choosing best match (%s) ... please interrupt if this is not the desired compound."
                   % search_result.name[0])
            self.__display_compound(search_result.name[0], search_result.InChI[0])
            return database.metabolites.get_by_id(search_result.index[0])

    @staticmethod
    def __display_product_search_result(search_result):
        if util.in_ipnb():
            Designer.__display_product_search_results_html(search_result)
        else:
            Designer.__display_product_search_results_cli(search_result)

    @staticmethod
    def __display_compound(name, inchi):
        if util.in_ipnb():
            Designer.__display_compound_html(name, inchi)
        else:
            Designer.__display_compound_cli(name, inchi)

    @staticmethod
    def __display_compound_html(name, inchi):
        svg = Designer.__generate_svg(inchi)
        display(HTML("""
        <p><strong>%s</strong><br/>
            %s
        </p>
        """ % (name, svg)))

    @staticmethod
    def __display_compound_cli(name, inchi):
        text = Designer.__generate_ascii(inchi)
        bold(name)
        print(text)

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
        # TODO: remove copy hack.
        model = load_model(original_model._id)
        with Grid(nrows=2, title="Production envelopes for %s (%s)" % (host.name, model.id)) as grid:
            for i, pathway in enumerate(predicted_pathways):
                pathway_id = "Pathway %i" % (i+1)
                bold(pathway_id)
                display((pathway.data_frame()))
                with TimeMachine() as tm:
                    pathway.plug_model(model, tm)
                    production_envelope = phenotypic_phase_plane(model, variables=[pathway.exchange])
                    production_envelope.plot(grid, title=pathway_id, width=320, height=300)


design = Designer()