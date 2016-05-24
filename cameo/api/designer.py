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

import re

import numpy as np

try:
    from IPython.core.display import display
    from IPython.core.display import HTML
except ImportError:
    def display(*args, **kwargs):
        print(*args, **kwargs)

    def HTML(*args, **kwargs):
        print(*args, **kwargs)

from pandas import DataFrame

from cameo import Metabolite, Model, fba
from cameo import config, util
from cameo.core.result import Result
from cameo.api.hosts import hosts, Host
from cameo.api.products import products
from cameo.exceptions import SolveError
from cameo.strain_design.heuristic import OptGene
from cameo.ui import notice, searching, stop_loader
from cameo.strain_design import pathway_prediction
from cameo.util import TimeMachine
from cameo.models import universal

from cameo.visualization import visualization

import logging

__all__ = ['design']

logger = logging.getLogger(__name__)


# TODO: implement cplex preference (if available)


class _OptimizationRunner(object):
    def __call__(self, strategy, *args, **kwargs):
        (model, pathway) = (strategy[1], strategy[2])
        with TimeMachine() as tm:
            pathway.plug_model(model, tm)
            opt_gene = OptGene(model=model, plot=False)
            opt_gene_designs = opt_gene.run(target=pathway.product.id, biomass=model.biomass,
                                            substrate=model.carbon_source, max_evaluations=10000)
            opt_gene_designs._process_solutions()

            # TODO: OptKnock is quite slow. Improve OptKnock (model simplification?)
            # opt_knock = OptKnock(model=model)
            # opt_knock_designs = opt_knock.run(max_knockouts=5, target=pathway.product.id,
            #                                   max_results=5, biomass=model.biomass, substrate=model.carbon_source)

            designs = opt_gene_designs  # + opt_knock_designs

            return designs, pathway


class DesignerResult(Result):
    def __init__(self, designs, *args, **kwargs):
        super(DesignerResult, self).__init__(*args, **kwargs)
        self.designs = designs

    def _repr_latex_(self):
        pass


class StrainDesings(Result):
    def __init__(self, organism, designs, *args, **kwargs):
        super(StrainDesings, self).__init__(*args, **kwargs)


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
            database = universal.metanetx_universal_model_bigg_rhea

        notice("Starting searching for compound %s" % product)
        product = self.__translate_product_to_universal_reactions_model_metabolite(product, database)
        pathways = self.predict_pathways(product, hosts=hosts, database=database)
        print("Optimizing %i pathways" % sum(len(p) for p in pathways.values()))
        optimization_reports = self.optimize_strains(pathways, view)
        return optimization_reports

    @staticmethod
    def optimize_strains(pathways, view):
        """
        """
        runner = _OptimizationRunner()
        designs = [(host, model, pathway) for (host, model) in pathways for pathway in pathways[host, model]
                   if pathway.needs_optimization(model, objective=model.biomass)]
        return view.map(runner, designs)

    def predict_pathways(self, product, hosts=None, database=None):  # TODO: make this work with a single host or model
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
            logging.debug('Processing host {}'.format(host.name))
            if isinstance(host, Model):
                host = Host(name='UNKNOWN_HOST', models=[host])
            for model in list(host.models):
                identifier = searching()
                logging.debug('Processing model {} for host {}'.format(model.id, host.name))
                notice('Predicting pathways for product %s in %s (using model %s).'
                       % (product.name, host, model.id))
                logging.debug('Predicting pathways for model {}'.format(model.id))
                pathway_predictor = pathway_prediction.PathwayPredictor(model,
                                                                        universal_model=database,
                                                                        compartment_regexp=re.compile(".*_c$"))
                # TODO adjust these numbers to something reasonable
                predicted_pathways = pathway_predictor.run(product, max_predictions=4, timeout=3 * 60, silent=True)
                stop_loader(identifier)
                pathways[(host, model)] = predicted_pathways
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
            self.__display_compound(search_result.InChI[0])
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
        else:
            Designer.__display_compound_cli(inchi)

    @staticmethod
    def __display_compound_html(inchi):
        svg = Designer.__generate_svg(inchi)
        display(HTML("""
        <p>
            %s
        </p>
        """ % svg))

    @staticmethod
    def __display_compound_cli(inchi):
        text = Designer.__generate_ascii(inchi)
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
        predicted_pathways.plot_production_envelopes(original_model,
                                                     title="Production envelopes for %s (%s)" % (
                                                         host.name, original_model.id),
                                                     objective=original_model.biomass)

    @staticmethod
    def calculate_yield(model, source, product):
        try:
            flux_dist = fba(model, objective=product)
            return flux_dist[product.id] / abs(flux_dist[source.id])
        except SolveError:
            return 0.0


design = Designer()
