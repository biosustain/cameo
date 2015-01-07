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

from types import StringType
from functools import partial
import progressbar
from cameo import Metabolite, Model
from cameo.api import _METANETX as METANETX
from cameo.api.hosts import hosts, Host
from cameo.api.products import products
from cameo.strain_design.pathway_prediction import PathwayPredictor
from cameo.util import TimeMachine, DisplayItemsWidget


class Designs(object):
    # TODO: implement results object for Designer
    pass

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

    def __call__(self, product='L-glutamate', hosts=hosts):
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
        product = self.__translate_product_to_universal_reactions_model_metabolite(product)
        pathways = self.predict_pathways(product, hosts=hosts)
        return pathways

    def predict_pathways(self, product, hosts=hosts):  #TODO: make this work with a single host or model
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
        product = self.__translate_product_to_universal_reactions_model_metabolite(product)
        pbar_hosts = progressbar.ProgressBar(widgets=["Processing", DisplayItemsWidget(list(hosts)), progressbar.Bar(), progressbar.ETA(), progressbar.Percentage()])
        for host in pbar_hosts(list(hosts)):
            if isinstance(host, Model):
                host = Host(name='UNKNOWN_HOST', models=[host])
            pbar_models = progressbar.ProgressBar(widgets=["Processing", DisplayItemsWidget([model.id for model in host.models]), progressbar.Bar(), progressbar.ETA(), progressbar.Percentage()])
            for model in pbar_models(list(host.models)):
                # TODO: Check if product is already part of model
                pathway_predictor = PathwayPredictor(model, universal_model=METANETX['universal_model'], mapping=METANETX['bigg2mnx'])
                predicted_pathways = pathway_predictor.run(product, max_predictions=5, timeout=15*60)  # TODO adjust these numbers to something reasonable
                pathways[(host, model)] = predicted_pathways
        return pathways

    def calculate_maximum_yields(self, pathways):
        """"""
        for (host, model), pathway in pathways.iteritems():
            tm = TimeMachine()
            tm(do=partial(model.add_reactions, pathway), undo=partial(model.remove_reactions, pathway))
            maximum_theoretical_yield()


    def __translate_product_to_universal_reactions_model_metabolite(self, product):
        if isinstance(product, Metabolite):
            return product
        elif isinstance(product, StringType):
            search_result = products.search(product)
            print "Found %d compounds that match query '%s'" % (len(search_result), product)
            print repr(search_result)
            print "Choosing best match (%s) ... please interrupt if this is not the desired compound." % search_result.name[0]
            return METANETX['universal_model'].metabolites.get_by_id(search_result.index[0])

design = Designer()


