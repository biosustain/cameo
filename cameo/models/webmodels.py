# Copyright 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
WebModels API
-------------

A high level API for retrieving models from the
http://darwin.di.uminho.pt/models and http://bigg.ucsd.edu databases
"""

from __future__ import absolute_import, print_function


import io
import tempfile
import logging

import requests

import lazy_object_proxy
import optlang

from functools import partial

from pandas import DataFrame

from cobra.io import load_json_model, read_sbml_model

from cameo.util import str_to_valid_variable_name
from cameo.core.solver_based_model import to_solver_based_model


__all__ = ['index_models_minho', 'index_models_bigg', 'bigg', 'minho']

logger = logging.getLogger(__name__)


class NotFoundException(Exception):
    def __init__(self, type, index, *args, **kwargs):
        message = "Could not retrieve %s for entry with index %i" % (type, index)
        Exception.__init__(self, message, *args, **kwargs)


def load_webmodel(query, solver_interface):
    logger.debug('Querying webmodels ... trying http://bigg.ucsd.edu first')
    try:
        model = get_model_from_bigg(query, solver_interface=solver_interface)
    except Exception:
        logger.debug('Querying webmodels ... trying minho next')
        try:
            df = index_models_minho()
        except requests.ConnectionError as e:
            logger.error("You need to be connected to the internet to load an online model.")
            raise e
        except Exception as e:
            logger.error("Something went wrong while looking up available webmodels.")
            raise e
        try:
            index = df.query('name == "%s"' % query).id.values[0]
            model = get_model_from_uminho(index)
            # handle = get_sbml_file(index)
            # path = handle.name
        except IndexError:
            raise ValueError("%s is neither a file nor a model ID." % query)
    return model


def index_models_minho(host="http://darwin.di.uminho.pt/models"):
    """
    Retrieves a summary of all models in the database.

    Parameters
    ----------
    host: the service host (optional, default: http://darwin.di.uminho.pt/models)

    Returns
    -------
    pandas.DataFrame
        summary of the models in the database
    """
    uri = host + "/models.json"
    try:
        response = requests.get(uri)
    except requests.ConnectionError as e:
        logger.error("Cannot reach %s. Are you sure that you are connected to the internet?" % host)
        raise e
    if response.ok:
        try:
            json = response.json()
        except Exception as e:
            logger.error('No json could be decoded from server response coming from {}.'.format(host))
            raise e
        else:
            index = DataFrame(json, columns=["id", "name", "doi", "author",
                                             "year", "formats", "organism",
                                             "taxonomy", "optflux_validated"])
            index.columns = ["id", "name", "doi", "author", "year", "formats", "organism", "taxonomy", "validated"]
            return index
    else:
        raise Exception("Could not index available models. %s returned status code %d" % (host, response.status_code))


def get_model_from_uminho(index, host="http://darwin.di.uminho.pt/models"):
    sbml_file = get_sbml_file(index, host)
    sbml_file.close()
    return to_solver_based_model(read_sbml_model(sbml_file.name))


def get_sbml_file(index, host="http://darwin.di.uminho.pt/models"):
    temp = tempfile.NamedTemporaryFile(delete=False)
    uri = host + "/models/%i.sbml" % index
    try:
        response = requests.get(uri)
    except requests.ConnectionError as e:
        logger.error("Cannot reach {}. Are you sure that you are connected to the internet?".format(host))
        raise e
    if response.ok:

        temp.write(response.text.encode('utf-8'))
        temp.flush()
        return temp
    else:
        raise NotFoundException("sbml", index)


def index_models_bigg():
    try:
        response = requests.get('http://bigg.ucsd.edu/api/v2/models')
    except requests.ConnectionError as e:
        logger.error("Cannot reach http://bigg.ucsd.edu. Are you sure that you are connected to the internet?")
        raise e
    if response.ok:
        try:
            json = response.json()
        except Exception as e:
            logger.error('No json could be decoded from server response coming from http://bigg.ucsd.edu.')
            raise e
        else:
            return DataFrame.from_dict(json['results'])
    else:
        raise Exception(
            "Could not index available models. bigg.ucsd.edu returned status code {}".format(response.status_code))


def get_model_from_bigg(id, solver_interface=optlang):
    try:
        response = requests.get('http://bigg.ucsd.edu/api/v2/models/{}/download'.format(id))
    except requests.ConnectionError as e:
        logger.error("Cannot reach http://bigg.ucsd.edu. Are you sure that you are connected to the internet?")
        raise e
    if response.ok:
        with io.StringIO(response.text) as f:
            return to_solver_based_model(load_json_model(f), solver_interface=solver_interface)
    else:
        raise Exception(
            "Could not download model {}. bigg.ucsd.edu returned status code {}".format(id, response.status_code))


class ModelDB(object):
    pass


bigg = ModelDB()
try:
    model_ids = index_models_bigg().bigg_id
except requests.ConnectionError:
    bigg.no_models_available = "Cameo couldn't reach http://bigg.ucsd.edu at initialization time." \
                               "Are you connected to the internet?"
except Exception as e:
    bigg.no_models_available = "Cameo could reach http://bigg.ucsd.edu at initialization time" \
                               "but something went wrong while decoding the server response."
    logger.debug(e)
else:
    for id in model_ids:
        setattr(bigg, str_to_valid_variable_name(id), lazy_object_proxy.Proxy(partial(get_model_from_bigg, id)))

minho = ModelDB()
try:
    minho_models = index_models_minho()
except requests.ConnectionError as e:
    minho.no_models_available = "Cameo couldn't reach http://darwin.di.uminho.pt/models at initialization time." \
                                "Are you connected to the internet?"
    logger.debug(e)
except Exception as e:
    minho.no_models_available = "Cameo could reach http://darwin.di.uminho.pt/models at initialization time" \
                                "but something went wrong while decoding the server response."
    logger.debug(e)
else:
    model_indices = minho_models.id
    model_ids = minho_models.name
    for index, id in zip(model_indices, model_ids):
        setattr(minho, str_to_valid_variable_name(id), lazy_object_proxy.Proxy(partial(get_model_from_uminho, index)))

    validated_models = minho_models[minho_models.validated]
    minho.validated = ModelDB()
    model_indices = validated_models.id
    model_ids = validated_models.name
    for index, id in zip(model_indices, model_ids):
        setattr(minho.validated, str_to_valid_variable_name(id),
                lazy_object_proxy.Proxy(partial(get_model_from_uminho, index)))


if __name__ == "__main__":
    print(index_models_minho())
    from cameo import load_model

    model = load_model(get_sbml_file(2))
    print(model.objective)
