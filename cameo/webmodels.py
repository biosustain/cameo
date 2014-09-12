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

"""
WebModels API
-------------

An high level API for retrieving models from the
http://darwin.di.uminho.pt/models database
"""

import json
import requests
from pandas import DataFrame
import tempfile


class NotFoundException(Exception):
    def __init__(self, type, index, *args, **kwargs):
        message = "Could not retrieve %s for entry with index %i" % (type, index)
        Exception.__init__(self, message, *args, **kwargs)


def index_models(host="http://darwin.di.uminho.pt/models"):
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
    response = json.loads(requests.get(uri).text)
    return DataFrame(response, columns=["name", "doi", "author", "year", "formats", "organism", "taxonomy"])


def get_sbml_file(index, host="http://darwin.di.uminho.pt/models"):
    temp = tempfile.NamedTemporaryFile()
    uri = host + "/models/%i.sbml" % index
    response = requests.get(uri)
    if response.status_code == 200:

        temp.write(response.text)
        temp.flush()
        return temp
    raise NotFoundException("sbml", index)

if __name__ == "__main__":
    print index_models()
    from cameo import load_model
    model = load_model(get_sbml_file(2))
    print model.objective
