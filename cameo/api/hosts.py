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


import os
import cameo
from cameo.util import IntelliContainer
from cameo import load_model

MODEL_DIRECTORY = os.path.join(os.path.join(os.path.split(cameo.__path__[0])[0]), 'tests/data')


class Host(object):

    def __init__(self, name='', models=[]):
        self.name = name
        self.models = IntelliContainer()
        for id in models:
            self.models[id] = ModelFacade(id)


class ModelFacade(object):

    def __init__(self, id):
        self.id = id
        self._model = None

    def __getattr__(self, value):
        if self._model is None:
            super(ModelFacade, self).__setattr__('_model', load_model(os.path.join(MODEL_DIRECTORY, self.id + '.xml')))
            return getattr(self._model, value)
        else:
            return getattr(self._model, value)

    def __dir__(self):
        if self._model is None:
            self._model = load_model(os.path.join(MODEL_DIRECTORY, self.id + '.xml'))
        return dir(self._model)


class Hosts(object):

    def __init__(self, host_spec):
        self.host_spec = host_spec
        for host_id, information in self.host_spec.iteritems():
            setattr(self, host_id, Host(**information))


    # def __getattr__(self, value):
    #     if getattr(self, 'host_spec').has_key(value):
    #         model load_model('../test/data/' + self.id + '.xml')
    #     else:
    #         raise AttributeError("Host %s is not available" % value)

    def __dir__(self):
        return self.host_spec.keys()


HOST_SPECS = {'ecoli': {'name': 'Escherichia coli', 'models': ('iAF1260', 'iJO1366', 'EcoliCore')},
         'scerevisiae': {'name': 'Saccharomyces cerevisiae', 'models': ('iND750', 'iMM904')}
}

hosts = Hosts(HOST_SPECS)
