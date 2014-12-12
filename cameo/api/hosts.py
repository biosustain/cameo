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

    def __str__(self):
        return self.name


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
        self._host_spec = host_spec
        self._hosts = list()
        for host_id, information in self._host_spec.iteritems():
            host = Host(**information)
            self._hosts.append(host)
            setattr(self, host_id, host)

    def __iter__(self):
        return iter(self._hosts)

    def __dir__(self):
        return self._host_spec.keys()


HOST_SPECS = {'ecoli': {'name': 'Escherichia coli', 'models': ('EcoliCore', 'iJO1366',)}, #  'iAF1260', 'iJO1366',
            'scerevisiae': {'name': 'Saccharomyces cerevisiae', 'models': ('iMM904', )} # 'iND750',
}

hosts = Hosts(HOST_SPECS)
