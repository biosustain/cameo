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

from cobra.core.DictList import DictList
from cameo import load_model

import sys
mod = sys.modules[__name__]

class ModelFacade(object):

    def __init__(self, id):
        self.id = id
        self.model = None

    def __getattr__(self, value):
        if getattr(self,) is None:
            self.model load_model('../test/data/' + self.id + '.xml')
        else:
            return getattr(self.model)


class Host(object):

    def __init__(self, name='', models=[]):
        self.name = name
        self.models = DictList()
        for id in models:
            blub = ModelFaced(id)
            self.models.append(blub)

hosts = {'ecoli': {'name': 'Escherichia coli', 'models': ('iAF1260', 'iJO1366')},
         'scerevisiae': {'name': 'Saccharomyces cerevisiae', 'models': ('iND750', 'iMM904')}
}

for host_id, information in hosts.iteritems():
    setattr(mod, host_id, Host(**information))