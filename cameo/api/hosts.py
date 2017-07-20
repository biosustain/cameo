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

"""This module provides reference models for a number of supported host organisms (currently Escherichia coli and
Sacharomyces cerevisiae)."""

from __future__ import absolute_import, print_function

import os
from functools import partial

import six
from lazy_object_proxy import Proxy

import cameo
from cameo import load_model
from cameo import util


__all__ = ['hosts']


MODEL_DIRECTORY = os.path.join(os.path.join(cameo.__path__[0]), 'models/json')


class Host(object):
    def __init__(self, name='', models=None, biomass=None, carbon_sources=None):
        models = models or []
        biomass = biomass or []
        carbon_sources = carbon_sources or []
        self.name = name
        self.models = util.IntelliContainer()
        for id, biomass, carbon_source in zip(models, biomass, carbon_sources):
            def lazy_model_init(path):
                model = load_model(path)
                setattr(model, "biomass", biomass)
                setattr(model, "carbon_source", carbon_source)
                return model
            model = Proxy(partial(lazy_model_init, os.path.join(MODEL_DIRECTORY, id + '.json')))
            self.models[id] = model

    def __str__(self):
        return self.name


class Hosts(object):
    def __init__(self, host_spec, aliases=None):
        self._host_spec = host_spec
        self._hosts = list()
        for host_id, information in six.iteritems(self._host_spec):
            host = Host(**information)
            self._hosts.append(host)
            setattr(self, host_id, host)
        if aliases and isinstance(aliases, list):
            for pair in aliases:
                setattr(self, pair[1], getattr(self, pair[0]))

    def __iter__(self):
        return iter(self._hosts)

    def __dir__(self):
        return list(self._host_spec.keys())


HOST_SPECS = {
    # 'iAF1260', 'iJO1366', 'EcoliCore'
    'ecoli': {
        'name': 'Escherichia coli',
        'models': ('iJO1366',),
        'biomass': ('BIOMASS_Ec_iJO1366_core_53p95M',),
        'carbon_sources': ('EX_glc__D_e',)
    },
    # 'iND750',
    'scerevisiae': {
        'name': 'Saccharomyces cerevisiae',
        'models': ('iMM904',),
        'biomass': ('BIOMASS_SC5_notrace',),
        'carbon_sources': ('EX_glc__D_e',)

    }
}

hosts = Hosts(HOST_SPECS, aliases=[('scerevisiae', 'yeast')])
