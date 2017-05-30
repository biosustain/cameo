# -*- coding: utf-8 -*-
# Copyright 2015 Novo Nordisk Foundation Center for Biosustainability, DTU.
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

import os
import pickle
import re

import pytest

from cameo import api, load_model
from cameo import models, config
from cameo.api.hosts import Host
from cameo.api.products import Compound

MODELS = os.path.dirname(models.__file__)

UNIVERSALMODEL = load_model(os.path.join(MODELS, 'json/iJO1366.json'))
UNIVERSALMODEL.remove_reactions(UNIVERSALMODEL.exchanges)


def test_api():
    mock_host = Host('core',
                     models=['e_coli_core'],
                     biomass=['BIOMASS_Ecoli_core_w_GAM'],
                     carbon_sources=['EX_glc__D_e'])

    api.design.debug = True
    pathways = api.design.predict_pathways(product=UNIVERSALMODEL.metabolites.ser__L_c, hosts=[mock_host],
                                           database=UNIVERSALMODEL, aerobic=True)
    optimization_reports = api.design.optimize_strains(pathways, config.default_view, aerobic=True)
    pickle.loads(pickle.dumps(optimization_reports))
    assert len(optimization_reports) > 0


def test_compound_repr():
    pytest.mark.skipif(not re.match('Open Babel.*', os.popen('obabel').read()),
                       reason='Skipping because OpenBabel is not installed.')
    compound = Compound('InChI=1S/H2O/h1H2')
    assert re.match(r"^<\?xml version=\"1\.0\"\?>.*</svg>$", compound._repr_svg_().replace('\n', ''))
    assert compound._repr_html_() == compound._repr_svg_()


def test_products():
    assert api.products.search('3-hydroxy propionate').index[0] == 'MNXM872'
    assert len(api.products.search('old spice')) == 0


def test_hosts():
    assert api.hosts.ecoli.models.iJO1366.id == 'iJO1366'
    assert api.hosts.scerevisiae.models.iMM904.id == 'iMM904'
