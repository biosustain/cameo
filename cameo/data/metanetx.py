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

from __future__ import absolute_import, print_function

__all__ = ['bigg2mnx', 'mnx2bigg', 'all2mnx', 'mnx2all']

import os
import json
import gzip

import pandas

import cameo


with gzip.open(os.path.join(cameo._cameo_data_path, 'metanetx.json.gz'),
               'rt') as f:
    _METANETX = json.load(f)

bigg2mnx = _METANETX['bigg2mnx']
mnx2bigg = _METANETX['mnx2bigg']
all2mnx = _METANETX['all2mnx']
mnx2all = {v: k for k, v in all2mnx.items()}

with gzip.open(os.path.join(cameo._cameo_data_path,
                            'metanetx_chem_prop.json.gz'), 'rt') as f:
    chem_prop = pandas.read_json(f)
