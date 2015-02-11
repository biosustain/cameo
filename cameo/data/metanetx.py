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

__all__ = ['bigg2mnx', 'mnx2bigg']

import os
import pickle
import gzip

import cameo

with open(os.path.join(cameo._cameo_data_path, 'metanetx.pickle')) as f:
    _METANETX = pickle.load(f)

bigg2mnx = _METANETX['bigg2mnx']
mnx2bigg = _METANETX['mnx2bigg']

with gzip.open(os.path.join(cameo._cameo_data_path, 'metanetx_chem_prop.pklz')) as f:
    chem_prop = pickle.load(f)

