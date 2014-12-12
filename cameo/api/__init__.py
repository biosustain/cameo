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

# import logging
# logging.basicConfig(level=logging.INFO)
# logger = logging.getLogger(__name__)
#
# logger.info('Initializing cameo advanced programming interface. Be patient this might take a while ...')

import os
import gzip
import cPickle as pickle
import cameo

# with gzip.open(os.path.join(cameo._cameo_data_path, 'metanetx.pgz')) as f:
#     _METANETX = pickle.load(f)

with open(os.path.join(cameo._cameo_data_path, 'metanetx.pickle')) as f:
    _METANETX = pickle.load(f)

from cameo.api.hosts import hosts
from cameo.api.designer import design
from cameo.api.products import products

