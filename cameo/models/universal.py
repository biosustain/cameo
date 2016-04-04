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

from __future__ import absolute_import

import sys
import os
import glob
import cameo
from cameo.io import load_model
from cameo import util

from functools import partial
from lazy_object_proxy import Proxy

__all__ = ['universal']


class ModelDB(object):
    pass


universal = ModelDB()

for file_path in glob.glob(os.path.join(os.path.dirname(cameo.__file__), 'models', 'universal_models', '*.json')):
    model_id = os.path.splitext(os.path.basename(file_path))[0]
    setattr(universal, util.str_to_valid_variable_name(model_id), Proxy(partial(load_model, file_path)))

sys.modules[__name__] = universal
