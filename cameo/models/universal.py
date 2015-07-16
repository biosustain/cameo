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

__all__ = ['universal']

import six


import os
import glob
import sys
import six.moves.cPickle as pickle
import cameo
from cameo import util


class ModelFacadeUniversal(util.ModelFacade):

    def _load_model(self):
        if six.PY2:
            return pickle.load(open(self._id))
        else:
            return pickle.load(open(self._id, 'rb'), encoding='bytes')

class ModelDB(object): pass

universal = ModelDB()

for file_path in glob.glob(os.path.join(os.path.dirname(cameo.__file__), 'data', 'universal_models', '*.pickle')):
    model_id = os.path.splitext(os.path.basename(file_path))[0]
    setattr(universal, util.str_to_valid_variable_name(model_id), ModelFacadeUniversal(file_path))
