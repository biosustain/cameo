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

import six


def _():

    import os
    import glob
    import sys
    import six.moves.cPickle as pickle
    from cameo import Model

    CURRENT_MODULE = sys.modules[__name__]
    CURRENT_PATH = os.path.dirname(os.path.realpath(__file__))

    class ModelFacade(Model):

        def __init__(self, id):
            self.id = id
            self._model = None

        def __getattr__(self, value):
            if self._model is None:
                self._load_lazily()
                return getattr(self._model, value)
            else:
                return getattr(self._model, value)

        def __dir__(self):
            if self._model is None:
                self._load_lazily()
            return dir(self._model)

        def _load_lazily(self):
            with open(os.path.join(CURRENT_PATH, self.id + '.pickle')) as f:
                if six.PY2:
                    self._model = pickle.load(f)
                else:
                    self._model = pickle.load(f, encoding='bytes')

    for file_path in glob.glob(os.path.join(CURRENT_PATH, '*.pickle')):
        model_id = os.path.splitext(os.path.basename(file_path))[0]
        setattr(CURRENT_MODULE, model_id, ModelFacade(model_id))

_()