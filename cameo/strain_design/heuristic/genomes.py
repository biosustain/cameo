# Copyright 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.
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

from ordered_set import OrderedSet


class MultipleChromosomeGenome(object):
    def __init__(self, keys=[], *args, **kwargs):
        super(MultipleChromosomeGenome, self).__init__(*args, **kwargs)
        self.chromosomes = {}
        self.keys = keys
        for key in keys:
            self.chromosomes[key] = OrderedSet()

    def __getitem__(self, key):
        return self.chromosomes[key]

    def __delitem__(self, key):
        del self.chromosomes[key]

    def __setitem__(self, key, value):
        self.chromosomes[key] = value
