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
import six


class MultipleChromosomeGenome(object):
    def __init__(self, keys=[], *args, **kwargs):
        super(MultipleChromosomeGenome, self).__init__(*args, **kwargs)
        self.chromosomes = {k: OrderedSet() for k in keys}
        self.keys = keys

    def __getitem__(self, key):
        return self.chromosomes[key]

    def __delitem__(self, key):
        del self.chromosomes[key]

    def __setitem__(self, key, value):
        self.chromosomes[key] = OrderedSet(value)

    def copy(self):
        new_genome = MultipleChromosomeGenome(self.keys)
        for key in self.keys:
            new_genome[key] = self[key]
        return new_genome

    def __repr__(self):
        return "| ".join(["%s: %s" % (k, list(v)) for k, v in six.iteritems(self.chromosomes)])
