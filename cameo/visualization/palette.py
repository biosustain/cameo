# Copyright 2016 Joao Cardoso
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import six

from palettable.colorbrewer.diverging import *
from palettable.colorbrewer.qualitative import *
from palettable.colorbrewer.sequential import *

from palettable.cubehelix.cubehelix import Cubehelix
from palettable.tableau.tableau import TableauMap
from palettable.wesanderson.wesanderson import WesAndersonMap

from palettable.palette import Palette


class PaletteMapper(object):
    def __init__(self):
        self._palettes = {}

    def add_palettalbe_palette(self, key, obj):
        if isinstance(obj, Palette):
            split = key.split("_")
            name = split[0]
            size = int(split[1])
            reverse = len(split) == 3
            if name not in self._palettes:
                self._palettes[name] = {}
            if size not in self._palettes[name]:
                self._palettes[name][size] = {}
            self._palettes[name][size][reverse] = obj

    def map_palette(self, name, n, reverse=False):
        palette = self._palettes[name]
        max_n = max(palette.keys())
        min_n = min(palette.keys())
        if n < min_n:
            n = min_n
        if n > max_n:
            n = max_n

        return palette[n][reverse]


mapper = PaletteMapper()
for key, obj in six.iteritems(dict(globals())):
    mapper.add_palettalbe_palette(key, obj)
