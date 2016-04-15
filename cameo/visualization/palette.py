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

from palettable.colorbrewer.colorbrewer import _load_maps_by_type as load_colorbrewer_palette  # noqa

from palettable.cubehelix.cubehelix import _get_all_maps as load_cubehelix_palettes
from palettable.tableau.tableau import _get_all_maps as load_tableau_palettes
from palettable.wesanderson.wesanderson import _get_all_maps as load_wesanderson_palettes

from palettable.palette import Palette


class PaletteMapper(object):
    def __init__(self):
        self._palettes = {}

    def add_palettalbe_palette(self, name, palette):
        split = name.split("_")
        name = split[0]
        reverse = False
        if len(split) == 2:
            size = int(split[1])
        else:
            if split[-1] == "r":
                reverse = True
                size = int(split[-2])
            else:
                size = int(split[-1])

        if name not in self._palettes:
            self._palettes[name] = {}
        if size not in self._palettes[name]:
            self._palettes[name][size] = {}
        self._palettes[name][size][reverse] = palette

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

for key, obj in six.iteritems(load_colorbrewer_palette("diverging")):
    mapper.add_palettalbe_palette(key, obj)

for key, obj in six.iteritems(load_colorbrewer_palette("sequential")):
    mapper.add_palettalbe_palette(key, obj)

for key, obj in six.iteritems(load_colorbrewer_palette("qualitative")):
    mapper.add_palettalbe_palette(key, obj)

for key, obj in six.iteritems(load_cubehelix_palettes()):
    mapper.add_palettalbe_palette(key, obj)

for key, obj in six.iteritems(load_tableau_palettes()):
    mapper.add_palettalbe_palette(key, obj)

for key, obj in six.iteritems(load_wesanderson_palettes()):
    mapper.add_palettalbe_palette(key, obj)
