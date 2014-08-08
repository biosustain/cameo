# Copyright 2013 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import json
import logging
import subprocess
import tempfile
from functools import partial
from escher import Builder

import os
from IPython.display import HTML, SVG
from cameo.flux_analysis.analysis import _ids_to_reactions
from cameo.util import TimeMachine


log = logging.getLogger(__name__)

import collections
import functools

from IPython.display import HTML, SVG


class memoized(object):
    '''Decorator. Caches a function's return value each time it is called.
    If called later with the same arguments, the cached value is returned
    (not reevaluated).
    '''

    def __init__(self, func):
        self.func = func
        self.cache = {}

    def __call__(self, *args):
        if not isinstance(args, collections.Hashable):
            # uncachen blow up.
            return self.func(*args)
        if args in self.cache:
            return self.cache[args]
        else:
            value = self.func(*args)
            self.cache[args] = value
            return value

    def __repr__(self):
        '''Return the function's docstring.'''
        return self.func.__doc__

    def __get__(self, obj, objtype):
        '''Support instance methods.'''
        return functools.partial(self.__call__, obj)


def pathviz_maps():
    """Return a list of maps available in pathviz.m"""
    output = subprocess.check_output(['pathviz.m', '-l'])
    return output.split('\n')


def pathviz_svg(map_id='EcoliCore_coreMap', **kwargs):
    config = {"map": map_id, "ImageSize": 800., "Boundary": False}
    for key, value in kwargs.iteritems():
        config[key] = value
    fd, tmp_pathviz_input = tempfile.mkstemp(prefix='pathviz_svg_', suffix='.json')
    with open(tmp_pathviz_input, 'w') as fhandle:
        json.dump(config, fhandle)
    fd, tmp_pathviz_output = tempfile.mkstemp(prefix='pathviz_svg_', suffix='.svg')
    output = subprocess.check_output(['pathviz.m', '-o', tmp_pathviz_output, '-i', tmp_pathviz_input])
    with open(tmp_pathviz_output, 'r') as fhandle:
        svg_map = fhandle.read()
    return SVG(svg_map)


def pathviz_cdf(map_id='EcoliCore_coreMap', **kwargs):
    """Draw pathway using pathviz.m

    Parameters
    ----------
    """
    html_template = """<!DOCTYPE html>
<html>
<head>
  <title>Map</title>
</head>
<body>

<script type="text/javascript" src="http://www.wolfram.com/cdf-player/plugin/v2.1/cdfplugin.js"></script>
<script type="text/javascript">
var cdf = new cdfplugin();
cdf.embed("%s", 942, 678);
</script>

</body>
</html>"""
    html_template2 = """
<iframe width=800 height=700 src="files/%s" id="CDF"></iframe>
"""
    config = {"map": map_id, "ImageSize": 800., "Boundary": False}
    for key, value in kwargs.iteritems():
        config[key] = value
    fd, tmp_pathviz_input = tempfile.mkstemp(prefix='pathviz_svg_', suffix='.json')
    with open(tmp_pathviz_input, 'w') as fhandle:
        json.dump(config, fhandle)
    fd, tmp_pathviz_output = tempfile.mkstemp(prefix='pathviz_cdf_', suffix='.cdf', dir='.')
    output = subprocess.check_output(['pathviz.m', '-o', tmp_pathviz_output, '-i', tmp_pathviz_input])
    log.debug(output)
    fd, tmp_html = tempfile.mkstemp(prefix='pathviz_cdf_embedded_', suffix='.html', dir='.')
    with open(tmp_html, 'w') as fhandle:
        fhandle.write(html_template % os.path.basename(tmp_pathviz_output))
    return HTML(html_template2 % os.path.basename(tmp_html))


def draw_knockout_result(model, map_name, simulation_method, knockouts, *args, **kwargs):
    tm = TimeMachine()

    try:
        for reaction in _ids_to_reactions(model, knockouts):
            tm(do=partial(setattr, reaction, 'lower_bound', 0),
               undo=partial(setattr, reaction, 'lower_bound', reaction.lower_bound))
            tm(do=partial(setattr, reaction, 'upper_bound', 0),
               undo=partial(setattr, reaction, 'upper_bound', reaction.upper_bound))

        solution = simulation_method(model, *args, **kwargs).x_dict
        tm.reset()

        return Builder(map_name, reaction_data=solution)

    except Exception as e:
        tm.reset()
        raise e
