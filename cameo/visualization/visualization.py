# Copyright 2013 Novo Nordisk Foundation Center for Biosustainability, DTU.
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

import os
import six
import json
import cameo
import logging
import tempfile
import subprocess

import networkx as nx

from functools import partial
from io import BytesIO
from escher import Builder
from cameo.util import TimeMachine

try:
    from IPython.display import HTML, SVG
except ImportError:
    pass

__all__ = ['graph_to_svg', 'draw_knockout_result', 'inchi_to_svg']

logger = logging.getLogger(__name__)


def pathviz_maps():
    """Return a list of maps available in pathviz.m"""
    output = subprocess.check_output(['pathviz.m', '-l'])
    return output.split('\n')


def pathviz_svg(map_id='EcoliCore_coreMap', **kwargs):
    config = {"map": map_id, "ImageSize": 800., "Boundary": False}
    for key, value in six.iteritems(kwargs):
        config[key] = value
    fd, tmp_pathviz_input = tempfile.mkstemp(prefix='pathviz_svg_', suffix='.json')
    with open(tmp_pathviz_input, 'w') as fhandle:
        json.dump(config, fhandle)
    fd, tmp_pathviz_output = tempfile.mkstemp(prefix='pathviz_svg_', suffix='.svg')
    output = subprocess.check_output(['pathviz.m', '-o', tmp_pathviz_output, '-i', tmp_pathviz_input])
    logger.debug(output)
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
    for key, value in six.iteritems(kwargs):
        config[key] = value
    fd, tmp_pathviz_input = tempfile.mkstemp(prefix='pathviz_svg_', suffix='.json')
    with open(tmp_pathviz_input, 'w') as fhandle:
        json.dump(config, fhandle)
    fd, tmp_pathviz_output = tempfile.mkstemp(prefix='pathviz_cdf_', suffix='.cdf', dir='.')
    output = subprocess.check_output(['pathviz.m', '-o', tmp_pathviz_output, '-i', tmp_pathviz_input])
    logger.debug(output)
    fd, tmp_html = tempfile.mkstemp(prefix='pathviz_cdf_embedded_', suffix='.html', dir='.')
    with open(tmp_html, 'w') as fhandle:
        fhandle.write(html_template % os.path.basename(tmp_pathviz_output))
    return HTML(html_template2 % os.path.basename(tmp_html))


def draw_knockout_result(model, map_name, simulation_method, knockouts, *args, **kwargs):
    tm = TimeMachine()

    try:
        for reaction in model._ids_to_reactions(knockouts):
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


def inchi_to_svg(inchi, file=None, debug=False, three_d=False):
    """Generate an SVG drawing from an InChI string.

    Parameters
    ----------
    inchi : str
        An InChI string.

    Returns
    -------
    str
        A vector graphics of the compound represented as SVG.

    Examples
    --------
    Draw water
    >>> inchi_to_svg('InChI=1S/H2O/h1H2')
    '<?xml version="1.0"?>\n<svg version="1.1" id="topsvg"\nxmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"\nxmlns:cml="http://www.xml-cml.org/schema" x="0" y="0" width="200px" height="200px" viewBox="0 0 100 100">\n<title>OBDepict</title>\n<rect x="0" y="0" width="100" height="100" fill="white"/>\n<text text-anchor="middle" font-size="6" fill ="black" font-family="sans-serif"\nx="50" y="98" ></text>\n<g transform="translate(0,0)">\n<svg width="100" height="100" x="0" y="0" viewBox="0 0 80 80"\nfont-family="sans-serif" stroke="rgb(0,0,0)" stroke-width="2"  stroke-linecap="round">\n<text x="36" y="48" fill="rgb(255,12,12)"  stroke="rgb(255,12,12)" stroke-width="1" font-size="16" >OH</text>\n<text x="60" y="51.68" fill="rgb(255,12,12)"  stroke="rgb(255,12,12)" stroke-width="1" font-size="13" >2</text>\n</svg>\n</g>\n</svg>\n\n'  # noqa
    """
    in_file = tempfile.NamedTemporaryFile()
    in_file.write(inchi.encode('utf-8'))
    in_file.flush()

    out_file = None
    gen = "--gen3d" if three_d else "--gen2d"
    error_level = 5 if debug else 1
    try:
        if file is not None:
            os.system("obabel -iinchi %s -osvg -O %s %s -xh 40 ---errorlevel %d" %
                      (in_file.name, file.name, gen, error_level))
            return file.name
        else:
            out_file = tempfile.NamedTemporaryFile("w+")
            os.system("obabel -iinchi %s -osvg -O %s %s -xh 40 ---errorlevel %d"
                      % (in_file.name, out_file.name, gen, error_level))
            return out_file.read()
    finally:
        in_file.close()
        if out_file is not None:
            out_file.close()


def inchi_to_ascii(inchi, file=None, debug=False):
    """Generate an ASCII drawing from an InChI string.

    Parameters
    ----------
    inchi : str
        An InChI string.

    Returns
    -------
    str
        A vector graphics of the compound represented as SVG.

    Examples
    --------
    Draw water
    >>> inchi_to_ascii('InChI=1S/H2O/h1H2')
    """

    in_file = tempfile.NamedTemporaryFile()
    in_file.write(inchi.encode('utf-8'))
    in_file.flush()

    out_file = None

    error_level = 5 if debug else 1
    try:
        if file is not None:
            os.system("obabel -iinchi %s -oascii -O %s --gen3d -xh 40 ---errorlevel %d"
                      % (in_file.name, file.name, error_level))
            return file.name
        else:
            out_file = tempfile.NamedTemporaryFile()
            os.system("obabel -iinchi %s -oascii -O %s --gen3d -xh 40 ---errorlevel %d"
                      % (in_file.name, out_file.name, error_level))
            return out_file.read()
    finally:
        in_file.close()
        if out_file is not None:
            out_file.close()


def graph_to_svg(g, layout=nx.spring_layout):
    """return the SVG of a matplotlib figure generated from a graph"""
    import matplotlib.pyplot as plt

    layout = layout(g)
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    # draw reaction nodes
    # rxn_nodes = [node for node in g.nodes() if isinstance(node, cameo.Reaction)]
    # draw metabolites
    met_nodes = [node for node in g.nodes() if isinstance(node, cameo.Metabolite)]
    nx.draw_networkx_edges(g, nodelist=met_nodes, pos=layout, ax=ax, edge_color='gray', arrows=False, node_color='b')
    labels = dict()
    for node in g.nodes():
        if isinstance(node, cameo.Reaction):
            labels[node] = node.name
        else:
            labels[node] = node.name
    nx.draw_networkx_labels(g, labels=labels, pos=layout, font_size=9)
    output = BytesIO()
    fig.savefig(output, format='svg')
    plt.close(fig)
    return output.getvalue()
