# -*- coding: utf-8 -*-
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

import os
import re

import nose
from nose.plugins.skip import SkipTest
from nose.tools import assert_equal, assert_raises_regexp

from cameo import api
from cameo.api.products import Compound


def test_compound_repr():
    if not re.match('Open Babel!.*', os.popen('obabel').read()):
        raise SkipTest('Skipping because OpenBabel is not installed.')
    compound = Compound('InChI=1S/H2O/h1H2')
    reference_output = '''<?xml version="1.0"?>
<svg xmlns="http://www.w3.org/2000/svg"
xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:cml="http://www.xml-cml.org/schema" width="80" height="80" x="0" y="0" font-family="sans-serif" stroke="rgb(0,0,0)" stroke-width="1"  stroke-linecap="round">
<rect x="0%" y="0%" width="100%" height="100%" stroke-width="0" fill="rgb(255,255,255)"  />
<text x="36" y="48" fill="rgb(255,12,12)"  stroke="rgb(255,12,12)" stroke-width="1" font-size="16" >OH</text>
<text x="60" y="51.68" fill="rgb(255,12,12)"  stroke="rgb(255,12,12)" stroke-width="1" font-size="13" >2</text>
<text font-size="18" fill ="black" font-family="sans-serif"
x="10" y="20" ></text>
<title> - OBDepict</title>
</svg>
'''
    assert_equal(compound._repr_svg_(), reference_output)
    assert_equal(compound._repr_html_(), compound._repr_svg_())


def test_products():
    assert_equal(api.products.search('3-hydroxy propionate').index[0], 'MNXM872')
    with assert_raises_regexp(Exception, "No compound matches found for query.*"):
        api.products.search('old spice')


def test_hosts():
    assert_equal(api.hosts.ecoli.models.iJO1366.id, 'iJO1366')
    assert_equal(api.hosts.scerevisiae.models.iMM904.id, 'iMM904')


if __name__ == '__main__':
    nose.runmodule()
