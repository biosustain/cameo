# Copyright 2015 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import absolute_import, print_function
from uuid import uuid4
from IPython.core.display import HTML, Javascript
from IPython.core.display import display

from cameo import util
import logging
import os


ASSETS = os.path.join(os.path.dirname(__file__), "assets")

SEARCHING_IMAGE_FILE = os.path.join(ASSETS, "searching.gif")
with open(SEARCHING_IMAGE_FILE, "rb") as f:
    SEARCHING_IMAGE = str(f.read()).encode('base64').replace('\n', '')


LOADING_IMAGE_FILE = os.path.join(ASSETS, "searching.gif")
with open(SEARCHING_IMAGE_FILE, "rb") as f:
    LOADING_IMAGE = str(f.read()).encode('base64').replace('\n', '')


logger = logging.getLogger(__name__)


def notice(message):
    if util.in_ipnb():
        display(HTML("<span>%s</span>" % message))
    else:
        print(message)


def bold(message):
    if util.in_ipnb():
        display(HTML("<strong>%s</strong>" % message))
    else:
        print("\033[1m" + message + "\033[0m")


def searching():
    if util.in_ipnb():
        identifier = str(uuid4())
        display(HTML("""
        <img class="loading" id="%s" style="margin:auto; text-align:center;" src="data:image/gif;base64,%s"/>
        """ % (identifier, SEARCHING_IMAGE)))
        return identifier
    else:
        logger.debug("loading only works on Jupyter notebooks")


def loading():
    if util.in_ipnb():
        identifier = str(uuid4())
        display(HTML("""
        <img class="loading" id="%s" style="margin:auto; text-align:center;" src="data:image/gif;base64,%s"/>
        """ % (identifier, LOADING_IMAGE)))
        return identifier
    else:
        logger.debug("loading only works on Jupyter notebooks")


def stop_loader(identifier):
    if util.in_ipnb():
        display(Javascript("""
        jQuery("#%s").remove();
        """ % identifier))
    else:
        logger.debug("loading only works on Jupyter notebooks")