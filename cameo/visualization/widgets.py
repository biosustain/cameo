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

from IPython.html import widgets
from IPython.utils.traitlets import Unicode


# documentation on widgets:
# http://nbviewer.ipython.org/github/ipython/ipython/blob/master/examples/Interactive%20Widgets/Custom%20Widget%20-%20Hello%20World.ipynb
class OptimizationResultWidget(widgets.DOMWidget):
    _view_name = Unicode("OptimizationResult", sync=True)

    # TODO: create html/js table with scroll. Maybe with http://www.datatables.net
    # TODO: create js for object (select column/value).


