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

non_zero_flux_threshold = 1e-6
ndecimals = 6

import logging
log = logging.getLogger(__name__)

# Determine available solver interfaces
solvers = {}

try:
    from optlang import glpk_interface

    solvers['glpk'] = glpk_interface
except ImportError:
    pass
try:
    from optlang import cplex_interface

    solvers['cplex'] = cplex_interface
except ImportError:
    pass

# Determine if bokeh is available
# TODO: This should also check if a bokeh server is actually running.
try:
    import bokeh
    use_bokeh = True
    bokeh_url = 'default'
except ImportError:
    use_bokeh = False
    bokeh_url = None

try:
    import matplotlib
    use_matplotlib = True
except ImportError:
    use_matplotlib = False

#Determine a default parallelization view


def in_ipnb():
    return False

try:
    from IPython import parallel, get_ipython
    from IPython.kernel.zmq import serialize

    def in_ipnb():
        try:
            return get_ipython().config.get('IPKernelApp').get('parent_appname') == 'ipython-notebook' or \
                get_ipython().config.get('KernelApp').get('parent_appname') == 'ipython-notebook'
        except Exception:
            return False

    client = parallel.Client()
    client.block = True
    default_view = client.direct_view()
except Exception as e:
    from .parallel import SequentialView
    default_view = SequentialView()
    # try:
    #     from .parallel import MultiprocessingView
    #     default_view = MultiprocessingView()
    # except ImportError:
    #     from .parallel import SequentialView
    #     default_view = SequentialView()


