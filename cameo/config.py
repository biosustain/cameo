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

non_zero_flux_threshold = 1e-6
ndecimals = 6

import logging
log = logging.getLogger(__name__)

try:
    import bokeh
    use_bokeh = True
except ImportError:
    use_bokeh = False

try:
    import matplotlib
    use_matplotlib = True
except ImportError:
    use_matplotlib = False

try:
    from IPython import parallel
    from IPython.kernel.zmq import serialize
    client = parallel.Client()
    client.block = True
    default_view = client.direct_view()
except Exception:
    try:
        from .parallel import MultiprocessingView
        default_view = MultiprocessingView()
    except ImportError:
        from .parallel import SequentialView
        default_view = SequentialView()
# except:
#     log.debug(
#         'IPython parallel not available ... using std lib multiprocessing')
#     try:
#         from .parallel import MultiprocessingView
#         default_view = MultiprocessingView(4)
#     except:
#         log.debug(
#             'multiprocessing not available ... using wrapper for standard map and apply (SequentialView)')
#         from .parallel import SequentialView
#         default_view = SequentialView()