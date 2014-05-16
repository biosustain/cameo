# Copyright 2013 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import logging
log = logging.getLogger(__name__)

try:
    from IPython import parallel
    client = parallel.Client()
    client.block = True
    default_view = client.direct_view()
except:
    log.debug(
        'IPython parallel not available ... using std lib multiprocessing')
    try:
        from .parallel import MultiprocessingView
        default_view = MultiprocessingView(4)
    except:
        log.debug(
            'multiprocessing not available ... using wrapper for standard map and apply (SequentialView)')
        from .parallel import SequentialView
        default_view = SequentialView()
