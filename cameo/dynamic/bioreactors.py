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
import numpy as np


def batch_reactor(t, volume, x, s, inflow_rate, outflow_rate, delta_x, x_feed, delta_s, s_feed):
    assert inflow_rate == 0 and outflow_rate == 0
    return 0, 0, np.zeros(len(delta_x)), np.zeros(len(x_feed)), np.zeros(len(delta_s)), np.zeros(len(s_feed))
