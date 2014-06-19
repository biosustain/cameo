# Copyright 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from functools import partial
from cameo.util import TimeMachine


def fba(model, objective=None):
    """Perform flux balance analysis."""
    tm = TimeMachine()
    if objective is not None:
        tm(do=partial(setattr, model, 'objective', objective),
           undo=partial(setattr, model, 'objective', model.objective))
    solution = model.solve()
    tm.reset()
    return solution


def pfba(model, objective=None):
    pass

def moma(model, objective=None):
    pass

def lmoma(model, objective=None):
    pass


if __name__ == '__main__':
    from cameo import load_model

    model = load_model('../../tests/data/EcoliCore.xml')
    print fba(model)
    print fba(model, 'TPI')
    print model.objective