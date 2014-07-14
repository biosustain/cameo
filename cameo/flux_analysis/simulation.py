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
from optlang import Variable, Constraint, Objective


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
    solution = fba(model, objective=objective)
    tm = TimeMachine()
    objective_terms = list()
    additional_constraints = list()
    for reaction in model.reactions:
        if reaction.reversibility:
            aux_variable = Variable(reaction.variable.id + '_aux', lb=0, ub=-1*reaction.lower_bound)
            tm(do=partial(setattr, reaction, 'lower_bound', 0), undo=partial(setattr, reaction, 'lower_bound', reaction.lower_bound))
            objective_terms.append(reaction.variable)
            objective_terms.append(aux_variable)
            additional_constraints.append(Constraint)
        else:
            objective_terms.append(reaction.variable)
    abs_objective = Objective(sum(objective_terms), name='Minimize total flux', direction='min')
    tm(do=partial(setattr, model, 'objective', abs_objective))

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