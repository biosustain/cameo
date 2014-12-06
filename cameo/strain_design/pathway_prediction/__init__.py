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

__all__ = ['PathwayPredictor']

from functools import partial
from progressbar import ProgressBar
from cameo.exceptions import SolveError
from cameo.util import TimeMachine
from sympy import Add

import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class PathwayPredictor(object):

    def __init__(self, model, universal_model):
        self.model = model.copy()
        for exchange in self.model.exchanges:
            if len(exchange.get_reactants()) > 0 and exchange.lower_bound <= 0:
                exchange.upper_bound = 999999.
        self._y_vars = list()
        switches = list()
        reactions = list()
        for reaction in universal_model.reactions:
            if reaction in self.model.reactions:
                continue
            reactions.append(reaction.copy())
        self.model.add_reactions(reactions)

        pbar = ProgressBar()
        for reaction in pbar(reactions):
            y = self.model.solver.interface.Variable('y_'+reaction.id, lb=0, ub=1, type='binary')
            self._y_vars.append(y)
            switch_lb = self.model.solver.interface.Constraint(y*reaction.lower_bound - reaction.variable, name='switch_lb_'+reaction.id, ub=0)
            switches.append(switch_lb)
            switch_ub = self.model.solver.interface.Constraint(y*reaction.upper_bound - reaction.variable, name='switch_ub_'+reaction.id, lb=0)
            switches.append(switch_ub)
        self.model.solver.add(switches)
        self.model.objective = self.model.solver.interface.Objective(Add(*self._y_vars), direction='min')

    def run(self, target=None, max_predictions=float("inf"), min_production=.1, timeout=None):
        pathways = list()
        with TimeMachine() as tm:
            tm(do=partial(setattr, self.model.solver.configuration, 'timeout', timeout),
               undo=partial(setattr, self.model.solver.configuration, 'timeout', self.model.solver.configuration.timeout))
            demand_reaction = self.model.add_demand(target)
            tm(do=str, undo=partial(self.model.remove_reactions, [demand_reaction]))
            demand_reaction.lower_bound = min_production
            self.model.solver.configuration.verbosity = 3
            counter = 0
            while counter <= max_predictions:
                counter += 1
                try:
                    solution = self.model.solve()
                except SolveError:
                    break
                vars_to_cut = list()
                for y_var in self._y_vars:
                    if y_var.primal == 1.0:
                        vars_to_cut.append(y_var)
                pathways.append([self.model.reactions.get_by_id(y_var.name[2:]) for y_var in vars_to_cut])
                integer_cut = self.model.solver.interface.Constraint(Add(*vars_to_cut), name="integer_cut_" + str(counter), ub=len(vars_to_cut)-1)
                # print "Adding as constraint", integer_cut  # TODO: add logger
                tm(do=partial(self.model.solver._add_constraint, integer_cut),
                   undo=partial(self.model.solver._remove_constraint, integer_cut))
            self.model.solver.configuration.verbosity = 0
        return pathways

if __name__ == '__main__':
    from cameo import api
    universal_model = api.hosts.ecoli.models.iJO1366
    model = api.hosts.ecoli.models.EcoliCore

    pathway_predictor = PathwayPredictor(model, universal_model)
    pathways = pathway_predictor.run(target=universal_model.metabolites.ser_dsh_L_c, max_predictions=3, timeout=None)
    print pathways
