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

import os
import gzip
import cPickle as pickle
from functools import partial

import cameo
from cameo.exceptions import SolveError
from cameo import Reaction
from cameo.util import TimeMachine
from sympy import Add

import logging
logger = logging.getLogger(__name__)


class PathwayPredictor(object):
    """Pathway predictions from a universal set of reaction.

    Parameters
    ----------
    model : SolverBasedModel
        The model that represents the host organism.
    universal_model : SolverBasedModel, optional
        The model that represents the universal set of reactions.
        A default model will be used if omitted.
    mapping : dict, optional
        A dictionary that contains a mapping between metabolite
        identifiers in `model` and `universal_model`

    Attributes
    ----------
    model : SolverBasedModel
        The provided model + universal_model + adapter reactions

    Examples
    --------
    Determine production pathways for propane-1,3-diol (MNXM2861 in the metanetx namespace)

    >>> from cameo.api import hosts
    >>> pathway_predictor = PathwayPredictor(hosts.ecoli.iJO1366)
    >>> pathway_predictor.run(product=pathway_predictor.model.metabolites.MNXM2861)
    """

    def __init__(self, model, universal_model=None, mapping=None):
        """"""
        if universal_model is None:
            logger.debug("Loading default universal model.")
            from cameo.api import _METANETX as metanetx
            universal_model = metanetx['universal_model']
        self.model = model.copy()
        for exchange in self.model.exchanges:
            if len(exchange.reactants) > 0 and exchange.lower_bound <= 0:
                exchange.upper_bound = 999999.

        logger.debug("Adding reactions from universal model to host model.")
        reactions = list()
        for reaction in universal_model.reactions:
            if reaction in self.model.reactions:
                continue
            reactions.append(reaction.copy())
        self.model.add_reactions(reactions)

        logger.debug("Adding adapter reactions to connect model with universal model.")
        self._adpater_reactions = list()
        if mapping is None:
            mapping = metanetx['bigg2mnx']
        for metabolite in model.metabolites:  # model is the original host model
            name = metabolite.id[0:-2]
            try:
                mnx_name = mapping[name]
            except KeyError:
                continue
                # print name, 'N/A'
            adapter_reaction = Reaction('adapter_' + metabolite.id + '_' + mnx_name)
            adapter_reaction.lower_bound = -1000
            try:
                adapter_reaction.add_metabolites({self.model.metabolites.get_by_id(metabolite.id): -1, self.model.metabolites.get_by_id(mnx_name): 1})
            except KeyError:
                pass
            else:
                self._adpater_reactions.append(adapter_reaction)
                self.model.add_reaction(adapter_reaction)

        logger.debug("Adding switches.")
        self._y_vars = list()
        switches = list()
        for reaction in reactions:
            if reaction.id.startswith('DM_'):
                continue  # demand reactions don't need integer switches
            y = self.model.solver.interface.Variable('y_'+reaction.id, lb=0, ub=1, type='binary')
            self._y_vars.append(y)
            switch_lb = self.model.solver.interface.Constraint(y*reaction.lower_bound - reaction.variable, name='switch_lb_'+reaction.id, ub=0)
            switches.append(switch_lb)
            switch_ub = self.model.solver.interface.Constraint(y*reaction.upper_bound - reaction.variable, name='switch_ub_'+reaction.id, lb=0)
            switches.append(switch_ub)
        self.model.solver.add(switches)
        logger.debug("Setting minimization of switch variables as objective.")
        self.model.objective = self.model.solver.interface.Objective(Add(*self._y_vars), direction='min')
        self.model.solver.configuration.verbosity = 3

    def run(self, product=None, max_predictions=float("inf"), min_production=.1, timeout=None):
        """Run pathway prediction for a desired product.

        Parameters
        ----------
        product : Metabolite
        max_predictions : int, optional
            The maximum number of predictions to compute.
        min_production : float
            The minimum acceptable production flux to product.
        timeout : int
            The time limit [seconds] per attempted prediction.

        Returns
        -------
        list
            A list of pathways (list of reactions)
        """
        pathways = list()
        with TimeMachine() as tm:
            tm(do=partial(setattr, self.model.solver.configuration, 'timeout', timeout),
               undo=partial(setattr, self.model.solver.configuration, 'timeout', self.model.solver.configuration.timeout))
            try:
                demand_reaction = self.model.reactions.get_by_id('DM_' + product.id)
            except KeyError:
                demand_reaction = self.model.add_demand(product)
            tm(do=str, undo=partial(self.model.remove_reactions, [demand_reaction]))
            demand_reaction.lower_bound = min_production
            # self.model.solver.configuration.verbosity = 3
            counter = 1
            while counter <= max_predictions:
                logger.debug('Predicting pathway No. %d' % counter)
                try:
                    self.model.solve()
                except SolveError:
                    logger.debug('No pathway could be predicted. Terminating pathway predictions.')
                    break
                vars_to_cut = list()
                for y_var in self._y_vars:
                    if y_var.primal == 1.0:
                        vars_to_cut.append(y_var)
                if len(vars_to_cut) == 0:  # no pathway found:
                    logger.debug("It seems %s is a native product in model %s. No further predictions are attempted." % (product, self.model))
                    break
                pathway = [self.model.reactions.get_by_id(y_var.name[2:]) for y_var in vars_to_cut]
                # Figure out adapter reactions to include

                for adapter in self._adpater_reactions:
                    if abs(adapter.variable.primal) > 1e-6:
                        logger.debug("Adapter %s added to pathway with primal %e" % (adapter.id, adapter.variable.primal))
                        pathway.append(adapter)
                pathways.append(pathway)
                integer_cut = self.model.solver.interface.Constraint(Add(*vars_to_cut), name="integer_cut_" + str(counter), ub=len(vars_to_cut)-1)
                logger.debug('Adding integer cut.')
                tm(do=partial(self.model.solver._add_constraint, integer_cut),
                   undo=partial(self.model.solver._remove_constraint, integer_cut))
                counter += 1
            # self.model.solver.configuration.verbosity = 0
        return pathways



if __name__ == '__main__':

    from cameo.api import hosts
    pathway_predictor = PathwayPredictor(hosts.ecoli.models.EcoliCore)
    print pathway_predictor.run(product=pathway_predictor.model.metabolites.MNXM53)  # MNXM53 = L-serine
