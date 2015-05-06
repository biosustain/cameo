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

__all__ = ['PathwayPredictor']

from functools import partial

from cameo.core.result import Result
from cameo.exceptions import SolveError
from cameo import Reaction, Model, Metabolite
from cameo.data import metanetx
from cameo.util import TimeMachine
from sympy import Add

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class PathwayPredictions(Result):

    def __init__(self, pathways, *args, **kwargs):
        super(PathwayPredictions, self).__init__(*args, **kwargs)
        # TODO: sort the pathways to make them easier to read
        self.pathways = pathways
        print(self.pathways)

    def _repr_html_(self):
        raise NotImplementedError

    def __str__(self):
        string = str()
        for i, pathway in enumerate(self.pathways):
            string += 'Pathway No. {}'.format(i + 1)
            for reaction in pathway:
                string += '{}, {}:'.format(reaction.id, reaction.name, reaction.build_reaction_string(use_metabolite_names=True))
        return string

    @property
    def plot(self):
        # TODO: small pathway visualizations would be great.
        raise NotImplementedError


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
            from cameo.data.universal_models import metanetx_universal_model_bigg_rhea
            self.universal_model = metanetx_universal_model_bigg_rhea
        elif isinstance(universal_model, Model):
            self.universal_model = universal_model
        else:
            raise ValueError('Provided universal_model %s is not a model.' % universal_model)
        self.products = self.universal_model.metabolites
        self.model = model.copy()
        try:
            logger.debug('Trying to set solver to cplex to speed up pathway predictions.')
            self.model.solver = 'cplex'
        except ValueError:
            logger.debug('cplex not available for pathway predictions.')
        for exchange in self.model.exchanges:
            if len(exchange.reactants) > 0 and exchange.lower_bound <= 0:
                exchange.upper_bound = 999999.

        logger.debug("Adding reactions from universal model to host model.")
        reactions = list()
        for reaction in self.universal_model.reactions:
            if reaction in self.model.reactions:
                continue
            reactions.append(reaction.copy())
        self.model.add_reactions(reactions)

        logger.debug("Adding adapter reactions to connect model with universal model.")
        self._adpater_reactions = list()
        if mapping is None:
            self.mapping = metanetx.all2mnx
        for metabolite in model.metabolites:  # model is the original host model
            name = metabolite.id[0:-2]
            try:
                mnx_name = self.mapping[name]
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
        y_vars = list()
        switches = list()
        for reaction in reactions:
            if reaction.id.startswith('DM_'):
                continue  # demand reactions don't need integer switches
            y = self.model.solver.interface.Variable('y_'+reaction.id, lb=0, ub=1, type='binary')
            y_vars.append(y)
            if reaction.reverse_variable is not None:
                flux_term = reaction.variable - reaction.reverse_variable
            else:
                flux_term = reaction.variable
            switch_lb = self.model.solver.interface.Constraint(y*reaction.lower_bound - flux_term, name='switch_lb_'+reaction.id, ub=0)
            switches.append(switch_lb)
            switch_ub = self.model.solver.interface.Constraint(y*reaction.upper_bound - flux_term, name='switch_ub_'+reaction.id, lb=0)
            switches.append(switch_ub)
        self.model.solver.add(switches)
        logger.debug("Setting minimization of switch variables as objective.")
        self.model.objective = self.model.solver.interface.Objective(Add(*y_vars), direction='min')
        self._y_vars_ids = [var.name for var in y_vars]

    def run(self, product=None, max_predictions=float("inf"), min_production=.1, timeout=None):
        """Run pathway prediction for a desired product.

        Parameters
        ----------
        product : Metabolite, str
            Metabolite or id or name of metabolite to find production pathways for.
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
        if isinstance(product, bytes):
            flag = False
            for metabolite in self.model.metabolites:
                if metabolite.id == product:
                    product = metabolite
                    flag = True
                    break
                if metabolite.name == product:
                    product = metabolite
                    flag = True
                    break
            if not flag:
                raise ValueError("Specified product '{product}' could not be found. Try searching pathway_predictor_obj.universal_metabolites.metabolites".format(product=product))
        elif isinstance(product, Metabolite):
            try:
                product = self.model.metabolites.get_by_id(product.id)
            except KeyError:
                raise ValueError('Provided product %s cannot be found in universal reaction database.' % product)
        else:
            raise ValueError('Provided product %s is neither a metabolite nor an ID or name.' % product)
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
            counter = 1
            while counter <= max_predictions:
                logger.debug('Predicting pathway No. %d' % counter)
                try:
                    solution = self.model.solve()
                except SolveError:
                    logger.debug('No pathway could be predicted. Terminating pathway predictions.')
                    break
                vars_to_cut = list()
                for i, y_var_id in enumerate(self._y_vars_ids):
                    y_var = self.model.solver.variables[y_var_id]
                    if y_var.primal == 1.0:
                        vars_to_cut.append(y_var)
                if len(vars_to_cut) == 0:  # no pathway found:
                    logger.info("It seems %s is a native product in model %s. Let's see if we can find better heterologous pathways." % (product, self.model))
                    df =  solution.data_frame
                    # knockout adapter with native product
                    for adapter in self._adpater_reactions:
                        if product in adapter.metabolites:
                            logger.debug('Knocking out adapter reaction %s containing native product.' % adapter)
                            adapter.knock_out(time_machine=tm)
                    continue
                pathway = [self.model.reactions.get_by_id(y_var.name[2:]) for y_var in vars_to_cut]
                logger.debug('Pathway predicted: %s' % '\t'.join([reaction.annotation['Description'] for reaction in pathway]))
                # Figure out adapter reactions to include
                for adapter in self._adpater_reactions:
                    if abs(adapter.variable.primal) > 1e-6:
                        # logger.debug('Adapter %s added to pathway with primal %e' % (adapter.id, adapter.variable.primal))
                        pathway.append(adapter)
                pathways.append(pathway)
                integer_cut = self.model.solver.interface.Constraint(Add(*vars_to_cut), name="integer_cut_" + str(counter), ub=len(vars_to_cut)-1)
                logger.debug('Adding integer cut.')
                tm(do=partial(self.model.solver._add_constraint, integer_cut),
                   undo=partial(self.model.solver._remove_constraint, integer_cut))
                counter += 1
            # self.model.solver.configuration.verbosity = 0
        return PathwayPredictions(pathways)

    def _predict_heterolgous_pathways_for_native_compound(self, product):
        # Determined reactions that produce product natively and knock them out
        pass

if __name__ == '__main__':

    from cameo.api import hosts
    pathway_predictor = PathwayPredictor(hosts.ecoli.models.EcoliCore)
    print(pathway_predictor.run(product=pathway_predictor.model.metabolites.MNXM53))  # MNXM53 = L-serine
