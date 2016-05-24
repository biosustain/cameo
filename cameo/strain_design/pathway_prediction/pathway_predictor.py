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
from math import ceil
from cameo.visualization.plotting import plotter

import six

import re
from functools import partial

from cameo.core.result import Result
from cameo.core.pathway import Pathway
from cameo import models, phenotypic_phase_plane
from cameo.exceptions import SolveError
from cameo import Model, Metabolite
from cameo.data import metanetx
from cameo.util import TimeMachine

from cameo.strain_design.strain_design import StrainDesignResult, StrainDesign
from cameo.strain_design.pathway_prediction import util

import sympy

from sympy import Add, Mul, RealNumber

import logging

__all__ = ['PathwayPredictor']

NegativeOne = sympy.singleton.S.NegativeOne


logger = logging.getLogger(__name__)


class PathwayResult(Pathway, Result, StrainDesign):
    def __init__(self, reactions, exchanges, adapters, product, *args, **kwargs):
        Result.__init__(self, *args, **kwargs)
        Pathway.__init__(self, reactions, *args, **kwargs)
        StrainDesign.__init__(self, knock_ins=[r.id for r in reactions], manipulation_type="reactions")
        self.exchanges = exchanges
        self.adapters = adapters
        self.product = product

    def plot(self, **kwargs):
        pass

    def needs_optimization(self, model, objective=None):
        area = self.production_envelope(model, objective).area
        return area > 1e-5

    def production_envelope(self, model, objective=None):
        with TimeMachine() as tm:
            self.plug_model(model, tm)
            return phenotypic_phase_plane(model, variables=[objective], objective=self.product)

    def plug_model(self, model, tm=None, adapters=True, exchanges=True):
        if tm is not None:
            tm(do=partial(model.add_reactions, self.reactions),
               undo=partial(model.remove_reactions, self.reactions, delete=False, remove_orphans=True))
            if adapters:
                tm(do=partial(model.add_reactions, self.adapters),
                   undo=partial(model.remove_reactions, self.adapters, delete=False, remove_orphans=True))
            if exchanges:
                tm(do=partial(model.add_reactions, self.exchanges),
                   undo=partial(model.remove_reactions, self.exchanges, delete=False, remove_orphans=True))
            tm(do=partial(setattr, self.product, "lower_bound", 0),
               undo=partial(setattr, self.product, "lower_bound", self.product.lower_bound))
            try:
                tm(do=partial(model.add_reaction, self.product),
                   undo=partial(model.remove_reactions, [self.product], delete=False, remove_orphans=True))
            except Exception:
                logger.warning("Exchange %s already in model" % self.product.id)
                pass
        else:
            model.add_reactions(self.reactions)
            if adapters:
                model.add_reactions(self.adapters)
            if exchanges:
                model.add_reactions(self.exchanges)

            self.product.lower_bound = 0
            try:
                model.add_reaction(self.product)
            except Exception:
                logger.warning("Exchange %s already in model" % self.product.id)
                pass


class PathwayPredictions(StrainDesignResult):
    __method_name__ = "PathwayPredictor"

    def __init__(self, pathways, *args, **kwargs):
        super(PathwayPredictions, self).__init__(*args, **kwargs)
        # TODO: sort the pathways to make them easier to read
        self.pathways = pathways

    def plug_model(self, model, index, tm=None):
        self.pathways[index].plug_model(model, tm)

    def __str__(self):
        string = str()
        for i, pathway in enumerate(self.pathways):
            string += 'Pathway No. {}'.format(i + 1)
            for reaction in pathway.reactions:
                string += '{}, {}:'.format(reaction.id, reaction.name,
                                           reaction.build_reaction_string(use_metabolite_names=True))
        return string

    def plot(self, grid=None, width=None, height=None, title=None):
        # TODO: small pathway visualizations would be great.
        raise NotImplementedError

    def plot_production_envelopes(self, model, objective=None, title=None):
        rows = int(ceil(len(self.pathways) / 2.0))
        title = "Production envelops for %s" % self.pathways[0].product.name if title is None else title
        grid = plotter.grid(n_rows=rows, title=title)
        with grid:
            for i, pathway in enumerate(self.pathways):
                ppp = pathway.production_envelope(model, objective=objective)
                ppp.plot(grid=grid, width=450, title="Pathway %i" % (i + 1))

    def __iter__(self):
        for p in self.pathways:
            yield p

    def __len__(self):
        return len(self.pathways)


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

    def __init__(self, model, universal_model=None, mapping=None, compartment_regexp=None):
        """"""
        self.original_model = model
        if compartment_regexp is None:
            compartment_regexp = re.compile(".*")

        if universal_model is None:
            logger.debug("Loading default universal model.")
            self.universal_model = models.universal.metanetx_universal_model_bigg_rhea
        elif isinstance(universal_model, Model):
            self.universal_model = universal_model
        else:
            raise ValueError('Provided universal_model %s is not a model.' % universal_model)
        self.products = self.universal_model.metabolites

        if mapping is None:
            self.mapping = metanetx.all2mnx
            for metabolite in self.original_model.metabolites:
                try:
                    self.universal_model.metabolites.get_by_id(metabolite.id)
                except KeyError:
                    pass
                else:
                    self.mapping[metabolite.id] = metabolite.id

        self.model = model.copy()

        try:
            logger.info('Trying to set solver to cplex to speed up pathway predictions.')
            self.model.solver = 'cplex'
        except ValueError:
            logger.info('cplex not available for pathway predictions.')

        self.new_reactions = self._extend_model(model.exchanges)

        logger.debug("Adding adapter reactions to connect model with universal model.")
        self.adpater_reactions = util.create_adapter_reactions(model.metabolites, self.universal_model,
                                                               self.mapping, compartment_regexp)
        self.model.add_reactions(self.adpater_reactions)
        self._add_switches(self.new_reactions)

    def run(self, product=None, max_predictions=float("inf"), min_production=.1, timeout=None, callback=None, silent=False):
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

        product = self._find_product(product)

        pathways = list()
        with TimeMachine() as tm:
            tm(do=partial(setattr, self.model.solver.configuration, 'timeout', timeout),
               undo=partial(setattr, self.model.solver.configuration, 'timeout',
                            self.model.solver.configuration.timeout))
            try:
                demand_reaction = self.model.reactions.get_by_id('DM_' + product.id)
            except KeyError:
                demand_reaction = self.model.add_demand(product)
                tm(do=str, undo=partial(self.model.remove_reactions, [demand_reaction], delete=False))
            demand_reaction.lower_bound = min_production
            counter = 1
            while counter <= max_predictions:
                logger.debug('Predicting pathway No. %d' % counter)
                try:
                    self.model.solve()
                except SolveError as e:
                    logger.error('No pathway could be predicted. Terminating pathway predictions.')
                    logger.error(e)
                    break

                vars_to_cut = list()
                for i, y_var_id in enumerate(self._y_vars_ids):
                    y_var = self.model.solver.variables[y_var_id]
                    if y_var.primal == 1.0:
                        vars_to_cut.append(y_var)
                logger.info(vars_to_cut)
                if len(vars_to_cut) == 0:
                    # no pathway found:
                    logger.info(
                        "It seems %s is a native product in model %s. Let's see if we can find better heterologous pathways." % (
                            product, self.model))
                    # knockout adapter with native product
                    for adapter in self.adpater_reactions:
                        if product in adapter.metabolites:
                            logger.info('Knocking out adapter reaction %s containing native product.' % adapter)
                            adapter.knock_out(time_machine=tm)
                    continue

                pathway = [self.model.reactions.get_by_id(y_var.name[2:]) for y_var in vars_to_cut]
                logger.info('Pathway predicted: %s' % '\t'.join(
                    [r.build_reaction_string(use_metabolite_names=True) for r in pathway]))
                # Figure out adapter reactions to include
                adapters = [adapter for adapter in self.adpater_reactions if abs(adapter.flux) != 0]

                # Figure out exchange reactions to include
                exchanges = [exchange for exchange in self._exchanges if abs(exchange.flux) != 0]

                pathway = PathwayResult(pathway, exchanges, adapters, demand_reaction)
                if not silent:
                    util.display_pathway(pathway, counter)

                pathways.append(pathway)
                if callback is not None:
                    callback(pathway)
                integer_cut = self.model.solver.interface.Constraint(Add(*vars_to_cut),
                                                                     name="integer_cut_" + str(counter),
                                                                     ub=len(vars_to_cut) - 1)
                logger.debug('Adding integer cut.')
                tm(do=partial(self.model.solver._add_constraint, integer_cut),
                   undo=partial(self.model.solver._remove_constraint, integer_cut))
                counter += 1
            # self.model.solver.configuration.verbosity = 0
            return PathwayPredictions(pathways)

    def _predict_heterolgous_pathways_for_native_compound(self, product):
        # Determined reactions that produce product natively and knock them out
        pass

    def _add_switches(self, reactions):
        logger.info("Adding switches.")
        y_vars = list()
        switches = list()
        self._exchanges = list()
        for reaction in reactions:
            if reaction.id.startswith('DM_'):
                # demand reactions don't need integer switches
                self._exchanges.append(reaction)
                continue

            y = self.model.solver.interface.Variable('y_' + reaction.id, lb=0, ub=1, type='binary')
            y_vars.append(y)
            # The following is a complicated but efficient way to write the following constraints

            # switch_lb = self.model.solver.interface.Constraint(y * reaction.lower_bound - reaction.flux_expression,
            #                                                    name='switch_lb_' + reaction.id, ub=0)
            # switch_ub = self.model.solver.interface.Constraint(y * reaction.upper_bound - reaction.flux_expression,
            #                                                    name='switch_ub_' + reaction.id, lb=0)
            forward_var_term = Mul._from_args((RealNumber(-1), reaction.forward_variable))
            reverse_var_term = Mul._from_args((RealNumber(-1), reaction.reverse_variable))
            switch_lb_y_term = Mul._from_args((RealNumber(reaction.lower_bound), y))
            switch_ub_y_term = Mul._from_args((RealNumber(reaction.upper_bound), y))
            switch_lb = self.model.solver.interface.Constraint(
                Add._from_args((switch_lb_y_term, forward_var_term, reverse_var_term)), name='switch_lb_' + reaction.id,
                ub=0, sloppy=True)
            switch_ub = self.model.solver.interface.Constraint(
                Add._from_args((switch_ub_y_term, forward_var_term, reverse_var_term)), name='switch_ub_' + reaction.id,
                lb=0, sloppy=True)
            switches.extend([switch_lb, switch_ub])
        self.model.solver.add(y_vars)
        self.model.solver.add(switches, sloppy=True)
        logger.info("Setting minimization of switch variables as objective.")
        self.model.objective = self.model.solver.interface.Objective(Add(*y_vars), direction='min')
        self._y_vars_ids = [var.name for var in y_vars]

    def _extend_model(self, original_exchanges):
        for exchange in self.model.exchanges:
            if len(exchange.reactants) > 0 >= exchange.lower_bound:
                exchange.upper_bound = 999999.

        logger.info("Adding reactions from universal model to host model.")
        new_reactions = list()
        original_model_metabolites = [self.mapping.get(m.id[0:-2], m.id) for
                                      r in original_exchanges for m, coeff in six.iteritems(r.metabolites)
                                      if len(r.metabolites) == 1 and coeff < 0 < r.upper_bound]

        universal_exchanges = self.universal_model.exchanges
        for reaction in self.universal_model.reactions:
            if reaction in self.model.reactions:
                continue
            if reaction in universal_exchanges:
                metabolite = list(reaction.metabolites.keys())[0]
                if metabolite.id in original_model_metabolites:
                    continue

            new_reactions.append(reaction.copy())
        self.model.add_reactions(new_reactions)

        return new_reactions

    def _find_product(self, product):
        if isinstance(product, str):
            for metabolite in self.model.metabolites:
                if metabolite.id == product:
                    return metabolite
                if metabolite.name == product:
                    return metabolite
            raise ValueError(
                "Specified product '{product}' could not be found. Try searching pathway_predictor_obj.universal_model.metabolites".format(
                    product=product))
        elif isinstance(product, Metabolite):
            try:
                return self.model.metabolites.get_by_id(product.id)
            except KeyError:
                raise ValueError('Provided product %s cannot be found in universal reaction database.' % product)
        else:
            raise ValueError('Provided product %s is neither a metabolite nor an ID or name.' % product)


if __name__ == '__main__':
    from cameo.api import hosts

    pathway_predictor = PathwayPredictor(hosts.ecoli.models.EcoliCore)
    print(pathway_predictor.run(product=pathway_predictor.model.metabolites.MNXM53))  # MNXM53 = L-serine
