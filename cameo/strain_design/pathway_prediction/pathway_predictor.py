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

import logging
import re
import warnings
from collections import Counter
from functools import partial
from math import ceil
from copy import copy

from cobra import DictList
from sympy import Add, Mul, RealNumber

from cobra import Model, Metabolite, Reaction
from cobra.util import SolverNotFound
from cobra.exceptions import OptimizationError

from cameo import fba
from cameo import models, phenotypic_phase_plane
from cameo.config import non_zero_flux_threshold
from cameo.core.pathway import Pathway
from cameo.core.result import Result, MetaInformation
from cameo.core.strain_design import StrainDesignMethodResult, StrainDesign, StrainDesignMethod
from cameo.core.target import ReactionKnockinTarget
from cameo.data import metanetx
from cameo.strain_design.pathway_prediction import util
from cameo.util import TimeMachine
from cameo.visualization.plotting import plotter

__all__ = ['PathwayPredictor']

logger = logging.getLogger(__name__)

add = Add._from_args
mul = Mul._from_args


class PathwayResult(Pathway, Result, StrainDesign):
    def __init__(self, reactions, exchanges, adapters, product, *args, **kwargs):
        self._meta_information = MetaInformation()
        self.reactions = reactions
        self.exchanges = exchanges
        self.adapters = adapters
        self.product = product
        self.targets = self._build_targets()

    def _replace_adapted_metabolites(self, reaction):
        """
        Replace adapted metabolites by model metabolites

        Parameters
        ----------
        reaction: cameo.core.reaction.Reaction

        Returns
        -------
        cameo.core.reaction.Reaction
        """
        stoichiometry = {}

        for metabolite, coefficient in reaction.metabolites.items():
            found = False
            for adapter in self.adapters:
                if metabolite == adapter.products[0]:
                    metabolite = adapter.reactants[0]
                    found = False
                    break
            if not found:
                metabolite = metabolite

            stoichiometry[metabolite] = coefficient

        reaction = Reaction(id=reaction.id,
                            name=reaction.name,
                            lower_bound=reaction.lower_bound,
                            upper_bound=reaction.upper_bound)
        reaction.add_metabolites(stoichiometry)

        return reaction

    def _build_targets(self):
        targets = DictList()
        for reaction in self.reactions:
            reaction = self._replace_adapted_metabolites(reaction)
            if reaction.id in metanetx.mnx2all:
                target = ReactionKnockinTarget(reaction.id, reaction, accession_id=reaction.id, accession_db='metanetx')
            else:
                target = ReactionKnockinTarget(reaction.id, reaction)
            targets.append(target)

        for reaction in self.exchanges:
            reaction = self._replace_adapted_metabolites(reaction)
            targets.append(ReactionKnockinTarget(reaction.id, reaction))

        product = self._replace_adapted_metabolites(self.product)
        product.lower_bound = 0
        targets.append(ReactionKnockinTarget(product.id, product))

        return targets

    def plot(self, **kwargs):
        pass

    def needs_optimization(self, model, objective=None):
        area = self.production_envelope(model, objective).area
        return area > 1e-5

    def production_envelope(self, model, objective=None):
        with model:
            self.apply(model)
            return phenotypic_phase_plane(model, variables=[objective], objective=self.product.id)

    def plug_model(self, model, adapters=True, exchanges=True):
        warnings.warn("The 'plug_model' method as been deprecated. Use apply instead.", DeprecationWarning)
        model.add_reactions(self.reactions)
        if adapters:
            model.add_reactions(self.adapters)
        if exchanges:
            model.add_reactions(self.exchanges)
        try:
            model.add_reaction(self.product)
        except Exception:
            logger.warning("Exchange %s already in model" % self.product.id)
            pass
        self.product.lower_bound = 0


class PathwayPredictions(StrainDesignMethodResult):
    __method_name__ = "PathwayPredictor"

    def __init__(self, pathways, *args, **kwargs):
        super(PathwayPredictions, self).__init__(pathways, *args, **kwargs)

    @property
    def pathways(self):
        return self._designs

    def plug_model(self, model, index):
        warnings.warn("The 'plug_model' method as been deprecated. You can use result[i].apply instead",
                      DeprecationWarning)
        self.pathways[index].plug_model(model)

    def __getitem__(self, item):
        return self.pathways[item]

    def __str__(self):
        string = str()
        for i, pathway in enumerate(self.pathways):
            string += 'Pathway No. {}'.format(i + 1)
            for reaction in pathway.reactions:
                string += '{}, {}: {}'.format(
                    reaction.id, reaction.name,
                    reaction.build_reaction_string(use_metabolite_names=True)
                )
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


class PathwayPredictor(StrainDesignMethod):
    """Pathway predictions from a universal set of reaction.

    Parameters
    ----------
    model : cobra.Model
        The model that represents the host organism.
    universal_model : cobra.Model, optional
        The model that represents the universal set of reactions.
        A default model will be used if omitted.
    mapping : dict, optional
        A dictionary that contains a mapping between metabolite
        identifiers in `model` and `universal_model`
    compartment_regexp : str, optional
        A regular expression that matches host metabolites' compartments
        that should be connected to the universal reaction model. If not
        provided, the compartment containing most metabolites will be
        chosen.

    Attributes
    ----------
    model : cobra.Model
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
            compartments_tally = Counter(metabolite.compartment for metabolite in self.original_model.metabolites)
            most_common_compartment = compartments_tally.most_common(n=1)[0][0]
            compartment_regexp = re.compile('^' + most_common_compartment + '$')
        else:
            compartment_regexp = re.compile(compartment_regexp)

        if universal_model is None:
            logger.debug("Loading default universal model.")
            self.universal_model = models.universal.metanetx_universal_model_bigg
        elif isinstance(universal_model, Model):
            self.universal_model = universal_model
        else:
            raise ValueError('Provided universal_model %s is not a model.' % universal_model)
        self.products = self.universal_model.metabolites

        if mapping is None:
            self.mapping = metanetx.all2mnx
        else:
            self.mapping = mapping

        self.model = model.copy()

        try:
            logger.info('Trying to set solver to cplex to speed up pathway predictions.')
            self.model.solver = 'cplex'
        except SolverNotFound:
            logger.info('cplex not available for pathway predictions.')

        self.new_reactions = self._extend_model(model.boundary)

        logger.debug("Adding adapter reactions to connect model with universal model.")
        self.adpater_reactions = util.create_adapter_reactions(model.metabolites, self.universal_model,
                                                               self.mapping, compartment_regexp)
        self.model.add_reactions(self.adpater_reactions)
        self._add_switches(self.new_reactions)

    def run(self, product=None, max_predictions=float("inf"), min_production=.1,
            timeout=None, callback=None, silent=False, allow_native_exchanges=False):
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
        callback : function
            A function that takes a successfully predicted pathway.
        silent : bool
            If True will print the pathways and max flux values.
        allow_native_exchanges: bool
            If True, exchange reactions for native metabolites will be allowed.

        Returns
        -------
        PathwayPredictions
            The predicted pathways.
        """

        product = self._find_product(product)

        pathways = list()
        with TimeMachine() as tm, self.model:
            tm(do=partial(setattr, self.model.solver.configuration, 'timeout', timeout),
               undo=partial(setattr, self.model.solver.configuration, 'timeout',
                            self.model.solver.configuration.timeout))
            try:
                product_reaction = self.model.reactions.get_by_id('DM_' + product.id)
            except KeyError:
                product_reaction = self.model.add_boundary(product, type='demand')

            product_reaction.lower_bound = min_production
            pathway_counter = 1
            integer_cut_counter = 1
            while pathway_counter <= max_predictions:
                logger.debug('Predicting pathway No. %d' % pathway_counter)
                try:
                    self.model.slim_optimize(error_value=None)
                except OptimizationError as err:
                    logger.error('No pathway could be predicted. Terminating pathway predictions.')
                    logger.error(err)
                    break

                vars_to_cut = list()
                for i, y_var_id in enumerate(self._y_vars_ids):
                    y_var = self.model.solver.variables[y_var_id]
                    if y_var.primal == 1.0:
                        vars_to_cut.append(y_var)
                logger.info(vars_to_cut)

                if len(vars_to_cut) == 0:
                    # no pathway found:
                    logger.info("It seems %s is a native product in model %s. "
                                "Let's see if we can find better heterologous pathways.", product, self.model)
                    # knockout adapter with native product
                    for adapter in self.adpater_reactions:
                        if product in adapter.metabolites:
                            logger.info('Knocking out adapter reaction %s '
                                        'containing native product.', adapter)
                            adapter.knock_out()
                    continue

                pathway = [self.model.reactions.get_by_id(y_var.name[2:]) for y_var in vars_to_cut]

                pathway_metabolites = set([m for pathway_reaction in pathway for m in pathway_reaction.metabolites])
                logger.info('Pathway predicted: %s', '\t'.join(
                    [r.build_reaction_string(use_metabolite_names=True) for r in pathway]))
                pathway_metabolites.add(product)

                # Figure out adapter reactions to include
                adapters = [adapter for adapter in self.adpater_reactions if adapter.products[0] in pathway_metabolites]

                # Figure out exchange reactions to include
                exchanges = [exchange for exchange in self._exchanges
                             if abs(exchange.flux) > non_zero_flux_threshold and exchange.id != product_reaction.id]

                if allow_native_exchanges:
                    exchanges = [exchange for exchange in exchanges
                                 if list(exchange.metabolites)[0] in pathway_metabolites]

                pathway = PathwayResult(pathway, exchanges, adapters, product_reaction)
                if not silent:
                    util.display_pathway(pathway, pathway_counter)

                integer_cut = self.model.solver.interface.Constraint(Add(*vars_to_cut),
                                                                     name="integer_cut_" + str(integer_cut_counter),
                                                                     ub=len(vars_to_cut) - 1)
                logger.debug('Adding integer cut.')
                tm(
                    do=partial(self.model.solver.add, integer_cut),
                    undo=partial(self.model.solver.remove, integer_cut))

                # Test pathway in the original model
                with self.original_model:
                    pathway.apply(self.original_model)
                    self.original_model.objective = pathway.product.id
                    try:
                        production_flux = self.original_model.slim_optimize(error_value=None)
                    except OptimizationError as err:
                        logger.error(err)
                        logger.error(
                            "Addition of pathway %r made the model unsolvable. "
                            "Skipping pathway.", pathway)
                        continue
                    else:
                        if production_flux > non_zero_flux_threshold:
                            pathways.append(pathway)
                            logger.info("Max flux: %.5G", production_flux)
                            pathway_counter += 1
                            if callback is not None:
                                callback(pathway)
                        else:
                            logger.warning(
                                "Pathway %r could not be verified. Production "
                                "flux %.5G is below the requirement %.5G. "
                                "Skipping.", pathway, production_flux,
                                non_zero_flux_threshold)
                    finally:
                        integer_cut_counter += 1

            return PathwayPredictions(pathways)

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
            forward_term = mul((RealNumber(-1), reaction.forward_variable))
            reverse_term = mul((RealNumber(-1), reaction.reverse_variable))
            switch_lb_term = mul((RealNumber(reaction.lower_bound), y))
            switch_ub_term = mul((RealNumber(reaction.upper_bound), y))
            switch_lb = self.model.solver.interface.Constraint(add((switch_lb_term, forward_term, reverse_term)),
                                                               name='switch_lb_' + reaction.id, ub=0, sloppy=True)
            switch_ub = self.model.solver.interface.Constraint(add((switch_ub_term, forward_term, reverse_term)),
                                                               name='switch_ub_' + reaction.id, lb=0, sloppy=True)
            switches.extend([switch_lb, switch_ub])

        self.model.solver.add(y_vars)
        self.model.solver.add(switches, sloppy=True)

        logger.info("Setting minimization of switch variables as objective.")
        self.model.objective = self.model.solver.interface.Objective(Add(*y_vars), direction='min')
        self._y_vars_ids = [var.name for var in y_vars]

    def _extend_model(self, original_exchanges):
        for exchange in self.model.boundary:
            if len(exchange.reactants) > 0 >= exchange.lower_bound:
                exchange.upper_bound = 999999.

        logger.info("Adding reactions from universal model to host model.")
        new_reactions = list()
        original_model_metabolites = [self.mapping.get('bigg:' + m.id[0:-2], m.id) for
                                      r in original_exchanges for m, coeff in r.metabolites.items()
                                      if len(r.metabolites) == 1 and coeff < 0 < r.upper_bound]

        universal_exchanges = self.universal_model.boundary
        for reaction in self.universal_model.reactions:
            if reaction in self.model.reactions:
                continue
            if reaction in universal_exchanges:
                metabolite = list(reaction.metabolites.keys())[0]
                if metabolite.id in original_model_metabolites:
                    continue

            new_reactions.append(copy(reaction))
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
                "Specified product '{product}' could not be found. "
                "Try searching pathway_predictor_obj.universal_model.metabolites".format(product=product))
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
