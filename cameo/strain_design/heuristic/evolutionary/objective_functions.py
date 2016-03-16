# -*- coding: utf-8 -*-
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

from functools import partial

from cobra import Reaction
from cameo import config, flux_variability_analysis
from cameo.util import TimeMachine

__all__ = ['biomass_product_coupled_yield', 'product_yield', 'number_of_knockouts']


class ObjectiveFunction(object):
    """
    Blueprint for objective function.

    All objective functions must extend this class and override the __call__ method.
    Because they implement the __call__ method, they will behave like functions.

    Methods
    -------
    __call__(model, solution, decoded_representation)
        Calculates the fitness of the solution

    """

    def __init__(self, *args, **kwargs):
        super(ObjectiveFunction, self).__init__(*args, **kwargs)

    def __call__(self, model, solution, decoded_representation):
        raise NotImplementedError

    @property
    def reactions(self):
        """
        Returns the reactions ids that require flux data for this evaluation function.
        """
        return []

    def _repr_latex_(self):
        return self.name

    @property
    def name(self):
        return self.__class__.__name__


class biomass_product_coupled_yield(ObjectiveFunction):
    """
    Biomass-Product Coupled Yield: (v[biomass] * v[product]) / v[substrate] [1]

    Parameters
    ----------
    biomass: str or Reaction
        biomass reaction identifier
    product: str or Reaction
        product reaction identifier
    substrate: str or Reaction
        substrate reaction identifier

    Returns
    -------
    float
        fitness value

    References
    ----------
    [1] Patil, K. R., Rocha, I., Förster, J., & Nielsen, J. (2005). "Evolutionary programming as a
    platform for in silico metabolic engineering". BMC Bioinformatics, 6, 308.
    doi:10.1186/1471-2105-6-308
    """

    def __init__(self, biomass, product, substrate, *args, **kwargs):
        super(biomass_product_coupled_yield, self).__init__(*args, **kwargs)
        if isinstance(biomass, Reaction):
            biomass = biomass.id
        self.biomass = biomass
        if isinstance(product, Reaction):
            product = product.id
        self.product = product
        if isinstance(substrate, Reaction):
            substrate = substrate.id
        self.substrate = substrate
        self.__name__ = self.__class__.__name__

    def __call__(self, model, solution, decoded_representation):
        try:
            biomass_flux = round(solution.fluxes[self.biomass], config.ndecimals)
            product_flux = round(solution.fluxes[self.product], config.ndecimals)
            substrate_flux = round(abs(solution.fluxes[self.substrate]), config.ndecimals)
            return round((biomass_flux * product_flux) / substrate_flux, config.ndecimals)

        except ZeroDivisionError:
            return 0.0

    def _repr_latex_(self):
        return "$$bpcy = \\frac{(%s * %s)}{%s}$$" % (
            self.biomass.replace("_", "\\_"), self.product.replace("_", "\\_"), self.substrate.replace("_", "\\_"))

    @property
    def name(self):
        return "bpcy = (%s * %s) / %s" % (self.biomass, self.product, self.substrate)

    @property
    def reactions(self):
        return [self.biomass, self.product, self.substrate]


class biomass_product_coupled_min_yield(ObjectiveFunction):
    """
    Biomass-Product Coupled Minimum Yield: (v[biomass] * min(v[product])) / v[substrate] [1]

    Parameters
    ----------
    biomass: str or Reaction
        biomass reaction identifier
    product: str or Reaction
        product reaction identifier
    substrate: str or Reaction
        substrate reaction identifier

    Returns
    -------
    float
        fitness value

    """

    def __init__(self, biomass, product, substrate, *args, **kwargs):
        super(biomass_product_coupled_min_yield, self).__init__(*args, **kwargs)
        if isinstance(biomass, Reaction):
            biomass = biomass.id
        self.biomass = biomass
        if isinstance(product, Reaction):
            product = product.id
        self.product = product
        if isinstance(substrate, Reaction):
            substrate = substrate.id
        self.substrate = substrate
        self.__name__ = self.__class__.__name__

    def __call__(self, model, solution, decoded_representation):
        try:
            biomass_flux = round(solution.fluxes[self.biomass], config.ndecimals)
            with TimeMachine() as tm:
                biomass_reaction = model.reactions.get_by_id(self.biomass)
                tm(do=partial(setattr, biomass_reaction, 'lower_bound', biomass_flux),
                   undo=partial(setattr, biomass_reaction, 'lower_bound', biomass_reaction.lower_bound))
                # for reaction_ko in decoded_representation[0]:
                #     reaction_ko.knock_out(tm)

            fva_res = flux_variability_analysis(model, reactions=[self.product])

            min_product_flux = fva_res["lower_bound"][self.product]
            substrate_flux = round(abs(solution.fluxes[self.substrate]), config.ndecimals)
            return round((biomass_flux * min_product_flux) / substrate_flux, config.ndecimals)

        except ZeroDivisionError:
            return 0.0

    def _repr_latex_(self):
        return "$$bpcy = \\frac{(%s * min(%s))}{%s}$$" % (
            self.biomass.replace("_", "\\_"), self.product.replace("_", "\\_"), self.substrate.replace("_", "\\_"))

    @property
    def name(self):
        return "bpcy = (%s * min(%s)) / %s" % (self.biomass, self.product, self.substrate)

    @property
    def reactions(self):
        return [self.biomass, self.product, self.substrate]


class product_yield(ObjectiveFunction):
    """
    Product Yield Objective function: v[product]/v[substrate]

    Parameters
    ----------
    product: str or Reaction
        product reaction identifier
    substrate: str or Reaction
        substrate reaction identifier

    Returns
    -------
    float
        fitness value
    """

    def __init__(self, product, substrate, *args, **kwargs):
        super(product_yield, self).__init__(*args, **kwargs)
        if isinstance(product, Reaction):
            product = product.id
        self.product = product
        if isinstance(substrate, Reaction):
            substrate = substrate.id
        self.substrate = substrate

    def __call__(self, model, solution, decoded_representation):
        try:
            product_flux = round(solution.fluxes[self.product], config.ndecimals)
            substrate_flux = round(abs(solution.fluxes[self.substrate]), config.ndecimals)
            return round(product_flux / substrate_flux, config.ndecimals)
        except ZeroDivisionError:
            return 0.0

    def _repr_latex_(self):
        return "$$yield = \\frac{%s}{%s}$$" % (self.product, self.substrate)

    @property
    def name(self):
        return "yield = (%s / %s)" % (self.product, self.substrate)

    @property
    def reactions(self):
        return [self.product, self.substrate]


class number_of_knockouts(ObjectiveFunction):
    """
    Number of Knockouts objective function.
    If sense is maximize then fitness is the number of knockouts, otherwise 1/#knockouts

    Parameters
    ----------
    sense: str
        'max' or 'min'

    Returns
    -------
    float
        fitness value
    """

    def __init__(self, sense='min', *args, **kwargs):
        super(number_of_knockouts, self).__init__(*args, **kwargs)
        self.sense = sense

    def __call__(self, model, solution, decoded_representation):
        if self.sense == 'max':
            return len(decoded_representation[1])
        else:
            return round(1.0 / len(decoded_representation[1]), config.ndecimals)

    def _repr_latex_(self):
        return "$$ %s\\:\\#knockouts $$" % self.sense

    @property
    def name(self):
        if self.sense == 'max':
            return "max knockouts"
        else:
            return "min knockouts"
