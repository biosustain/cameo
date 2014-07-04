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
from cameo import config


class biomass_product_coupled_yield():
    """
    Biomass-Product Coupled Yield: (v[biomass] * v[product]) / v[substrate]

    :param biomass: biomass reaction identifier
    :param product: product reaction identifier
    :param substrate: substrate reaction identifier
    :return: fitness value
    """

    def __init__(self, biomass, product, substrate):
        self.biomass = biomass
        self.product = product
        self.substrate = substrate
        self.name = "bpcy = (%s * %s) / %s" % (biomass, product, substrate)
        self.__name__ = self.__class__.__name__

    def __call__(self, model, solution, decoded_representation):
        try:
            biomass_flux = round(solution.get_primal_by_id(self.biomass), config.ndecimals)
            product_flux = round(solution.get_primal_by_id(self.product), config.ndecimals)
            substrate_flux = round(abs(solution.get_primal_by_id(self.substrate)), config.ndecimals)
            return round((biomass_flux * product_flux) / substrate_flux, config.ndecimals)

        except ZeroDivisionError:
            return 0

def product_yield(product, substrate):
    """
    Product Yield Objective function: v[product]/v[substrate]
    :param product: product reaction identifier
    :param substrate: substrate reaction identifier
    :return: fitness value
    """
    def f(model, solution, decoded_representation):
        try:
            product_flux = round(solution.get_primal_by_id(product), config.ndecimals)
            substrate_flux = round(abs(solution.get_primal_by_id(substrate)), config.ndecimals)
            return round(product_flux / substrate_flux, config.ndecimals)
        except ZeroDivisionError:
            return 0

    f.name = "yield = (%s / %s)" % (product, substrate)

    return f


def number_of_knockouts(sense='min'):
    """
    Number of Knockouts objective function.
    If sense is maximize then fitness is the number of knockouts, otherwise 1/#knockouts

    :param sense: 'max' or 'min'
    :return: fitness value
    """
    def f(model, solution, decoded_representation):
        if sense == 'max':
            return len(decoded_representation[1])
        else:
            return round(1.0 / len(decoded_representation[1]), config.ndecimals)

    if sense == max:
        f.name = "max #knockouts"
    else:
        f.name = "min #knockouts"

    return f
