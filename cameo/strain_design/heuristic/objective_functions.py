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
    Biomass-Product Coupled Yield: (v[biomass] * v[product]) / v[substrate] [1]

    Parameters
    ----------
    biomass: str
        biomass reaction identifier
    product: str
        product reaction identifier
    substrate: str
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
            return 0.0

    def _repr_latex_(self):
        return "$$bpcy = \\frac{(%s * %s)}{%s}$$" % (self.biomass.replace("_", "\\_"), self.product.replace("_", "\\_"), self.substrate.replace("_", "\\_"))


class product_yield():
    """
    Product Yield Objective function: v[product]/v[substrate]

    Parameters
    ----------
    product: str
        product reaction identifier
    substrate: str
        substrate reaction identifier

    Returns
    -------
    float
        fitness value
    """
    def __init__(self, product, substrate):
        self.product = product
        self.substrate = substrate
        self.name = "yield = (%s / %s)" % (product, substrate)

    def __call__(self, model, solution, decoded_representation):
        try:
            product_flux = round(solution.get_primal_by_id(self.product), config.ndecimals)
            substrate_flux = round(abs(solution.get_primal_by_id(self.substrate)), config.ndecimals)
            return round(product_flux / substrate_flux, config.ndecimals)
        except ZeroDivisionError:
            return 0.0

    def _repr_latex_(self):
        return "$$yield = \\frac{%s}{%s}$$" % (self.product, self.substrate)


class number_of_knockouts():
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

    def __init__(self, sense='min'):
        self.sense = sense
        if sense == 'max':
            self.name = "max #knockouts"
        else:
            self.name = "min #knockouts"

    def __call__(self, model, solution, decoded_representation):
        if self.sense == 'max':
            return len(decoded_representation[1])
        else:
            return round(1.0 / len(decoded_representation[1]), config.ndecimals)

    def _repr_latex_(self):
        return "$$ %s\\:\\#knockouts $$" % self.sense