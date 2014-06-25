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


def bpcy(biomass, product, substrate):
    def f(model, solution, decoded_representation):
        try:
            biomass_flux = solution.get_primal_by_id(biomass)
            product_flux = solution.get_primal_by_id(product)
            substrate_flux = abs(solution.get_primal_by_id(substrate))
            return (biomass_flux * product_flux) / substrate_flux

        except ZeroDivisionError, e:
            return 0

    f.__name__ = "bpcy = (%s * %s) / %s" % (biomass, product, substrate)

    return f


def product_yield(product, substrate):
    def f(model, solution, decoded_representation):
        try:
            return solution.get_primal_by_id(product) / -solution.get_primal_by_id(substrate)
        except ZeroDivisionError:
            return 0

    f.__name__ = "yield = (%s / %s)" % (product, substrate)

    return f


def number_of_knockouts(sense='min'):
    def f(model, solution, decoded_representation):
        if sense == 'max':
            return len(decoded_representation[1])
        else:
            return 1 / len(decoded_representation[1])

    if sense == max:
        f.__name__ = "max #knockouts"
    else:
        f.__name__ = "min #knockouts"

    return f
