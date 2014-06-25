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
    def f(solution):
        try:
            return (solution.x_dict[biomass] * solution.x_dict[product]) / -solution.x_dict[substrate]

        except ZeroDivisionError:
            return 0

    f.__name__ = "bpcy = (%s * %s) / %s" % (biomass, product, substrate)

    return f


def product_yield(product, substrate):
    def f(solution):
        try:
            return solution[product] / -solution[substrate]
        except ZeroDivisionError:
            return 0

    f.__name__ = "yield = (%s / %s)" % (product, substrate)

    return f


def number_of_knockouts(sense='min'):
    def f(solution):
        if sense == 'max':
            return solution.knockouts
        else:
            return 1/solution.knockouts

    if sense == max:
        f.__name__ = "max #knockouts"
    else:
        f.__name__ = "min #knockouts"

    return f
