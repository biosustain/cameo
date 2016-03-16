# Copyright 2015 Novo Nordisk Foundation Center for Biosustainability, DTU.
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

__all__ = ['distance_based_on_molecular_formula']


def distance_based_on_molecular_formula(metabolite1, metabolite2, normalize=True):
    """Calculate the distance of two metabolites bases on the molecular formula

    Arguments
    ---------
    metabolite1 : Metabolite
        The first metabolite.
    metabolite2 : Metabolite
        The second metabolite.
    normalize : bool, optional
        If the distance should be normalized by the total number of elements in both metabolites (defaults to True).

    Returns
    -------
    float
        The distance between metabolite1 and metabolite2.
    """
    if len(metabolite1.formula.elements) == 0 or len(metabolite2.formula.elements) == 0:
        raise ValueError('Cannot calculate distance between metabolites %s and %s' % (metabolite1, metabolite2))
    elements = set(list(metabolite1.formula.elements.keys()) + list(metabolite2.formula.elements.keys()))
    distance = 0.
    for element in elements:
        distance += abs(metabolite1.formula.elements.get(element, 0) - metabolite2.formula.elements.get(element, 0))
    if normalize:
        return distance / sum(list(metabolite1.formula.elements.values()) + list(metabolite2.formula.elements.values()))
    else:
        return distance
