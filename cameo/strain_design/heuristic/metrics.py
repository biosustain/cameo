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

from __future__ import absolute_import

from math import sqrt


def euclidean_distance(wt, mutant):
    return sqrt(sum([(wt[r] - mutant[r]) ** 2 for r in list(wt.keys())]))


def manhattan_distance(wt, mutant):
    return sum([(wt[r] - mutant[r]) for r in list(wt.keys())])
