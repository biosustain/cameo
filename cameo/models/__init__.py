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

"""This sub-package provides access (an function internet connection is needed
to models available on the BiGG database (http://bigg.ucsd.edu) and a model repository
hosted by the University of Minho (http://optflux.org/models).
Furthermore, it provides universal reaction database models compiled from different
sources (KEGG, RHEA, BIGG, BRENDA) and that have been obtained from http://metanetx.org"""

from .webmodels import *
from .universal import *
