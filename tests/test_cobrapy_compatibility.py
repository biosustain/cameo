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

import types

from cobra.test import create_test_model
from cobra.test.flux_analysis import TestCobraFluxAnalysis
from cobra.test.unit_tests import CobraTestCase, TestReactions

from cameo.core.solver_based_model import to_solver_based_model, SolverBasedModel


def setUp(self):
    # Make Model pickable and then load a solver based version of test_pickle
    self.model = to_solver_based_model(create_test_model())
    self.model_class = SolverBasedModel


for cls in (CobraTestCase, TestReactions, TestCobraFluxAnalysis):
    cls.setUp = types.MethodType(setUp, cls)

del TestCobraFluxAnalysis.test_single_gene_deletion
del TestCobraFluxAnalysis.test_phenotype_phase_plane  # Avoid bug in cobra tests (AttributeError on self.skip() )
