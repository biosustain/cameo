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

import os
import unittest

from cameo import load_model
from cameo.core.pathway import Pathway
from cameo.flux_analysis.analysis import PhenotypicPhasePlaneResult
from cameo.strain_design.pathway_prediction import PathwayPredictor
from cameo.strain_design.pathway_prediction.pathway_predictor import PathwayResult
from cameo.util import TimeMachine

TESTDIR = os.path.dirname(__file__)
TESTMODEL = load_model(os.path.join(TESTDIR, 'data/EcoliCore.xml'))
UNIVERSALMODEL = load_model(os.path.join(TESTDIR, 'data/iJO1366.xml'))
UNIVERSALMODEL.remove_reactions(UNIVERSALMODEL.exchanges)

TRAVIS = os.getenv('TRAVIS', False)

PATHWAYPREDICTOR = PathwayPredictor(TESTMODEL, universal_model=UNIVERSALMODEL)


class Wrapper:
    class AbstractPathwayPredictorTestCase(unittest.TestCase):

        def test_setting_incorrect_universal_model_raises(self):
            with self.assertRaisesRegexp(ValueError, 'Provided universal_model.*'):
                PathwayPredictor(TESTMODEL, universal_model='Mickey_Mouse')

        # def test_predict_native_compound_returns_shorter_alternatives(self):
        #     result = self.pathway_predictor.run(product='Phosphoenolpyruvate', max_predictions=1)
        #     self.assertTrue(len(result.pathways) == 1)
        #     self.assertTrue(len(result.pathways[0].pathway) == 3)
        #     self.assertTrue(len(result.pathways[0].adapters) == 0)

        def test_predict_non_native_compound(self):
            result = self.pathway_predictor.run(product='L-Serine', max_predictions=1)
            self.assertTrue(len(result) == 1)
            self.assertTrue(len(result.pathways) == 1)
            self.assertTrue(len(result.pathways[0].reactions) == 3)
            self.assertTrue(len(result.pathways[0].adapters) == 0)

        def test_contains_right_adapters_and_exchanges(self):
            result = self.pathway_predictor.run(product='L-Serine', max_predictions=1)
            for pathway in result:
                for reaction in pathway.reactions:
                    for metabolite in reaction.metabolites:
                        try:
                            met = TESTMODEL.metabolites.get_by_id(metabolite.id)
                        except KeyError:
                            metabolite_ids = [met.id in adapter for adapter in pathway.adapters
                                              for met in adapter.metabolites]

                            metabolite_ids += [met.id in exchange for exchange in pathway.exchanges
                                               for met in exchange.metabolites]
                            for r in pathway.reactions:
                                if r != reaction:
                                    metabolite_ids += [met.id for met in r.metabolites]

                            metabolite_ids += [met.id for met in pathway.product.metabolites]

                            self.assertTrue(metabolite.id in metabolite_ids)


class PathwayPredictorCPLEXTestCase(Wrapper.AbstractPathwayPredictorTestCase):
    def setUp(self):
            TESTMODEL.solver = "cplex"
            self.pathway_predictor = PathwayPredictor(TESTMODEL, universal_model=UNIVERSALMODEL)


class PathwayPredictorGLPKTestCase(Wrapper.AbstractPathwayPredictorTestCase):
    def setUp(self):
            TESTMODEL.solver = "glpk"
            self.pathway_predictor = PathwayPredictor(TESTMODEL, universal_model=UNIVERSALMODEL)


class PathwayPredictionsTestCase(unittest.TestCase):
    def setUp(self):
        self.result = PATHWAYPREDICTOR.run(product='L-Serine', max_predictions=1)

    def test_pathway(self):
        model = TESTMODEL.copy()
        pathway = self.result[0]
        biomass = 'Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2'
        self.assertIsInstance(pathway, PathwayResult)
        self.assertIsInstance(pathway, Pathway)
        self.assertIsInstance(pathway.production_envelope(model, objective=biomass), PhenotypicPhasePlaneResult)
        self.assertTrue(pathway.needs_optimization(model, objective=biomass))

    def test_plug_model_without_time_machine(self):
        model = TESTMODEL.copy()
        self.result[0].plug_model(model)
        for reaction in self.result[0].reactions:
            self.assertIn(reaction, model.reactions)

        for reaction in self.result[0].exchanges:
            self.assertIn(reaction, model.reactions)

        for reaction in self.result[0].adapters:
            self.assertIn(reaction, model.reactions)

    def test_plug_model_with_time_machine(self):
        model = TESTMODEL.copy()
        with TimeMachine() as tm:
            self.result[0].plug_model(model, tm=tm)

        for reaction in self.result[0].reactions:
            self.assertNotIn(reaction, model.reactions)

        for reaction in self.result[0].exchanges:
            self.assertNotIn(reaction, model.reactions)

        for reaction in self.result[0].adapters:
            self.assertNotIn(reaction, model.reactions)