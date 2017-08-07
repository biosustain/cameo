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

import pickle
import re
from os.path import join

import pytest

from cameo import load_model
from cameo.config import solvers
from cameo.core.pathway import Pathway
from cameo.flux_analysis.analysis import PhenotypicPhasePlaneResult
from cameo.strain_design.pathway_prediction import PathwayPredictor
from cameo.strain_design.pathway_prediction.pathway_predictor import \
    PathwayResult
from cameo.util import TimeMachine


@pytest.fixture(scope="module", params=list(solvers))
def pathway_predictor(request, data_directory, universal_model):
    core_model = load_model(join(data_directory, 'EcoliCore.xml'), sanitize=False)
    core_model.solver = request.param
    predictor = PathwayPredictor(core_model, universal_model=universal_model)
    return core_model, predictor


@pytest.fixture(scope="function")
def pathway_predictor_result(pathway_predictor):
    core_model, predictor = pathway_predictor
    return core_model, predictor.run(product='L-Serine', max_predictions=1)


class TestPathwayPredictor:
    def test_incorrect_arguments_raises(self, pathway_predictor):
        model, _ = pathway_predictor
        with pytest.raises(ValueError) as excinfo:
            PathwayPredictor(model, universal_model='Mickey_Mouse')
        assert re.search(r'Provided universal_model.*', str(excinfo.value))
        with pytest.raises(ValueError) as excinfo:
            PathwayPredictor(model, compartment_regexp='Mickey_Mouse')

    def test_predict_non_native_compound(self, pathway_predictor):
        model, predictor = pathway_predictor
        result = predictor.run(product='L-Serine', max_predictions=1)
        assert len(result) == 1
        assert len(result.pathways) == 1
        assert len(result.pathways[0].reactions) == 3
        assert len(result.pathways[0].adapters) == 0

    def test_pathway_predictor_benchmark(self, benchmark, pathway_predictor):
        model, predictor = pathway_predictor
        benchmark(predictor.run, product='L-Serine', max_predictions=1)

    def test_contains_right_adapters_and_exchanges(self, pathway_predictor):
        model, predictor = pathway_predictor
        result = predictor.run(product='L-Serine', max_predictions=1)
        for pathway in result:
            for reaction in pathway.reactions:
                for metabolite in reaction.metabolites:
                    try:
                        model.metabolites.get_by_id(metabolite.id)
                    except KeyError:
                        metabolite_ids = [met.id in adapter for adapter in pathway.adapters
                                          for met in adapter.metabolites]

                        metabolite_ids += [met.id in exchange for exchange in pathway.exchanges
                                           for met in exchange.metabolites]
                        for r in pathway.reactions:
                            if r != reaction:
                                metabolite_ids += [met.id for met in r.metabolites]

                        metabolite_ids += [met.id for met in pathway.product.metabolites]

                        assert metabolite.id in metabolite_ids


class TestPathwayResult:
    def test_pathway(self, pathway_predictor_result):
        model, result = pathway_predictor_result
        pathway = result[0]
        biomass = 'Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2'
        assert isinstance(pathway, PathwayResult)
        assert isinstance(pathway, Pathway)
        assert isinstance(pathway.production_envelope(model, objective=biomass), PhenotypicPhasePlaneResult)
        assert pathway.needs_optimization(model, objective=biomass)

    def test_pickle_pathway(self, pathway_predictor_result):
        model, result = pathway_predictor_result
        dump = pickle.dumps(result[0])
        result_recovered = pickle.loads(dump)

        assert set(r.id for r in result[0].reactions) == set(r.id for r in result_recovered.reactions)
        assert set(r.id for r in result[0].targets) == set(r.id for r in result_recovered.targets)
        assert set(r.id for r in result[0].adapters) == set(r.id for r in result_recovered.adapters)
        assert set(r.id for r in result[0].exchanges) == set(r.id for r in result_recovered.exchanges)
        assert result[0].product.id == result_recovered.product.id

    def test_plug_model_without_context(self, pathway_predictor_result):
        model, result = pathway_predictor_result
        model = model.copy()
        result[0].plug_model(model)
        for reaction in result[0].reactions:
            assert reaction in model.reactions

        for reaction in result[0].exchanges:
            assert reaction in model.reactions

        for reaction in result[0].adapters:
            assert reaction in model.reactions

    def test_plug_model_with_context(self, pathway_predictor_result):
        model, result = pathway_predictor_result
        with model:
            result[0].plug_model(model)

        for reaction in result[0].reactions:
            assert reaction not in model.reactions

        for reaction in result[0].exchanges:
            assert reaction not in model.reactions

        for reaction in result[0].adapters:
            assert reaction not in model.reactions
