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

import pandas
import pytest
from pandas import DataFrame
from pandas.util.testing import assert_frame_equal

from cobra.exceptions import Infeasible

import cameo
from cameo import fba
from cameo.config import solvers
from cameo.strain_design.deterministic.flux_variability_based import (FSEOF,
                                                                      DifferentialFVA,
                                                                      FSEOFResult)
from cameo.strain_design.deterministic.linear_programming import OptKnock

TRAVIS = bool(os.getenv('TRAVIS', False))
TESTDIR = os.path.dirname(__file__)


@pytest.fixture(scope='module')
def cplex_optknock(model):
    cplex_core = model.copy()
    cplex_core.reactions.Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2.lower_bound = 0.1
    cplex_core.solver = "cplex"
    return cplex_core, OptKnock(cplex_core)


@pytest.fixture(scope='module')
def diff_fva(model):
    return DifferentialFVA(model, model.reactions.EX_succ_lp_e_rp_, points=5)


class TestFSEOF:
    def test_fseof(self, model):
        objective = model.objective
        fseof = FSEOF(model)
        fseof_result = fseof.run(target="EX_succ_lp_e_rp_")
        assert isinstance(fseof_result, FSEOFResult)
        assert objective.expression == model.objective.expression

    def test_fseof_result(self, model):
        fseof = FSEOF(model)
        fseof_result = fseof.run(target=model.reactions.EX_ac_lp_e_rp_)
        assert isinstance(fseof_result.data_frame, DataFrame)
        assert fseof_result.target is model.reactions.EX_ac_lp_e_rp_
        assert fseof_result.model is model


class TestDifferentialFVA:
    def test_minimal_input(self, diff_fva):
        result = diff_fva.run()
        ref_df = pandas.read_csv(os.path.join(TESTDIR, 'data/REFERENCE_DiffFVA1.csv'), index_col=0)
        this_df = result.nth_panel(0)
        this_df.index.name = None
        pandas.util.testing.assert_frame_equal(this_df[ref_df.columns], ref_df)

    def test_apply_designs(self, model, diff_fva):
        result = diff_fva.run()
        works = []
        for strain_design in result:
            with model:
                strain_design.apply(model)
                try:
                    solution = fba(model, objective="Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2")
                    works.append(solution["EX_succ_lp_e_rp_"] > 1e-6 and solution.objective_value > 1e-6)
                except Infeasible:
                    works.append(False)
        assert any(works)

    def test_diff_fva_benchmark(self, diff_fva, benchmark):
        benchmark(diff_fva.run)

    def test_with_reference_model(self, model):
        reference_model = model.copy()
        biomass_rxn = reference_model.reactions.Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2
        biomass_rxn.lower_bound = 0.3
        target = reference_model.reactions.EX_succ_lp_e_rp_
        target.lower_bound = 2
        result = DifferentialFVA(model, target, reference_model=reference_model, points=5).run()
        ref_df = pandas.read_csv(os.path.join(TESTDIR, 'data/REFERENCE_DiffFVA2.csv'), index_col=0)
        this_df = result.nth_panel(0)
        this_df.index.name = None
        pandas.util.testing.assert_frame_equal(this_df[ref_df.columns], ref_df)


@pytest.mark.skipif('cplex' not in solvers, reason="No cplex interface available")
class TestOptKnock:
    def test_optknock_runs(self, cplex_optknock):
        _, optknock = cplex_optknock
        result = optknock.run(max_knockouts=0, target="EX_ac_lp_e_rp_",
                              biomass="Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2", max_results=1)
        assert len(result) == 1
        assert len(result.knockouts[0]) == 0
        assert len(list(result)) == 1
        assert isinstance(result.data_frame, DataFrame)

    def test_optknock_benchmark(self, cplex_optknock, benchmark):
        _, optknock = cplex_optknock
        benchmark(optknock.run, max_knockouts=2, target="EX_ac_lp_e_rp_",
                  biomass="Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2", max_results=1)

    def test_result_is_correct(self, cplex_optknock):
        model, optknock = cplex_optknock
        result = optknock.run(max_knockouts=1, target="EX_ac_lp_e_rp_",
                              biomass="Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2", max_results=1)
        production = result.production[0]
        knockouts = result.knockouts[0]
        for knockout in knockouts:
            model.reactions.get_by_id(knockout).knock_out()
        fva = cameo.flux_variability_analysis(model, fraction_of_optimum=1, remove_cycles=False,
                                              reactions=["EX_ac_lp_e_rp_"])
        assert abs(fva["upper_bound"][0] - production) < 1e-6

    def test_invalid_input(self, cplex_optknock):
        _, optknock = cplex_optknock
        with pytest.raises(ValueError):
            optknock.run(target="EX_ac_lp_e_rp_")
        with pytest.raises(ValueError):
            optknock.run(biomass="Biomass_Ecoli_core_N_lp_w_fsh_GAM_rp__Nmet2")
