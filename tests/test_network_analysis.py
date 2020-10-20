# -*- coding: utf-8 -*-
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

from cameo.network_analysis import (model_to_network, reactions_to_network,
                                    remove_highly_connected_nodes)
from cameo.network_analysis.util import \
    distance_based_on_molecular_formula as dbmf


def test_distance_based_on_molecular_formula(salmonella):
    assert dbmf(salmonella.metabolites[0], salmonella.metabolites[0], normalize=False) == 0
    assert dbmf(salmonella.metabolites[0], salmonella.metabolites[0], normalize=True) == 0
    assert dbmf(salmonella.metabolites[0], salmonella.metabolites[1], normalize=False) == 58.0
    assert dbmf(salmonella.metabolites[0], salmonella.metabolites[1], normalize=True) == 0.6590909090909091


def test_model_to_network(salmonella):
    assert model_to_network(salmonella).edges == reactions_to_network(salmonella.reactions).edges
    network = model_to_network(salmonella)
    assert len(network.nodes) == 1761
    assert len(network.edges) == 4924
    network = model_to_network(salmonella, max_distance=1.)
    # nodes = network.nodes()
    # print(set(nodes).difference(set(core_model_one.metabolites)))
    # print(set(core_model_one.metabolites).difference(set(nodes)))
    assert len(network.nodes) == 1800
    assert len(network.edges) == 12853


def test_remove_highly_connected_nodes(salmonella):
    network = model_to_network(salmonella)
    assert salmonella.metabolites.atp_c in network.nodes()
    assert salmonella.metabolites.adp_c in network.nodes()
    remove_highly_connected_nodes(network, max_degree=10, ignore=[salmonella.metabolites.atp_c])
    assert len(network.nodes) == 1671
    assert len(network.edges) == 2342
    assert salmonella.metabolites.atp_c in network.nodes()
    assert salmonella.metabolites.adp_c not in network.nodes()
