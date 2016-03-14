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

from cobra.test import create_test_model
from nose.tools import assert_equal, assert_in, assert_not_in

from cameo.network_analysis import model_to_network, reactions_to_network, remove_highly_connected_nodes
from cameo.network_analysis.util import distance_based_on_molecular_formula

TEST_MODEL = create_test_model()


def test_distance_based_on_molecular_formula():
    assert_equal(
        distance_based_on_molecular_formula(TEST_MODEL.metabolites[0], TEST_MODEL.metabolites[0], normalize=False), 0)
    assert_equal(
        distance_based_on_molecular_formula(TEST_MODEL.metabolites[0], TEST_MODEL.metabolites[0], normalize=True), 0)
    assert_equal(
        distance_based_on_molecular_formula(TEST_MODEL.metabolites[0], TEST_MODEL.metabolites[1], normalize=False),
        58.0)
    assert_equal(
        distance_based_on_molecular_formula(TEST_MODEL.metabolites[0], TEST_MODEL.metabolites[1], normalize=True),
        0.6590909090909091)


def test_model_to_network():
    assert_equal(model_to_network(TEST_MODEL).edges(), reactions_to_network(TEST_MODEL.reactions).edges())
    network = model_to_network(TEST_MODEL)
    assert_equal(len(network.nodes()), 1761)
    assert_equal(len(network.edges()), 4924)
    network = model_to_network(TEST_MODEL, max_distance=1.)
    # nodes = network.nodes()
    # print(set(nodes).difference(set(TEST_MODEL.metabolites)))
    # print(set(TEST_MODEL.metabolites).difference(set(nodes)))
    assert_equal(len(network.nodes()), 1800)
    assert_equal(len(network.edges()), 12853)


def test_remove_highly_connected_nodes():
    network = model_to_network(TEST_MODEL)
    assert_in(TEST_MODEL.metabolites.atp_c, network.nodes())
    assert_in(TEST_MODEL.metabolites.adp_c, network.nodes())
    remove_highly_connected_nodes(network, max_degree=10, ignore=[TEST_MODEL.metabolites.atp_c])
    assert_equal(len(network.nodes()), 1671)
    assert_equal(len(network.edges()), 2342)
    assert_in(TEST_MODEL.metabolites.atp_c, network.nodes())
    assert_not_in(TEST_MODEL.metabolites.adp_c, network.nodes())


if __name__ == '__main__':
    import nose

    nose.runmodule()
