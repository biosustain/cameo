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

from cobra import Metabolite
import six

__all__ = ['model_to_network', 'reactions_to_network', 'remove_highly_connected_nodes']

import networkx as nx
from cameo.network_analysis.util import distance_based_on_molecular_formula


def model_to_network(model, distance_function=distance_based_on_molecular_formula, exchanges=False):
    """Convert a model into a networkx graph.

    Parameters
    ----------
    model: SolverBasedModel
        The model.

    Returns
    -------
    networkx.MultiDiGraph
    """
    return reactions_to_network(model.reactions, distance_function=distance_function, exchanges=exchanges)


def reactions_to_network(reactions, distance_function=distance_based_on_molecular_formula, exchanges=False):
    """Convert a list of reactions into a networkx graph.

    Parameters
    ----------
    reactions: list
        The list of reactions.

    Returns
    -------
    networkx.MultiDiGraph
    """
    r_edges = list()
    metabolite_nodes = set()

    e_edges = list()
    exchange_nodes = set()

    for reaction in reactions:
        for substrate in reaction.reactants:
            for product in reaction.products:
                try:
                    distance = distance_function(substrate, product)
                except ValueError:
                    distance = 0.
                if distance <= 0.3:
                    r_edges.append((substrate, product, dict(reaction=reaction)))
                    if reaction.reversibility:
                        r_edges.append((product, substrate, dict(reaction=reaction)))
                metabolite_nodes.add(product)
            metabolite_nodes.add(substrate)
        if len(reaction.metabolites) == 1 and exchanges:
            met = six.next(six.iterkeys(reaction.metabolites))
            met_ex = Metabolite(id=met.id, formula=met.formula)
            e_edges.append((met, met_ex, dict(reaction=reaction)))
            exchange_nodes.add(met_ex)

    multi_graph = nx.MultiDiGraph()

    multi_graph.add_nodes_from([(m, dict(id=m.id)) for m in metabolite_nodes], color="#088da5", type="metabolite")
    multi_graph.add_nodes_from([(m, dict(id=m.id)) for m in exchange_nodes], color="#fdeae0", type="exchange")

    multi_graph.add_edges_from(r_edges)
    multi_graph.add_edges_from(e_edges)
    return multi_graph


def remove_highly_connected_nodes(network, max_degree=10, ignore=[]):
    """Remove highly connected nodes.

    Parameters
    ----------
    network: networkx graph
    max_degree: int (default 10)
        Remove nodes with degree > max_degree
    ignore: list
        List of nodes to ignore.
    """
    to_remove = [node for node, degree in network.degree_iter() if degree > max_degree and node not in ignore]
    network.remove_nodes_from(to_remove)