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


__all__ = ['model_to_network', 'reactions_to_network', 'remove_highly_connected_nodes']

import networkx as nx


def model_to_network(model, network_type=nx.Graph):
    """Convert a model into a networkx graph.

    Parameters
    ----------
    model: SolverBasedModel
        The model.
    network_type: networkx.Graph or or other networkx graph types, optional (default networkx.Graph)
        The type of networkx graph that should be returned.

    Returns
    -------
    networkx.Graph (default)
        Depends on network_type parameter.
    """
    return reactions_to_network(model.reactions, network_type=network_type)

def reactions_to_network(reactions, network_type=nx.Graph):
    """Convert a list of reactions into a networkx graph.

    Parameters
    ----------
    reactions: list
        The list of reactions.
    network_type: networkx.Graph or or other networkx graph types, optional (default networkx.Graph)
        The type of networkx graph that should be returned.

    Returns
    -------
    networkx.Graph (default)
        Depends on network_type parameter.
    """
    edges = list()
    for reaction in reactions:
        for substrate in reaction.reactants:
            edges.append((substrate, reaction))
        for product in reaction.products:
            edges.append((reaction, product))
    if reaction.reversibility:
        edges.extend([(node2, node1) for node1, node2 in edges])
    return network_type(edges)

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