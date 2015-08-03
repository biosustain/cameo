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

__all__ = ['model_to_network', 'reactions_to_network', 'remove_highly_connected_nodes']

import networkx as nx
from cameo.network_analysis.util import distance_based_on_molecular_formula


def model_to_network(model, *args, **kwargs):
    """Convert a model into a networkx graph.

    Calls reactions_to_network with model.reactions.

    Parameters
    ----------
    model : SolverBasedModel
        The model.

    Returns
    -------
    networkx.MultiDiGraph

    See Also
    --------
    reactions_to_network
    """
    return reactions_to_network(model.reactions, *args, **kwargs)


def reactions_to_network(reactions, max_distance=0.3):
    """Convert a list of reactions into a networkx graph.

    Parameters
    ----------
    reactions : list
        The list of reactions.
    max_distance : float, optional
        A threshold on the normalized distance between two compounds. If distance is above this threshold, no edge between
        those compounds will be created.

    Returns
    -------
    networkx.MultiDiGraph

    See Also
    --------
    distance_based_on_molecular_formula
    """
    edges = list()
    for reaction in reactions:
        substrates = list(reaction.reactants)
        for substrate in substrates:
            products = list(reaction.products)
            for product in products:
                try:
                    distance = distance_based_on_molecular_formula(substrate, product, normalize=True)
                except ValueError:
                    distance = 0.
                if distance <= max_distance:
                    if reaction.reversibility:
                        edges.append((product, substrate, dict(reaction=reaction)))
                        edges.append((substrate, product, dict(reaction=reaction)))
                    elif reaction.lower_bound > 0:
                        edges.append((substrate, product, dict(reaction=reaction)))
                    else:
                        edges.append((product, substrate, dict(reaction=reaction)))

    multi_graph = nx.MultiDiGraph(edges)
    return multi_graph


def remove_highly_connected_nodes(network, max_degree=10, ignore=[]):
    """Remove highly connected nodes (changes network in place).

    Parameters
    ----------
    network : networkx graph
    max_degree : int (default 10)
        Remove nodes with degree > max_degree
    ignore : list
        List of nodes to ignore.

    Returns
    -------
    None
    """
    to_remove = [node for node, degree in network.degree_iter() if degree > max_degree and node not in ignore]
    network.remove_nodes_from(to_remove)
