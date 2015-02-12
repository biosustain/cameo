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


__all__ = ['reactions_to_network']

import networkx as nx


def model_to_network(model, network_type=nx.DiGraph):
    return reactions_to_network(model.reactions, network_type=nx.DiGraph)

def reactions_to_network(reactions, network_type=nx.DiGraph):
    edges = list()
    for reaction in reactions:
        for substrate in reaction.reactants:
            edges.append((substrate, reaction))
        for product in reaction.products:
            edges.append((reaction, product))
    if reaction.reversibility:
        edges.extend([(node2, node1) for node1, node2 in edges])
    return nx.DiGraph(edges)