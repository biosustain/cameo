# Copyright 2015 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""This module implements a pathway data structure that can be added to a model for example."""

from functools import partial, reduce

from pandas import DataFrame

from cobra import Metabolite, Reaction


# TODO: Load pathways from SBML and JSON
# TODO: Define the product
# TODO: Visualization with ESCHER
class Pathway(object):
    """
    Representation of a pathway (a set of reactions)

    Attributes
    ----------

    reactions: list o Reaction
        The list of reactions in the pathway.

    """

    def __init__(self, reactions, *args, **kwargs):
        super(Pathway, self).__init__(*args, **kwargs)
        self.reactions = reactions

    @property
    def data_frame(self):
        return DataFrame([[r.build_reaction_string(True), r.lower_bound, r.upper_bound] for r in self.reactions],
                         columns=["equation", "lower_bound", "upper_bound"], index=[r.id for r in self.reactions])

    @classmethod
    def from_file(cls, file_path, sep="\t"):
        """
        Read a pathway from a file.

        The file format is:
        reaction_id<sep>equation<sep>lower_limit<sep>upper_limit<sep>name<sep>comments\n

        The equation is defined by:
        coefficient * substrate_name#substrate_id + ... <=> coefficient * product_name#product_id

        Parameters
        ----------
        file_path: str
            The path to the file containing the pathway
        sep: str
            The separator between elements in the file (default: "\t")

        Returns
        -------
        Pathway

        """
        reactions = []
        metabolites = {}
        with open(file_path, "r") as pathway:
            for line in pathway:
                if line.startswith("#"):
                    continue
                line = line.strip("\n")

                identifier, equation, lower, upper, name, comments = line.split(sep)
                stoichiometry = _parse_equation(equation, metabolites)
                reaction = _build_reaction(identifier, stoichiometry, float(upper), float(lower), name, comments)
                reactions.append(reaction)

        return cls(reactions)

    def to_file(self, file_path, sep="\t"):
        """
        Writes the pathway to a file.

        Parameters
        ----------
        file_path: str
            The path to the file where the pathway will be written
        sep: str
            The separator between elements in the file (default: "\t")

        See Also
        --------
        Pathway.from_file

        """
        with open(file_path, "w") as output_file:
            for reaction in self.reactions:
                equation = _build_equation(reaction.metabolites)
                output_file.write(sep.join(map(str, [
                    reaction.id,
                    equation,
                    reaction.lower_bound,
                    reaction.upper_bound,
                    reaction.name,
                    reaction.notes.get("pathway_note", "")
                ])) + "\n")

    def plug_model(self, model):
        """
        Plug the pathway to a model.

        Metabolites are matched in the model by id. For metabolites with no ID in the model, an exchange reaction
        is added to the model

        Parameters
        ----------
        model: cobra.Model
            The model to plug in the pathway

        """
        model.add_reactions(self.reactions)
        metabolites = set(reduce(lambda x, y: x + y, [list(r.metabolites.keys()) for r in self.reactions], []))
        exchanges = [model.add_boundary(m) for m in metabolites if len(m.reactions) == 1]
        for exchange in exchanges:
            exchange.lower_bound = 0


def _build_reaction(identifier, stoichiometry, upper, lower, name, comments):
    reaction = Reaction(identifier)
    reaction.id = identifier
    reaction.name = name
    reaction.add_metabolites(stoichiometry)
    reaction.upper_bound = upper
    reaction.lower_bound = lower
    reaction.notes["pathway_note"] = comments
    return reaction


def _parse_equation(equation, metabolites):
    stoichiometry = {}
    reactants, products = equation.split(" <=> ")
    for reactant in reactants.split(" + "):
        coeff, metabolite = reactant.split(" * ")
        name, metabolite_id = metabolite.split("#")
        if metabolite_id in metabolites:
            met = metabolites[metabolite_id]
        else:
            metabolites[metabolite_id] = met = Metabolite(id=metabolite_id, name=name)
        stoichiometry[met] = -float(coeff)

    for product in products.split(" + "):
        coeff, metabolite = product.split(" * ")
        name, metabolite_id = metabolite.split("#")
        if metabolite_id in metabolites:
            met = metabolites[metabolite_id]
        else:
            metabolites[metabolite_id] = met = Metabolite(id=metabolite_id, name=name)
        stoichiometry[met] = float(coeff)

    return stoichiometry


def _build_equation(stoichiometry):
    products = {m: v for m, v in stoichiometry.items() if v > 0}
    reactants = {m: v for m, v in stoichiometry.items() if v < 0}

    products = " + ".join(["%f * %s#%s" % (v, m.name, m.id) for m, v in products.items()])
    reactants = " + ".join(["%f * %s#%s" % (-v, m.name, m.id) for m, v in reactants.items()])
    return reactants + " <=> " + products
