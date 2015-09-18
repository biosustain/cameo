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

from functools import partial, reduce
from cobra.core import Metabolite
from pandas import DataFrame
import six
from cameo.core.reaction import Reaction


class Pathway(object):
    def __init__(self, reactions, *args, **kwargs):
        super(Pathway, self).__init__(*args, **kwargs)
        self.reactions = reactions

    @property
    def data_frame(self):
        return DataFrame([[r.build_reaction_string(True), r.lower_bound, r.upper_bound] for r in self.reactions],
                         columns=["equation", "lower_bound", "upper_bound"], index=[r.id for r in self.reactions])

    @classmethod
    def from_file(cls, file_path, sep="\t"):
        reactions = []
        with open(file_path, "r") as pathway:
            for line in pathway:
                if line.startswith("#"):
                    continue
                line = line.strip("\n")

                identifier, equation, lower, upper, name, comments = line.split(sep)
                stoichiometry = _parse_equation(equation)
                reaction = _build_reaction(identifier, stoichiometry, float(upper), float(lower), name, comments)
                reactions.append(reaction)

        return cls(reactions)

    def to_file(self, file_path, sep="\t"):
        with open(file_path, "w") as output_file:
            for reaction in self.reactions:
                equation = _build_equation(reaction.metabolites)
                output_file.write(reaction.id + sep +
                                  equation + sep +
                                  reaction.lower_bound + sep +
                                  reaction.upper_bound + sep +
                                  reaction.name + sep +
                                  reaction.notes.get("pathway_note", "") + "\n")

    def plug_model(self, model, tm=None):
        metabolites = reduce(lambda x, y: x+y, [list(r.metabolites.keys()) for r in self.reactions], [])
        exchanges = [model.add_demand(m, prefix="EX_", time_machine=tm) for m in metabolites if m not in model.metabolites]
        for exchange in exchanges:
            exchange.lower_bound = 0
        if tm is not None:
            tm(do=partial(model.add_reactions, self.reactions),
               undo=partial(model.remove_reactions, self.reactions, delete=False))
        else:
            model.add_reactions(self.reactions)


def _build_reaction(identifier, stoichiometry, upper, lower, name, comments):
    reaction = Reaction(identifier)
    reaction.name = name
    reaction.add_metabolites(stoichiometry)
    reaction.upper_bound = upper
    reaction.lower_bound = lower
    reaction.notes["pathway_note"] = comments
    return reaction


def _parse_equation(equation):
    stoichiometry = {}
    reactants, products = equation.split(" <=> ")
    for reactant in reactants.split(" + "):
        coeff, metabolite = reactant.split(" * ")
        name, metabolite_id = metabolite.split("#")
        met = Metabolite(id=metabolite_id, name=name)
        stoichiometry[met] = -float(coeff)

    for product in products.split(" + "):
        coeff, metabolite = product.split(" * ")
        name, metabolite_id = metabolite.split("#")
        met = Metabolite(id=metabolite_id, name=name)
        stoichiometry[met] = float(coeff)

    return stoichiometry


def _build_equation(stoichiometry):
    products = {m: v for m, v in six.iteritems(stoichiometry) if v > 0}
    reactants = {m: v for m, v in six.iteritems(stoichiometry) if v > 0}

    products = " + ".join(["%f * %s#%s" % (v, m.name, m.id) for m, v in products])
    reactants = " + ".join(["%f * %s#%s" % (-v, m.name, m.id) for m, v in reactants])
    return reactants + " <=> " + products