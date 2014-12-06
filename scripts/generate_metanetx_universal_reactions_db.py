# Copyright 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.
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

from pandas import read_table
from cobra.io import write_sbml_model
from cameo.io import _apply_sanitize_rules, ID_SANITIZE_RULES_TAB_COMPLETION
from cameo import Reaction, Metabolite, Model

def parse_reaction(formula, irrev_arrow='-->', rev_arrow='<=>'):

    def parse_rhs_side(string):
        return parse_one_side(string, factor=1.)

    def parse_lhs_side(string):
        return parse_one_side(string, factor=-1.)

    def parse_one_side(string, factor=1.):

        def parse_coeff_and_metabolites(term):
            try:
                coeff, metabolite_id = term.strip().split(' ')
            except:
                raise ValueError('Something is fishy with the provided term %s' % term)
            return Metabolite(metabolite_id), factor*float(coeff)

        terms = string.split('+')
        return dict([parse_coeff_and_metabolites(term) for term in terms])

    try:
        (lhs, rhs) = formula.split(irrev_arrow)
    except ValueError:
        try:
            (lhs, rhs) = formula.split(rev_arrow)
        except ValueError:
            raise ValueError('Something is wrong with the provided reaction formula %s' % formula)
    stoichiometry = parse_lhs_side(lhs)
    stoichiometry.update(parse_rhs_side(rhs))
    return stoichiometry

chem_xref = read_table('../data/metanetx/chem_xref.tsv.gz', skiprows=124, compression='gzip')
chem_xref.columns = [name.replace('#', '') for name in chem_xref.columns]
reac_xref = read_table('../data/metanetx/reac_xref.tsv.gz', skiprows=107, compression='gzip')
reac_xref.columns = [name.replace('#', '') for name in reac_xref.columns]
reac_prop = read_table('../data/metanetx/reac_prop.tsv.gz', skiprows=107, compression='gzip')
reac_prop.columns = [name.replace('#', '') for name in reac_prop.columns]

bigg_selection = chem_xref[['bigg' in blub for blub in chem_xref.XREF]]
sanitized_XREF = [_apply_sanitize_rules(id, ID_SANITIZE_RULES_TAB_COMPLETION) for id in bigg_selection.XREF]
bigg2mnx = dict(zip(sanitized_XREF, bigg_selection.MNX_ID))
mnx2bigg = dict(zip(bigg_selection.MNX_ID, sanitized_XREF))

reaction_selection = reac_prop[[('bigg:' in source) or ('rhea:' in source) for source in reac_prop.Source]]
reactions = list()
for i, row in reaction_selection.iterrows():
    try:
        stoichiometry = parse_reaction(row.Equation, rev_arrow='=')
    except ValueError:
        continue
    else:
#         for k, v in stoichiometry.items():
#             try:
#                 bigg_id = mnx2bigg[k.id]
#             except KeyError:
#                 continue
#             stoichiometry.pop(k)
#             stoichiometry[Metabolite(bigg_id.replace('bigg:', ''))] = v
        mets = [met.id for met in stoichiometry.keys()]
        if len(mets) != len(set(mets)):
            continue
        reaction = Reaction(row.MNX_ID)
        reaction.add_metabolites(stoichiometry)
        if row.Balance:
            reaction.lower_bound = -1*reaction.upper_bound
        reactions.append(reaction)

metanetx_model = Model('metanetx_universal_model_bigg_rhea')
metanetx_model.add_reactions(reactions)
write_sbml_model(metanetx_model, '../data/universal_reaction_database_metanetx_bigg_rhea.xml')