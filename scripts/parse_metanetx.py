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

import cPickle as pickle
import re
import optlang
from pandas import read_table
from cameo.io import _apply_sanitize_rules, ID_SANITIZE_RULES_TAB_COMPLETION, ID_SANITIZE_RULES_SIMPHENY
from cameo import Reaction, Metabolite, Model
from cobra.core.Formula import Formula


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
reac_prop = read_table('../data/metanetx/reac_prop.tsv.gz', skiprows=107, compression='gzip', index_col=0)
reac_prop.columns = [name.replace('#', '') for name in reac_prop.columns]
chem_prop = read_table('../data/metanetx/chem_prop.tsv.gz', skiprows=125, compression='gzip', index_col=0, names=['name', 'formula', 'charge', 'mass', 'InChI', 'SMILES', 'source'])


REVERSE_ID_SANITIZE_RULES_SIMPHENY = [(value, key) for key, value in ID_SANITIZE_RULES_SIMPHENY]

# Metabolites
bigg_selection = chem_xref[['bigg' in blub for blub in chem_xref.XREF]]
sanitized_XREF = [_apply_sanitize_rules(_apply_sanitize_rules(id.replace('bigg:', ''), REVERSE_ID_SANITIZE_RULES_SIMPHENY), ID_SANITIZE_RULES_TAB_COMPLETION) for id in bigg_selection.XREF]
bigg2mnx = dict(zip(sanitized_XREF, bigg_selection.MNX_ID))
mnx2bigg = dict(zip(bigg_selection.MNX_ID, sanitized_XREF))

# Reactions
bigg_selection = reac_xref[['bigg' in blub for blub in reac_xref.XREF]]
sanitized_XREF = [_apply_sanitize_rules(_apply_sanitize_rules(id.replace('bigg:', ''), REVERSE_ID_SANITIZE_RULES_SIMPHENY), ID_SANITIZE_RULES_TAB_COMPLETION) for id in bigg_selection.XREF]
bigg2mnx.update(dict(zip(sanitized_XREF, bigg_selection.MNX_ID)))
mnx2bigg.update(dict(zip(bigg_selection.MNX_ID, sanitized_XREF)))

# Select which reactions to include in universal reaction database
reaction_selection = reac_prop[[(('bigg:' in source) or ('rhea:' in source)) and re.match('.*biomass.*', source, re.I) is None for source in reac_prop.Source]]
reactions = list()
for index, row in reaction_selection.iterrows():
    try:
        stoichiometry = parse_reaction(row.Equation, rev_arrow='=')
    except ValueError:
        continue
    else:
        for met, coeff in stoichiometry.iteritems():
            met.name = chem_prop.loc[met.id]['name']
            try:
                met.formula = Formula(chem_prop.loc[met.id].formula)
            except:
                pass
            try:
                met.charge = int(chem_prop.loc[met.id].charge)
            except:
                pass
            rest = chem_prop.loc[met.id].to_dict()
            met.annotation = dict((key, rest[key]) for key in rest if key in ('mass', 'InChI', 'source'))
            # print met.id
            # print met.name
            # print met.annotation
        mets = [met.id for met in stoichiometry.keys()]
        if len(mets) != len(set(mets)):
            continue
        reaction = Reaction(index)
        reaction.add_metabolites(stoichiometry)
        if reaction.check_mass_balance() != []:
            continue
        if row.Balance:
            reaction.lower_bound = -1*reaction.upper_bound
        print row['Source']
        reaction.name = row['Source']
        rest = row.to_dict()
        reaction.annotation = dict((key, rest[key]) for key in rest if key in ('EC', 'Description'))
        reactions.append(reaction)

metanetx_model = Model('metanetx_universal_model_bigg_rhea', solver_interface=optlang.interface)
metanetx_model.add_reactions(reactions)
# Add sinks for all metabolites
for metabolite in metanetx_model.metabolites:
    metanetx_model.add_demand(metabolite)

# Select which reactions to include in universal reaction database
reaction_selection = reac_prop[[(('bigg:' in source) or ('rhea:' in source) or ('kegg:' in source)) and re.match('.*biomass.*', source, re.I) is None for source in reac_prop.Source]]
reactions = list()
for index, row in reaction_selection.iterrows():
    try:
        stoichiometry = parse_reaction(row.Equation, rev_arrow='=')
    except ValueError:
        continue
    else:
        for met, coeff in stoichiometry.iteritems():
            met.name = chem_prop.loc[met.id]['name']
            try:
                met.formula = Formula(chem_prop.loc[met.id].formula)
            except:
                pass
            try:
                met.charge = int(chem_prop.loc[met.id].charge)
            except:
                pass
            rest = chem_prop.loc[met.id].to_dict()
            met.annotation = dict((key, rest[key]) for key in rest if key in ('mass', 'InChI', 'source'))
            # print met.id
            # print met.name
            # print met.annotation
        mets = [met.id for met in stoichiometry.keys()]
        if len(mets) != len(set(mets)):
            continue
        reaction = Reaction(index)
        reaction.add_metabolites(stoichiometry)
        if reaction.check_mass_balance() != []:
            continue
        if row.Balance:
            reaction.lower_bound = -1*reaction.upper_bound
        print row['Source']
        reaction.name = row['Source']
        rest = row.to_dict()
        reaction.annotation = dict((key, rest[key]) for key in rest if key in ('EC', 'Description'))
        reactions.append(reaction)

metanetx_model2 = Model('metanetx_universal_model_bigg_rhea_kegg', solver_interface=optlang.interface)
metanetx_model2.add_reactions(reactions)
# Add sinks for all metabolites
for metabolite in metanetx_model2.metabolites:
    metanetx_model2.add_demand(metabolite)

# Store all relevant metanetx data in a pickled dictionary
metanetx = dict()
metanetx['universal_model'] = metanetx_model
metanetx['universal_model_low_quality'] = metanetx_model2
metabolite_ids = [metabolite.id for metabolite in metanetx_model.metabolites]
metanetx['chem_prop'] = chem_prop.loc[metabolite_ids]
metanetx['bigg2mnx'] = bigg2mnx
metanetx['mnx2bigg'] = mnx2bigg

# with gzip.open('../data/metanetx.pgz', 'w') as f:
#     pickle.dump(metanetx, f)

with open('../cameo/data/metanetx.pickle', 'w') as f:
    pickle.dump(metanetx, f)
