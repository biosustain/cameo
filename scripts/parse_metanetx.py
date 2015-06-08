# Copyright 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import cPickle as pickle
import re
import gzip
import optlang
from pandas import read_table
from cameo.io import _apply_sanitize_rules, ID_SANITIZE_RULES_TAB_COMPLETION, ID_SANITIZE_RULES_SIMPHENY
from cameo import Reaction, Metabolite, Model
from cobra.core.Formula import Formula

import logging

logger = logging.getLogger('parse_metanetx')
# logger.setLevel(logging.DEBUG)


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
            return Metabolite(metabolite_id), factor * float(coeff)

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


def construct_universal_model(list_of_db_prefixes):
    # Select which reactions to include in universal reaction database

    reaction_selection = reac_prop[[
        any([source.startswith(db_prefix) for db_prefix in list_of_db_prefixes]) and re.match('.*biomass.*', source,
                                                                                              re.I) is None for source
        in reac_prop.Source]]
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
                    logger.debug('Cannot parse formula %s. Skipping formula' % chem_prop.loc[met.id].formula)
                    continue
                # if met.formula.weight is None:
                #     logger.debug('Cannot calculate weight for formula %s. Skipping reaction %s' % (met.formula, row.Equation))
                #     # print('Cannot calculate weight for formula %s. Skipping reaction %s' % (met.formula, row.Equation))
                #     continue
                try:
                    met.charge = int(chem_prop.loc[met.id].charge)
                except ValueError:
                    logger.debug('Cannot parse charge %s. Skipping charge' % chem_prop.loc[met.id].charge)
                    pass
                rest = chem_prop.loc[met.id].to_dict()
                met.annotation = dict((key, rest[key]) for key in rest if key in ('mass', 'InChI', 'source'))
            mets = [met.id for met in stoichiometry.keys()]
            if len(mets) != len(set(mets)):
                continue
            reaction = Reaction(index)
            reaction.add_metabolites(stoichiometry)
            if reaction.check_mass_balance() != []:
                continue
            if row.Balance:
                reaction.lower_bound = -1 * reaction.upper_bound
            reaction.name = row['Source']
            rest = row.to_dict()
            reaction.annotation = dict((key, rest[key]) for key in rest if key in ('EC', 'Description'))
            reactions.append(reaction)

    model = Model('metanetx_universal_model_' + '_'.join(list_of_db_prefixes), solver_interface=optlang.interface)
    model.add_reactions(reactions)
    # Add sinks for all metabolites
    for metabolite in model.metabolites:
        model.add_demand(metabolite)
    return model


if __name__ == '__main__':

    import logging

    logging.basicConfig(level='INFO')

    # load metanetx data
    chem_xref = read_table('../data/metanetx/chem_xref.tsv.gz', skiprows=124, compression='gzip')
    chem_xref.columns = [name.replace('#', '') for name in chem_xref.columns]
    reac_xref = read_table('../data/metanetx/reac_xref.tsv.gz', skiprows=107, compression='gzip')
    reac_xref.columns = [name.replace('#', '') for name in reac_xref.columns]
    reac_prop = read_table('../data/metanetx/reac_prop.tsv.gz', skiprows=107, compression='gzip', index_col=0)
    reac_prop.columns = [name.replace('#', '') for name in reac_prop.columns]
    chem_prop = read_table('../data/metanetx/chem_prop.tsv.gz', skiprows=125, compression='gzip', index_col=0,
                           names=['name', 'formula', 'charge', 'mass', 'InChI', 'SMILES', 'source'])

    REVERSE_ID_SANITIZE_RULES_SIMPHENY = [(value, key) for key, value in ID_SANITIZE_RULES_SIMPHENY]

    metanetx = dict()
    # Metabolites
    bigg_selection = chem_xref[['bigg' in blub for blub in chem_xref.XREF]]
    sanitized_XREF = [
        _apply_sanitize_rules(_apply_sanitize_rules(id.replace('bigg:', ''), REVERSE_ID_SANITIZE_RULES_SIMPHENY),
                              ID_SANITIZE_RULES_TAB_COMPLETION) for id in bigg_selection.XREF]
    bigg2mnx = dict(zip(sanitized_XREF, bigg_selection.MNX_ID))
    mnx2bigg = dict(zip(bigg_selection.MNX_ID, sanitized_XREF))

    # Reactions
    bigg_selection = reac_xref[['bigg' in blub for blub in reac_xref.XREF]]
    sanitized_XREF = [
        _apply_sanitize_rules(_apply_sanitize_rules(id.replace('bigg:', ''), REVERSE_ID_SANITIZE_RULES_SIMPHENY),
                              ID_SANITIZE_RULES_TAB_COMPLETION) for id in bigg_selection.XREF]
    bigg2mnx.update(dict(zip(sanitized_XREF, bigg_selection.MNX_ID)))
    mnx2bigg.update(dict(zip(bigg_selection.MNX_ID, sanitized_XREF)))

    # put into final result dict
    metanetx['bigg2mnx'] = bigg2mnx
    metanetx['mnx2bigg'] = mnx2bigg

    all2mnx = dict()
    for other_id, mnx_id in chem_xref[['XREF', 'MNX_ID']].values:
        cleaned_key = _apply_sanitize_rules(
            _apply_sanitize_rules(other_id.split(':')[1], REVERSE_ID_SANITIZE_RULES_SIMPHENY),
            ID_SANITIZE_RULES_TAB_COMPLETION)
        all2mnx[cleaned_key] = mnx_id
    for other_id, mnx_id in reac_xref[['XREF', 'MNX_ID']].values:
        cleaned_key = _apply_sanitize_rules(
            _apply_sanitize_rules(other_id.split(':')[1], REVERSE_ID_SANITIZE_RULES_SIMPHENY),
            ID_SANITIZE_RULES_TAB_COMPLETION)
        all2mnx[cleaned_key] = mnx_id

    metanetx['all2mnx'] = all2mnx
    with open('../cameo/data/metanetx.pickle', 'wb') as f:
        pickle.dump(metanetx, f)

    # generate universal reaction models
    db_combinations = [('bigg',), ('rhea',) , ('bigg', 'rhea'), ('bigg', 'rhea', 'kegg'), ('bigg', 'rhea', 'kegg', 'brenda')]
    for db_combination in db_combinations:
        universal_model = construct_universal_model(db_combination)
        with open('../cameo/data/universal_models/{model_name}.pickle'.format(model_name=universal_model.id) , 'wb') as f:
            pickle.dump(universal_model, f)

    chem_prop_filtered = chem_prop[[any([source.startswith(db) for db in ('bigg', 'rhea', 'kegg', 'brenda', 'chebi')]) for source in chem_prop.source]]
    chem_prop_filtered = chem_prop_filtered.dropna(subset=['name'])
    with gzip.open('../cameo/data/metanetx_chem_prop.pklz','wb') as f:
        pickle.dump(chem_prop_filtered, f)