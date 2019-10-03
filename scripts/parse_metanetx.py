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

import sys

import logging
import gzip
import json
import sys

import requests
from cobra.core.formula import Formula
from cobra.io.json import save_json_model
from pandas import read_table, notnull

from cobra import Reaction, Metabolite, Model
from cameo.io import _apply_sanitize_rules, ID_SANITIZE_RULES_TAB_COMPLETION, ID_SANITIZE_RULES_SIMPHENY

logger = logging.getLogger('parse_metanetx')


# logger.setLevel(logging.DEBUG)


def parse_reaction(formula, irrev_arrow='-->', rev_arrow='<=>'):
    """Parse a metanetx reaction formula."""
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


def construct_universal_model(list_of_db_prefixes, reac_xref, reac_prop, chem_prop):
    """"Construct a universal model based on metanetx.

    Parameters
    ----------
    list_of_db_prefixes : list
        A list of database prefixes, e.g., ['bigg', 'rhea']
    reac_xref : pandas.DataFrame
        A dataframe of http://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_xref.tsv
    reac_prop : pandas.DataFrame
        A dataframe of http://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_prop.tsv
    chem_prop : pandas.DataFrame
        A dataframe of http://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_prop.tsv

    """
    # Select which reactions to include in universal reaction database

    mnx_reaction_id_selection = set()
    for db_prefix in list_of_db_prefixes:
        mnx_reaction_id_selection.update(reac_xref[reac_xref.XREF.str.startswith(db_prefix)].MNX_ID)

    reaction_selection = reac_prop.loc[mnx_reaction_id_selection]
    reactions = list()
    for index, row in reaction_selection.iterrows():
        try:
            stoichiometry = parse_reaction(row.Equation, rev_arrow='=')
        except ValueError:
            continue
        else:
            for met, coeff in stoichiometry.items():
                met.name = chem_prop.loc[met.id]['name']
                try:
                    met.formula = Formula(chem_prop.loc[met.id].formula)
                except:
                    logger.debug('Cannot parse formula %s. Skipping formula' % chem_prop.loc[met.id].formula)
                    continue
                try:
                    met.charge = int(chem_prop.loc[met.id].charge)
                except (ValueError, TypeError):
                    logger.debug('Cannot parse charge %s. Skipping charge' % chem_prop.loc[met.id].charge)
                    pass
                rest = chem_prop.loc[met.id].to_dict()
                met.annotation = dict((key, rest[key]) for key in rest if key in ('mass', 'InChI', 'source'))
            mets = [met.id for met in stoichiometry.keys()]
            if len(mets) != len(set(mets)):
                continue
            reaction = Reaction(index)
            reaction.add_metabolites(stoichiometry)
            try:
                if len(reaction.check_mass_balance()) != 0:
                    continue
            except (AttributeError, ValueError) as e:
                logger.debug(str(e))
                continue
            if row.Balance:
                reaction.lower_bound = -1 * reaction.upper_bound
            reaction.name = row['Source']
            row = row.fillna("")
            rest = row.to_dict()
            reaction.annotation = dict((key, rest[key]) for key in rest if key in ('EC', 'Description'))
            reactions.append(reaction)

    model = Model('metanetx_universal_model_' + '_'.join(list_of_db_prefixes))
    model.add_reactions(reactions)
    # Add demands for all metabolites
    for metabolite in model.metabolites:
        model.add_boundary(metabolite, type='demand')
    return model


def load_metanetx_files():
    """"Update metanetx data."""
    BASE_URL = 'http://www.metanetx.org/cgi-bin/mnxget/mnxref/{}.tsv'
    for filename in ['chem_prop', 'chem_xref', 'reac_prop', 'reac_xref', 'comp_prop', 'comp_xref']:
        response = requests.get(BASE_URL.format(filename))
        filepath = '../data/metanetx/{}.tsv.gz'.format(filename)
        compress_by_lines(response, filepath)


def compress_by_lines(response, filepath):
    prev_line = next(response.iter_lines())
    with gzip.open(filepath, 'wb') as f:
        for line in response.iter_lines(decode_unicode=response.encoding):
            if line.startswith('#'):
                prev_line = line
                continue
            if prev_line:
                f.write(str.encode(prev_line + '\n'))
                prev_line = None
            f.write(str.encode(line + '\n'))


def add_to_all_mapping(dataframe, mapping):
    for other_id, mnx_id in dataframe[['XREF', 'MNX_ID']].values:
        cleaned_key = _apply_sanitize_rules(
            _apply_sanitize_rules(other_id, REVERSE_ID_SANITIZE_RULES_SIMPHENY),
            ID_SANITIZE_RULES_TAB_COMPLETION)
        mapping[cleaned_key] = mnx_id


def add_to_bigg_mapping(xref, bigg2mnx, mnx2bigg):
    bigg_selection = xref[['bigg' in blub for blub in xref.XREF]]
    sanitized_XREF = [
        _apply_sanitize_rules(_apply_sanitize_rules(id, REVERSE_ID_SANITIZE_RULES_SIMPHENY),
                              ID_SANITIZE_RULES_TAB_COMPLETION) for id in bigg_selection.XREF]
    bigg2mnx.update(dict(zip(sanitized_XREF, bigg_selection.MNX_ID)))
    mnx2bigg.update(dict(zip(bigg_selection.MNX_ID, sanitized_XREF)))


if __name__ == '__main__':

    import logging

    logging.basicConfig(level='INFO')

    if len(sys.argv) > 1 and sys.argv[1] == '--load':
        load_metanetx_files()

    # load metanetx data
    chem_xref = read_table('../data/metanetx/chem_xref.tsv.gz', compression='gzip')
    chem_xref.columns = [name.replace('#', '') for name in chem_xref.columns]
    reac_xref = read_table('../data/metanetx/reac_xref.tsv.gz', compression='gzip')
    reac_xref.columns = [name.replace('#', '') for name in reac_xref.columns]
    reac_prop = read_table('../data/metanetx/reac_prop.tsv.gz', compression='gzip', index_col=0)
    reac_prop.columns = [name.replace('#', '') for name in reac_prop.columns]
    chem_prop = read_table('../data/metanetx/chem_prop.tsv.gz', compression='gzip', index_col=0,
                           names=['name', 'formula', 'charge', 'mass', 'InChI', 'SMILES', 'source'])

    # replace NaN with None
    chem_prop = chem_prop.where((notnull(chem_prop)), "")

    REVERSE_ID_SANITIZE_RULES_SIMPHENY = [(value, key) for key, value in ID_SANITIZE_RULES_SIMPHENY]

    metanetx = dict()
    metanetx['all2mnx'] = dict()
    metanetx['bigg2mnx'] = dict()
    metanetx['mnx2bigg'] = dict()
    # Metabolites
    for xref in [chem_xref, reac_xref]:
        add_to_bigg_mapping(xref, metanetx['bigg2mnx'], metanetx['mnx2bigg'])
        add_to_all_mapping(xref, metanetx['all2mnx'])

    with gzip.open('../cameo/data/metanetx.json.gz', 'wt') as f:
        json.dump(metanetx, f)

    # generate universal reaction models
    db_combinations = [('bigg',), ('rhea',), ('bigg', 'rhea'), ('bigg', 'rhea', 'kegg'),
                       ('bigg', 'rhea', 'kegg', 'brenda')]
    largest_universal_model = None
    for db_combination in db_combinations:
        universal_model = construct_universal_model(db_combination, reac_xref, reac_prop, chem_prop)
        if largest_universal_model is None:
            largest_universal_model = universal_model
        else:
            if len(set(largest_universal_model.metabolites)) < len(set(universal_model.metabolites)):
                largest_universal_model = universal_model
        # The following is a hack; uncomment the following
        from cobra.io.json import _REQUIRED_REACTION_ATTRIBUTES

        _REQUIRED_REACTION_ATTRIBUTES.add('annotation')
        with open('../cameo/models/universal_models/{model_name}.json'.format(model_name=universal_model.id), 'w') as f:
            save_json_model(universal_model, f)

    metabolite_ids_in_largest_model = [m.id for m in largest_universal_model.metabolites]
    chem_prop_filtered = chem_prop.dropna(subset=['name'])
    chem_prop_filtered = chem_prop_filtered[chem_prop_filtered.index.isin(metabolite_ids_in_largest_model)]
    with gzip.open('../cameo/data/metanetx_chem_prop.json.gz', 'wt') as f:
        chem_prop_filtered.to_json(f, force_ascii=True)
