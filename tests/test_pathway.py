import filecmp
import os
import pytest
import tempfile

import pandas as pd

from cobra import Metabolite, Reaction
from cameo.core.pathway import Equation, Pathway


@pytest.fixture(scope="module")
def but_csv(data_directory):
    return os.path.join(data_directory, '1-butanol.csv')


@pytest.fixture(scope="module")
def but_xml(data_directory):
    return os.path.join(data_directory, '1-butanol.xml')


@pytest.fixture(scope='module')
def but_df(but_csv):
    return pd.read_csv(but_csv)


@pytest.fixture(scope='module')
def but_met():
    mets = {}
    mets['b2coa_c'] = Metabolite(
        id='b2coa_c',
        name='Crotonoyl-CoA',
        compartment='c'
    )
    mets['nadh_c'] = Metabolite(
        id='nadh_c',
        name='Nicotinamide adenine dinucleotide - reduced',
        compartment='c'
    )
    mets['h_c'] = Metabolite(
        id='h_c',
        name='H+',
        compartment='c'
    )
    mets['btcoa_c'] = Metabolite(
        id='btcoa_c',
        name='Gamma-butyrobetainyl-CoA',
        compartment='c'
    )
    mets['nad_c'] = Metabolite(
        id='nad_c',
        name='Deamino-NAD+',
        compartment='c'
    )
    return mets



@pytest.fixture(scope='module')
def but_st(but_met):
    st = {}
    st[but_met['b2coa_c']] = -1
    st[but_met['nadh_c']] = -1
    st[but_met['h_c']] = -1
    st[but_met['btcoa_c']] = 1
    st[but_met['nad_c']] = 1
    return st


class TestEquation:
    def test_contains_arrow(self):
        segments = []
        segments.append(('aa=>aa', True))
        segments.append(('aa->a', True))
        segments.append(('aa->', True))
        segments.append(('aa ->a', True))
        segments.append(('aa -> ', True))
        segments.append(('-> a', True))
        segments.append(('a>a', False))
        segments.append(('a', False))
        segments.append(('-a', False))
        segments.append(('', False))
        for tes in segments:
            assert Equation.contains_arrow(tes[0]) == tes[1]

    def test_tokenize(self, but_df, but_met):
        equations = but_df['equation'].to_list()
        tokens = [x for x in Equation.tokenize(equations[0])]

        # Theorical list
        ths = []
        ths.append(Equation.Token(kind='COEFF', value=1, group='reactant', target=False))
        ths.append(Equation.Token(kind='MET', value=but_met['b2coa_c'], group='reactant', target=False))
        ths.append(Equation.Token(kind='COEFF', value=1, group='reactant', target=False))
        ths.append(Equation.Token(kind='MET', value=but_met['nadh_c'], group='reactant', target=False))
        ths.append(Equation.Token(kind='COEFF', value=1, group='reactant', target=False))
        ths.append(Equation.Token(kind='MET', value=but_met['h_c'], group='reactant', target=False))
        ths.append(Equation.Token(kind='ARROW', value='=> ', group='product', target=False))
        ths.append(Equation.Token(kind='COEFF', value=1, group='product', target=False))
        ths.append(Equation.Token(kind='MET', value=but_met['btcoa_c'], group='product', target=False))
        ths.append(Equation.Token(kind='COEFF', value=1, group='product', target=False))
        ths.append(Equation.Token(kind='MET', value=but_met['nad_c'], group='product', target=False))

        assert len(tokens) == len(ths)
        for token, th in zip(tokens, ths):
            # cobra.Metabolite does not implement yet __eq__()
            if isinstance(token.value, Metabolite) and isinstance(th.value, Metabolite):
                assert token.value.id == th.value.id
                assert token.value.name == th.value.name
                assert token.value.compartment == th.value.compartment
            else:
                assert token == th

    def test_from_string(self, but_df, but_met):
        equations = but_df['equation'].to_list()
        equation = Equation.from_string(equations[0])
        stochiometry = {}
        for k, v in equation.stochiometry.items():
            stochiometry[k.id] = v

        # Theorical value
        ths = {}
        ths['b2coa_c'] = -1.0
        ths['nadh_c'] = -1.0
        ths['h_c'] = -1.0
        ths['btcoa_c'] = 1.0
        ths['nad_c'] = 1.0

        assert stochiometry == ths

    def test_to_string(self, but_df, but_st):
        equation = Equation(but_st)

        # Theorical value
        equations = but_df['equation'].to_list()

        assert equation.to_string() == equations[0]

    def test_to_reaction(self, but_df, but_st):
        sequation = but_df.loc[0]
        equation = Equation.from_string(sequation['equation'])

        reaction = equation.to_reaction(
            identifier=sequation['id'],
            name=sequation['name'],
            lower=sequation['lower_bound'],
            upper=sequation['upper_bound'],
            comments=sequation['comments']
        )

        # Theorical value
        th = Reaction(id='ButCoaDeh', name='Butyryl-CoA dehydrogenase', lower_bound=0, upper_bound=1000)
        th.add_metabolites(but_st)
        th.notes["pathway_note"] = 'metanetx.reaction:undefined'

        # cobra.Reaction does not implement yet __eq__()
        assert reaction.id == th.id
        assert reaction.name == th.name
        assert reaction.bounds == th.bounds
        assert reaction.notes == th.notes


class TestPathway:
    def test_from_model(self, but_xml, but_csv):
        pathway_model = Pathway.from_model(but_xml)
        pathway_file = Pathway.from_file(but_csv)

        df_model = pathway_model.data_frame.sort_index()
        df_file = pathway_file.data_frame.sort_index()
        assert df_model.equals(df_file)

    def test_from_file(self, but_csv):
        pathway = Pathway.from_file(but_csv)

        assert len(pathway.reactions) == 6

    def test_to_file(self, but_csv, but_df):
        pathway = Pathway.from_file(but_csv)
        tmp = ''
        with tempfile.NamedTemporaryFile(delete=False) as fod:
            tmp = fod.name
            pathway.to_file(tmp, sep=",")
            lines = fod.read().splitlines()

        assert len(lines) == but_df.shape[0] + 1

        os.remove(fod.name)

    def test_plug_model(self, salmonella, but_csv):
        pathway = Pathway.from_file(but_csv)

        assert len([x for x in salmonella.reactions if x.id == 'ButOlTrE']) == 0
        pathway.plug_model(salmonella)
        assert len([x for x in salmonella.reactions if x.id == 'ButOlTrE']) == 1
