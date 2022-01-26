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
import csv
import re

import pandas as pd

from collections import namedtuple
from io import StringIO
from functools import partial, reduce
from pandas import DataFrame, read_csv
from typing import Dict, List
from cobra import Metabolite, Model, Reaction

from cameo.io import load_model


# TODO: Define the product
# TODO: Visualization with ESCHER
class MetaboliteError(Exception):
    """Raised when something went wrong when metabolite string is malformated."""
    def __init__(self, msg):
        self._msg = msg

    def __str__(self):
        return 'Metabolite is malformated: {}'.format(self._msg)


class EquationError(Exception):
    """Raised when something went wrong when equation string is malformated."""
    def __init__(self, msg, equation=''):
        self._msg = msg
        self._equation = equation

    def __str__(self):
        return 'Equation is malformated: {}\n{}'.format(self._msg, self._equation)


class PathwayFileError(Exception):
    """Raised when something went wrong when pathway file is malformated."""
    def __init__(self, msg):
        self._msg = msg

    def __str__(self):
        return 'Pathway file is malformated: {}'.format(self._msg)


class Equation(object):
    """
    Representation of an equation (a set of metabolites)

    Attributes
    ----------

    stoichiometry: dict of Metabolite, coefficient

    """

    ReArrow = r'\s*[-=]+>\s*'
    Token = namedtuple('Token', ['kind', 'value', 'group', 'target'])

    def __init__(
        self,
        stochiometry: Dict[Metabolite, float],
        *args,
        **kwargs
    ) -> None:
        """Initialize an Equation.

        Parameters
        ----------
        stochiometry: dict
            A dictionary of:
            - key: cobra.Metabolite
            - value: int

        Returns
        -------
        None
        """
        super(Equation, self).__init__(*args, **kwargs)
        self.stochiometry = stochiometry

    @staticmethod
    def contains_arrow(
        segment: str
    ) -> bool:
        """Check if a string contains an arrow defined to represent an Equation.

        Parameters
        ----------
        segment: str
            A string to check

        Returns
        -------
        bool
        """
        if re.search(Equation.ReArrow, segment):
            return True
        return False

    @staticmethod
    def tokenize(
        segment: str
    ) -> 'Equation.Token':
        """Parse a String to return differents tokens describing Metabolite & coefficient \
        in the same order of the String.

        Parameters
        ----------
        segment: str
            A string to parse

        Returns
        -------
        TokenEquation
            namedtuple: kind(COEFF|MET), value(int|cobra.Metabolite), group(product|reactant), target(bool)
        """
        token_specification = [
            ('ARROW', Equation.ReArrow),
            ('COEFFSEP', r'[\*]{1}'),
            ('MET', r'(#{1}[\d\w\+\-\s]+#?){3,4}'),
            ('COEFF', r'\-?(\d+)'),
            ('AND', r'[\+]{1}'),
            ('SKIP', r'[ \t]+'),
            ('END', r'\Z')
        ]
        tok_regex = '|'.join('(?P<%s>%s)' % pair for pair in token_specification)

        is_product, is_target = False, False
        next_met = False
        for mo in re.finditer(tok_regex, segment):
            kind = mo.lastgroup
            value = mo.group()
            if kind == 'ARROW':
                is_product = True
            elif kind == 'AND':
                continue
            elif kind == 'COEFF' and not next_met:
                value = int(value)
            elif kind == 'COEFFSEP':
                next_met = True
                continue
            elif kind == 'MET' and next_met:
                values = value.split('#')
                values = [x for x in values if x != '']
                if len(values) < 3 or len(values) > 4:
                    raise MetaboliteError(
                        'When a metabolite is provided, it needs to be \
                        formated as #<name>#<id>#<compartement_id>#[<target>#]'
                    )
                elif len(values) == 4:
                    m = re.search(r'target', values[3], flags=re.I)
                    if m is None:
                        raise MetaboliteError('When four # is filled, a target keyword is attempted')
                    else:
                        is_target = True
                value = Metabolite(
                    id=values[1],
                    name=values[0],
                    compartment=values[2]
                )
                next_met = False
            elif kind == 'SKIP':
                continue
            elif kind == 'END':
                continue
            group = 'product' if is_product else 'reactant'
            token = Equation.Token(
                kind,
                value,
                group,
                is_target
            )
            is_target = False

            yield token

    @staticmethod
    def from_string(
        equation: str
    ) -> 'Equation':
        """Construct an Equation object from a string.

        Parameters
        ----------
        equation: str
            An equation to parse.
            Must have a specific syntax:
              - contains an arrow like '->' to separate product/reactant
              - describe product and reactant like:
                '<coefficient:int> * <metabolite:str>'
              - product and reactant are combined with '+' operator
              - a metabolite is described like:
                '#<id:str>#<name:str>#<compartment_id:str#'
                Optionaly, you could specify if it's a target by adding to the current metabolite:
                'target#'

        Returns
        -------
        Equation
            A novel object.
        """
        stoichiometry = {}

        nb_coef, nb_met = 0, 0
        nb_target, nb_arrow = 0, 0
        coef = None
        for token in Equation.tokenize(equation):
            if token.kind == 'ARROW':
                nb_arrow += 1
            elif token.kind == 'COEFF':
                nb_coef += 1
                if token.group == 'product':
                    coef = abs(float(token.value))
                else:
                    coef = -abs(float(token.value))
            elif token.kind == 'MET':
                nb_met += 1
                if coef is None:
                    raise MetaboliteError('A coefficient must be set with each metabolite')
                stoichiometry[token.value] = coef
                coef = None
                if token.target:
                    nb_target += 1

        if nb_arrow != 1:
            raise EquationError('An arrow must be indicated in equation string or must be surrounded by spaces')
        elif nb_coef != nb_met:
            raise EquationError('The number of coefficient must be the same as the number of metabolites')
        elif nb_target > 1:
            raise EquationError('A maximum of one metabolite target must be indicated')

        return Equation(stoichiometry)

    def to_string(self) -> str:
        """Represent an Equation to a string.

        Returns
        -------
        String
            A representation of an Equation.
        """

        products = {m: v for m, v in self.stochiometry.items() if v > 0}
        reactants = {m: v for m, v in self.stochiometry.items() if v < 0}

        products = " + ".join(['%s * #%s#%s#%s#' % (int(v), m.name, m.id, m.compartment) for m, v in products.items()])
        reactants = " + ".join(['%s * #%s#%s#%s#' % (-int(v), m.name, m.id, m.compartment) for m, v in reactants.items()])
        return reactants + " <=> " + products

    def to_reaction(
        self,
        identifier: str,
        name: str,
        lower: float = 0.0,
        upper: float = 1000.0,
        comments: str = ''
    ) -> Reaction:
        """Construct an Equation object from a string.

        Parameters
        ----------
        identifier: str
            An id for the reaction
        name: str
            A name for the reaction
        lower: float
            A lower bound for the reaction (Default: 0.0)
        upper: float
            An upper bound for the reaction (Default: 1000.0)
        comments: str
            A comment to associate for the reaction (Default: '')

        Returns
        -------
        Reaction
            A novel object.
        """
        reaction = Reaction(identifier)
        reaction.id = identifier
        reaction.name = name
        reaction.add_metabolites(self.stochiometry)
        reaction.upper_bound = float(upper)
        reaction.lower_bound = float(lower)
        reaction.notes["pathway_note"] = comments
        return reaction


class Pathway(object):
    """Representation of a pathway (a set of reactions)

    Attributes
    ----------
    reactions: list of cobra.Reaction
        The list of reactions in the pathway.
    """

    FileHeader = ['id', 'name', 'equation', 'lower_bound', 'upper_bound', 'comments']

    def __init__(
        self,
        reactions: List[Reaction],
        *args,
        **kwargs
    ) -> None:
        """Initialize a Pathway.

        Attributes
        ----------
        reactions: list of cobra.Reaction
            The list of reactions in the pathway.

        Returns
        -------
        None
        """
        super(Pathway, self).__init__(*args, **kwargs)
        self.reactions = reactions

    @property
    def data_frame(self) -> DataFrame:
        """Represents the reactions with a dataframe.

        Returns
        -------
        pd.DataFrame
            A dataframe with :
              - row: a reaction
              - col: equation, lower bound, upper bound
        """

        return DataFrame([[r.build_reaction_string(True), r.lower_bound, r.upper_bound] for r in self.reactions],
                         columns=["equation", "lower_bound", "upper_bound"], index=[r.id for r in self.reactions])

    @classmethod
    def from_model(
        cls,
        path_or_handle,
        *args,
        **kwargs
    ) -> 'Pathway':
        """Read a pathway from a model.
        Call cameo.io.load_model and use the same arguments.

        Returns
        -------
        Pathway
            A novel object

        See Also
        --------
        cameo.io.load_model
        Pathway.from_file
        """

        model = load_model(path_or_handle, *args, **kwargs)
        return Pathway(model.reactions)

    @classmethod
    def from_file(
        cls,
        file_path: str,
        sep: str = 'auto'
    ) -> 'Pathway':
        """Read a pathway from a file.

        The file format is:
        id<sep>name<sep>equation<sep>lower_bound<sep>upper_bound<sep>comments\n
        This header must be provided in the file.

        The equation is defined by:
        coefficient * #substrate_name#substrate_id#substrate_compartment# + ... => \
            coefficient * #product_name#product_id#product_compartment#

        Parameters
        ----------
        file_path: str
            The path to the file containing the pathway
        sep: str
            The separator between elements in the file (default: "auto")

        Returns
        -------
        Pathway
            A novel object
        """
        # Sniff if file is well formated.
        if sep == 'auto':
            with open(file_path) as fid:
                dialect = csv.Sniffer().sniff(fid.readline())
            sep = dialect.delimiter

        # Check header.
        lines = []
        with open(file_path) as fid:
            csv_reader = csv.reader(fid, delimiter=sep)
            lines = [x for x in csv_reader if not x[0].startswith('#')]

        if len(lines) < 1:
            raise PathwayFileError('You must provide at lear one entry except the header.')
        if len(lines[0]) != 6:
            raise PathwayFileError('Too much fields are in the file, you must respect the format attempted.')

        header = lines[0]
        if Equation.contains_arrow(sep.join(header)):
            raise PathwayFileError('Pathway file must contain an header.')

        # Buffering.
        lines = [sep.join(x) for x in lines]
        lines = '\n'.join(lines)

        # Load.
        df = read_csv(StringIO(lines), sep=sep, na_values='', skiprows=1, names=Pathway.FileHeader)

        # Build reactions
        def _build_reaction(x):
            eq = Equation.from_string(x['equation'])
            return eq.to_reaction(x['id'], x['name'], x['lower_bound'], x['upper_bound'], x['comments'])
        reactions = [x for x in df.apply(_build_reaction, axis=1)]

        return cls(reactions)

    def to_file(
        self,
        file_path: str,
        sep: str = '\t'
    ) -> None:
        """Writes the pathway to a file.

        Parameters
        ----------
        file_path: str
            The path to the file where the pathway will be written
        sep: str
            The separator between elements in the file (default: "\t")

        Returns
        -------
        None

        See Also
        --------
        Pathway.from_file
        """
        with open(file_path, 'w') as output_file:
            output_file.write(sep.join(Pathway.FileHeader) + '\n')
            for reaction in self.reactions:
                equation = Equation(reaction.metabolites)
                comments = reaction.notes.get('pathway_note', '')
                if pd.isna(comments):
                    comments = ''
                output_file.write(sep.join(map(str, [
                    reaction.id,
                    equation.to_string(),
                    reaction.lower_bound,
                    reaction.upper_bound,
                    reaction.name,
                    comments
                ])) + '\n')

    def plug_model(
        self,
        model: Model
    ) -> Model:
        """Plug the pathway to a model.

        Metabolites are matched in the model by id.\
        For metabolites with no ID in the model, an exchange reaction is added to the model

        Parameters
        ----------
        model: cobra.Model
            The model to plug in the pathway

        Returns
        -------
        cobra.Model
            A model with pathway reactions added.
        """
        model.add_reactions(self.reactions)
        metabolites = set(reduce(lambda x, y: x + y, [list(r.metabolites.keys()) for r in self.reactions], []))
        exchanges = [model.add_boundary(m) for m in metabolites if len(m.reactions) == 1]
        for exchange in exchanges:
            exchange.lower_bound = 0.0
