# Copyright 2016 The Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
Manage manipulations such as swapping reaction cofactors, over-express or down-regulate genes and reactions.

"""

from functools import partial

import six

from cameo.core import Reaction, Gene


def over_express(gene_or_reaction, ref_value, value, time_machine=None):
    if isinstance(gene_or_reaction, Reaction):
        over_express_reaction(gene_or_reaction, ref_value, value, time_machine=time_machine)
    elif isinstance(gene_or_reaction, Gene):
        over_express_gene(gene_or_reaction, ref_value, value, time_machine=time_machine)

    else:
        raise ValueError("'gene_or_reaction' is neither instance of Gene or Reaction")


def down_regulate(gene_or_reaction, ref_value, value, time_machine=None):
    if isinstance(gene_or_reaction, Reaction):
        down_regulate_reaction(gene_or_reaction, ref_value, value, time_machine=time_machine)
    elif isinstance(gene_or_reaction, Gene):
        down_regulate_gene(gene_or_reaction, ref_value, value, time_machine=time_machine)

    else:
        raise ValueError("'gene_or_reaction' is neither instance of Gene or Reaction")


def over_express_reaction(reaction, ref_value, value, time_machine=None):
    """
    lb                           0                           ub
    |--------------------------- ' ---------------------------|
                <- - -|----------'
                                 '----------|- - - ->

    reaction: cameo.Reaction
        The reaction to over-express.
    ref_value: float
        The flux value to come from.
    value: float
        The flux value to achieve.
    time_machine: TimeMachine
        The action stack manager.

    """

    if abs(value) < abs(ref_value):
        raise ValueError("'value' is lower than 'ref_value', this is down_regulation (%f < %f)" % (value, ref_value))

    if value > 0:
        reaction.change_bounds(lb=value, time_machine=time_machine)
    elif value < 0:
        reaction.change_bounds(ub=value, time_machine=time_machine)
    else:
        reaction.knock_out(time_machine=time_machine)


def down_regulate_reaction(reaction, ref_value, value, time_machine=None):
    """
    lb                           0                           ub
    |--------------------------- ' ---------------------------|
                 |- - >----------'
                                 '----------<- - - -|

    reaction: cameo.Reaction
        The reaction to down_regulate.
    ref_value: float
        The flux value to come from.
    value: float
        The flux value to achieve.
    time_machine: TimeMachine
        The action stack manager.

    """
    if abs(value) > abs(ref_value):
        raise ValueError("'value' is higher than 'ref_value', this is down_regulation (%f < %f)" % (value, ref_value))

    if value > 0:
        reaction.change_bounds(ub=value, time_machine=time_machine)
    elif value < 0:
        reaction.change_bounds(lb=value, time_machine=time_machine)
    else:
        reaction.knock_out(time_machine=time_machine)


def over_express_gene(gene, ref_value, value, time_machine=None):
    """
    0                       expression
    |---------------------------|
    |---------|- - ->

    gene: cameo.Gene
        The gene to over-express.
    ref_value: float
        The expression value to come from.
    value: float
        The expression value to achieve.
    time_machine: TimeMachine
        The action stack manager.

    """

    if value < 0:
        raise ValueError("'value' expression value is always positive")

    if value < ref_value:
        raise ValueError("'value' is lower than 'ref_value', this is down_regulation (%f < %f)" % (value, ref_value))

    raise NotImplementedError


def down_regulate_gene(gene, ref_value, value, time_machine=None):
    """
    0                       expression
    |---------------------------|
    |---------<- - -|

    gene: cameo.Gene
        The reaction to down_regulate.
    ref_value: float
        The flux value to come from.
    value: float
        The flux value to achieve.
    time_machine: TimeMachine
        The action stack manager.

    """

    if value < 0:
        raise ValueError("'value' expression value is always positive")

    if value > ref_value:
        raise ValueError("'value' is higher than 'ref_value', this is down_regulation (%f < %f)" % (value, ref_value))

    raise NotImplementedError


def swap_cofactors(reaction, model, swap_pairs, inplace=True, time_machine=None):
    """
    Swaps the cofactors of a reaction. For speed, it can be done inplace which just changes the coefficients.
    If not done inplace, it will create a new Reaction, add it to the model, and knockout the original reaction.

    reaction: cameo.Reaction
        The reaction to swap.
    model: cameo.SolverBasedModel
        A constraint-based model.
    swap_pairs: tuple
        A tuple of (cofactors, equivalent_cofactors)
    inplace: bool
        If replace is done inplace, it changes the coefficients in the matrix. Otherwise, it creates a new reaction
        with the other cofactors and adds it to the model.
    time_machine: : TimeMachine
        The action stack manager.

    Returns
    -------
        Reaction
            A reaction with swapped cofactors (the same if inplace).
    """

    if all(reaction.metabolites.get(met, False) for met in swap_pairs[0]):
        new_coefficients = {met: -reaction.metabolites[met] for met in swap_pairs[0]}
        new_coefficients.update({new_met: reaction.metabolites[met] for met, new_met in zip(*swap_pairs)})
        revert_coefficients = {met: -coeff for met, coeff in six.iteritems(new_coefficients)}
    elif all(reaction.metabolites.get(met, False) for met in swap_pairs[1]):
        new_coefficients = {met: -reaction.metabolites[met] for met in swap_pairs[1]}
        new_coefficients.update({new_met: reaction.metabolites[met] for new_met, met in zip(*swap_pairs)})
        revert_coefficients = {met: -coeff for met, coeff in six.iteritems(new_coefficients)}
    else:
        raise ValueError("%s: Invalid swap pairs %s (%s)" % (reaction.id, str(swap_pairs), reaction.reaction))

    def _inplace(reaction, stoichiometry):
        reaction.add_metabolites(stoichiometry, combine=True)

    def _replace(reaction, stoichiometry):
        new_reaction = Reaction(id="%s_swap" % reaction.id, name=reaction.name,
                                lower_bound=reaction.lower_bound, upper_bound=reaction.upper_bound)
        new_reaction.stoichiometry = reaction
        new_reaction.add_metabolites(stoichiometry)
        return new_reaction

    if inplace:
        if time_machine:
            time_machine(do=partial(_inplace, reaction, new_coefficients),
                         undo=partial(_inplace, reaction, revert_coefficients))
        else:
            _inplace(reaction, new_coefficients)

        return reaction
    else:
        new_reaction = _replace(reaction, new_coefficients)
        if time_machine:
            time_machine(do=partial(model.add_reactions, [new_reaction]),
                         undo=partial(model.remove_reactions, [new_reaction], delete=False))
            time_machine(do=partial(model.remove_reactions, [reaction], delete=False),
                         undo=partial(model.add_reactions, [new_reaction]))
            reaction.knock_out(time_machine)
        else:
            model.add_reactions([new_reaction])
            reaction.knock_out()

        return new_reaction
