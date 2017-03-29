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

from cobra import Reaction


def increase_flux(reaction, ref_value, value):
    """
    lb                           0                           ub
    |--------------------------- ' ---------------------------|
                <- - -|----------'
                                 '----------|- - - ->

    Parameters
    ----------
    reaction: cobra.Reaction
        The reaction to over-express.
    ref_value: float
        The flux value to come from.
    value: float
        The flux value to achieve.

    """
    if abs(value) < abs(ref_value):
        raise ValueError("'value' is lower than 'ref_value', this is increase_flux (%f < %f)" % (value, ref_value))

    if value > 0:
        reaction.lower_bound = value
    elif value < 0:
        reaction.upper_bound = value
    else:
        reaction.knock_out()


def decrease_flux(reaction, ref_value, value):
    """
    lb                           0                           ub
    |--------------------------- ' ---------------------------|
                 |- - >----------'
                                 '----------<- - - -|

    Parameters
    ----------
    reaction: cobra.Reaction
        The reaction to down_regulate.
    ref_value: float
        The flux value to come from.
    value: float
        The flux value to achieve.

    """
    if abs(value) > abs(ref_value):
        raise ValueError("'value' is higher than 'ref_value', this is decrease_flux (%f < %f)" % (value, ref_value))

    if value > 0:
        reaction.upper_bound = value
    elif value < 0:
        reaction.lower_bound = value
    else:
        reaction.knock_out()


def reverse_flux(reaction, ref_value, value):
    """

    Forces a reaction to have a minimum flux level in the opposite direction of a reference state.

    lb                           0                           ub
    |--------------------------- ' ---------------------------|
                      <----------'- - - - - - - ->

    Parameters
    ----------
    reaction: cobra.Reaction
        The reaction that will be inverted.
    ref_value: float
        The flux value to come from.
    value: float
        The flux value to achieve.

    """
    if (value >= 0) == (ref_value >= 0):
        raise ValueError("'value' and 'ref_value' cannot have the same sign (%.5f, %.5f)" % (value, ref_value))

    if value > 0:
        reaction.upper_bound = value
    elif value < 0:
        reaction.lower_bound = value
    else:
        reaction.knock_out()


def swap_cofactors(reaction, model, swap_pairs, inplace=True):
    """
    Swaps the cofactors of a reaction. For speed, it can be done inplace which just changes the coefficients.
    If not done inplace, it will create a new Reaction, add it to the model, and knockout the original reaction.

    Parameters
    ----------
    reaction: cobra.Reaction
        The reaction to swap.
    model: cameo.cobra.Model
        A constraint-based model.
    swap_pairs: tuple
        A tuple of (cofactors, equivalent_cofactors)
    inplace: bool
        If replace is done inplace, it changes the coefficients in the matrix. Otherwise, it creates a new reaction
        with the other cofactors and adds it to the model.

    Returns
    -------
        Reaction
            A reaction with swapped cofactors (the same if inplace).
    """
    if all(reaction.metabolites.get(met, False) for met in swap_pairs[0]):
        new_coefficients = {met: -reaction.metabolites[met] for met in swap_pairs[0]}
        new_coefficients.update({new_met: reaction.metabolites[met] for met, new_met in zip(*swap_pairs)})
    elif all(reaction.metabolites.get(met, False) for met in swap_pairs[1]):
        new_coefficients = {met: -reaction.metabolites[met] for met in swap_pairs[1]}
        new_coefficients.update({new_met: reaction.metabolites[met] for new_met, met in zip(*swap_pairs)})
    else:
        raise ValueError("%s: Invalid swap pairs %s (%s)" % (reaction.id, str(swap_pairs), reaction.reaction))

    def _inplace(rxn, stoichiometry):
        rxn.add_metabolites(stoichiometry, combine=True)

    def _replace(rxn, stoichiometry):
        new_reaction = Reaction(id="%s_swap" % rxn.id, name=rxn.name,
                                lower_bound=rxn.lower_bound, upper_bound=rxn.upper_bound)
        new_reaction.stoichiometry = rxn
        new_reaction.add_metabolites(stoichiometry)
        return new_reaction

    if inplace:
        _inplace(reaction, new_coefficients)
        return reaction
    else:
        new_reaction = _replace(reaction, new_coefficients)
        model.add_reactions([new_reaction])
        reaction.knock_out()
        return new_reaction
