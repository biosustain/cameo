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

from functools import partial

import six

from cameo import Reaction


def swap_cofactors(reaction, model, swap_pairs, inplace=True, time_machine=None):
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
        else:
            model.add_reactions([new_reaction])

        return new_reaction
