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

from __future__ import absolute_import, print_function

from functools import partial
from cobra.manipulation.delete import find_gene_knockout_reactions
from cobra.core import Reaction
from cameo.flux_analysis.simulation import fba
from cameo.exceptions import SolveError
from cameo.util import TimeMachine


def gene_knockout_growth(gene_id, model, threshold=10 ** -6, simulation_method=fba,
                         normalize=True, biomass=None, biomass_flux=None, *args, **kwargs):
    if biomass_flux is None:
        s = model.solve()
        biomass_flux = s.f
    if 'reference' not in kwargs:
        kwargs['reference'] = s.x_dict
    gene = model.genes.get_by_id(gene_id)
    knockouts = find_gene_knockout_reactions(model, [gene])
    tm = TimeMachine()

    for reaction in knockouts:
        tm(do=partial(setattr, reaction, 'lower_bound', 0),
           undo=partial(setattr, reaction, 'lower_bound', reaction.lower_bound))
        tm(do=partial(setattr, reaction, 'upper_bound', 0),
           undo=partial(setattr, reaction, 'upper_bound', reaction.upper_bound))

    try:
        s = simulation_method(model, *args, **kwargs)
        f = s.get_primal_by_id(biomass)
        if f >= threshold:
            if normalize:
                f = f / biomass_flux
        else:
            f = 0
    except SolveError:
        f = float('nan')
    finally:
        tm.reset()

    return f


def reaction_component_production(model, reaction):
    tm = TimeMachine()
    for metabolite in reaction.metabolites:
        test = Reaction("EX_%s_temp" % metabolite.id)
        test._metabolites[metabolite] = -1
        # hack frozen set from cobrapy to be able to add a reaction
        metabolite._reaction = set(metabolite._reaction)
        tm(do=partial(model.add_reactions, [test]), undo=partial(model.remove_reactions, [test]))
        tm(do=partial(setattr, model, 'objective', test.id), undo=partial(setattr, model, 'objective', model.objective))
        try:
            print(metabolite.id, "= ", model.solve().f)
        except SolveError:
            print(metabolite, " cannot be produced (reactions: %s)" % metabolite.reactions)
        finally:
            tm.reset()
