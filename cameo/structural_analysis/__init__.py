# Copyright 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
from itertools import combinations, chain
from ordered_set import OrderedSet
import sympy
import logging

logger = logging.getLogger("cameo")

mul = sympy.Mul._from_args
add = sympy.Add._from_args
RealNumber = sympy.RealNumber
NegativeOne = sympy.singleton.S.NegativeOne

N = "-"

def identify_currency_metabolites_by_pattern(model, top=20, min_combination=3, max_combination=5):
    """
    Identify currency metabolites and their patterns in reactions.

    Parameters
    ---------
    model : cobra.Model
    top : int
        Top max connected metabolites in the network, default: 20
    min_combination : int
        Size of the smallest pattern, default: 3
    max_combination : int
        Size of the longest pattern, default: 5

    Returns
    -------
    possible_motifs : dict
        The identified motifs and their rank
    currency_metabolites : list
        The ids of the most connected metabolites
    """
    metabolites_degree = dict([[m.id, len(m.reactions)] for m in model.metabolites])
    metabolite_ids = metabolites_degree.keys()
    metabolite_ids = sorted(metabolite_ids, key=metabolites_degree.get)
    metabolite_ids.reverse()
    currency_metabolites = metabolite_ids[:top]

    possible_motifs = {}

    for reaction in model.reactions:
        metabolites = [m.id for m in reaction.metabolites.keys() if m.id in currency_metabolites]
        r = xrange(min_combination, min(len(metabolites), max_combination)+1)
        motifs = list(set(chain(*[[frozenset(OrderedSet(c)) for c in combinations(metabolites, n)] for n in r])))
        if len(motifs) > 0:
            while len(motifs) > 0:
                motif = motifs.pop()
                print motifs
                print motif
                add = True
                for other_motif in motifs:
                    if motif.issubset(other_motif):
                        add = False
                    if other_motif.issubset(motif):
                        motifs.remove(other_motif)

                if add:
                    if motif in possible_motifs:
                        possible_motifs[motif] += 1
                    else:
                        possible_motifs[motif] = 1
    motifs = possible_motifs.keys()
    motifs = sorted(motifs, key=possible_motifs.get)
    motifs.reverse()

    for m in motifs[:10]:
        print m, ": ", possible_motifs[m]

    return possible_motifs, currency_metabolites


def shortest_elementary_flux_modes(model=None, k=5, M=100000):
    """
    Computing the shortest elementary flux modes in genome-scale metabolic networks.[1]

    Parameters
    ----------

    model : SolverBaseModel
        A COBRA model.
    k : int
        The number of iterations
    M : int
        A large number, default: 100000

    """
    efm_model = model.solver.interface.Model()

    z_map = dict()
    t_map = dict()
    constraints = list()
    exchanges = model.exchanges
    for reaction in model.reactions:
        if not reaction in exchanges:
            z = model.solver.interface.Variable("z_%s" % reaction.id, type='binary')
            z_map[reaction.id] = z
            t = model.solver.interface.Variable("t_%s" % reaction.id, lb=0, type='integer')
            t_map[reaction.id] = t

            expression_1 = add([t, mul([RealNumber(-M), z])])
            constraint_1 = model.solver.interface.Constraint(expression_1, name="c1_%s" % reaction.id, ub=0)
            constraints.append(constraint_1)

            expression_2 = add([z, mul([NegativeOne, t])])
            constraint_2 = model.solver.interface.Constraint(expression_2, name="c2_%s" % reaction.id, ub=0)
            constraints.append(constraint_2)


            if reaction.reversibility:
                z_rev = model.solver.interface.Variable("z_%s" % reaction._get_reverse_id(), type='binary')
                z_map[reaction._get_reverse_id()] = z_rev
                t_rev = model.solver.interface.Variable("t_%s" % reaction._get_reverse_id(), lb=0, type='integer')
                t_map[reaction._get_reverse_id()] = t_rev
                exp = add([z, z_rev])
                rev_constraint = model.solver.interface.Constraint(exp, ub=1)
                constraints.append(rev_constraint)

                expression_1_rev = add([t_rev, mul([RealNumber(-M), z_rev])])
                constraint_1_rev = model.solver.interface.Constraint(expression_1_rev,
                                                                     name="c1_%s_rev" % reaction.id,
                                                                     ub=0)
                constraints.append(constraint_1_rev)

                expression_2_rev = add([z_rev, mul([NegativeOne, t_rev])])
                constraint_2_rev = model.solver.interface.Constraint(expression_2_rev,
                                                                     name="c2_%s_rev" % reaction.id,
                                                                     ub=0)
                constraints.append(constraint_2_rev)

    for met in model.metabolites:
        steady_state_expression = sum([mul([RealNumber(r.metabolites[met]), t_map[r.id]]) for r in met.reactions if not r in exchanges])
        steady_state_constraint = model.solver.interface.Constraint(steady_state_expression,
                                                                    name="ss_%s" % met.id,
                                                                    lb=0,
                                                                    ub=0)
        constraints.append(steady_state_constraint)

    trivial_solutions_expression = sum(z_map.values())
    trivial_solutions_constraint = model.solver.interface.Constraint(trivial_solutions_expression, lb=1)
    constraints.append(trivial_solutions_constraint)

    objective = model.solver.interface.Objective(add(z_map.values()), direction='min')
    efm_model.add(z_map.values())
    efm_model.add(t_map.values())
    efm_model.add(constraints)
    efm_model.objective = objective
    #iteration 1
    status = efm_model.optimize()
    logger.debug("Iteration 1: %s" % status)
    [logger.debug("%s: %f" % (z, z.primal)) for z in z_map.values() if z.primal > 0]
    for i in xrange(k):
        previous_solution_part = list()
        sum_part = 0
        for reaction in model.reactions:
            if not reaction in exchanges:
                z = z_map[reaction.id]
                previous_solution_part.append(mul([RealNumber(z.primal), z]))
                sum_part += z.primal
                if reaction.reversibility:
                    z = z_map[reaction._get_reverse_id()]
                    previous_solution_part.append(mul([RealNumber(z.primal), z]))
                    sum_part += z.primal
        iteration_expression = add(previous_solution_part)
        iteration_constraint = model.solver.interface.Constraint(iteration_expression, ub=sum_part-1)
        efm_model.add(iteration_constraint)

        #interation n
        status = efm_model.optimize()
        logger.debug("Iteration %i: %s" % (i+2, status))
        [logger.debug("%s: %f" % (z, z.primal)) for z in z_map.values() if z.primal > 0]
        efm_model.remove(iteration_constraint)