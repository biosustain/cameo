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
from functools import partial
from itertools import combinations, chain
from ordered_set import OrderedSet
import sympy
import logging
from cameo.util import TimeMachine

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


class EFMModel(object):
    """
    Implementation of optimization model used in Figueiredo et al 2009[1] and Kamp et al 2014 [2].

    .. [1] Luis F. de Figueiredo, Adam Podhorski, Angel Rubio, Christoph Kaleta, John E. Beasley,
       Stefan Schuster and Francisco J. Planes "Computing the shortest elementary flux modes in
       genome-scale metabolic networks" Bioinformatics, vol 25, issue 23, pp. 3158-3165, 2009.
    .. [2] Axel von Kamp, Steffen Klamt "Enumeration of Smallest Intervention Strategies in Genome-Scale
       Metabolic Networks" PLOS Computational Biology, vol 10, issue 01, pp. e1003378, 2014.
    """
    def __init__(self, model, M=100000, matrix=None, include_exchanges=False, flux_type="integer"):
        self.z_map = dict()
        self.t_map = dict()
        self.M = M
        if matrix is None:
            matrix = model.S
        self.matrix = matrix
        self.model = model.solver.interface.Model()
        self._populate_model(model, include_exchanges, flux_type)

    def add(self, obj):
        self.model.add(obj)

    def remove(self, obj):
        self.model.remove(obj)

    @property
    def z_vars(self):
        return self.z_map.values()

    def z_var(self, reaction_id):
        return self.z_map[reaction_id]

    @property
    def t_vars(self):
        return self.t_map.values()

    def t_var(self, reaction_id):
        return self.t_map[reaction_id]

    def _populate_model(self, model, include_exchanges, flux_type):
        constraints = list()
        exchanges = model.exchanges
        for reaction in model.reactions:
            if not reaction in exchanges or include_exchanges:
                z, t = self.add_reaction(reaction.id, model, constraints, flux_type)
                if reaction.reversibility:
                    z_rev, t_rev = self.add_reaction(reaction._get_reverse_id(), model, constraints, flux_type)
                    exp = add([z, z_rev])
                    rev_constraint = model.solver.interface.Constraint(exp, ub=1, name="rev_%s" % reaction.id)
                    constraints.append(rev_constraint)

        for met in model.metabolites:
            aux = list()
            for r in met.reactions:
                if not r in exchanges or include_exchanges:
                    i = model.metabolites.index(met.id)
                    j = model.reactions.index(r.id)
                    coeff = self.matrix[i, j]

                    aux.append(mul([RealNumber(coeff), self.t_map[r.id]]))
                    if r.reversibility:
                        aux.append(mul([RealNumber(-coeff), self.t_map[r._get_reverse_id()]]))

            steady_state_expression = add(aux)
            steady_state_constraint = model.solver.interface.Constraint(steady_state_expression,
                                                                        name="ss_%s" % met.id,
                                                                        lb=0,
                                                                        ub=0)
            constraints.append(steady_state_constraint)

        trivial_solutions_expression = add(self.z_vars)
        trivial_solutions_constraint = model.solver.interface.Constraint(trivial_solutions_expression, lb=1)
        constraints.append(trivial_solutions_constraint)

        objective = model.solver.interface.Objective(add(self.z_map.values()), direction='min')
        self.model.add(self.z_vars)
        self.model.add(self.t_vars)
        self.model.add(constraints)
        self.model.objective = objective

    def add_reaction(self, r_id, model, constraints, flux_type):
        z = model.solver.interface.Variable("z_%s" % r_id, type='binary')
        self.z_map[r_id] = z
        t = model.solver.interface.Variable("t_%s" % r_id, lb=0, ub=self.M, type=flux_type)
        self.t_map[r_id] = t

        expression_1 = add([mul([RealNumber(-self.M), z]), t])
        constraint_1 = model.solver.interface.Constraint(expression_1, name="c1_%s" % r_id, ub=0)
        constraints.append(constraint_1)

        expression_2 = add([z, mul([NegativeOne, t])])
        constraint_2 = model.solver.interface.Constraint(expression_2, name="c2_%s" % r_id, ub=0)
        constraints.append(constraint_2)

        return z, t

    def optimize(self):
        return self.model.optimize()


#FIXME: the matrix must be converted into integer values to work!!!!
def shortest_elementary_flux_modes(model=None, k=5, M=100000, efm_model=None, matrix=None):
    """
    Calculates the shortest elementary flux modes using MILP.[1]
    This requires that the stoichiometric matrix is converted into a integer matrix.

    If you the efm model is provided, it will not be rebuild.
    If the matrix is provided, it will not be recomputed.

    Parameters
    ----------

    model : SolverBaseModel
        A COBRA model.
    k : int
        The number of iterations
    M : int
        A large number, default: 100000
    efm_model : EFMModel
        A previous built model, optional
    matrix : numpy.Array
        A numpy.Array like object that represents the integer version of the stoichiometric matrix, optional

    .. [1] Luis F. de Figueiredo, Adam Podhorski, Angel Rubio, Christoph Kaleta, John E. Beasley,
       Stefan Schuster and Francisco J. Planes, "Computing the shortest elementary flux modes in
       genome-scale metabolic networks" Bioinformatics, vol 25, issue 23, pp. 3158-3165, 2009.
    """
    if matrix is None:
        logger.debug("Computing integer matrix...")
        matrix = model.integerS

    if efm_model is None:
        efm_model = EFMModel(model, M=M, matrix=matrix, include_exchanges=False)
    else:
        assert isinstance(efm_model, EFMModel)

    #iteration 1
    logger.debug("Iteration 1:")
    status = efm_model.optimize()
    logger.debug("Iteration 1: %s" % status)
    [logger.debug("%s: %f" % (z, z.primal)) for z in efm_model.z_vars if z.primal > 0]
    [logger.debug("%s: %f" % (t, t.primal)) for t in efm_model.t_vars if t.primal > 0]
    for i in xrange(k-1):
        previous_solution_part = list()
        sum_part = 0
        for reaction in model.reactions:
            if not reaction in model.exchanges:
                z = efm_model.z_var(reaction.id)
                previous_solution_part.append(mul([RealNumber(z.primal), z]))
                sum_part += z.primal
                if reaction.reversibility:
                    z = efm_model.z_var(reaction._get_reverse_id())
                    previous_solution_part.append(mul([RealNumber(z.primal), z]))
                    sum_part += z.primal
        iteration_expression = add(previous_solution_part)
        iteration_constraint = model.solver.interface.Constraint(iteration_expression, ub=sum_part-1, name="iter")

        efm_model.add(iteration_constraint)

        #interation n
        logger.debug("Iteration %i:" % (i+2))
        status = efm_model.optimize()
        logger.debug("Iteration %i: %s" % (i+2, status))
        [logger.debug("%s: %f" % (z, z.primal)) for z in efm_model.z_vars if z.primal > 0]
        efm_model.remove(iteration_constraint)

    return efm_model


def fixed_size_elementary_modes(model=None, c=1, efm_model=None, size=20):
    """
    Calculates the shortest elementary flux modes using MILP with explicit indicator variables.[1]

    Parameters
    model : SolverBasedModel
        A COBRA model.
    c : int
        Threshold to determine if the variable is greater than zero, default is 1.
    em_model : EFMModel
        A prebuilt efm optimization model.
    """

    if efm_model is None:
        efm_model = EFMModel(model, M=c, include_exchanges=True, matrix=model.S, flux_type="continuous")

    with TimeMachine() as tm:
        expression = sum(efm_model.z_vars)
        size_constraint = model.solver.interface.Constraint(expression, name="fixed_size", lb=size, ub=size)

        tm(do=partial(efm_model.add, size_constraint), undo=partial(efm_model.remove, size_constraint))
        logger.debug("Start optimization")
        status = efm_model.optimize()
        logger.debug("Optimization: %s" % status)
        [logger.debug("%s: %f" % (z, z.primal)) for z in efm_model.z_vars if z.primal > 0]
        [logger.debug("%s: %f" % (t, t.primal)) for t in efm_model.t_vars if t.primal > 0]