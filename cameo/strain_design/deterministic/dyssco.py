# Copyright 2015 Novo Nordisk Foundation Center for Biosustainability,
# Technical University of Denmark.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

""" This module implements classes and methods useful for the Dynamic Strain Scanning Optimization strategy

Cite:
Kai Zhuang et al. 2013. Dynamic strain scanning optimization: an efficient strain design strategy for balanced yield,
titer, and productivity.
"""
from functools import partial

import cobra
import math
from pandas import DataFrame

from cameo import phenotypic_phase_plane, config
from cameo.dynamic.bioreactor import Organism, IdealFedBatch, BioReactor
from cameo.dynamic import dfba
from cameo.strain_design import StrainDesignMethod
from cameo.config import in_ipnb
from cameo import plot_utils
from ipython_notebook_utils import ProgressBar as IPythonNBPorgressBar


TITER = 'product_titer'
YIELD = 'product_yield'
PRODUCTIVITY = 'productivity'


class _DFBAPerformanceEvaluator(object):
    def __init__(self, reactor, substrate, product, tf, dt, solver):
        assert isinstance(reactor, BioReactor)
        assert isinstance(tf, (float, int))
        assert isinstance(dt, (float, int))
        self.reactor = reactor
        self.product = product
        self.substrate = substrate
        self.tf = tf
        self.dt = dt
        self.solver = solver

    def __call__(self, strain, reporter, *args, **kwargs):
        self.reactor.organisms = [strain]
        if reporter is not None:
            reporters = [reporter]
        else:
            reporters = []
        try:
            solution = dfba(self.reactor, 0, self.tf, self.dt, solver=self.solver, reporters=reporters)

            return {
                TITER: self.reactor.calculate_titer_from_dfba(solution, self.product.id),
                YIELD: self.reactor.calculate_yield_from_dfba(solution, self.substrate.id, self.product.id),
                PRODUCTIVITY: self.reactor.calculate_productivity_from_dfba(solution, self.product.id)
            }
        except:
            return {
                TITER: 0,
                YIELD: 0,
                PRODUCTIVITY: 0
            }


class DySScO(StrainDesignMethod):

    def __init__(self, organism=None, product=None, substrate=None, reactor=None,
                 dt=1, tf=10, update_function=None, dfba_solver='dopri5', *args, **kwargs):
        super(DySScO, self).__init__(*args, **kwargs)
        assert isinstance(organism, Organism), "organism must be instance of Organism"
        if isinstance(product, str):
            product = organism.model.reactions.get_by_id(product)
        if isinstance(substrate, str):
            substrate = organism.model.reactions.get_by_id(substrate)
        assert isinstance(product, cobra.Reaction), "product must be a reaction id or an instance of (cobra.Reaction)"
        assert isinstance(substrate, cobra.Reaction), "substrate must be a reaction id or an instance of (cobra.Reaction)"
        assert isinstance(reactor, BioReactor)
        self.organism = organism
        self.product = product
        self.substrate = substrate
        self.reactor = reactor
        self.tf = tf
        self.dt = dt
        self.dfba_solver = dfba_solver
        self.update_function=update_function

    def _generate_strains(self, number_of_strains, view=config.default_view):
        envelope = phenotypic_phase_plane(self.organism.model,
                                          variables=[self.product],
                                          points=number_of_strains,
                                          view=view)

        plot_utils.plot_production_envelope(envelope, self.product.id, highligt=range(1, len(envelope) - 1))

        for i in range(1, len(envelope) - 1):
            row = envelope.loc[i]
            strain_id = "_%i" % i
            strain = Organism(self.organism.model, strain_id, self.organism.objective, dict(self.organism.constraints))
            strain.update = partial(self.update_function, strain)
            strain.constraints[self.product.id] = (math.floor(row[self.product.id]), math.ceil(row[self.product.id]))
            yield strain

    def run(self, ode_solver="dopri5", number_of_strains=10, view=config.default_view):
        strains = self._generate_strains(number_of_strains+2, view=view)
        evaluator = _DFBAPerformanceEvaluator(self.reactor, self.product, self.substrate,
                                              self.tf, self.dt, self.dfba_solver)

        reporters = []
        if in_ipnb():
            for i in range(1, number_of_strains + 1):
                progress_bar = IPythonNBPorgressBar(size=self.tf, label="Strain %i" % i)
                progress_bar.start()
                reporters.append(lambda t, t0, tf, y: progress_bar.set(t))

        performance = DataFrame.from_records(list(view.map(evaluator, strains, reporters)))

        plot_utils.plot_dfba_performance(performance, [s.id for s in strains])
        return performance


if __name__ == "__main__":
    from cameo import load_model
    import os

    path = os.path.abspath(os.path.dirname(__file__))
    model = load_model(os.path.join(path, "../../../../tests/data/ecoli_core_model.xml"))

    organism = Organism(model=model)

    def update(self, volume, growth_rate, substrates):
        env = self.environment

        index = env.metabolites.index('EX_glc_e')
        vlb_glc = float(-10 * substrates[index] / (substrates[index] + 1))
        self.constraints['EX_glc_e'] = (vlb_glc, 0)

        self.constraints['EX_ac_e'] = (0, 10000)
        self.constraints['EX_o2_e'] = (-10, 10000)

    model.update = update

    product = 'EX_ac_e'
    substrate = 'EX_glc_e'
    reactor = IdealFedBatch(primary_substrate=substrate,
                            metabolites=['EX_glc_e', 'EX_ac_e', 'EX_o2_e'],
                            initial_conditions=[1.0, 0.1, 0.1])

    dissco = DySScO(organism=organism, product=product, substrate=substrate, reactor=reactor, update_function=update)



#
#
# def make_envelope_strains(base_organism, r_substrate, r_target, N=10, constraints=None):
#     """
#     Create N strains along the product envelope.
#         (Used for Steps 1 and 2 of DySScO strategy)
#
#     Arguments:
#         base_organism: Organism -- the host organism used to product the target product
#         r_substrate: str -- the rxn id of the substrate
#         r_target: str -- the rxn id of the target product
#         N: int -- the number of strains to be generated along the production envelope
#         constraints: dict -- custom constraints
#
#     Returns:
#         strains: list of Organism -- N strains that are fixed to the production envelope
#     """
#     strains = []
#
#     # add custom constraints to base_organism
#     if constraints:
#         base_organism.fba_constraints.update(constraints)
#
#     # create the product envelope
#     xvals, ymins, ymaxs = production_envelope(base_organism.model, r_target, steps=N, constraints=constraints)
#
#     # finding the maximum r_substrate uptake rate
#     if r_substrate in base_organism.fba_constraints:
#         vSmax = base_organism.fba_constraints[r_substrate][0]
#     else:
#         vSmax = base_organism.model.bounds[r_substrate][0]
#
#     # create new strains along the production envelope
#     for i, mu in enumerate(xvals):
#         # creating a new strain
#         strain = deepcopy(base_organism)                        # create a deepcopy of the base_organism
#         strain.fba_constraints[r_target] = (ymaxs[i], ymaxs[i])   # fix target production at ymax[i]
#         #strain.Y = float(-ymaxs[i]/vSmax)                       # store the yield of the strain
#         #strain.mu = mu                                          # store the growth rate of the strain
#         strain.id = base_organism.id + '_mu_' + str(round(mu, 3))
#         strains.append(strain)
#
#     return strains
#
#
# def calculate_performances(strains, bioreactor, r_substrate, r_target, t0, tf, dt, initial_conditions=[],
#                           dfba_solver='dopri5', additional_yields=[], verbose=False, save_dfba_solution=False,
#                           func_dfba2yield=None):
#     """
#     calculates the performances of a list of strains in a given bioreactor
#
#     Arguments:
#         strains: list (of Organism)
#         bioreactor: Bioreactor
#         r_substrate: str -- reaction id of the substrate uptake reaction
#         r_target: str -- reaction id of the target metabolite exchange reaction
#         t0: float -- initial time
#         tf: float -- final time
#         dt: float -- time step
#         initial_conditions: list (of float) -- the initial conditions in the order of V0, X0, S0 (default: None)
#         dfba_solver: str -- ODE solver.  (default: 'dopri5')
#         additional_yields: list (of str) -- the reaction ids of the additional yields (yields other than target yield)
#                                             to be calculated.
#         verbose: bool -- Verbosity control.  (default: False).
#         save_dfba_solution: bool -- controls whether dfba solutions are returned (default: False)
#         func_dfba2yield: function -- if None, yield is calculated from FBA solutions.
#                                      if a function is passed in, yield is calculated from dFBA solutions.
#     Returns:
#         performances: list (of Dict) -- a list of dictionaries, each entry contains the calculated performance metrics
#         of a strain
#     """
#
#     performances = []
#
#     for strain in strains:
#         performance = calculate_performance(strain, bioreactor, r_substrate, r_target, t0, tf, dt,
#                                                 initial_conditions, dfba_solver, additional_yields, verbose,
#                                                 save_dfba_solution)
#         performances.append(performance)
#
#     return performances
#
#
# def calculate_performance(strain, bioreactor, r_substrate, r_target, t0, tf, dt, initial_conditions=[],
#                           dfba_solver='vode', additional_yields=[], verbose=False, save_dfba_solution=False):
#     """
#     calculates the performances of a list of strains in a given bioreactor
#
#     Arguments:
#         strain: Organism
#         bioreactor: Bioreactor
#         r_substrate: str -- reaction id of the substrate uptake reaction
#         r_target: str -- reaction id of the target metabolite exchange reaction
#         t0: float -- initial time
#         tf: float -- final time
#         dt: float -- time step
#         initial_conditions: list (of float) -- the initial conditions in the order of V0, X0, S0 (default: None)
#         dfba_solver: str -- ODE solver.  (default: 'dopri5')
#         additional_yields: list (of str) -- the reaction ids of the additional yields (yields other than target yield)
#                                             to be calculated.
#         verbose: bool -- Verbosity control.  (default: False).
#         save_dfba_solution: bool -- controls whether dfba solutions are returned (default: False)
#
#     Returns:
#         performance: Dict -- contains the calculated performance metrics of a strain
#     """
#     performance = {'strain_id': strain.id}
#     r_biomass = strain.model.detect_biomass_reaction()
#
#     # perform FBA simulation
#     if verbose:
#         print 'Performing FBA simulation.'
#     if hasattr(strain, 'solver'):
#         fba_solution = strain.solver.solve_lp(strain.fba_objective, constraints=strain.fba_constraints)
#     else:
#         strain.solver = solver_instance()
#         strain.solver.build_problem(strain.model)
#         fba_solution = strain.solver.solve_lp(strain.fba_objective, constraints=strain.fba_constraints)
#
#     # growth, substrate uptake, and target production rates from FBA solution
#     v_biomass = fba_solution.values[r_biomass]
#     v_target = fba_solution.values[r_target]
#     v_substrate = fba_solution.values[r_substrate]
#
#     # if the strain does not grow, set growth, titer, productivity to zero, and calculate yield from FBA
#     if v_biomass <= 10**-6:
#         performance['growth_rate'] = 0
#         performance['biomass_yield'] = 0
#         performance['product_titer'] = 0
#         performance['productivity'] = 0
#         performance['product_yield'] = - v_target/v_substrate
#
#         for r_id in additional_yields:
#             pid = 'yield_' + r_id.lstrip('R_EX_').rstrip('_e')
#             performance[pid] = - fba_solution.values[r_id] / v_substrate
#         dfba_solution = None
#
#         if verbose:
#             print 'none growing'
#
#     # if the strain grows but does not produce, set yield, titer, productivity to zero, and calculate growth from FBA
#     elif v_target <= 10**-6:
#         performance['growth_rate'] = v_biomass
#         performance['biomass_yield'] = - v_biomass/v_substrate
#         performance['product_titer'] = 0
#         performance['productivity'] = 0
#         performance['product_yield'] = 0
#
#         for r_id in additional_yields:
#             pid = 'yield_' + r_id.lstrip('R_EX_').rstrip('_e')
#             performance[pid] = 0
#         dfba_solution = None
#
#         if verbose:
#             print 'none producing'
#
#     # if the strain both grows and produces, perform dFBA simulation, and calculate yield, titer, productivity
#     else:
#         # perform dFBA simulation
#         bioreactor.organisms([strain])
#
#         if verbose:
#             print 'Performing dFBA simulation.'
#         dfba_solution = dfba(bioreactor, t0, tf, dt, initial_conditions, solver=dfba_solver, verbose=verbose)
#
#         # calculate yield using dFBA solution if the method is known,
#         # otherwise calculate yield using FBA solution
#         try:
#             performance['product_yield'] = bioreactor.calculate_yield_from_dfba()
#         except NotImplementedError:
#             performance['product_yield'] = - v_target / v_substrate
#
#         # calculate titer and productivity using dFBA solution.
#             # if specific calculation functions exist for the bioreactor, use them
#             # otherwise, assume ideal fedbatch reactor.
#
#         try:
#             performance['product_titer'] = bioreactor.calculate_titer_from_dfba()
#         except NotImplementedError:
#             performance['product_titer'] = dfba_solution[r_target].max()
#
#         try:
#             performance['productivity'] = bioreactor.calculate_productivity_from_dfba()
#         except NotImplementedError:
#              # the index at which the production is finished.  round off to 4 decimal places
#             index = numpy.around(dfba_solution[r_target]).argmax()
#             performance['productivity'] = performance['product_titer']/dfba_solution['time'][index]
#
#         performance['growth_rate'] = v_biomass
#         performance['biomass_yield'] = - v_biomass/v_substrate
#
#         # calculate additional yields
#         for r_id in additional_yields:
#             id = 'yield_' + r_id.lstrip('R_EX_').rstrip('_e')
#             performance[id] = - fba_solution.values[r_id] / v_substrate
#
#     if save_dfba_solution:
#         performance['dfba_solution'] = dfba_solution
#
#     return performance
