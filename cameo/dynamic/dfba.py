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

""" This module implements the Dynamic Flux Balance Analysis with support for multi-organism communities.
This multi-organism version of the dFBA is derived from the Dynamic Multi-species Metabolic Modeling (DyMMM) framework.

For dFBA, please cite:
    Mahadevan et al. 2002. Dynamic flux balance analysis of diauxic growth in Escherichia coli.

For DyMMM, please cite:
    Zhuang et al. 2011. Genome-scale dynamic modeling of the competition between Rhodoferax and Geobacter in anoxic
        subsurface environments.
    Zhuang et al. 2012. The design of long-term effective uranium bioremediation strategy using a community metabolic
        model.
"""

from collections import OrderedDict
from cameo import plot_utils




def dfba(bioreactor, t0, tf, dt, initial_conditions=None, solver='dopri5', reporters=[]):
    """
    Dynamic Flux Balance Analysis with Multi-organism support

    Parameters
    ----------
        bioreactor: BioReactor
            A BioReactor instance with defined organisms and metabolites.
        t0: float
            Initial time.
        tf: float
            Final time.
        dt: float
            Time step.
        initial_conditions: list
            the initial conditions in the order of V0, X0, S0 (default: None).
        solver: str
            ODE solver (default: 'dopri5').
        verbose: bool
            Verbosity control (default: False).

    Returns
    -------
        results: OrderedDict -- simulation results
    """
    t, y = bioreactor.integrate(t0, tf, dt, initial_conditions, solver, reporters=reporters)
    i = 0

    organisms = {}

    for organism in bioreactor.organisms:
        i += 1
        organisms[organism.id] = y[:, i]

    metabolites = {}
    for metabolite in bioreactor.metabolites:
        i += 1
        metabolites[metabolite] = y[:, i]

    return DynamicFBAResult(bioreactor, t, y[:, 0], organisms, metabolites)


def combinatorial_dfba(organisms, bioreactors, t0, tf, dt, initial_conditions=None, solver='dopri5'):
    """
    Run dFBA for all possible combinations of the given organisms and reactors.
    For example,
        given two organisms "ecoli" and "scerevisiae", and two reactors "batch" and "fedbatch",
        the call combinatorial_dfba([ecoli, scerevisiae], [batch, fedbatch], t0, ft, dt] will perform four simulations:
            1. ecoli in batch
            2. ecoli in fedbatch
            3. scerevisiae in batch
            4. scerevisiae in fedbtach

    Arguments:
        organisms: list of Organism
        bioreactors: list of Bioreactor
        t0: float -- initial time
        tf: float -- final time
        dt: float -- time step
        initial_conditions: list of float -- the initial conditions in the order of V0, X0, S0 (default: None)
        solver: str -- ODE solver.  (default: 'dopri5')
        verbose: bool -- Verbosity control.  (default: False).

    Returns:
        result: OrderedDict -- a dictionary of dfba results
    """

    result = OrderedDict()

    for organism in organisms:
        for bioreactor in bioreactors:
            bioreactor.organisms = [organism]
            dfba_result = dfba(bioreactor, t0, tf, dt, initial_conditions, solver=solver)
            result[organism.id, bioreactor.id] = dfba_result

    return result