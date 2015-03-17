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

"""
Ported from FRAMED. Original author: Kai Zhuang.
"""

from functools import partial

from scipy.integrate import ode
import numpy
import collections
import logging
from cameo.exceptions import SolveError
from cameo.util import TimeMachine

logger = logging.getLogger(__name__)


class Organism(object):
    """
    Organism describes a generic biological organism.

    Attributes
    ----------
    model : SolverBasedModel
        A constraint-based model.
    objective : see flux_analysis.simulation.fba
        The objective function of the simulation.
    constrains : a map of constraints to be applied on each simulation

    Methods
    -------
    update()

    """
    def __init__(self, model, id=None, objective=None, constraints={}):

        self.model = model

        if id:
            self.id = id
        else:
            self.id = model.id

        self.objective = objective
        self.constraints = constraints
        self.solution = None

        self.environment = None  # upon initiation, the organism is not placed in any environment

    def update(self, *args, **kw):
        """
        The update() method is used to change the internal state of the organism.
        - This method is called at each integration step.
        - One usage of this method is to describe how the FBA uptake constraint changes in response to the changes in
            the metabolite concentrations.

        ** this is an abstract method, must be implemented in strain specific subclasses **
         """

        raise NotImplementedError


class Environment(object):
    """
    This class describes a generic environment that contains a number of organisms and metabolites


    Properties
    ----------
    organisms: list of organisms
    metabolites: list of metabolites

    Methods
    -------
    update()

    """

    def __init__(self, organisms=[], metabolites=[]):
        assert isinstance(organisms, list)
        assert isinstance(metabolites, list)
        self.organisms = organisms
        self._metabolites = metabolites

    def update(self, *args, **kwargs):
        """
        The update() method is used to change the internal state of the environment.
        This method is called at each integration step.

        ** this is an abstract method, must be implemented for specific environments **
        """
        raise NotImplementedError("update() method must be implemented for the specific environment")

    @property
    def organisms(self):
        return self._organisms

    @organisms.setter
    def organisms(self, organisms):
        self._organisms = []
        for organism in organisms:
            self._add_organism(organism)

    @property
    def metabolites(self):
        return self._metabolites

    @metabolites.setter
    def metabolites(self, metabolites):
        self._metabolites = metabolites

    def _add_organism(self, organism):
        organism.environment = self
        self._organisms.append(organism)


class DynamicSystem(object):
    """
    This class describes a generic dynamic system.

    """

    def __init__(self, initial_conditions):
        self._initial_conditions = initial_conditions

    @property
    def initial_conditions(self):
        return self._initial_conditions

    @initial_conditions.setter
    def initial_conditions(self, initial_conditions):
        self._initial_conditions = initial_conditions

    def _ode_rhs(self, y, t):
        """
        This is the Right Hand Side of the system of ODE that describe the dynamic multi-species system.

        Parameters
        ----------
        y: list
            The state variables: liquid volume(y[0]), biomass concentrations, and metabolite concentrations
        t: float
            Time

        """
        raise NotImplementedError("the RHS of the ODE must be described for the each specific environment")

    def integrate(self, t0, tf, dt, initial_conditions, solver, reporters=[]):
        """
        the integrate() solves the ODE of the dynamic system using the designated solver
        :param t0: initial time
        :param tf: final time
        :param dt: time step
        :param initial_conditions: array-like -- initial conditions of the ODE system
        :param solver: the designated solver
        :return:
        """
        if solver == 'analytical':
            try:
                t, y = self.analytical_integrator(t0, tf, dt, initial_conditions, solver, reporters=reporters)
                return t, y
            except NotImplementedError:
                logger.warn('analytical solver have no been implemented yet. It will use numerical solver dopri5.')

        t, y = self.numerical_integrator(t0, tf, dt, initial_conditions, solver)
        return t, y

    def numerical_integrator(self, t0, tf, dt, initial_conditions=None, solver='dopri5', reporters=[]):
        """
        The numerical_integrator() method integrates the ODE of the dynamic system using a numerical solver
        """
        if initial_conditions:
            y0 = initial_conditions
        else:
            y0 = self._initial_conditions

        m_dfba_ode = ode(self._ode_rhs).set_integrator(solver)
        m_dfba_ode.set_initial_value(y0, t0)

        t = [t0]
        y = [y0]

        while m_dfba_ode.successful() and m_dfba_ode.t < tf:
            m_dfba_ode.integrate(m_dfba_ode.t + dt)

            t.append(m_dfba_ode.t)
            y.append(m_dfba_ode.y)
            for reporter in reporters:
                print "Reporting", t0, m_dfba_ode.t, tf, y
                reporter(t0, m_dfba_ode.t, tf, y)


        t = numpy.array(t)
        y = numpy.array(y)
        return t, y

    def analytical_integrator(self, t0, tf, dt, initial_conditions, solver, reporters=[]):
        """
        Integrates the ODE of the dynamic system using a user-defined analytical method.

        ** this is an abstract method, must be implemented for specific dynamic systems **
        """
        raise NotImplementedError


class BioReactor(Environment, DynamicSystem):
    """
    This class describes a generic bio-reactor with one influent (feed) stream and one effluent stream

    Attributes
    ----------

    organisms: list
        A list of Organism.
    metabolites: list
        A list of metabolite names.
    inflow_rate: float
        Volume inflow rate
    outflow_rate: float
        Volume outflow rate
    max_volume: float
        The maximum liquid capacity of the bioreactor.
    x_feed: float
        The concentration of organisms in the feed stream [g/L].
    s_feed: float
        The concentration of metabolites in the feed stream [mmol/L].
    delta_x: float
        Custom defined terms to dX/dt [g/L/hr].
    delta_s: list
        Custom defined terms to dX/dt [mmol/L/hr]
    initial_conditions: list
    """

    def __init__(self, organisms=[], metabolites=[], id='Generic Bioreactor', inflow_rate=0, outflow_rate=0,
                 max_volume=None, x_feed=0.0, s_feed=0.0, delta_x=None, delta_s=None, initial_conditions=[]):

        if not isinstance(organisms, collections.Iterable):
            organisms = [organisms]

        if not isinstance(metabolites, collections.Iterable):
            metabolites = [metabolites]

        self.id = id

        self.inflow_rate = inflow_rate
        self.outflow_rate = outflow_rate
        self._max_volume = max_volume

        self._s_feed = s_feed or numpy.zeros(len(metabolites))
        self._x_feed = x_feed or numpy.zeros(len(organisms))

        self._delta_x = delta_x or numpy.zeros(len(organisms))
        self._delta_s = delta_s or numpy.zeros(len(metabolites))

        Environment.__init__(self, organisms, metabolites)
        DynamicSystem.__init__(self, initial_conditions)

    @property
    def metabolites(self):
        return Environment.metabolites.fget(self)

    @metabolites.setter
    def metabolites(self, metabolites):
        Environment.metabolites.fset(self, metabolites)
        if len(metabolites) != len(self._s_feed):
            self._s_feed = numpy.zeros(len(self.metabolites))
            self._delta_s = numpy.zeros(len(self.metabolites))

    @property
    def organisms(self):
        return Environment.organisms.fget(self)

    @organisms.setter
    def organisms(self, organisms):
        Environment.organisms.fset(self, organisms)
        if len(organisms) != len(self._x_feed):
            self._x_feed = numpy.zeros(len(self.organisms))
            self._delta_x = numpy.zeros(len(self.organisms))


    @property
    def max_volume(self):
        return self._max_volume

    @property
    def x_feed(self):
        return self._x_feed

    @property
    def s_feed(self):
        return self._s_feed

    @property
    def delta_x(self):
        return self._delta_x

    @property
    def delta_s(self):
        return self._delta_s

    def update(self, time, volume, s, x):
        if self.max_volume:
            if volume > self.max_volume:
                raise ValueError('Volume %f exceeds limit %f' % (volume, self.max_volume))

    def _ode_rhs(self, t, y):
        """
        the RHS of the ODE that describe the bioreactor system
        :param y:
            y[0]: volume
            y[1] to y[number_of_organisms]: biomass of the organisms
            y[number_of_organisms + 1] to y[-1] concentration of metabolites

        :param t: time
        :return: dy
        """
        number_of_organisms = len(self._organisms)
        number_of_metabolites = len(self._metabolites)
        assert (len(y) == 1 + number_of_organisms + number_of_metabolites), \
            "%s == t + [%s, %s]" % (y, self.organisms, self.metabolites)

        dy = numpy.zeros(len(y))

        # creating class variables V, X, S, time from y and t.
        # making them class variables so that class methods like update() can access them
        volume = y[0]
        x = y[1:number_of_organisms + 1]
        s = y[number_of_organisms + 1:]


        # assigning growth rates and metabolic production/consumption rates here
        # in this method, these rates are calculated using FBA

        # fluxes through metabolites
        vs = numpy.zeros([number_of_organisms, number_of_metabolites])
        # growth rates of organisms
        mu = numpy.zeros([number_of_organisms])

        for i, organism in enumerate(self._organisms):
            organism.update(volume, x[i], s)   # updating the internal states of the organism
            # eg. updating the uptake constraints based on metabolite concentrations
            with TimeMachine() as tm:
                for reaction_id, constraints in organism.constraints.items():
                    reaction = organism.model.reactions.get_by_id(reaction_id)
                    tm(do=partial(setattr, reaction, 'lower_bound', constraints[0]),
                       undo=partial(setattr, reaction, 'lower_bound', reaction.lower_bound))
                    tm(do=partial(setattr, reaction, 'upper_bound', constraints[1]),
                       undo=partial(setattr, reaction, 'upper_bound', reaction.upper_bound))
                if organism.objective:
                    tm(do=partial(setattr, organism.model, 'objective', organism.objective),
                       undo=partial(setattr, organism.model, 'objective', organism.model.objective))

                try:
                    organism.solution = organism.model.solve()
                    mu[i] = organism.solution.f

                    for j, metabolite in enumerate(self._metabolites):
                        vs[i, j] = organism.solution.get_primal_by_id(metabolite)
                except SolveError:
                    organism.solution = None
                    mu[i] = 0
                    for j, metabolite in enumerate(self._metabolites):
                        vs[i, j] = 0

        # updating the internal states of the bioreactor
        # eg. flow rates, feed concentrations, and custom defined dX/dt and dS/dt terms
        self.update(t, volume, x, s)

        # calculating the rates of change of reactor volume[L], biomass [g/L] and metabolite [mmol/L]
        # dV/dt [L/hr]
        dy[0] = self.inflow_rate - self.outflow_rate
        # dX/dt [g/L/hr]
        dy[1:number_of_organisms + 1] = mu * x + self.inflow_rate / volume * (self.x_feed - x) + self.delta_x
        # dS/dt [mmol/L/hr]
        dy[number_of_organisms + 1:] = numpy.dot(x, vs) + self.inflow_rate / volume * (self.s_feed - s) + self.delta_s

        return dy

    def calculate_yield_from_dfba(self, dfba_solution, r_substrate, r_product):
        """
        Abstract used for calculating product yield from dFBA solution.
        This is useful for certain analysis methods (eg. DySScO).

        This should be implemented for each specific bioreactor.
        """
        raise NotImplementedError

    def calculate_titer_from_dfba(self, dfba_solution, r_target):
        """
        Abstract used for calculating product titer from dFBA solution.
        This is useful for certain analysis methods (eg. DySScO).

        This should be implemented for each specific bioreactor.
        """
        raise NotImplementedError

    def calculate_productivity_from_dfba(self, dfba_solution, r_target):
        """
        Abstract used for calculating productivity from dFBA solution.
        This is useful for certain analysis methods (eg. DySScO).

        This should be implemented for specific each bioreactor.
        """
        raise NotImplementedError
