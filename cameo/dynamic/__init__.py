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
from functools import partial
import numpy as np
from pandas import DataFrame

from scipy.integrate import ode
from cameo import plot_utils, fba
from cameo.core.result import Result
from cameo.dynamic.bioreactors import batch_reactor
from cameo.exceptions import SolveError
from cameo.util import TimeMachine

import logging

logger = logging.getLogger(__name__)


class DynamicFBAResult(Result):
    def data_frame(self):
        d = {
            "Time": self.time,
            "Volume": self.volume,
        }
        d.update(self.biomass)
        d.update(self.metabolite_concentrations)
        return DataFrame()

    def __init__(self, reactor, time, volume, biomass, metabolite_concentrations, *args, **kwargs):
        super(DynamicFBAResult, self).__init__(*args, **kwargs)
        self.time = time
        self.volume = volume
        self.biomass = biomass
        self.metabolite_concentrations = metabolite_concentrations
        self.reactor = reactor

    def __getitem__(self, item):
        try:
            return self.metabolite_concentrations[item]
        except KeyError:
            return self.biomass[item]

    def plot(self):
        plot_utils.plot_dfba_solution(self)

    def __len__(self):
        return len(self.biomass) + len(self.metabolite_concentrations)

    def __repr__(self):
        return "DFBA Result (t=%s, v=%s)\nMetabolites:\n%s\nStrains:%s" % (
            self.time,
            self.volume,
            "".join(["\n\t%s: %s (mmol/L)" % (k, v) for k, v in self.metabolite_concentrations.iteritems()]),
            "".join(["\n\t%s: %s (g/L)" % (k, v) for k, v in self.biomass.iteritems()]),
        )


def _bioreactor_ode(t, y, metabolites, models, models_dynamics, reactor, simulation_method, simulation_kwargs):
        """
        Description the reactor system.

        Parameters
        ----------
        y: np.array
            y[0]: volume;
            y[1; number_of_organisms]: biomass of the organisms;
            y[number_of_organisms + 1:]  concentration of metabolites.
        t: float
            time

        Returns
        -------

        dy: np.array
            The variation of volume, metabolite concentrations and biomass in the reactor.
        """

        n = len(models)
        m = len(metabolites)

        dy = np.zeros(len(y))

        # creating class variables from y and t.
        volume, inflow_rate, outflow_rate = y[0], y[1], y[2]
        logger.debug("Volume: %f, In: %f, Out: %f" % (volume, inflow_rate, outflow_rate))

        x, delta_x, x_feed = y[3:n+3], y[n+3: 2*n+3], y[2*n+3: 3*n+3]
        for i, model in enumerate(models):
            logger.debug("Model %s => Cdw: %f (g/L), DeltaX: %f, XFeed: %f" % (model.id, x[i], delta_x[i], x_feed[i]))

        s, delta_s, s_feed = y[3*n+3:3*n+m+3], y[3*n+m+3:3*n+2*m+3], y[3*n+2*m+3:3*n+3*m+3]
        for i, met in enumerate(metabolites):
            logger.debug("%s => C: %f (mmol/L), DeltaS: %f, SFeed: %f" % (met, s[i], delta_s[i], s_feed[i]))

        # fluxes through metabolites
        vs = np.zeros([len(models), len(metabolites)])
        # growth rates of organisms
        mu = np.zeros([len(models)])

        for i, model in enumerate(models):
            # updating the internal states of the organism
            # eg. updating the uptake constraints based on metabolite concentrations
            constraints, objective = models_dynamics[i](volume, x[i], metabolites, s)

            with TimeMachine() as tm:
                for reaction_id, c in constraints.items():
                    logger.debug("Reaction: %s (lb=%f, ub=%f)" % (reaction_id, c[0], c[1]))
                    r = models[i].reactions.get_by_id(reaction_id)
                    tm(do=partial(setattr, r, 'lower_bound', c[0]), undo=partial(setattr, r, 'lower_bound', r.lower_bound))
                    tm(do=partial(setattr, r, 'upper_bound', c[1]), undo=partial(setattr, r, 'upper_bound', r.upper_bound))

                try:
                    solution = simulation_method(model, objective=objective, raw=True, **simulation_kwargs)
                    mu[i] = solution[model.objective]

                    for j, metabolite in enumerate(metabolites):
                        vs[i, j] = solution[metabolite]
                except SolveError:
                    mu[i] = 0
                    for j, metabolite in enumerate(metabolites):
                        vs[i, j] = 0

        # updating the internal states of the bioreactor
        # eg. flow rates, feed concentrations, and custom defined dX/dt and dS/dt terms
        d_inflow_rate, d_outflow_rate, d_delta_x, d_x_feed, d_delta_s, d_s_feed = reactor(t, volume, x, s,
                                                                                          inflow_rate, outflow_rate,
                                                                                          delta_x, x_feed,
                                                                                          delta_s, s_feed)

        # calculating the rates of change of reactor volume[L], biomass [g/L] and metabolite [mmol/L]
        # dV/dt [L/hr]
        dy[0:3] = [inflow_rate - outflow_rate, d_inflow_rate, d_outflow_rate]
        # dX/dt [g/L/hr]
        dx = mu * x + inflow_rate / volume * (x_feed - x) + delta_x
        # dS/dt [mmol/L/hr]
        ds = np.dot(x, vs) + inflow_rate / volume * (s_feed - s) + delta_s

        dy[3:3*n+3] = np.append(dx, [d_delta_x , d_x_feed])
        dy[3*n+3:3*n+3*m+3] = np.append(ds, [d_delta_s, d_s_feed])

        return dy


def _extend_list(alist, value, message):
    if isinstance(value, list):
        assert len(value) > 0, message
        value = [value[i] if i < len(value) else 0 for i, _ in enumerate(alist)]
    elif isinstance(value, (float, int)):
        value = [value for _ in alist]

    return value


def _ensure_metabolite_params(metabolites, initial, delta_s, s_feed):
    initial = _extend_list(metabolites, initial, "At least one initial concentration must be provided")
    delta_s = _extend_list(metabolites, delta_s, "At least one variation rate must be provided")
    s_feed = _extend_list(metabolites, s_feed, "At least one feed rate must be provided")

    return initial, delta_s, s_feed


def _ensure_biomass_params(models, initial, delta_x, x_feed):
    initial = _extend_list(models, initial, "At least one initial biomass must be provided")
    delta_x = _extend_list(models, delta_x, "At least one variation rate must be provided")
    x_feed = _extend_list(models, x_feed, "At least one feed rate must be provided")

    return initial, delta_x, x_feed


def batch_dfba(metabolites=None, initial_concentrations=None, models=None, initial_biomass=None, models_dynamics=None,
               initial_volume=None, t0=0, tf=10, dt=1, simulation_method=fba, simulation_kwargs=None,
               ode_solver='dopri5', ode_kwargs=None, reporters=None):
    return dfba(batch_reactor,
                metabolites=metabolites,
                initial_concentrations=initial_concentrations,
                delta_s=0.,
                s_feed=0.,
                models=models,
                initial_biomass=initial_biomass,
                models_dynamics=models_dynamics,
                delta_x=0,
                x_feed=0,
                initial_volume=initial_volume,
                inflow_rate=0.,
                outflow_rate=.0,
                t0=t0,
                tf=tf,
                dt=dt,
                simulation_method=simulation_method,
                simulation_kwargs=simulation_kwargs,
                ode_solver=ode_solver,
                ode_kwargs=ode_kwargs,
                reporters=reporters)


def dfba(reactor, metabolites=None, initial_concentrations=None, delta_s=None, s_feed=None,
         models=None, initial_biomass=None, models_dynamics=None, delta_x=None, x_feed=None,
         initial_volume=None, inflow_rate=None, outflow_rate=None, t0=0, tf=10, dt=1,
         simulation_method=fba, simulation_kwargs=None, ode_solver='dopri5', ode_kwargs=None,
         reporters=None):
    initial_conc, delta_s, s_feed = _ensure_metabolite_params(metabolites, initial_concentrations, delta_s, s_feed)
    initial_biom, delta_x, x_feed = _ensure_metabolite_params(models, initial_biomass, delta_x, x_feed)
    if reporters is None:
        reporters = []
    if simulation_kwargs is None:
        simulation_kwargs = {}
    if ode_kwargs is None:
        ode_kwargs = {}

    return _dfba(reactor,
                 metabolites=metabolites,
                 initial_concentrations=initial_conc,
                 delta_s=delta_s,
                 s_feed=s_feed,
                 models=models,
                 initial_biomass=initial_biom,
                 models_dynamics=models_dynamics,
                 delta_x=delta_x,
                 x_feed=x_feed,
                 initial_volume=initial_volume,
                 inflow_rate=inflow_rate,
                 outflow_rate=outflow_rate,
                 t0=t0,
                 tf=tf,
                 dt=dt,
                 simulation_method=simulation_method,
                 simulation_kwargs=simulation_kwargs,
                 ode_solver=ode_solver,
                 ode_kwargs=ode_kwargs,
                 reporters=reporters)


def _dfba(reactor, metabolites=None, initial_concentrations=None, delta_s=None, s_feed=None,
          models=None, initial_biomass=None, models_dynamics=None, delta_x=None, x_feed=None,
          initial_volume=None, inflow_rate=None, outflow_rate=None, t0=0, tf=10, dt=1,
          simulation_method=None, simulation_kwargs=None, ode_solver='dopri5', ode_kwargs=None,
          reporters=None):

    assert initial_volume > 0, "Cannot starta fermentation without volume!"
    assert len(metabolites) == len(initial_concentrations), \
        "Initial concentration and metabolites must have the same size."
    assert len(metabolites) == len(delta_s), \
        "Metabolite variation rates (delta_s) and metabolites must have the same size."
    assert len(metabolites) == len(s_feed), "Metabolite feed rates (s_feed) and metabolites must have the same size."
    assert len(models) == len(models_dynamics), "Dynamic behavior must be defined for every model"
    assert len(models) == len(initial_biomass), "Initial biomass and models must have the same size."
    assert len(models) == len(delta_x), "Biomass variation rates (delta_x) and models must have the same size."
    assert len(models) == len(x_feed), "Biomass feed rates (x_feed) and models must have the same size."

    y0 = [initial_volume, inflow_rate, outflow_rate] +\
         initial_biomass + delta_x + x_feed + initial_concentrations + delta_s + s_feed

    dfba_ode = ode(_bioreactor_ode).set_integrator(ode_solver, **ode_kwargs)
    dfba_ode.set_initial_value(y0, t0)
    dfba_ode.set_f_params(metabolites, models, models_dynamics, reactor, simulation_method, simulation_kwargs)

    t = [t0]
    y = [y0]

    while dfba_ode.successful() and dfba_ode.t < tf:
        dfba_ode.integrate(dfba_ode.t + dt)
        t.append(dfba_ode.t)
        y.append(dfba_ode.y)
        for reporter in reporters:
            print "Reporting", t0, dfba_ode.t, tf, y
            reporter(t0, dfba_ode.t, tf, y)

    t = np.array(t)
    y = np.array(y)

    biomass = {}
    i = 2
    for model in models:
        i += 1
        biomass[model.id] = y[:, i]

    metabolite_concentrations = {}
    i = 3 * len(models) + 2
    for metabolite in metabolites:
        i += 1
        metabolite_concentrations[metabolite] = y[:, i]

    return DynamicFBAResult(reactor, t, y[:, 0], biomass, metabolite_concentrations)