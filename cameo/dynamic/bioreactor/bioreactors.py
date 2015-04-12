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

from cameo.dynamic.bioreactor.base import BioReactor

# value definition for oxygen_availability flag
ANAEROBIC = 0
AEROBIC = 1
MICRO_AEROBIC = 2


class BioReactorOX(BioReactor):
    """
    Bioreactor class with oxygen_availability flag
    """

    def __init__(self, oxygen_availability=None, **kwargs):
        super(BioReactorOX, self).__init__(**kwargs)
        self.oxygen_availability = oxygen_availability


class IdealBatch(BioReactorOX):
    """
    This class describes an ideal batch reactor.
        - flow_rate_in, flow_rate_out, Xfeed, Sfeed are all set to zero (no feeding in batch reactor).
    """

    def __init__(self, id="IdealBatch", **kwargs):
        """
        Arguments:
            organisms: list of Organism
            metabolites: list of string
            volume_max: float -- liquid capacity of the bioreactor
            deltaX: custom defined terms to dX/dt [g/L/hr]
            deltaS: list of float -- special custom defined terms to dX/dt [mmol/L/hr]
            initial_conditions: list of float
        """
        kwargs['id'] = kwargs.get('id', id)
        kwargs['inflow_rate'] = 0.0
        kwargs['outflow_rate'] = 0.0
        super(IdealBatch, self).__init__(**kwargs)

    def update(self, time, volume, s, x):
        initial_volume = self.initial_conditions[0]
        if volume != initial_volume:
            raise ValueError("Batch volume cannot be changed (initial=%f, current=%f)" % (initial_volume, volume))

    def calculate_yield_from_dfba(self, dfba_solution, r_substrate, r_product):
        """
        calculates the product yield from dFBA solution
        """
        s_f = dfba_solution[r_substrate][-1]
        s0 = dfba_solution[r_substrate][0]
        p_f = dfba_solution[r_product][-1]
        product_yield = p_f / (s_f - s0)

        return product_yield

        # TODO: Check whether is max or last time point
    def calculate_titer_from_dfba(self, dfba_solution, r_product):
        return dfba_solution[r_product][-1]

    def calculate_productivity_from_dfba(self, dfba_solution, r_product):
        return dfba_solution[r_product][-1]/dfba_solution.time[-1]


class IdealFedBatch(BioReactorOX):
    """
    This class describes an ideal fed batch reactor with a single primary substrate.
        - The flow_rate_in is automatically adjusted using the following rules:
            - otherwise, calculates flow_rate_in so that substrate concentration is maintained (d_substrate/dt = 0)
        - The primary substrate (usually the carbon & energy source) can be specified in the __init__() method.
          If it is not specified, the first element of metabolites is assumed to be the substrate
    """

    def __init__(self, id='IdealFedBatch', primary_substrate=None, **kwargs):
        kwargs['id'] = kwargs.get('id', id)
        super(IdealFedBatch, self).__init__(**kwargs)

        # if the substrate is unspecified, it is assumed to be metabolites[0]
        if primary_substrate:
            assert (primary_substrate in self.metabolites)
            self.primary_substrate = primary_substrate
        else:
            self.primary_substrate = self.metabolites[0]

    def update(self, time, volume, x, s):
        """
        calculates the flow_rate_in of the fed batch reactor is calculated here.
            - if liquid volume >= volume_max, then the tank is full, set flow_rate_in to zero
            - otherwise, calculate the flow rate so that d_substrate/dt = 0

        :param time: float -- the simulation time.  can be used in subclasses to trigger time-specific events
        """
        if self.max_volume and (volume >= self.max_volume):
            self.inflow_rate = 0
        else:
            met_id = self._metabolites.index(self.primary_substrate)
            self.inflow_rate = 0
            for org_id, organism in enumerate(self._organisms):
                if organism.solution:
                    vs = organism.solution.get_primal_by_id(self.primary_substrate)
                    self.inflow_rate -= vs * x[org_id] * volume / (self.s_feed[met_id] - s[met_id])
