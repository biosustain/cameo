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
from itertools import izip

import logging
from pandas import DataFrame
from progressbar import ProgressBar
from cameo import config, flux_variability_analysis
from cameo.parallel import SequentialView, MultiprocessingView
from cameo.strain_design import StrainDesignMethod
from cameo.flux_analysis.analysis import phenotypic_phase_plane


logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


class DifferentialFVA(StrainDesignMethod):
    """Differential flux variability analysis.__init__.py

    Compare flux ranges of a reference model to and a set of models that
    are parameterized to lie on a grid of evenly space points in the
    production envelope of a model.

    ...
    """

    @staticmethod
    def _interval_overlap(interval1, interval2):
        return min(interval1[1] - interval2[0], interval2[1] - interval1[0])

    @classmethod
    def _interval_gap(cls, interval1, interval2):
        overlap = cls._interval_overlap(interval1, interval2)
        if overlap >= 0:
            return 0
        else:
            if abs(interval1[1]) > abs(interval2[1]):
                return overlap
            else:
                return -1 * overlap

    def __init__(self, design_space_model=None, reference_model=None, target=None, variables=[],
                 normalize_fluxes_by=None, points=20):
        super(DifferentialFVA, self).__init__()
        self.design_space_model = design_space_model
        self.reference_model = reference_model
        self.target = target
        self.variables = variables
        self.points = points
        self.envelope = None
        self.grid = None
        self.normalize_fluxes_by = normalize_fluxes_by

    def _init_search_grid(self):
        self.envelope = phenotypic_phase_plane(
            self.design_space_model, self.variables, objective=self.target, points=self.points)
        intervals = self.envelope[['objective_lower_bound', 'objective_upper_bound']]
        max_distance = 0.
        max_interval = None
        for i, (lb, ub) in intervals.iterrows():
            distance = abs(ub - lb)
            if distance > max_distance:
                max_distance = distance
                max_interval = (lb, ub)
        step_size = (max_interval[1] - max_interval[0]) / (self.points - 1)
        print step_size
        grid = list()
        for i, row in self.envelope.iterrows():
            variables = row[self.variables]
            lb = row.objective_lower_bound
            ub = row.objective_upper_bound
            coordinate = lb
            while coordinate < ub:
                grid.append(list(variables.values) + [coordinate])
                coordinate += step_size
            grid.append(list(variables.values) + [ub])
        columns = self.variables + [self.target]
        self.grid = DataFrame(grid, columns=columns)

    def _set_bounds(self, point):
        for variable in self.variables:
            reaction = self.design_space_model.reactions.get_by_id(variable)
            bound = point[variable]
            reaction.lower_bound, reaction.upper_bound = bound, bound
        target_reaction = self.design_space_model.reactions.get_by_id(self.target)
        target_bound = point[self.target]
        target_reaction.lower_bound, target_reaction.upper_bound = target_bound, target_bound

    def run(self, view=None):
        if view is None:
            view = config.default_view
        else:
            view = view
        reference_flux_ranges = flux_variability_analysis(self.reference_model, view=view)
        self.reference_flux_ranges = reference_flux_ranges
        self._init_search_grid()
        fva_solutions = list()
        progress = ProgressBar(len(self.grid))
        progress.start()
        for i, point in self.grid.iterrows():
            progress.update(i + 1)
            # print "Processing point:", point
            self._set_bounds(point)
            fva_sol = flux_variability_analysis(self.design_space_model, view=view)
            fva_solutions.append(fva_sol)
        reference_intervals = reference_flux_ranges[['lower_bound', 'upper_bound']].values
        for sol in fva_solutions:
            intervals = sol[['lower_bound', 'upper_bound']].values
            gaps = [self._interval_gap(interval1, interval2) for interval1, interval2 in
                    izip(reference_intervals, intervals)]
            sol['gaps'] = gaps
        if self.normalize_fluxes_by is not None:
            for sol in fva_solutions:
                normalized_intervals = sol[['lower_bound', 'upper_bound']].values / sol.lower_bound[
                    self.normalize_fluxes_by]
                normalized_gaps = [self._interval_gap(interval1, interval2) for interval1, interval2 in
                                   izip(reference_intervals, normalized_intervals)]
                sol['normalized_gaps'] = normalized_gaps

        return fva_solutions


if __name__ == '__main__':
    from cameo.io import load_model

    model = load_model(
        '/Users/niko/Arbejder/Dev/cameo/tests/data/EcoliCore.xml')

    solution = model.solve()
    max_growth = solution.f

    reference_model = model.copy()
    biomass_rxn = model.reactions.get_by_id('Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2')
    reference_model.reactions.get_by_id(
        'Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2').lower_bound = max_growth

    diffFVA = DifferentialFVA(design_space_model=model,
                              reference_model=reference_model,
                              target='EX_succ_LPAREN_e_RPAREN_',
                              variables=['Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2',
                                         'EX_o2_LPAREN_e_RPAREN_'],
                              normalize_fluxes_by='Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2',
                              points=10
    )
    # result = diffFVA.run(view=SequentialView())
    result = diffFVA.run(view=MultiprocessingView())

