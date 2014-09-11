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
from itertools import izip
import logging
from pandas import DataFrame, pandas
from progressbar import ProgressBar
from progressbar.widgets import ETA, Bar
from cameo import config, flux_variability_analysis
from cameo.parallel import SequentialView, MultiprocessingView
from cameo.solver_based_model import Reaction
from cameo.strain_design import StrainDesignMethod
from cameo.flux_analysis.analysis import phenotypic_phase_plane
from cameo.util import TimeMachine

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


class DifferentialFVA(StrainDesignMethod):
    """Differential flux variability analysis.

    Compares flux ranges of a reference model to a set of models that
    have been parameterized to lie on a grid of evenly spaced points in the
    n-dimensional production envelope (n being the number of reaction bounds
    to be varied).
    ::
        production
        ^
        |---------.          * reference_model
        | . . . . .\         . design_space_model
        | . . . . . \\
        | . . . . . .\\
        | . . . . . . \\
        o--------------*- >
                     growth

    Overexpression, downregulation, knockout, flux-reversal and other
    strain engineering targets can be inferred from the resulting comparison.

    Parameters
    ----------
    design_space_model : SolverBasedModel
        A model whose flux ranges will be scanned.
    reference_model : SolverBasedModel
        A model whose flux ranges represent the reference state and all calculates
        flux ranges will be compared to.
    objective : str or Reaction
        The reaction to be maximized.
    variables : iterable
        A iterable of n reactions (or IDs) to be scanned.
    exclude : iterable
        An iterable of reactions (or IDs) to be excluded in the analysis (exchange
        reactions will not be analyzed automatically).
    normalize_ranges_by : str or Reaction, optional
        A reaction ID that specifies a flux by whom all calculated flux ranges
        will be normalized by.
    points : int, optional
        Number of points to lay on the surface of the n-dimensional production envelope.
    """

    def __init__(self, design_space_model, reference_model, objective, variables=[],
                 exclude=[], normalize_ranges_by=None, points=20):
        super(DifferentialFVA, self).__init__()

        self.design_space_model = design_space_model
        self.reference_model = reference_model

        if isinstance(objective, Reaction):
            self.objective = objective.id
        else:
            self.objective = objective

        self.variables = list()
        for variable in variables:
            if isinstance(variable, Reaction):
                self.variables.append(variable.id)
            else:
                self.variables.append(variable)

        self.exclude = list()
        for elem in exclude:
            if isinstance(elem, Reaction):
                self.exclude.append(elem.id)
            else:
                self.exclude.append(elem)
        self.exclude += [reaction.id for reaction in design_space_model.exchanges]
        self.exclude += [reaction.id for reaction in reference_model.exchanges]
        self.exclude = set(self.exclude).difference(set([self.objective] + self.variables))

        self.points = points
        self.envelope = None
        self.grid = None

        if isinstance(normalize_ranges_by, Reaction):
            self.normalize_ranges_by = normalize_ranges_by.id
        else:
            self.normalize_ranges_by = normalize_ranges_by

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

    def _init_search_grid(self, surface_only=False, improvements_only=True):
        """Initialize the grid of points to be scanned within the production envelope."""
        self.envelope = phenotypic_phase_plane(
            self.design_space_model, self.variables, objective=self.objective, points=self.points)
        intervals = self.envelope[['objective_lower_bound', 'objective_upper_bound']]
        max_distance = 0.
        max_interval = None
        for i, (lb, ub) in intervals.iterrows():
            distance = abs(ub - lb)
            if distance > max_distance:
                max_distance = distance
                max_interval = (lb, ub)
        step_size = (max_interval[1] - max_interval[0]) / (self.points - 1)
        grid = list()
        minimal_reference_production = self.reference_flux_ranges['lower_bound'][self.objective]
        for i, row in self.envelope.iterrows():
            variables = row[self.variables]
            lb = row.objective_lower_bound
            if improvements_only:
                lb = max(lb, minimal_reference_production) + step_size
            ub = row.objective_upper_bound
            if not surface_only:
                coordinate = lb
                while coordinate < ub:
                    grid.append(list(variables.values) + [coordinate])
                    coordinate += step_size
            if improvements_only and ub <= minimal_reference_production:
                continue
            else:
                grid.append(list(variables.values) + [ub])
        columns = self.variables + [self.objective]
        self.grid = DataFrame(grid, columns=columns)

    def run(self, surface_only=False, improvements_only=True, view=None):
        """Run the differential flux variability analysis.

        Parameters
        ----------
        surface_only : bool, optional
            If only the optimal surface should be scanned.
        view : SequentialView or MultiprocessingView or ipython.cluster.DirectView, optional
            A parallelization view.
        surface_only : bool, optional
            If only the surface of the n-dimensional production envelope should be scanned.
        improvements_only : bool, optional
            If only grid points should should be scanned that constitute and improvement in production
            over the reference state.

        Returns
        -------
        pandas.Panel
            A pandas Panel containing a results DataFrame for every grid point scanned.
        """
        with TimeMachine() as tm:
            # Make sure that the design_space_model is initialized to its original state later
            for variable in self.variables:
                reaction = self.design_space_model.reactions.get_by_id(variable)
                tm(do=int, undo=partial(setattr, reaction, 'lower_bound', reaction.lower_bound))
                tm(do=int, undo=partial(setattr, reaction, 'upper_bound', reaction.upper_bound))
            target_reaction = self.design_space_model.reactions.get_by_id(self.objective)
            tm(do=int, undo=partial(setattr, target_reaction, 'lower_bound', reaction.lower_bound))
            tm(do=int, undo=partial(setattr, target_reaction, 'upper_bound', reaction.upper_bound))

            if view is None:
                view = config.default_view
            else:
                view = view

            included_reactions = [reaction.id for reaction in self.reference_model.reactions if
                                  not reaction.id in self.exclude]
            self.reference_flux_ranges = flux_variability_analysis(self.reference_model, reactions=included_reactions,
                                                                   view=view, remove_cycles=False)
            self._init_search_grid(surface_only=surface_only, improvements_only=improvements_only)

            progress = ProgressBar(len(self.grid), widgets=['Scanning grid points ', Bar(), ' ', ETA()])
            func_obj = _DifferentialFvaEvaluator(self.design_space_model, self.variables, self.objective,
                                                 included_reactions)
            results = list(progress(view.imap(func_obj, self.grid.iterrows())))

        solutions = dict((tuple(point.to_dict().items()), fva_result) for (point, fva_result) in results)
        reference_intervals = self.reference_flux_ranges[['lower_bound', 'upper_bound']].values
        for sol in solutions.itervalues():
            intervals = sol[['lower_bound', 'upper_bound']].values
            gaps = [self._interval_gap(interval1, interval2) for interval1, interval2 in
                    izip(reference_intervals, intervals)]
            sol['gaps'] = gaps
        if self.normalize_ranges_by is not None:
            for sol in solutions.itervalues():
                normalized_intervals = sol[['lower_bound', 'upper_bound']].values / sol.lower_bound[
                    self.normalize_ranges_by]
                normalized_gaps = [self._interval_gap(interval1, interval2) for interval1, interval2 in
                                   izip(reference_intervals, normalized_intervals)]
                sol['normalized_gaps'] = normalized_gaps
        for df in solutions.itervalues():
            ko_selection = df[(df.lower_bound == 0) &
                              (df.upper_bound == 0) &
                              (self.reference_flux_ranges.lower_bound != 0) &
                              self.reference_flux_ranges.upper_bound != 0]
            df['KO'] = False
            df['KO'][ko_selection.index] = True

        for df in solutions.itervalues():
            flux_reversal_selection = df[((self.reference_flux_ranges.upper_bound < 0) & (df.lower_bound > 0) |
                                          ((self.reference_flux_ranges.lower_bound > 0) & (df.upper_bound < 0)))]
            df['flux_reversal'] = False
            df['flux_reversal'][flux_reversal_selection.index] = True

        for df in solutions.itervalues():
            flux_reversal_selection = df[((self.reference_flux_ranges.lower_bound <= 0) & (df.lower_bound > 0)) | ((self.reference_flux_ranges.upper_bound >= 0) & (df.upper_bound <= 0))]
            df['suddenly_essential'] = False
            df['suddenly_essential'][flux_reversal_selection.index] = True

        # solutions['reference_flux_ranges'] = self.reference_flux_ranges
        return pandas.Panel(solutions)


class _DifferentialFvaEvaluator(object):
    def __init__(self, model, variables, objective, included_reactions):
        self.model = model
        self.variables = variables
        self.objective = objective
        self.included_reactions = included_reactions

    def __call__(self, point):
        self._set_bounds(point[1])
        return (point[1], flux_variability_analysis(self.model, reactions=self.included_reactions, remove_cycles=False,
                                                    view=SequentialView()))

    def _set_bounds(self, point):
        for variable in self.variables:
            reaction = self.model.reactions.get_by_id(variable)
            bound = point[variable]
            reaction.lower_bound, reaction.upper_bound = bound, bound
        target_reaction = self.model.reactions.get_by_id(self.objective)
        target_bound = point[self.objective]
        target_reaction.lower_bound, target_reaction.upper_bound = target_bound, target_bound


if __name__ == '__main__':
    from cameo.io import load_model
    from cameo.util import Timer

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
                              objective='EX_succ_LPAREN_e_RPAREN_',
                              variables=['Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2',
                                         'EX_o2_LPAREN_e_RPAREN_'],
                              normalize_ranges_by='Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2',
                              points=10
    )
    result = diffFVA.run(surface_only=True, view=SequentialView())
    print result.to_json()

    with Timer('Sequential'):
        result = diffFVA.run(surface_only=True, view=SequentialView())
    with Timer('Multiprocessing'):
        result = diffFVA.run(surface_only=True, view=MultiprocessingView())
        # try:
        # from IPython.parallel import Client
        #     client = Client()
        #     view = client.load_balanced_view()
        #     view.block = True
        # except:
        #     pass
        # else:
        #     with Timer('IPython'):
        #         result = diffFVA.run(surface_only=False, view=view)

        # model = load_model(
        #     '/Users/niko/Arbejder/Dev/cameo/tests/data/iJO1366.xml')
        #
        # reference_model = model.copy()
        # biomass_rxn = reference_model.reactions.get_by_id('Ec_biomass_iJO1366_core_53p95M')
        # biomass_rxn.lower_bound = .9 * reference_model.solve().f
        #
        #
        # diffFVA = DifferentialFVA(model, reference_model, 'EX_trp_DASH_L_LPAREN_e_RPAREN_', ['Ec_biomass_iJO1366_core_53p95M'],
        #                       normalize_ranges_by='Ec_biomass_iJO1366_core_53p95M', points=10)
        # with Timer('Sequential'):
        #     result = diffFVA.run(surface_only=True, view=SequentialView())
        # with Timer('Multiprocessing'):
        #     result = diffFVA.run(surface_only=True, view=MultiprocessingView())
        # try:
        #     from IPython.parallel import Client
        #     client = Client()
        #     view = client.load_balanced_view()
        #     with Timer('IPython'):
        #         result = diffFVA.run(surface_only=True, view=())
        # except:
        #     pass

