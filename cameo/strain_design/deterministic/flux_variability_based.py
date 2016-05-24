# -*- coding: utf-8 -*-
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


import os
import re
import six
import warnings
import logging

import numpy

from functools import partial
from uuid import uuid4


try:
    from IPython.core.display import display, HTML, Javascript
except ImportError:
    pass

from IProgress import ProgressBar
from pandas import DataFrame, pandas
from cameo.visualization.plotting import plotter

from cameo import config

from cameo.ui import notice
from cameo.util import TimeMachine, in_ipnb
from cameo.config import non_zero_flux_threshold
from cameo.parallel import SequentialView

from cameo.core.reaction import Reaction
from cameo.core.metabolite import Metabolite

from cameo.visualization.escher_ext import NotebookBuilder
from cameo.visualization.palette import mapper, Palette

from cameo.flux_analysis.analysis import flux_variability_analysis, phenotypic_phase_plane, PhenotypicPhasePlaneResult
from cameo.flux_analysis.simulation import pfba, fba

from cameo.strain_design.strain_design import StrainDesignMethod, StrainDesignResult, StrainDesign




with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from IPython.html.widgets import interact, IntSlider
    except ImportError:
        try:
            from ipywidgets import interact, IntSlider
        except ImportError:
            pass

zip = my_zip = six.moves.zip

__all__ = ['DifferentialFVA', 'FSEOF']

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
        | . . . . . \
        | . . . . . .\
        | . . . . . . \
        o--------------*- >
                     growth

    Overexpression, downregulation, knockout, flux-reversal and other
    strain engineering targets can be inferred from the resulting comparison.

    Parameters
    ----------
    design_space_model : SolverBasedModel
        A model whose flux ranges will be scanned.
    objective : str or Reaction or Metabolite
        A reaction whose flux or a metabolite whose production should be maximized.
    variables : iterable, optional
        A iterable of n reactions (or IDs) to be scanned (defaults to current objective in design_space_model).
    reference_model : SolverBasedModel, optional
        A model whose flux ranges represent the reference state and all calculated
        flux ranges will be compared to. Defaults to design_space_model constrained
        to its maximum objective value.
    exclude : iterable
        An iterable of reactions (or IDs) to be excluded in the analysis (exchange
        reactions will not be analyzed automatically).
    normalize_ranges_by : str or Reaction, optional
        A reaction ID that specifies a flux by whom all calculated flux ranges
        will be normalized by.
    points : int, optional
        Number of points to lay on the surface of the n-dimensional production envelope (defaults to 10).

    Examples
    --------
    >>> from cameo import models
    >>> from cameo.strain_design.deterministic import DifferentialFVA
    >>> model = models.bigg.e_coli_core
    >>> reference_model = model.copy()
    >>> reference_model.reactions.Biomass_Ecoli_core_w_GAM.lower_bound = reference_model.solve().objective_value
    >>> diffFVA = DifferentialFVA(design_space_model=model,
                          reference_model=reference_model,
                          objective=model.reactions.EX_succ_e,
                          variables=[model.reactions.Biomass_Ecoli_core_w_GAM],
                          normalize_ranges_by=model.reactions.Biomass_Ecoli_core_w_GAM,
                          points=10)
    >>> result = diffFVA.run(surface_only=True)
    >>> result.plot()
    """

    def __init__(self, design_space_model, objective, variables=None, reference_model=None,
                 exclude=(), normalize_ranges_by=None, points=10):
        super(DifferentialFVA, self).__init__()

        self.design_space_model = design_space_model
        if reference_model is None:
            self.reference_model = self.design_space_model.copy()
            self.reference_model.fix_objective_as_constraint()
        else:
            self.reference_model = reference_model

        if isinstance(objective, Reaction):
            self.objective = objective.id
        elif isinstance(objective, Metabolite):
            try:
                self.reference_model.add_demand(objective)
            except ValueError:
                pass
            try:
                self.objective = self.design_space_model.add_demand(objective).id
            except ValueError:
                self.objective = self.design_space_model.reactions.get_by_id("DM_" + objective.id).id
        elif isinstance(objective, str):
            self.objective = objective
        else:
            raise ValueError('You need to provide an objective as a Reaction, Metabolite or a reaction id')

        if variables is None:
            # try to establish the current objective reaction
            obj_var_ids = [variable.name for variable in self.design_space_model.objective.expression.free_symbols]
            obj_var_ids = [re.sub('_reverse.*', '', id) for id in obj_var_ids]
            if len(set(obj_var_ids)) != 1:
                raise ValueError(
                    "The current objective in design_space_model is not a single reaction objective. DifferentialFVA does not support composite objectives.")
            else:
                self.variables = [self.design_space_model.reactions.get_by_id(obj_var_ids[0]).id]
        else:
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
        self.exclude += [reaction.id for reaction in self.design_space_model.exchanges]
        self.exclude += [reaction.id for reaction in self.reference_model.exchanges]
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

    def run(self, surface_only=True, improvements_only=True, view=None):
        """Run the differential flux variability analysis.

        Parameters
        ----------
        surface_only : bool, optional
            If only the surface of the n-dimensional production envelope should be scanned (defaults to True).
        improvements_only : bool, optional
            If only grid points should should be scanned that constitute and improvement in production
            over the reference state (defaults to True).
        view : SequentialView or MultiprocessingView or ipython.cluster.DirectView, optional
            A parallelization view (defaults to SequentialView).

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
            tm(do=int, undo=partial(setattr, target_reaction, 'lower_bound', target_reaction.lower_bound))
            tm(do=int, undo=partial(setattr, target_reaction, 'upper_bound', target_reaction.upper_bound))

            if view is None:
                view = config.default_view
            else:
                view = view

            included_reactions = [reaction.id for reaction in self.reference_model.reactions if
                                  reaction.id not in self.exclude]
            # FIXME: reference flux ranges should be for a fraction (say 0.75) of the original objective
            self.reference_flux_ranges = flux_variability_analysis(self.reference_model, reactions=included_reactions,
                                                                   view=view, remove_cycles=False).data_frame
            self._init_search_grid(surface_only=surface_only, improvements_only=improvements_only)

            progress = ProgressBar(len(self.grid))
            func_obj = _DifferentialFvaEvaluator(self.design_space_model, self.variables, self.objective,
                                                 included_reactions)
            results = list(progress(view.imap(func_obj, self.grid.iterrows())))

        solutions = dict((tuple(point.iteritems()), fva_result.data_frame) for (point, fva_result) in results)
        reference_intervals = self.reference_flux_ranges[['lower_bound', 'upper_bound']].values
        for sol in six.itervalues(solutions):
            intervals = sol[['lower_bound', 'upper_bound']].values
            gaps = [self._interval_gap(interval1, interval2) for interval1, interval2 in
                    my_zip(reference_intervals, intervals)]
            sol['gaps'] = gaps
        if self.normalize_ranges_by is not None:
            for sol in six.itervalues(solutions):
                normalized_intervals = sol[['lower_bound', 'upper_bound']].values / sol.lower_bound[
                    self.normalize_ranges_by]
                normalized_gaps = [self._interval_gap(interval1, interval2) for interval1, interval2 in
                                   my_zip(reference_intervals, normalized_intervals)]
                sol['normalized_gaps'] = normalized_gaps
        chopped_ref_lower_bounds = self.reference_flux_ranges.lower_bound.apply(lambda x: 0 if abs(x) < non_zero_flux_threshold else x)
        chopped_ref_upper_bounds = self.reference_flux_ranges.upper_bound.apply(lambda x: 0 if abs(x) < non_zero_flux_threshold else x)
        for df in six.itervalues(solutions):
            df_chopped = df.applymap(lambda x: 0 if abs(x) < non_zero_flux_threshold else x)
            ko_selection = df[(df_chopped.lower_bound == 0) &
                              (df_chopped.upper_bound == 0) &
                              (chopped_ref_lower_bounds != 0) &
                              chopped_ref_upper_bounds != 0]
            df['KO'] = False
            df.loc[ko_selection.index, 'KO'] = True

        for df in six.itervalues(solutions):
            df_chopped = df.applymap(lambda x: 0 if abs(x) < non_zero_flux_threshold else x)
            flux_reversal_selection = df[((chopped_ref_upper_bounds < 0) & (df_chopped.lower_bound > 0) |
                                          ((chopped_ref_lower_bounds > 0) & (df_chopped.upper_bound < 0)))]
            df['flux_reversal'] = False
            df.loc[flux_reversal_selection.index, 'flux_reversal'] = True

        for df in six.itervalues(solutions):
            df_chopped = df.applymap(lambda x: 0 if abs(x) < non_zero_flux_threshold else x)
            suddenly_essential_selection = df[((df_chopped.lower_bound <= 0) & (df_chopped.lower_bound > 0)) | (
                (chopped_ref_upper_bounds >= 0) & (df_chopped.upper_bound <= 0))]
            df['suddenly_essential'] = False
            df.loc[suddenly_essential_selection.index, 'suddenly_essential'] = True

        if self.objective is None:
            objective = self.reference_model.objective
        else:
            objective = self.objective

        if isinstance(objective, Reaction):
            if hasattr(self.objective, 'nice_id'):
                nice_objective_id = objective.nice_id
                objective = objective.id
            else:
                objective = objective.id
                nice_objective_id = objective
        else:
            objective = str(self.objective)
            nice_objective_id = str(objective)

        return DifferentialFVAResult(pandas.Panel(solutions), self.envelope, self.reference_flux_ranges,
                                     self.variables, objective, nice_objective_id=nice_objective_id,
                                     nice_variable_ids=self.variables)


class DifferentialFVAResult(PhenotypicPhasePlaneResult):
    def __init__(self, solutions, phase_plane, reference_fva, variables_ids, objective, *args, **kwargs):
        if isinstance(phase_plane, PhenotypicPhasePlaneResult):
            phase_plane = phase_plane._phase_plane
        super(DifferentialFVAResult, self).__init__(phase_plane, variables_ids, objective, *args, **kwargs)
        self.reference_fva = reference_fva
        self.solutions = solutions

    def __getitem__(self, item):
        columns = ["lower_bound", "upper_bound", "gaps", "normalized_gaps", "KO", "flux_reversal", "suddenly_essential"]
        rows = list(range(len(self.solutions)))
        values = numpy.ndarray((len(rows), len(columns)))
        for i in rows:
            values[i] = self.solutions.iloc[i].loc[item].values

        data = DataFrame(values, index=rows, columns=columns)
        data["KO"] = data["KO"].values.astype(numpy.bool)
        data["flux_reversal"] = data["flux_reversal"].values.astype(numpy.bool)
        data["suddenly_essential"] = data["suddenly_essential"].values.astype(numpy.bool)
        return data

    def plot(self, index=None, variables=None, grid=None, width=None, height=None, title=None, palette=None, **kwargs):
        if len(self.variable_ids) > 1:
            notice("Multi-dimensional plotting is not supported")
            return
        if index is not None:
            self._plot_flux_variability_analysis(index, variables=variables, width=width, grid=grid, palette=palette)
        else:
            self._plot_production_envelope(title=title, grid=grid, width=width, height=height)

    def _plot_flux_variability_analysis(self, index, variables=None, title=None,
                                        width=None, height=None, palette=None, grid=None):
        if variables is None:
            variables = self.reference_fva.index[0:10]

        title = "Compare WT solution %i" % index if title is None else title

        wt_fva_res = self.reference_fva.loc[variables]
        strain_fva_res = self.solutions.iloc[index].loc[variables]
        dataframe = pandas.DataFrame(columns=["lb", "ub", "strain", "reaction"])
        for reaction_id, row in wt_fva_res.iterrows():
            _df = pandas.DataFrame([[row['lower_bound'], row['upper_bound'], "WT", reaction_id]],
                                   columns=dataframe.columns)
            dataframe = dataframe.append(_df)

        for reaction_id, row in strain_fva_res.iterrows():
            _df = pandas.DataFrame([[row['lower_bound'], row['upper_bound'], "Strain %i" % index, reaction_id]],
                                   columns=dataframe.columns)
            dataframe = dataframe.append(_df)

        plot = plotter.flux_variability_analysis(dataframe, grid=grid, width=width, height=height,
                                                 title=title, x_axis_label="Reactions", y_axis_label="Flux limits",
                                                 palette=palette)

        plotter.display(plot)

    def _plot_production_envelope(self, title=None, width=None, height=None, grid=None):
        title = "DifferentialFVA Result" if title is None else title
        x = [elem[0][1] for elem in list(self.solutions.items)]
        y = [elem[1][1] for elem in list(self.solutions.items)]
        colors = ["red" for _ in x]
        points = zip(x, y)
        super(DifferentialFVAResult, self).plot(title=title, grid=grid, width=width, heigth=height,
                                                points=points, points_colors=colors)

    def _repr_html_(self):
        def _data_frame(solution):
            notice("%s: %f" % self.solutions.axes[0][solution - 1][0])
            notice("%s: %f" % self.solutions.axes[0][solution - 1][1])
            display(self.solutions.iloc[solution - 1])

        interact(_data_frame, solution=[1, len(self.solutions)])
        return ""

    @property
    def data_frame(self):
        return self.solutions

    @property
    def normalized_gaps(self):
        return numpy.concatenate(self.solutions.iloc[:, :, 3].values).astype(float)

    def display_on_map(self, index=0, map_name=None, palette="RdYlBu", **kwargs):
        # TODO: hack escher to use iterative maps
        self._display_on_map_static(index, map_name, palette=palette, **kwargs)

    def plot_scale(self, palette="YlGnBu"):
        """
        Generates a color scale based on the flux distribution.
        It makes an array containing the absolute values and minus absolute values.

        The colors set as follows (p standsfor palette colors array):
        min   -2*std  -std      0      std      2*std   max
        |-------|-------|-------|-------|-------|-------|
        p[0] p[0] ..  p[1] .. p[2] ..  p[3] .. p[-1] p[-1]


        Arguments
        ---------
        palette: Palette, list, str
            A Palette from palettable of equivalent, a list of colors (size 5) or a palette name

        Returns
        -------
        tuple:
            ((-2*std, color), (-std, color) (0 color) (std, color) (2*std, color))

        """
        if isinstance(palette, str):
            palette = mapper.map_palette(palette, 5)
            palette = palette.hex_colors

        elif isinstance(palette, Palette):
            palette = palette.hex_colors

        values = [abs(v) for v in self.normalized_gaps if not (numpy.isnan(v) or numpy.isinf(v))]
        values += [-v for v in values]

        std = numpy.std(values)

        return (-2 * std, palette[0]), (-std, palette[1]), (0, palette[2]), (std, palette[3]), (2 * std, palette[4])

    def _display_on_map_static(self, index, map_name, palette="RdYlBu", **kwargs):
        try:
            import escher
            if os.path.exists(map_name):
                map_json = map_name
                map_name = None
            else:
                map_json = None

            values = self.normalized_gaps

            values = values[~numpy.isnan(values)]
            values = values[~numpy.isinf(values)]

            data = self.solutions.iloc[index]
            # Find values above decimal precision
            data = data[numpy.abs(data.normalized_gaps.astype(float)) > non_zero_flux_threshold]
            # Remove NaN rows
            data = data[~numpy.isnan(data.normalized_gaps.astype(float))]

            reaction_data = dict(data.normalized_gaps)
            for rid, gap in six.iteritems(reaction_data):
                if numpy.isposinf(gap):
                    gap = numpy.max(values)
                elif numpy.isneginf(gap):
                    gap = numpy.min(values)

                reaction_data[rid] = gap

            scale = self.plot_scale(palette)
            reaction_data['min'] = min(numpy.abs(values) * -1)
            reaction_data['max'] = max(numpy.abs(values))

            reaction_scale = [dict(type='min', color=scale[0][1], size=24),
                              dict(type='value', value=scale[0][0], color=scale[0][1], size=21),
                              dict(type='value', value=scale[1][0], color=scale[1][1], size=16),
                              dict(type='value', value=scale[2][0], color=scale[2][1], size=8),
                              dict(type='value', value=scale[3][0], color=scale[3][1], size=16),
                              dict(type='value', value=scale[4][0], color=scale[4][1], size=21),
                              dict(type='max', color=scale[4][1], size=24)]

            builder = escher.Builder(map_name=map_name, map_json=map_json, reaction_data=reaction_data,
                                     reaction_scale=reaction_scale)

            if in_ipnb():
                from IPython.display import display
                display(builder.display_in_notebook())
            else:
                builder.display_in_browser()

        except ImportError:
            print("Escher must be installed in order to visualize maps")

    def _display_on_map_iteractive(self, index, map_name, **kwargs):
        view = _MapView(self.solutions, map_name, **kwargs)
        slider = IntSlider(min=1, max=len(self.solutions), value=index + 1)
        slider.on_trait_change(lambda x: view(slider.get_state("value")["value"]))
        display(slider)
        view(1)


class _MapView(object):
    def __init__(self, solutions, map_name, **kwargs):
        self.solutions = solutions
        self.map_name = map_name
        self.builder = None
        self.kwargs_for_escher = kwargs

    def __call__(self, index):
        reaction_data = dict(self.solutions.iloc[index - 1].gaps)
        axis = self.solutions.axes[0]
        if self.builder is None:
            self._init_builder(reaction_data, axis[index - 1][0], axis[index - 1][1])
        else:
            self.builder.update(reaction_data)
            self.update_header(axis[index - 1][0], axis[index - 1][1])

    def update_header(self, objective, variable):
        display(Javascript("""
            jQuery("#objective-%s").text("%f");\n
            jQuery("#variable-%s").text("%f");
        """ % (self.header_id, objective[1], self.header_id, variable[1])))

    def _init_builder(self, reaction_data, objective, variable):
        self.header_id = str(uuid4()).replace("-", "_")
        display(HTML("""
        <p>
            %s&nbsp;<span id="objective-%s">%f</span></br>
            %s&nbsp;<span id="variable-%s">%f</span>
        </p>
        """ % (objective[0], self.header_id, objective[1],
               variable[0], self.header_id, variable[1])))

        self.builder = NotebookBuilder(map_name=self.map_name,
                                       reaction_data=reaction_data,
                                       reaction_scale=[
                                           dict(type='min', color="red", size=20),
                                           dict(type='median', color="grey", size=7),
                                           dict(type='max', color='green', size=20)], **self.kwargs_for_escher)
        display(self.builder.display_in_notebook())


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


class FSEOF(StrainDesignMethod):
    """
    Performs a Flux Scanning based on Enforced Objective Flux (FSEOF) analysis.

    Parameters
    ----------
    model : SolverBasedModel
    enforced_reaction : Reaction
        The flux that will be enforced. Reaction object or reaction id string.
    primary_objective : Reaction
        The primary objective flux (defaults to model.objective).

    References
    ----------
    .. [1] H. S. Choi, S. Y. Lee, T. Y. Kim, and H. M. Woo, 'In silico identification of gene amplification targets
    for improvement of lycopene production.,' Appl Environ Microbiol, vol. 76, no. 10, pp. 3097–3105, May 2010.

    """

    def __init__(self, model, primary_objective=None, *args, **kwargs):
        super(FSEOF, self).__init__(*args, **kwargs)
        self.model = model

        if primary_objective is None:
            self.primary_objective = model.objective
        elif isinstance(primary_objective, Reaction):
            if primary_objective in model.reactions:
                self.primary_objective = primary_objective
            else:
                raise ValueError("The reaction " + primary_objective.id + " does not belong to the model")
        elif isinstance(primary_objective, str):
            if primary_objective in model.reactions:
                self.primary_objective = model.reactions.get_by_id(primary_objective)
            else:
                raise ValueError("No reaction " + primary_objective + " found in the model")
        elif isinstance(primary_objective, type(model.objective)):
            self.primary_objective = primary_objective
        else:
            raise TypeError("Primary objective must be an Objective, Reaction or a string")

    def run(self, target=None, max_enforced_flux=0.9, number_of_results=10, exclude=(), simulation_method=fba,
            simulation_kwargs=None):
        """
        Performs a Flux Scanning based on Enforced Objective Flux (FSEOF) analysis.

        Parameters
        ----------
        target: str, Reaction, Metabolite
            The target for optimization.
        max_enforced_flux : float, optional
            The maximal flux of secondary_objective that will be enforced, relative to the theoretical maximum (defaults to 0.9).
        number_of_results : int, optional
            The number of enforced flux levels (defaults to 10).
        exclude : Iterable of reactions or reaction ids that will not be included in the output.

        Returns
        -------
        FseofResult
            An object containing the identified reactions and the used parameters.

        References
        ----------
        .. [1] H. S. Choi, S. Y. Lee, T. Y. Kim, and H. M. Woo, 'In silico identification of gene amplification targets
        for improvement of lycopene production.,' Appl Environ Microbiol, vol. 76, no. 10, pp. 3097–3105, May 2010.

        """
        model = self.model
        target = model.reaction_for(target)

        simulation_kwargs = simulation_kwargs if simulation_kwargs is not None else {}
        simulation_kwargs['objective'] = self.primary_objective

        if 'reference' not in simulation_kwargs:
            reference = simulation_kwargs['reference'] = pfba(model, **simulation_kwargs)
        else:
            reference = simulation_kwargs['reference']

        ndecimals = config.ndecimals

        # Exclude list
        exclude = list(exclude) + model.exchanges
        exclude_ids = [target.id]
        for reaction in exclude:
            if isinstance(reaction, Reaction):
                exclude_ids.append(reaction.id)
            else:
                exclude_ids.append(reaction)

        with TimeMachine() as tm:

            tm(do=int, undo=partial(setattr, model, "objective", model.objective))
            tm(do=int, undo=partial(setattr, target, "lower_bound", target.lower_bound))
            tm(do=int, undo=partial(setattr, target, "upper_bound", target.upper_bound))

            # Find initial flux of enforced reaction
            initial_fluxes = reference.fluxes
            initial_flux = round(initial_fluxes[target.id], ndecimals)

            # Find theoretical maximum of enforced reaction
            max_theoretical_flux = round(fba(model, objective=target.id, reactions=[target.id]).fluxes[target.id],
                                         ndecimals)

            max_flux = max_theoretical_flux * max_enforced_flux

            # Calculate enforcement levels
            levels = [initial_flux + (i + 1) * (max_flux - initial_flux) / number_of_results for i in
                      range(number_of_results)]

            # FSEOF results
            results = {reaction.id: [] for reaction in model.reactions}

            for level in levels:
                target.lower_bound = level
                target.upper_bound = level
                solution = simulation_method(model, **simulation_kwargs)
                for reaction_id, flux in solution.fluxes.items():
                    results[reaction_id].append(round(flux, ndecimals))

        # Test each reaction
        fseof_reactions = []
        for reaction_id, fluxes in results.items():
            if reaction_id not in exclude_ids \
                    and max(abs(max(fluxes)), abs(min(fluxes))) > abs(reference[reaction_id]) \
                    and min(fluxes) * max(fluxes) >= 0:
                fseof_reactions.append(model.reactions.get_by_id(reaction_id))

        results = {rea.id: results[rea.id] for rea in fseof_reactions}
        run_args = dict(max_enforced_flux=max_enforced_flux,
                        number_of_results=number_of_results,
                        solution_method=simulation_method,
                        simulation_kwargs=simulation_kwargs,
                        exclude=exclude)

        return FSEOFResult(fseof_reactions, target, model, self.primary_objective, levels, results, run_args, reference)


class FSEOFResult(StrainDesignResult):
    """
    Object for storing a FSEOF result.

    Attributes:
    -----------
    reactions: list
        A list of the reactions that are found to increase with product formation.
    enforced_levels: list
        A list of the fluxes that the enforced reaction was constrained to.
    data_frame: DataFrame
        A pandas DataFrame containing the fluxes for every reaction for each enforced flux.
    run_args: dict
        The arguments that the analysis was run with. To repeat do 'FSEOF.run(**FSEOFResult.run_args)'.

    """

    __method_name__ = "FSEOF"

    def plot(self, grid=None, width=None, height=None, title=None):
        pass

    def __init__(self, reactions, target, model, primary_objective, enforced_levels, reaction_results,
                 run_args, reference, *args, **kwargs):
        super(FSEOFResult, self).__init__(*args, **kwargs)
        self._reactions = reactions
        self._target = target
        self._model = model
        self._primary_objective = primary_objective
        self._run_args = run_args
        self._enforced_levels = enforced_levels
        self._reaction_results = reaction_results
        self._reference_fluxes = {r: reference.fluxes[r.id] for r in reactions}

    def __len__(self):
        return len(self.reactions)

    def __iter__(self):
        ref_fluxes = self._reference_fluxes
        for i, level in enumerate(self.enforced_levels):
            knockouts = [r for r, v in six.iteritems(self._reaction_results) if ref_fluxes[r.id] > 0 and v[i] == 0]
            over_expression = {r: v for r, v in six.iteritems(self._reaction_results) if v[i] > ref_fluxes[r.id]}
            down_regulation = {r: v for r, v in six.iteritems(self._reaction_results) if v[i] < ref_fluxes[r.id]}
            yield StrainDesign(knockouts=knockouts, over_expression=over_expression,
                               down_regulation=down_regulation, manipulation_type="reactions")

    def __eq__(self, other):
        return isinstance(other, self.__class__) and self.target == other.target and self.reactions == other.reactions

    @property
    def reactions(self):
        return self._reactions

    @property
    def model(self):
        return self._model

    @property
    def target(self):
        return self._target

    @property
    def primary_objective(self):
        return self._primary_objective

    @property
    def run_args(self):
        return self._run_args

    @property
    def enforced_levels(self):
        return self._enforced_levels

    def _repr_html_(self):
        template = """
<strong>Model:</strong> %(model)s</br>
<strong>Enforced objective:</strong> %(objective)s</br>
<strong>Primary objective:</strong> %(primary)s</br>
<br>
<strong>Reaction fluxes</strong><br><br>
%(df)s
"""
        return template % {'objective': self.target.id,
                           'reactions': "<br>".join(reaction.id for reaction in self.reactions),
                           'model': self.model.id,
                           'primary': str(self._primary_objective),
                           'df': self.data_frame._repr_html_()}

    @property
    def data_frame(self):
        df = pandas.DataFrame(self._reaction_results).transpose()
        df.columns = (i + 1 for i in range(len(self._enforced_levels)))
        df.loc[self.target.id] = self._enforced_levels
        return df

        # if __name__ == '__main__':
        #     from cameo.io import load_model
        #     from cameo.util import Timer
        #
        #     model = load_model(
        #         '/Users/niko/Arbejder/Dev/cameo/tests/data/EcoliCore.xml')
        #
        #     solution = model.solve()
        #     max_growth = solution.f
        #
        #     reference_model = model.copy()
        #     biomass_rxn = model.reactions.get_by_id('Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2')
        #     reference_model.reactions.get_by_id(
        #         'Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2').lower_bound = max_growth
        #
        #     diffFVA = DifferentialFVA(design_space_model=model,
        #                               reference_model=reference_model,
        #                               objective='EX_succ_LPAREN_e_RPAREN_',
        #                               variables=['Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2',
        #                                          'EX_o2_LPAREN_e_RPAREN_'],
        #                               normalize_ranges_by='Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2',
        #                               points=10
        #                               )
        #     result = diffFVA.run(surface_only=True, view=SequentialView())
        #
        #     with Timer('Sequential'):
        #         result = diffFVA.run(surface_only=True, view=SequentialView())
        #     with Timer('Multiprocessing'):
        #         result = diffFVA.run(surface_only=True, view=MultiprocessingView())
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
