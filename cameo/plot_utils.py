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
from uuid import uuid1

from cameo import config, util
from cameo.config import in_ipnb

try:
    # import matplotlib.pyplot as plt

    def plot_production_envelope_ipython_matplotlib(envelope, key, highligt=[]):
        pass
        # plt.plot(envelope["objective_upper_bound"], envelope[key], title="Production envelop")
        # plt.xlabel("growth")
        # plt.ylabel(key)

    def plot_dfba_performance_ipython_matplotlib(performance, ids):
        pass

    def plot_dfba_solution_ipython_matplotlib(dfba_solution):
        pass

except ImportError:

    def plot_production_envelope_ipython_matplotlib(envelope, key, highligt=[]):
        pass

    def plot_dfba_performance_ipython_matplotlib(performance, ids):
        pass

    def plot_dfba_solution_ipython_matplotlib(dfba_solution):
        pass

try:
    from bokeh.plotting import *
    from bokeh.charts import Bar
    from IPython.core.display import display_html

    def plot_production_envelope_ipython_bokeh(envelope, key, highligt=[]):
        output_notebook(url=config.bokeh_url, docname=str(uuid1()))
        p = figure(title="Production envelope")
        p.xaxis.axis_label = 'growth'
        p.yaxis.axis_label = key
        colors = ["green" if i in highligt else "red" for i in envelope.index]
        p.circle(envelope["objective_upper_bound"], envelope[key], color=colors, size=12)
        p.line(envelope["objective_upper_bound"], envelope[key], color="blue")
        show(p)

    def plot_dfba_performance_ipython_bokeh(performance, ids):
        output_notebook(url=config.bokeh_url, docname=str(uuid1()))
        bar = Bar(performance, ids, title="Performance", legend=True, stacked=True)
        display_html(bar.html)

    def plot_dfba_solution_ipython_bokeh(dfba_solution):
        output_notebook(url=config.bokeh_url, docname=str(uuid1()))
        p = figure(title="%s" % dfba_solution.reactor.id)
        p.xaxis.axis_label = "time"
        p.yaxis.axis_label = ""
        colors = util.generate_colors(len(dfba_solution))
        i = 0
        for metabolite in dfba_solution.metabolite_rates:
            p.line(dfba_solution.time, dfba_solution[metabolite],
                   color=colors[i],
                   legend="%s (mmol/L/h)" % metabolite,
                   size=12)
            i += 1
        for organism in dfba_solution.growth_rates:
            p.line(dfba_solution.time, dfba_solution[organism],
                   color=colors[i],
                   legend="%s (g/L/h)" % organism,
                   size=12)
            i += 1
        show(p)

except ImportError:

    def plot_production_envelope_ipython_bokeh(envelope, key, highligt=[]):
        pass

    def plot_dfba_performance_ipython_bokeh(performance, ids):
        pass

    def plot_dfba_solution_ipython_bokeh(dfba_solution):
        pass

try:
    from bashplotlib import scatterplot

    def plot_production_envelope_cli(envelope, key, highligt=[]):
        scatterplot.plot_scatter(None, envelope["objective_upper_bound"], envelope[key], "*")

    def plot_dfba_performance_cli(performance, ids):
        pass

    def plot_dfba_solution_cli(dfba_solution):
        pass

except ImportError:
    def plot_production_envelope_cli(envelope, key, highligt=[]):
        pass

    def plot_dfba_performance_cli(performance, ids):
        pass

    def plot_dfba_solution_cli(dfba_solution):
        pass


def plot_production_envelope(envelope, key, highligt=[]):
    if in_ipnb():
        if config.use_bokeh:
            plot_production_envelope_ipython_bokeh(envelope, key, highligt=highligt)
        elif config.use_matplotlib:
            plot_production_envelope_ipython_matplotlib(envelope, key, highligt=highligt)
    else:
        plot_production_envelope_cli(envelope, key, highligt=highligt)


def plot_dfba_performance(performance, ids):
    if in_ipnb():
        if config.use_bokeh:
            plot_dfba_performance_ipython_bokeh(performance, ids)
        elif config.use_matplotlib:
            plot_dfba_performance_ipython_matplotlib(performance, ids)
    else:
        plot_dfba_performance_cli(performance, ids)


def plot_dfba_solution(dfba_solution):
    if in_ipnb():
        if config.use_bokeh:
            plot_dfba_solution_ipython_bokeh(dfba_solution)
        elif config.use_matplotlib:
            plot_dfba_performance_ipython_matplotlib(dfba_solution)
    else:
        plot_dfba_solution_cli(dfba_solution)