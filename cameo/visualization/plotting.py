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

from uuid import uuid1

from cameo import config, util

try:
    # import matplotlib.pyplot as plt

    def plot_production_envelope_ipython_matplotlib(envelope, key, highligt=[]):
        pass
        # plt.plot(envelope["objective_upper_bound"], envelope[key], title="Production envelop")
        # plt.xlabel("growth")
        # plt.ylabel(key)

except ImportError:

    def plot_production_envelope_ipython_matplotlib(envelope, key, highligt=[]):
        pass

try:
    from bokeh.plotting import *

    def plot_production_envelope_ipython_bokeh(envelope, key, highligt=[]):
        output_notebook(url=config.bokeh_url, docname=str(uuid1()))
        p = figure(title="Production envelope")
        p.xaxis.axis_label = 'growth'
        p.yaxis.axis_label = key
        colors = ["green" if i in highligt else "red" for i in envelope.index]
        if "label" in envelope.columns:
            p.text(envelope["label"], envelope["objective_upper_bound"], envelope[key])
        p.circle(envelope["objective_upper_bound"], envelope[key], color=colors, size=12)
        p.line(envelope["objective_upper_bound"], envelope[key], color="blue")
        show(p)

except ImportError:

    def plot_production_envelope_ipython_bokeh(envelope, key, highligt=[]):
        pass

try:
    from bashplotlib import scatterplot

    def plot_production_envelope_cli(envelope, key, highligt=[]):
        scatterplot.plot_scatter(None, envelope["objective_upper_bound"], envelope[key], "*")

except ImportError:
    def plot_production_envelope_cli(envelope, key, highligt=[]):
        pass


def plot_production_envelope(envelope, key, highligt=[]):
    if util.in_ipnb():
        if config.use_bokeh:
            plot_production_envelope_ipython_bokeh(envelope, key, highligt=highligt)
        elif config.use_matplotlib:
            plot_production_envelope_ipython_matplotlib(envelope, key, highligt=highligt)
    else:
        plot_production_envelope_cli(envelope, key, highligt=highligt)