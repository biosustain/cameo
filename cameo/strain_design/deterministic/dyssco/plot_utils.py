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
from cameo import config

try:
    import matplotlib.pyplot as plt

    def plot_envelope_ipython_matplotlib(envelope, key):
        plt.plot(envelope["objective_upper_bound"], envelope[key], title="Production envelop")
        plt.xlabel("growth")
        plt.ylabel(key)

except ImportError:

    def plot_envelope_ipython_matplotlib(envelope, key):
        pass

try:
    from bokeh.plotting import *

    def plot_envelope_ipython_bokeh(envelope, key):
        output_notebook(url=config.bokeh_url, docname=str(uuid1()))
        hold()
        p = figure(title="Production envelope")
        p.xaxis.axis_label = 'growth'
        p.yaxis.axis_label = key
        p.circle(envelope["objective_upper_bound"], envelope[key], color=["blue"])
        show(p)

except ImportError:

    def plot_envelope_ipython_bokeh(envelope, key):
        pass


try:
    from bashplotlib import scatterplot

    def plot_envelope_cli(envelope, key):
        scatterplot.plot_scatter(None, envelope["objective_upper_bound"], envelope[key], "*")

except ImportError:
    def plot_envelope_cli(envelope, key):
        pass