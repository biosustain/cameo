## Cameo—Computer Aided Metabolic Engineering and Optimization

[![PyPI](https://img.shields.io/pypi/v/cameo.svg)](https://pypi.python.org/pypi/cameo)
[![License](http://img.shields.io/badge/license-APACHE2-blue.svg)](http://img.shields.io/badge/license-APACHE2-blue.svg)
[![Build Status](https://travis-ci.org/biosustain/cameo.svg?branch=master)](https://travis-ci.org/biosustain/cameo)
[![Coverage Status](https://coveralls.io/repos/biosustain/cameo/badge.svg?branch=devel)](https://coveralls.io/r/biosustain/cameo?branch=devel)
[![DOI](https://zenodo.org/badge/5031/biosustain/cameo.svg)](https://zenodo.org/badge/latestdoi/5031/biosustain/cameo)



### What is Cameo?
**Cameo** is a high-level python library developed to aid the strain design process in metabolic engineering projects. The library provides a modular framework of simulation methods, strain design methods, access to models, that targets developers that want  custom analysis workflows. 

Computationally heavy methods have been parallelized and can be run on a clusters using the IPython parallelization framework (see example and documentation for more details). The default fallback is python's multiprocessing library.

Furthermore, it exposes a high-level API to users that just want to compute promising strain designs. 

You got curious? Head over to [try.cameo.bio](http://try.cameo.bio) and give it a try.

### Installation
Use pip to install Cameo from [PyPI](https://pypi.python.org/pypi/cameo) (we recommend doing this inside a [virtual environment](http://docs.python-guide.org/en/latest/dev/virtualenvs/)).

    pip install cameo

We highly recommend updating `pip` beforehand (`pip install pip --upgrade`).

In case you downloaded the source code, run

	pip install -e .  # recommended

while you are in the top level directory. You might need to run these commands with administrative privileges if you're not using a virtual environment (using `sudo` for example).


### Examples

A number of examples are available as static ([nbviewer.ipython.org](http://nbviewer.ipython.org/github/biosustain/cameo-notebooks/tree/master/)) or executable Jupyter (née IPython) notebooks ([try.cameo.bio](http://try.cameo.bio)).

#### High-level API (for users)
Compute strain engineering strategies for a desired product in a number of host organisms using the high-level interface.

	from cameo.api import design
	design(product='L-Serine')

[Output](http://nbviewer.ipython.org/github/biosustain/cameo-notebooks/blob/master/8-high-level-API.ipynb)

#### Low-level API (for developers)

Find gene knockout targets using evolutionary computation.

	from cameo import models
	from cameo.strain_design.heuristic import GeneKnockoutOptimization
	from cameo.strain_design.heuristic.objective_functions import biomass_product_coupled_yield
	
	model = models.bigg.e_coli_core
	obj = biomass_product_coupled_yield(
	    model.reactions.Biomass_Ecoli_core_w_GAM,
	    model.reactions.EX_succ_e,
	    model.reactions.EX_glc_e)
	ko = GeneKnockoutOptimization(model=model, objective_function=obj)
	ko.run(max_evaluations=50000, n=1, mutation_rate=0.15, indel_rate=0.185)

[Output](http://nbviewer.ipython.org/github/biosustain/cameo-notebooks/blob/master/6-predict-gene-knockout-strategies.ipynb)

Predict heterologous pathways for a desired chemical.

	from cameo.strain_design import pathway_prediction
	predictor = pathway_prediction.PathwayPredictor(model)
	pathways = predictor.run(product="vanillin")

[Output](http://nbviewer.ipython.org/github/biosustain/cameo-notebooks/blob/master/7-predict-heterologous-pathways.ipynb)


### Dependencies
This library depends on:

- [cobrapy](https://github.com/opencobra/cobrapy) for constraint-based modeling
- [optlang](https://github.com/biosustain/optlang) for heuristic optimization and mathematical programming

Furthermore, the following dependencies are needed: 

- [numpy](http://www.numpy.org/) and [scipy](http://www.scipy.org/) for obvious reasons.
- [IPython](http://ipython.org/) is needed for parallel computations and notebook interface.
- [bokeh](http://bokeh.pydata.org/) is needed for reporting progress and plotting in the IPython notebook interface.
- [pandas](http://pandas.pydata.org/) is needed because most functions returns results as pandas DataFrames.
- [inspyred](https://pypi.python.org/pypi/inspyred) for evolutionary computations.

