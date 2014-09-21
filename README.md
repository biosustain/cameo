## cameo - computer assisted metabolic engineering & optimization

[![Build Status](https://travis-ci.org/biosustain/cameo.svg?branch=devel)](https://travis-ci.org/biosustain/cameo)
[![Coverage Status](https://coveralls.io/repos/biosustain/cameo/badge.png?branch=devel)](https://coveralls.io/r/biosustain/cameo?branch=devel)
[![Documentation Status](https://readthedocs.org/projects/cameo/badge/?version=latest)](https://readthedocs.org/projects/cameo/?badge=latest)

### Vision
Provide a high-level python library to aid _in silico_ strain design process in metabolic engineering projects. The library provides a modular architecture that enables the efficient construction of custom analysis workflows.

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

Computationally heavy methods have been parallelized and can be run on a clusters using the IPython parallelization framework (see example and documentation for more details). The default fallback is python's multiprocessing library.

### Installation
Run
    python setup.py install
to install cameo. Installation is still a little bit shaky, so if it fails due to version mismatches, try to install the appropriate version manually and then retry `python setup.py install`. For example:
pip install pytz==2013b --upgrade
