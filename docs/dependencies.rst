Dependencies
============

Cameo has the following hard dependencies:


* [optlang](https://pypi.python.org/pypi/optlang) (for defining optimization problems)
* [numpy](http://www.numpy.org/) and [scipy](http://www.scipy.org/) (for obvious reasons)
* [inspyred](https://pypi.python.org/pypi/inspyred) (for heuristic optimizations)
* [escher](https://pypi.python.org/pypi/Escher) (for pathway visualizations)
* [blessings](https://pypi.python.org/pypi/blessings) and [IProgress](https://pypi.python.org/pypi/IProgress) (for displaying progress bars)
* [lazy-proxy-object](https://pypi.python.org/pypi/lazy-object-proxy) (for keeping models unevaluated)


Optionally, the following soft dependencies can be installed

* Jupyter notebook (cameo is tightly integrated with the notebook interface)
* bokeh (for plotting and showing progress)

The following dependencies are needed for development:

* [nose]() (for running unit tests)
* [rednose]() (for running unit tests)
* [coverage]() (for determining test coverage)

The following dependencies are needed for generating documentation:

* sphinx
* mock
* numpydoc
