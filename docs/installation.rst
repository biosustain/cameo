============
Installation
============

Setting up a virtual environment first
======================================

We highly recommended installing cameo inside a virtual environment (virtualenv_).
virtualenvwrapper_ tremendously simplifies using virtualenv_ and can easily
be installed using virtualenv-burrito_. Once you installed virtualenv_ and virtualenvwrapper_, run

.. code-block:: bash

    $ mkvirtualenv cameo  # or whatever you'd like to call your virtual environment
    $ workon cameo

and then continue with the installation instructions described below.

Alternatively you could use `conda` if you're already an anaconda user (there is no conda recipe for cameo though so you'll
still need install using `pip`).

.. code-block:: bash

    $ mkvirtualenv cameo  # or whatever you'd like to call your virtual environment
    $ workon cameo

Non-python dependencies
=======================

cameo relies on optlang_ to solve optimization problems. Currently, optlang supports either glpk_ (open source) or cplex_
(academic licenses available), which are not python tools. At least one of them has to be installed before one can proceed
with the cameo installation.

GLPK
----

Using cameo with glpk_ also requires swig_ to be installed (in order to generate python bindings).
On ubuntu (or other similar linux platforms) we recommend using :code:`apt-get`:

.. code-block:: bash

    $ sudo apt-get install libglpk-dev glpk-utils swig

On macs we recommend using homebrew_.

.. code-block:: bash

    $ brew install swig
    $ brew install glpk

CPLEX
-----

The cplex_ contains a python directory (similar to :code:`IBM/ILOG/CPLEX_Studio1251/cplex/python/x86-64_osx`). Inside
this directory run

.. code-block:: bash

    $ python setup.py install

to install the python bindings.

Hard dependencies
=================

Cameo has the following hard dependencies.

* [optlang](https://pypi.python.org/pypi/optlang) (for defining optimization problems)
* [numpy](http://www.numpy.org/) and [scipy](http://www.scipy.org/) (for obvious reasons)
* [inspyred](https://pypi.python.org/pypi/inspyred) (for heuristic optimizations)
* [escher](https://pypi.python.org/pypi/Escher) (for pathway visualizations)
* [blessings](https://pypi.python.org/pypi/blessings) and [IProgress](https://pypi.python.org/pypi/IProgress) (for displaying progress bars)
* [lazy-proxy-object](https://pypi.python.org/pypi/lazy-object-proxy) (for keeping models unevaluated)


Normal installation
===================

.. warning::
    cameo is still under heavy development. We recommend installing the development version (see below)
    if you would like to stay up-to-date with the latest changes.

cameo can be installed using `pip`.

.. code-block:: bash

    $ pip install cameo

This will also install the hard dependencies mentioned above.

Soft dependencies
=================

We highly recommend installing the following soft dependencies to provide a richer user experience when working with cameo.

- [Jupyter notebook](https://pypi.python.org/pypi/jupyter)) >= 1.0.0 (for interactive modeling environment)
- [bokeh](https://pypi.python.org/pypi/bokeh)) >= 0.11.0 (for plotting)
- [escher](https://pypi.python.org/pypi/escher)) >= 1.2.1 (for pathway visualizations)

All of these soft dependencies can be installed using pip.

.. code-block:: bash

    $ pip install jupyter bokeh escher

Furthermore, the following dependencies are needed for developing and contributing to cameo.

-  `sphinx`_) (for generating documentation)
-  `numpydoc`_ (for using numpy doc strings)
-  `nose`_ (for running unit tests)
-  `rednose`_ (for running unit tests)
-  `codecov`_ (for determining test coverage)


Development setup
=================

`pip` can also be used to install cameo directly from the `github repository <https://github.com/biosustain/cameo>`_.

.. code-block:: bash

    $ pip install -e git+https://github.com/biosustain/cameo.git@devel#egg=cameo

Alternatively, you can clone the repository (or your fork) and then run

.. code-block:: bash

    $ pip install -e .

within the cameo directory.

.. _homebrew: http://brew.sh/
.. _swig: http://www.swig.org/
.. _glpk: https://www.gnu.org/software/glpk/
.. _cplex: http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/
.. _optlang: https://github.com/biosustain/optlang
.. _virtualenv-burrito: https://github.com/brainsik/virtualenv-burrito
.. _virtualenv: https://pypi.python.org/pypi/virtualenv
.. _virtualenvwrapper: https://pypi.python.org/pypi/virtualenvwrapper

.. _sphinx: https://pypi.python.org/pypi/sphinx
.. _numpydoc: https://pypi.python.org/pypi/numpydoc
.. _nose: https://pypi.python.org/pypi/nose/
.. _rednose: https://pypi.python.org/pypi/rednose
.. _codecov: https://pypi.python.org/pypi/codecov
