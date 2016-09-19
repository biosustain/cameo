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

Alternatively you can use ``conda`` if you're already an anaconda user (there is no conda recipe for cameo though so you'll
still need to install it using ``pip``). Do the following to create a virtual environment and get some of the heavier dependencies out of the way.

.. code-block:: bash

    $ conda create -y -n cameo3.4 python=3.4 scipy numpy pandas numexpr matplotlib

Non-python dependencies
=======================

Cameo relies on optlang_ to solve optimization problems. Currently, optlang supports either glpk_ (open source) or cplex_
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

Installation
============

Cameo can be installed using ``pip`` (don't forget to activate your virtual environment in case you created one).

.. code-block:: bash

    $ pip install cameo


Soft dependencies
=================

The following soft dependencies can be installed all at once using ``pip install cameo[all]`` or individually
by specifying individual categories of dependencies (for example ``pip install cameo[swiglpk, sbml, ...]``).
The following categories are available::

    'docs': ['Sphinx>=1.3.5', 'numpydoc>=0.5'],
    'swiglpk': ['swiglpk>=1.2.14'],
    'plotly': ['plotly>=1.9.6'],
    'bokeh': ['bokeh>=0.11.1'],
    'jupyter': ['jupyter>=1.0.0', 'ipywidgets>=4.1.1'],
    'test': ['nose>=1.3.7', 'rednose>=0.4.3', 'coverage>=4.0.3'],
    'parallel': ['redis>=2.10.5', 'ipyparallel>=5.0.1'],
    'sbml': ['python-libsbml>=5.13.0', 'lxml>=3.6.0']


Development setup
=================

``pip`` can also be used to install cameo directly from the `github repository <https://github.com/biosustain/cameo>`_.

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
