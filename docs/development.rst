=================
Design principles
=================

The following documentation is intended for people that want to use cameo to develop new methods
and/or that would like to contribute to its development.


cameo vs. cobrapy
=================

Importing a model
-----------------

cobrapy (load a model in SBML format):

.. code-block:: python

    from cobra.io import read_sbml_model
    model = read_sbml_model('path/to/model.xml')

cameo (load models from different formats):

.. code-block:: python

    from cameo import load_model
    # read SBML model
    model = load_model('path/to/model.xml')
    # ... or read a pickled model
    model = load_model('path/to/model.pickle')
    # ... or just import a model by ID from http://darwin.di.uminho.pt/models
    iAF1260 = load_model('iAF1260')

Solving models
--------------

cobrapy:

.. code-block:: python

    solution = model.optimize()
    if solution.status == 'optimal':
        # proceed

.. code-block:: python

    try:
        solution = model.solve()
    except cameo.exceptions.SolverError:
        print "A non-optimal solution was returned by the solver"
    else:
        # proceed

It is important to note that cameo models maintain `optimize` to maintain
compatibility with cobrapy but we discourage its use.

    optlang
    copy_vs_time_machine


The optlang solver interface
============================

`optlang`_ is a generic interface solver interface to . It is based on the

.. _optlang: http://biosustain.github.io/optlang/

For efficiency reasons, cameo does not utilize the cobrapy's interface to LP and MILP solver. Instead it utilizes
optlang.


Changing the solver
-------------------

The LP/MILP solver can be changed in the following way.

.. code-block:: python

    model.solver = 'cplex'

Currently `'cplex'` and `'glpk' are supported.

cameo (load models from different formats):

.. code-block:: python

    from cameo import load_model
    # read SBML model
    model = load_model('path/to/model.xml')
    # ... or read a pickled model
    model = load_model('path/to/model.pickle')
    # ... or just import a model by ID from http://darwin.di.uminho.pt/models
    iAF1260 = load_model('iAF1260')


Manipulating the solver object
------------------------------

...


Avoiding copies
===============

Cameo users are encouraged to avoid making expensive copies of models (on the order seconds) and other data structures. Instead, we put forward a design pattern based on transactions (see ...)

.. code-block:: python

    from cameo.util import TimeMachine
    with TimeMachine() as tm:
        model.reactions.knock_out(time_machine=tm)
