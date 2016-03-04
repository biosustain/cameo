=================
Design principles
=================

The following documentation is intended for people that want to use cameo to develop new methods
and/or that would like to contribute to its development.


cameo vs. cobrapy
=================

The following provides a comprehensive side-by-side comparison of cobrapy and cameo aiming to make it easier for users
who are already familiar with cobrapy to get started with cameo. While cameo uses and extends the same data structures
as provided by cobrapy (`Reaction`, `SolverBasedModel`, etc.) and is thus backwards-compatible to it, it deviates

Solver interface
----------------

Cameo deviates from cobrapy in the way optimization problems are solved by using a separate solver interface provided by
the optlang package (see also below). The following benefits ....

* Methods that require solving multiple succession will run a lot faster since previously found solution will be reused.
* Implementation of novel or published becomes a lot easier since optlang (based on the very popular symbolic math library
sympy) facilitates the formulation of constraints and objectives using equations (similar to GAMS) instead of matrix formalism.
* Adding of additional constraints (even non-metabolic ), is straight forwards and eliminates the problem in cobrapy of having
to define opti (check out the ice cream sandwich ...)
* The optimization problem is always accessible and has a one-to-one correspondence to the model.

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

Furthermore, cameo provides direct access to a number of other models through `cameo.models`.

.. code-block:: python

    from cameo import models
    models.e_coli_core.solve().f

Solving models
--------------

When models are optimized with cobrapy, `model.optimize()` returns the status code of the solver and leaves it to the
user to determine if the problem was successfully solved.

.. code-block:: python

    solution = model.optimize()
    if solution.status == 'optimal':
        # proceed

In our personal opinion, we believe that the more pythonic way is to raise an Exception if the problem could not be solved.

.. code-block:: python

    try:
        solution = model.solve()
    except cameo.exceptions.SolverError:
        print "A non-optimal solution was returned by the solver"
    else:
        # proceed

It is important to note that cameo models still provide the `optimize` method to maintain backwards
compatibility with cobrapy but we discourage its use.

    optlang
    copy_vs_time_machine

Convenience functions
---------------------

Cameo implements a number of convenience functions that are (currently) not available in cobrapy. For example, instead of
running flux variability analysis, one can quickly obtain the effective lower and upper bound

.. code-block:: python

    model.reaction.PGK.effective_lower_bound


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
