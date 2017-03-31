===========
Development
===========

The following documentation is intended for people that want to use cameo to develop new methods
and/or that would like to contribute to its development.


cameo vs. cobrapy
~~~~~~~~~~~~~~~~~

While cameo uses and extends the same data structures as cobrapy, there exist a few notable differences.
The following provides a comprehensive side-by-side comparison of cobrapy and cameo aiming to make it easier for users
who are already familiar with cobrapy to get started with cameo.

Solver interface
----------------

Cameo deviates from cobrapy in the way optimization problems are solved by using a separate solver interface provided by
the optlang package (see :ref:`optlang_interface`), which has the following benefits:

* Methods that require solving a model multiple times will run faster since previously found solutions will be
  automatically re-used by the solvers to warm-start the next optimization.
* Implementation of novel or published methods becomes easier since optlang (based on the popular symbolic math
  library sympy) facilitates the formulation of constraints and objectives using equations (similar to GAMS)
  instead of matrix formalism.
* Adding additional constraints (even non-metabolic) is straight forwards and eliminates the problem in cobrapy
  of having to define dummy metabolites etc.

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

In order to avoid users accidentally working with non-optimal solutions, cameo will raise an exception instead.

.. code-block:: python

    try:
        solution = model.solve()
    except cameo.exceptions.SolverError:
        print "A non-optimal solution was returned by the solver"
    else:
        # proceed

It is important to note that cameo models still provide the `~cameo.core.SolverBasedModel.optimize` method to maintain backwards
compatibility with cobrapy but we discourage its use.

    optlang
    copy_vs_time_machine


The optlang solver interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For efficiency reasons, cameo does not utilize cobrapy's interfaces to LP and MILP solvers.
Instead it utilizes optlang_, which is a generic interface to a number of free and commercial optimization solvers.
It is based on the popular symbolic math library sympy_ and thus enables the formulation of optimization problems
using equations instead of matrix formalism.

Changing the solver
-------------------

The LP/MILP solver can be changed in the following way.

.. code-block:: python

    model.solver = 'cplex'

Currently `cplex`, `glpk`, and `gurobi` are supported.

Manipulating the solver object
------------------------------

The solver object in cameo is always accessible through `~SolverBased.solver`.
For example, one can inspect the optimization problem in CPLEX LP format by printing the solver object.

.. code-block:: python

    print(model.solver)

Having access to the `optlang`_ solver object provides for a very convenient way for manipulating the optimization problem.
For example, it is straightforward to add additional constraints, for example, a flux ratio constraint.

.. code-block:: python

    reaction1 = model.reactions.PGI
    reaction2 = model.reactions.G6PDH2r
    ratio = 5
    flux_ratio_constraint = model.solver.interface.Constraint(
        reaction1.flux_expression - ratio * reaction2.flux_expression,
        lb=0,
        ub=0)
    model.solver.add(flux_ratio_constraint)

This will constrain the flux split between glycolysis and pentose phosphate patwhay to 20.
`model.solver.interface` hereby provides access to


Good coding practices
=====================

Cameo developers and users are encouraged to avoid making copies of models and other data structures. Instead, we put
 forward a design pattern based on transactions.

.. code-block:: python

    from cameo.util import TimeMachine
    with TimeMachine() as tm:
        model.reactions.knock_out(time_machine=tm)


.. _optlang: http://biosustain.github.io/optlang/
.. _sympy: http://www.sympy.org/en/index.html

