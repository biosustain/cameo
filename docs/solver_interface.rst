============
Dependencies
============

For efficiency reasons, cameo does not utilize the cobrapy's interface to LP and MILP solver. Instead it utilizes
optlang.


Changing the solver
==================

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


~~~~~~~
