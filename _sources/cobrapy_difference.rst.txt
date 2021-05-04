=================
cameo vs. cobrapy
=================

-----------------
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

--------------
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
compatibility with cobrapy.
