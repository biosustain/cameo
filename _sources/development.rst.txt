===========
Development
===========

The following documentation is intended for people that want to use cameo to develop new methods
and/or that would like to contribute to its development.

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
