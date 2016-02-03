.. warning::
    These pages are under construction. Feel free to look around ...

Welcome to cameo!
=================

|Documentation Status| |Build Status| |Coverage Status| |DOI|

**Cameo** is a high-level python library developed to aid the strain
design process in metabolic engineering projects. The library provides a
modular framework of simulation methods, strain design methods, access
to models, that targets developers that want custom analysis workflows.

Computationally heavy methods have been parallelized and can be run on a
clusters using the IPython parallelization framework (see example and
documentation for more details). The default fallback is python's
multiprocessing library.

Furthermore, it exposes a high-level API to users that just want to compute
promising strain designs.

::

    from cameo.api import design
    design(product='L-Serine')

You got curious? Head over to `try.cameo.bio <http://try.cameo.bio>`__
and give it a try.

Table of Contents
-----------------

.. toctree::
    :maxdepth: 3

    dependencies
    installation
    01-quick-start
    02-import-models
    03-simulate-models
    04-analyze-models
    05-predict-gene-knockout-strategies
    06-predict-expression-modulation-targets
    07-predict-heterologous-pathways
    08-high-level-API
    09-vanillin-production
    11-multiprocess
    cobrapy_difference
    API

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. |Build Status| image:: https://travis-ci.org/biosustain/cameo.svg?branch=master
   :target: https://travis-ci.org/biosustain/cameo
.. |Coverage Status| .. image:: https://codecov.io/github/biosustain/cameo/coverage.svg?branch=devel
    :target: https://codecov.io/github/biosustain/cameo?branch=devel
.. |DOI| image:: https://zenodo.org/badge/doi/10.5281/zenodo.19827.svg
   :target: http://dx.doi.org/10.5281/zenodo.19827
