=================
Welcome to cameo!
=================

|Build Status| |Coverage Status| |DOI|

**Cameo** is a high-level python library developed to aid the strain
design process in metabolic engineering projects. The library provides a
modular framework of simulation methods, strain design methods, access
to models, that targets developers that want custom analysis workflows.

Computationally heavy methods have been parallelized and can be run on a
clusters using the IPython parallelization framework (see example and
documentation for more details). The default fallback is python's
multiprocessing library.

Furthermore, it exposes a high-level API to users that simply want to compute
promising strain designs.

::

    from cameo.api import design
    design(product='L-Serine')

You got curious? Head over to `try.cameo.bio <http://try.cameo.bio>`__
and give it a try.

User's guide
============

.. toctree::
    :maxdepth: 2

    installation
    FAQ <FAQ>
    tutorials

Developers's guide
==================

.. toctree::
    :maxdepth: 2

    development
    contributing

API
===

.. toctree::
    :maxdepth: 2

    API


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. |Build Status| image:: https://travis-ci.org/biosustain/cameo.svg?branch=master
   :target: https://travis-ci.org/biosustain/cameo
.. |Coverage Status| image:: https://codecov.io/github/biosustain/cameo/coverage.svg?branch=devel
   :target: https://codecov.io/github/biosustain/cameo?branch=devel
.. |DOI| image:: https://zenodo.org/badge/doi/10.5281/zenodo.19827.svg
   :target: http://dx.doi.org/10.5281/zenodo.19827
