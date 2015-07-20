.. warning::
    These pages are under construction. Feel free to look around ...

cameo
=====

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
    :maxdepth: 2

    installation
    1-quick-start
    2-import-models
    3-simulate-models
    4-analyze-models
    5-predict-expression-modulation-targets
    6-predict-gene-knockout-strategies
    7-predict-heterologous-pathways
    8-high-level-AP
    parallelization
    cobrapy_difference
    How_to
    API

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. |Documentation Status| image:: https://readthedocs.org/projects/cameo/badge/?version=devel
   :target: https://readthedocs.org/projects/cameo/?badge=devel
.. |Build Status| image:: https://travis-ci.org/biosustain/cameo.svg?branch=devel
   :target: https://travis-ci.org/biosustain/cameo
.. |Coverage Status| image:: https://coveralls.io/repos/biosustain/cameo/badge.svg?branch=devel
   :target: https://coveralls.io/r/biosustain/cameo?branch=devel
.. |DOI| image:: https://zenodo.org/badge/doi/10.5281/zenodo.19827.svg
   :target: http://dx.doi.org/10.5281/zenodo.19827
