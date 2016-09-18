=================
Welcome to cameo!
=================

|PyPI| |License| |Build Status| |Coverage Status| |DOI|

**Cameo** is a high-level python library developed to aid the strain
design process in metabolic engineering projects. The library provides a
modular framework of simulation methods, strain design methods, access
to models, that targets developers that want custom analysis workflows.

Computationally heavy methods have been parallelized and can be run on a
clusters using the IPython parallelization framework (see example and
documentation for more details). The default fallback is python's
multiprocessing library.

Furthermore, it will expose (in the near future) a high-level API to users that simply want to compute
promising strain designs (work in progress ...).

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


.. |PyPI| image:: https://img.shields.io/pypi/v/cameo.svg
   :target: https://pypi.python.org/pypi/cameo
.. |License| image:: http://img.shields.io/badge/license-APACHE2-blue.svg
   :target: http://img.shields.io/badge/license-APACHE2-blue.svg
.. |Build Status| image:: https://travis-ci.org/biosustain/cameo.svg?branch=master
   :target: https://travis-ci.org/biosustain/cameo
.. |Coverage Status| image:: https://coveralls.io/repos/biosustain/cameo/badge.svg?branch=devel
   :target: https://coveralls.io/r/biosustain/cameo?branch=devel
.. |DOI| image:: https://zenodo.org/badge/5031/biosustain/cameo.svg
   :target: https://zenodo.org/badge/latestdoi/5031/biosustain/cameo
