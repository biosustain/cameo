cameo - computer aided metabolic engineering & optimization
-----------------------------------------------------------

|Documentation Status| |Build Status| |Coverage Status| |DOI|

Vision
~~~~~~

Cameo is a high-level python library developed to aid the *in silico*
strain design process in metabolic engineering projects. The library
provides a modular architecture that enables the efficient construction
of custom analysis workflows.

Dependencies
~~~~~~~~~~~~

This library depends on:

-  `cobrapy <https://github.com/opencobra/cobrapy>`__ for
   constraint-based modeling
-  `optlang <https://github.com/biosustain/optlang>`__ for heuristic
   optimization and mathematical programming

Furthermore, the following dependencies are needed:

-  `numpy <http://www.numpy.org/>`__ and
   `scipy <http://www.scipy.org/>`__ for obvious reasons.
-  `IPython <http://ipython.org/>`__ is needed for parallel computations
   and notebook interface.
-  `bokeh <http://bokeh.pydata.org/>`__ is needed for reporting progress
   and plotting in the IPython notebook interface.
-  `pandas <http://pandas.pydata.org/>`__ is needed because most
   functions returns results as pandas DataFrames.
-  `inspyred <https://pypi.python.org/pypi/inspyred>`__ for evolutionary
   computations.

Computationally heavy methods have been parallelized and can be run on a
clusters using the IPython parallelization framework (see example and
documentation for more details). The default fallback is python's
multiprocessing library.

Installation
~~~~~~~~~~~~

Run python setup.py install to install cameo. Installation is still a
little bit shaky, so if it fails due to version mismatches, try to
install the appropriate version manually and then retry
``python setup.py install``. For example: pip install pytz==2013b
--upgrade

.. |Documentation Status| image:: https://readthedocs.org/projects/cameo/badge/?version=devel
   :target: https://readthedocs.org/projects/cameo/?badge=devel
.. |Build Status| image:: https://travis-ci.org/biosustain/cameo.svg?branch=devel
   :target: https://travis-ci.org/biosustain/cameo
.. |Coverage Status| image:: https://coveralls.io/repos/biosustain/cameo/badge.svg?branch=devel
   :target: https://coveralls.io/r/biosustain/cameo?branch=devel
.. |DOI| image:: https://zenodo.org/badge/doi/10.5281/zenodo.19827.svg
   :target: http://dx.doi.org/10.5281/zenodo.19827
