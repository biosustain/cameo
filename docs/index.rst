.. warning::
    These pages are under construction. Feel free to look around ...

cameo
=====

|BuildStatus| |CoverageStatus|

**cameo** is a python library for computer-aided metabolic engineering and optimization
of microbes. It caters to different audiences as it provides a high-level interface
that allows computer savvy bench biologists to predict heterologous pathways and enumerate engineering strategies
and a modular framework for computational modelers that facilitates methods development and the construction
of custom analysis workflows and prediction tools.

cameo is based on `cobrapy <http://opencobra.github.io/cobrapy/>`_, the community standard for
constraint-based modeling of genome-scale metabolic models, and uses and extends the same data structures. For
efficiency reasons, however, it differs quite significantly from cobrapy in the way optimization
problems are solved and defined. Nevertheless, cameo models remain 100% compatible with cobrapy.
If you are already a cobrapy user you might want to start reading :doc:`here <cobrapy_difference>` to get an overview of the main differences.


Table of Contents
-----------------

.. toctree::
    :maxdepth: 2

    installation
    quickstart
    cameo_high_level_interface
    cameo_low_level_interface
    parallelization
    cobrapy_difference
    How_to
    API

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. |BuildStatus| image:: https://travis-ci.org/biosustain/cameo.svg?branch=devel
    :target: https://travis-ci.org/biosustain/cameo
.. |CoverageStatus| image:: https://coveralls.io/repos/biosustain/cameo/badge.png?branch=devel
    :target: https://coveralls.io/r/biosustain/cameo?branch=devel