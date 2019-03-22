Cameo—Computer Aided Metabolic Engineering and Optimization
-----------------------------------------------------------

.. summary-start

|Join the chat at https://gitter.im/biosustain/cameo| |PyPI| |License|
|Build Status| |Coverage Status| |DOI| |zenhub|

What is cameo?
~~~~~~~~~~~~~~

**Cameo** is a high-level python library developed to aid the strain
design process in metabolic engineering projects. The library provides a
modular framework of simulation and strain design methods that targets
developers that want to develop new design algorithms and custom analysis workflows.
Furthermore, it exposes a high-level API to users that just want to
compute promising strain designs.

Curious? Head over to `try.cameo.bio <http://try.cameo.bio>`__
and give it a try.

Please cite https://doi.org/10.1021/acssynbio.7b00423 if you've used cameo in a scientific publication.

.. summary-end

Installation
~~~~~~~~~~~~

.. installation-start

Use pip to install cameo from `PyPI <https://pypi.python.org/pypi/cameo>`__.

::

    $ pip install cameo


In case you downloaded or cloned the source code from `GitHub <https://github.com/biosustain/cameo>`__
or your own fork, you can run the following to install cameo for development.

::

    $ pip install -e <path-to-cameo-repo>  # recommended


You might need to run these commands with administrative
privileges if you're not using a virtual environment (using ``sudo`` for example).
Please check the `documentation <http://cameo.bio/installation.html>`__
for further details.

.. installation-end

Documentation and Examples
~~~~~~~~~~~~~~~~~~~~~~~~~~

Documentation is available on `cameo.bio <http://cameo.bio>`__. Numerous `Jupyter notebooks <http://nbviewer.ipython.org/github/biosustain/cameo-notebooks/tree/master/>`__
provide examples and tutorials and also form part of the documentation. They are also availabe in executable form on (`try.cameo.bio <http://try.cameo.bio>`__).
Furthermore, course materials for a two day cell factory engineering course are available `here <https://biosustain.github.io/cell-factory-design-course/>`__.

.. showcase-start

High-level API (for users)
^^^^^^^^^^^^^^^^^^^^^^^^^^

Compute strain engineering strategies for a desired product in a number
of host organisms using the high-level interface (runtime is on the order of hours).

::

    from cameo.api import design
    design(product='L-Serine')

`Output <http://nbviewer.ipython.org/github/biosustain/cameo-notebooks/blob/master/08-high-level-API.ipynb>`__


The high-level API can also be called from the command line.

::

    $ cameo design vanillin

For more information run

::

    $ cameo --help

Low-level API (for developers)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Find gene knockout targets using evolutionary computation.

::

    from cameo import models
    from cameo.strain_design.heuristic import GeneKnockoutOptimization
    from cameo.strain_design.heuristic.objective_functions import biomass_product_coupled_yield

    model = models.bigg.e_coli_core
    obj = biomass_product_coupled_yield(
        model.reactions.Biomass_Ecoli_core_w_GAM,
        model.reactions.EX_succ_e,
        model.reactions.EX_glc_e)
    ko = GeneKnockoutOptimization(model=model, objective_function=obj)
    ko.run(max_evaluations=50000, n=1, mutation_rate=0.15, indel_rate=0.185)

`Output <http://nbviewer.ipython.org/github/biosustain/cameo-notebooks/blob/master/05-predict-gene-knockout-strategies.ipynb>`__

Predict heterologous pathways for a desired chemical.

::

    from cameo.strain_design import pathway_prediction
    predictor = pathway_prediction.PathwayPredictor(model)
    pathways = predictor.run(product="vanillin")

`Output <http://nbviewer.ipython.org/github/biosustain/cameo-notebooks/blob/master/07-predict-heterologous-pathways.ipynb>`__

.. showcase-end


Contributions
~~~~~~~~~~~~~

... are very welcome! Please read the `guideline <CONTRIBUTING.rst>`__ for instructions how to contribute.


.. url-marker

.. |Join the chat at https://gitter.im/biosustain/cameo| image:: https://badges.gitter.im/biosustain/cameo.svg
   :target: https://gitter.im/biosustain/cameo?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge
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
.. |zenhub| image:: https://img.shields.io/badge/Shipping_faster_with-ZenHub-5e60ba.svg?style=flat-square
   :target: https://zenhub.com
