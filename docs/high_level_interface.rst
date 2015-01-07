==============================
High-level interface for users
==============================

.. warning::
    The high-level interface is not fully implemented yet.

Users primarily interested in using cameo as a tool
for enumerating metabolic engineering strategies have access to cameo's advanced programming interface via :mod:`cameo.api`
that provides access to potential products (:mod:`cameo.api.products`), host organisms (:mod:`cameo.api.hosts`) and
a configurable design function (:attr:`cameo.api.design`).

Running :func:`cameo.api.design` requires only minimal input.

.. code-block:: python

    from cobra import api
    report = api.design(product='L-serine')


Get started with `IPython notebook <http://nbviewer.ipython.org/github/biosustain/cameo/blob/devel/docs/high_level_interface.ipynb>`_
