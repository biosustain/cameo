
Low-level interface for developers
==================================

.. code:: python

    from cameo import load_model
    model = load_model('iJO1366')
Structure and design
--------------------

Relation to cobrapy
~~~~~~~~~~~~~~~~~~~

Formulating and solving optimization problems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

http://optlang.readthedocs.org/en/latest/

Best practices
--------------

Avoid copies at all costs
~~~~~~~~~~~~~~~~~~~~~~~~~

Copying models is expensive!

.. code:: python

    %time model_copy = model.copy()

.. parsed-literal::

    CPU times: user 3.17 s, sys: 84.7 ms, total: 3.26 s
    Wall time: 3.25 s


To avoid excessive copying of models, we suggest the following utility
included in cameo. For example, temporary reaction deletion, solve the
model and magically returns to it original state.

.. code:: python

    from functools import partial
    from cameo.util import TimeMachine
    with TimeMachine() as tm:
        tm
TimeMachine performs ``do`` steps immediately and stores all
corresponding ``undo`` actions in order. Running all ``undo`` steps can
be achieved with ``TimeMachine.reset()`` (wich is run) We recommend
using ``TimeMachine`` in conjunction with ``with``, as in case of an
unforseen exception all ``undo`` steps will be run and manipulated
models will return to their original states. This is more elegant than
using explicit ``try-except-else-finally`` statements.

IPython notebook
~~~~~~~~~~~~~~~~

Click
`here <http://nbviewer.ipython.org/github/biosustain/cameo/blob/devel/docs/cameo_low_level_interface.ipynb>`__
to download this page as an IPython notebook.
