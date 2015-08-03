
.. code:: python

    from pandas import options
    options.display.max_rows = 8
    from cameo import load_model
    model = load_model("iJO1366")

Simulating models with
======================

**c**\ omputer **a**\ ided **m**\ etabolic **e**\ ngineering and
**o**\ ptimization

**cameo** uses and extends the model data structures defined by
`cobrapy <https://opencobra.github.io/cobrapy/>`__, our favorite
**CO**\ nstraints-\ **B**\ ased **R**\ econstruction and **A**\ nalysis
tool for **Py**\ thon. **cameo** is thus 100% compatible with
**cobrapy**. For efficiency reasons, however, **cameo** implements its
own simulation methods that take advantage of a more advanced solver
interface.

Primer: Constraint-Based Modeling
---------------------------------

Constraint-based modeling is a powerful modeling framework for analyzing
metabolism on the genome scale (`McCloskey et al.,
2013 <http://www.ncbi.nlm.nih.gov/pubmed/23632383>`__). For a model that
encompasses :math:`n` reactions that involve :math:`m` metabolites,
:math:`\mathbf{S}` is a matrix of dimension :math:`m \times n` that
encodes the stoichiometry of the metabolic reaction system; it is
usually referred to as stoichiometric matrix. Assuming that the system
is in a steady state—the concentration of metabolites are constant—the
system of flux-balances can be formulated as

.. math::


   \begin{align}
   \mathbf{S} \mathbf{v} = 0\,,
   \end{align}

where :math:`\mathbf{v}` is the vector of flux rates. With the addition
of a biologically meaningful objective, flux capacity constraints,
information about the reversibility of reactions under physiological
conditions, an optimization problem can be formulated that can easily be
solved using `linear
programming <https://en.wikipedia.org/wiki/Linear_programming>`__.

, e.g., maximimization of biomass production,Given the maximization of
growth rate as one potential biological objective :math:`v_{biomass}`,
i.e., the flux of an artificial reaction that consumes biomass
components in empirically determined proportions, and assuming that the
cell is evolutionary optimized to achieve that objective, and
incorporating knowledge about reaction reversibility, uptake and
secretion rates, and maximum flux capacities in the form of lower and
uppers bounds (:math:`\mathbf{v}_{lb}` and :math:`\mathbf{v}_{ub}`) on
the flux variables :math:`\mathbf{v}`, one can formulate and solve an
optimization problem to identify an optimal set of flux rates using flux
balance analysis (FBA):

.. math::


   \begin{align}
    Max ~ & ~ Z_{obj} = \mathbf{c}^{T} \mathbf{v}\\
    \text{s.t.}~ & ~ \mathbf{S} \mathbf{v} = 0 \\
    ~ & ~ \mathbf{v}_{lb} \leq \mathbf{v} \leq \mathbf{v}_{ub} \,.
   \end{align}

Flux Balance Analysis
---------------------

In **cameo**, flux balance analysis can be performed with the function
``fba``.

.. code:: python

    from cameo import fba
    fba_result = fba(model)

Basically, ``fba`` calls ``model.solve()`` and wraps the optimization
solution in a ``FluxDistributionResult`` object. The maximum objective
values (corresponding to a maximum growth rate) can obtained throug
``result.objective_value``.

.. code:: python

    fba_result.objective_value




.. parsed-literal::

    0.9823718127269799



Parsimonious Flux Balance Analysis
----------------------------------

Parsimonious flux balance analysis (`Lewis et al.,
2010 <http://www.ncbi.nlm.nih.gov/pubmed/20664636>`__), a variant of
FBA, performs FBA in in a first step to determine the maximum objective
value :math:`Z_{obj}`, fixes it in form of an additional model
constraint (:math:`\mathbf{c}^{T} \mathbf{v} \ge Z_{obj}`), and then
minimizes in a second optimization the :math:`L_1` norm of
:math:`\mathbf{v}`. The assumption behind the pFBA is that cells try to
minimize flux magnitude as well in order to keep the costs of protein
low.

.. math::


   \begin{align}
    Max ~ & ~ \lvert \mathbf{v} \rvert\\
    \text{s.t.}~ & ~ \mathbf{S} \mathbf{v} = 0 \\
    & ~ \mathbf{c}^{T} \mathbf{v} \ge Z_{obj} \\
    ~ & ~ \mathbf{v}_{lb} \leq \mathbf{v} \leq \mathbf{v}_{ub} \,.
   \end{align}

In **cameo**, pFBA can be performed with the function ``pfba``.

.. code:: python

    from cameo import pfba
    pfba_result = pfba(model)

The ``objective_function`` value is :math:`\lvert \mathbf{v} \rvert` ...

.. code:: python

    pfba_result.objective_value




.. parsed-literal::

    699.0222751839377



... whis is significantly smaller than flux vector of the original FBA
solution.

.. code:: python

    abs(fba_result.data_frame.flux).sum()




.. parsed-literal::

    764.91487969777245



Setp 2: Simulate knockouts phenotypes
-------------------------------------

Although PFBA and FBA can be used to simulate the effect of knockouts,
other methods have been proven more valuable for that task: MOMA and
ROOM. In *cameo* we implement a linear version of MOMA.

--------------

Simulating knockouts:

-  Manipulate the bounds of the reaction (or use the shorthand method
   knock\_out)

.. code:: python

    model.reactions.PGI




.. raw:: html

    
            <table>
                <tr>
                    <td><strong>Id</strong></td><td>PGI</td>
                </tr>
                <tr>
                    <td><strong>Stoichiometry</strong></td><td>g6p_c <=> f6p_c</td>
                </tr>
                <tr>
                    <td><strong>Lower bound</strong></td><td>-999999.000000</td>
                </tr>
                <tr>
                    <td><strong>Upper bound</strong></td><td>999999.000000</td>
                </tr>
            </table>
            



.. code:: python

    model.reactions.PGI.knock_out()
    model.reactions.PGI




.. raw:: html

    
            <table>
                <tr>
                    <td><strong>Id</strong></td><td>PGI</td>
                </tr>
                <tr>
                    <td><strong>Stoichiometry</strong></td><td>g6p_c --> f6p_c</td>
                </tr>
                <tr>
                    <td><strong>Lower bound</strong></td><td>0.000000</td>
                </tr>
                <tr>
                    <td><strong>Upper bound</strong></td><td>0.000000</td>
                </tr>
            </table>
            



-  Simulate using different methods:

.. code:: python

    %time
    fba_knockout_result = simulation.fba(model)
    fba_knockout_result[model.objective]


.. parsed-literal::

    CPU times: user 2 µs, sys: 0 ns, total: 2 µs
    Wall time: 5.01 µs




.. parsed-literal::

    0.905983



.. code:: python

    pfba_knockout_result = simulation.pfba(model)
    pfba_knockout_result[model.objective]




.. parsed-literal::

    0.905983



MOMA and ROOM relly on a reference (wild-type) flux distribution and we
can use the one previously computed.

**Parsimonious FBA references seem to produce better results using this
methods**

.. code:: python

    lmoma_result["2 * EX_glc_lp_e_rp_"]




.. parsed-literal::

    -18.7358



.. code:: python

    %time
    lmoma_result = simulation.lmoma(model, reference=pfba_result.fluxes)
    lmoma_result[model.objective]


.. parsed-literal::

    CPU times: user 2 µs, sys: 1 µs, total: 3 µs
    Wall time: 5.01 µs




.. parsed-literal::

    0.791393



.. code:: python

    %time
    room_result = simulation.room(model, reference=pfba_result.fluxes)
    room_result[model.objective]


.. parsed-literal::

    CPU times: user 2 µs, sys: 1 µs, total: 3 µs
    Wall time: 5.01 µs




.. parsed-literal::

    0.887440



.. code:: python

    room_result




.. parsed-literal::

    <cameo.core.result.FluxDistributionResult at 0x10aa75b50>



