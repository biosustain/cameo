
Simulate models
===============

cameo uses and extends the model data structures defined by
`cobrapy <https://opencobra.github.io/cobrapy/>`__, our favorite
COnstraints-Based Reconstruction and Analysis tool for Python. cameo is
thus 100% compatible with cobrapy. For efficiency reasons, however,
cameo implements its own simulation methods that take advantage of a
more advanced solver interface.

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

Load a model.

.. code:: python

    from cameo import load_model
    model = load_model('iJO1366')

In cameo, flux balance analysis can be performed with the function
`fba`.

.. code:: python

    from cameo import fba
    %time fba_result = fba(model)


.. parsed-literal::

    CPU times: user 410 ms, sys: 6.43 ms, total: 416 ms
    Wall time: 421 ms


Basically, `fba` calls `model.solve()` and wraps the optimization
solution in a `FluxDistributionResult` object. The maximum objective
values (corresponding to a maximum growth rate) can obtained throug
`result.objective_value`.

.. code:: python

    fba_result.data_frame




.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>flux</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>DM_4crsol_c</th>
          <td>0.000219</td>
        </tr>
        <tr>
          <th>DM_5drib_c</th>
          <td>0.000221</td>
        </tr>
        <tr>
          <th>DM_aacald_c</th>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>DM_amob_c</th>
          <td>0.000002</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
        </tr>
        <tr>
          <th>ZN2t3pp</th>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>ZN2tpp</th>
          <td>0.000335</td>
        </tr>
        <tr>
          <th>ZNabcpp</th>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>Zn2tex</th>
          <td>0.000335</td>
        </tr>
      </tbody>
    </table>
    <p>2583 rows × 1 columns</p>
    </div>



Parsimonious Flux Balance Analysis
----------------------------------

Parsimonious flux balance analysis (`Lewis et al.,
2010 <http://www.ncbi.nlm.nih.gov/pubmed/20664636>`__), a variant of
FBA, performs FBA in in a first step to determine the maximum objective
value :math:`Z_{obj}`, fixes it in form of an additional model
constraint (:math:`\mathbf{c}^{T} \mathbf{v} \ge Z_{obj}`), and then
minimizes in a second optimization the :math:`L_1` norm of
:math:`\mathbf{v}`. The assumption behind pFBA is that cells try to
minimize flux magnitude as well in order to keep protein costs low.

.. math::


   \begin{align}
    Max ~ & ~ \lvert \mathbf{v} \rvert\\
    \text{s.t.}~ & ~ \mathbf{S} \mathbf{v} = 0 \\
    & ~ \mathbf{c}^{T} \mathbf{v} \ge Z_{obj} \\
    ~ & ~ \mathbf{v}_{lb} \leq \mathbf{v} \leq \mathbf{v}_{ub} \,.
   \end{align}

In cameo, pFBA can be performed with the function `pfba`.

.. code:: python

    from cameo import pfba
    %time pfba_result = pfba(model)


.. parsed-literal::

    CPU times: user 429 ms, sys: 19 ms, total: 448 ms
    Wall time: 460 ms


The `objective_function` value is :math:`\lvert \mathbf{v} \rvert` ...

.. code:: python

    pfba_result.objective_value




.. parsed-literal::

    699.0222751839472



... whis is significantly smaller than flux vector of the original FBA
solution.

.. code:: python

    abs(fba_result.data_frame.flux).sum()




.. parsed-literal::

    702.81946544045991



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
                    <td><strong>Name</strong></td><td>Glucose-6-phosphate isomerase</td>
                </tr>
                <tr>
                    <td><strong>Stoichiometry</strong></td><td>g6p_c <=> f6p_c</td>
                </tr>
                <tr>
                    <td><strong>Lower bound</strong></td><td>-1000.000000</td>
                </tr>
                <tr>
                    <td><strong>Upper bound</strong></td><td>1000.000000</td>
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
                    <td><strong>Name</strong></td><td>Glucose-6-phosphate isomerase</td>
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

    %time fba_knockout_result = fba(model)
    fba_knockout_result[model.reactions.BIOMASS_Ec_iJO1366_core_53p95M]


.. parsed-literal::

    CPU times: user 266 ms, sys: 4.96 ms, total: 271 ms
    Wall time: 272 ms




.. parsed-literal::

    0.9761293262947403



.. code:: python

    %time pfba_knockout_result = pfba(model)
    pfba_knockout_result[model.reactions.BIOMASS_Ec_iJO1366_core_53p95M]


.. parsed-literal::

    CPU times: user 374 ms, sys: 3.5 ms, total: 378 ms
    Wall time: 379 ms




.. parsed-literal::

    0.9761293262947374



MOMA and ROOM relly on a reference (wild-type) flux distribution and we
can use the one previously computed.

**Parsimonious FBA references seem to produce better results using this
methods**

.. code:: python

    from cameo.flux_analysis.simulation import room, lmoma

.. code:: python

    %time lmoma_result = lmoma(model, reference=pfba_result.fluxes)
    lmoma_result[model.reactions.BIOMASS_Ec_iJO1366_core_53p95M]


.. parsed-literal::

    CPU times: user 16.3 s, sys: 246 ms, total: 16.6 s
    Wall time: 16.8 s




.. parsed-literal::

    0.8724093536243601



.. code:: python

    %time room_result = room(model, reference=pfba_result.fluxes)
    room_result[model.reactions.BIOMASS_Ec_iJO1366_core_53p95M]

