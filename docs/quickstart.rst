
Quickstart
==========

.. code:: python

    import pandas
    pandas.options.display.max_rows = 15
    from cameo import load_model

.. code:: python

    # model = load_model('Ecoli core Model')
    model = load_model('iJO1366')

.. code:: python

    solution = model.solve()

.. code:: python

    solution.to_frame()




.. raw:: html

    <div style="max-height:1000px;max-width:1500px;overflow:auto;">
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>fluxes</th>
          <th>reduced_costs</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>DM_4CRSOL</th>
          <td> 0.000219</td>
          <td> 0.000000</td>
        </tr>
        <tr>
          <th>DM_5DRIB</th>
          <td> 0.000221</td>
          <td> 0.000000</td>
        </tr>
        <tr>
          <th>DM_AACALD</th>
          <td> 0.000000</td>
          <td> 0.000000</td>
        </tr>
        <tr>
          <th>DM_AMOB</th>
          <td> 0.000002</td>
          <td> 0.000000</td>
        </tr>
        <tr>
          <th>DM_MTHTHF</th>
          <td> 0.000440</td>
          <td> 0.000000</td>
        </tr>
        <tr>
          <th>DM_OXAM</th>
          <td> 0.000000</td>
          <td> 0.000000</td>
        </tr>
        <tr>
          <th>Ec_biomass_iJO1366_WT_53p95M</th>
          <td> 0.000000</td>
          <td>-0.963591</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>XYLt2pp</th>
          <td> 0.000000</td>
          <td> 0.000000</td>
        </tr>
        <tr>
          <th>XYLtex</th>
          <td> 0.000000</td>
          <td> 0.000000</td>
        </tr>
        <tr>
          <th>ZN2abcpp</th>
          <td> 0.000000</td>
          <td>-0.004148</td>
        </tr>
        <tr>
          <th>ZN2t3pp</th>
          <td> 0.000000</td>
          <td>-0.001037</td>
        </tr>
        <tr>
          <th>ZN2tpp</th>
          <td> 0.000335</td>
          <td> 0.000000</td>
        </tr>
        <tr>
          <th>ZNabcpp</th>
          <td> 0.000000</td>
          <td>-0.004148</td>
        </tr>
        <tr>
          <th>Zn2tex</th>
          <td> 0.000335</td>
          <td> 0.000000</td>
        </tr>
      </tbody>
    </table>
    <p>2583 rows × 2 columns</p>
    </div>



