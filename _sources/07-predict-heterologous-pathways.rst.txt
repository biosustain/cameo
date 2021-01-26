
.. code:: ipython3

    from IPython.display import display
    import re

Predict heterologous pathways
=============================

Predicting heterologous pathways is an important strategy to generate
new viable strains. Because portfolio of available reactions is very
large, computer assisted pathway design becomes essential. **Cameo**
implements a shortest pathways search algorithm using an universal
biochemical reaction database.

.. raw:: html

   <div class="alert alert-warning">

If you’re running this notebook on
`try.cameo.bio <http://try.cameo.bio>`__, things might run very slow due
to our inability to provide access to the proprietary
`CPLEX <https://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/>`__
solver on a public webserver. Furthermore, Jupyter kernels might crash
and restart due to memory limitations on the server.

.. raw:: html

   </div>

.. code:: ipython3

    from cameo import models
    from cameo.strain_design import pathway_prediction

.. code:: ipython3

    model = models.bigg.iMM904

.. code:: ipython3

    predictor = pathway_prediction.PathwayPredictor(model)

.. code:: ipython3

    pathways = predictor.run(product="vanillin", max_predictions=4)



.. raw:: html

    <span>Pathway 1</span>



.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>equation</th>
          <th>lower_bound</th>
          <th>upper_bound</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>MNXR5340</th>
          <td>H(+) + NADH + O2 + vanillate &lt;=&gt; H2O + 3,4-dih...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR5336</th>
          <td>2.0 H(+) + NADH + vanillate &lt;=&gt; H2O + vanillin...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR68718</th>
          <td>H2O + 3,4-dihydroxybenzoate &lt;=&gt; 3-dehydroshiki...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
      </tbody>
    </table>
    </div>


.. parsed-literal::

    Max flux: 3.36842



.. raw:: html

    <span>Pathway 2</span>



.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>equation</th>
          <th>lower_bound</th>
          <th>upper_bound</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>MNXR5340</th>
          <td>H(+) + NADH + O2 + vanillate &lt;=&gt; H2O + 3,4-dih...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR5336</th>
          <td>2.0 H(+) + NADH + vanillate &lt;=&gt; H2O + vanillin...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR230</th>
          <td>H(+) + 4-hydroxybenzoate + O2 + NADPH &lt;=&gt; H2O ...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
      </tbody>
    </table>
    </div>


.. parsed-literal::

    Max flux: 1.90533



.. raw:: html

    <span>Pathway 3</span>



.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>equation</th>
          <th>lower_bound</th>
          <th>upper_bound</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>MNXR4008</th>
          <td>H(+) + 3-oxoadipate &lt;=&gt; H2O + 5-oxo-4,5-dihydr...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR184</th>
          <td>3-oxoadipyl-CoA + succinate &lt;=&gt; 3-oxoadipate +...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR5340</th>
          <td>H(+) + NADH + O2 + vanillate &lt;=&gt; H2O + 3,4-dih...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR5336</th>
          <td>2.0 H(+) + NADH + vanillate &lt;=&gt; H2O + vanillin...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR228</th>
          <td>CO2 + 5-oxo-4,5-dihydro-2-furylacetate &lt;=&gt; H(+...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR4119</th>
          <td>2.0 H(+) + 3-carboxy-cis,cis-muconate &lt;=&gt; 3,4-...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR209</th>
          <td>CoA + 3-oxoadipyl-CoA &lt;=&gt; acetyl-CoA + succiny...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR3655</th>
          <td>2-(carboxymethyl)-5-oxo-2,5-dihydro-2-furoate ...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
      </tbody>
    </table>
    </div>


.. parsed-literal::

    Max flux: 5.59223



.. raw:: html

    <span>Pathway 4</span>



.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>equation</th>
          <th>lower_bound</th>
          <th>upper_bound</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>MNXR5338</th>
          <td>2.0 H(+) + NADH + 3,4-dihydroxybenzoate &lt;=&gt; H2...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR1041</th>
          <td>diphosphate + AMP + caffeoyl-CoA &lt;=&gt; CoA + ATP...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR4974</th>
          <td>O2 + 2.0 trans-4-coumarate &lt;=&gt; 2.0 trans-caffeate</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR227</th>
          <td>diphosphate + AMP + 4-coumaroyl-CoA &lt;=&gt; CoA + ...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR5340</th>
          <td>H(+) + NADH + O2 + vanillate &lt;=&gt; H2O + 3,4-dih...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR5336</th>
          <td>2.0 H(+) + NADH + vanillate &lt;=&gt; H2O + vanillin...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR18369</th>
          <td>CoA + H2O + 4-coumaroyl-CoA + NAD(+) &lt;=&gt; H(+) ...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR232</th>
          <td>H(+) + CoA + 4-hydroxybenzoate &lt;=&gt; H2O + 4-hyd...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR1039</th>
          <td>acetyl-CoA + 3,4-dihydroxybenzaldehyde &lt;=&gt; H2O...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
      </tbody>
    </table>
    </div>


.. parsed-literal::

    Max flux: 2.24390

