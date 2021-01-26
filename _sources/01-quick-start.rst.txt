
Getting started with cameo
==========================

**cameo** reuses and extends model data structures defined by
`cobrapy <https://opencobra.github.io/cobrapy/>`__
(**CO**\ nstraints-\ **B**\ ased **R**\ econstruction and **A**\ nalysis
tool for **Py**\ thon). So, in addition to following this quick start
guide and other **cameo** tutorials, we encourage you to explore
cobrapy’s
`documentation <https://cobrapy.readthedocs.org/en/latest/cobra.core.html>`__
as well.

Step 1: Load a model
--------------------

Loading a model is easy. Just import the `~cameo.io.load_model`
function.

.. code:: ipython3

    from cameo import load_model

For example, load a genome-scale metabolic reconstruction of
*Escherichia coli*.

.. code:: ipython3

    model = load_model("iJO1366")

Models, reactions, metabolites, etc., return HTML when evaluated in
Jupyter notebooks and can be easily inspected.

.. code:: ipython3

    model




.. raw:: html

    
            <table>
                <tr>
                    <td><strong>Name</strong></td>
                    <td>iJO1366</td>
                </tr><tr>
                    <td><strong>Memory address</strong></td>
                    <td>0x01120756d8</td>
                </tr><tr>
                    <td><strong>Number of metabolites</strong></td>
                    <td>1805</td>
                </tr><tr>
                    <td><strong>Number of reactions</strong></td>
                    <td>2583</td>
                </tr><tr>
                    <td><strong>Objective expression</strong></td>
                    <td>-1.0*BIOMASS_Ec_iJO1366_core_53p95M_reverse_5c8b1 + 1.0*BIOMASS_Ec_iJO1366_core_53p95M</td>
                </tr><tr>
                    <td><strong>Compartments</strong></td>
                    <td>extracellular space, cytosol, periplasm</td>
                </tr>
              </table>



Step 2: Simulate a model
------------------------

The model can be simulated by executing `optimize`.

.. code:: ipython3

    solution = model.optimize()

A quick overview of the solution can be obtained in inspecting it.

.. code:: ipython3

    solution




.. raw:: html

    <strong><em>Optimal</em> solution with objective value 0.982</strong><br><div>
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
          <th>fluxes</th>
          <th>reduced_costs</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>12DGR120tipp</th>
          <td>0.000000</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>12DGR140tipp</th>
          <td>0.000000</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>12DGR141tipp</th>
          <td>0.000000</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>12DGR160tipp</th>
          <td>0.000000</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>12DGR161tipp</th>
          <td>0.000000</td>
          <td>-0.008295</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>ZN2abcpp</th>
          <td>0.000000</td>
          <td>-0.008295</td>
        </tr>
        <tr>
          <th>ZN2t3pp</th>
          <td>0.000000</td>
          <td>-0.002074</td>
        </tr>
        <tr>
          <th>ZN2tpp</th>
          <td>0.000335</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>ZNabcpp</th>
          <td>0.000000</td>
          <td>-0.008295</td>
        </tr>
        <tr>
          <th>Zn2tex</th>
          <td>0.000335</td>
          <td>-0.000000</td>
        </tr>
      </tbody>
    </table>
    <p>2583 rows × 2 columns</p>
    </div>



A data frame representation of the solution is accessible via
`solution.to_frame()`.

.. code:: ipython3

    solution.to_frame()




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
          <th>fluxes</th>
          <th>reduced_costs</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>12DGR120tipp</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>12DGR140tipp</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>12DGR141tipp</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>12DGR160tipp</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>12DGR161tipp</th>
          <td>0.000000</td>
          <td>-8.295308e-03</td>
        </tr>
        <tr>
          <th>12DGR180tipp</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>12DGR181tipp</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>12PPDRtex</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>12PPDRtpp</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>12PPDStex</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>12PPDStpp</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>14GLUCANabcpp</th>
          <td>0.000000</td>
          <td>-5.551115e-17</td>
        </tr>
        <tr>
          <th>14GLUCANtexi</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>23CAMPtex</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>23CCMPtex</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>23CGMPtex</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>23CUMPtex</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>23DAPPAt2pp</th>
          <td>0.000000</td>
          <td>1.517883e-18</td>
        </tr>
        <tr>
          <th>23DAPPAtex</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>23PDE2pp</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>23PDE4pp</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>23PDE7pp</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>23PDE9pp</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>26DAHtex</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>2AGPA120tipp</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>2AGPA140tipp</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>2AGPA141tipp</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>2AGPA160tipp</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>2AGPA161tipp</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>2AGPA180tipp</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>VALTRS</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>VALabcpp</th>
          <td>0.000000</td>
          <td>-6.221481e-03</td>
        </tr>
        <tr>
          <th>VALt2rpp</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>VALtex</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>VPAMTr</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>WCOS</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>X5PL3E</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>XAND</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>XANt2pp</th>
          <td>0.000000</td>
          <td>-2.073827e-03</td>
        </tr>
        <tr>
          <th>XANtex</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>XANtpp</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>XMPtex</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>XPPT</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>XTSNH</th>
          <td>0.000000</td>
          <td>-8.295308e-03</td>
        </tr>
        <tr>
          <th>XTSNt2rpp</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>XTSNtex</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>XYLI1</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>XYLI2</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>XYLK</th>
          <td>0.000000</td>
          <td>-1.387779e-17</td>
        </tr>
        <tr>
          <th>XYLK2</th>
          <td>0.000000</td>
          <td>-1.387779e-17</td>
        </tr>
        <tr>
          <th>XYLUt2pp</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>XYLUtex</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>XYLabcpp</th>
          <td>0.000000</td>
          <td>-6.221481e-03</td>
        </tr>
        <tr>
          <th>XYLt2pp</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>XYLtex</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>ZN2abcpp</th>
          <td>0.000000</td>
          <td>-8.295308e-03</td>
        </tr>
        <tr>
          <th>ZN2t3pp</th>
          <td>0.000000</td>
          <td>-2.073827e-03</td>
        </tr>
        <tr>
          <th>ZN2tpp</th>
          <td>0.000335</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>ZNabcpp</th>
          <td>0.000000</td>
          <td>-8.295308e-03</td>
        </tr>
        <tr>
          <th>Zn2tex</th>
          <td>0.000335</td>
          <td>-0.000000e+00</td>
        </tr>
      </tbody>
    </table>
    <p>2583 rows × 2 columns</p>
    </div>



Data frames make it very easy to process results. For example, let’s
take a look at reactions with flux != 0

.. code:: ipython3

    solution.to_frame().query('fluxes != 0')




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
          <th>fluxes</th>
          <th>reduced_costs</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>3OAR140</th>
          <td>0.076452</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>3OAS140</th>
          <td>0.076452</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>5DOAN</th>
          <td>0.000221</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>A5PISO</th>
          <td>0.038226</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>AACPS3</th>
          <td>0.125378</td>
          <td>-1.387779e-16</td>
        </tr>
        <tr>
          <th>AACPS4</th>
          <td>0.147776</td>
          <td>-2.775558e-17</td>
        </tr>
        <tr>
          <th>AACPS7</th>
          <td>0.076452</td>
          <td>-2.775558e-17</td>
        </tr>
        <tr>
          <th>ACACT1r</th>
          <td>0.349606</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>ACACT2r</th>
          <td>0.349606</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>ACACT3r</th>
          <td>0.349606</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>ACACT4r</th>
          <td>0.349606</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>ACACT5r</th>
          <td>0.349606</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>ACACT6r</th>
          <td>0.273154</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>ACACT7r</th>
          <td>0.273154</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>ACCOAC</th>
          <td>0.076458</td>
          <td>-7.090682e-17</td>
        </tr>
        <tr>
          <th>ACGK</th>
          <td>0.290578</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>ACGS</th>
          <td>0.290578</td>
          <td>2.927346e-17</td>
        </tr>
        <tr>
          <th>ACHBS</th>
          <td>0.285408</td>
          <td>-6.938894e-18</td>
        </tr>
        <tr>
          <th>ACLS</th>
          <td>0.858857</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>ACOAD1f</th>
          <td>-0.349606</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>ACOAD2f</th>
          <td>-0.349606</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>ACOAD3f</th>
          <td>-0.349606</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>ACOAD4f</th>
          <td>-0.349606</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>ACOAD5f</th>
          <td>-0.349606</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>ACOAD6f</th>
          <td>-0.273154</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>ACOAD7f</th>
          <td>-0.125378</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>ACODA</th>
          <td>0.290578</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>ACONTa</th>
          <td>4.857777</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>ACONTb</th>
          <td>4.857777</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>ACOTA</th>
          <td>-0.290578</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>TMDS</th>
          <td>0.025705</td>
          <td>-5.551115e-17</td>
        </tr>
        <tr>
          <th>TMPK</th>
          <td>0.000219</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>TMPPP</th>
          <td>0.000219</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>TPI</th>
          <td>7.645371</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>TRDR</th>
          <td>0.243502</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>TRPAS2</th>
          <td>-0.055841</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>TRPS3</th>
          <td>0.055841</td>
          <td>2.775558e-17</td>
        </tr>
        <tr>
          <th>TYRL</th>
          <td>0.000219</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>TYRTA</th>
          <td>-0.135684</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>U23GAAT</th>
          <td>0.038226</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>UAAGDS</th>
          <td>0.027298</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>UAGAAT</th>
          <td>0.038226</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>UAGCVT</th>
          <td>0.027298</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>UAGDP</th>
          <td>0.092822</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>UAGPT3</th>
          <td>0.027298</td>
          <td>-5.551115e-17</td>
        </tr>
        <tr>
          <th>UAMAGS</th>
          <td>0.027298</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>UAMAS</th>
          <td>0.027298</td>
          <td>5.551115e-17</td>
        </tr>
        <tr>
          <th>UAPGR</th>
          <td>0.027298</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>UDCPDP</th>
          <td>0.027298</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>UDCPDPS</th>
          <td>0.000054</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>UGMDDS</th>
          <td>0.027298</td>
          <td>1.110223e-16</td>
        </tr>
        <tr>
          <th>UHGADA</th>
          <td>0.038226</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>UMPK</th>
          <td>0.371375</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>UPP3MT</th>
          <td>0.000219</td>
          <td>1.110223e-16</td>
        </tr>
        <tr>
          <th>UPP3S</th>
          <td>0.000438</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>UPPDC1</th>
          <td>0.000219</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>USHD</th>
          <td>0.019113</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>VALTA</th>
          <td>-0.415702</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>ZN2tpp</th>
          <td>0.000335</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>Zn2tex</th>
          <td>0.000335</td>
          <td>-0.000000e+00</td>
        </tr>
      </tbody>
    </table>
    <p>436 rows × 2 columns</p>
    </div>



Step 3: Exploring a model
-------------------------

Objects—models, reactions, metabolites, genes—can easily be explored in
the Jupyter notebook, taking advantage of tab completion. For example,
place your cursor after the period in `model.reactions.` and press the
TAB key. A dialog will appear that allows you to navigate the list of
reactions encoded in the model.

.. code:: ipython3

    model.reactions.PGK # delete PGK, place your cursor after the period and press the TAB key.




.. raw:: html

    
            <table>
                <tr>
                    <td><strong>Reaction identifier</strong></td><td>PGK</td>
                </tr><tr>
                    <td><strong>Name</strong></td><td>Phosphoglycerate kinase</td>
                </tr><tr>
                    <td><strong>Memory address</strong></td>
                    <td>0x01129829b0</td>
                </tr><tr>
                    <td><strong>Stoichiometry</strong></td>
                    <td>
                        <p style='text-align:right'>3pg_c + atp_c <=> 13dpg_c + adp_c</p>
                        <p style='text-align:right'>3-Phospho-D-glycerate + ATP C10H12N5O13P3 <=> 3-Phospho-D-glyceroyl phosphate + ADP C10H12N5O10P2</p>
                    </td>
                </tr><tr>
                    <td><strong>GPR</strong></td><td>b2926</td>
                </tr><tr>
                    <td><strong>Lower bound</strong></td><td>-1000.0</td>
                </tr><tr>
                    <td><strong>Upper bound</strong></td><td>1000.0</td>
                </tr>
            </table>
            



For example, you can access the E4PD (*Erythrose 4-phosphate
dehydrogenase*) reaction in the model.

.. code:: ipython3

    model.reactions.E4PD




.. raw:: html

    
            <table>
                <tr>
                    <td><strong>Reaction identifier</strong></td><td>E4PD</td>
                </tr><tr>
                    <td><strong>Name</strong></td><td>Erythrose 4-phosphate dehydrogenase</td>
                </tr><tr>
                    <td><strong>Memory address</strong></td>
                    <td>0x0112606160</td>
                </tr><tr>
                    <td><strong>Stoichiometry</strong></td>
                    <td>
                        <p style='text-align:right'>e4p_c + h2o_c + nad_c <=> 4per_c + 2.0 h_c + nadh_c</p>
                        <p style='text-align:right'>D-Erythrose 4-phosphate + H2O H2O + Nicotinamide adenine dinucleotide <=> 4-Phospho-D-erythronate + 2.0 H+ + Nicotinamide adenine dinucleotide - reduced</p>
                    </td>
                </tr><tr>
                    <td><strong>GPR</strong></td><td>b2927 or b1779</td>
                </tr><tr>
                    <td><strong>Lower bound</strong></td><td>-1000.0</td>
                </tr><tr>
                    <td><strong>Upper bound</strong></td><td>1000.0</td>
                </tr>
            </table>
            



Be aware that, due to variable naming restrictions in Python, dot
notation access to reactions (and other objects) might not work in some
cases.

.. code:: ipython3

    # model.reactions.12DGR120tipp  # uncommenting and running this cell will produce a syntax error

In those cases you need to use the `model.reactions.get_by_id`.

.. code:: ipython3

    model.reactions.get_by_id('12DGR120tipp')




.. raw:: html

    
            <table>
                <tr>
                    <td><strong>Reaction identifier</strong></td><td>12DGR120tipp</td>
                </tr><tr>
                    <td><strong>Name</strong></td><td>1,2 diacylglycerol transport via flipping (periplasm to cytoplasm, n-C12:0)</td>
                </tr><tr>
                    <td><strong>Memory address</strong></td>
                    <td>0x0112506ba8</td>
                </tr><tr>
                    <td><strong>Stoichiometry</strong></td>
                    <td>
                        <p style='text-align:right'>12dgr120_p --> 12dgr120_c</p>
                        <p style='text-align:right'>1,2-Diacyl-sn-glycerol (didodecanoyl, n-C12:0) --> 1,2-Diacyl-sn-glycerol (didodecanoyl, n-C12:0)</p>
                    </td>
                </tr><tr>
                    <td><strong>GPR</strong></td><td></td>
                </tr><tr>
                    <td><strong>Lower bound</strong></td><td>0.0</td>
                </tr><tr>
                    <td><strong>Upper bound</strong></td><td>1000.0</td>
                </tr>
            </table>
            



Metabolites are accessible through `model.metabolites`. For example,
D-glucose in the cytosol compartment.

.. code:: ipython3

    model.metabolites.glc__D_c




.. raw:: html

    
            <table>
                <tr>
                    <td><strong>Metabolite identifier</strong></td><td>glc__D_c</td>
                </tr><tr>
                    <td><strong>Name</strong></td><td>D-Glucose</td>
                </tr><tr>
                    <td><strong>Memory address</strong></td>
                    <td>0x01120db4a8</td>
                </tr><tr>
                    <td><strong>Formula</strong></td><td>C6H12O6</td>
                </tr><tr>
                    <td><strong>Compartment</strong></td><td>c</td>
                </tr><tr>
                    <td><strong>In 19 reaction(s)</strong></td><td>
                        GLCt2pp, GLCATr, TRE6PH, MLTG2, G6PP, LACZ, MLTG3, AMALT1, HEX1, XYLI2, MLTG4, GALS3, GLCabcpp, AMALT2, MLTG5, AMALT3, AMALT4, TREH, MLTG1</td>
                </tr>
            </table>



And it is easy to find the associated reactions

.. code:: ipython3

    model.metabolites.glc__D_c.reactions




.. parsed-literal::

    frozenset({<Reaction AMALT1 at 0x11257a978>,
               <Reaction AMALT2 at 0x11257aa58>,
               <Reaction AMALT3 at 0x11257aa90>,
               <Reaction AMALT4 at 0x11257ab00>,
               <Reaction G6PP at 0x112671940>,
               <Reaction GALS3 at 0x112679630>,
               <Reaction GLCATr at 0x1126864e0>,
               <Reaction GLCabcpp at 0x11268d240>,
               <Reaction GLCt2pp at 0x11268d438>,
               <Reaction HEX1 at 0x1126bf588>,
               <Reaction LACZ at 0x1126e2940>,
               <Reaction MLTG1 at 0x1129253c8>,
               <Reaction MLTG2 at 0x112925518>,
               <Reaction MLTG3 at 0x112925550>,
               <Reaction MLTG4 at 0x1129255c0>,
               <Reaction MLTG5 at 0x112925630>,
               <Reaction TRE6PH at 0x112a0b4e0>,
               <Reaction TREH at 0x112a0b748>,
               <Reaction XYLI2 at 0x112a30940>})



A list of the genes encoded in the model can be accessed via
`model.genes`.

.. code:: ipython3

    model.genes[0:10]




.. parsed-literal::

    [<Gene b0002 at 0x103f25748>,
     <Gene b0003 at 0x1034ec2e8>,
     <Gene b0004 at 0x1034fe3c8>,
     <Gene b0007 at 0x1034fe438>,
     <Gene b0008 at 0x1034f4780>,
     <Gene b0009 at 0x1034f4c88>,
     <Gene b0019 at 0x11237cd30>,
     <Gene b0025 at 0x11237cd68>,
     <Gene b0026 at 0x11237cda0>,
     <Gene b0029 at 0x11237cdd8>]



Other additional attributes can be accessed to explore the model. For
example, exchange reactions that allow certain metabolites to enter or
leave the model can be accessed through `model.exchanges`.

.. code:: ipython3

    model.exchanges[0:10]




.. parsed-literal::

    [<Reaction DM_4crsol_c at 0x1125f7160>,
     <Reaction DM_5drib_c at 0x1125f7470>,
     <Reaction DM_aacald_c at 0x1125f74a8>,
     <Reaction DM_amob_c at 0x1125f75c0>,
     <Reaction DM_mththf_c at 0x1125f75f8>,
     <Reaction DM_oxam_c at 0x1125f7630>,
     <Reaction EX_12ppd__R_e at 0x112613438>,
     <Reaction EX_12ppd__S_e at 0x112613470>,
     <Reaction EX_14glucan_e at 0x112613668>,
     <Reaction EX_15dap_e at 0x112613780>]



Or, the current medium can be accessed through `model.medium`.

.. code:: ipython3

    model.medium




.. parsed-literal::

    {'EX_ca2_e': 1000.0,
     'EX_cbl1_e': 0.01,
     'EX_cl_e': 1000.0,
     'EX_co2_e': 1000.0,
     'EX_cobalt2_e': 1000.0,
     'EX_cu2_e': 1000.0,
     'EX_fe2_e': 1000.0,
     'EX_fe3_e': 1000.0,
     'EX_glc__D_e': 10.0,
     'EX_h2o_e': 1000.0,
     'EX_h_e': 1000.0,
     'EX_k_e': 1000.0,
     'EX_mg2_e': 1000.0,
     'EX_mn2_e': 1000.0,
     'EX_mobd_e': 1000.0,
     'EX_na1_e': 1000.0,
     'EX_nh4_e': 1000.0,
     'EX_ni2_e': 1000.0,
     'EX_o2_e': 1000.0,
     'EX_pi_e': 1000.0,
     'EX_sel_e': 1000.0,
     'EX_slnt_e': 1000.0,
     'EX_so4_e': 1000.0,
     'EX_tungs_e': 1000.0,
     'EX_zn2_e': 1000.0}



It is also possible to get a list of essential reactions …

.. code:: ipython3

    from cobra.flux_analysis import find_essential_reactions
    list(find_essential_reactions(model))[0:10]




.. parsed-literal::

    [<Reaction THDPS at 0x1129fc080>,
     <Reaction UAAGDS at 0x112a1a080>,
     <Reaction PPNDH at 0x1129a60b8>,
     <Reaction CYSTL at 0x1125d40f0>,
     <Reaction E4PD at 0x112606160>,
     <Reaction KDOPP at 0x1126e2198>,
     <Reaction DHAD2 at 0x1125ea198>,
     <Reaction APRAUR at 0x1125881d0>,
     <Reaction METS at 0x11291e1d0>,
     <Reaction PAPPT3 at 0x11297c1d0>]



… and essential genes.

.. code:: ipython3

    from cobra.flux_analysis import find_essential_genes
    list(find_essential_genes(model))[0:10]




.. parsed-literal::

    [<Gene b1662 at 0x1123b6048>,
     <Gene b4245 at 0x112506128>,
     <Gene b3633 at 0x1123f0198>,
     <Gene b3634 at 0x1123f01d0>,
     <Gene b3639 at 0x1123f0208>,
     <Gene b4261 at 0x112506208>,
     <Gene b1415 at 0x1123ae240>,
     <Gene b4262 at 0x112506240>,
     <Gene b3642 at 0x1123f0278>,
     <Gene b1693 at 0x1123b6278>]


