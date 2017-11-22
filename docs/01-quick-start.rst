
Getting started with cameo
==========================

**cameo** reuses and extends model data structures defined by
`cobrapy <https://opencobra.github.io/cobrapy/>`__
(**CO**\ nstraints-\ **B**\ ased **R**\ econstruction and **A**\ nalysis
tool for **Py**\ thon). So, in addition to following this quick start
guide and other **cameo** tutorials, we encourage you to explore
cobrapy's
`documentation <https://cobrapy.readthedocs.org/en/latest/cobra.core.html>`__
as well.

Step 1: Load a model
--------------------

Loading a model is easy. Just import the `~cameo.io.load_model`
function.

.. code:: ipython3

    from cameo import load_model

For example, load the most current genome-scale metabolic reconstruction
of *Escherichia coli*.

.. code:: ipython3

    model = load_model("iJO1366")

Models, reactions, metabolites, etc., provide return HTML when evaluated
in Jupyter notebooks and can thus be easily inspected.

.. code:: ipython3

    model




.. raw:: html

    
            <table>
                <tr>
                    <td><strong>Name</strong></td>
                    <td>iJO1366</td>
                </tr><tr>
                    <td><strong>Memory address</strong></td>
                    <td>0x0115c8ceb8</td>
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
                    <td>periplasm, cytosol, extracellular space</td>
                </tr>
              </table>



Step 2: Simulate a model
------------------------

The model can be simulated by executing
`~cameo.core.solver_based_model.SolverBasedModel.solve`.

.. code:: ipython3

    solution = model.optimize()

A quick overview of the solution can be obtained in form of a pandas
`DataFrame <http://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.html>`__
(all solution objects in cameo provide access to data frames through a
`data_frame` attribute).

.. code:: ipython3

    solution




.. raw:: html

    <h3>Optimal solution with objective value 0.982</h3><br><div>
    <style>
        .dataframe thead tr:only-child th {
            text-align: right;
        }
    
        .dataframe thead th {
            text-align: left;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
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
          <th>DM_4crsol_c</th>
          <td>2.1907e-04</td>
          <td>0.0000</td>
        </tr>
        <tr>
          <th>DM_5drib_c</th>
          <td>2.2103e-04</td>
          <td>0.0000</td>
        </tr>
        <tr>
          <th>DM_aacald_c</th>
          <td>-0.0000e+00</td>
          <td>0.0000</td>
        </tr>
        <tr>
          <th>DM_amob_c</th>
          <td>1.9647e-06</td>
          <td>0.0000</td>
        </tr>
        <tr>
          <th>DM_mththf_c</th>
          <td>4.4010e-04</td>
          <td>0.0000</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>ZN2abcpp</th>
          <td>0.0000e+00</td>
          <td>-0.0083</td>
        </tr>
        <tr>
          <th>ZN2t3pp</th>
          <td>0.0000e+00</td>
          <td>-0.0021</td>
        </tr>
        <tr>
          <th>ZN2tpp</th>
          <td>3.3499e-04</td>
          <td>0.0000</td>
        </tr>
        <tr>
          <th>ZNabcpp</th>
          <td>0.0000e+00</td>
          <td>-0.0083</td>
        </tr>
        <tr>
          <th>Zn2tex</th>
          <td>3.3499e-04</td>
          <td>-0.0000</td>
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
    <style>
        .dataframe thead tr:only-child th {
            text-align: right;
        }
    
        .dataframe thead th {
            text-align: left;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
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
          <th>DM_4crsol_c</th>
          <td>2.1907e-04</td>
          <td>0.0000</td>
        </tr>
        <tr>
          <th>DM_5drib_c</th>
          <td>2.2103e-04</td>
          <td>0.0000</td>
        </tr>
        <tr>
          <th>DM_aacald_c</th>
          <td>-0.0000e+00</td>
          <td>0.0000</td>
        </tr>
        <tr>
          <th>DM_amob_c</th>
          <td>1.9647e-06</td>
          <td>0.0000</td>
        </tr>
        <tr>
          <th>DM_mththf_c</th>
          <td>4.4010e-04</td>
          <td>0.0000</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>ZN2abcpp</th>
          <td>0.0000e+00</td>
          <td>-0.0083</td>
        </tr>
        <tr>
          <th>ZN2t3pp</th>
          <td>0.0000e+00</td>
          <td>-0.0021</td>
        </tr>
        <tr>
          <th>ZN2tpp</th>
          <td>3.3499e-04</td>
          <td>0.0000</td>
        </tr>
        <tr>
          <th>ZNabcpp</th>
          <td>0.0000e+00</td>
          <td>-0.0083</td>
        </tr>
        <tr>
          <th>Zn2tex</th>
          <td>3.3499e-04</td>
          <td>-0.0000</td>
        </tr>
      </tbody>
    </table>
    <p>2583 rows × 2 columns</p>
    </div>



Data frames make it very easy to process results. For example, let's
take a look at reactions with flux != 0

.. code:: ipython3

    solution.to_frame().query('fluxes != 0')




.. raw:: html

    <div>
    <style>
        .dataframe thead tr:only-child th {
            text-align: right;
        }
    
        .dataframe thead th {
            text-align: left;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
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
          <th>DM_4crsol_c</th>
          <td>2.1907e-04</td>
          <td>0.0000e+00</td>
        </tr>
        <tr>
          <th>DM_5drib_c</th>
          <td>2.2103e-04</td>
          <td>0.0000e+00</td>
        </tr>
        <tr>
          <th>DM_amob_c</th>
          <td>1.9647e-06</td>
          <td>0.0000e+00</td>
        </tr>
        <tr>
          <th>DM_mththf_c</th>
          <td>4.4010e-04</td>
          <td>0.0000e+00</td>
        </tr>
        <tr>
          <th>BIOMASS_Ec_iJO1366_core_53p95M</th>
          <td>9.8237e-01</td>
          <td>1.8492e-15</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>UPPDC1</th>
          <td>2.1907e-04</td>
          <td>0.0000e+00</td>
        </tr>
        <tr>
          <th>USHD</th>
          <td>1.9113e-02</td>
          <td>0.0000e+00</td>
        </tr>
        <tr>
          <th>VALTA</th>
          <td>-4.1570e-01</td>
          <td>0.0000e+00</td>
        </tr>
        <tr>
          <th>ZN2tpp</th>
          <td>3.3499e-04</td>
          <td>0.0000e+00</td>
        </tr>
        <tr>
          <th>Zn2tex</th>
          <td>3.3499e-04</td>
          <td>-0.0000e+00</td>
        </tr>
      </tbody>
    </table>
    <p>437 rows × 2 columns</p>
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
                    <td><strong>Stoichiometry</strong></td>
                    <td>
                        <p style='text-align:right'>3pg_c + atp_c <=> 13dpg_c + adp_c</p>
                        <p style='text-align:right'>3-Phospho-D-glycerate + ATP <=> 3-Phospho-D-glyceroyl phosphate + ADP</p>
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
                    <td><strong>Stoichiometry</strong></td>
                    <td>
                        <p style='text-align:right'>e4p_c + h2o_c + nad_c <=> 4per_c + 2.0 h_c + nadh_c</p>
                        <p style='text-align:right'>D-Erythrose 4-phosphate + H2O + Nicotinamide adenine dinucleotide <=> 4-Phospho-D-erythronate + 2.0 H+ + Nicotinamide adenine dinucleotide - reduced</p>
                    </td>
                </tr><tr>
                    <td><strong>GPR</strong></td><td>b2927 or b1779</td>
                </tr><tr>
                    <td><strong>Lower bound</strong></td><td>-1000.0</td>
                </tr><tr>
                    <td><strong>Upper bound</strong></td><td>1000.0</td>
                </tr>
            </table>
            



Be aware though that due variable naming restrictions in Python dot
notation access to reactions (and other objects) might not work in some
cases.

.. code:: ipython3

    # model.reactions.12DGR120tipp  # uncommenting and running this cell will produce a syntax error

In these cases you need to use the `model.reactions.get_by_id`.

.. code:: ipython3

    model.reactions.get_by_id('12DGR120tipp')




.. raw:: html

    
            <table>
                <tr>
                    <td><strong>Reaction identifier</strong></td><td>12DGR120tipp</td>
                </tr><tr>
                    <td><strong>Name</strong></td><td>1,2 diacylglycerol transport via flipping (periplasm to cytoplasm, n-C12:0)</td>
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
D-glucose in the cytosolic compartment.

.. code:: ipython3

    model.metabolites.glc__D_c




.. raw:: html

    
            <table>
                <tr>
                    <td><strong>Metabolite identifier</strong></td><td>glc__D_c</td>
                </tr>
                <tr>
                    <td><strong>Name</strong></td><td>D-Glucose</td>
                </tr>
                <tr>
                    <td><strong>Formula</strong></td><td>C6H12O6</td>
                </tr>
            </table>



And it is easy to find the associated reactions

.. code:: ipython3

    model.metabolites.glc__D_c.reactions




.. parsed-literal::

    frozenset({<Reaction MLTG1 at 0x1163779b0>,
               <Reaction MLTG2 at 0x1163779e8>,
               <Reaction TRE6PH at 0x116557ac8>,
               <Reaction G6PP at 0x116206b00>,
               <Reaction MLTG3 at 0x116377ba8>,
               <Reaction GLCabcpp at 0x11623f3c8>,
               <Reaction MLTG4 at 0x116377c18>,
               <Reaction AMALT2 at 0x11605b438>,
               <Reaction GLCt2pp at 0x11623f470>,
               <Reaction MLTG5 at 0x116377c88>,
               <Reaction TREH at 0x116557d30>,
               <Reaction AMALT1 at 0x11605b588>,
               <Reaction AMALT3 at 0x11605b6a0>,
               <Reaction GLCATr at 0x1162336a0>,
               <Reaction AMALT4 at 0x11605b710>,
               <Reaction XYLI2 at 0x1165a1f28>,
               <Reaction HEX1 at 0x1162a7748>,
               <Reaction LACZ at 0x1162f0f98>,
               <Reaction GALS3 at 0x1162157f0>})



A list of the genes encoded in the model can be accessed via
`model.genes`.

.. code:: ipython3

    model.genes[0:10]




.. parsed-literal::

    [<Gene b2215 at 0x10c6f3780>,
     <Gene b1377 at 0x10950b4e0>,
     <Gene b0241 at 0x109351be0>,
     <Gene b0929 at 0x109351048>,
     <Gene b4035 at 0x109351d68>,
     <Gene b4033 at 0x109344b38>,
     <Gene b4034 at 0x115e60518>,
     <Gene b4032 at 0x115e60550>,
     <Gene b4036 at 0x115e60588>,
     <Gene b4213 at 0x115e605c0>]



A few additional attributes have been added that are not available in a
`cobrapy <https://opencobra.github.io/cobrapy/>`__ model. For example,
exchange reactions that allow certain metabolites to enter or leave the
model can be accessed through `model.exchanges`.

.. code:: ipython3

    model.exchanges[0:10]




.. parsed-literal::

    [<Reaction DM_4crsol_c at 0x115f44390>,
     <Reaction DM_5drib_c at 0x115f443c8>,
     <Reaction DM_aacald_c at 0x115f44400>,
     <Reaction DM_amob_c at 0x115f44438>,
     <Reaction DM_mththf_c at 0x115f44470>,
     <Reaction DM_oxam_c at 0x115f444a8>,
     <Reaction EX_12ppd__R_e at 0x115f44550>,
     <Reaction EX_12ppd__S_e at 0x115f44588>,
     <Reaction EX_14glucan_e at 0x115f445c0>,
     <Reaction EX_15dap_e at 0x115f445f8>]



Or, the current medium can be accessed through `model.medium`.

.. code:: ipython3

    model.medium.T




.. raw:: html

    <div>
    <style>
        .dataframe thead tr:only-child th {
            text-align: right;
        }
    
        .dataframe thead th {
            text-align: left;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>bound</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>EX_ca2_e</th>
          <td>1000.00</td>
        </tr>
        <tr>
          <th>EX_cbl1_e</th>
          <td>0.01</td>
        </tr>
        <tr>
          <th>EX_cl_e</th>
          <td>1000.00</td>
        </tr>
        <tr>
          <th>EX_co2_e</th>
          <td>1000.00</td>
        </tr>
        <tr>
          <th>EX_cobalt2_e</th>
          <td>1000.00</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
        </tr>
        <tr>
          <th>EX_sel_e</th>
          <td>1000.00</td>
        </tr>
        <tr>
          <th>EX_slnt_e</th>
          <td>1000.00</td>
        </tr>
        <tr>
          <th>EX_so4_e</th>
          <td>1000.00</td>
        </tr>
        <tr>
          <th>EX_tungs_e</th>
          <td>1000.00</td>
        </tr>
        <tr>
          <th>EX_zn2_e</th>
          <td>1000.00</td>
        </tr>
      </tbody>
    </table>
    <p>25 rows × 1 columns</p>
    </div>



It is also possible to get a list of essential reactions ...

.. code:: ipython3

    from cameo.flux_analysis.analysis import find_essential_reactions
    find_essential_reactions(model)[0:10]




.. parsed-literal::

    [<Reaction DM_4crsol_c at 0x115f44390>,
     <Reaction DM_5drib_c at 0x115f443c8>,
     <Reaction DM_amob_c at 0x115f44438>,
     <Reaction DM_mththf_c at 0x115f44470>,
     <Reaction BIOMASS_Ec_iJO1366_core_53p95M at 0x115f44518>,
     <Reaction EX_ca2_e at 0x115f5c3c8>,
     <Reaction EX_cl_e at 0x115f5c588>,
     <Reaction EX_cobalt2_e at 0x115f5c668>,
     <Reaction EX_cu2_e at 0x115f5c860>,
     <Reaction EX_glc__D_e at 0x115f697b8>]



... and essential genes.

.. code:: ipython3

    from cameo.flux_analysis.analysis import find_essential_genes
    find_essential_genes(model)[0:10]




.. parsed-literal::

    [<Gene b4245 at 0x115e90048>,
     <Gene b0109 at 0x115f08080>,
     <Gene b2838 at 0x115ea80f0>,
     <Gene b0423 at 0x115f380f0>,
     <Gene b2574 at 0x115e90128>,
     <Gene b3809 at 0x115ea8128>,
     <Gene b4407 at 0x115f38128>,
     <Gene b0175 at 0x115ea8160>,
     <Gene b3992 at 0x115f38160>,
     <Gene b0928 at 0x115e90198>]


