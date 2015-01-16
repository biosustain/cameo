
High-level interface for users
==============================

Users primarily interested in using cameo as a tool for enumerating
metabolic engineering strategies have access to cameo's advanced
programming interface via :mod:``cameo.api`` that provides access to
potential products (:mod:``cameo.api.products``), host organisms
(:mod:``cameo.api.hosts``) and a configurable design function
(:attr:``cameo.api.design``). Running :func:``cameo.api.design``
requires only minimal input.

.. code:: python

    from cameo import api
    report = api.design(product='L-serine')

.. parsed-literal::

    Found 5 compounds that match query 'L-serine'
                         name     formula charge     mass  \
    MNXM53           L-serine     C3H7NO3      0  105.093   
    MNXM114         L-proline     C5H9NO2      0   115.13   
    MNXM3635       L-mimosine   C8H10N2O4      0  198.176   
    MNXM89905      L-allysine    C6H11NO3      0  145.156   
    MNXM384    L-saccharopine  C11H19N2O6     -1  275.278   
    
                                                           InChI  \
    MNXM53     InChI=1S/C3H7NO3/c4-2(1-5)3(6)7/h2,5H,1,4H2,(H...   
    MNXM114    InChI=1S/C5H9NO2/c7-5(8)4-2-1-3-6-4/h4,6H,1-3H...   
    MNXM3635   InChI=1S/C8H10N2O4/c9-5(8(13)14)3-10-2-1-6(11)...   
    MNXM89905  InChI=1S/C6H11NO3/c7-5(6(9)10)3-1-2-4-8/h4-5H,...   
    MNXM384    InChI=1S/C11H20N2O6/c12-7(10(16)17)3-1-2-6-13-...   
    
                                                          SMILES       source  \
    MNXM53                             [NH3+][C@@H](CO)C([O-])=O  chebi:33384   
    MNXM114                           [O-]C(=O)[C@@H]1CCC[NH2+]1  chebi:60039   
    MNXM3635           [NH3+][C@@H](CN1C=CC(=O)C(O)=C1)C([O-])=O  chebi:29063   
    MNXM89905                  [H]C(=O)CCC[C@H]([NH3+])C([O-])=O  chebi:58321   
    MNXM384    [NH3+][C@@H](CCCC[NH2+][C@@H](CCC([O-])=O)C([O...  chebi:57951   
    
               search_rank  
    MNXM53               0  
    MNXM114              1  
    MNXM3635             2  
    MNXM89905            3  
    MNXM384              4  
    Choosing best match (L-serine) ... please interrupt if this is not the desired compound.


.. parsed-literal::

    /Users/niko/Arbejder/Dev/cameo/cameo/api/products.py:75 [1;31mSettingWithCopyWarning[0m: 
    A value is trying to be set on a copy of a slice from a DataFrame.
    Try using .loc[row_indexer,col_indexer] = value instead
    
    See the the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy



.. parsed-literal::

    <IPython.core.display.Javascript at 0x114c52910>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x114c52950>



.. raw:: html

    <div class="pb" id="1acca05b-d708-447a-94e3-9eaa83e2e0ed"><table class="pb ui-widget"><tr>
    <td class="pb_widget">Processing</td>
    <td class="pb_widget"><div id="2a06b7ac-2238-478b-9760-047eaeeeae44">Escherichia coli</div></td>
    <td class="pb_widget_fill">
            <div class="pb_bar" id="449fe787-96dd-4956-8dce-2e41aa6993f7"></div>
            <script type="text/javascript">
                $("div#449fe787-96dd-4956-8dce-2e41aa6993f7").progressbar({value: 0, max: 2});
            </script>
            </td>
    <td class="pb_widget"><div id="3377f440-b83b-46a4-a8e4-86060a64f33a">ETA:  --:--:--</div></td>
    <td class="pb_widget"><div id="5e845f92-cc4a-4eb0-9bcd-40baaad86a94">  0%</div></td>
    </tr></table><div>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x114c52d10>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x114c52e10>



.. raw:: html

    <div class="pb" id="3f37040f-cb32-40ed-83d2-bef6286147e3"><table class="pb ui-widget"><tr>
    <td class="pb_widget">Processing</td>
    <td class="pb_widget"><div id="a788981e-381e-4f45-874b-42973f136864">iJO1366</div></td>
    <td class="pb_widget_fill">
            <div class="pb_bar" id="d61130ff-1f04-465f-a783-276b043f3cc1"></div>
            <script type="text/javascript">
                $("div#d61130ff-1f04-465f-a783-276b043f3cc1").progressbar({value: 0, max: 1});
            </script>
            </td>
    <td class="pb_widget"><div id="668a1c24-71c8-4533-afda-0dd7d7b5f15c">ETA:  --:--:--</div></td>
    <td class="pb_widget"><div id="030427bf-502f-4b09-8782-4aa484eee32e">  0%</div></td>
    </tr></table><div>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x11bbcaa90>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x11bbcaa50>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x1048f3e50>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x11aa711d0>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x11c280cd0>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x11aa71150>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x11aa71090>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x11aa71110>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x11aa5eed0>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x11aa5ee90>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x11aa5ef10>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x11aa5ef50>



.. raw:: html

    <div class="pb" id="047f9345-b4d4-4a35-92dd-176b6e9ecf2a"><table class="pb ui-widget"><tr>
    <td class="pb_widget">Processing</td>
    <td class="pb_widget"><div id="a78483e6-40a3-4970-9b1a-2a4343ad97cf">iMM904</div></td>
    <td class="pb_widget_fill">
            <div class="pb_bar" id="611cf549-2d62-4c76-9b78-3b9064116faa"></div>
            <script type="text/javascript">
                $("div#611cf549-2d62-4c76-9b78-3b9064116faa").progressbar({value: 0, max: 1});
            </script>
            </td>
    <td class="pb_widget"><div id="eed5af8c-c2aa-477f-a95f-0222357104a0">ETA:  --:--:--</div></td>
    <td class="pb_widget"><div id="95a2b3c2-9f75-4897-986d-6e20b9545061">  0%</div></td>
    </tr></table><div>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x114c52e90>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x11934d650>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x125185510>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x1251852d0>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x126889f10>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x1251851d0>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x125185250>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x125185190>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x125172e10>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x125172050>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x125185210>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x125172d90>


IPython notebook
~~~~~~~~~~~~~~~~

Click
`here <http://nbviewer.ipython.org/github/biosustain/cameo/blob/devel/docs/cameo_high_level_interface.ipynb>`__
to download this page as an IPython notebook.
