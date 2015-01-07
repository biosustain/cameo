
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

    <IPython.core.display.Javascript at 0x11802d450>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x11802d490>



.. raw:: html

    <div class="pb" id="38aca156-f895-4256-9f08-6a157e7bbd5a"><table class="pb ui-widget"><tr>
    <td class="pb_widget">Processing</td>
    <td class="pb_widget"><div id="3d1b0bcd-2100-4795-b8ff-ff3ee0fe0d7a">Escherichia coli</div></td>
    <td class="pb_widget_fill">
            <div class="pb_bar" id="d9b360f3-1c38-42a7-86e8-66ee69e3e955"></div>
            <script type="text/javascript">
                $("div#d9b360f3-1c38-42a7-86e8-66ee69e3e955").progressbar({value: 0, max: 2});
            </script>
            </td>
    <td class="pb_widget"><div id="23fcfffa-4199-4ad2-9524-e438f188391b">ETA:  --:--:--</div></td>
    <td class="pb_widget"><div id="738534f1-e60d-43cc-b5eb-c8204f605ee9">  0%</div></td>
    </tr></table><div>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x11802d850>



.. parsed-literal::

    <IPython.core.display.Javascript at 0x11802d950>



.. raw:: html

    <div class="pb" id="c5c4033d-17a4-4eb7-84fe-406bcd3cce84"><table class="pb ui-widget"><tr>
    <td class="pb_widget">Processing</td>
    <td class="pb_widget"><div id="21417de5-31a1-4796-b831-32bba54d9b17">iJO1366</div></td>
    <td class="pb_widget_fill">
            <div class="pb_bar" id="10f5171b-b208-44a1-91d2-b71d5e552c33"></div>
            <script type="text/javascript">
                $("div#10f5171b-b208-44a1-91d2-b71d5e552c33").progressbar({value: 0, max: 1});
            </script>
            </td>
    <td class="pb_widget"><div id="51357a8d-03f0-497d-a573-5c1d589bbe82">ETA:  --:--:--</div></td>
    <td class="pb_widget"><div id="1a9188ef-ef47-4005-a90b-0e112b3fe35f">  0%</div></td>
    </tr></table><div>


.. parsed-literal::

    CPLEX Error  1561: Not enough data in SAV file.
    CPLEX Error  1561: Not enough data in SAV file.


::


    ---------------------------------------------------------------------------
    CplexSolverError                          Traceback (most recent call last)

    <ipython-input-3-7607d1c24e0c> in <module>()
          1 from cameo import api
    ----> 2 report = api.design(product='L-serine')
    

    /Users/niko/Arbejder/Dev/cameo/cameo/api/designer.pyc in __call__(self, product, hosts)
         67         """
         68         product = self.__translate_product_to_universal_reactions_model_metabolite(product)
    ---> 69         pathways = self.predict_pathways(product, hosts=hosts)
         70         return pathways
         71 


    /Users/niko/Arbejder/Dev/cameo/cameo/api/designer.pyc in predict_pathways(self, product, hosts)
         93             for model in pbar_models(list(host.models)):
         94                 # TODO: Check if product is already part of model
    ---> 95                 pathway_predictor = PathwayPredictor(model, universal_model=METANETX['universal_model'], mapping=METANETX['bigg2mnx'])
         96                 predicted_pathways = pathway_predictor.run(product, max_predictions=5, timeout=15*60)  # TODO adjust these numbers to something reasonable
         97                 pathways[(host, model)] = predicted_pathways


    /Users/niko/Arbejder/Dev/cameo/cameo/strain_design/pathway_prediction/__init__.pyc in __init__(self, model, universal_model, mapping)
         66             from cameo.api import _METANETX as metanetx
         67             universal_model = metanetx['universal_model']
    ---> 68         self.model = model.copy()
         69         for exchange in self.model.exchanges:
         70             if len(exchange.reactants) > 0 and exchange.lower_bound <= 0:


    /Users/niko/Arbejder/Dev/cameo/cameo/solver_based_model.pyc in copy(self)
        589             model_copy._solver = deepcopy(self.solver)
        590         except:  # Cplex has an issue with deep copies
    --> 591             model_copy._solver = copy(self.solver)
        592         return model_copy
        593 


    /usr/local/Cellar/python/2.7.8_2/Frameworks/Python.framework/Versions/2.7/lib/python2.7/copy.pyc in copy(x)
         94                 raise Error("un(shallow)copyable object of type %s" % cls)
         95 
    ---> 96     return _reconstruct(x, rv, 0)
         97 
         98 


    /usr/local/Cellar/python/2.7.8_2/Frameworks/Python.framework/Versions/2.7/lib/python2.7/copy.pyc in _reconstruct(x, info, deep, memo)
        334             state = deepcopy(state, memo)
        335         if hasattr(y, '__setstate__'):
    --> 336             y.__setstate__(state)
        337         else:
        338             if isinstance(state, tuple) and len(state) == 2:


    /Users/niko/Arbejder/Dev/optlang/optlang/cplex_interface.pyc in __setstate__(self, repr_dict)
        479         tmp_file = tempfile.mktemp(suffix=".sav")
        480         open(tmp_file, 'w').write(repr_dict['cplex_repr'])
    --> 481         problem = cplex.Cplex(tmp_file)
        482         self.__init__(problem=problem)
        483 


    /Users/niko/.virtualenvs/cameo2.7.8/lib/python2.7/site-packages/cplex/__init__.pyc in __init__(self, *args)
        713                     filetype = args[1]
        714                 self._lp = _internal._procedural.createprob(env._e, filename, env.parameters.read.apiencoding.get())
    --> 715                 _internal._procedural.readcopyprob(env._e, self._lp, filename, filetype)
        716             else:
        717                 self._lp = _internal._procedural.createprob(env._e, "", env.parameters.read.apiencoding.get())


    /Users/niko/.virtualenvs/cameo2.7.8/lib/python2.7/site-packages/cplex/_internal/_procedural.pyc in readcopyprob(env, lp, filename, filetype)
        399     else:
        400         status = CR.CPXXreadcopyprob(env, lp, filename, filetype)
    --> 401     check_status(env, status)
        402     return
        403 


    /Users/niko/.virtualenvs/cameo2.7.8/lib/python2.7/site-packages/cplex/_internal/_procedural.pyc in __call__(self, env, status, from_cb)
        115                 else:
        116                     error_string = geterrorstring(env, status)
    --> 117             raise CplexSolverError(error_string, env, status)
        118 
        119 check_status = StatusChecker()


    CplexSolverError: CPLEX Error  1561: Not enough data in SAV file.



IPython notebook
~~~~~~~~~~~~~~~~

Click
`here <http://nbviewer.ipython.org/github/biosustain/cameo/blob/devel/docs/cameo_high_level_interface.ipynb>`__
to download this page as an IPython notebook.
