
High-level interface for users
==============================

Users primarily interested in using cameo as a tool for enumerating
metabolic engineering strategies have access to cameo's advanced
programming interface via ``cameo.api`` that provides access to
potential products (``cameo.api.products``), host organisms
(``cameo.api.hosts``) and a configurable design function
(``cameo.api.design``). Running ``cameo.api.design`` requires only
minimal input.

.. code:: python

    import logging
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

.. code:: python

    from cameo import api

.. code:: python

    report = api.design(product='3-hydroxy propionate')


.. parsed-literal::

    Found 5 compounds that match query '3-hydroxy propionate'
                                        name    formula charge     mass  \
    MNXM872              3-hydroxypropanoate     C3H5O3     -1    89.07   
    MNXM1526               3-hydroxypropanal     C3H6O2      0   74.079   
    MNXM10011   3-hydroxy-3-phenylpropionate     C9H9O3     -1  165.166   
    MNXM2955              6-hydroxyprotopine  C20H19NO6      0  369.368   
    MNXM101843    3-hydroxyphenyl propanoate    C9H10O3      0  166.174   
    
                                                            InChI  \
    MNXM872     InChI=1S/C3H6O3/c4-2-1-3(5)6/h4H,1-2H2,(H,5,6)...   
    MNXM1526               InChI=1S/C3H6O2/c4-2-1-3-5/h2,5H,1,3H2   
    MNXM10011   InChI=1S/C9H10O3/c10-8(6-9(11)12)7-4-2-1-3-5-7...   
    MNXM2955    InChI=1S/C20H19NO6/c1-21-8-14-11(2-3-16-20(14)...   
    MNXM101843  InChI=1S/C9H10O3/c1-2-9(11)12-8-5-3-4-7(10)6-8...   
    
                                                          SMILES       source  \
    MNXM872                                         OCCC([O-])=O  chebi:16510   
    MNXM1526                                         [H]C(=O)CCO  chebi:17871   
    MNXM10011                          OC(CC([O-])=O)C1=CC=CC=C1  chebi:63469   
    MNXM2955    CN1CC2=C3OCOC3=CC=C2CC(=O)C2=C(CC1O)C=C1OCOC1=C2  chebi:17104   
    MNXM101843                            CCC(=O)OC1=CC=CC(O)=C1  chebi:65228   
    
                search_rank  
    MNXM872               0  
    MNXM1526              1  
    MNXM10011             2  
    MNXM2955              3  
    MNXM101843            4  
    Choosing best match (3-hydroxypropanoate) ... please interrupt if this is not the desired compound.
    Predicting pathways for product 3-hydroxypropanoate and host Escherichia coli using model iJO1366.
    [[<Reaction MNXR161 at 0x1305fe210>, <Reaction adapter_h_c_MNXM1 at 0x135aa6c10>, <Reaction adapter_h2o2_c_MNXM22 at 0x135aa6850>, <Reaction adapter_hco3_c_MNXM60 at 0x135aa6750>, <Reaction adapter_msa_c_MNXM244 at 0x135a5aed0>, <Reaction adapter_nadph_c_MNXM6 at 0x135a5a650>, <Reaction adapter_h2o_p_MNXM2 at 0x13586c350>], [<Reaction MNXR8993 at 0x131ee0110>, <Reaction adapter_h_c_MNXM1 at 0x135aa6c10>, <Reaction adapter_h2o2_c_MNXM22 at 0x135aa6850>, <Reaction adapter_hco3_c_MNXM60 at 0x135aa6750>, <Reaction adapter_msa_c_MNXM244 at 0x135a5aed0>, <Reaction adapter_nadph_c_MNXM6 at 0x135a5a650>, <Reaction adapter_h2o_p_MNXM2 at 0x13586c350>], [<Reaction MNXR4458 at 0x130ad8990>, <Reaction MNXR8996 at 0x131ee0a50>, <Reaction MNXR9356 at 0x132047d10>, <Reaction adapter_2h3oppan_c_MNXM475 at 0x135bd6710>, <Reaction adapter_adp_c_MNXM7 at 0x135b65b90>, <Reaction adapter_db4p_c_MNXM2887 at 0x135b178d0>, <Reaction adapter_h2o_c_MNXM2 at 0x135aa6b10>, <Reaction adapter_h2o2_c_MNXM22 at 0x135aa6850>, <Reaction adapter_hco3_c_MNXM60 at 0x135aa6750>, <Reaction adapter_nadp_c_MNXM5 at 0x135a5a750>, <Reaction adapter_ppcoa_c_MNXM86 at 0x135a34290>, <Reaction adapter_pppi_c_MNXM332 at 0x135a21dd0>, <Reaction adapter_nh4_e_MNXM15 at 0x13592b210>, <Reaction adapter_pi_e_MNXM9 at 0x135918550>], [<Reaction MNXR4458 at 0x130ad8990>, <Reaction MNXR9031 at 0x131f2c410>, <Reaction MNXR9356 at 0x132047d10>, <Reaction adapter_2h3oppan_c_MNXM475 at 0x135bd6710>, <Reaction adapter_adp_c_MNXM7 at 0x135b65b90>, <Reaction adapter_db4p_c_MNXM2887 at 0x135b178d0>, <Reaction adapter_h2o_c_MNXM2 at 0x135aa6b10>, <Reaction adapter_h2o2_c_MNXM22 at 0x135aa6850>, <Reaction adapter_hco3_c_MNXM60 at 0x135aa6750>, <Reaction adapter_nadp_c_MNXM5 at 0x135a5a750>, <Reaction adapter_ppcoa_c_MNXM86 at 0x135a34290>, <Reaction adapter_pppi_c_MNXM332 at 0x135a21dd0>, <Reaction adapter_nh4_e_MNXM15 at 0x13592b210>, <Reaction adapter_pi_e_MNXM9 at 0x135918550>], [<Reaction MNXR4458 at 0x130ad8990>, <Reaction MNXR9031 at 0x131f2c410>, <Reaction MNXR15262 at 0x1328bdf50>, <Reaction adapter_2h3oppan_c_MNXM475 at 0x135bd6710>, <Reaction adapter_adp_c_MNXM7 at 0x135b65b90>, <Reaction adapter_db4p_c_MNXM2887 at 0x135b178d0>, <Reaction adapter_h2o_c_MNXM2 at 0x135aa6b10>, <Reaction adapter_h2o2_c_MNXM22 at 0x135aa6850>, <Reaction adapter_hco3_c_MNXM60 at 0x135aa6750>, <Reaction adapter_nadp_c_MNXM5 at 0x135a5a750>, <Reaction adapter_ppcoa_c_MNXM86 at 0x135a34290>, <Reaction adapter_pppi_c_MNXM332 at 0x135a21dd0>, <Reaction adapter_nh4_e_MNXM15 at 0x13592b210>, <Reaction adapter_pi_e_MNXM9 at 0x135918550>]]
    Predicting pathways for product 3-hydroxypropanoate and host Saccharomyces cerevisiae using model iMM904.
    [[<Reaction MNXR8993 at 0x1473d1950>, <Reaction MNXR32735 at 0x147fd1350>, <Reaction adapter_ala_dsh_B_c_MNXM144 at 0x14a3921d0>, <Reaction adapter_cer3_26_c_MNXM63157 at 0x14a36c850>, <Reaction adapter_pyr_c_MNXM23 at 0x14a29d4d0>, <Reaction adapter_atp_g_MNXM3 at 0x14a1f05d0>, <Reaction adapter_gdpmann_g_MNXM82 at 0x14a1f0290>, <Reaction adapter_h_g_MNXM1 at 0x14a1f0250>, <Reaction adapter_nadh_m_MNXM10 at 0x14a17f410>, <Reaction adapter_adp_x_MNXM7 at 0x14a10dad0>, <Reaction adapter_pmtcoa_x_MNXM88 at 0x14a0e7d50>, <Reaction adapter_ttdca_x_MNXM314 at 0x14a0e7550>], [<Reaction MNXR738 at 0x132553b10>, <Reaction MNXR8993 at 0x1473d1950>, <Reaction adapter_ala_dsh_B_c_MNXM144 at 0x14a3921d0>, <Reaction adapter_cer3_26_c_MNXM63157 at 0x14a36c850>, <Reaction adapter_pyr_c_MNXM23 at 0x14a29d4d0>, <Reaction adapter_atp_g_MNXM3 at 0x14a1f05d0>, <Reaction adapter_gdpmann_g_MNXM82 at 0x14a1f0290>, <Reaction adapter_h_g_MNXM1 at 0x14a1f0250>, <Reaction adapter_nadh_m_MNXM10 at 0x14a17f410>, <Reaction adapter_adp_x_MNXM7 at 0x14a10dad0>, <Reaction adapter_pmtcoa_x_MNXM88 at 0x14a0e7d50>, <Reaction adapter_ttdca_x_MNXM314 at 0x14a0e7550>], [<Reaction MNXR3607 at 0x131b16050>, <Reaction MNXR8993 at 0x1473d1950>, <Reaction adapter_ala_dsh_B_c_MNXM144 at 0x14a3921d0>, <Reaction adapter_cer3_26_c_MNXM63157 at 0x14a36c850>, <Reaction adapter_pyr_c_MNXM23 at 0x14a29d4d0>, <Reaction adapter_atp_g_MNXM3 at 0x14a1f05d0>, <Reaction adapter_gdpmann_g_MNXM82 at 0x14a1f0290>, <Reaction adapter_h_g_MNXM1 at 0x14a1f0250>, <Reaction adapter_nadh_m_MNXM10 at 0x14a17f410>, <Reaction adapter_adp_x_MNXM7 at 0x14a10dad0>, <Reaction adapter_pmtcoa_x_MNXM88 at 0x14a0e7d50>, <Reaction adapter_ttdca_x_MNXM314 at 0x14a0e7550>], [<Reaction MNXR5704 at 0x130ff7390>, <Reaction MNXR8993 at 0x1473d1950>, <Reaction adapter_ala_dsh_B_c_MNXM144 at 0x14a3921d0>, <Reaction adapter_cer3_26_c_MNXM63157 at 0x14a36c850>, <Reaction adapter_pyr_c_MNXM23 at 0x14a29d4d0>, <Reaction adapter_atp_g_MNXM3 at 0x14a1f05d0>, <Reaction adapter_gdpmann_g_MNXM82 at 0x14a1f0290>, <Reaction adapter_h_g_MNXM1 at 0x14a1f0250>, <Reaction adapter_nadh_m_MNXM10 at 0x14a17f410>, <Reaction adapter_adp_x_MNXM7 at 0x14a10dad0>, <Reaction adapter_pmtcoa_x_MNXM88 at 0x14a0e7d50>, <Reaction adapter_ttdca_x_MNXM314 at 0x14a0e7550>], [<Reaction MNXR8993 at 0x1473d1950>, <Reaction MNXR8994 at 0x1473d1c10>, <Reaction adapter_ala_dsh_B_c_MNXM144 at 0x14a3921d0>, <Reaction adapter_cer3_26_c_MNXM63157 at 0x14a36c850>, <Reaction adapter_pyr_c_MNXM23 at 0x14a29d4d0>, <Reaction adapter_atp_g_MNXM3 at 0x14a1f05d0>, <Reaction adapter_gdpmann_g_MNXM82 at 0x14a1f0290>, <Reaction adapter_h_g_MNXM1 at 0x14a1f0250>, <Reaction adapter_nadh_m_MNXM10 at 0x14a17f410>, <Reaction adapter_adp_x_MNXM7 at 0x14a10dad0>, <Reaction adapter_pmtcoa_x_MNXM88 at 0x14a0e7d50>, <Reaction adapter_ttdca_x_MNXM314 at 0x14a0e7550>]]


.. parsed-literal::

    /Users/niko/Dev/cameo/cameo/api/products.py:76 [1;31mSettingWithCopyWarning[0m: 
    A value is trying to be set on a copy of a slice from a DataFrame.
    Try using .loc[row_indexer,col_indexer] = value instead
    
    See the the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy


.. code:: python

    report




.. parsed-literal::

    {(<cameo.api.hosts.Host at 0x117484410>,
      <cameo.api.hosts.ModelFacade at 0x117484550>): [],
     (<cameo.api.hosts.Host at 0x117484510>,
      <cameo.api.hosts.ModelFacade at 0x117484610>): []}



IPython notebook
~~~~~~~~~~~~~~~~

Click
`here <http://nbviewer.ipython.org/github/biosustain/cameo/blob/devel/docs/cameo_high_level_interface.ipynb>`__
to download this page as an IPython notebook.
