
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

Loading a model is easy. Just import the ``~cameo.io.load_model``
function.

.. code:: python

    from cameo import load_model



.. raw:: html

    




For example, load the most current genome-scale metabolic reconstruction
of *Escherichia coli*.

.. code:: python

    model = load_model("iJO1366")

Models, reactions, metabolites, etc., provide return HTML when evaluated
in Jupyter notebooks and can thus be easily inspected.

.. code:: python

    model




.. raw:: html

    <table>
    <tr>
    <td>Name</td>
    <td>iJO1366</td>
    </tr>
    <tr>
    <td>Number of metabolites</td>
    <td>1805</td>
    </tr>
    <tr>
    <td>Number of reactions</td>
    <td>2583</td>
    </tr>
    <tr>
    <td>Reactions</td>
    <td><div style="width:100%; max-height:300px; overflow:auto">4crsol_c --> <br>5drib_c --> <br>aacald_c --> <br>amob_c --> <br>mththf_c --> <br>oxam_c --> <br>0.012379 nh4_c + 0.000223 thf_c + 0.499149 ala__L_c + 0.000674 cu2_c + 0.001345 murein3p3p_p + 0.000223 ribflv_c + 0.000223 enter_c + 0.140101 utp_c + 0.000279 accoa_c + 0.000223 fad_c + 0.008151 colipa_e + 0.000223 gthrd_c + 0.214798 pro__L_c + 0.234232 asn__L_c + 0.000223 mql8_c + 2e-06 btn_c + 0.00118 clpn181_p + 0.000112 nadp_c + 0.000223 hemeO_c + 0.025612 dctp_c + 0.154187 glycogen_c + 0.001961 pg181_p + 7e-06 mocogdp_c + 7e-06 mococdp_c + 0.000673 murein4px4px4p_p + 0.003805 pg161_p + 0.006744 spmd_c + 0.129799 ctp_c + 0.282306 ile__L_c + 9.8e-05 succoa_c + 0.000324 zn2_c + 0.005381 murein4p4p_p + 0.005448 murein4px4p_p + 0.000223 5mthf_c + 0.000307 ni2_c + 0.595297 gly_c + 0.004952 cl_c + 0.000223 q8h2_c + 0.000223 mlthf_c + 0.024805 dttp_c + 0.004892 pg160_p + 0.234232 asp__L_c + 48.752916 h2o_c + 0.004126 so4_c + 0.007428 fe3_c + 0.006388 fe2_c + 0.000605 murein3px4p_p + 0.180021 phe__L_c + 0.255712 gln__L_c + 0.000658 mn2_c + 0.024732 pe161_p + 0.000223 chor_c + 0.000223 pydx5p_c + 0.000223 10fthf_c + 0.209121 gtp_c + 0.000223 adocbl_c + 0.012747 pe181_p + 0.002944 clpn160_p + 0.333448 lys__L_c + 0.092056 his__L_c + 0.000116 bmocogdp_c + 0.004952 ca2_c + 0.03327 ptrc_c + 0.18569 k_c + 0.031798 pe160_p + 0.000223 sheme_c + 0.000223 pheme_c + 0.002288 pg181_c + 3e-06 lipopb_c + 0.255712 glu__L_c + 2.4e-05 cobalt2_c + 2.5e-05 2fe2s_c + 0.149336 met__L_c + 0.012366 pe160_c + 0.088988 cys__L_c + 0.000223 thmpp_c + 0.001787 nad_c + 0.000168 coa_c + 0.209684 ser__L_c + 5.5e-05 udcpdp_c + 0.005707 pg160_c + 0.004439 pg161_c + 0.000223 2dmmql8_c + 0.004957 pe181_c + 0.008253 mg2_c + 0.009618 pe161_c + 0.28742 arg__L_c + 54.119975 atp_c + 0.055234 trp__L_c + 3.1e-05 malcoa_c + 0.246506 thr__L_c + 0.024805 datp_c + 4.5e-05 nadh_c + 0.025612 dgtp_c + 0.411184 val__L_c + 0.000335 nadph_c + 0.00229 clpn161_p + 0.000248 4fe4s_c + 0.133993 tyr__L_c + 0.437778 leu__L_c + 0.000223 amet_c + 7e-06 mobd_c --> 0.749831 ppi_c + 53.95 h_c + 53.95 adp_c + 53.945874 pi_c<br>0.000691 mn2_c + 0.000223 pydx5p_c + 0.006715 fe2_c + 0.513689 ala__L_c + 0.027017 dctp_c + 0.137896 tyr__L_c + 0.000223 ribflv_c + 0.144104 utp_c + 0.094738 his__L_c + 0.215096 gtp_c + 0.000223 fad_c + 0.000122 bmocogdp_c + 0.423162 val__L_c + 0.26316 glu__L_c + 2.5e-05 cobalt2_c + 0.09158 cys__L_c + 7e-06 mobd_c + 2e-06 btn_c + 0.195193 k_c + 0.045946 pe160_p + 0.241055 asn__L_c + 0.000223 sheme_c + 0.00026 4fe4s_c + 0.013013 nh4_c + 0.005205 ca2_c + 2.6e-05 2fe2s_c + 0.000447 nadp_c + 0.253687 thr__L_c + 0.133508 ctp_c + 0.001831 nad_c + 0.056843 trp__L_c + 0.215792 ser__L_c + 0.013894 murein5px4p_p + 0.000709 cu2_c + 0.343161 lys__L_c + 0.153686 met__L_c + 0.000576 coa_c + 0.019456 kdo2lipid4_e + 0.290529 ile__L_c + 0.008675 mg2_c + 54.124831 atp_c + 0.054154 pe161_c + 0.000223 2ohph_c + 0.000223 thmpp_c + 0.000223 thf_c + 5.5e-05 udcpdp_c + 0.000323 ni2_c + 0.612638 gly_c + 0.005205 cl_c + 0.000223 mlthf_c + 0.026166 dttp_c + 0.026166 datp_c + 0.241055 asp__L_c + 0.26316 gln__L_c + 0.004338 so4_c + 0.000223 10fthf_c + 0.027017 dgtp_c + 0.007808 fe3_c + 48.601527 h2o_c + 0.000341 zn2_c + 0.000223 pheme_c + 0.185265 phe__L_c + 0.017868 pe160_c + 0.221055 pro__L_c + 0.295792 arg__L_c + 0.02106 pe161_p + 0.450531 leu__L_c + 0.000223 amet_c --> 0.773903 ppi_c + 53.945662 pi_c + 53.95 h_c + 53.95 adp_c<br>12ppd__R_e --> <br>12ppd__S_e --> <br>14glucan_e --> <br>15dap_e --> <br>23camp_e --> <br>23ccmp_e --> <br>23cgmp_e --> <br>23cump_e --> <br>23dappa_e --> <br>26dap__M_e --> <br>2ddglcn_e --> <br>34dhpac_e --> <br>3amp_e --> <br>3cmp_e --> <br>3gmp_e --> <br>3hcinnm_e --> <br>3hpp_e --> <br>3hpppn_e --> <br>3ump_e --> <br>4abut_e --> <br>4hoxpacd_e --> <br>5dglcn_e --> <br>5mtr_e --> <br>LalaDglu_e --> <br>LalaDgluMdap_e --> <br>LalaDgluMdapDala_e --> <br>LalaLglu_e --> <br>ac_e --> <br>acac_e --> <br>acald_e --> <br>acgal_e --> <br>acgal1p_e --> <br>acgam_e --> <br>acgam1p_e --> <br>acmana_e --> <br>acmum_e --> <br>acnam_e --> <br>acolipa_e --> <br>acser_e --> <br>ade_e --> <br>adn_e --> <br>adocbl_e --> <br>ag_e --> <br>agm_e --> <br>akg_e --> <br>ala_B_e --> <br>ala__D_e --> <br>ala__L_e --> <br>alaala_e --> <br>all__D_e --> <br>alltn_e --> <br>amp_e --> <br>anhgm_e --> <br>arab__L_e --> <br>arbt_e --> <br>arbtn_e --> <br>arbtn_fe3_e --> <br>arg__L_e --> <br>ascb__L_e --> <br>asn__L_e --> <br>aso3_e --> <br>asp__L_e --> <br>btn_e --> <br>but_e --> <br>butso3_e --> <br>ca2_e <=> <br>cbi_e --> <br>cbl1_e <=> <br>cd2_e --> <br>cgly_e --> <br>chol_e --> <br>chtbs_e --> <br>cit_e --> <br>cl_e <=> <br>cm_e --> <br>cmp_e --> <br>co2_e <=> <br>cobalt2_e <=> <br>colipa_e --> <br>colipap_e --> <br>cpgn_e --> <br>cpgn_un_e --> <br>crn_e --> <br>crn__D_e --> <br>csn_e --> <br>cu_e --> <br>cu2_e <=> <br>cyan_e --> <br>cynt_e --> <br>cys__D_e --> <br>cys__L_e --> <br>cytd_e --> <br>dad_2_e --> <br>damp_e --> <br>dca_e --> <br>dcmp_e --> <br>dcyt_e --> <br>ddca_e --> <br>dgmp_e --> <br>dgsn_e --> <br>dha_e --> <br>dimp_e --> <br>din_e --> <br>dms_e --> <br>dmso_e --> <br>dopa_e --> <br>doxrbcn_e --> <br>dtmp_e --> <br>dump_e --> <br>duri_e --> <br>eca4colipa_e --> <br>enlipa_e --> <br>enter_e --> <br>etha_e --> <br>ethso3_e --> <br>etoh_e --> <br>f6p_e --> <br>fald_e --> <br>fe2_e <=> <br>fe3_e <=> <br>fe3dcit_e --> <br>fe3dhbzs_e --> <br>fe3hox_e --> <br>fe3hox_un_e --> <br>fecrm_e --> <br>fecrm_un_e --> <br>feenter_e --> <br>feoxam_e --> <br>feoxam_un_e --> <br>for_e --> <br>fru_e --> <br>frulys_e --> <br>fruur_e --> <br>fuc__L_e --> <br>fum_e --> <br>fusa_e --> <br>g1p_e --> <br>g3pc_e --> <br>g3pe_e --> <br>g3pg_e --> <br>g3pi_e --> <br>g3ps_e --> <br>g6p_e --> <br>gal_e --> <br>gal_bD_e --> <br>gal1p_e --> <br>galct__D_e --> <br>galctn__D_e --> <br>galctn__L_e --> <br>galt_e --> <br>galur_e --> <br>gam_e --> <br>gam6p_e --> <br>gbbtn_e --> <br>gdp_e --> <br>glc__D_e <=> <br>glcn_e --> <br>glcr_e --> <br>glcur_e --> <br>glcur1p_e --> <br>gln__L_e --> <br>glu__L_e --> <br>gly_e --> <br>glyald_e --> <br>glyb_e --> <br>glyc_e --> <br>glyc__R_e --> <br>glyc2p_e --> <br>glyc3p_e --> <br>glyclt_e --> <br>gmp_e --> <br>gsn_e --> <br>gthox_e --> <br>gthrd_e --> <br>gtp_e --> <br>gua_e --> <br>h_e <=> <br>h2_e --> <br>h2o_e <=> <br>h2o2_e --> <br>h2s_e --> <br>hacolipa_e --> <br>halipa_e --> <br>hdca_e --> <br>hdcea_e --> <br>hg2_e --> <br>his__L_e --> <br>hom__L_e --> <br>hxa_e --> <br>hxan_e --> <br>idon__L_e --> <br>ile__L_e --> <br>imp_e --> <br>indole_e --> <br>inost_e --> <br>ins_e --> <br>isetac_e --> <br>k_e <=> <br>kdo2lipid4_e --> <br>lac__D_e --> <br>lac__L_e --> <br>lcts_e --> <br>leu__L_e --> <br>lipa_e --> <br>lipa_cold_e --> <br>lipoate_e --> <br>lys__L_e --> <br>lyx__L_e --> <br>mal__D_e --> <br>mal__L_e --> <br>malt_e --> <br>malthx_e --> <br>maltpt_e --> <br>malttr_e --> <br>maltttr_e --> <br>man_e --> <br>man6p_e --> <br>manglyc_e --> <br>melib_e --> <br>meoh_e --> <br>met__D_e --> <br>met__L_e --> <br>metsox_R__L_e --> <br>metsox_S__L_e --> <br>mg2_e <=> <br>mincyc_e --> <br>minohp_e --> <br>mmet_e --> <br>mn2_e <=> <br>mnl_e --> <br>mobd_e <=> <br>mso3_e --> <br>n2o_e --> <br>na1_e <=> <br>nac_e --> <br>nh4_e <=> <br>ni2_e <=> <br>nmn_e --> <br>no_e --> <br>no2_e --> <br>no3_e --> <br>novbcn_e --> <br>o16a4colipa_e --> <br>o2_e <=> <br>o2s_e --> <br>ocdca_e --> <br>ocdcea_e --> <br>octa_e --> <br>orn_e --> <br>orot_e --> <br>pacald_e --> <br>peamn_e --> <br>phe__L_e --> <br>pheme_e --> <br>pi_e <=> <br>pnto__R_e --> <br>ppa_e --> <br>ppal_e --> <br>pppn_e --> <br>ppt_e --> <br>pro__L_e --> <br>progly_e --> <br>psclys_e --> <br>pser__L_e --> <br>ptrc_e --> <br>pydam_e --> <br>pydx_e --> <br>pydxn_e --> <br>pyr_e --> <br>quin_e --> <br>r5p_e --> <br>rfamp_e --> <br>rib__D_e --> <br>rmn_e --> <br>sbt__D_e --> <br>sel_e <=> <br>ser__D_e --> <br>ser__L_e --> <br>skm_e --> <br>slnt_e <=> <br>so2_e --> <br>so3_e --> <br>so4_e <=> <br>spmd_e --> <br>succ_e --> <br>sucr_e --> <br>sulfac_e --> <br>tartr__D_e --> <br>tartr__L_e --> <br>taur_e --> <br>tcynt_e --> <br>thm_e --> <br>thr__L_e --> <br>thrp_e --> <br>thym_e --> <br>thymd_e --> <br>tma_e --> <br>tmao_e --> <br>tre_e --> <br>trp__L_e --> <br>tsul_e --> <br>ttdca_e --> <br>ttdcea_e --> <br>ttrcyc_e --> <br>tungs_e <=> <br>tym_e --> <br>tyr__L_e --> <br>tyrp_e --> <br>uacgam_e --> <br>udpacgal_e --> <br>udpg_e --> <br>udpgal_e --> <br>udpglcur_e --> <br>ump_e --> <br>ura_e --> <br>urea_e --> <br>uri_e --> <br>val__L_e --> <br>xan_e --> <br>xmp_e --> <br>xtsn_e --> <br>xyl__D_e --> <br>xylu__L_e --> <br>zn2_e <=> <br>12dgr120_p --> 12dgr120_c<br>12dgr140_p --> 12dgr140_c<br>12dgr141_p --> 12dgr141_c<br>12dgr160_p --> 12dgr160_c<br>12dgr161_p --> 12dgr161_c<br>12dgr180_p --> 12dgr180_c<br>12dgr181_p --> 12dgr181_c<br>12ppd__R_e <=> 12ppd__R_p<br>12ppd__R_p <=> 12ppd__R_c<br>12ppd__S_e <=> 12ppd__S_p<br>12ppd__S_p <=> 12ppd__S_c<br>14glucan_p + atp_c + h2o_c --> 14glucan_c + pi_c + h_c + adp_c<br>14glucan_e --> 14glucan_p<br>23camp_e <=> 23camp_p<br>23ccmp_e <=> 23ccmp_p<br>23cgmp_e <=> 23cgmp_p<br>23cump_e <=> 23cump_p<br>23dappa_p + h_p --> 23dappa_c + h_c<br>23dappa_e <=> 23dappa_p<br>23cump_p + h2o_p --> h_p + 3ump_p<br>23ccmp_p + h2o_p --> 3cmp_p + h_p<br>h2o_p + 23camp_p --> h_p + 3amp_p<br>23cgmp_p + h2o_p --> h_p + 3gmp_p<br>26dap__M_e <=> 26dap__M_p<br>2ddecg3p_p --> 2ddecg3p_c<br>2tdecg3p_p --> 2tdecg3p_c<br>2tdec7eg3p_p --> 2tdec7eg3p_c<br>2hdecg3p_p --> 2hdecg3p_c<br>2hdec9eg3p_p --> 2hdec9eg3p_c<br>2odecg3p_p --> 2odecg3p_c<br>2odec11eg3p_p --> 2odec11eg3p_c<br>2agpe120_p --> 2agpe120_c<br>2agpe140_p --> 2agpe140_c<br>2agpe141_p --> 2agpe141_c<br>2agpe160_p --> 2agpe160_c<br>2agpe161_p --> 2agpe161_c<br>2agpe180_p --> 2agpe180_c<br>2agpe181_p --> 2agpe181_c<br>ddca_c + 2agpe120_c + atp_c --> ppi_c + pe120_c + amp_c<br>2agpe140_c + ttdca_c + atp_c --> pe140_c + ppi_c + amp_c<br>ttdcea_c + 2agpe141_c + atp_c --> pe141_c + ppi_c + amp_c<br>hdca_c + 2agpe160_c + atp_c --> pe160_c + ppi_c + amp_c<br>hdcea_c + 2agpe161_c + atp_c --> ppi_c + pe161_c + amp_c<br>2agpe180_c + ocdca_c + atp_c --> ppi_c + pe180_c + amp_c<br>2agpe181_c + ocdcea_c + atp_c --> ppi_c + pe181_c + amp_c<br>2agpg120_p --> 2agpg120_c<br>2agpg140_p --> 2agpg140_c<br>2agpg141_p --> 2agpg141_c<br>2agpg160_p --> 2agpg160_c<br>2agpg161_p --> 2agpg161_c<br>2agpg180_p --> 2agpg180_c<br>2agpg181_p --> 2agpg181_c<br>ddca_c + 2agpg120_c + atp_c --> pg120_c + ppi_c + amp_c<br>2agpg140_c + ttdca_c + atp_c --> ppi_c + amp_c + pg140_c<br>ttdcea_c + 2agpg141_c + atp_c --> ppi_c + pg141_c + amp_c<br>hdca_c + 2agpg160_c + atp_c --> pg160_c + ppi_c + amp_c<br>2agpg161_c + hdcea_c + atp_c --> pg161_c + ppi_c + amp_c<br>ocdca_c + 2agpg180_c + atp_c --> pg180_c + ppi_c + amp_c<br>2agpg181_c + ocdcea_c + atp_c --> ppi_c + pg181_c + amp_c<br>2dhguln_c + nadh_c + h_c --> glcn_c + nad_c<br>2dhguln_c + nadph_c + h_c --> glcn_c + nadp_c<br>2dhguln_c + nadh_c + h_c --> nad_c + idon__L_c<br>2dhguln_c + nadph_c + h_c --> nadp_c + idon__L_c<br>2mahmp_c + h2o_c --> pi_c + 4ampm_c + h_c<br>34dhpac_e <=> 34dhpac_p<br>3amac_c + h2o_c + h_c --> nh4_c + msa_c<br>3amp_e <=> 3amp_p<br>3cmp_e <=> 3cmp_p<br>3gmp_e <=> 3gmp_p<br>3hdecACP_c --> h2o_c + tdec2eACP_c<br>3hddecACP_c --> h2o_c + tddec2eACP_c<br>3hcddec5eACP_c --> t3c5ddeceACP_c + h2o_c<br>3hmrsACP_c --> h2o_c + tmrs2eACP_c<br>3hcmrs7eACP_c --> t3c7mrseACP_c + h2o_c<br>3hpalmACP_c --> h2o_c + tpalm2eACP_c<br>3hcpalm9eACP_c --> h2o_c + t3c9palmeACP_c<br>3hoctaACP_c --> toctd2eACP_c + h2o_c<br>3hcvac11eACP_c --> t3c11vaceACP_c + h2o_c<br>3haACP_c --> but2eACP_c + h2o_c<br>3hhexACP_c --> thex2eACP_c + h2o_c<br>3hoctACP_c --> toct2eACP_c + h2o_c<br>3hcinnm_c + h_c + nadh_c + o2_c --> dhcinnm_c + nad_c + h2o_c<br>3hpppn_c + nadh_c + h_c + o2_c --> nad_c + dhpppn_c + h2o_c<br>3hpp_e <=> 3hpp_p<br>3hpp_c + h_c --> 3hpp_p + h_p<br>3dhguln_c + atp_c --> adp_c + 3dhgulnp_c + h_c<br>3ump_p + h2o_p --> pi_p + uri_p<br>3cmp_p + h2o_p --> cytd_p + pi_p<br>h2o_p + 3amp_p --> pi_p + adn_p<br>h2o_p + 3gmp_p --> gsn_p + pi_p<br>3odecACP_c + nadph_c + h_c <=> nadp_c + 3hdecACP_c<br>3oddecACP_c + h_c + nadph_c <=> nadp_c + 3hddecACP_c<br>nadph_c + 3ocddec5eACP_c + h_c --> nadp_c + 3hcddec5eACP_c<br>nadph_c + 3omrsACP_c + h_c <=> 3hmrsACP_c + nadp_c<br>3ocmrs7eACP_c + nadph_c + h_c --> 3hcmrs7eACP_c + nadp_c<br>3opalmACP_c + nadph_c + h_c <=> nadp_c + 3hpalmACP_c<br>nadph_c + h_c + 3ocpalm9eACP_c --> nadp_c + 3hcpalm9eACP_c<br>nadph_c + 3ooctdACP_c + h_c <=> nadp_c + 3hoctaACP_c<br>nadph_c + 3ocvac11eACP_c + h_c --> 3hcvac11eACP_c + nadp_c<br>actACP_c + nadph_c + h_c <=> nadp_c + 3haACP_c<br>nadph_c + h_c + 3ohexACP_c <=> nadp_c + 3hhexACP_c<br>nadph_c + 3ooctACP_c + h_c <=> nadp_c + 3hoctACP_c<br>ocACP_c + h_c + malACP_c --> co2_c + 3odecACP_c + ACP_c<br>h_c + dcaACP_c + malACP_c --> co2_c + 3oddecACP_c + ACP_c<br>h_c + malACP_c + cdec3eACP_c --> co2_c + ACP_c + 3ocddec5eACP_c<br>ddcaACP_c + h_c + malACP_c --> 3omrsACP_c + co2_c + ACP_c<br>cddec5eACP_c + h_c + malACP_c --> co2_c + 3ocmrs7eACP_c + ACP_c<br>malACP_c + myrsACP_c + h_c --> co2_c + 3opalmACP_c + ACP_c<br>tdeACP_c + h_c + malACP_c --> co2_c + ACP_c + 3ocpalm9eACP_c<br>palmACP_c + h_c + malACP_c --> co2_c + 3ooctdACP_c + ACP_c<br>h_c + hdeACP_c + malACP_c --> co2_c + ACP_c + 3ocvac11eACP_c<br>h_c + butACP_c + malACP_c --> 3ohexACP_c + co2_c + ACP_c<br>h_c + hexACP_c + malACP_c --> co2_c + ACP_c + 3ooctACP_c<br>coa_c + oxadpcoa_c --> accoa_c + succoa_c<br>LalaDgluMdap_p + h2o_c + atp_c --> h_c + pi_c + LalaDgluMdap_c + adp_c<br>LalaDgluMdap_e <=> LalaDgluMdap_p<br>3ump_e <=> 3ump_p<br>dopa_p + h2o_p + o2_p --> nh4_p + 34dhpac_p + h2o2_p<br>4hoxpacd_e <=> 4hoxpacd_p<br>h2o_c + phthr_c --> 4hthr_c + pi_c<br>h2o_c + LalaDgluMdapDala_c --> ala__D_c + LalaDgluMdap_c<br>LalaDgluMdapDala_p + h2o_p --> LalaDgluMdap_p + ala__D_p<br>LalaDgluMdapDala_p + atp_c + h2o_c --> LalaDgluMdapDala_c + pi_c + h_c + adp_c<br>LalaDgluMdapDala_e <=> LalaDgluMdapDala_p<br>nadph_c + 5dglcn_c + h_c <=> glcn_c + nadp_c<br>5dglcn_p + h_p <=> 5dglcn_c + h_c<br>5dglcn_e <=> 5dglcn_p<br>dad_5_c + h2o_c --> 5drib_c + ade_c<br>5mtr_e <=> 5mtr_p<br>5mtr_c + h_c --> h_p + 5mtr_p<br>ru5p__D_c <=> ara5p_c<br>ACP_c + ttdca_c + atp_c --> ppi_c + myrsACP_c + amp_c<br>ttdcea_c + ACP_c + atp_c --> tdeACP_c + ppi_c + amp_c<br>hdca_c + ACP_c + atp_c --> palmACP_c + ppi_c + amp_c<br>hdcea_c + ACP_c + atp_c --> ppi_c + hdeACP_c + amp_c<br>ACP_c + ocdcea_c + atp_c --> ppi_c + octeACP_c + amp_c<br>ACP_c + ocdca_c + atp_c --> ocdcaACP_c + ppi_c + amp_c<br>ddca_c + ACP_c + atp_c --> ddcaACP_c + ppi_c + amp_c<br>ACP_c + dca_c + atp_c --> ppi_c + dcaACP_c + amp_c<br>ACP_c + atp_c + octa_c --> ocACP_c + ppi_c + amp_c<br>aact_c + o2_c + h2o_c --> mthgxl_c + h2o2_c + nh4_c<br>unagamu_c + dtdp4aaddg_c --> unagamuf_c + dtdp_c + h_c<br>14glucan_c --> malthx_c<br>14glucan_p --> malthx_p<br>arbt6p_c + h2o_c --> g6p_c + hqn_c<br>4abut_c + akg_c --> glu__L_c + sucsal_c<br>nad_c + 4abutn_c + h2o_c --> 4abut_c + 2.0 h_c + nadh_c<br>4abut_p + h_p --> 4abut_c + h_c<br>4abut_e <=> 4abut_p<br>accoa_c + acac_c --> ac_c + aacoa_c<br>2.0 accoa_c <=> coa_c + aacoa_c<br>accoa_c + btcoa_c <=> coa_c + 3ohcoa_c<br>hxcoa_c + accoa_c <=> coa_c + 3oocoa_c<br>accoa_c + occoa_c <=> 3odcoa_c + coa_c<br>dcacoa_c + accoa_c <=> coa_c + 3oddcoa_c<br>accoa_c + ddcacoa_c <=> coa_c + 3otdcoa_c<br>tdcoa_c + accoa_c <=> 3ohdcoa_c + coa_c<br>3oodcoa_c + coa_c <=> accoa_c + pmtcoa_c<br>h_p + acac_p <=> acac_c + h_c<br>acac_e <=> acac_p<br>acald_c + coa_c + nad_c <=> h_c + nadh_c + accoa_c<br>acald_e <=> acald_p<br>acald_p <=> acald_c<br>anth_c + accoa_c --> coa_c + acanth_c<br>adocbip_c + gtp_c + h_c --> ppi_c + agdpcbi_c<br>hco3_c + accoa_c + atp_c --> malcoa_c + h_c + pi_c + adp_c<br>ppa_c + coa_c + atp_c --> ppcoa_c + pi_c + adp_c<br>acgal1p_p + h2o_p --> acgal_p + pi_p<br>acgal1p_e <=> acgal1p_p<br>acgal_e <=> acgal_p<br>acgam1p_p + h2o_p --> pi_p + acgam_p<br>acgam1p_e <=> acgam1p_p<br>atp_c + acgam_c --> acgam6p_c + adp_c + h_c<br>udcpp_c + uacgam_c --> ump_c + unaga_c<br>pep_c + acgam_p --> acgam6p_c + pyr_c<br>acgam_e <=> acgam_p<br>acglu_c + atp_c --> acg5p_c + adp_c<br>glu__L_c + accoa_c --> coa_c + acglu_c + h_c<br>2obut_c + pyr_c + h_c --> 2ahbut_c + co2_c<br>ac_c + atp_c <=> actp_c + adp_c<br>2.0 pyr_c + h_c --> co2_c + alac__S_c<br>h2o_c + acmum6p_c --> acgam6p_c + lac__D_c<br>uacmamu_c + unaga_c --> unagamu_c + udp_c + h_c<br>pep_c + acmana_p --> acmanap_c + pyr_c<br>acmana_e <=> acmana_p<br>pep_c + acmum_p --> pyr_c + acmum6p_c<br>acmum_e --> acmum_p<br>h_p + acnam_p --> acnam_c + h_c<br>acnam_e <=> acnam_p<br>acnam_c --> acmana_c + pyr_c<br>fad_c + btcoa_c <=> b2coa_c + fadh2_c<br>hxcoa_c + fad_c <=> hx2coa_c + fadh2_c<br>fad_c + occoa_c <=> fadh2_c + oc2coa_c<br>dcacoa_c + fad_c <=> dc2coa_c + fadh2_c<br>fad_c + ddcacoa_c <=> fadh2_c + dd2coa_c<br>tdcoa_c + fad_c <=> fadh2_c + td2coa_c<br>fad_c + pmtcoa_c <=> hdd2coa_c + fadh2_c<br>stcoa_c + fad_c <=> fadh2_c + od2coa_c<br>ACP_c + accoa_c <=> coa_c + acACP_c<br>acorn_c + h2o_c --> orn_c + ac_c<br>acolipa_p + atp_c + h2o_c --> acolipa_e + h_c + pi_c + adp_c<br>acon_T_c <=> acon_C_c<br>acon_T_c + amet_c --> ahcys_c + aconm_c<br>cit_c <=> acon_C_c + h2o_c<br>acon_C_c + h2o_c <=> icit_c<br>akg_c + acorn_c <=> acg5sa_c + glu__L_c<br>pi_c + ddcaACP_c + h_c --> ACP_c + ddcap_c<br>myrsACP_c + pi_c + h_c --> ACP_c + ttdcap_c<br>pi_c + tdeACP_c + h_c --> ACP_c + ttdceap_c<br>palmACP_c + pi_c + h_c --> ACP_c + hdcap_c<br>pi_c + h_c + hdeACP_c --> ACP_c + hdceap_c<br>pi_c + ocdcaACP_c + h_c --> ACP_c + ocdcap_c<br>pi_c + octeACP_c + h_c --> ocdceap_c + ACP_c<br>apoACP_c + coa_c --> pap_c + ACP_c + h_c<br>coa_c + ac_c + atp_c --> ppi_c + accoa_c + amp_c<br>acser_e <=> acser_p<br>acser_c --> acser_p<br>h_p + ac_p <=> ac_c + h_c<br>ac_p + na1_p --> ac_c + na1_c<br>ac_e <=> ac_p<br>adn_c + h_c + h2o_c --> ins_c + nh4_c<br>4adcho_c --> pyr_c + 4abz_c + h_c<br>chor_c + gln__L_c --> 4adcho_c + glu__L_c<br>ade_c + h2o_c + h_c --> hxan_c + nh4_c<br>h_p + ade_p <=> ade_c + h_c<br>ade_e <=> ade_p<br>amp_c + atp_c <=> 2.0 adp_c<br>amp_c + gtp_c <=> gdp_c + adp_c<br>itp_c + amp_c <=> idp_c + adp_c<br>amet_c + h_c --> ametam_c + co2_c<br>atp_c --> ppi_c + camp_c<br>adn_c + atp_c --> amp_c + adp_c + h_c<br>adn_c + h2o_c --> rib__D_c + ade_c<br>h_p + adn_p --> adn_c + h_c<br>h_p + adn_p <=> adn_c + h_c<br>adn_e <=> adn_p<br>adocbi_c + atp_c --> adp_c + adocbip_c + h_c<br>rdmbzi_c + agdpcbi_c --> adocbl_c + gmp_c + h_c<br>adocbl_p + atp_c + h2o_c --> h_c + adocbl_c + pi_c + adp_c<br>adocbl_e + h_p --> adocbl_p + h_c<br>adprib_c + h2o_c --> amp_c + 2.0 h_c + r5p_c<br>ade_c + prpp_c --> ppi_c + amp_c<br>aps_c + atp_c --> paps_c + adp_c + h_c<br>dcamp_c <=> fum_c + amp_c<br>25aics_c <=> fum_c + aicar_c<br>imp_c + asp__L_c + gtp_c --> gdp_c + 2.0 h_c + pi_c + dcamp_c<br>acgam6p_c + h2o_c --> gam6p_c + ac_c<br>anhgm3p_c + h2o_c --> anhgm_c + LalaDgluMdap_c<br>anhgm3p_p + h2o_p --> LalaDgluMdap_p + anhgm_p<br>anhgm3p_c + h2o_c --> anhm3p_c + acgam_c<br>h_p + anhgm3p_p --> anhgm3p_c + h_c<br>h2o_c + anhgm4p_c --> anhgm_c + LalaDgluMdapDala_c<br>h2o_p + anhgm4p_p --> LalaDgluMdapDala_p + anhgm_p<br>h2o_c + anhgm4p_c --> anhgm3p_c + ala__D_c<br>h2o_p + anhgm4p_p --> anhgm3p_p + ala__D_p<br>anhgm4p_c + h2o_c --> anhm4p_c + acgam_c<br>h_p + anhgm4p_p --> anhgm4p_c + h_c<br>anhgm_c + h2o_c --> anhm_c + acgam_c<br>adphep_DD_c --> adphep_LD_c<br>agm_c + h2o_c --> urea_c + ptrc_c<br>h_p + anhgm_p --> anhgm_c + h_c<br>agm_e <=> agm_p<br>ddcaACP_c + 1ddecg3p_c --> ACP_c + pa120_c<br>myrsACP_c + 1tdecg3p_c --> ACP_c + pa140_c<br>1tdec7eg3p_c + tdeACP_c --> ACP_c + pa141_c<br>palmACP_c + 1hdecg3p_c --> pa160_c + ACP_c<br>hdeACP_c + 1hdec9eg3p_c --> ACP_c + pa161_c<br>ocdcaACP_c + 1odecg3p_c --> pa180_c + ACP_c<br>1odec11eg3p_c + octeACP_c --> pa181_c + ACP_c<br>acg5sa_c + pi_c + nadp_c <=> nadph_c + acg5p_c + h_c<br>h_e + ag_c --> ag_e + h_c<br>ahcys_c + h2o_c --> rhcys_c + ade_c<br>10fthf_c + aicar_c <=> fprica_c + thf_c<br>air_c + hco3_c + atp_c --> h_c + pi_c + 5caiz_c + adp_c<br>5aizc_c <=> 5caiz_c<br>akg_c + coa_c + nad_c --> co2_c + nadh_c + succoa_c<br>h_p + akg_p <=> akg_c + h_c<br>akg_e <=> akg_p<br>alaala_c + h2o_c --> 2.0 ala__D_c<br>alaala_p + atp_c + h2o_c --> h_c + pi_c + alaala_c + adp_c<br>2.0 ala__D_c + atp_c <=> h_c + pi_c + alaala_c + adp_c<br>alaala_e <=> alaala_p<br>LalaDglu_c <=> LalaLglu_c<br>ala__L_c <=> ala__D_c<br>pydx5p_c + ala__D_c --> pyr_c + pyam5p_c<br>ala__L_c + akg_c <=> pyr_c + glu__L_c<br>ala__L_c + pydx5p_c --> pyr_c + pyam5p_c<br>trnaala_c + ala__L_c + atp_c --> ppi_c + amp_c + alatrna_c<br>ala__L_p + atp_c + h2o_c --> h_c + ala__L_c + pi_c + adp_c<br>ala__L_p + h_p --> ala__L_c + h_c<br>ala__L_p + h_p <=> ala__L_c + h_c<br>ala__L_p + na1_p --> ala__L_c + na1_c<br>ala__L_e <=> ala__L_p<br>nadh_c + glyald_c + h_c <=> nad_c + glyc_c<br>etoh_c + nad_c <=> nadh_c + acald_c + h_c<br>nad_c + pacald_c + h2o_c <=> pac_c + 2.0 h_c + nadh_c<br>acald_c + nad_c + h2o_c --> nadh_c + 2.0 h_c + ac_c<br>acald_c + nadp_c + h2o_c --> nadph_c + 2.0 h_c + ac_c<br>nadp_c + ppal_c + h2o_c --> 2.0 h_c + ppa_c + nadph_c<br>nad_c + btal_c + h2o_c --> but_c + 2.0 h_c + nadh_c<br>all__D_c + atp_c --> adp_c + all6p_c + h_c<br>all6p_c <=> allul6p_c<br>alltt_c + 2.0 h_c + 2.0 h2o_c --> co2_c + urdglyc_c + 2.0 nh4_c<br>alltn_c + h2o_c --> alltt_c + h_c<br>h_p + alltn_p <=> alltn_c + h_c<br>alltn_e <=> alltn_p<br>allul6p_c <=> f6p_c<br>all__D_p + atp_c + h2o_c --> all__D_c + h_c + pi_c + adp_c<br>all__D_e <=> all__D_p<br>pe160_p + alpp_p --> lpp_p + 2agpe160_p<br>alpp_p + pg160_p --> lpp_p + 2agpg160_p<br>mthgxl_c + nadph_c + h_c --> acetol_c + nadp_c<br>h_c + acetol_c + nadh_c --> nad_c + 12ppd__R_c<br>altrn_c --> 2ddglcn_c + h2o_c<br>anhm3p_c + h2o_c --> LalaDgluMdap_c + anhm_c<br>anhm4p_c + h2o_c --> anhm_c + LalaDgluMdapDala_c<br>anhm4p_c + h2o_c --> anhm3p_c + ala__D_c<br>malttr_c + malt_c --> glc__D_c + maltttr_c<br>malt_c + maltttr_c --> glc__D_c + maltpt_c<br>malt_c + maltpt_c --> glc__D_c + malthx_c<br>malt_c + malthx_c --> glc__D_c + malthp_c<br>acmanap_c <=> acgam6p_c<br>acmana_c + atp_c --> acmanap_c + adp_c + h_c<br>8aonn_c + amet_c <=> amob_c + dann_c<br>2dmmql8_c + amet_c --> ahcys_c + mql8_c + h_c<br>air_c + nad_c + h2o_c --> nadh_c + 2.0 for_c + 3.0 h_c + 4ampm_c<br>amp_c + h2o_c --> ade_c + r5p_c<br>h2o_c + cgly_c --> gly_c + cys__L_c<br>h2o_c + progly_c --> pro__L_c + gly_c<br>amp_e <=> amp_p<br>anhgm_e <=> anhgm_p<br>h2o_c + anhm_c + atp_c --> adp_c + h_c + acmum6p_c<br>anth_c + prpp_c --> pran_c + ppi_c<br>chor_c + gln__L_c --> anth_c + h_c + pyr_c + glu__L_c<br>2aobut_c + h_c --> aact_c + co2_c<br>ala__L_c + pimACP_c --> ACP_c + 8aonn_c + co2_c<br>ap4a_c + h2o_c --> 2.0 h_c + 2.0 adp_c<br>2.0 atp_c + h_c --> ap4a_c + ppi_c<br>ap5a_c + h2o_c --> adp_c + atp_c + 2.0 h_c<br>ametam_c + 15dap_c --> 5mta_c + na15dap_c + h_c<br>glyc3p_c + ddcap_c --> pi_c + 1ddecg3p_c + h_c<br>glyc3p_c + ttdcap_c --> pi_c + 1tdecg3p_c + h_c<br>glyc3p_c + ttdceap_c --> 1tdec7eg3p_c + pi_c + h_c<br>glyc3p_c + hdcap_c --> 1hdecg3p_c + pi_c + h_c<br>glyc3p_c + hdceap_c --> h_c + 1hdec9eg3p_c + pi_c<br>glyc3p_c + ocdcap_c --> pi_c + 1odecg3p_c + h_c<br>ocdceap_c + glyc3p_c --> 1odec11eg3p_c + pi_c + h_c<br>h2o_c + ddcap_c --> pi_c + ddca_c + 2.0 h_c<br>h2o_c + ttdcap_c --> pi_c + ttdca_c + 2.0 h_c<br>ttdceap_c + h2o_c --> ttdcea_c + pi_c + 2.0 h_c<br>hdcap_c + h2o_c --> pi_c + 2.0 h_c + hdca_c<br>h2o_c + hdceap_c --> pi_c + 2.0 h_c + hdcea_c<br>h2o_c + ocdcap_c --> ocdca_c + pi_c + 2.0 h_c<br>ocdceap_c + h2o_c --> pi_c + ocdcea_c + 2.0 h_c<br>h_c + nadh_c + aact_c <=> nad_c + appl_c<br>nadph_c + 5apru_c + h_c --> 5aprbu_c + nadp_c<br>arab__L_c <=> rbl__L_c<br>2.0 arbtn_fe3_c + fadh2_c --> 2.0 arbtn_c + 2.0 h_c + fad_c + 2.0 fe2_c<br>2.0 arbtn_fe3_c + fmnh2_c --> 2.0 arbtn_c + 2.0 h_c + 2.0 fe2_c + fmn_c<br>2.0 arbtn_fe3_c + rbflvrd_c --> 2.0 arbtn_c + 2.0 h_c + 2.0 fe2_c + ribflv_c<br>arbtn_fe3_p + atp_c + h2o_c --> arbtn_fe3_c + h_c + pi_c + adp_c<br>arbtn_e + fe3_e --> arbtn_fe3_e<br>arbtn_p + h_p --> arbtn_e + h_c<br>h_p + arbtn_fe3_e --> arbtn_fe3_p + h_c<br>arbtn_c + h_p --> arbtn_p + h_c<br>pep_c + arbt_p --> pyr_c + arbt6p_c<br>arbt_e --> arbt_p<br>arab__L_p + atp_c + h2o_c --> h_c + pi_c + arab__L_c + adp_c<br>arab__L_p + h_p <=> arab__L_c + h_c<br>arab__L_c + h_p --> arab__L_p + h_c<br>arab__L_e <=> arab__L_p<br>agm_c + arg__L_p <=> arg__L_c + agm_p<br>arg__L_c + h_c --> agm_c + co2_c<br>h_p + arg__L_p --> agm_p + co2_p<br>orn_c + arg__L_p <=> arg__L_c + orn_p<br>argsuc_c <=> arg__L_c + fum_c<br>asp__L_c + atp_c + citr__L_c --> argsuc_c + ppi_c + h_c + amp_c<br>arg__L_c + trnaarg_c + atp_c --> ppi_c + argtrna_c + amp_c<br>arg__L_p + atp_c + h2o_c --> arg__L_c + pi_c + h_c + adp_c<br>arg__L_c + h_p --> h_c + arg__L_p<br>arg__L_e <=> arg__L_p<br>pi_c + aspsa_c + nadp_c <=> nadph_c + h_c + 4pasp_c<br>h2o_c + ascb6p_c --> h_c + 3dhgulnp_c<br>ascb__L_p + pep_c --> pyr_c + ascb6p_c<br>ascb__L_e <=> ascb__L_p<br>asn__L_c + h2o_c --> nh4_c + asp__L_c<br>asn__L_p + h2o_p --> nh4_p + asp__L_p<br>atp_c + asp__L_c + h2o_c + gln__L_c --> asn__L_c + h_c + ppi_c + glu__L_c + amp_c<br>asp__L_c + nh4_c + atp_c --> asn__L_c + ppi_c + h_c + amp_c<br>asn__L_c + trnaasn_c + atp_c --> ppi_c + asntrna_c + amp_c<br>asn__L_p + atp_c + h2o_c --> asn__L_c + h_c + adp_c + pi_c<br>h_p + asn__L_p <=> asn__L_c + h_c<br>asn__L_e <=> asn__L_p<br>aso3_c + atp_c + h2o_c --> h_c + aso3_p + pi_c + adp_c<br>aso3_e <=> aso3_p<br>asp__L_c + h_c --> co2_c + ala_B_c<br>cbp_c + asp__L_c --> pi_c + cbasp_c + h_c<br>atp_c + asp__L_c <=> 4pasp_c + adp_c<br>q8_c + asp__L_c --> iasp_c + h_c + q8h2_c<br>mqn8_c + asp__L_c --> mql8_c + iasp_c + h_c<br>fum_c + asp__L_c --> iasp_c + succ_c + h_c<br>o2_c + asp__L_c --> h2o2_c + iasp_c + h_c<br>asp__L_c --> nh4_c + fum_c<br>akg_c + asp__L_c <=> oaa_c + glu__L_c<br>asp__L_c + trnaasp_c + atp_c --> asptrna_c + ppi_c + amp_c<br>asp__L_p + atp_c + h2o_c --> asp__L_c + h_c + pi_c + adp_c<br>2.0 h_p + asp__L_p --> asp__L_c + 2.0 h_c<br>3.0 h_p + asp__L_p --> asp__L_c + 3.0 h_c<br>h_p + asp__L_p --> asp__L_c + h_c<br>h_p + asp__L_p <=> asp__L_c + h_c<br>asp__L_e <=> asp__L_p<br>2.0 gthrd_c + aso4_c --> aso3_c + gthox_c + h2o_c<br>arg__L_c + succoa_c --> sucarg_c + coa_c + h_c<br>nadp_c + athr__L_c <=> nadph_c + 2aobut_c + h_c<br>h2o_c + atp_c + h_c --> itp_c + nh4_c<br>h2o_c + atp_c --> pi_c + adp_c + h_c<br>prpp_c + atp_c --> ppi_c + prbatp_c<br>pi_c + 4.0 h_p + adp_c <=> 3.0 h_c + h2o_c + atp_c<br>ala_B_p + h_p --> ala_B_c + h_c<br>ala_B_e <=> ala_B_p<br>betald_c + nad_c + h2o_c --> glyb_c + 2.0 h_c + nadh_c<br>nadp_c + betald_c + h2o_c --> nadph_c + glyb_c + 2.0 h_c<br>mptamp_c + moco_c --> amp_c + cu2_c + bmoco_c<br>bmoco_c + gtp_c + h_c --> ppi_c + bmoco1gdp_c<br>gtp_c + bmoco1gdp_c + h_c --> bmocogdp_c + ppi_c<br>pap_c + h2o_c --> pi_c + amp_c<br>h_c + nadh_c + btnso_c --> btn_c + nad_c + h2o_c<br>nadph_c + h_c + btnso_c --> btn_c + h2o_c + nadp_c<br>h_p + btn_p --> btn_c + h_c<br>btn_e <=> btn_p<br>2fe2s_c + dtbt_c + amet_c --> dad_5_c + 2fe1s_c + met__L_c + h_c + btn_c<br>but_c + accoa_c --> ac_c + btcoa_c<br>butso3_p + atp_c + h2o_c --> h_c + butso3_c + adp_c + pi_c<br>butso3_e <=> butso3_p<br>but_p + h_p <=> but_c + h_c<br>but_e <=> but_p<br>bwco_c + gtp_c + h_c --> ppi_c + bwco1gdp_c<br>bwco1gdp_c + gtp_c + h_c --> bwcogdp_c + ppi_c<br>mptamp_c + wco_c --> cu2_c + bwco_c + amp_c<br>h_p + ca2_c --> ca2_p + h_c<br>ca2_e <=> ca2_p<br>15dap_c + lys__L_p + h_p --> lys__L_c + h_c + 15dap_p<br>2.0 h2o2_c --> o2_c + 2.0 h2o_c<br>ca2_c + na1_p <=> ca2_p + na1_c<br>cbi_c + atp_c + h_c <=> pppi_c + adocbi_c<br>cbi_e + h_p --> cbi_p + h_c<br>cbi_p + atp_c + h2o_c --> cbi_c + h_c + pi_c + adp_c<br>cbl1_p + atp_c + h2o_c --> cbl1_c + h_c + pi_c + adp_c<br>h_p + cbl1_e --> cbl1_p + h_c<br>cbl1_c + atp_c + h_c <=> pppi_c + adocbl_c<br>cbm_c + 2.0 h_c --> co2_c + nh4_c<br>co2_c + nh4_c + atp_c <=> 2.0 h_c + cbp_c + adp_c<br>hco3_c + h2o_c + 2.0 atp_c + gln__L_c --> 2.0 h_c + cbp_c + pi_c + glu__L_c + 2.0 adp_c<br>cdg_c + nh4_c + atp_c --> h_c + h2o_c + preq0_c + pi_c + adp_c<br>cd2_c + atp_c + h2o_c --> cd2_p + h_c + pi_c + adp_c<br>h_p + cd2_c --> cd2_p + h_c<br>cd2_e <=> cd2_p<br>cd2_p --> cd2_c<br>h2o_c + cdpdddecg_c --> cmp_c + pa120_c + 2.0 h_c<br>cdpdtdecg_c + h2o_c --> pa140_c + cmp_c + 2.0 h_c<br>h2o_c + cdpdtdec7eg_c --> 2.0 h_c + pa141_c + cmp_c<br>h2o_c + cdpdhdecg_c --> pa160_c + cmp_c + 2.0 h_c<br>h2o_c + cdpdhdec9eg_c --> 2.0 h_c + cmp_c + pa161_c<br>cdpdodecg_c + h2o_c --> pa180_c + cmp_c + 2.0 h_c<br>h2o_c + cdpdodec11eg_c --> pa181_c + cmp_c + 2.0 h_c<br>preq0_c + 2.0 nadph_c + 3.0 h_c --> 2.0 nadp_c + preq1_c<br>cph4_c + h_c --> cdg_c + nh4_c<br>4c2me_c + atp_c --> 2p4c2me_c + adp_c + h_c<br>2.0 amet_c + pe161_c --> 2.0 ahcys_c + 2.0 h_c + cpe160_c<br>2.0 amet_c + pg161_c --> cpg160_c + 2.0 ahcys_c + 2.0 h_c<br>pe181_c + 2.0 amet_c --> 2.0 ahcys_c + 2.0 h_c + cpe180_c<br>pg181_c + 2.0 amet_c --> 2.0 ahcys_c + cpg180_c + 2.0 h_c<br>h2o_c + cgly_p + atp_c --> h_c + pi_c + adp_c + cgly_c<br>cgly_e <=> cgly_p<br>chol_p + atp_c + h2o_c --> chol_c + h_c + pi_c + adp_c<br>h_p + chol_p --> chol_c + h_c<br>chol_e <=> chol_p<br>nad_c + chol_c --> nadh_c + betald_c + h_c<br>chor_c --> pphn_c<br>3psme_c --> pi_c + chor_c<br>chor_c --> pyr_c + 4hbz_c<br>pep_c + chtbs_p --> pyr_c + chtbs6p_c<br>chtbs_e <=> chtbs_p<br>cinnm_c + h_c + nadh_c + o2_c --> cenchddd_c + nad_c<br>cit_c --> oaa_c + ac_c<br>h_p + cit_c --> cit_p + h_c<br>cit_p + succ_c --> succ_p + cit_c<br>cit_e <=> cit_p<br>lipa_cold_p + atp_c + h2o_c --> lipa_cold_e + h_c + pi_c + adp_c<br>clpn120_p + h2o_p --> pa120_p + pg120_p + h_p<br>h2o_p + clpn140_p --> h_p + pa140_p + pg140_p<br>clpn141_p + h2o_p --> pg141_p + h_p + pa141_p<br>clpn160_p + h2o_p --> pa160_p + h_p + pg160_p<br>clpn161_p + h2o_p --> pg161_p + h_p + pa161_p<br>clpn180_p + h2o_p --> pa180_p + h_p + pg180_p<br>clpn181_p + h2o_p --> h_p + pa181_p + pg181_p<br>2.0 pg120_p <=> clpn120_p + glyc_p<br>2.0 pg140_p <=> glyc_p + clpn140_p<br>2.0 pg141_p <=> clpn141_p + glyc_p<br>2.0 pg160_p <=> clpn160_p + glyc_p<br>2.0 pg161_p <=> glyc_p + clpn161_p<br>2.0 pg180_p <=> clpn180_p + glyc_p<br>2.0 pg181_p <=> clpn181_p + glyc_p<br>2.0 cl_p + h_c --> 2.0 cl_c + h_p<br>cl_e <=> cl_p<br>cmp_c + h2o_c --> csn_c + r5p_c<br>cmp_e <=> cmp_p<br>cm_e <=> cm_p<br>cm_p + h_p --> cm_e + h_c<br>co2_e <=> co2_p<br>co2_p <=> co2_c<br>cobalt2_c + atp_c + h2o_c --> cobalt2_p + h_c + pi_c + adp_c<br>h_p + cobalt2_c --> cobalt2_p + h_c<br>cobalt2_e <=> cobalt2_p<br>cobalt2_p --> cobalt2_c<br>colipa_p + udcpdp_p --> colipap_p + udcpp_p<br>colipap_p + atp_c + h2o_c --> pi_c + h_c + colipap_e + adp_c<br>colipa_c + atp_c + h2o_c --> h_c + pi_c + adp_c + colipa_p<br>colipa_p + h2o_c + atp_c --> h_c + pi_c + colipa_e + adp_c<br>2.0 cpgn_c + fadh2_c --> 2.0 h_c + fad_c + 2.0 fe2_c + 2.0 cpgn_un_c<br>2.0 cpgn_c + fmnh2_c --> 2.0 h_c + 2.0 fe2_c + fmn_c + 2.0 cpgn_un_c<br>2.0 cpgn_c + rbflvrd_c --> 2.0 h_c + 2.0 fe2_c + 2.0 cpgn_un_c + ribflv_c<br>cpgn_un_p + h_p --> cpgn_un_e + h_c<br>h_p + cpgn_un_c --> cpgn_un_p + h_c<br>cpgn_p + atp_c + h2o_c --> h_c + pi_c + cpgn_c + adp_c<br>fe3_e + cpgn_un_e --> cpgn_e<br>cpgn_e + h_p --> cpgn_p + h_c<br>ahdt_c + h2o_c --> acald_c + h_c + pppi_c + cph4_c<br>gtp_c + h2o_c --> cpmp_c + ppi_c<br>cpppg3_c + 2.0 h_c + o2_c --> 2.0 co2_c + pppg9_c + 2.0 h2o_c<br>cpppg3_c + 2.0 amet_c --> 2.0 met__L_c + 2.0 co2_c + 2.0 dad_5_c + pppg9_c<br>bbtcoa_c + crn_c <=> gbbtn_c + crncoa_c<br>crn_c + coa_c + atp_c --> crncoa_c + pi_c + adp_c<br>crncoa_c <=> crnDcoa_c<br>ctbtcoa_c + crn_c <=> crncoa_c + ctbt_c<br>crncoa_c <=> ctbtcoa_c + h2o_c<br>crn__D_c + coa_c + atp_c --> crnDcoa_c + pi_c + adp_c<br>crn__D_p + atp_c + h2o_c --> crn__D_c + h_c + pi_c + adp_c<br>h_p + crn__D_p <=> h_c + crn__D_c<br>crn__D_e <=> crn__D_p<br>crn_p + atp_c + h2o_c --> crn_c + h_c + pi_c + adp_c<br>crn_p + h_p <=> crn_c + h_c<br>crn_p + gbbtn_c --> crn_c + gbbtn_p<br>crn_p + crn__D_c --> crn_c + crn__D_p<br>crn_e <=> crn_p<br>oaa_c + accoa_c + h2o_c --> cit_c + h_c + coa_c<br>h_c + csn_c + h2o_c --> nh4_c + ura_c<br>csn_p + h_p --> csn_c + h_c<br>csn_e <=> csn_p<br>ctbt_c + coa_c + atp_c --> ctbtcoa_c + pi_c + adp_c<br>ctbt_p + atp_c + h2o_c --> ctbt_c + h_c + pi_c + adp_c<br>ctbt_p + h_p <=> ctbt_c + h_c<br>tdecoa_c <=> td2coa_c<br>hdcoa_c <=> hdd2coa_c<br>odecoa_c <=> od2coa_c<br>utp_c + atp_c + h2o_c + gln__L_c --> ctp_c + pi_c + glu__L_c + 2.0 h_c + adp_c<br>4.0 h_p + 4.0 cu_p + o2_p --> 4.0 cu2_p + 2.0 h2o_p<br>cu_c + atp_c + h2o_c --> cu_p + h_c + pi_c + adp_c<br>cu2_c + atp_c + h2o_c --> cu2_p + h_c + pi_c + adp_c<br>cu2_e <=> cu2_p<br>cu2_p --> cu2_c<br>h_e + cu_c --> cu_e + h_c<br>cu_e <=> cu_p<br>tsul_c + cyan_c --> tcynt_c + so3_c + h_c<br>cyan_p + tsul_p --> so3_p + h_p + tcynt_p<br>cyan_e <=> cyan_p<br>cynt_c + hco3_c + 3.0 h_c --> nh4_c + 2.0 co2_c<br>cynt_p + h_p --> cynt_c + h_c<br>cynt_e <=> cynt_p<br>cys__D_c + h2o_c --> pyr_c + nh4_c + h2s_c<br>cys__L_c + h2o_c --> pyr_c + nh4_c + h2s_c<br>h2o_c + atp_c + cys__D_p --> h_c + pi_c + cys__D_c + adp_c<br>cys__D_e <=> cys__D_p<br>acser_c + h2s_c --> ac_c + cys__L_c + h_c<br>3sala_c + 2.0 h_c --> ala__L_c + so2_c<br>cyst__L_c + h2o_c --> pyr_c + nh4_c + hcys__L_c<br>cys__L_c + trnacys_c + atp_c --> cystrna_c + ppi_c + amp_c<br>cys__L_c + atp_c + h2o_c --> cys__L_p + h_c + pi_c + adp_c<br>cys__L_p + atp_c + h2o_c --> cys__L_c + h_c + pi_c + adp_c<br>cys__L_e <=> cys__L_p<br>cys__L_c --> cys__L_p<br>2.0 h_c + 0.5 o2_c + mql8_c --> mqn8_c + 2.0 h_p + h2o_c<br>2.0 h_c + 0.5 o2_c + q8h2_c --> 2.0 h_p + q8_c + h2o_c<br>4.0 h_c + 0.5 o2_c + q8h2_c --> 4.0 h_p + q8_c + h2o_c<br>h2o_c + cytd_c + h_c --> nh4_c + uri_c<br>cytd_c + h2o_c --> rib__D_c + csn_c<br>gtp_c + cytd_c --> cmp_c + gdp_c + h_c<br>h_p + cytd_p --> cytd_c + h_c<br>h_p + cytd_p <=> cytd_c + h_c<br>cytd_e <=> cytd_p<br>atp_c + cmp_c <=> cdp_c + adp_c<br>dcmp_c + atp_c <=> adp_c + dcdp_c<br>lac__D_p + h_p <=> lac__D_c + h_c<br>lac__D_e <=> lac__D_p<br>ala__D_c + fad_c + h2o_c --> pyr_c + nh4_c + fadh2_c<br>h_c + h2o_c + dad_2_c --> din_c + nh4_c<br>atp_c + damp_c <=> dadp_c + adp_c<br>dad_2_p + h_p --> dad_2_c + h_c<br>dad_2_e <=> dad_2_p<br>12dgr120_c + atp_c --> pa120_c + adp_c + h_c<br>12dgr140_c + atp_c --> pa140_c + adp_c + h_c<br>12dgr141_c + atp_c --> adp_c + pa141_c + h_c<br>12dgr160_c + atp_c --> adp_c + pa160_c + h_c<br>12dgr161_c + atp_c --> h_c + adp_c + pa161_c<br>12dgr180_c + atp_c --> pa180_c + adp_c + h_c<br>12dgr181_c + atp_c --> pa181_c + adp_c + h_c<br>h_p + ala__D_p --> ala__D_c + h_c<br>ala__D_e <=> ala__D_p<br>damp_e <=> damp_p<br>23dappa_c + h2o_c --> pyr_c + 2.0 nh4_c<br>h_c + 26dap__M_c --> lys__L_c + co2_c<br>26dap_LL_c <=> 26dap__M_c<br>h2o_c + 26dap__M_p + atp_c --> adp_c + h_c + pi_c + 26dap__M_c<br>15dap_e <=> 15dap_p<br>h_c + pa120_c + ctp_c --> ppi_c + cdpdddecg_c<br>pa140_c + h_c + ctp_c --> cdpdtdecg_c + ppi_c<br>h_c + ctp_c + pa141_c --> ppi_c + cdpdtdec7eg_c<br>h_c + pa160_c + ctp_c --> cdpdhdecg_c + ppi_c<br>h_c + pa161_c + ctp_c --> ppi_c + cdpdhdec9eg_c<br>pa180_c + h_c + ctp_c --> cdpdodecg_c + ppi_c<br>pa181_c + ctp_c + h_c --> cdpdodec11eg_c + ppi_c<br>h2o_c + datp_c + h_c --> nh4_c + ditp_c<br>ru5p__D_c --> h_c + db4p_c + for_c<br>dann_c + co2_c + atp_c --> dtbt_c + 3.0 h_c + pi_c + adp_c<br>chtbs6p_c + h2o_c --> acgam6p_c + acgam_c<br>dca_e <=> dca_p<br>dcmp_e <=> dcmp_p<br>dctp_c + h2o_c + h_c --> nh4_c + dutp_c<br>dcyt_c + h2o_c + h_c --> duri_c + nh4_c<br>h_p + dcyt_p --> dcyt_c + h_c<br>dcyt_e <=> dcyt_p<br>ddca_e --> ddca_p<br>2dh3dgal_c + atp_c --> adp_c + 2dh3dgal6p_c + h_c<br>h_p + 2ddglcn_p <=> 2ddglcn_c + h_c<br>2ddglcn_e <=> 2ddglcn_p<br>2ddglcn_c + atp_c --> 2ddg6p_c + adp_c + h_c<br>e4p_c + pep_c + h2o_c --> pi_c + 2dda7p_c<br>2dh3dgal6p_c <=> pyr_c + g3p_c<br>dgmp_c + atp_c <=> dgdp_c + adp_c<br>dgmp_e <=> dgmp_p<br>h_p + dgsn_p --> dgsn_c + h_c<br>dgsn_e <=> dgsn_p<br>23dhacoa_c + h2o_c <=> 3hadpcoa_c<br>23dhmb_c --> 3mob_c + h2o_c<br>23dhmp_c --> h2o_c + 3mop_c<br>pep_c + dha_c --> dhap_c + pyr_c<br>dha_e <=> dha_p<br>dha_p <=> dha_c<br>nad_c + 23ddhb_c <=> 23dhb_c + nadh_c + h_c<br>23dhb_c + atp_c + h_c --> 23dhba_c + ppi_c<br>23dhbzs_c + h2o_c --> 23dhb_c + ser__L_c<br>nad_c + cenchddd_c --> dhcinnm_c + h_c + nadh_c<br>dhcinnm_c + o2_c --> hkntd_c + h_c<br>23dhdp_c + nadph_c + h_c --> nadp_c + thdp_c<br>pyr_c + aspsa_c --> 23dhdp_c + 2.0 h2o_c + h_c<br>dhf_c + nadph_c + h_c <=> nadp_c + thf_c<br>dhpt_c + glu__L_c + atp_c --> h_c + dhf_c + pi_c + adp_c<br>dhmpt_c + nadph_c + h_c --> nadp_c + thmnp_c<br>octdp_c + dhna_c + h_c --> ppi_c + co2_c + 2dmmql8_c<br>sbzcoa_c + h_c --> h2o_c + 14dhncoa_c<br>14dhncoa_c + h2o_c --> coa_c + dhna_c + h_c<br>dhnpt_c <=> gcald_c + 6hmhpt_c<br>dhnpt_c <=> dhmpt_c<br>q8_c + dhor__S_c --> orot_c + q8h2_c<br>mqn8_c + dhor__S_c --> orot_c + mql8_c<br>fum_c + dhor__S_c --> orot_c + succ_c<br>h2o_c + dhor__S_c <=> cbasp_c + h_c<br>nad_c + cechddd_c --> nadh_c + h_c + dhpppn_c<br>25drapp_c + h2o_c + h_c --> 5apru_c + nh4_c<br>6hmhptpp_c + 4abz_c --> dhpt_c + ppi_c<br>dhptd_c --> mdhdhf_c<br>dhptdn_c + nadph_c + 3.0 h_c --> thptdn_c + nadp_c<br>nadh_c + dhptdn_c + 3.0 h_c --> thptdn_c + nad_c<br>ahdt_c <=> dhmptp_c<br>2dda7p_c --> pi_c + 3dhq_c<br>3dhq_c --> 3dhsk_c + h2o_c<br>dimp_e <=> dimp_p<br>din_p + h_p --> h_c + din_c<br>din_e <=> din_p<br>25dkglcn_c + nadph_c + h_c --> 2dhguln_c + nadp_c<br>h_c + 25dkglcn_c + nadh_c --> nad_c + 5dglcn_c<br>nadph_c + 25dkglcn_c + h_c --> nadp_c + 5dglcn_c<br>dmpp_c + ipdp_c --> grdp_c + ppi_c<br>h2mb4p_c + h_c + nadh_c --> dmpp_c + nad_c + h2o_c<br>2omhmbl_c + amet_c --> h_c + ahcys_c + q8h2_c<br>dmso_c + mql8_c --> mqn8_c + dms_c + h2o_c<br>mql8_c + dmso_p --> dms_p + mqn8_c + h2o_p<br>2dmmql8_c + dmso_c --> dms_c + 2dmmq8_c + h2o_c<br>2dmmql8_c + dmso_p --> dms_p + h2o_p + 2dmmq8_c<br>dmso_e <=> dmso_p<br>dmso_p <=> dmso_c<br>dms_e <=> dms_p<br>dhpmp_c + h2o_c --> pi_c + dhnpt_c<br>ahdt_c + h2o_c --> ppi_c + dhpmp_c + h_c<br>23doguln_c + nadh_c + h_c --> nad_c + 3dhguln_c<br>dopa_e <=> dopa_p<br>doxrbcn_e <=> doxrbcn_p<br>doxrbcn_p + h_p --> doxrbcn_e + h_c<br>dpcoa_c + atp_c --> coa_c + adp_c + h_c<br>2dhp_c + nadph_c + h_c --> nadp_c + pant__R_c<br>2dr5p_c --> acald_c + g3p_c<br>q8_c + dsbard_p --> dsbaox_p + q8h2_c<br>dsbard_p + mqn8_c --> dsbaox_p + mql8_c<br>2.0 gthrd_p + dsbcox_p --> dsbcrd_p + gthox_p<br>trdrd_c + dsbdox_c --> dsbdrd_c + trdox_c<br>dsbgox_p + 2.0 gthrd_p --> dsbgrd_p + gthox_p<br>nadp_c + ser__D_c <=> 2amsa_c + nadph_c + h_c<br>ser__D_p + h_p --> ser__D_c + h_c<br>ser__D_e <=> ser__D_p<br>tartr__D_c --> oaa_c + h2o_c<br>dtmp_c + atp_c <=> dtdp_c + adp_c<br>dtmp_e <=> dtmp_p<br>dump_e <=> dump_p<br>56dura_c + nad_c <=> nadh_c + ura_c + h_c<br>duri_c + atp_c --> dump_c + adp_c + h_c<br>pi_c + duri_c <=> 2dr1p_c + ura_c<br>duri_p + h_p --> duri_c + h_c<br>duri_e <=> duri_p<br>h2o_c + dutp_c --> dump_c + ppi_c + h_c<br>h_c + nadph_c + dxyl5p_c --> nadp_c + 2me4p_c<br>g3p_c + pyr_c + h_c --> co2_c + dxyl5p_c<br>dxyl_c + atp_c --> h_c + adp_c + dxyl5p_c<br>e4p_c + nad_c + h2o_c <=> nadh_c + 2.0 h_c + 4per_c<br>nadh_c + h_c + tdec2eACP_c --> dcaACP_c + nad_c<br>h_c + nadph_c + tdec2eACP_c --> dcaACP_c + nadp_c<br>nadh_c + tddec2eACP_c + h_c --> nad_c + ddcaACP_c<br>tddec2eACP_c + nadph_c + h_c --> nadp_c + ddcaACP_c<br>nadh_c + t3c5ddeceACP_c + h_c --> nad_c + cddec5eACP_c<br>t3c5ddeceACP_c + nadph_c + h_c --> nadp_c + cddec5eACP_c<br>tmrs2eACP_c + nadh_c + h_c --> nad_c + myrsACP_c<br>tmrs2eACP_c + nadph_c + h_c --> myrsACP_c + nadp_c<br>nadh_c + t3c7mrseACP_c + h_c --> tdeACP_c + nad_c<br>nadph_c + t3c7mrseACP_c + h_c --> nadp_c + tdeACP_c<br>nadh_c + tpalm2eACP_c + h_c --> palmACP_c + nad_c<br>nadph_c + tpalm2eACP_c + h_c --> palmACP_c + nadp_c<br>nadh_c + t3c9palmeACP_c + h_c --> nad_c + hdeACP_c<br>t3c9palmeACP_c + nadph_c + h_c --> nadp_c + hdeACP_c<br>toctd2eACP_c + nadh_c + h_c --> nad_c + ocdcaACP_c<br>toctd2eACP_c + nadph_c + h_c --> nadp_c + ocdcaACP_c<br>t3c11vaceACP_c + h_c + nadh_c --> nad_c + octeACP_c<br>t3c11vaceACP_c + nadph_c + h_c --> nadp_c + octeACP_c<br>nadh_c + but2eACP_c + h_c --> butACP_c + nad_c<br>but2eACP_c + nadph_c + h_c --> butACP_c + nadp_c<br>nadh_c + thex2eACP_c + h_c --> hexACP_c + nad_c<br>thex2eACP_c + nadph_c + h_c --> hexACP_c + nadp_c<br>toct2eACP_c + h_c + nadh_c --> nad_c + ocACP_c<br>toct2eACP_c + nadph_c + h_c --> nadp_c + ocACP_c<br>eca4colipa_p + atp_c + h2o_c --> eca4colipa_e + h_c + pi_c + adp_c<br>eca4und_p + colipa_p --> udcpdp_p + h_p + eca4colipa_p<br>2.0 unagamuf_p --> h_p + eca2und_p + udcpdp_p<br>eca2und_p + unagamuf_p --> eca3und_p + h_p + udcpdp_p<br>eca3und_p + unagamuf_p --> eca4und_p + h_p + udcpdp_p<br>unagamuf_c --> unagamuf_p<br>3hbcoa_c <=> b2coa_c + h2o_c<br>3hhcoa_c <=> hx2coa_c + h2o_c<br>3hocoa_c <=> h2o_c + oc2coa_c<br>3hdcoa_c <=> dc2coa_c + h2o_c<br>3hddcoa_c <=> h2o_c + dd2coa_c<br>3htdcoa_c <=> td2coa_c + h2o_c<br>3hhdcoa_c <=> hdd2coa_c + h2o_c<br>3hodcoa_c <=> od2coa_c + h2o_c<br>2ddg6p_c --> g3p_c + pyr_c<br>6pgc_c --> 2ddg6p_c + h2o_c<br>kdo2lipid4_c + ddcaACP_c --> ACP_c + kdo2lipid4L_c<br>myrsACP_c + kdo2lipid4L_c --> ACP_c + lipa_c<br>hdeACP_c + kdo2lipid4_c --> kdo2lipid4p_c + ACP_c<br>kdo2lipid4p_c + myrsACP_c --> ACP_c + lipa_cold_c<br>egmeACP_c + nadph_c + h_c --> nadp_c + gmeACP_c<br>enlipa_p + atp_c + h2o_c --> enlipa_e + h_c + pi_c + adp_c<br>2pg_c <=> pep_c + h2o_c<br>3.0 seramp_c + 3.0 23dhba_c --> enter_c + 6.0 amp_c + 9.0 h_c<br>enter_c + 3.0 h2o_c --> 3.0 23dhbzs_c + 3.0 h_c<br>3.0 h2o_c + feenter_c --> fe3_c + 3.0 23dhbzs_c + 3.0 h_c<br>epmeACP_c + nadph_c + h_c --> nadp_c + pmeACP_c<br>etha_c --> acald_c + nh4_c<br>h_p + etha_p --> etha_c + h_c<br>etha_e <=> etha_p<br>ethso3_p + atp_c + h2o_c --> h_c + pi_c + adp_c + ethso3_c<br>ethso3_e <=> ethso3_p<br>etoh_e <=> etoh_p<br>etoh_p <=> etoh_c<br>f6p_c <=> g3p_c + dha_c<br>h2o_c + f6p_c --> pi_c + fru_c<br>f6p_p + 2.0 pi_c --> 2.0 pi_p + f6p_c<br>f6p_e <=> f6p_p<br>dcaACP_c + h2o_c --> ACP_c + dca_c + h_c<br>h2o_c + ddcaACP_c --> ACP_c + ddca_c + h_c<br>myrsACP_c + h2o_c --> ACP_c + ttdca_c + h_c<br>tdeACP_c + h2o_c --> ttdcea_c + ACP_c + h_c<br>palmACP_c + h2o_c --> ACP_c + hdca_c + h_c<br>hdeACP_c + h2o_c --> ACP_c + hdcea_c + h_c<br>ocACP_c + h2o_c --> ACP_c + octa_c + h_c<br>dcacoa_c + h2o_c --> dca_c + coa_c + h_c<br>ddcacoa_c + h2o_c --> coa_c + ddca_c + h_c<br>tdcoa_c + h2o_c --> coa_c + ttdca_c + h_c<br>tdecoa_c + h2o_c --> ttdcea_c + coa_c + h_c<br>h2o_c + pmtcoa_c --> coa_c + hdca_c + h_c<br>hdcoa_c + h2o_c --> coa_c + hdcea_c + h_c<br>stcoa_c + h2o_c --> ocdca_c + coa_c + h_c<br>h2o_c + odecoa_c --> coa_c + ocdcea_c + h_c<br>hxcoa_c + h2o_c --> coa_c + hxa_c + h_c<br>occoa_c + h2o_c --> coa_c + octa_c + h_c<br>dca_p + coa_c + h_p + atp_c --> ppi_c + h_c + dcacoa_c + amp_c<br>ddca_p + coa_c + h_p + atp_c --> ddcacoa_c + ppi_c + h_c + amp_c<br>coa_c + ttdca_p + h_p + atp_c --> tdcoa_c + ppi_c + h_c + amp_c<br>ttdcea_p + coa_c + h_p + atp_c --> ppi_c + h_c + tdecoa_c + amp_c<br>coa_c + h_p + hdca_p + atp_c --> pmtcoa_c + ppi_c + h_c + amp_c<br>atp_c + coa_c + h_p + hdcea_p --> ppi_c + h_c + hdcoa_c + amp_c<br>coa_c + h_p + ocdca_p + atp_c --> stcoa_c + ppi_c + h_c + amp_c<br>ocdcea_p + coa_c + h_p + atp_c --> ppi_c + h_c + odecoa_c + amp_c<br>hxa_p + coa_c + h_p + atp_c --> ppi_c + hxcoa_c + h_c + amp_c<br>octa_p + h_p + coa_c + atp_c --> ppi_c + h_c + occoa_c + amp_c<br>fad_c + nadh_c + h_c --> nad_c + fadh2_c<br>nadph_c + fad_c + h_c --> nadp_c + fadh2_c<br>nad_c + hmgth_c <=> Sfglutth_c + nadh_c + h_c<br>fald_e <=> fald_p<br>fald_p <=> fald_c<br>gthrd_c + fald_c <=> hmgth_c<br>fdp_c <=> dhap_c + g3p_c<br>s17bp_c <=> e4p_c + dhap_c<br>fdp_c + h2o_c --> pi_c + f6p_c<br>fuc__L_c <=> fcl__L_c<br>fcl__L_c + atp_c --> fc1p_c + adp_c + h_c<br>fc1p_c <=> dhap_c + lald__L_c<br>fe2_c + ppp9_c --> pheme_c + 2.0 h_c<br>2.0 h_c + for_p + q8_c --> h_p + co2_p + q8h2_c<br>mqn8_c + 2.0 h_c + for_p --> mql8_c + co2_p + h_p<br>fmnh2_c + o2_c + isetac_c --> so3_c + gcald_c + h_c + fmn_c + h2o_c<br>mso3_c + fmnh2_c + o2_c --> fald_c + so3_c + h_c + fmn_c + h2o_c<br>fmnh2_c + ethso3_c + o2_c --> acald_c + so3_c + h_c + fmn_c + h2o_c<br>fmnh2_c + butso3_c + o2_c --> so3_c + h_c + fmn_c + btal_c + h2o_c<br>sulfac_c + fmnh2_c + o2_c --> so3_c + h_c + fmn_c + glx_c + h2o_c<br>fe2_p + atp_c + h2o_c --> h_c + pi_c + fe2_c + adp_c<br>fe2_p + h_p --> fe2_c + h_c<br>fe2_c + h_p --> fe2_p + h_c<br>fe2_e <=> fe2_p<br>fe2_p --> fe2_c<br>fe3dcit_p + h2o_c + atp_c --> fe3_c + 2.0 cit_c + h_c + pi_c + adp_c<br>h_p + fe3dcit_e --> fe3dcit_p + h_c<br>fe3dhbzs_c --> fe3_c + 23dhbzs_c<br>fe3dhbzs_p + atp_c + h2o_c --> fe3dhbzs_c + h_c + pi_c + adp_c<br>fe3dhbzs_e + h_p --> fe3dhbzs_p + h_c<br>2.0 fe3hox_c + fadh2_c --> fad_c + 2.0 h_c + 2.0 fe3hox_un_c + 2.0 fe2_c<br>2.0 fe3hox_c + fmnh2_c --> 2.0 h_c + 2.0 fe3hox_un_c + 2.0 fe2_c + fmn_c<br>2.0 fe3hox_c + rbflvrd_c --> 2.0 h_c + 2.0 fe3hox_un_c + 2.0 fe2_c + ribflv_c<br>fe3hox_un_c + h_p --> fe3hox_un_p + h_c<br>h_p + fe3hox_un_p --> fe3hox_un_e + h_c<br>fe3hox_p + atp_c + h2o_c --> fe3hox_c + pi_c + h_c + adp_c<br>fe3hox_un_e + fe3_e --> fe3hox_e<br>fe3hox_e + h_p --> fe3hox_p + h_c<br>2.0 fe3_c + fadh2_c --> fad_c + 2.0 fe2_c + 2.0 h_c<br>fe3_p + atp_c + h2o_c --> fe3_c + h_c + pi_c + adp_c<br>fe3_e <=> fe3_p<br>2.0 fecrm_c + fadh2_c --> 2.0 h_c + fad_c + 2.0 fe2_c + 2.0 fecrm_un_c<br>2.0 fecrm_c + fmnh2_c --> 2.0 h_c + 2.0 fe2_c + fmn_c + 2.0 fecrm_un_c<br>2.0 fecrm_c + rbflvrd_c --> 2.0 h_c + 2.0 fe2_c + 2.0 fecrm_un_c + ribflv_c<br>h_p + fecrm_un_p --> fecrm_un_e + h_c<br>h_p + fecrm_un_c --> fecrm_un_p + h_c<br>fecrm_p + atp_c + h2o_c --> fecrm_c + h_c + pi_c + adp_c<br>fe3_e + fecrm_un_e --> fecrm_e<br>fecrm_e + h_p --> fecrm_p + h_c<br>2.0 feenter_c + fadh2_c --> 2.0 h_c + 2.0 enter_c + fad_c + 2.0 fe2_c<br>fmnh2_c + 2.0 feenter_c --> 2.0 h_c + 2.0 enter_c + 2.0 fe2_c + fmn_c<br>rbflvrd_c + 2.0 feenter_c --> 2.0 h_c + 2.0 enter_c + 2.0 fe2_c + ribflv_c<br>feenter_p + atp_c + h2o_c --> h_c + pi_c + adp_c + feenter_c<br>fe3_e + enter_e --> feenter_e<br>h_p + enter_p --> h_c + enter_e<br>feenter_e + h_p --> feenter_p + h_c<br>enter_c + h_p --> h_c + enter_p<br>2.0 feoxam_c + fadh2_c --> 2.0 feoxam_un_c + 2.0 h_c + fad_c + 2.0 fe2_c<br>2.0 feoxam_c + fmnh2_c --> 2.0 feoxam_un_c + 2.0 h_c + 2.0 fe2_c + fmn_c<br>2.0 feoxam_c + rbflvrd_c --> 2.0 feoxam_un_c + 2.0 h_c + 2.0 fe2_c + ribflv_c<br>feoxam_un_p + h_p --> feoxam_un_e + h_c<br>feoxam_un_c + h_p --> feoxam_un_p + h_c<br>feoxam_p + atp_c + h2o_c --> feoxam_c + h_c + pi_c + adp_c<br>feoxam_un_e + fe3_e --> feoxam_e<br>h_p + feoxam_e --> feoxam_p + h_c<br>4.0 fe2_p + o2_p + 4.0 h_p --> 2.0 h2o_p + 4.0 fe3_p<br>2.0 4fe4s_c + 2.0 h_c + h2o2_c --> 2.0 fe3_c + 2.0 3fe4s_c + 2.0 h2o_c<br>2.0 4fe4s_c + 2.0 no_c + 2.0 h_c --> 2.0 fe3_c + n2o_c + 2.0 3fe4s_c + h2o_c<br>fe2_c + 3fe4s_c --> 4fe4s_c<br>suc6p_c + h2o_c --> g6p_c + fru_c<br>h_c + for_c --> co2_c + h2_c<br>nadph_c + 2.0 flxso_c --> nadp_c + 2.0 flxr_c + h_c<br>ribflv_c + nadph_c + h_c --> rbflvrd_c + nadp_c<br>ribflv_c + nadh_c + h_c --> rbflvrd_c + nad_c<br>mettrna_c + 10fthf_c --> fmettrna_c + thf_c + h_c<br>fmn_c + atp_c + h_c --> fad_c + ppi_c<br>nadh_c + fmn_c + h_c --> fmnh2_c + nad_c<br>fmn_c + nadph_c + h_c --> fmnh2_c + nadp_c<br>5fthf_c + h_c --> methf_c + h2o_c<br>oxa_c + forcoa_c <=> oxalcoa_c + for_c<br>for_p + h_p --> h_c + for_c<br>for_e <=> for_p<br>for_c --> for_p<br>mql8_c + fum_c --> mqn8_c + succ_c<br>2dmmql8_c + fum_c --> 2dmmq8_c + succ_c<br>f1p_c + atp_c --> fdp_c + adp_c + h_c<br>h2o_c + frulysp_c <=> g6p_c + lys__L_c<br>psclys_c <=> frulys_c<br>frulys_c + atp_c --> h_c + adp_c + frulysp_c<br>frulys_p + h_p --> frulys_c + h_c<br>frulys_e <=> frulys_p<br>fruur_p + h_p <=> fruur_c + h_c<br>fruur_e <=> fruur_p<br>fru_p + pep_c --> pyr_c + f6p_c<br>fru_p + pep_c --> pyr_c + f1p_c<br>fru_e <=> fru_p<br>h2o_c + 10fthf_c --> h_c + thf_c + for_c<br>thf_c + for_c + atp_c --> pi_c + 10fthf_c + adp_c<br>fuc__L_e <=> fuc__L_p<br>fuc__L_p + h_p <=> fuc__L_c + h_c<br>fum_c + h2o_c <=> mal__L_c<br>2.0 h_p + fum_p --> fum_c + 2.0 h_c<br>3.0 h_p + fum_p --> fum_c + 3.0 h_c<br>fum_e <=> fum_p<br>fusa_e <=> fusa_p<br>h_p + fusa_p --> fusa_e + h_c<br>accoa_c + gam1p_c --> acgam1p_c + coa_c + h_c<br>g1p_p + h2o_p --> glc__D_p + pi_p<br>g1p_c + dttp_c + h_c --> dtdpglu_c + ppi_c<br>g1p_e <=> g1p_p<br>glu1sa_c <=> 5aop_c<br>glyc2p_c + h2o_c --> glyc_c + pi_c<br>glyc2p_p + h2o_p --> pi_p + glyc_p<br>ddcaACP_c + glyc3p_c --> ACP_c + 1ddecg3p_c<br>glyc3p_c + myrsACP_c --> 1tdecg3p_c + ACP_c<br>glyc3p_c + tdeACP_c --> 1tdec7eg3p_c + ACP_c<br>palmACP_c + glyc3p_c --> 1hdecg3p_c + ACP_c<br>glyc3p_c + hdeACP_c --> 1hdec9eg3p_c + ACP_c<br>glyc3p_c + ocdcaACP_c --> 1odecg3p_c + ACP_c<br>glyc3p_c + octeACP_c --> 1odec11eg3p_c + ACP_c<br>g3pc_p + atp_c + h2o_c --> h_c + pi_c + adp_c + g3pc_c<br>g3pc_e <=> g3pc_p<br>glyc3p_c + nadp_c <=> dhap_c + nadph_c + h_c<br>glyc3p_c + q8_c --> dhap_c + q8h2_c<br>glyc3p_c + mqn8_c --> dhap_c + mql8_c<br>glyc3p_c + 2dmmq8_c --> 2dmmql8_c + dhap_c<br>g3pe_p + atp_c + h2o_c --> h_c + g3pe_c + pi_c + adp_c<br>g3pe_e <=> g3pe_p<br>g3pg_p + atp_c + h2o_c --> g3pg_c + h_c + pi_c + adp_c<br>g3pg_e <=> g3pg_p<br>g3pi_p + atp_c + h2o_c --> g3pi_c + h_c + pi_c + adp_c<br>g3pi_e <=> g3pi_p<br>g3ps_p + atp_c + h2o_c --> g3ps_c + pi_c + h_c + adp_c<br>g3ps_e <=> g3ps_p<br>glyc3p_c + h2o_c --> glyc_c + pi_c<br>glu5sa_c --> h_c + 1pyr5c_c + h2o_c<br>nadph_c + h_c + glu5p_c --> pi_c + glu5sa_c + nadp_c<br>gam6p_c + h2o_c --> nh4_c + f6p_c<br>nadp_c + g6p_c <=> nadph_c + 6pgl_c + h_c<br>g6p_c + h2o_c --> glc__D_c + pi_c<br>2.0 pi_c + g6p_p --> g6p_c + 2.0 pi_p<br>g6p_e <=> g6p_p<br>h2o_p + gal1p_p --> gal_p + pi_p<br>gal1p_e <=> gal1p_p<br>gal_bD_e <=> gal_bD_p<br>galct__D_c --> 5dh4dglc_c + h2o_c<br>nad_c + galctn__L_c --> nadh_c + h_c + tagur_c<br>galctn__D_c --> 2dh3dgal_c + h2o_c<br>h_p + galctn__L_p --> galctn__L_c + h_c<br>galctn__L_e <=> galctn__L_p<br>galctn__D_p + h_p --> galctn__D_c + h_c<br>galctn__D_e <=> galctn__D_p<br>h_p + galct__D_p <=> galct__D_c + h_c<br>galct__D_e <=> galct__D_p<br>gal_c + atp_c <=> adp_c + gal1p_c + h_c<br>gal_bD_p --> gal_p<br>h2o_c + melib_c --> glc__D_c + gal_c<br>gicolipa_c + udpg_c --> h_c + udp_c + gagicolipa_c<br>pep_c + galt_p --> galt1p_c + pyr_c<br>galt_e <=> galt_p<br>galur_p + h_p <=> h_c + galur_c<br>galur_e <=> galur_p<br>utp_c + g1p_c + h_c --> udpg_c + ppi_c<br>gal_p + atp_c + h2o_c --> h_c + pi_c + adp_c + gal_c<br>gal_p + h_p --> gal_c + h_c<br>gal_e <=> gal_p<br>gam6p_p + 2.0 pi_c --> 2.0 pi_p + gam6p_c<br>gam6p_e <=> gam6p_p<br>pep_c + gam_p --> gam6p_c + pyr_c<br>gam_e <=> gam_p<br>g3p_c + pi_c + nad_c <=> h_c + 13dpg_c + nadh_c<br>10fthf_c + gar_c <=> fgam_c + thf_c + h_c<br>gar_c + for_c + atp_c --> fgam_c + h_c + pi_c + adp_c<br>gbbtn_e <=> gbbtn_p<br>gcald_c + nad_c + h2o_c --> 2.0 h_c + nadh_c + glyclt_c<br>gdpddman_c --> gdpofuc_c<br>gdp_c + atp_c --> amp_c + h_c + ppgpp_c<br>h2o_c + gdpmann_c --> gdp_c + man_c + h_c<br>h2o_c + gdpmann_c --> gmp_c + man1p_c + 2.0 h_c<br>gdptp_c + h2o_c --> ppi_c + gtp_c<br>gdp_e <=> gdp_p<br>f6p_c + gln__L_c --> gam6p_c + glu__L_c<br>nadp_c + ggbutal_c + h2o_c <=> nadph_c + 2.0 h_c + gg4abut_c<br>gg4abut_c + h2o_c --> glu__L_c + 4abut_c<br>ggptrc_c + o2_c + h2o_c --> h2o2_c + nh4_c + ggbutal_c<br>ptrc_c + glu__L_c + atp_c --> ggptrc_c + h_c + pi_c + adp_c<br>nadh_c + h_c + sucsal_c <=> nad_c + ghb_c<br>thf_c + ser__L_c <=> gly_c + h2o_c + mlthf_c<br>gmp_c + atp_c <=> gdp_c + adp_c<br>glycogen_c --> bglycogen_c<br>glc__D_c + accoa_c <=> acglc__D_c + coa_c<br>glc__D_p + h2o_p + q8_c --> glcn_p + h_p + q8h2_c<br>glcn_p + h_p <=> glcn_c + h_c<br>glcn_e <=> glcn_p<br>glycogen_c + pi_c --> g1p_c<br>bglycogen_c + pi_c --> g1p_c<br>5dh4dglc_c --> pyr_c + 2h3oppan_c<br>glcr_c --> 5dh4dglc_c + h2o_c<br>h_p + glcr_p <=> glcr_c + h_c<br>glcr_e <=> glcr_p<br>adpglc_c --> glycogen_c + adp_c + h_c<br>icolipa_c + udpg_c --> gicolipa_c + udp_c + h_c<br>udpg_c + gagicolipa_c --> udp_c + h_c + ggagicolipa_c<br>udpg_c + ggagicolipa_c --> gggagicolipa_c + udp_c + h_c<br>glcur1p_e <=> glcur1p_p<br>glcur_p + h_p <=> glcur_c + h_c<br>glcur_e <=> glcur_p<br>glc__D_p + atp_c + h2o_c --> h_c + glc__D_c + pi_c + adp_c<br>pep_c + glc__D_p --> g6p_c + pyr_c<br>glc__D_p + h_p --> glc__D_c + h_c<br>glc__D_e <=> glc__D_p<br>glc__D_e --> glc__D_p<br>bglycogen_c --> glycogen_c<br>g1p_c + atp_c + h_c --> adpglc_c + ppi_c<br>nh4_c + glu__L_c + atp_c --> gln__L_c + h_c + pi_c + adp_c<br>trnagln_c + gln__L_c + atp_c --> glntrna_c + ppi_c + amp_c<br>gln__L_p + atp_c + h2o_c --> gln__L_c + h_c + pi_c + adp_c<br>gln__L_e <=> gln__L_p<br>galt1p_c + nad_c <=> tag6p__D_c + nadh_c + h_c<br>glu__L_c + atp_c --> glu5p_c + adp_c<br>glu__L_p + 4abut_c <=> 4abut_p + glu__L_c<br>cys__L_c + glu__L_c + atp_c --> h_c + pi_c + glucys_c + adp_c<br>glu__L_c + h_c --> 4abut_c + co2_c<br>nadp_c + glu__L_c + h2o_c <=> h_c + akg_c + nh4_c + nadph_c<br>gln__L_c + h2o_c --> nh4_c + glu__L_c<br>h2o_p + gln__L_p --> nh4_p + glu__L_p<br>prpp_c + gln__L_c + h2o_c --> pram_c + ppi_c + glu__L_c<br>glu__D_c <=> glu__L_c<br>h_c + akg_c + nadph_c + gln__L_c --> nadp_c + 2.0 glu__L_c<br>nadph_c + glutrna_c + h_c --> nadp_c + glu1sa_c + trnaglu_c<br>trnaglu_c + glu__L_c + atp_c --> glutrna_c + ppi_c + amp_c<br>glu__L_p + atp_c + h2o_c --> h_c + pi_c + adp_c + glu__L_c<br>h_p + glu__L_p <=> glu__L_c + h_c<br>glu__L_p + na1_p --> glu__L_c + na1_c<br>glu__L_e <=> glu__L_p<br>2.0 glx_c + h_c --> 2h3oppan_c + co2_c<br>glyald_e <=> glyald_p<br>glyald_p <=> glyald_c<br>accoa_c + gly_c <=> 2aobut_c + coa_c<br>glyb_p + atp_c + h2o_c --> glyb_c + h_c + pi_c + adp_c<br>glyb_p + h_p --> glyb_c + h_c<br>glyb_e <=> glyb_p<br>glyc2p_p + atp_c + h2o_c --> h_c + glyc2p_c + pi_c + adp_c<br>glyc2p_e <=> glyc2p_p<br>glyc3p_p + atp_c + h2o_c --> glyc3p_c + h_c + pi_c + adp_c<br>pi_c + glyc3p_p --> glyc3p_c + pi_p<br>glyc3p_e <=> glyc3p_p<br>h_p + glyc__R_p <=> h_c + glyc__R_c<br>glyc__R_e <=> glyc__R_p<br>glyc_c + nad_c --> h_c + nadh_c + dha_c<br>atp_c + glyc__R_c --> adp_c + 3pg_c + h_c<br>glyc__R_c + atp_c --> adp_c + 2pg_c + h_c<br>thf_c + gly_c + nad_c --> nadh_c + co2_c + nh4_c + mlthf_c<br>nadh_c + glx_c + h_c --> nad_c + glyclt_c<br>nadph_c + glx_c + h_c --> glyclt_c + nadp_c<br>glyclt_p + h_p <=> glyclt_c + h_c<br>glyclt_p + na1_p --> glyclt_c + na1_c<br>glyclt_e <=> glyclt_p<br>q8_c + glyclt_c --> glx_c + q8h2_c<br>mqn8_c + glyclt_c --> mql8_c + glx_c<br>2dmmq8_c + glyclt_c --> 2dmmql8_c + glx_c<br>glyc_e <=> glyc_p<br>glyc_c <=> glyc_p<br>glyc_c + atp_c --> glyc3p_c + adp_c + h_c<br>lgt__S_c + h2o_c --> gthrd_c + h_c + lac__D_c<br>mthgxl_c + h2o_c --> lac__D_c + h_c<br>trnagly_c + gly_c + atp_c --> ppi_c + amp_c + glytrna_c<br>gly_p + h_p --> gly_c + h_c<br>gly_p + h_p <=> gly_c + h_c<br>gly_p + na1_p --> gly_c + na1_c<br>gly_e <=> gly_p<br>gdpmann_c --> gdpddman_c + h2o_c<br>gmhep1p_c + atp_c + h_c --> adphep_DD_c + ppi_c<br>gmhep7p_c + atp_c --> adp_c + gmhep17bp_c + h_c<br>gmhep17bp_c + h2o_c --> gmhep1p_c + pi_c<br>nadph_c + 2.0 h_c + gmp_c --> imp_c + nh4_c + nadp_c<br>atp_c + h2o_c + xmp_c + gln__L_c --> 2.0 h_c + ppi_c + gmp_c + glu__L_c + amp_c<br>gmp_e <=> gmp_p<br>nadp_c + 6pgc_c --> ru5p__D_c + nadph_c + co2_c<br>glcn_c + atp_c --> 6pgc_c + adp_c + h_c<br>gdpofuc_c + nadph_c + h_c --> nadp_c + gdpfuc_c<br>gp4g_c + h2o_c --> 2.0 gdp_c + 2.0 h_c<br>g3pc_c + h2o_c --> glyc3p_c + chol_c + h_c<br>g3pc_p + h2o_p --> h_p + glyc3p_p + chol_p<br>g3pe_c + h2o_c --> glyc3p_c + etha_c + h_c<br>g3pe_p + h2o_p --> etha_p + h_p + glyc3p_p<br>h2o_c + g3ps_c --> glyc3p_c + h_c + ser__L_c<br>h2o_p + g3ps_p --> h_p + glyc3p_p + ser__L_p<br>h2o_c + g3pg_c --> glyc3p_c + glyc_c + h_c<br>g3pg_p + h2o_p --> h_p + glyc3p_p + glyc_p<br>h2o_c + g3pi_c --> inost_c + glyc3p_c + h_c<br>g3pi_p + h2o_p --> glyc3p_p + h_p + inost_p<br>grdp_c + ipdp_c --> frdp_c + ppi_c<br>2.0 gthrd_c + grxox_c --> grxrd_c + gthox_c<br>gsn_c + atp_c --> gmp_c + adp_c + h_c<br>h_p + gsn_p --> gsn_c + h_c<br>gsn_e <=> gsn_p<br>gtspmd_c + h2o_c --> gthrd_c + spmd_c<br>gthrd_c + spmd_c + atp_c --> h_c + gtspmd_c + pi_c + adp_c<br>gthox_e <=> gthox_p<br>gthox_c + nadph_c + h_c <=> 2.0 gthrd_c + nadp_c<br>h2o2_c + 2.0 gthrd_c --> gthox_c + 2.0 h2o_c<br>gthrd_p + h2o_p --> cgly_p + glu__L_p<br>gthrd_c + atp_c + h2o_c --> gthrd_p + h_c + pi_c + adp_c<br>gthrd_p + atp_c + h2o_c --> gthrd_c + h_c + pi_c + adp_c<br>gthrd_e <=> gthrd_p<br>gly_c + glucys_c + atp_c --> gthrd_c + h_c + pi_c + adp_c<br>h2o_c + gtp_c --> ahdt_c + h_c + for_c<br>gtp_c + 3.0 h2o_c --> for_c + 2.0 h_c + ppi_c + 25drapp_c<br>gdptp_c + h2o_c --> pi_c + h_c + ppgpp_c<br>atp_c + gtp_c --> gdptp_c + amp_c + h_c<br>h2o_c + gtp_c + h_c --> nh4_c + xtp_c<br>gtp_e <=> gtp_p<br>gtp_c --> 35cgmp_c + ppi_c<br>h_c + h2o_c + gua_c --> xan_c + nh4_c<br>prpp_c + gua_c --> gmp_c + ppi_c<br>gua_p + h_p --> h_c + gua_c<br>gua_e <=> gua_p<br>gua_p <=> gua_c<br>glcur_c <=> fruur_c<br>galur_c <=> tagur_c<br>glcur1p_p + h2o_p --> glcur_p + pi_p<br>h2o2_e <=> h2o2_p<br>h2o_e <=> h2o_p<br>h2o_p <=> h2o_c<br>h2s_c + 2.0 o2_c --> so4_c + 2.0 h_c<br>h2s_c --> h2s_p<br>h2s_e <=> h2s_p<br>h2_e <=> h2_p<br>h2_p <=> h2_c<br>nadh_c + h_c + aacoa_c <=> nad_c + 3hbcoa_c<br>3ohcoa_c + h_c + nadh_c <=> nad_c + 3hhcoa_c<br>3oocoa_c + h_c + nadh_c <=> 3hocoa_c + nad_c<br>3odcoa_c + nadh_c + h_c <=> nad_c + 3hdcoa_c<br>nadh_c + 3oddcoa_c + h_c <=> nad_c + 3hddcoa_c<br>nadh_c + 3otdcoa_c + h_c <=> nad_c + 3htdcoa_c<br>h_c + nadh_c + 3ohdcoa_c <=> nad_c + 3hhdcoa_c<br>3oodcoa_c + nadh_c + h_c <=> 3hodcoa_c + nad_c<br>nad_c + 3hadpcoa_c <=> oxadpcoa_c + nadh_c + h_c<br>octdp_c + 4hbz_c --> 3ophb_c + ppi_c<br>3hcinnm_p + h_p <=> h_c + 3hcinnm_c<br>3hcinnm_e <=> 3hcinnm_p<br>co2_c + h2o_c <=> hco3_c + h_c<br>hcys__L_c + amet_c --> met__L_c + ahcys_c + h_c<br>hcys__L_c + mmet_c --> 2.0 met__L_c + h_c<br>hdca_e --> hdca_p<br>hdcea_e --> hdcea_p<br>pheme_c + frdp_c + h2o_c --> hemeO_c + ppi_c<br>hhlipa_c + atp_c --> phhlipa_c + adp_c + h_c<br>hphhlipa_c + atp_c --> adp_c + phphhlipa_c + h_c<br>adphep_LD_c + lipa_c --> adp_c + hlipa_c + h_c<br>adphep_LD_c + hlipa_c --> hhlipa_c + adp_c + h_c<br>adphep_LD_c + phhlipa_c --> h_c + hphhlipa_c + adp_c<br>gggagicolipa_c + adphep_LD_c --> colipa_c + adp_c + h_c<br>4mhetz_c + atp_c --> adp_c + 4mpetz_c + h_c<br>glc__D_c + atp_c --> g6p_c + adp_c + h_c<br>atp_c + man_c --> man6p_c + adp_c + h_c<br>atp_c + fru_c --> adp_c + h_c + f6p_c<br>h_p + hxa_p <=> hxa_c + h_c<br>hg2_c + atp_c + h2o_c --> hg2_p + h_c + pi_c + adp_c<br>h_p + hg2_c --> hg2_p + h_c<br>hg2_e <=> hg2_p<br>histd_c + 2.0 nad_c + h2o_c --> his__L_c + 3.0 h_c + 2.0 nadh_c<br>hisp_c + h2o_c --> pi_c + histd_c<br>his__L_c + trnahis_c + atp_c --> histrna_c + ppi_c + amp_c<br>his__L_p + atp_c + h2o_c --> his__L_c + h_c + pi_c + adp_c<br>h_p + his__L_p <=> his__L_c + h_c<br>his__L_e <=> his__L_p<br>h2o_c + hkndd_c --> h_c + op4en_c + succ_c<br>hkntd_c + h2o_c --> fum_c + op4en_c + h_c<br>4.0 ppbng_c + h2o_c --> hmbil_c + 4.0 nh4_c<br>4ahmmp_c + atp_c --> adp_c + 4ampm_c + h_c<br>hom__L_c + h_p --> hom__L_p + h_c<br>hom__L_e <=> hom__L_p<br>4h2opntn_c --> pyr_c + acald_c<br>atp_c + 6hmhpt_c --> amp_c + 6hmhptpp_c + h_c<br>dhpppn_c + o2_c --> h_c + hkndd_c<br>h_p + 3hpppn_p <=> 3hpppn_c + h_c<br>3hpppn_e <=> 3hpppn_p<br>hpyr_c <=> 2h3oppan_c<br>hpyr_c + nadh_c + h_c --> nad_c + glyc__R_c<br>hpyr_c + nadph_c + h_c --> nadp_c + glyc__R_c<br>hom__L_c + nadp_c <=> aspsa_c + nadph_c + h_c<br>hom__L_c + atp_c --> adp_c + phom_c + h_c<br>hom__L_c + succoa_c --> suchms_c + coa_c<br>glu__L_c + imacp_c --> akg_c + hisp_c<br>hxan_c + nad_c + h2o_c --> h_c + xan_c + nadh_c<br>hxa_e <=> hxa_p<br>accoa_c + hxa_c --> hxcoa_c + ac_c<br>hxan_c + prpp_c --> imp_c + ppi_c<br>q8_c + 2.0 h_c + h2_c --> 2.0 h_p + q8h2_c<br>mqn8_c + h2_c + 2.0 h_c --> 2.0 h_p + mql8_c<br>2dmmq8_c + h2_c + 2.0 h_c --> 2dmmql8_c + 2.0 h_p<br>h2o_c + pyam5p_c --> pi_c + pydam_c<br>hxan_e <=> hxan_p<br>hxan_p <=> hxan_c<br>h_e <=> h_p<br>iscssh_c + 2fe1s_c + iscu_c --> iscs_c + 4.0 h_c + iscu_2fe2s_c<br>2.0 iscssh_c + 2.0 fe2_c + iscu_c + fadh2_c --> 2.0 iscs_c + 6.0 h_c + fad_c + iscu_2fe2s_c<br>2.0 iscssh_c + iscu_2fe2s_c + 2.0 fe2_c + fadh2_c --> iscu_2fe2s2_c + 2.0 iscs_c + 6.0 h_c + fad_c<br>iscu_2fe2s_c + 4.0 h_c --> 2fe2s_c + iscu_c<br>iscu_2fe2s2_c + fadh2_c + 2.0 h_c --> fad_c + iscu_4fe4s_c<br>iscu_4fe4s_c + 4.0 h_c --> 4fe4s_c + iscu_c<br>nadp_c + icit_c <=> nadph_c + akg_c + co2_c<br>chor_c <=> ichor_c<br>chor_c --> ichor_c<br>h2o_c + ichor_c --> pyr_c + 23ddhb_c<br>icit_c --> glx_c + succ_c<br>iscs_c + cys__L_c --> ala__L_c + iscssh_c<br>nadh_c + 5dglcn_c + h_c <=> nad_c + idon__L_c<br>nadph_c + 5dglcn_c + h_c --> nadp_c + idon__L_c<br>idon__L_p + h_p <=> idon__L_c + h_c<br>idon__L_e <=> idon__L_p<br>prlp_c + gln__L_c --> h_c + aicar_c + glu__L_c + eig3p_c<br>eig3p_c --> h2o_c + imacp_c<br>2cpr5p_c + h_c --> co2_c + 3ig3p_c + h2o_c<br>ile__L_c + akg_c <=> glu__L_c + 3mop_c<br>ile__L_c + trnaile_c + atp_c --> ppi_c + iletrna_c + amp_c<br>ile__L_p + atp_c + h2o_c --> h_c + ile__L_c + pi_c + adp_c<br>h_p + ile__L_p <=> ile__L_c + h_c<br>ile__L_e <=> ile__L_p<br>imp_c + h2o_c <=> fprica_c<br>imp_c + nad_c + h2o_c --> nadh_c + h_c + xmp_c<br>imp_e <=> imp_p<br>h_c + indole_c --> h_p + indole_p<br>h_p + indole_p <=> h_c + indole_c<br>indole_e <=> indole_p<br>inost_p + na1_p --> inost_c + na1_c<br>ins_c + h2o_c --> hxan_c + rib__D_c<br>ins_c + atp_c --> imp_c + adp_c + h_c<br>inost_e <=> inost_p<br>ins_p + h_p --> ins_c + h_c<br>ins_p + h_p <=> ins_c + h_c<br>ins_e <=> ins_p<br>ipdp_c <=> dmpp_c<br>h2mb4p_c + h_c + nadh_c --> nad_c + ipdp_c + h2o_c<br>nad_c + 3c2hmp_c --> nadh_c + h_c + 3c4mop_c<br>3c2hmp_c <=> 2ippm_c + h2o_c<br>h2o_c + 2ippm_c <=> 3c3hmp_c<br>3mob_c + accoa_c + h2o_c --> 3c3hmp_c + h_c + coa_c<br>isetac_p + atp_c + h2o_c --> h_c + pi_c + adp_c + isetac_c<br>isetac_e <=> isetac_p<br>kdo2lipid4_c + atp_c + h2o_c --> kdo2lipid4_p + h_c + pi_c + adp_c<br>kdo2lipid4_p + atp_c + h2o_c --> h_c + kdo2lipid4_e + pi_c + adp_c<br>nadp_c + 23dhmb_c <=> alac__S_c + nadph_c + h_c<br>2ahbut_c + nadph_c + h_c <=> nadp_c + 23dhmp_c<br>h_c + acACP_c + malACP_c --> co2_c + actACP_c + ACP_c<br>h_c + malACP_c + accoa_c --> co2_c + actACP_c + coa_c<br>kdo_c + ctp_c --> ckdo_c + ppi_c<br>h2o_c + kdo8p_c --> pi_c + kdo_c<br>ara5p_c + pep_c + h2o_c --> pi_c + kdo8p_c<br>3dhgulnp_c + h_c --> xu5p__L_c + co2_c<br>k_p + atp_c + h2o_c --> k_c + h_c + pi_c + adp_c<br>h_p + k_p --> k_c + h_c<br>k_c + h_p --> k_p + h_c<br>k_e <=> k_p<br>lac__L_c + q8_c --> pyr_c + q8h2_c<br>lac__L_c + mqn8_c --> pyr_c + mql8_c<br>h_p + lac__L_p <=> lac__L_c + h_c<br>lac__L_e <=> lac__L_p<br>uLa4n_p + colipa_p --> acolipa_p + udcpp_p<br>lcts_c + h2o_c --> glc__D_c + gal_c<br>lcts_p + h2o_p --> glc__D_p + gal_p<br>h2o_c + LalaDgluMdap_c --> LalaDglu_c + 26dap__M_c<br>LalaDglu_e <=> LalaDglu_p<br>h_p + LalaDglu_p --> LalaDglu_c + h_c<br>LalaLglu_e <=> LalaLglu_p<br>LalaLglu_p + h_p --> LalaLglu_c + h_c<br>nadh_c + h_c + mthgxl_c --> lald__D_c + nad_c<br>mthgxl_c + nadph_c + h_c --> lald__L_c + nadp_c<br>LalaLglu_c + h2o_c --> ala__L_c + glu__L_c<br>lald__L_c + nad_c + h2o_c --> lac__L_c + 2.0 h_c + nadh_c<br>lald__D_c + h_c + nadh_c <=> nad_c + 12ppd__R_c<br>lald__L_c + nadh_c + h_c <=> nad_c + 12ppd__S_c<br>lcts_c + h_p --> lcts_p + h_c<br>lcts_e <=> lcts_p<br>lcts_p + h_p <=> lcts_c + h_c<br>nad_c + lac__D_c <=> h_c + nadh_c + pyr_c<br>lac__D_c + q8_c --> pyr_c + q8h2_c<br>glu__L_c + 4mop_c --> akg_c + leu__L_c<br>trnaleu_c + atp_c + leu__L_c --> leutrna_c + ppi_c + amp_c<br>leu__L_p + atp_c + h2o_c --> pi_c + h_c + leu__L_c + adp_c<br>leu__L_p + h_p <=> h_c + leu__L_c<br>leu__L_e <=> leu__L_p<br>gthrd_c + mthgxl_c --> lgt__S_c<br>lipa_cold_c + atp_c + h2o_c --> lipa_cold_p + h_c + pi_c + adp_c<br>h_e + hdca_e + colipa_e --> h2o_e + hacolipa_e<br>h_e + lipa_e + hdca_e --> halipa_e + h2o_e<br>lipoamp_c --> amp_c + lipopb_c<br>lipoate_c + atp_c --> ppi_c + lipoamp_c<br>lipa_c + atp_c + h2o_c --> lipa_p + h_c + pi_c + adp_c<br>lipa_p + atp_c + h2o_c --> lipa_e + h_c + pi_c + adp_c<br>h_c + ocACP_c --> ACP_c + octapb_c<br>h_c + 4fe4s_c + nad_c + 2.0 amet_c + octapb_c --> 2fe2s_c + 2.0 met__L_c + nadh_c + lipopb_c + 2.0 dad_5_c + 2.0 fe2_c<br>h_p + lipoate_p --> lipoate_c + h_c<br>lipoate_e <=> lipoate_p<br>lipidX_c + u23ga_c --> udp_c + lipidAds_c + h_c<br>h2o_p + 1ddecg3p_p --> ddca_p + h_p + glyc3p_p<br>1tdecg3p_p + h2o_p --> h_p + glyc3p_p + ttdca_p<br>1tdec7eg3p_p + h2o_p --> ttdcea_p + h_p + glyc3p_p<br>1hdecg3p_p + h2o_p --> glyc3p_p + hdca_p + h_p<br>1hdec9eg3p_p + h2o_p --> h_p + glyc3p_p + hdcea_p<br>1odecg3p_p + h2o_p --> h_p + glyc3p_p + ocdca_p<br>1odec11eg3p_p + h2o_p --> glyc3p_p + h_p + ocdcea_p<br>1agpe120_p + h2o_p --> ddca_p + g3pe_p + h_p<br>h2o_p + 1agpe140_p --> g3pe_p + h_p + ttdca_p<br>1agpe141_p + h2o_p --> ttdcea_p + h_p + g3pe_p<br>h2o_p + 1agpe160_p --> g3pe_p + hdca_p + h_p<br>1agpe161_p + h2o_p --> g3pe_p + h_p + hdcea_p<br>1agpe180_p + h2o_p --> g3pe_p + h_p + ocdca_p<br>h2o_p + 1agpe181_p --> g3pe_p + h_p + ocdcea_p<br>1agpg120_p + h2o_p --> ddca_p + h_p + g3pg_p<br>1agpg140_p + h2o_p --> h_p + g3pg_p + ttdca_p<br>1agpg141_p + h2o_p --> ttdcea_p + h_p + g3pg_p<br>1agpg160_p + h2o_p --> hdca_p + h_p + g3pg_p<br>1agpg161_p + h2o_p --> h_p + g3pg_p + hdcea_p<br>1agpg180_p + h2o_p --> h_p + g3pg_p + ocdca_p<br>1agpg181_p + h2o_p --> g3pg_p + h_p + ocdcea_p<br>2ddecg3p_c + h2o_c --> glyc3p_c + ddca_c + 2.0 h_c<br>h2o_c + 2tdecg3p_c --> glyc3p_c + 2.0 h_c + ttdca_c<br>h2o_c + 2tdec7eg3p_c --> ttdcea_c + glyc3p_c + 2.0 h_c<br>2hdecg3p_c + h2o_c --> 2.0 h_c + hdca_c + glyc3p_c<br>h2o_c + 2hdec9eg3p_c --> glyc3p_c + hdcea_c + 2.0 h_c<br>h2o_c + 2odecg3p_c --> glyc3p_c + 2.0 h_c + ocdca_c<br>2odec11eg3p_c + h2o_c --> glyc3p_c + ocdcea_c + 2.0 h_c<br>2agpe120_c + pg120_c --> g3pe_c + apg120_c<br>2agpe140_c + pg140_c --> apg140_c + g3pe_c<br>pg141_c + 2agpe141_c --> g3pe_c + apg141_c<br>pg160_c + 2agpe160_c --> g3pe_c + apg160_c<br>2agpe161_c + pg161_c --> g3pe_c + apg161_c<br>2agpe180_c + pg180_c --> g3pe_c + apg180_c<br>pg181_c + 2agpe181_c --> g3pe_c + apg181_c<br>pg120_c + 2agpg120_c --> apg120_c + g3pg_c<br>2agpg140_c + pg140_c --> apg140_c + g3pg_c<br>pg141_c + 2agpg141_c --> g3pg_c + apg141_c<br>pg160_c + 2agpg160_c --> g3pg_c + apg160_c<br>pg161_c + 2agpg161_c --> apg161_c + g3pg_c<br>2agpg180_c + pg180_c --> apg180_c + g3pg_c<br>pg181_c + 2agpg181_c --> g3pg_c + apg181_c<br>2agpe120_c + h2o_c --> g3pe_c + ddca_c + h_c<br>2agpe140_c + h2o_c --> g3pe_c + h_c + ttdca_c<br>h2o_c + 2agpe141_c --> g3pe_c + ttdcea_c + h_c<br>h2o_c + 2agpe160_c --> g3pe_c + hdca_c + h_c<br>2agpe161_c + h2o_c --> g3pe_c + hdcea_c + h_c<br>h2o_c + 2agpe180_c --> g3pe_c + ocdca_c + h_c<br>h2o_c + 2agpe181_c --> g3pe_c + ocdcea_c + h_c<br>h2o_c + 2agpg120_c --> g3pg_c + ddca_c + h_c<br>h2o_c + 2agpg140_c --> g3pg_c + ttdca_c + h_c<br>h2o_c + 2agpg141_c --> ttdcea_c + g3pg_c + h_c<br>h2o_c + 2agpg160_c --> g3pg_c + hdca_c + h_c<br>h2o_c + 2agpg161_c --> h_c + g3pg_c + hdcea_c<br>2agpg180_c + h2o_c --> ocdca_c + g3pg_c + h_c<br>2agpg181_c + h2o_c --> ocdcea_c + g3pg_c + h_c<br>nadp_c + ser__L_c <=> 2amsa_c + nadph_c + h_c<br>lys__L_c + h_c --> 15dap_c + co2_c<br>lys__L_c + atp_c + trnalys_c --> ppi_c + amp_c + lystrna_c<br>lys__L_p + atp_c + h2o_c --> lys__L_c + h_c + pi_c + adp_c<br>lys__L_p + h_p --> lys__L_c + h_c<br>lys__L_c + h_p --> lys__L_p + h_c<br>lys__L_e <=> lys__L_p<br>lyx__L_c --> xylu__L_c<br>h_p + lyx__L_p --> lyx__L_c + h_c<br>lyx__L_e <=> lyx__L_p<br>mnl1p_c + nad_c <=> nadh_c + f6p_c + h_c<br>h_c + malACP_c --> acACP_c + co2_c<br>malcoa_c + amet_c --> ahcys_c + malcoame_c<br>nad_c + mal__D_c --> nadh_c + pyr_c + co2_c<br>mal__D_p + 2.0 h_p --> mal__D_c + 2.0 h_c<br>mal__D_e <=> mal__D_p<br>accoa_c + glx_c + h2o_c --> mal__L_c + h_c + coa_c<br>malt_c + accoa_c <=> coa_c + acmalt_c<br>malthx_p + atp_c + h2o_c --> h_c + malthx_c + pi_c + adp_c<br>malthx_e --> malthx_p<br>maltpt_p + atp_c + h2o_c --> maltpt_c + h_c + pi_c + adp_c<br>maltpt_e --> maltpt_p<br>malttr_p + atp_c + h2o_c --> h_c + malttr_c + pi_c + adp_c<br>malttr_e --> malttr_p<br>maltttr_p + atp_c + h2o_c --> h_c + maltttr_c + adp_c + pi_c<br>maltttr_e --> maltttr_p<br>malt_p + atp_c + h2o_c --> h_c + malt_c + pi_c + adp_c<br>pep_c + malt_p --> pyr_c + malt6p_c<br>malt_e --> malt_p<br>mal__L_p + 2.0 h_p --> mal__L_c + 2.0 h_c<br>mal__L_p + 3.0 h_p --> mal__L_c + 3.0 h_c<br>mal__L_c + h_p --> mal__L_p + h_c<br>mal__L_e <=> mal__L_p<br>man1p_c + gdp_c + h_c --> pi_c + gdpmann_c<br>man6p_c <=> f6p_c<br>2.0 pi_c + man6p_p --> man6p_c + 2.0 pi_p<br>man6p_e <=> man6p_p<br>mana_c + nad_c <=> nadh_c + fruur_c + h_c<br>pep_c + manglyc_p --> pyr_c + man6pglyc_c<br>manglyc_e <=> manglyc_p<br>man6pglyc_c + h2o_c --> man6p_c + glyc__R_c<br>pep_c + man_p --> pyr_c + man6p_c<br>man_e <=> man_p<br>2mcit_c --> 2mcacn_c + h2o_c<br>micit_c <=> pyr_c + succ_c<br>oaa_c + ppcoa_c + h2o_c --> h_c + coa_c + 2mcit_c<br>ACP_c + malcoa_c <=> coa_c + malACP_c<br>mercppyr_c + cyan_c --> tcynt_c + pyr_c + h_c<br>murein5p5p_p --> murein5px4p_p + ala__D_p<br>murein5p5p_p --> alaala_p + murein5px3p_p<br>murein5p5p5p_p --> murein5px4px4p_p + 2.0 ala__D_p<br>murein5px4p_p + h2o_p --> murein4px4p_p + ala__D_p<br>murein5px4px4p_p + h2o_p --> murein4px4px4p_p + ala__D_p<br>h2o_p + murein5p5p_p --> murein5p4p_p + ala__D_p<br>h2o_p + murein5p4p_p --> murein4p4p_p + ala__D_p<br>murein5p3p_p + h2o_p --> murein4p3p_p + ala__D_p<br>murein4px4p_p + h2o_p --> murein4p4p_p<br>murein3px4p_p + h2o_p --> murein4p3p_p<br>murein5px4p_p + h2o_p --> murein5p4p_p<br>murein4px4px4p_p + h2o_p --> murein4px4p4p_p<br>nad_c + mal__L_c <=> oaa_c + nadh_c + h_c<br>mal__L_c + q8_c --> oaa_c + q8h2_c<br>mal__L_c + mqn8_c --> oaa_c + mql8_c<br>nad_c + mal__L_c --> nadh_c + pyr_c + co2_c<br>nadp_c + mal__L_c --> nadph_c + pyr_c + co2_c<br>2mecdp_c + h_c + 2.0 flxr_c --> 2.0 flxso_c + h2mb4p_c + h2o_c<br>2p4c2me_c --> 2mecdp_c + cmp_c<br>melib_p + h_p --> melib_c + h_c<br>h_p + melib_c --> melib_p + h_c<br>melib_e <=> melib_p<br>meoh_e <=> meoh_p<br>meoh_p <=> meoh_c<br>h_c + 2me4p_c + ctp_c --> 4c2me_c + ppi_c<br>met__L_c + atp_c + h2o_c --> ppi_c + pi_c + amet_c<br>h2o_c + met__D_p + atp_c --> adp_c + h_c + pi_c + met__D_c<br>met__D_e <=> met__D_p<br>h2o2_c + met__L_c --> metsox_S__L_c + h2o_c<br>h2o2_c + met__L_c --> metsox_R__L_c + h2o_c<br>5mthf_c + hcys__L_c --> met__L_c + thf_c + h_c<br>metsox_S__L_p + atp_c + h2o_c --> metsox_S__L_c + h_c + pi_c + adp_c<br>metsox_S__L_e <=> metsox_S__L_p<br>metsox_R__L_p + atp_c + h2o_c --> h_c + pi_c + adp_c + metsox_R__L_c<br>metsox_R__L_e <=> metsox_R__L_p<br>trdrd_c + metsox_S__L_c --> met__L_c + trdox_c + h2o_c<br>trdrd_c + metsox_R__L_c --> met__L_c + trdox_c + h2o_c<br>trnamet_c + met__L_c + atp_c --> ppi_c + mettrna_c + amp_c<br>met__L_p + atp_c + h2o_c --> h_c + met__L_c + pi_c + adp_c<br>met__L_e <=> met__L_p<br>mg2_p + 2.0 h_c <=> mg2_c + 2.0 h_p<br>mg2_e <=> mg2_p<br>mg2_p --> mg2_c<br>mg2_p + atp_c + h2o_c --> mg2_c + h_c + pi_c + adp_c<br>dhap_c --> pi_c + mthgxl_c<br>mi1p__D_c + h2o_c --> inost_c + pi_c<br>2mcacn_c + h2o_c <=> micit_c<br>mincyc_e <=> mincyc_p<br>mincyc_p + h_p --> mincyc_e + h_c<br>minohp_e --> minohp_p<br>murein5px4p_p + h2o_p --> murein3px4p_p + alaala_p<br>murein4p4p_p + h2o_p --> murein4p3p_p + ala__D_p<br>h2o_p + murein5p5p_p --> murein5p3p_p + alaala_p<br>murein4p3p_p + h2o_p --> murein3p3p_p + ala__D_p<br>murein5px3p_p + h2o_p --> alaala_p + murein3px3p_p<br>murein3px3p_p + h2o_p --> murein3p3p_p<br>murein5px3p_p + h2o_p --> murein5p3p_p<br>malttr_c + h2o_c --> glc__D_c + malt_c<br>maltttr_c + h2o_c --> glc__D_c + malttr_c<br>maltpt_c + h2o_c --> glc__D_c + maltttr_c<br>h2o_c + malthx_c --> glc__D_c + maltpt_c<br>h2o_c + malthp_c --> glc__D_c + malthx_c<br>murein4p4p_p --> 2.0 anhgm4p_p<br>murein4p3p_p --> anhgm3p_p + anhgm4p_p<br>murein3p3p_p --> 2.0 anhgm3p_p<br>murein4px4p4p_p --> murein4px4p_p + anhgm4p_p<br>pi_c + maltpt_c <=> maltttr_c + g1p_c<br>pi_c + malthx_c <=> maltpt_c + g1p_c<br>malthp_c + pi_c <=> g1p_c + malthx_c<br>mmcoa__S_c + h_c --> co2_c + ppcoa_c<br>h_p + mmet_p --> h_c + mmet_c<br>mmet_e <=> mmet_p<br>succoa_c --> mmcoa__S_c<br>mn2_c + h_p --> mn2_p + h_c<br>mn2_p --> mn2_c<br>man6p_c + h2o_c --> pi_c + man_c<br>pep_c + mnl_p --> mnl1p_c + pyr_c<br>mnl_e <=> mnl_p<br>mana_c --> h2o_c + 2ddglcn_c<br>h_p + mn2_p --> mn2_c + h_c<br>mn2_e <=> mn2_p<br>iscssh_c + moadamp_c + nadh_c --> iscs_c + moadcosh_c + nad_c + amp_c<br>lipidA_c + ckdo_c --> kdolipid4_c + cmp_c + h_c<br>ckdo_c + kdolipid4_c --> kdo2lipid4_c + cmp_c + h_c<br>phphhlipa_c + ckdo_c --> cmp_c + kphphhlipa_c + h_c<br>mobd_p + atp_c + h2o_c --> mobd_c + pi_c + h_c + adp_c<br>mobd_e <=> mobd_p<br>h_c + moco_c + ctp_c --> mococdp_c + ppi_c<br>mptamp_c + mobd_c + 2.0 h_c --> amp_c + cu2_c + moco_c + h2o_c<br>gtp_c + moco_c + h_c --> mocogdp_c + ppi_c<br>3mob_c + mlthf_c + h2o_c --> 2dhp_c + thf_c<br>mal__L_c + o2_c <=> h2o2_c + oaa_c<br>mpt_c + atp_c + h_c --> mptamp_c + ppi_c<br>2.0 uaagmda_c --> 2.0 udcpdp_c + murein5p5p_p + 2.0 h_c<br>uaagmda_c + murein5p5p_p --> murein5p5p5p_p + udcpdp_c + h_c<br>cpmp_c + 2.0 moadcosh_c + cu2_c --> 2.0 moadcoo_c + 5.0 h_c + mpt_c<br>moadcoo_c + atp_c + h_c --> moadamp_c + ppi_c<br>msa_c + nadph_c + h_c --> 3hpp_c + nadp_c<br>mso3_p + atp_c + h2o_c --> mso3_c + h_c + pi_c + adp_c<br>mso3_e <=> mso3_p<br>h2o_c + 5mta_c --> ade_c + 5mtr_c<br>methf_c + h2o_c <=> 10fthf_c + h_c<br>nadp_c + mlthf_c <=> methf_c + nadph_c<br>nadh_c + mlthf_c + 2.0 h_c --> 5mthf_c + nad_c<br>mdhdhf_c + h2o_c --> mththf_c<br>Nmtrp_c + o2_c + h2o_c --> trp__L_c + h2o2_c + fald_c<br>n2o_e <=> n2o_p<br>n2o_p <=> n2o_c<br>acg5sa_c + h2o_c --> ac_c + glu5sa_c<br>nac_e <=> nac_p<br>nac_p --> nac_c<br>nad_c + h2o_c --> amp_c + 2.0 h_c + nmn_c<br>nadh_c + mqn8_c + h_c --> mql8_c + nad_c<br>4.0 h_c + nadh_c + q8_c --> 3.0 h_p + nad_c + q8h2_c<br>mqn8_c + 4.0 h_c + nadh_c --> mql8_c + nad_c + 3.0 h_p<br>4.0 h_c + nadh_c + 2dmmq8_c --> 2dmmql8_c + nad_c + 3.0 h_p<br>nadh_c + q8_c + h_c --> nad_c + q8h2_c<br>nadh_c + 2dmmq8_c + h_c --> 2dmmql8_c + nad_c<br>nad_c + atp_c --> nadp_c + adp_c + h_c<br>nad_c + h2o_c --> ncam_c + adprib_c + h_c<br>h_c + q8_c + nadph_c --> nadp_c + q8h2_c<br>mqn8_c + nadph_c + h_c --> nadp_c + mql8_c<br>nadph_c + 2dmmq8_c + h_c --> 2dmmql8_c + nadp_c<br>nadp_c + h2o_c --> pi_c + nad_c<br>dnad_c + nh4_c + atp_c --> ppi_c + h_c + nad_c + amp_c<br>nadph_c + nad_c --> nadp_c + nadh_c<br>nac_c + h2o_c + prpp_c + atp_c --> ppi_c + nicrnt_c + pi_c + adp_c<br>3.0 h_p + 2.0 na1_c --> 3.0 h_c + 2.0 na1_p<br>2.0 h_p + na1_c --> 2.0 h_c + na1_p<br>h_p + na1_c --> h_c + na1_p<br>na1_e <=> na1_p<br>gdp_c + atp_c <=> gtp_c + adp_c<br>udp_c + atp_c <=> utp_c + adp_c<br>cdp_c + atp_c <=> adp_c + ctp_c<br>dtdp_c + atp_c <=> dttp_c + adp_c<br>dgdp_c + atp_c <=> dgtp_c + adp_c<br>dudp_c + atp_c <=> dutp_c + adp_c<br>atp_c + dcdp_c <=> dctp_c + adp_c<br>dadp_c + atp_c <=> datp_c + adp_c<br>nh4_e <=> nh4_p<br>nh4_p <=> nh4_c<br>h_c + 2.0 no_c + nadh_c --> n2o_c + nad_c + h2o_c<br>ni2_c + atp_c + h2o_c --> ni2_p + h_c + pi_c + adp_c<br>ni2_c + h_p --> ni2_p + h_c<br>ni2_e <=> ni2_p<br>ni2_p --> ni2_c<br>ni2_p + atp_c + h2o_c --> ni2_c + h_c + pi_c + adp_c<br>nmn_c + atp_c + h_c --> nad_c + ppi_c<br>nmn_c + h2o_c --> nicrnt_c + nh4_c<br>nmn_c + h2o_c --> h_c + ncam_c + r5p_c<br>nmn_p --> nmn_c<br>h2o_c + nmn_p --> h_c + r5p_c + ncam_c<br>nmn_e <=> nmn_p<br>ncam_c + h2o_c --> nh4_c + nac_c<br>nicrnt_c + atp_c + h_c <=> dnad_c + ppi_c<br>dmbzid_c + nicrnt_c --> nac_c + 5prdmbz_c + h_c<br>quln_c + 2.0 h_c + prpp_c --> ppi_c + nicrnt_c + co2_c<br>no2_p + h_p <=> no2_c + h_c<br>no2_e <=> no2_p<br>no3_p + q8h2_c --> h2o_p + no2_p + q8_c<br>no3_c + 2.0 h_c + q8h2_c --> no2_c + 2.0 h_p + h2o_c + q8_c<br>mql8_c + no3_p --> no2_p + mqn8_c + h2o_p<br>no3_c + 2.0 h_c + mql8_c --> no2_c + mqn8_c + 2.0 h_p + h2o_c<br>no2_c + no3_p --> no3_c + no2_p<br>no3_e <=> no3_p<br>2.0 no_c + nadh_c + 2.0 o2_c --> 2.0 no3_c + h_c + nad_c<br>2.0 no_c + 2.0 o2_c + nadph_c --> 2.0 no3_c + nadp_c + h_c<br>novbcn_e <=> novbcn_p<br>novbcn_p + h_p --> novbcn_e + h_c<br>no_e <=> no_p<br>no_p <=> no_c<br>dump_c + h2o_c --> duri_c + pi_c<br>xmp_c + h2o_c --> pi_c + xtsn_c<br>xmp_p + h2o_p --> pi_p + xtsn_p<br>imp_c + h2o_c --> pi_c + ins_c<br>imp_p + h2o_p --> ins_p + pi_p<br>dimp_c + h2o_c --> pi_c + din_c<br>dimp_p + h2o_p --> din_p + pi_p<br>dump_p + h2o_p --> duri_p + pi_p<br>ump_c + h2o_c --> pi_c + uri_c<br>h2o_p + ump_p --> pi_p + uri_p<br>dcmp_c + h2o_c --> dcyt_c + pi_c<br>h2o_p + dcmp_p --> pi_p + dcyt_p<br>cmp_c + h2o_c --> pi_c + cytd_c<br>h2o_p + cmp_p --> cytd_p + pi_p<br>dtmp_c + h2o_c --> pi_c + thymd_c<br>dtmp_p + h2o_p --> thymd_p + pi_p<br>damp_c + h2o_c --> pi_c + dad_2_c<br>damp_p + h2o_p --> dad_2_p + pi_p<br>h2o_c + amp_c --> pi_c + adn_c<br>amp_p + h2o_p --> pi_p + adn_p<br>dgmp_c + h2o_c --> dgsn_c + pi_c<br>dgmp_p + h2o_p --> dgsn_p + pi_p<br>gmp_c + h2o_c --> pi_c + gsn_c<br>gmp_p + h2o_p --> pi_p + gsn_p<br>h2o_c + atp_c --> pi_c + adp_c + h_c<br>itp_c + h2o_c --> pi_c + idp_c + h_c<br>h2o_c + ditp_c --> pi_c + didp_c + h_c<br>h2o_c + xtp_c --> xdp_c + pi_c + h_c<br>h2o_c + gtp_c --> pi_c + gdp_c + h_c<br>gtp_p + h2o_p --> pi_p + h_p + gdp_p<br>ctp_c + h2o_c --> pi_c + cdp_c + h_c<br>dgtp_c + h2o_c --> dgmp_c + ppi_c + h_c<br>h2o_c + ditp_c --> dimp_c + ppi_c + h_c<br>h2o_c + xtp_c --> xmp_c + ppi_c + h_c<br>h2o_c + gtp_c --> gmp_c + ppi_c + h_c<br>dctp_c + h2o_c --> dcmp_c + ppi_c + h_c<br>h2o_c + ctp_c --> h_c + ppi_c + cmp_c<br>h2o_c + datp_c --> ppi_c + damp_c + h_c<br>atp_c + h2o_c --> amp_c + h_c + ppi_c<br>h2o_c + dttp_c --> dtmp_c + ppi_c + h_c<br>h2o_c + utp_c --> ump_c + ppi_c + h_c<br>itp_c + h2o_c --> imp_c + ppi_c + h_c<br>dgtp_c + h2o_c --> pppi_c + dgsn_c<br>gtp_c + h2o_c --> pppi_c + gsn_c<br>no2_c + 5.0 h_c + 3.0 nadh_c --> 3.0 nad_c + nh4_c + 2.0 h2o_c<br>no2_p + 2.0 h_p + 3.0 q8h2_c --> nh4_p + 2.0 h2o_p + 3.0 q8_c<br>no2_p + 3.0 mql8_c + 2.0 h_p --> 3.0 mqn8_c + 2.0 h2o_p + nh4_p<br>o16a4colipa_p + atp_c + h2o_c --> o16a4colipa_e + h_c + pi_c + adp_c<br>o16a4und_p + colipa_p --> h_p + o16a4colipa_p + udcpdp_p<br>2.0 o16aund_p --> h_p + o16a2und_p + udcpdp_p<br>o16a2und_p + o16aund_p --> h_p + udcpdp_p + o16a3und_p<br>o16aund_p + o16a3und_p --> o16a4und_p + h_p + udcpdp_p<br>accoa_c + ragund_c --> coa_c + aragund_c<br>o16aund_c --> o16aund_p<br>garagund_c + udpgalfur_c --> gfgaragund_c + udp_c + h_c<br>aragund_c + udpg_c --> garagund_c + udp_c + h_c<br>gfgaragund_c + udpg_c --> h_c + udp_c + o16aund_c<br>o2s_e <=> o2s_p<br>o2_e <=> o2_p<br>o2_p <=> o2_c<br>oaa_c + h_c --> pyr_c + co2_c<br>2obut_c + coa_c --> ppcoa_c + for_c<br>cbp_c + orn_c <=> pi_c + h_c + citr__L_c<br>ocdca_e --> ocdca_p<br>ocdcea_e --> ocdcea_p<br>octa_e <=> octa_p<br>frdp_c + 5.0 ipdp_c --> octdp_c + 5.0 ppi_c<br>h_c + atp_c + octa_c --> ppi_c + amp_c + octapb_c<br>hgmeACP_c --> egmeACP_c + h2o_c<br>ogmeACP_c + nadph_c + h_c --> nadp_c + hgmeACP_c<br>h_c + malACP_c + malcoame_c --> ogmeACP_c + coa_c + co2_c<br>ohpb_c + glu__L_c <=> akg_c + phthr_c<br>amet_c + 2ohph_c --> 2omph_c + ahcys_c + h_c<br>2ombzl_c + amet_c --> 2ommbl_c + ahcys_c + h_c<br>h_c + 3c4mop_c --> co2_c + 4mop_c<br>0.5 o2_c + 2ommbl_c --> 2omhmbl_c<br>3.0 h2o_c + nad_c + 2ommbl_c + 2.0 atp_c --> 3.0 h_c + 2.0 pi_c + 2omhmbl_c + nadh_c + 2.0 adp_c<br>orot5p_c + h_c --> co2_c + ump_c<br>2omph_c + 0.5 o2_c --> 2ombzl_c<br>2omph_c + 3.0 h2o_c + nad_c + 2.0 atp_c --> 3.0 h_c + 2ombzl_c + 2.0 pi_c + nadh_c + 2.0 adp_c<br>op4en_c + h2o_c --> 4h2opntn_c<br>3ophb_c + h_c --> co2_c + 2oph_c<br>0.5 o2_c + 2oph_c --> 2ohph_c<br>2oph_c + 3.0 h2o_c + nad_c + 2.0 atp_c --> 2ohph_c + 3.0 h_c + 2.0 pi_c + nadh_c + 2.0 adp_c<br>hpmeACP_c --> epmeACP_c + h2o_c<br>nadph_c + opmeACP_c + h_c --> nadp_c + hpmeACP_c<br>gmeACP_c + h_c + malACP_c --> co2_c + opmeACP_c + ACP_c<br>orn_c + h_c --> co2_c + ptrc_c<br>orn_p + atp_c + h2o_c --> h_c + orn_c + pi_c + adp_c<br>orn_e <=> orn_p<br>orot_p + 2.0 h_p --> orot_c + 2.0 h_c<br>orot_e <=> orot_p<br>ppi_c + orot5p_c <=> prpp_c + orot_c<br>pi_c + oxur_c --> cbp_c + oxam_c<br>oxalcoa_c + h_c --> forcoa_c + co2_c<br>nadp_c + 2oxpaccoa_c + 2.0 h2o_c --> 3oxdhscoa_c + 2.0 h_c + nadph_c<br>coa_c + 3oxdhscoa_c --> accoa_c + 23dhacoa_c<br>nad_c + 1pyr5c_c + 2.0 h2o_c --> h_c + nadh_c + glu__L_c<br>1pyr5c_c + nadph_c + 2.0 h_c --> pro__L_c + nadp_c<br>pa120_c + atp_c + h2o_c --> h_c + pi_c + adp_c + pa120_p<br>pa140_c + atp_c + h2o_c --> h_c + pi_c + adp_c + pa140_p<br>pa141_c + h2o_c + atp_c --> h_c + pi_c + adp_c + pa141_p<br>pa160_c + atp_c + h2o_c --> pa160_p + h_c + pi_c + adp_c<br>h2o_c + atp_c + pa161_c --> pa161_p + h_c + pi_c + adp_c<br>pa180_c + atp_c + h2o_c --> h_c + pi_c + adp_c + pa180_p<br>pa181_c + atp_c + h2o_c --> h_c + pi_c + adp_c + pa181_p<br>pacald_p + h_p <=> pacald_c + h_c<br>pacald_e <=> pacald_p<br>phaccoa_c + h_c + o2_c + nadph_c --> nadp_c + rephaccoa_c + h2o_c<br>pac_c + coa_c + atp_c --> ppi_c + phaccoa_c + amp_c<br>ala_B_c + pant__R_c + atp_c --> pnto__R_c + ppi_c + h_c + amp_c<br>pa120_c + h2o_c --> pi_c + 12dgr120_c<br>pa120_p + h2o_p --> 12dgr120_p + pi_p<br>h2o_c + pa140_c --> 12dgr140_c + pi_c<br>pa140_p + h2o_p --> 12dgr140_p + pi_p<br>pa141_c + h2o_c --> pi_c + 12dgr141_c<br>h2o_p + pa141_p --> 12dgr141_p + pi_p<br>pa160_c + h2o_c --> pi_c + 12dgr160_c<br>pa160_p + h2o_p --> 12dgr160_p + pi_p<br>h2o_c + pa161_c --> pi_c + 12dgr161_c<br>pa161_p + h2o_p --> 12dgr161_p + pi_p<br>pa180_c + h2o_c --> 12dgr180_c + pi_c<br>h2o_p + pa180_p --> 12dgr180_p + pi_p<br>pa181_c + h2o_c --> pi_c + 12dgr181_c<br>pa181_p + h2o_p --> pi_p + 12dgr181_p<br>udcpp_c + ugmda_c --> ump_c + uagmda_c<br>paps_c + trdrd_c --> trdox_c + so3_c + 2.0 h_c + pap_c<br>paps_c + grxrd_c --> so3_c + grxox_c + 2.0 h_c + pap_c<br>h2o_c + camp_c --> amp_c + h_c<br>35cgmp_c + h2o_c --> gmp_c + h_c<br>pyr_c + coa_c + nad_c --> co2_c + nadh_c + accoa_c<br>nad_c + pdx5p_c --> nadh_c + pydx5p_c + h_c<br>o2_c + pdx5p_c --> h2o2_c + pydx5p_c<br>phthr_c + nad_c + dxyl5p_c --> co2_c + pdx5p_c + h_c + nadh_c + pi_c + 2.0 h2o_c<br>pdx5p_c + h2o_c --> pi_c + pydxn_c<br>pe120_c + atp_c + h2o_c --> pe120_p + h_c + pi_c + adp_c<br>pe140_c + atp_c + h2o_c --> pe140_p + h_c + pi_c + adp_c<br>pe141_c + atp_c + h2o_c --> pe141_p + pi_c + h_c + adp_c<br>pe160_c + atp_c + h2o_c --> h_c + pi_c + pe160_p + adp_c<br>pe161_c + atp_c + h2o_c --> h_c + pe161_p + pi_c + adp_c<br>pe180_c + atp_c + h2o_c --> pe180_p + h_c + pi_c + adp_c<br>pe181_c + atp_c + h2o_c --> h_c + pi_c + adp_c + pe181_p<br>h2o_p + peamn_p + o2_p --> pacald_p + nh4_p + h2o2_p<br>peamn_e <=> peamn_p<br>nad_c + 4per_c <=> nadh_c + ohpb_c + h_c<br>pe161_p + lipa_p --> 12dgr161_p + enlipa_p<br>lipa_p + pe181_p --> enlipa_p + 12dgr181_p<br>f6p_c + atp_c --> fdp_c + adp_c + h_c<br>tag6p__D_c + atp_c --> tagdp__D_c + adp_c + h_c<br>s7p_c + atp_c --> adp_c + s17bp_c + h_c<br>pyr_c + coa_c --> accoa_c + for_c<br>pg120_c + atp_c + h2o_c --> pg120_p + h_c + pi_c + adp_c<br>pg140_c + atp_c + h2o_c --> pg140_p + h_c + pi_c + adp_c<br>pg141_c + atp_c + h2o_c --> h_c + pi_c + adp_c + pg141_p<br>pg160_c + atp_c + h2o_c --> pg160_p + h_c + pi_c + adp_c<br>pg161_c + h2o_c + atp_c --> h_c + pi_c + pg161_p + adp_c<br>pg180_c + atp_c + h2o_c --> h_c + pi_c + adp_c + pg180_p<br>pg181_c + atp_c + h2o_c --> h_c + pi_c + pg181_p + adp_c<br>gam1p_c <=> gam6p_c<br>nad_c + 3pg_c --> nadh_c + 3php_c + h_c<br>g6p_c <=> f6p_c<br>3pg_c + atp_c <=> 13dpg_c + adp_c<br>h2o_c + 6pgl_c --> 6pgc_c + h_c<br>h2o_c + 2pglyc_c --> pi_c + glyclt_c<br>2pg_c <=> 3pg_c<br>g1p_c <=> g6p_c<br>pgp120_c + atp_c + h2o_c --> h_c + pi_c + adp_c + pgp120_p<br>pgp140_c + atp_c + h2o_c --> h_c + pgp140_p + pi_c + adp_c<br>pgp141_c + atp_c + h2o_c --> h_c + pi_c + pgp141_p + adp_c<br>pgp160_c + atp_c + h2o_c --> h_c + pi_c + adp_c + pgp160_p<br>pgp161_c + atp_c + h2o_c --> pgp161_p + h_c + pi_c + adp_c<br>pgp180_c + atp_c + h2o_c --> pgp180_p + h_c + pi_c + adp_c<br>pgp181_c + atp_c + h2o_c --> h_c + pgp181_p + pi_c + adp_c<br>pgp120_c + h2o_c --> pg120_c + pi_c<br>pgp120_p + h2o_p --> pg120_p + pi_p<br>pgp140_c + h2o_c --> pi_c + pg140_c<br>pgp140_p + h2o_p --> pg140_p + pi_p<br>h2o_c + pgp141_c --> pg141_c + pi_c<br>pgp141_p + h2o_p --> pi_p + pg141_p<br>pgp160_c + h2o_c --> pg160_c + pi_c<br>pgp160_p + h2o_p --> pg160_p + pi_p<br>pgp161_c + h2o_c --> pi_c + pg161_c<br>h2o_p + pgp161_p --> pg161_p + pi_p<br>pgp180_c + h2o_c --> pi_c + pg180_c<br>pgp180_p + h2o_p --> pi_p + pg180_p<br>pgp181_c + h2o_c --> pg181_c + pi_c<br>pgp181_p + h2o_p --> pi_p + pg181_p<br>glyc3p_c + cdpdddecg_c --> pgp120_c + cmp_c + h_c<br>cdpdtdecg_c + glyc3p_c --> pgp140_c + cmp_c + h_c<br>glyc3p_c + cdpdtdec7eg_c --> h_c + cmp_c + pgp141_c<br>glyc3p_c + cdpdhdecg_c --> pgp160_c + cmp_c + h_c<br>glyc3p_c + cdpdhdec9eg_c --> pgp161_c + cmp_c + h_c<br>cdpdodecg_c + glyc3p_c --> pgp180_c + cmp_c + h_c<br>glyc3p_c + cdpdodec11eg_c --> pgp181_c + cmp_c + h_c<br>pheme_c + atp_c + h2o_c --> h_c + pi_c + pheme_p + adp_c<br>pheme_p --> pheme_e<br>phe__L_c + akg_c <=> phpyr_c + glu__L_c<br>phe__L_c + trnaphe_c + atp_c --> ppi_c + amp_c + phetrna_c<br>h_p + phe__L_p <=> phe__L_c + h_c<br>phe__L_e <=> phe__L_p<br>minohp_p + 6.0 h2o_p --> 6.0 pi_p + inost_p<br>h_p + pi_p <=> pi_c + h_c<br>pi_e <=> pi_p<br>pi_p + atp_c + h2o_c --> h_c + 2.0 pi_c + adp_c<br>pa120_p + h2o_p --> ddca_p + 2ddecg3p_p<br>pa140_p + h2o_p --> ttdca_p + 2tdecg3p_p<br>h2o_p + pa141_p --> ttdcea_p + 2tdec7eg3p_p<br>h2o_p + pa160_p --> 2hdecg3p_p + hdca_p<br>pa161_p + h2o_p --> 2hdec9eg3p_p + hdcea_p<br>h2o_p + pa180_p --> ocdca_p + 2odecg3p_p<br>pa181_p + h2o_p --> 2odec11eg3p_p + ocdcea_p<br>pe120_p + h2o_p --> ddca_p + h_p + 2agpe120_p<br>h2o_p + pe140_p --> h_p + 2agpe140_p + ttdca_p<br>h2o_p + pe141_p --> ttdcea_p + 2agpe141_p + h_p<br>pe160_p + h2o_p --> hdca_p + h_p + 2agpe160_p<br>pe161_p + h2o_p --> 2agpe161_p + h_p + hdcea_p<br>pe180_p + h2o_p --> 2agpe180_p + h_p + ocdca_p<br>pe181_p + h2o_p --> ocdcea_p + h_p + 2agpe181_p<br>pg120_p + h2o_p --> ddca_p + 2agpg120_p + h_p<br>pg140_p + h2o_p --> h_p + 2agpg140_p + ttdca_p<br>h2o_p + pg141_p --> 2agpg141_p + ttdcea_p + h_p<br>pg160_p + h2o_p --> hdca_p + h_p + 2agpg160_p<br>pg161_p + h2o_p --> 2agpg161_p + h_p + hdcea_p<br>pg180_p + h2o_p --> 2agpg180_p + h_p + ocdca_p<br>h2o_p + pg181_p --> ocdcea_p + 2agpg181_p + h_p<br>pa120_p + h2o_p --> ddca_p + h_p + 1ddecg3p_p<br>pa140_p + h2o_p --> 1tdecg3p_p + h_p + ttdca_p<br>h2o_p + pa141_p --> ttdcea_p + 1tdec7eg3p_p + h_p<br>pa160_p + h2o_p --> hdca_p + 1hdecg3p_p + h_p<br>pa161_p + h2o_p --> 1hdec9eg3p_p + h_p + hdcea_p<br>pa180_p + h2o_p --> h_p + 1odecg3p_p + ocdca_p<br>pa181_p + h2o_p --> 1odec11eg3p_p + h_p + ocdcea_p<br>pe120_p + h2o_p --> ddca_p + h_p + 1agpe120_p<br>h2o_p + pe140_p --> h_p + 1agpe140_p + ttdca_p<br>h2o_p + pe141_p --> ttdcea_p + 1agpe141_p + h_p<br>pe160_p + h2o_p --> hdca_p + h_p + 1agpe160_p<br>pe161_p + h2o_p --> 1agpe161_p + h_p + hdcea_p<br>pe180_p + h2o_p --> h_p + 1agpe180_p + ocdca_p<br>pe181_p + h2o_p --> h_p + 1agpe181_p + ocdcea_p<br>pg120_p + h2o_p --> ddca_p + 1agpg120_p + h_p<br>pg140_p + h2o_p --> h_p + 1agpg140_p + ttdca_p<br>pg141_p + h2o_p --> ttdcea_p + h_p + 1agpg141_p<br>h2o_p + pg160_p --> h_p + hdca_p + 1agpg160_p<br>pg161_p + h2o_p --> h_p + hdcea_p + 1agpg161_p<br>pg180_p + h2o_p --> 1agpg180_p + h_p + ocdca_p<br>h2o_p + pg181_p --> ocdcea_p + 1agpg181_p + h_p<br>man1p_c <=> man6p_c<br>5aprbu_c + h2o_c --> 4r5au_c + pi_c<br>pmeACP_c + h2o_c --> meoh_c + pimACP_c<br>4ampm_c + atp_c --> 2mahmp_c + adp_c<br>pnto__R_c + atp_c --> adp_c + 4ppan_c + h_c<br>pnto__R_p + na1_p --> pnto__R_c + na1_c<br>pnto__R_e <=> pnto__R_p<br>nadh_c + poaac_c --> 3amac_c + nad_c + h2o_c<br>2.0 flxso_c + pyr_c + coa_c <=> h_c + accoa_c + 2.0 flxr_c + co2_c<br>h2o_c + pyr_c + q8_c --> co2_c + ac_c + q8h2_c<br>h2o_c + ppi_c --> 2.0 pi_c + h_c<br>pppi_c + h2o_c --> pi_c + ppi_c + h_c<br>ppap_c + adp_c <=> ppa_c + atp_c<br>ppal_e <=> ppal_p<br>ppal_c <=> ppal_p<br>ppa_p + na1_p --> ppa_c + na1_c<br>ppa_e <=> ppa_p<br>2.0 5aop_c --> h_c + ppbng_c + 2.0 h2o_c<br>pep_c + co2_c + h2o_c --> oaa_c + h_c + pi_c<br>4ppcys_c + h_c --> pan4p_c + co2_c<br>oaa_c + atp_c --> co2_c + pep_c + adp_c<br>ppcoa_c + succ_c --> ppa_c + succoa_c<br>h2o_c + ppgpp_c --> gdp_c + ppi_c<br>ppi_c + atp_c <=> pppi_c + adp_c<br>pi_c + atp_c <=> ppi_c + adp_c<br>r1p_c <=> r5p_c<br>2dr1p_c <=> 2dr5p_c<br>4ppan_c + ctp_c + cys__L_c --> cmp_c + ppi_c + h_c + 4ppcys_c<br>nad_c + pphn_c --> 34hpp_c + nadh_c + co2_c<br>pphn_c + h_c --> phpyr_c + co2_c + h2o_c<br>pppg9_c + 1.5 o2_c --> ppp9_c + 3.0 h2o_c<br>3.0 fum_c + pppg9_c --> ppp9_c + 3.0 succ_c<br>nadh_c + h_c + o2_c + pppn_c --> nad_c + cechddd_c<br>pppn_p + h_p <=> pppn_c + h_c<br>pppn_e <=> pppn_p<br>pyr_c + atp_c + h2o_c --> 2.0 h_c + amp_c + pi_c + pep_c<br>ppt_p + h2o_p --> h2_p + pi_p<br>ppt_e <=> ppt_p<br>pram_c + gly_c + atp_c <=> gar_c + h_c + pi_c + adp_c<br>fpram_c + atp_c --> air_c + 2.0 h_c + pi_c + adp_c<br>pran_c --> 2cpr5p_c<br>prbamp_c + h2o_c --> prfp_c<br>5aizc_c + asp__L_c + atp_c --> h_c + pi_c + 25aics_c + adp_c<br>h2o_c + prbatp_c --> h_c + ppi_c + prbamp_c<br>fgam_c + atp_c + h2o_c + gln__L_c --> fpram_c + h_c + pi_c + glu__L_c + adp_c<br>prfp_c <=> prlp_c<br>pro__L_c + fad_c --> 1pyr5c_c + fadh2_c + h_c<br>progly_p + atp_c + h2o_c --> h_c + pi_c + adp_c + progly_c<br>progly_e <=> progly_p<br>trnapro_c + pro__L_c + atp_c --> protrna_c + ppi_c + amp_c<br>pro__L_p + atp_c + h2o_c --> h_c + pro__L_c + pi_c + adp_c<br>h_p + pro__L_p <=> pro__L_c + h_c<br>pro__L_p + na1_p --> pro__L_c + na1_c<br>pro__L_e <=> pro__L_p<br>r5p_c + atp_c <=> amp_c + prpp_c + h_c<br>h_p + psclys_p --> psclys_c + h_c<br>psclys_e <=> psclys_p<br>pep_c + skm5p_c <=> 3psme_c + pi_c<br>ps120_c + h_c --> co2_c + pe120_c<br>ps140_c + h_c --> pe140_c + co2_c<br>ps141_c + h_c --> co2_c + pe141_c<br>ps160_c + h_c --> co2_c + pe160_c<br>h_c + ps161_c --> co2_c + pe161_c<br>ps180_c + h_c --> pe180_c + co2_c<br>ps181_c + h_c --> co2_c + pe181_c<br>glu__L_c + 3php_c --> akg_c + pser__L_c<br>pser__L_e <=> pser__L_p<br>pser__L_c + h2o_c --> pi_c + ser__L_c<br>pser__L_p + h2o_p --> ser__L_p + pi_p<br>ser__L_c + cdpdddecg_c --> ps120_c + cmp_c + h_c<br>cdpdtdecg_c + ser__L_c --> ps140_c + cmp_c + h_c<br>ser__L_c + cdpdtdec7eg_c --> h_c + ps141_c + cmp_c<br>cdpdhdecg_c + ser__L_c --> ps160_c + cmp_c + h_c<br>cdpdhdec9eg_c + ser__L_c --> cmp_c + ps161_c + h_c<br>cdpdodecg_c + ser__L_c --> ps180_c + cmp_c + h_c<br>ser__L_c + cdpdodec11eg_c --> ps181_c + cmp_c + h_c<br>pi_c + ppcoa_c --> ppap_c + coa_c<br>pi_c + accoa_c <=> coa_c + actp_c<br>thrp_p + h2o_p --> pi_p + thr__L_p<br>pan4p_c + atp_c + h_c --> dpcoa_c + ppi_c<br>orn_c + ptrc_p <=> ptrc_c + orn_p<br>akg_c + ptrc_c --> glu__L_c + 4abutn_c<br>ptrc_p + atp_c + h2o_c --> ptrc_c + h_c + pi_c + adp_c<br>h_p + ptrc_p --> ptrc_c + h_c<br>ptrc_e <=> ptrc_p<br>pi_c + adn_c <=> ade_c + r1p_c<br>pi_c + dad_2_c <=> 2dr1p_c + ade_c<br>gsn_c + pi_c <=> r1p_c + gua_c<br>dgsn_c + pi_c <=> 2dr1p_c + gua_c<br>pi_c + ins_c <=> hxan_c + r1p_c<br>pi_c + din_c <=> hxan_c + 2dr1p_c<br>pi_c + xtsn_c <=> xan_c + r1p_c<br>pyam5p_c + o2_c + h2o_c --> pydx5p_c + h2o2_c + nh4_c<br>pydam_c + atp_c --> pyam5p_c + adp_c + h_c<br>pydam_e <=> pydam_p<br>pydam_p --> pydam_c<br>pydx_c + atp_c --> pydx5p_c + adp_c + h_c<br>pydxn_c + atp_c --> h_c + pdx5p_c + adp_c<br>pydxn_e <=> pydxn_p<br>pydxn_p --> pydxn_c<br>pydx5p_c + h2o_c --> pydx_c + pi_c<br>pydx_e <=> pydx_p<br>pydx_p --> pydx_c<br>pep_c + adp_c + h_c --> pyr_c + atp_c<br>pi_c + uri_c <=> ura_c + r1p_c<br>ura_c + h_c + nadh_c + o2_c --> nad_c + uracp_c<br>h_p + pyr_p <=> pyr_c + h_c<br>pyr_e <=> pyr_p<br>2.0 o2_c + q8h2_c --> q8_c + 2.0 h_c + 2.0 o2s_c<br>2.0 o2_c + mql8_c --> mqn8_c + 2.0 o2s_c + 2.0 h_c<br>quin_e <=> quin_p<br>quin_p <=> quin_c<br>quin_c + nad_c <=> 3dhq_c + 2.0 h_c + nadh_c<br>dhap_c + iasp_c --> pi_c + quln_c + 2.0 h2o_c<br>atp_c + r15bp_c --> adp_c + prpp_c<br>r1p_c + atp_c --> r15bp_c + adp_c + h_c<br>r5p_c + h2o_c --> rib__D_c + pi_c<br>r5p_p + h2o_p --> rib__D_p + pi_p<br>r5p_e <=> r5p_p<br>ribflv_c + atp_c --> h_c + fmn_c + adp_c<br>4r5au_c + db4p_c --> pi_c + 2.0 h2o_c + dmlz_c<br>2.0 dmlz_c --> 4r5au_c + ribflv_c<br>rib__D_c + atp_c --> adp_c + r5p_c + h_c<br>rbl__L_c + atp_c --> adp_c + ru5p__L_c + h_c<br>ru5p__L_c <=> xu5p__D_c<br>rephaccoa_c <=> 2oxpaccoa_c<br>rfamp_e <=> rfamp_p<br>rfamp_p + h_p --> rfamp_e + h_c<br>dtdprmn_c + kphphhlipa_c --> dtdp_c + icolipa_c + h_c<br>rhcys_c --> dhptd_c + hcys__L_c<br>rib__D_p + atp_c + h2o_c --> rib__D_c + pi_c + h_c + adp_c<br>rib__D_e <=> rib__D_p<br>rmn_c <=> rml_c<br>rml_c + atp_c --> rml1p_c + adp_c + h_c<br>rmn_e <=> rmn_p<br>h_p + rmn_p --> rmn_c + h_c<br>rml1p_c <=> dhap_c + lald__L_c<br>trdrd_c + adp_c --> trdox_c + dadp_c + h2o_c<br>grxrd_c + adp_c --> dadp_c + grxox_c + h2o_c<br>trdrd_c + gdp_c --> trdox_c + dgdp_c + h2o_c<br>grxrd_c + gdp_c --> h2o_c + grxox_c + dgdp_c<br>trdrd_c + cdp_c --> trdox_c + h2o_c + dcdp_c<br>grxrd_c + cdp_c --> h2o_c + grxox_c + dcdp_c<br>trdrd_c + udp_c --> trdox_c + dudp_c + h2o_c<br>grxrd_c + udp_c --> dudp_c + grxox_c + h2o_c<br>2.0 h_c + 2.0 flxr_c + atp_c --> 2.0 flxso_c + datp_c + h2o_c<br>gtp_c + 2.0 h_c + 2.0 flxr_c --> 2.0 flxso_c + dgtp_c + h2o_c<br>ctp_c + 2.0 flxr_c + 2.0 h_c --> 2.0 flxso_c + dctp_c + h2o_c<br>utp_c + 2.0 h_c + 2.0 flxr_c --> 2.0 flxso_c + dutp_c + h2o_c<br>ru5p__D_c <=> xu5p__D_c<br>r5p_c <=> ru5p__D_c<br>5prdmbz_c + h2o_c --> pi_c + rdmbzi_c<br>sufbcd_c + h2o_c + 2fe1s_c + sufsesh_c + atp_c --> sufbcd_2fe2s_c + 5.0 h_c + sufse_c + pi_c + adp_c<br>2.0 fe2_c + 2.0 sufsesh_c + atp_c + h2o_c + sufbcd_c + fadh2_c --> 2.0 sufse_c + 7.0 h_c + pi_c + sufbcd_2fe2s_c + fad_c + adp_c<br>2.0 sufsesh_c + 2.0 fe2_c + atp_c + h2o_c + sufbcd_2fe2s_c + fadh2_c --> 2.0 sufse_c + 7.0 h_c + pi_c + sufbcd_2fe2s2_c + fad_c + adp_c<br>sufbcd_2fe2s_c + 4.0 h_c --> 2fe2s_c + sufbcd_c<br>sufbcd_2fe2s2_c + fadh2_c + 2.0 h_c --> sufbcd_4fe4s_c + fad_c<br>sufbcd_4fe4s_c + 4.0 h_c --> 4fe4s_c + sufbcd_c<br>s7p_c --> gmhep7p_c<br>2.0 h_c + sucarg_c + 2.0 h2o_c --> co2_c + 2.0 nh4_c + sucorn_c<br>gtp_c + so4_c + atp_c + h2o_c --> gdp_c + ppi_c + pi_c + aps_c<br>sarcs_c + o2_c + h2o_c --> fald_c + h2o2_c + gly_c<br>nad_c + sbt6p_c <=> nadh_c + f6p_c + h_c<br>sbt__D_p + pep_c --> pyr_c + sbt6p_c<br>sbt__D_e <=> sbt__D_p<br>sufse_c + cys__L_c --> ala__L_c + sufsesh_c<br>h2o_c + sl26da_c --> 26dap_LL_c + succ_c<br>akg_c + sl26da_c <=> sl2a6o_c + glu__L_c<br>sertrna_sec_c + selnp_c --> pi_c + sectrna_c + h_c<br>4.0 gthrd_c + 2.0 h_c + slnt_c --> gthox_c + dgslnt_c + 3.0 h2o_c<br>dgslnt_c + h_c + nadph_c --> gthrd_c + nadp_c + gslnt_c<br>nadph_c + gslnt_c --> gthrd_c + seln_c + nadp_c<br>seln_c + atp_c + h2o_c --> selnp_c + pi_c + amp_c<br>mql8_c + sel_c --> mqn8_c + slnt_c + h2o_c<br>sel_e <=> sel_p<br>sel_p + 2.0 h_p --> sel_c + 2.0 h_c<br>ichor_c + akg_c + h_c --> co2_c + 2sephchc_c<br>atp_c + ser__L_c + h_c <=> ppi_c + seramp_c<br>accoa_c + ser__L_c <=> acser_c + coa_c<br>ser__D_c --> pyr_c + nh4_c<br>ser__L_c --> pyr_c + nh4_c<br>ser__L_c + trnaser_c + atp_c --> ppi_c + sertrna_c + amp_c<br>ser__L_c + trnasecys_c + atp_c --> ppi_c + sertrna_sec_c + amp_c<br>ser__L_p + h_p <=> ser__L_c + h_c<br>ser__L_p + na1_p --> na1_c + ser__L_c<br>ser__L_e <=> ser__L_p<br>Sfglutth_c + h2o_c --> gthrd_c + h_c + for_c<br>sucglu_c + h2o_c --> glu__L_c + succ_c<br>sucgsa_c + nad_c + h2o_c --> 2.0 h_c + sucglu_c + nadh_c<br>2sephchc_c --> pyr_c + 2shchc_c<br>nad_c + dscl_c --> nadh_c + scl_c + h_c<br>scl_c + fe2_c --> sheme_c + 3.0 h_c<br>3dhsk_c + h_c + nadph_c <=> nadp_c + skm_c<br>skm_c + atp_c --> adp_c + skm5p_c + h_c<br>suchms_c + cys__L_c --> cyst__L_c + succ_c + h_c<br>skm_p + h_p --> skm_c + h_c<br>skm_e <=> skm_p<br>slnt_e <=> slnt_p<br>slnt_p + h_p --> slnt_c + h_c<br>so2_e <=> so2_p<br>so2_p <=> so2_c<br>so3_e <=> so3_p<br>h_p + so4_p --> h_c + so4_c<br>so4_e <=> so4_p<br>sucorn_c + akg_c --> glu__L_c + sucgsa_c<br>spmd_c + accoa_c --> N1aspmd_c + coa_c + h_c<br>spmd_c + accoa_c --> n8aspmd_c + coa_c + h_c<br>spmd_p + atp_c + h2o_c --> h_c + spmd_c + adp_c + pi_c<br>spmd_c + h_p --> spmd_p + h_c<br>spmd_e <=> spmd_p<br>ametam_c + ptrc_c --> spmd_c + 5mta_c + h_c<br>2.0 o2s_c + 2.0 h_c --> h2o2_c + o2_c<br>2.0 h_p + 2.0 o2s_p --> h2o2_p + o2_p<br>nad_c + sucsal_c + h2o_c --> nadh_c + succ_c + 2.0 h_c<br>nadp_c + sucsal_c + h2o_c --> succ_c + 2.0 h_c + nadph_c<br>asp__L_p + succ_c --> succ_p + asp__L_c<br>sucbz_c + coa_c + atp_c --> ppi_c + sbzcoa_c + amp_c<br>2shchc_c --> sucbz_c + h2o_c<br>succ_p + 2.0 h_p --> succ_c + 2.0 h_c<br>succ_p + 3.0 h_p --> succ_c + 3.0 h_c<br>h_p + succ_c --> succ_p + h_c<br>succ_e <=> succ_p<br>q8_c + succ_c --> fum_c + q8h2_c<br>fum_p + succ_c --> fum_c + succ_p<br>succ_c + mal__L_p --> mal__L_c + succ_p<br>succ_c + coa_c + atp_c <=> succoa_c + pi_c + adp_c<br>sucr_e <=> sucr_p<br>tartr__D_p + succ_c --> succ_p + tartr__D_c<br>pep_c + sucr_p --> suc6p_c + pyr_c<br>sulfac_p + h2o_c + atp_c --> sulfac_c + h_c + pi_c + adp_c<br>sulfac_e <=> sulfac_p<br>so3_c + 5.0 h_c + 3.0 nadph_c --> h2s_c + 3.0 nadp_c + 3.0 h2o_c<br>so4_p + atp_c + h2o_c --> so4_c + pi_c + h_c + adp_c<br>tdec2eACP_c <=> cdec3eACP_c<br>nad_c + altrn_c <=> h_c + nadh_c + tagur_c<br>s7p_c + g3p_c <=> e4p_c + f6p_c<br>tartr__L_c --> oaa_c + h2o_c<br>tartr__D_e <=> tartr__D_p<br>tartr__L_p + succ_c <=> succ_p + tartr__L_c<br>tartr__L_e <=> tartr__L_p<br>3.0 h_p + tartr__D_p --> tartr__D_c + 3.0 h_c<br>taur_c + akg_c + o2_c --> succ_c + co2_c + so3_c + h_c + aacald_c<br>taur_p + atp_c + h2o_c --> taur_c + h_c + pi_c + adp_c<br>taur_e <=> taur_p<br>tcynt_e <=> tcynt_p<br>h2o_c + thmpp_c --> pi_c + h_c + thmmp_c<br>dtdp4addg_c + accoa_c --> dtdp4aaddg_c + coa_c + h_c<br>dtdp4d6dg_c + glu__L_c --> dtdp4addg_c + akg_c<br>dtdp4d6dg_c --> dtdp4d6dm_c<br>dtdp4d6dm_c + nadph_c + h_c --> dtdprmn_c + nadp_c<br>dtdpglu_c --> dtdp4d6dg_c + h2o_c<br>lipidAds_c + atp_c --> h_c + adp_c + lipidA_c<br>dsbdrd_c + dsbcox_p --> dsbcrd_p + dsbdox_c<br>dsbdrd_c + dsbgox_p --> dsbgrd_p + dsbdox_c<br>tagdp__D_c <=> dhap_c + g3p_c<br>nadp_c + nadh_c + 2.0 h_p --> 2.0 h_c + nad_c + nadph_c<br>thdp_c + succoa_c + h2o_c --> sl2a6o_c + coa_c<br>methf_c + h2o_c --> 5fthf_c + h_c<br>h2o2_c + trdrd_c --> trdox_c + 2.0 h2o_c<br>thymd_p + h_p --> thymd_c + h_c<br>thymd_p + h_p <=> thymd_c + h_c<br>thymd_e <=> thymd_p<br>h2o_c + atp_c + thm_p --> thm_c + h_c + pi_c + adp_c<br>thm_e <=> thm_p<br>athr__L_c --> gly_c + acald_c<br>thr__L_c --> gly_c + acald_c<br>nad_c + thr__L_c --> 2aobut_c + nadh_c + h_c<br>thr__L_c --> 2obut_c + nh4_c<br>thrp_e <=> thrp_p<br>phom_c + h2o_c --> pi_c + thr__L_c<br>trnathr_c + thr__L_c + atp_c --> ppi_c + thrtrna_c + amp_c<br>thr__L_p + atp_c + h2o_c --> thr__L_c + h_c + pi_c + adp_c<br>thr__L_c + h_p --> thr__L_p + h_c<br>h_p + thr__L_p <=> thr__L_c + h_c<br>thr__L_p + na1_p --> thr__L_c + na1_c<br>thr__L_e <=> thr__L_p<br>thym_c + h_p --> thym_p + h_c<br>thym_e <=> thym_p<br>nadph_c + dxyl5p_c + iscssh_c + atp_c + dhgly_c + h_c --> ppi_c + co2_c + amp_c + nadp_c + iscs_c + 2.0 h2o_c + 4mpetz_c<br>xu5p__D_c + r5p_c <=> s7p_c + g3p_c<br>xu5p__D_c + e4p_c <=> g3p_c + f6p_c<br>tmao_c + h_c + mql8_c --> mqn8_c + h2o_c + tma_c<br>tmao_p + mql8_c + h_p --> tma_p + h2o_p + mqn8_c<br>tmao_c + h_c + 2dmmql8_c --> h2o_c + 2dmmq8_c + tma_c<br>2dmmql8_c + tmao_p + h_p --> tma_p + h2o_p + 2dmmq8_c<br>tmao_e <=> tmao_p<br>tma_e <=> tma_p<br>thymd_c + atp_c --> adp_c + dtmp_c + h_c<br>pi_c + thymd_c <=> 2dr1p_c + thym_c<br>dump_c + mlthf_c --> dhf_c + dtmp_c<br>thm_c + atp_c --> adp_c + thmmp_c + h_c<br>thmmp_c + atp_c --> thmpp_c + adp_c<br>2mahmp_c + 4mpetz_c + h_c --> thmmp_c + ppi_c<br>dhap_c <=> g3p_c<br>dpcoa_c + atp_c --> 2tpr3dpcoa_c + ade_c<br>trdox_c + nadph_c + h_c --> trdrd_c + nadp_c<br>tre6p_c + h2o_c --> glc__D_c + g6p_c<br>tre6p_c + h2o_c --> pi_c + tre_c<br>g6p_c + udpg_c --> tre6p_c + udp_c + h_c<br>h2o_c + tre_c --> 2.0 glc__D_c<br>tre_p + h2o_p --> 2.0 glc__D_p<br>tre_p + pep_c --> tre6p_c + pyr_c<br>tre_e <=> tre_p<br>trp__L_c + h2o_c <=> indole_c + pyr_c + nh4_c<br>3ig3p_c + ser__L_c --> g3p_c + trp__L_c + h2o_c<br>indole_c + ser__L_c --> trp__L_c + h2o_c<br>3ig3p_c --> indole_c + g3p_c<br>trp__L_c + trnatrp_c + atp_c --> trptrna_c + ppi_c + amp_c<br>trp__L_p + h_p <=> trp__L_c + h_c<br>trp__L_e <=> trp__L_p<br>nadh_c + 2h3oppan_c + h_c <=> nad_c + glyc__R_c<br>tsul_p + atp_c + h2o_c --> h_c + pi_c + adp_c + tsul_c<br>tsul_e <=> tsul_p<br>ttdca_e --> ttdca_p<br>ttdcea_e --> ttdcea_p<br>ttrcyc_e <=> ttrcyc_p<br>h_p + ttrcyc_p --> ttrcyc_e + h_c<br>tungs_p + atp_c + h2o_c --> h_c + tungs_c + pi_c + adp_c<br>tungs_e <=> tungs_p<br>tym_e <=> tym_p<br>amet_c + nadph_c + tyr__L_c --> met__L_c + nadp_c + dhgly_c + 4crsol_c + dad_5_c + h_c<br>tym_p + h2o_p + o2_p --> 4hoxpacd_p + nh4_p + h2o2_p<br>h2o_p + tyrp_p --> tyr__L_p + pi_p<br>tyrp_e <=> tyrp_p<br>akg_c + tyr__L_c <=> 34hpp_c + glu__L_c<br>trnatyr_c + tyr__L_c + atp_c --> tyrtrna_c + ppi_c + amp_c<br>tyr__L_p + h_p <=> tyr__L_c + h_c<br>tyr__L_e <=> tyr__L_p<br>thmpp_c + adp_c + h_c --> pi_c + athtp_c<br>u3hga_c + 3hmrsACP_c --> ACP_c + u23ga_c + h_c<br>uamag_c + atp_c + 26dap__M_c --> ugmd_c + h_c + pi_c + adp_c<br>h2o_p + udpacgal_p --> acgal1p_p + 2.0 h_p + ump_p<br>uacgam_p + h2o_p --> acgam1p_p + 2.0 h_p + ump_p<br>uacgam_e <=> uacgam_p<br>2.0 nad_c + uacmam_c + h2o_c --> 2.0 nadh_c + uacmamu_c + 3.0 h_c<br>uacgam_c <=> uacmam_c<br>3hmrsACP_c + uacgam_c <=> u3aga_c + ACP_c<br>pep_c + uacgam_c --> pi_c + uaccg_c<br>acgam1p_c + utp_c + h_c --> uacgam_c + ppi_c<br>uacgam_c + uagmda_c --> uaagmda_c + udp_c + h_c<br>glu__D_c + uama_c + atp_c --> uamag_c + h_c + pi_c + adp_c<br>ala__L_c + uamr_c + atp_c --> h_c + pi_c + uama_c + adp_c<br>uaccg_c + nadph_c + h_c --> nadp_c + uamr_c<br>udcpdp_c + h2o_c --> udcpp_c + pi_c + h_c<br>frdp_c + 8.0 ipdp_c --> udcpdp_c + 8.0 ppi_c<br>udcpdp_p + h2o_p --> pi_p + udcpp_p + h_p<br>udcpp_p --> udcpp_c<br>udpacgal_e <=> udpacgal_p<br>udpg_c <=> udpgal_c<br>udpgal_c --> udpgalfur_c<br>udpgal_p + h2o_p --> 2.0 h_p + ump_p + gal1p_p<br>udpgal_e <=> udpgal_p<br>2.0 nad_c + udpg_c + h2o_c --> 3.0 h_c + 2.0 nadh_c + udpglcur_c<br>nad_c + udpglcur_c --> udpLa4o_c + nadh_c + co2_c<br>udpglcur_e <=> udpglcur_p<br>udpg_p + h2o_p --> g1p_p + 2.0 h_p + ump_p<br>udpg_e <=> udpg_p<br>udpLa4o_c + glu__L_c <=> udpLa4n_c + akg_c<br>udpglcur_p + h2o_p --> glcur1p_p + 2.0 h_p + ump_p<br>gal1p_c + udpg_c <=> g1p_c + udpgal_c<br>2.0 h_c + urdglyc_c + h2o_c --> co2_c + 2.0 nh4_c + glx_c<br>ugmd_c + alaala_c + atp_c --> h_c + ugmda_c + pi_c + adp_c<br>u3aga_c + h2o_c --> u3hga_c + ac_c<br>10fthf_c + udpLa4n_c --> udpLa4fn_c + thf_c + h_c<br>uLa4n_c --> uLa4n_p<br>uamr_c + LalaDgluMdap_c + atp_c --> ugmd_c + h_c + pi_c + adp_c<br>h2o_c + um4p_c --> ugmd_c + ala__D_c<br>uamr_c + LalaDgluMdapDala_c + atp_c --> um4p_c + h_c + pi_c + adp_c<br>ump_c + atp_c <=> udp_c + adp_c<br>ump_e <=> ump_p<br>uLa4fn_c + h2o_c --> uLa4n_c + for_c<br>udpLa4fn_c + udcpp_c --> udp_c + uLa4fn_c<br>2.0 amet_c + uppg3_c --> 2.0 ahcys_c + h_c + dscl_c<br>hmbil_c --> h2o_c + uppg3_c<br>uppg3_c + 4.0 h_c --> cpppg3_c + 4.0 co2_c<br>ura_c + prpp_c --> ump_c + ppi_c<br>uracp_c + h2o_c --> h_c + cbm_c + poaac_c<br>ura_p + h_p --> ura_c + h_c<br>ura_p + h_p <=> ura_c + h_c<br>ura_e <=> ura_p<br>urdglyc_c + nad_c --> nadh_c + oxur_c + h_c<br>urea_e <=> urea_p<br>urea_p <=> urea_c<br>urate_c + o2_c + 2.0 h2o_c --> alltn_c + co2_c + h2o2_c<br>uri_c + h2o_c --> rib__D_c + ura_c<br>uri_c + gtp_c --> ump_c + gdp_c + h_c<br>h_p + uri_p --> h_c + uri_c<br>h_p + uri_p <=> h_c + uri_c<br>uri_e <=> uri_p<br>u23ga_c + h2o_c --> lipidX_c + ump_c + 2.0 h_c<br>akg_c + val__L_c <=> 3mob_c + glu__L_c<br>trnaval_c + atp_c + val__L_c --> ppi_c + amp_c + valtrna_c<br>val__L_p + h2o_c + atp_c --> h_c + pi_c + adp_c + val__L_c<br>h_p + val__L_p <=> h_c + val__L_c<br>val__L_e <=> val__L_p<br>ala__L_c + 3mob_c <=> pyr_c + val__L_c<br>mptamp_c + 2.0 h_c + tungs_c --> wco_c + amp_c + cu2_c + h2o_c<br>xu5p__L_c --> ru5p__L_c<br>xan_c + nad_c + h2o_c --> nadh_c + h_c + urate_c<br>h_p + xan_p --> xan_c + h_c<br>xan_e <=> xan_p<br>xan_p <=> xan_c<br>xmp_e <=> xmp_p<br>xan_c + prpp_c --> xmp_c + ppi_c<br>xtsn_c + h2o_c --> xan_c + rib__D_c<br>h_p + xtsn_p <=> h_c + xtsn_c<br>xtsn_e <=> xtsn_p<br>xyl__D_c <=> xylu__D_c<br>glc__D_c <=> fru_c<br>xylu__D_c + atp_c --> xu5p__D_c + adp_c + h_c<br>xylu__L_c + atp_c --> xu5p__L_c + adp_c + h_c<br>h_p + xylu__L_p --> xylu__L_c + h_c<br>xylu__L_e <=> xylu__L_p<br>xyl__D_p + atp_c + h2o_c --> xyl__D_c + pi_c + h_c + adp_c<br>xyl__D_p + h_p --> h_c + xyl__D_c<br>xyl__D_e <=> xyl__D_p<br>h2o_c + zn2_c + atp_c --> adp_c + h_c + pi_c + zn2_p<br>h_p + zn2_c --> h_c + zn2_p<br>zn2_p --> zn2_c<br>h2o_c + atp_c + zn2_p --> adp_c + h_c + pi_c + zn2_c<br>zn2_e <=> zn2_p</div></td>
    </tr>
    </table>



Step 2: Simulate a model
------------------------

The model can be simulated by executing
``~cameo.core.solver_based_model.SolverBasedModel.solve``.

.. code:: python

    solution = model.solve()

A quick overview of the solution can be obtained in form of a pandas
`DataFrame <http://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.html>`__
(all solution objects in cameo provide access to data frames through a
``data_frame`` attribute).

.. code:: python

    solution




.. raw:: html

    <div>
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
          <td>0.000219</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>DM_5drib_c</th>
          <td>0.000221</td>
          <td>-3.775930e-24</td>
        </tr>
        <tr>
          <th>DM_aacald_c</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>DM_amob_c</th>
          <td>0.000002</td>
          <td>4.867441e-23</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
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
          <td>0.000000e+00</td>
        </tr>
      </tbody>
    </table>
    <p>2583 rows × 2 columns</p>
    </div>



The data frame is accessible through ``solution.data_frame``.

.. code:: python

    solution.data_frame




.. raw:: html

    <div>
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
          <td>0.000219</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>DM_5drib_c</th>
          <td>0.000221</td>
          <td>-3.775930e-24</td>
        </tr>
        <tr>
          <th>DM_aacald_c</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>DM_amob_c</th>
          <td>0.000002</td>
          <td>4.867441e-23</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
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
          <td>0.000000e+00</td>
        </tr>
      </tbody>
    </table>
    <p>2583 rows × 2 columns</p>
    </div>



Data frames make it very easy to process results. For example, let's
take a look at reactions with flux != 0

.. code:: python

    solution.data_frame.query('fluxes != 0')




.. raw:: html

    <div>
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
          <td>0.000219</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>DM_5drib_c</th>
          <td>0.000221</td>
          <td>-3.775930e-24</td>
        </tr>
        <tr>
          <th>DM_amob_c</th>
          <td>0.000002</td>
          <td>4.867441e-23</td>
        </tr>
        <tr>
          <th>DM_mththf_c</th>
          <td>0.000440</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>USHD</th>
          <td>0.019113</td>
          <td>7.372575e-18</td>
        </tr>
        <tr>
          <th>VALTA</th>
          <td>-0.415702</td>
          <td>-0.000000e+00</td>
        </tr>
        <tr>
          <th>ZN2tpp</th>
          <td>0.000335</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>Zn2tex</th>
          <td>0.000335</td>
          <td>0.000000e+00</td>
        </tr>
      </tbody>
    </table>
    <p>461 rows × 2 columns</p>
    </div>



Step 3: Exploring a model
-------------------------

Objects—models, reactions, metabolites, genes—can easily be explored in
the Jupyter notebook, taking advantage of tab completion. For example,
place your cursor after the period in ``model.reactions.`` and press the
TAB key. A dialog will appear that allows you to navigate the list of
reactions encoded in the model.

.. code:: python

    model.reactions.PGK # delete PGK, place your cursor after the period and press the TAB key.




.. raw:: html

    
            <table>
                <tr>
                    <td><strong>Id</strong></td><td>PGK</td>
                </tr>
                <tr>
                    <td><strong>Name</strong></td><td>Phosphoglycerate kinase</td>
                </tr>
                <tr>
                    <td><strong>Stoichiometry</strong></td><td>3pg_c + atp_c <=> 13dpg_c + adp_c</td>
                </tr>
                <tr>
                    <td><strong>Lower bound</strong></td><td>-1000.000000</td>
                </tr>
                <tr>
                    <td><strong>Upper bound</strong></td><td>1000.000000</td>
                </tr>
            </table>
            



For example, you can access the E4PD (*Erythrose 4-phosphate
dehydrogenase*) reaction in the model.

.. code:: python

    model.reactions.E4PD




.. raw:: html

    
            <table>
                <tr>
                    <td><strong>Id</strong></td><td>E4PD</td>
                </tr>
                <tr>
                    <td><strong>Name</strong></td><td>Erythrose 4-phosphate dehydrogenase</td>
                </tr>
                <tr>
                    <td><strong>Stoichiometry</strong></td><td>e4p_c + nad_c + h2o_c <=> nadh_c + 2.0 h_c + 4per_c</td>
                </tr>
                <tr>
                    <td><strong>Lower bound</strong></td><td>-1000.000000</td>
                </tr>
                <tr>
                    <td><strong>Upper bound</strong></td><td>1000.000000</td>
                </tr>
            </table>
            



Be aware though that due variable naming restrictions in Python dot
notation access to reactions (and other objects) might not work in some
cases.

.. code:: python

    # model.reactions.12DGR120tipp  # uncommenting and running this cell will produce a syntax error

In these cases you need to use the ``model.reactions.get_by_id``.

.. code:: python

    model.reactions.get_by_id('12DGR120tipp')




.. raw:: html

    
            <table>
                <tr>
                    <td><strong>Id</strong></td><td>12DGR120tipp</td>
                </tr>
                <tr>
                    <td><strong>Name</strong></td><td>1,2 diacylglycerol transport via flipping (periplasm to cytoplasm, n-C12:0)</td>
                </tr>
                <tr>
                    <td><strong>Stoichiometry</strong></td><td>12dgr120_p --> 12dgr120_c</td>
                </tr>
                <tr>
                    <td><strong>Lower bound</strong></td><td>0.000000</td>
                </tr>
                <tr>
                    <td><strong>Upper bound</strong></td><td>1000.000000</td>
                </tr>
            </table>
            



Metabolites are accessible through ``model.metabolites``. For example,
D-glucose in the cytosolic compartment.

.. code:: python

    model.metabolites.glc__D_c




.. parsed-literal::

    <Metabolite glc__D_c at 0x110493e80>



And it is easy to find the associated reactions

.. code:: python

    model.metabolites.glc__D_c.reactions




.. parsed-literal::

    frozenset({<Reaction TRE6PH at 0x1105aac18>,
               <Reaction GLCATr at 0x110652048>,
               <Reaction G6PP at 0x110652080>,
               <Reaction GALS3 at 0x110652898>,
               <Reaction MLTG5 at 0x11060d0b8>,
               <Reaction MLTG3 at 0x11060d550>,
               <Reaction LACZ at 0x110623d68>,
               <Reaction MLTG4 at 0x11060d588>,
               <Reaction AMALT1 at 0x10f193dd8>,
               <Reaction TREH at 0x1105aa9e8>,
               <Reaction MLTG1 at 0x11060d5f8>,
               <Reaction MLTG2 at 0x11060d630>,
               <Reaction AMALT2 at 0x10f193e80>,
               <Reaction HEX1 at 0x110635ef0>,
               <Reaction AMALT3 at 0x10f193f28>,
               <Reaction AMALT4 at 0x10f193f60>,
               <Reaction GLCabcpp at 0x110643ba8>,
               <Reaction XYLI2 at 0x11059c7b8>,
               <Reaction GLCt2pp at 0x110643be0>})



A list of the genes encoded in the model can be accessed via
``model.genes``.

.. code:: python

    model.genes[0:10]




.. parsed-literal::

    [<Gene b0241 at 0x11059c6d8>,
     <Gene b1377 at 0x11059c6a0>,
     <Gene b2215 at 0x11059c5f8>,
     <Gene b0929 at 0x11059c630>,
     <Gene b4035 at 0x11059c550>,
     <Gene b4033 at 0x11059c588>,
     <Gene b4034 at 0x11059c4e0>,
     <Gene b4032 at 0x11059c4a8>,
     <Gene b4036 at 0x11059c438>,
     <Gene b4213 at 0x11059c400>]



A few additional attributes have been added that are not available in a
`cobrapy <https://opencobra.github.io/cobrapy/>`__ model. For example,
exchange reactions that allow certain metabolites to enter or leave the
model can be accessed through ``model.exchanges``.

.. code:: python

    model.exchanges[0:10]




.. parsed-literal::

    [<Reaction DM_4crsol_c at 0x109adde48>,
     <Reaction DM_5drib_c at 0x11020e2b0>,
     <Reaction DM_aacald_c at 0x10a959860>,
     <Reaction DM_amob_c at 0x10a959b70>,
     <Reaction DM_mththf_c at 0x10a961780>,
     <Reaction DM_oxam_c at 0x10a961828>,
     <Reaction EX_12ppd__R_e at 0x10fb1d860>,
     <Reaction EX_12ppd__S_e at 0x10fb1d940>,
     <Reaction EX_14glucan_e at 0x10fb1dd30>,
     <Reaction EX_15dap_e at 0x10fb1d978>]



Or, the current medium can be accessed through ``model.medium``.

.. code:: python

    model.medium




.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>reaction_id</th>
          <th>reaction_name</th>
          <th>lower_bound</th>
          <th>upper_bound</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>EX_ca2_e</td>
          <td>Calcium exchange</td>
          <td>-1000.00</td>
          <td>1000.0</td>
        </tr>
        <tr>
          <th>1</th>
          <td>EX_cbl1_e</td>
          <td>Cob(I)alamin exchange</td>
          <td>-0.01</td>
          <td>1000.0</td>
        </tr>
        <tr>
          <th>2</th>
          <td>EX_cl_e</td>
          <td>Chloride exchange</td>
          <td>-1000.00</td>
          <td>1000.0</td>
        </tr>
        <tr>
          <th>3</th>
          <td>EX_co2_e</td>
          <td>CO2 exchange</td>
          <td>-1000.00</td>
          <td>1000.0</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>21</th>
          <td>EX_slnt_e</td>
          <td>Selenite exchange</td>
          <td>-1000.00</td>
          <td>1000.0</td>
        </tr>
        <tr>
          <th>22</th>
          <td>EX_so4_e</td>
          <td>Sulfate exchange</td>
          <td>-1000.00</td>
          <td>1000.0</td>
        </tr>
        <tr>
          <th>23</th>
          <td>EX_tungs_e</td>
          <td>Tungstate exchange</td>
          <td>-1000.00</td>
          <td>1000.0</td>
        </tr>
        <tr>
          <th>24</th>
          <td>EX_zn2_e</td>
          <td>Zinc exchange</td>
          <td>-1000.00</td>
          <td>1000.0</td>
        </tr>
      </tbody>
    </table>
    <p>25 rows × 4 columns</p>
    </div>



It is also possible to get a list of essential reactions ...

.. code:: python

    model.essential_reactions()[0:10]




.. parsed-literal::

    [<Reaction DM_4crsol_c at 0x109adde48>,
     <Reaction DM_5drib_c at 0x11020e2b0>,
     <Reaction DM_amob_c at 0x10a959b70>,
     <Reaction DM_mththf_c at 0x10a961780>,
     <Reaction BIOMASS_Ec_iJO1366_core_53p95M at 0x10fb1dcf8>,
     <Reaction EX_ca2_e at 0x10fb20ef0>,
     <Reaction EX_cl_e at 0x10fb250b8>,
     <Reaction EX_cobalt2_e at 0x10fb251d0>,
     <Reaction EX_cu2_e at 0x10fb25400>,
     <Reaction EX_glc__D_e at 0x10fb323c8>]



... and essential genes.

.. code:: python

    model.essential_genes()[0:10]




.. parsed-literal::

    [<Gene b2311 at 0x11050c0b8>,
     <Gene b0907 at 0x11050c160>,
     <Gene b3608 at 0x110544198>,
     <Gene b0109 at 0x1105141d0>,
     <Gene b0180 at 0x11059c240>,
     <Gene b1281 at 0x11050c2b0>,
     <Gene b0639 at 0x1105142e8>,
     <Gene b0154 at 0x110544470>,
     <Gene b3939 at 0x1104ec4a8>,
     <Gene b3187 at 0x11050c5c0>]


