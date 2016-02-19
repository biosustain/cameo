
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

.. code:: python

    from cameo import load_model



.. raw:: html

    
    
        <script type="text/javascript">
          
          (function(global) {
            function now() {
              return new Date();
            }
          
            if (typeof (window._bokeh_onload_callbacks) === "undefined") {
              window._bokeh_onload_callbacks = [];
            }
          
            function run_callbacks() {
              window._bokeh_onload_callbacks.forEach(function(callback) { callback() });
              delete window._bokeh_onload_callbacks
              console.info("Bokeh: all callbacks have finished");
            }
          
            function load_libs(js_urls, callback) {
              window._bokeh_onload_callbacks.push(callback);
              if (window._bokeh_is_loading > 0) {
                console.log("Bokeh: BokehJS is being loaded, scheduling callback at", now());
                return null;
              }
              if (js_urls == null || js_urls.length === 0) {
                run_callbacks();
                return null;
              }
              console.log("Bokeh: BokehJS not loaded, scheduling load and callback at", now());
              window._bokeh_is_loading = js_urls.length;
              for (var i = 0; i < js_urls.length; i++) {
                var url = js_urls[i];
                var s = document.createElement('script');
                s.src = url;
                s.async = false;
                s.onreadystatechange = s.onload = function() {
                  window._bokeh_is_loading--;
                  if (window._bokeh_is_loading === 0) {
                    console.log("Bokeh: all BokehJS libraries loaded");
                    run_callbacks()
                  }
                };
                s.onerror = function() {
                  console.warn("failed to load library " + url);
                };
                console.log("Bokeh: injecting script tag for BokehJS library: ", url);
                document.getElementsByTagName("head")[0].appendChild(s);
              }
            };var js_urls = ['https://cdn.pydata.org/bokeh/release/bokeh-0.11.0.min.js', 'https://cdn.pydata.org/bokeh/release/bokeh-widgets-0.11.0.min.js', 'https://cdn.pydata.org/bokeh/release/bokeh-compiler-0.11.0.min.js'];
          
            var inline_js = [
              function(Bokeh) {
                Bokeh.set_log_level("info");
              },
              function(Bokeh) {
                console.log("Bokeh: injecting CSS: https://cdn.pydata.org/bokeh/release/bokeh-0.11.0.min.css");
                Bokeh.embed.inject_css("https://cdn.pydata.org/bokeh/release/bokeh-0.11.0.min.css");
                console.log("Bokeh: injecting CSS: https://cdn.pydata.org/bokeh/release/bokeh-widgets-0.11.0.min.css");
                Bokeh.embed.inject_css("https://cdn.pydata.org/bokeh/release/bokeh-widgets-0.11.0.min.css");
              }
            ];
          
            function run_inline_js() {
              for (var i = 0; i < inline_js.length; i++) {
                inline_js[i](window.Bokeh);
              }
            }
          
            if (window._bokeh_is_loading === 0) {
              console.log("Bokeh: BokehJS loaded, going straight to plotting");
              run_inline_js();
            } else {
              load_libs(js_urls, function() {
                console.log("Bokeh: BokehJS plotting callback run at", now());
                run_inline_js();
              });
            }
          }(this));
        </script>


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
    <td><div style="width:100%; max-height:300px; overflow:auto">4crsol_c --> <br>5drib_c --> <br>aacald_c --> <br>amob_c --> <br>mththf_c --> <br>oxam_c --> <br>0.437778 leu__L_c + 0.000335 nadph_c + 0.000223 thmpp_c + 0.004952 ca2_c + 0.012379 nh4_c + 0.28742 arg__L_c + 0.000673 murein4px4px4p_p + 0.001345 murein3p3p_p + 0.180021 phe__L_c + 0.000223 pheme_c + 0.234232 asn__L_c + 0.005381 murein4p4p_p + 3e-06 lipopb_c + 0.333448 lys__L_c + 0.234232 asp__L_c + 0.000307 ni2_c + 0.002944 clpn160_p + 9.8e-05 succoa_c + 0.000279 accoa_c + 0.129799 ctp_c + 0.411184 val__L_c + 0.000223 hemeO_c + 3.1e-05 malcoa_c + 0.000223 ribflv_c + 0.209684 ser__L_c + 0.088988 cys__L_c + 0.025612 dctp_c + 0.000223 sheme_c + 0.031798 pe160_p + 0.000116 bmocogdp_c + 2e-06 btn_c + 0.006388 fe2_c + 0.000223 enter_c + 0.000324 zn2_c + 0.024805 datp_c + 0.004126 so4_c + 0.006744 spmd_c + 0.00118 clpn181_p + 0.000223 gthrd_c + 0.000674 cu2_c + 0.007428 fe3_c + 0.000223 fad_c + 0.008151 colipa_e + 0.255712 gln__L_c + 0.024805 dttp_c + 0.149336 met__L_c + 0.005707 pg160_c + 0.499149 ala__L_c + 0.024732 pe161_p + 0.255712 glu__L_c + 0.140101 utp_c + 0.055234 trp__L_c + 0.000112 nadp_c + 0.214798 pro__L_c + 0.008253 mg2_c + 0.282306 ile__L_c + 0.000223 mlthf_c + 0.03327 ptrc_c + 0.595297 gly_c + 0.000223 adocbl_c + 0.004957 pe181_c + 0.000223 amet_c + 0.092056 his__L_c + 0.000168 coa_c + 0.000605 murein3px4p_p + 0.000658 mn2_c + 0.025612 dgtp_c + 0.000223 5mthf_c + 0.00229 clpn161_p + 7e-06 mobd_c + 0.154187 glycogen_c + 7e-06 mococdp_c + 0.009618 pe161_c + 7e-06 mocogdp_c + 0.012747 pe181_p + 0.000223 pydx5p_c + 0.012366 pe160_c + 0.000223 chor_c + 0.000248 4fe4s_c + 5.5e-05 udcpdp_c + 0.001787 nad_c + 0.004892 pg160_p + 0.000223 mql8_c + 0.003805 pg161_p + 2.5e-05 2fe2s_c + 0.000223 2dmmql8_c + 0.246506 thr__L_c + 0.001961 pg181_p + 54.119975 atp_c + 0.000223 q8h2_c + 0.005448 murein4px4p_p + 0.004439 pg161_c + 48.752916 h2o_c + 0.002288 pg181_c + 0.209121 gtp_c + 0.18569 k_c + 0.004952 cl_c + 2.4e-05 cobalt2_c + 0.133993 tyr__L_c + 0.000223 thf_c + 0.000223 10fthf_c + 4.5e-05 nadh_c --> 53.95 h_c + 53.945874 pi_c + 53.95 adp_c + 0.749831 ppi_c<br>54.124831 atp_c + 0.153686 met__L_c + 0.343161 lys__L_c + 0.450531 leu__L_c + 0.513689 ala__L_c + 0.26316 glu__L_c + 0.000223 ribflv_c + 0.005205 ca2_c + 0.013013 nh4_c + 0.006715 fe2_c + 0.000323 ni2_c + 0.013894 murein5px4p_p + 0.007808 fe3_c + 0.241055 asp__L_c + 0.185265 phe__L_c + 0.000447 nadp_c + 0.290529 ile__L_c + 0.000223 pheme_c + 0.000223 mlthf_c + 0.000223 2ohph_c + 0.612638 gly_c + 0.054154 pe161_c + 0.000223 amet_c + 0.215792 ser__L_c + 0.000576 coa_c + 0.019456 kdo2lipid4_e + 0.027017 dgtp_c + 0.008675 mg2_c + 0.144104 utp_c + 7e-06 mobd_c + 0.423162 val__L_c + 0.045946 pe160_p + 0.026166 dttp_c + 0.02106 pe161_p + 0.000223 thmpp_c + 0.094738 his__L_c + 0.000223 thf_c + 0.000223 pydx5p_c + 0.09158 cys__L_c + 0.027017 dctp_c + 0.056843 trp__L_c + 0.195193 k_c + 0.000223 sheme_c + 5.5e-05 udcpdp_c + 0.133508 ctp_c + 0.000122 bmocogdp_c + 0.295792 arg__L_c + 0.017868 pe160_c + 0.221055 pro__L_c + 0.253687 thr__L_c + 2e-06 btn_c + 0.005205 cl_c + 0.241055 asn__L_c + 0.000341 zn2_c + 0.00026 4fe4s_c + 0.026166 datp_c + 2.6e-05 2fe2s_c + 48.601527 h2o_c + 0.004338 so4_c + 0.215096 gtp_c + 2.5e-05 cobalt2_c + 0.000691 mn2_c + 0.000709 cu2_c + 0.001831 nad_c + 0.000223 fad_c + 0.137896 tyr__L_c + 0.000223 10fthf_c + 0.26316 gln__L_c --> 53.945662 pi_c + 53.95 adp_c + 53.95 h_c + 0.773903 ppi_c<br>12ppd__R_e --> <br>12ppd__S_e --> <br>14glucan_e --> <br>15dap_e --> <br>23camp_e --> <br>23ccmp_e --> <br>23cgmp_e --> <br>23cump_e --> <br>23dappa_e --> <br>26dap__M_e --> <br>2ddglcn_e --> <br>34dhpac_e --> <br>3amp_e --> <br>3cmp_e --> <br>3gmp_e --> <br>3hcinnm_e --> <br>3hpp_e --> <br>3hpppn_e --> <br>3ump_e --> <br>4abut_e --> <br>4hoxpacd_e --> <br>5dglcn_e --> <br>5mtr_e --> <br>LalaDglu_e --> <br>LalaDgluMdap_e --> <br>LalaDgluMdapDala_e --> <br>LalaLglu_e --> <br>ac_e --> <br>acac_e --> <br>acald_e --> <br>acgal_e --> <br>acgal1p_e --> <br>acgam_e --> <br>acgam1p_e --> <br>acmana_e --> <br>acmum_e --> <br>acnam_e --> <br>acolipa_e --> <br>acser_e --> <br>ade_e --> <br>adn_e --> <br>adocbl_e --> <br>ag_e --> <br>agm_e --> <br>akg_e --> <br>ala_B_e --> <br>ala__D_e --> <br>ala__L_e --> <br>alaala_e --> <br>all__D_e --> <br>alltn_e --> <br>amp_e --> <br>anhgm_e --> <br>arab__L_e --> <br>arbt_e --> <br>arbtn_e --> <br>arbtn_fe3_e --> <br>arg__L_e --> <br>ascb__L_e --> <br>asn__L_e --> <br>aso3_e --> <br>asp__L_e --> <br>btn_e --> <br>but_e --> <br>butso3_e --> <br>ca2_e <=> <br>cbi_e --> <br>cbl1_e <=> <br>cd2_e --> <br>cgly_e --> <br>chol_e --> <br>chtbs_e --> <br>cit_e --> <br>cl_e <=> <br>cm_e --> <br>cmp_e --> <br>co2_e <=> <br>cobalt2_e <=> <br>colipa_e --> <br>colipap_e --> <br>cpgn_e --> <br>cpgn_un_e --> <br>crn_e --> <br>crn__D_e --> <br>csn_e --> <br>cu_e --> <br>cu2_e <=> <br>cyan_e --> <br>cynt_e --> <br>cys__D_e --> <br>cys__L_e --> <br>cytd_e --> <br>dad_2_e --> <br>damp_e --> <br>dca_e --> <br>dcmp_e --> <br>dcyt_e --> <br>ddca_e --> <br>dgmp_e --> <br>dgsn_e --> <br>dha_e --> <br>dimp_e --> <br>din_e --> <br>dms_e --> <br>dmso_e --> <br>dopa_e --> <br>doxrbcn_e --> <br>dtmp_e --> <br>dump_e --> <br>duri_e --> <br>eca4colipa_e --> <br>enlipa_e --> <br>enter_e --> <br>etha_e --> <br>ethso3_e --> <br>etoh_e --> <br>f6p_e --> <br>fald_e --> <br>fe2_e <=> <br>fe3_e <=> <br>fe3dcit_e --> <br>fe3dhbzs_e --> <br>fe3hox_e --> <br>fe3hox_un_e --> <br>fecrm_e --> <br>fecrm_un_e --> <br>feenter_e --> <br>feoxam_e --> <br>feoxam_un_e --> <br>for_e --> <br>fru_e --> <br>frulys_e --> <br>fruur_e --> <br>fuc__L_e --> <br>fum_e --> <br>fusa_e --> <br>g1p_e --> <br>g3pc_e --> <br>g3pe_e --> <br>g3pg_e --> <br>g3pi_e --> <br>g3ps_e --> <br>g6p_e --> <br>gal_e --> <br>gal_bD_e --> <br>gal1p_e --> <br>galct__D_e --> <br>galctn__D_e --> <br>galctn__L_e --> <br>galt_e --> <br>galur_e --> <br>gam_e --> <br>gam6p_e --> <br>gbbtn_e --> <br>gdp_e --> <br>glc__D_e <=> <br>glcn_e --> <br>glcr_e --> <br>glcur_e --> <br>glcur1p_e --> <br>gln__L_e --> <br>glu__L_e --> <br>gly_e --> <br>glyald_e --> <br>glyb_e --> <br>glyc_e --> <br>glyc__R_e --> <br>glyc2p_e --> <br>glyc3p_e --> <br>glyclt_e --> <br>gmp_e --> <br>gsn_e --> <br>gthox_e --> <br>gthrd_e --> <br>gtp_e --> <br>gua_e --> <br>h_e <=> <br>h2_e --> <br>h2o_e <=> <br>h2o2_e --> <br>h2s_e --> <br>hacolipa_e --> <br>halipa_e --> <br>hdca_e --> <br>hdcea_e --> <br>hg2_e --> <br>his__L_e --> <br>hom__L_e --> <br>hxa_e --> <br>hxan_e --> <br>idon__L_e --> <br>ile__L_e --> <br>imp_e --> <br>indole_e --> <br>inost_e --> <br>ins_e --> <br>isetac_e --> <br>k_e <=> <br>kdo2lipid4_e --> <br>lac__D_e --> <br>lac__L_e --> <br>lcts_e --> <br>leu__L_e --> <br>lipa_e --> <br>lipa_cold_e --> <br>lipoate_e --> <br>lys__L_e --> <br>lyx__L_e --> <br>mal__D_e --> <br>mal__L_e --> <br>malt_e --> <br>malthx_e --> <br>maltpt_e --> <br>malttr_e --> <br>maltttr_e --> <br>man_e --> <br>man6p_e --> <br>manglyc_e --> <br>melib_e --> <br>meoh_e --> <br>met__D_e --> <br>met__L_e --> <br>metsox_R__L_e --> <br>metsox_S__L_e --> <br>mg2_e <=> <br>mincyc_e --> <br>minohp_e --> <br>mmet_e --> <br>mn2_e <=> <br>mnl_e --> <br>mobd_e <=> <br>mso3_e --> <br>n2o_e --> <br>na1_e <=> <br>nac_e --> <br>nh4_e <=> <br>ni2_e <=> <br>nmn_e --> <br>no_e --> <br>no2_e --> <br>no3_e --> <br>novbcn_e --> <br>o16a4colipa_e --> <br>o2_e <=> <br>o2s_e --> <br>ocdca_e --> <br>ocdcea_e --> <br>octa_e --> <br>orn_e --> <br>orot_e --> <br>pacald_e --> <br>peamn_e --> <br>phe__L_e --> <br>pheme_e --> <br>pi_e <=> <br>pnto__R_e --> <br>ppa_e --> <br>ppal_e --> <br>pppn_e --> <br>ppt_e --> <br>pro__L_e --> <br>progly_e --> <br>psclys_e --> <br>pser__L_e --> <br>ptrc_e --> <br>pydam_e --> <br>pydx_e --> <br>pydxn_e --> <br>pyr_e --> <br>quin_e --> <br>r5p_e --> <br>rfamp_e --> <br>rib__D_e --> <br>rmn_e --> <br>sbt__D_e --> <br>sel_e <=> <br>ser__D_e --> <br>ser__L_e --> <br>skm_e --> <br>slnt_e <=> <br>so2_e --> <br>so3_e --> <br>so4_e <=> <br>spmd_e --> <br>succ_e --> <br>sucr_e --> <br>sulfac_e --> <br>tartr__D_e --> <br>tartr__L_e --> <br>taur_e --> <br>tcynt_e --> <br>thm_e --> <br>thr__L_e --> <br>thrp_e --> <br>thym_e --> <br>thymd_e --> <br>tma_e --> <br>tmao_e --> <br>tre_e --> <br>trp__L_e --> <br>tsul_e --> <br>ttdca_e --> <br>ttdcea_e --> <br>ttrcyc_e --> <br>tungs_e <=> <br>tym_e --> <br>tyr__L_e --> <br>tyrp_e --> <br>uacgam_e --> <br>udpacgal_e --> <br>udpg_e --> <br>udpgal_e --> <br>udpglcur_e --> <br>ump_e --> <br>ura_e --> <br>urea_e --> <br>uri_e --> <br>val__L_e --> <br>xan_e --> <br>xmp_e --> <br>xtsn_e --> <br>xyl__D_e --> <br>xylu__L_e --> <br>zn2_e <=> <br>12dgr120_p --> 12dgr120_c<br>12dgr140_p --> 12dgr140_c<br>12dgr141_p --> 12dgr141_c<br>12dgr160_p --> 12dgr160_c<br>12dgr161_p --> 12dgr161_c<br>12dgr180_p --> 12dgr180_c<br>12dgr181_p --> 12dgr181_c<br>12ppd__R_e <=> 12ppd__R_p<br>12ppd__R_p <=> 12ppd__R_c<br>12ppd__S_e <=> 12ppd__S_p<br>12ppd__S_p <=> 12ppd__S_c<br>atp_c + h2o_c + 14glucan_p --> adp_c + 14glucan_c + pi_c + h_c<br>14glucan_e --> 14glucan_p<br>23camp_e <=> 23camp_p<br>23ccmp_e <=> 23ccmp_p<br>23cgmp_e <=> 23cgmp_p<br>23cump_e <=> 23cump_p<br>h_p + 23dappa_p --> 23dappa_c + h_c<br>23dappa_e <=> 23dappa_p<br>h2o_p + 23cump_p --> h_p + 3ump_p<br>h2o_p + 23ccmp_p --> 3cmp_p + h_p<br>h2o_p + 23camp_p --> h_p + 3amp_p<br>h2o_p + 23cgmp_p --> h_p + 3gmp_p<br>26dap__M_e <=> 26dap__M_p<br>2ddecg3p_p --> 2ddecg3p_c<br>2tdecg3p_p --> 2tdecg3p_c<br>2tdec7eg3p_p --> 2tdec7eg3p_c<br>2hdecg3p_p --> 2hdecg3p_c<br>2hdec9eg3p_p --> 2hdec9eg3p_c<br>2odecg3p_p --> 2odecg3p_c<br>2odec11eg3p_p --> 2odec11eg3p_c<br>2agpe120_p --> 2agpe120_c<br>2agpe140_p --> 2agpe140_c<br>2agpe141_p --> 2agpe141_c<br>2agpe160_p --> 2agpe160_c<br>2agpe161_p --> 2agpe161_c<br>2agpe180_p --> 2agpe180_c<br>2agpe181_p --> 2agpe181_c<br>atp_c + 2agpe120_c + ddca_c --> ppi_c + pe120_c + amp_c<br>atp_c + 2agpe140_c + ttdca_c --> ppi_c + pe140_c + amp_c<br>atp_c + 2agpe141_c + ttdcea_c --> ppi_c + pe141_c + amp_c<br>atp_c + hdca_c + 2agpe160_c --> ppi_c + pe160_c + amp_c<br>atp_c + 2agpe161_c + hdcea_c --> ppi_c + pe161_c + amp_c<br>ocdca_c + atp_c + 2agpe180_c --> ppi_c + pe180_c + amp_c<br>atp_c + 2agpe181_c + ocdcea_c --> amp_c + ppi_c + pe181_c<br>2agpg120_p --> 2agpg120_c<br>2agpg140_p --> 2agpg140_c<br>2agpg141_p --> 2agpg141_c<br>2agpg160_p --> 2agpg160_c<br>2agpg161_p --> 2agpg161_c<br>2agpg180_p --> 2agpg180_c<br>2agpg181_p --> 2agpg181_c<br>atp_c + ddca_c + 2agpg120_c --> pg120_c + ppi_c + amp_c<br>atp_c + ttdca_c + 2agpg140_c --> ppi_c + pg140_c + amp_c<br>atp_c + 2agpg141_c + ttdcea_c --> ppi_c + pg141_c + amp_c<br>atp_c + 2agpg160_c + hdca_c --> amp_c + ppi_c + pg160_c<br>atp_c + hdcea_c + 2agpg161_c --> ppi_c + amp_c + pg161_c<br>ocdca_c + atp_c + 2agpg180_c --> pg180_c + ppi_c + amp_c<br>2agpg181_c + atp_c + ocdcea_c --> ppi_c + pg181_c + amp_c<br>nadh_c + 2dhguln_c + h_c --> nad_c + glcn_c<br>nadph_c + h_c + 2dhguln_c --> glcn_c + nadp_c<br>nadh_c + 2dhguln_c + h_c --> nad_c + idon__L_c<br>nadph_c + h_c + 2dhguln_c --> idon__L_c + nadp_c<br>h2o_c + 2mahmp_c --> 4ampm_c + pi_c + h_c<br>34dhpac_e <=> 34dhpac_p<br>h2o_c + 3amac_c + h_c --> msa_c + nh4_c<br>3amp_e <=> 3amp_p<br>3cmp_e <=> 3cmp_p<br>3gmp_e <=> 3gmp_p<br>3hdecACP_c --> h2o_c + tdec2eACP_c<br>3hddecACP_c --> h2o_c + tddec2eACP_c<br>3hcddec5eACP_c --> t3c5ddeceACP_c + h2o_c<br>3hmrsACP_c --> h2o_c + tmrs2eACP_c<br>3hcmrs7eACP_c --> h2o_c + t3c7mrseACP_c<br>3hpalmACP_c --> tpalm2eACP_c + h2o_c<br>3hcpalm9eACP_c --> t3c9palmeACP_c + h2o_c<br>3hoctaACP_c --> h2o_c + toctd2eACP_c<br>3hcvac11eACP_c --> h2o_c + t3c11vaceACP_c<br>3haACP_c --> h2o_c + but2eACP_c<br>3hhexACP_c --> h2o_c + thex2eACP_c<br>3hoctACP_c --> toct2eACP_c + h2o_c<br>o2_c + 3hcinnm_c + nadh_c + h_c --> nad_c + h2o_c + dhcinnm_c<br>o2_c + 3hpppn_c + h_c + nadh_c --> nad_c + h2o_c + dhpppn_c<br>3hpp_e <=> 3hpp_p<br>3hpp_c + h_c --> h_p + 3hpp_p<br>atp_c + 3dhguln_c --> adp_c + h_c + 3dhgulnp_c<br>h2o_p + 3ump_p --> pi_p + uri_p<br>h2o_p + 3cmp_p --> pi_p + cytd_p<br>h2o_p + 3amp_p --> pi_p + adn_p<br>h2o_p + 3gmp_p --> gsn_p + pi_p<br>nadph_c + 3odecACP_c + h_c <=> nadp_c + 3hdecACP_c<br>nadph_c + h_c + 3oddecACP_c <=> nadp_c + 3hddecACP_c<br>nadph_c + 3ocddec5eACP_c + h_c --> nadp_c + 3hcddec5eACP_c<br>nadph_c + 3omrsACP_c + h_c <=> 3hmrsACP_c + nadp_c<br>nadph_c + h_c + 3ocmrs7eACP_c --> 3hcmrs7eACP_c + nadp_c<br>nadph_c + h_c + 3opalmACP_c <=> 3hpalmACP_c + nadp_c<br>nadph_c + 3ocpalm9eACP_c + h_c --> 3hcpalm9eACP_c + nadp_c<br>nadph_c + h_c + 3ooctdACP_c <=> 3hoctaACP_c + nadp_c<br>nadph_c + 3ocvac11eACP_c + h_c --> 3hcvac11eACP_c + nadp_c<br>nadph_c + actACP_c + h_c <=> 3haACP_c + nadp_c<br>nadph_c + h_c + 3ohexACP_c <=> nadp_c + 3hhexACP_c<br>nadph_c + 3ooctACP_c + h_c <=> nadp_c + 3hoctACP_c<br>ocACP_c + malACP_c + h_c --> 3odecACP_c + ACP_c + co2_c<br>dcaACP_c + malACP_c + h_c --> 3oddecACP_c + ACP_c + co2_c<br>cdec3eACP_c + malACP_c + h_c --> 3ocddec5eACP_c + co2_c + ACP_c<br>ddcaACP_c + malACP_c + h_c --> 3omrsACP_c + ACP_c + co2_c<br>malACP_c + h_c + cddec5eACP_c --> ACP_c + co2_c + 3ocmrs7eACP_c<br>myrsACP_c + malACP_c + h_c --> 3opalmACP_c + ACP_c + co2_c<br>tdeACP_c + malACP_c + h_c --> 3ocpalm9eACP_c + ACP_c + co2_c<br>palmACP_c + malACP_c + h_c --> ACP_c + 3ooctdACP_c + co2_c<br>hdeACP_c + malACP_c + h_c --> 3ocvac11eACP_c + ACP_c + co2_c<br>butACP_c + malACP_c + h_c --> ACP_c + co2_c + 3ohexACP_c<br>hexACP_c + malACP_c + h_c --> ACP_c + 3ooctACP_c + co2_c<br>coa_c + oxadpcoa_c --> succoa_c + accoa_c<br>atp_c + h2o_c + LalaDgluMdap_p --> adp_c + pi_c + LalaDgluMdap_c + h_c<br>LalaDgluMdap_e <=> LalaDgluMdap_p<br>3ump_e <=> 3ump_p<br>dopa_p + h2o_p + o2_p --> 34dhpac_p + h2o2_p + nh4_p<br>4hoxpacd_e <=> 4hoxpacd_p<br>phthr_c + h2o_c --> 4hthr_c + pi_c<br>h2o_c + LalaDgluMdapDala_c --> LalaDgluMdap_c + ala__D_c<br>h2o_p + LalaDgluMdapDala_p --> ala__D_p + LalaDgluMdap_p<br>atp_c + LalaDgluMdapDala_p + h2o_c --> adp_c + pi_c + h_c + LalaDgluMdapDala_c<br>LalaDgluMdapDala_e <=> LalaDgluMdapDala_p<br>nadph_c + 5dglcn_c + h_c <=> glcn_c + nadp_c<br>h_p + 5dglcn_p <=> 5dglcn_c + h_c<br>5dglcn_e <=> 5dglcn_p<br>h2o_c + dad_5_c --> ade_c + 5drib_c<br>5mtr_e <=> 5mtr_p<br>5mtr_c + h_c --> 5mtr_p + h_p<br>ru5p__D_c <=> ara5p_c<br>atp_c + ttdca_c + ACP_c --> ppi_c + myrsACP_c + amp_c<br>atp_c + ACP_c + ttdcea_c --> ppi_c + tdeACP_c + amp_c<br>atp_c + ACP_c + hdca_c --> palmACP_c + ppi_c + amp_c<br>atp_c + hdcea_c + ACP_c --> hdeACP_c + ppi_c + amp_c<br>atp_c + ACP_c + ocdcea_c --> amp_c + ppi_c + octeACP_c<br>ocdca_c + atp_c + ACP_c --> ocdcaACP_c + ppi_c + amp_c<br>atp_c + ddca_c + ACP_c --> ppi_c + ddcaACP_c + amp_c<br>atp_c + dca_c + ACP_c --> ppi_c + dcaACP_c + amp_c<br>atp_c + octa_c + ACP_c --> ppi_c + ocACP_c + amp_c<br>o2_c + h2o_c + aact_c --> nh4_c + h2o2_c + mthgxl_c<br>unagamu_c + dtdp4aaddg_c --> dtdp_c + h_c + unagamuf_c<br>14glucan_c --> malthx_c<br>14glucan_p --> malthx_p<br>arbt6p_c + h2o_c --> g6p_c + hqn_c<br>4abut_c + akg_c --> glu__L_c + sucsal_c<br>4abutn_c + nad_c + h2o_c --> 4abut_c + nadh_c + 2.0 h_c<br>h_p + 4abut_p --> 4abut_c + h_c<br>4abut_e <=> 4abut_p<br>acac_c + accoa_c --> ac_c + aacoa_c<br>2.0 accoa_c <=> coa_c + aacoa_c<br>accoa_c + btcoa_c <=> coa_c + 3ohcoa_c<br>accoa_c + hxcoa_c <=> 3oocoa_c + coa_c<br>occoa_c + accoa_c <=> 3odcoa_c + coa_c<br>accoa_c + dcacoa_c <=> coa_c + 3oddcoa_c<br>accoa_c + ddcacoa_c <=> coa_c + 3otdcoa_c<br>accoa_c + tdcoa_c <=> 3ohdcoa_c + coa_c<br>3oodcoa_c + coa_c <=> pmtcoa_c + accoa_c<br>h_p + acac_p <=> acac_c + h_c<br>acac_e <=> acac_p<br>nad_c + acald_c + coa_c <=> accoa_c + nadh_c + h_c<br>acald_e <=> acald_p<br>acald_p <=> acald_c<br>anth_c + accoa_c --> coa_c + acanth_c<br>gtp_c + adocbip_c + h_c --> ppi_c + agdpcbi_c<br>atp_c + accoa_c + hco3_c --> adp_c + pi_c + malcoa_c + h_c<br>atp_c + coa_c + ppa_c --> adp_c + pi_c + ppcoa_c<br>h2o_p + acgal1p_p --> acgal_p + pi_p<br>acgal1p_e <=> acgal1p_p<br>acgal_e <=> acgal_p<br>h2o_p + acgam1p_p --> pi_p + acgam_p<br>acgam1p_e <=> acgam1p_p<br>atp_c + acgam_c --> adp_c + acgam6p_c + h_c<br>uacgam_c + udcpp_c --> ump_c + unaga_c<br>pep_c + acgam_p --> acgam6p_c + pyr_c<br>acgam_e <=> acgam_p<br>atp_c + acglu_c --> acg5p_c + adp_c<br>glu__L_c + accoa_c --> coa_c + acglu_c + h_c<br>pyr_c + h_c + 2obut_c --> 2ahbut_c + co2_c<br>atp_c + ac_c <=> actp_c + adp_c<br>2.0 pyr_c + h_c --> co2_c + alac__S_c<br>acmum6p_c + h2o_c --> acgam6p_c + lac__D_c<br>unaga_c + uacmamu_c --> udp_c + unagamu_c + h_c<br>pep_c + acmana_p --> pyr_c + acmanap_c<br>acmana_e <=> acmana_p<br>acmum_p + pep_c --> acmum6p_c + pyr_c<br>acmum_e --> acmum_p<br>h_p + acnam_p --> acnam_c + h_c<br>acnam_e <=> acnam_p<br>acnam_c --> acmana_c + pyr_c<br>fad_c + btcoa_c <=> b2coa_c + fadh2_c<br>hxcoa_c + fad_c <=> hx2coa_c + fadh2_c<br>fad_c + occoa_c <=> oc2coa_c + fadh2_c<br>fad_c + dcacoa_c <=> fadh2_c + dc2coa_c<br>fad_c + ddcacoa_c <=> dd2coa_c + fadh2_c<br>tdcoa_c + fad_c <=> td2coa_c + fadh2_c<br>pmtcoa_c + fad_c <=> hdd2coa_c + fadh2_c<br>stcoa_c + fad_c <=> fadh2_c + od2coa_c<br>accoa_c + ACP_c <=> coa_c + acACP_c<br>acorn_c + h2o_c --> ac_c + orn_c<br>acolipa_p + h2o_c + atp_c --> adp_c + pi_c + h_c + acolipa_e<br>acon_T_c <=> acon_C_c<br>amet_c + acon_T_c --> ahcys_c + aconm_c<br>cit_c <=> h2o_c + acon_C_c<br>h2o_c + acon_C_c <=> icit_c<br>acorn_c + akg_c <=> glu__L_c + acg5sa_c<br>h_c + ddcaACP_c + pi_c --> ddcap_c + ACP_c<br>myrsACP_c + pi_c + h_c --> ACP_c + ttdcap_c<br>tdeACP_c + h_c + pi_c --> ttdceap_c + ACP_c<br>palmACP_c + pi_c + h_c --> hdcap_c + ACP_c<br>hdeACP_c + h_c + pi_c --> hdceap_c + ACP_c<br>ocdcaACP_c + h_c + pi_c --> ocdcap_c + ACP_c<br>octeACP_c + pi_c + h_c --> ACP_c + ocdceap_c<br>coa_c + apoACP_c --> h_c + ACP_c + pap_c<br>atp_c + ac_c + coa_c --> ppi_c + amp_c + accoa_c<br>acser_e <=> acser_p<br>acser_c --> acser_p<br>ac_p + h_p <=> h_c + ac_c<br>ac_p + na1_p --> ac_c + na1_c<br>ac_e <=> ac_p<br>adn_c + h_c + h2o_c --> ins_c + nh4_c<br>4adcho_c --> 4abz_c + pyr_c + h_c<br>chor_c + gln__L_c --> glu__L_c + 4adcho_c<br>ade_c + h2o_c + h_c --> nh4_c + hxan_c<br>h_p + ade_p <=> ade_c + h_c<br>ade_e <=> ade_p<br>atp_c + amp_c <=> 2.0 adp_c<br>gtp_c + amp_c <=> gdp_c + adp_c<br>amp_c + itp_c <=> adp_c + idp_c<br>amet_c + h_c --> ametam_c + co2_c<br>atp_c --> ppi_c + camp_c<br>atp_c + adn_c --> adp_c + amp_c + h_c<br>adn_c + h2o_c --> ade_c + rib__D_c<br>h_p + adn_p --> adn_c + h_c<br>h_p + adn_p <=> adn_c + h_c<br>adn_e <=> adn_p<br>atp_c + adocbi_c --> adp_c + adocbip_c + h_c<br>rdmbzi_c + agdpcbi_c --> gmp_c + adocbl_c + h_c<br>atp_c + adocbl_p + h2o_c --> adp_c + pi_c + h_c + adocbl_c<br>h_p + adocbl_e --> adocbl_p + h_c<br>adprib_c + h2o_c --> amp_c + r5p_c + 2.0 h_c<br>ade_c + prpp_c --> ppi_c + amp_c<br>atp_c + aps_c --> adp_c + paps_c + h_c<br>dcamp_c <=> amp_c + fum_c<br>25aics_c <=> fum_c + aicar_c<br>gtp_c + imp_c + asp__L_c --> gdp_c + pi_c + dcamp_c + 2.0 h_c<br>acgam6p_c + h2o_c --> ac_c + gam6p_c<br>anhgm3p_c + h2o_c --> LalaDgluMdap_c + anhgm_c<br>h2o_p + anhgm3p_p --> anhgm_p + LalaDgluMdap_p<br>anhgm3p_c + h2o_c --> acgam_c + anhm3p_c<br>anhgm3p_p + h_p --> anhgm3p_c + h_c<br>h2o_c + anhgm4p_c --> LalaDgluMdapDala_c + anhgm_c<br>h2o_p + anhgm4p_p --> LalaDgluMdapDala_p + anhgm_p<br>h2o_c + anhgm4p_c --> ala__D_c + anhgm3p_c<br>h2o_p + anhgm4p_p --> anhgm3p_p + ala__D_p<br>h2o_c + anhgm4p_c --> anhm4p_c + acgam_c<br>h_p + anhgm4p_p --> h_c + anhgm4p_c<br>h2o_c + anhgm_c --> anhm_c + acgam_c<br>adphep_DD_c --> adphep_LD_c<br>h2o_c + agm_c --> urea_c + ptrc_c<br>h_p + anhgm_p --> anhgm_c + h_c<br>agm_e <=> agm_p<br>1ddecg3p_c + ddcaACP_c --> pa120_c + ACP_c<br>myrsACP_c + 1tdecg3p_c --> pa140_c + ACP_c<br>tdeACP_c + 1tdec7eg3p_c --> ACP_c + pa141_c<br>palmACP_c + 1hdecg3p_c --> ACP_c + pa160_c<br>hdeACP_c + 1hdec9eg3p_c --> pa161_c + ACP_c<br>1odecg3p_c + ocdcaACP_c --> ACP_c + pa180_c<br>1odec11eg3p_c + octeACP_c --> ACP_c + pa181_c<br>nadp_c + pi_c + acg5sa_c <=> nadph_c + acg5p_c + h_c<br>ag_c + h_e --> ag_e + h_c<br>h2o_c + ahcys_c --> ade_c + rhcys_c<br>aicar_c + 10fthf_c <=> thf_c + fprica_c<br>atp_c + hco3_c + air_c --> adp_c + pi_c + 5caiz_c + h_c<br>5aizc_c <=> 5caiz_c<br>coa_c + nad_c + akg_c --> succoa_c + nadh_c + co2_c<br>h_p + akg_p <=> akg_c + h_c<br>akg_e <=> akg_p<br>alaala_c + h2o_c --> 2.0 ala__D_c<br>atp_c + alaala_p + h2o_c --> adp_c + pi_c + alaala_c + h_c<br>atp_c + 2.0 ala__D_c <=> adp_c + pi_c + alaala_c + h_c<br>alaala_e <=> alaala_p<br>LalaDglu_c <=> LalaLglu_c<br>ala__L_c <=> ala__D_c<br>ala__D_c + pydx5p_c --> pyam5p_c + pyr_c<br>akg_c + ala__L_c <=> glu__L_c + pyr_c<br>pydx5p_c + ala__L_c --> pyam5p_c + pyr_c<br>atp_c + trnaala_c + ala__L_c --> alatrna_c + ppi_c + amp_c<br>ala__L_p + atp_c + h2o_c --> pi_c + ala__L_c + adp_c + h_c<br>ala__L_p + h_p --> h_c + ala__L_c<br>ala__L_p + h_p <=> h_c + ala__L_c<br>ala__L_p + na1_p --> ala__L_c + na1_c<br>ala__L_e <=> ala__L_p<br>h_c + glyald_c + nadh_c <=> nad_c + glyc_c<br>nad_c + etoh_c <=> acald_c + h_c + nadh_c<br>nad_c + h2o_c + pacald_c <=> pac_c + nadh_c + 2.0 h_c<br>acald_c + nad_c + h2o_c --> ac_c + nadh_c + 2.0 h_c<br>acald_c + nadp_c + h2o_c --> ac_c + nadph_c + 2.0 h_c<br>nadp_c + h2o_c + ppal_c --> nadph_c + ppa_c + 2.0 h_c<br>nad_c + h2o_c + btal_c --> but_c + nadh_c + 2.0 h_c<br>atp_c + all__D_c --> adp_c + h_c + all6p_c<br>all6p_c <=> allul6p_c<br>alltt_c + 2.0 h2o_c + 2.0 h_c --> co2_c + urdglyc_c + 2.0 nh4_c<br>h2o_c + alltn_c --> h_c + alltt_c<br>h_p + alltn_p <=> alltn_c + h_c<br>alltn_e <=> alltn_p<br>allul6p_c <=> f6p_c<br>atp_c + all__D_p + h2o_c --> adp_c + all__D_c + pi_c + h_c<br>all__D_e <=> all__D_p<br>alpp_p + pe160_p --> lpp_p + 2agpe160_p<br>pg160_p + alpp_p --> lpp_p + 2agpg160_p<br>nadph_c + mthgxl_c + h_c --> nadp_c + acetol_c<br>h_c + acetol_c + nadh_c --> nad_c + 12ppd__R_c<br>altrn_c --> h2o_c + 2ddglcn_c<br>h2o_c + anhm3p_c --> LalaDgluMdap_c + anhm_c<br>anhm4p_c + h2o_c --> anhm_c + LalaDgluMdapDala_c<br>anhm4p_c + h2o_c --> ala__D_c + anhm3p_c<br>malttr_c + malt_c --> glc__D_c + maltttr_c<br>maltttr_c + malt_c --> glc__D_c + maltpt_c<br>maltpt_c + malt_c --> malthx_c + glc__D_c<br>malthx_c + malt_c --> glc__D_c + malthp_c<br>acmanap_c <=> acgam6p_c<br>atp_c + acmana_c --> adp_c + h_c + acmanap_c<br>amet_c + 8aonn_c <=> amob_c + dann_c<br>2dmmql8_c + amet_c --> mql8_c + h_c + ahcys_c<br>air_c + nad_c + h2o_c --> 2.0 for_c + 4ampm_c + nadh_c + 3.0 h_c<br>amp_c + h2o_c --> ade_c + r5p_c<br>h2o_c + cgly_c --> gly_c + cys__L_c<br>h2o_c + progly_c --> pro__L_c + gly_c<br>amp_e <=> amp_p<br>anhgm_e <=> anhgm_p<br>anhm_c + h2o_c + atp_c --> acmum6p_c + adp_c + h_c<br>anth_c + prpp_c --> pran_c + ppi_c<br>chor_c + gln__L_c --> pyr_c + glu__L_c + anth_c + h_c<br>2aobut_c + h_c --> co2_c + aact_c<br>pimACP_c + ala__L_c --> ACP_c + 8aonn_c + co2_c<br>h2o_c + ap4a_c --> 2.0 adp_c + 2.0 h_c<br>2.0 atp_c + h_c --> ppi_c + ap4a_c<br>ap5a_c + h2o_c --> atp_c + adp_c + 2.0 h_c<br>ametam_c + 15dap_c --> na15dap_c + 5mta_c + h_c<br>ddcap_c + glyc3p_c --> h_c + pi_c + 1ddecg3p_c<br>glyc3p_c + ttdcap_c --> 1tdecg3p_c + h_c + pi_c<br>ttdceap_c + glyc3p_c --> 1tdec7eg3p_c + pi_c + h_c<br>hdcap_c + glyc3p_c --> 1hdecg3p_c + h_c + pi_c<br>hdceap_c + glyc3p_c --> pi_c + 1hdec9eg3p_c + h_c<br>ocdcap_c + glyc3p_c --> 1odecg3p_c + pi_c + h_c<br>glyc3p_c + ocdceap_c --> 1odec11eg3p_c + h_c + pi_c<br>ddcap_c + h2o_c --> ddca_c + 2.0 h_c + pi_c<br>h2o_c + ttdcap_c --> ttdca_c + 2.0 h_c + pi_c<br>ttdceap_c + h2o_c --> pi_c + ttdcea_c + 2.0 h_c<br>hdcap_c + h2o_c --> hdca_c + 2.0 h_c + pi_c<br>hdceap_c + h2o_c --> hdcea_c + 2.0 h_c + pi_c<br>h2o_c + ocdcap_c --> ocdca_c + 2.0 h_c + pi_c<br>h2o_c + ocdceap_c --> ocdcea_c + pi_c + 2.0 h_c<br>aact_c + h_c + nadh_c <=> nad_c + appl_c<br>nadph_c + h_c + 5apru_c --> 5aprbu_c + nadp_c<br>arab__L_c <=> rbl__L_c<br>2.0 arbtn_fe3_c + fadh2_c --> 2.0 fe2_c + 2.0 arbtn_c + fad_c + 2.0 h_c<br>fmnh2_c + 2.0 arbtn_fe3_c --> 2.0 fe2_c + fmn_c + 2.0 arbtn_c + 2.0 h_c<br>rbflvrd_c + 2.0 arbtn_fe3_c --> 2.0 fe2_c + ribflv_c + 2.0 arbtn_c + 2.0 h_c<br>atp_c + arbtn_fe3_p + h2o_c --> adp_c + pi_c + h_c + arbtn_fe3_c<br>arbtn_e + fe3_e --> arbtn_fe3_e<br>h_p + arbtn_p --> arbtn_e + h_c<br>h_p + arbtn_fe3_e --> arbtn_fe3_p + h_c<br>h_p + arbtn_c --> h_c + arbtn_p<br>arbt_p + pep_c --> arbt6p_c + pyr_c<br>arbt_e --> arbt_p<br>atp_c + arab__L_p + h2o_c --> adp_c + pi_c + h_c + arab__L_c<br>h_p + arab__L_p <=> h_c + arab__L_c<br>h_p + arab__L_c --> arab__L_p + h_c<br>arab__L_e <=> arab__L_p<br>agm_c + arg__L_p <=> arg__L_c + agm_p<br>arg__L_c + h_c --> co2_c + agm_c<br>h_p + arg__L_p --> co2_p + agm_p<br>orn_c + arg__L_p <=> orn_p + arg__L_c<br>argsuc_c <=> arg__L_c + fum_c<br>atp_c + citr__L_c + asp__L_c --> ppi_c + argsuc_c + amp_c + h_c<br>atp_c + arg__L_c + trnaarg_c --> ppi_c + argtrna_c + amp_c<br>atp_c + h2o_c + arg__L_p --> adp_c + arg__L_c + pi_c + h_c<br>h_p + arg__L_c --> h_c + arg__L_p<br>arg__L_e <=> arg__L_p<br>aspsa_c + nadp_c + pi_c <=> 4pasp_c + nadph_c + h_c<br>h2o_c + ascb6p_c --> 3dhgulnp_c + h_c<br>ascb__L_p + pep_c --> pyr_c + ascb6p_c<br>ascb__L_e <=> ascb__L_p<br>asn__L_c + h2o_c --> asp__L_c + nh4_c<br>h2o_p + asn__L_p --> asp__L_p + nh4_p<br>atp_c + h2o_c + asp__L_c + gln__L_c --> asn__L_c + ppi_c + glu__L_c + amp_c + h_c<br>atp_c + asp__L_c + nh4_c --> asn__L_c + amp_c + ppi_c + h_c<br>asn__L_c + atp_c + trnaasn_c --> ppi_c + asntrna_c + amp_c<br>asn__L_p + h2o_c + atp_c --> asn__L_c + adp_c + pi_c + h_c<br>h_p + asn__L_p <=> asn__L_c + h_c<br>asn__L_e <=> asn__L_p<br>atp_c + aso3_c + h2o_c --> adp_c + pi_c + aso3_p + h_c<br>aso3_e <=> aso3_p<br>h_c + asp__L_c --> co2_c + ala_B_c<br>cbp_c + asp__L_c --> cbasp_c + h_c + pi_c<br>atp_c + asp__L_c <=> adp_c + 4pasp_c<br>asp__L_c + q8_c --> q8h2_c + iasp_c + h_c<br>mqn8_c + asp__L_c --> iasp_c + h_c + mql8_c<br>fum_c + asp__L_c --> h_c + iasp_c + succ_c<br>o2_c + asp__L_c --> iasp_c + h_c + h2o2_c<br>asp__L_c --> fum_c + nh4_c<br>akg_c + asp__L_c <=> glu__L_c + oaa_c<br>atp_c + trnaasp_c + asp__L_c --> amp_c + asptrna_c + ppi_c<br>atp_c + h2o_c + asp__L_p --> adp_c + pi_c + h_c + asp__L_c<br>2.0 h_p + asp__L_p --> 2.0 h_c + asp__L_c<br>3.0 h_p + asp__L_p --> 3.0 h_c + asp__L_c<br>h_p + asp__L_p --> h_c + asp__L_c<br>h_p + asp__L_p <=> h_c + asp__L_c<br>asp__L_e <=> asp__L_p<br>aso4_c + 2.0 gthrd_c --> h2o_c + gthox_c + aso3_c<br>succoa_c + arg__L_c --> sucarg_c + coa_c + h_c<br>athr__L_c + nadp_c <=> nadph_c + 2aobut_c + h_c<br>atp_c + h2o_c + h_c --> nh4_c + itp_c<br>atp_c + h2o_c --> adp_c + pi_c + h_c<br>atp_c + prpp_c --> ppi_c + prbatp_c<br>adp_c + pi_c + 4.0 h_p <=> atp_c + h2o_c + 3.0 h_c<br>h_p + ala_B_p --> ala_B_c + h_c<br>ala_B_e <=> ala_B_p<br>nad_c + betald_c + h2o_c --> glyb_c + nadh_c + 2.0 h_c<br>nadp_c + betald_c + h2o_c --> nadph_c + glyb_c + 2.0 h_c<br>mptamp_c + moco_c --> cu2_c + bmoco_c + amp_c<br>bmoco_c + gtp_c + h_c --> ppi_c + bmoco1gdp_c<br>gtp_c + h_c + bmoco1gdp_c --> bmocogdp_c + ppi_c<br>pap_c + h2o_c --> amp_c + pi_c<br>btnso_c + nadh_c + h_c --> btn_c + nad_c + h2o_c<br>btnso_c + nadph_c + h_c --> btn_c + nadp_c + h2o_c<br>h_p + btn_p --> btn_c + h_c<br>btn_e <=> btn_p<br>dtbt_c + amet_c + 2fe2s_c --> met__L_c + btn_c + dad_5_c + 2fe1s_c + h_c<br>but_c + accoa_c --> ac_c + btcoa_c<br>atp_c + butso3_p + h2o_c --> adp_c + butso3_c + pi_c + h_c<br>butso3_e <=> butso3_p<br>h_p + but_p <=> but_c + h_c<br>but_e <=> but_p<br>bwco_c + gtp_c + h_c --> ppi_c + bwco1gdp_c<br>gtp_c + bwco1gdp_c + h_c --> ppi_c + bwcogdp_c<br>wco_c + mptamp_c --> cu2_c + bwco_c + amp_c<br>ca2_c + h_p --> ca2_p + h_c<br>ca2_e <=> ca2_p<br>h_p + lys__L_p + 15dap_c --> 15dap_p + lys__L_c + h_c<br>2.0 h2o2_c --> o2_c + 2.0 h2o_c<br>ca2_c + na1_p <=> ca2_p + na1_c<br>atp_c + cbi_c + h_c <=> adocbi_c + pppi_c<br>cbi_e + h_p --> cbi_p + h_c<br>atp_c + h2o_c + cbi_p --> adp_c + pi_c + cbi_c + h_c<br>cbl1_p + h2o_c + atp_c --> cbl1_c + pi_c + adp_c + h_c<br>h_p + cbl1_e --> cbl1_p + h_c<br>atp_c + cbl1_c + h_c <=> adocbl_c + pppi_c<br>cbm_c + 2.0 h_c --> co2_c + nh4_c<br>atp_c + co2_c + nh4_c <=> adp_c + cbp_c + 2.0 h_c<br>2.0 atp_c + hco3_c + h2o_c + gln__L_c --> 2.0 adp_c + pi_c + cbp_c + glu__L_c + 2.0 h_c<br>atp_c + cdg_c + nh4_c --> adp_c + preq0_c + pi_c + h2o_c + h_c<br>atp_c + cd2_c + h2o_c --> adp_c + cd2_p + pi_c + h_c<br>h_p + cd2_c --> cd2_p + h_c<br>cd2_e <=> cd2_p<br>cd2_p --> cd2_c<br>h2o_c + cdpdddecg_c --> pa120_c + 2.0 h_c + cmp_c<br>cdpdtdecg_c + h2o_c --> pa140_c + 2.0 h_c + cmp_c<br>h2o_c + cdpdtdec7eg_c --> pa141_c + 2.0 h_c + cmp_c<br>h2o_c + cdpdhdecg_c --> cmp_c + pa160_c + 2.0 h_c<br>cdpdhdec9eg_c + h2o_c --> cmp_c + 2.0 h_c + pa161_c<br>cdpdodecg_c + h2o_c --> 2.0 h_c + pa180_c + cmp_c<br>h2o_c + cdpdodec11eg_c --> pa181_c + 2.0 h_c + cmp_c<br>2.0 nadph_c + preq0_c + 3.0 h_c --> 2.0 nadp_c + preq1_c<br>h_c + cph4_c --> nh4_c + cdg_c<br>atp_c + 4c2me_c --> 2p4c2me_c + h_c + adp_c<br>2.0 amet_c + pe161_c --> 2.0 ahcys_c + 2.0 h_c + cpe160_c<br>2.0 amet_c + pg161_c --> cpg160_c + 2.0 h_c + 2.0 ahcys_c<br>2.0 amet_c + pe181_c --> cpe180_c + 2.0 ahcys_c + 2.0 h_c<br>2.0 amet_c + pg181_c --> cpg180_c + 2.0 ahcys_c + 2.0 h_c<br>atp_c + cgly_p + h2o_c --> adp_c + pi_c + cgly_c + h_c<br>cgly_e <=> cgly_p<br>atp_c + chol_p + h2o_c --> adp_c + pi_c + h_c + chol_c<br>h_p + chol_p --> h_c + chol_c<br>chol_e <=> chol_p<br>nad_c + chol_c --> betald_c + h_c + nadh_c<br>chor_c --> pphn_c<br>3psme_c --> pi_c + chor_c<br>chor_c --> 4hbz_c + pyr_c<br>chtbs_p + pep_c --> pyr_c + chtbs6p_c<br>chtbs_e <=> chtbs_p<br>cinnm_c + o2_c + nadh_c + h_c --> cenchddd_c + nad_c<br>cit_c --> ac_c + oaa_c<br>h_p + cit_c --> cit_p + h_c<br>cit_p + succ_c --> succ_p + cit_c<br>cit_e <=> cit_p<br>atp_c + h2o_c + lipa_cold_p --> adp_c + lipa_cold_e + pi_c + h_c<br>h2o_p + clpn120_p --> h_p + pg120_p + pa120_p<br>h2o_p + clpn140_p --> h_p + pg140_p + pa140_p<br>h2o_p + clpn141_p --> h_p + pa141_p + pg141_p<br>h2o_p + clpn160_p --> h_p + pa160_p + pg160_p<br>h2o_p + clpn161_p --> h_p + pg161_p + pa161_p<br>h2o_p + clpn180_p --> h_p + pg180_p + pa180_p<br>h2o_p + clpn181_p --> h_p + pa181_p + pg181_p<br>2.0 pg120_p <=> clpn120_p + glyc_p<br>2.0 pg140_p <=> glyc_p + clpn140_p<br>2.0 pg141_p <=> clpn141_p + glyc_p<br>2.0 pg160_p <=> glyc_p + clpn160_p<br>2.0 pg161_p <=> glyc_p + clpn161_p<br>2.0 pg180_p <=> clpn180_p + glyc_p<br>2.0 pg181_p <=> glyc_p + clpn181_p<br>2.0 cl_p + h_c --> 2.0 cl_c + h_p<br>cl_e <=> cl_p<br>h2o_c + cmp_c --> r5p_c + csn_c<br>cmp_e <=> cmp_p<br>cm_e <=> cm_p<br>cm_p + h_p --> cm_e + h_c<br>co2_e <=> co2_p<br>co2_p <=> co2_c<br>atp_c + cobalt2_c + h2o_c --> adp_c + pi_c + cobalt2_p + h_c<br>h_p + cobalt2_c --> cobalt2_p + h_c<br>cobalt2_e <=> cobalt2_p<br>cobalt2_p --> cobalt2_c<br>udcpdp_p + colipa_p --> colipap_p + udcpp_p<br>atp_c + colipap_p + h2o_c --> adp_c + pi_c + colipap_e + h_c<br>atp_c + colipa_c + h2o_c --> adp_c + pi_c + h_c + colipa_p<br>atp_c + h2o_c + colipa_p --> adp_c + colipa_e + pi_c + h_c<br>2.0 cpgn_c + fadh2_c --> 2.0 fe2_c + 2.0 cpgn_un_c + fad_c + 2.0 h_c<br>fmnh2_c + 2.0 cpgn_c --> 2.0 fe2_c + fmn_c + 2.0 cpgn_un_c + 2.0 h_c<br>rbflvrd_c + 2.0 cpgn_c --> 2.0 fe2_c + ribflv_c + 2.0 cpgn_un_c + 2.0 h_c<br>h_p + cpgn_un_p --> h_c + cpgn_un_e<br>h_p + cpgn_un_c --> cpgn_un_p + h_c<br>atp_c + cpgn_p + h2o_c --> adp_c + pi_c + h_c + cpgn_c<br>cpgn_un_e + fe3_e --> cpgn_e<br>h_p + cpgn_e --> h_c + cpgn_p<br>ahdt_c + h2o_c --> acald_c + cph4_c + pppi_c + h_c<br>gtp_c + h2o_c --> ppi_c + cpmp_c<br>o2_c + 2.0 h_c + cpppg3_c --> 2.0 co2_c + pppg9_c + 2.0 h2o_c<br>2.0 amet_c + cpppg3_c --> 2.0 met__L_c + pppg9_c + 2.0 co2_c + 2.0 dad_5_c<br>crn_c + bbtcoa_c <=> gbbtn_c + crncoa_c<br>crn_c + coa_c + atp_c --> adp_c + pi_c + crncoa_c<br>crncoa_c <=> crnDcoa_c<br>crn_c + ctbtcoa_c <=> ctbt_c + crncoa_c<br>crncoa_c <=> h2o_c + ctbtcoa_c<br>atp_c + coa_c + crn__D_c --> adp_c + pi_c + crnDcoa_c<br>atp_c + crn__D_p + h2o_c --> adp_c + crn__D_c + pi_c + h_c<br>crn__D_p + h_p <=> crn__D_c + h_c<br>crn__D_e <=> crn__D_p<br>h2o_c + crn_p + atp_c --> crn_c + adp_c + pi_c + h_c<br>h_p + crn_p <=> crn_c + h_c<br>gbbtn_c + crn_p --> crn_c + gbbtn_p<br>crn__D_c + crn_p --> crn_c + crn__D_p<br>crn_e <=> crn_p<br>accoa_c + oaa_c + h2o_c --> cit_c + coa_c + h_c<br>h2o_c + h_c + csn_c --> ura_c + nh4_c<br>h_p + csn_p --> h_c + csn_c<br>csn_e <=> csn_p<br>atp_c + coa_c + ctbt_c --> adp_c + pi_c + ctbtcoa_c<br>atp_c + ctbt_p + h2o_c --> adp_c + ctbt_c + pi_c + h_c<br>ctbt_p + h_p <=> ctbt_c + h_c<br>tdecoa_c <=> td2coa_c<br>hdcoa_c <=> hdd2coa_c<br>odecoa_c <=> od2coa_c<br>atp_c + utp_c + h2o_c + gln__L_c --> adp_c + glu__L_c + pi_c + ctp_c + 2.0 h_c<br>4.0 h_p + 4.0 cu_p + o2_p --> 2.0 h2o_p + 4.0 cu2_p<br>atp_c + h2o_c + cu_c --> adp_c + cu_p + pi_c + h_c<br>cu2_c + h2o_c + atp_c --> adp_c + pi_c + h_c + cu2_p<br>cu2_e <=> cu2_p<br>cu2_p --> cu2_c<br>cu_c + h_e --> cu_e + h_c<br>cu_e <=> cu_p<br>cyan_c + tsul_c --> tcynt_c + so3_c + h_c<br>tsul_p + cyan_p --> h_p + so3_p + tcynt_p<br>cyan_e <=> cyan_p<br>3.0 h_c + hco3_c + cynt_c --> 2.0 co2_c + nh4_c<br>h_p + cynt_p --> h_c + cynt_c<br>cynt_e <=> cynt_p<br>cys__D_c + h2o_c --> h2s_c + pyr_c + nh4_c<br>h2o_c + cys__L_c --> h2s_c + pyr_c + nh4_c<br>atp_c + cys__D_p + h2o_c --> adp_c + pi_c + cys__D_c + h_c<br>cys__D_e <=> cys__D_p<br>h2s_c + acser_c --> ac_c + h_c + cys__L_c<br>3sala_c + 2.0 h_c --> so2_c + ala__L_c<br>cyst__L_c + h2o_c --> nh4_c + pyr_c + hcys__L_c<br>atp_c + trnacys_c + cys__L_c --> amp_c + cystrna_c + ppi_c<br>atp_c + h2o_c + cys__L_c --> adp_c + pi_c + cys__L_p + h_c<br>atp_c + h2o_c + cys__L_p --> adp_c + pi_c + cys__L_c + h_c<br>cys__L_e <=> cys__L_p<br>cys__L_c --> cys__L_p<br>0.5 o2_c + 2.0 h_c + mql8_c --> 2.0 h_p + h2o_c + mqn8_c<br>0.5 o2_c + q8h2_c + 2.0 h_c --> 2.0 h_p + q8_c + h2o_c<br>0.5 o2_c + q8h2_c + 4.0 h_c --> 4.0 h_p + q8_c + h2o_c<br>cytd_c + h2o_c + h_c --> nh4_c + uri_c<br>cytd_c + h2o_c --> rib__D_c + csn_c<br>cytd_c + gtp_c --> gdp_c + cmp_c + h_c<br>h_p + cytd_p --> cytd_c + h_c<br>h_p + cytd_p <=> cytd_c + h_c<br>cytd_e <=> cytd_p<br>atp_c + cmp_c <=> cdp_c + adp_c<br>atp_c + dcmp_c <=> dcdp_c + adp_c<br>h_p + lac__D_p <=> h_c + lac__D_c<br>lac__D_e <=> lac__D_p<br>ala__D_c + h2o_c + fad_c --> fadh2_c + pyr_c + nh4_c<br>h2o_c + dad_2_c + h_c --> din_c + nh4_c<br>atp_c + damp_c <=> adp_c + dadp_c<br>dad_2_p + h_p --> dad_2_c + h_c<br>dad_2_e <=> dad_2_p<br>atp_c + 12dgr120_c --> adp_c + pa120_c + h_c<br>atp_c + 12dgr140_c --> pa140_c + h_c + adp_c<br>atp_c + 12dgr141_c --> adp_c + h_c + pa141_c<br>atp_c + 12dgr160_c --> adp_c + pa160_c + h_c<br>atp_c + 12dgr161_c --> adp_c + h_c + pa161_c<br>atp_c + 12dgr180_c --> adp_c + pa180_c + h_c<br>atp_c + 12dgr181_c --> adp_c + h_c + pa181_c<br>h_p + ala__D_p --> ala__D_c + h_c<br>ala__D_e <=> ala__D_p<br>damp_e <=> damp_p<br>h2o_c + 23dappa_c --> pyr_c + 2.0 nh4_c<br>h_c + 26dap__M_c --> co2_c + lys__L_c<br>26dap_LL_c <=> 26dap__M_c<br>atp_c + 26dap__M_p + h2o_c --> adp_c + pi_c + h_c + 26dap__M_c<br>15dap_e <=> 15dap_p<br>pa120_c + ctp_c + h_c --> ppi_c + cdpdddecg_c<br>pa140_c + h_c + ctp_c --> cdpdtdecg_c + ppi_c<br>ctp_c + h_c + pa141_c --> ppi_c + cdpdtdec7eg_c<br>ctp_c + pa160_c + h_c --> ppi_c + cdpdhdecg_c<br>ctp_c + h_c + pa161_c --> cdpdhdec9eg_c + ppi_c<br>ctp_c + pa180_c + h_c --> cdpdodecg_c + ppi_c<br>ctp_c + h_c + pa181_c --> ppi_c + cdpdodec11eg_c<br>h2o_c + h_c + datp_c --> nh4_c + ditp_c<br>ru5p__D_c --> for_c + h_c + db4p_c<br>atp_c + dann_c + co2_c --> dtbt_c + pi_c + adp_c + 3.0 h_c<br>h2o_c + chtbs6p_c --> acgam6p_c + acgam_c<br>dca_e <=> dca_p<br>dcmp_e <=> dcmp_p<br>h2o_c + dctp_c + h_c --> nh4_c + dutp_c<br>dcyt_c + h2o_c + h_c --> duri_c + nh4_c<br>h_p + dcyt_p --> dcyt_c + h_c<br>dcyt_e <=> dcyt_p<br>ddca_e --> ddca_p<br>atp_c + 2dh3dgal_c --> adp_c + 2dh3dgal6p_c + h_c<br>h_p + 2ddglcn_p <=> 2ddglcn_c + h_c<br>2ddglcn_e <=> 2ddglcn_p<br>atp_c + 2ddglcn_c --> 2ddg6p_c + h_c + adp_c<br>e4p_c + h2o_c + pep_c --> 2dda7p_c + pi_c<br>2dh3dgal6p_c <=> g3p_c + pyr_c<br>atp_c + dgmp_c <=> dgdp_c + adp_c<br>dgmp_e <=> dgmp_p<br>dgsn_p + h_p --> dgsn_c + h_c<br>dgsn_e <=> dgsn_p<br>h2o_c + 23dhacoa_c <=> 3hadpcoa_c<br>23dhmb_c --> h2o_c + 3mob_c<br>23dhmp_c --> h2o_c + 3mop_c<br>pep_c + dha_c --> dhap_c + pyr_c<br>dha_e <=> dha_p<br>dha_p <=> dha_c<br>nad_c + 23ddhb_c <=> 23dhb_c + h_c + nadh_c<br>atp_c + 23dhb_c + h_c --> ppi_c + 23dhba_c<br>h2o_c + 23dhbzs_c --> 23dhb_c + ser__L_c<br>cenchddd_c + nad_c --> nadh_c + dhcinnm_c + h_c<br>o2_c + dhcinnm_c --> hkntd_c + h_c<br>nadph_c + 23dhdp_c + h_c --> thdp_c + nadp_c<br>aspsa_c + pyr_c --> 2.0 h2o_c + h_c + 23dhdp_c<br>nadph_c + h_c + dhf_c <=> thf_c + nadp_c<br>atp_c + dhpt_c + glu__L_c --> pi_c + adp_c + dhf_c + h_c<br>nadph_c + dhmpt_c + h_c --> thmnp_c + nadp_c<br>octdp_c + h_c + dhna_c --> 2dmmql8_c + ppi_c + co2_c<br>h_c + sbzcoa_c --> h2o_c + 14dhncoa_c<br>14dhncoa_c + h2o_c --> coa_c + h_c + dhna_c<br>dhnpt_c <=> 6hmhpt_c + gcald_c<br>dhnpt_c <=> dhmpt_c<br>dhor__S_c + q8_c --> orot_c + q8h2_c<br>mqn8_c + dhor__S_c --> orot_c + mql8_c<br>fum_c + dhor__S_c --> orot_c + succ_c<br>h2o_c + dhor__S_c <=> cbasp_c + h_c<br>nad_c + cechddd_c --> nadh_c + dhpppn_c + h_c<br>25drapp_c + h2o_c + h_c --> 5apru_c + nh4_c<br>6hmhptpp_c + 4abz_c --> dhpt_c + ppi_c<br>dhptd_c --> mdhdhf_c<br>nadph_c + dhptdn_c + 3.0 h_c --> nadp_c + thptdn_c<br>dhptdn_c + 3.0 h_c + nadh_c --> nad_c + thptdn_c<br>ahdt_c <=> dhmptp_c<br>2dda7p_c --> pi_c + 3dhq_c<br>3dhq_c --> h2o_c + 3dhsk_c<br>dimp_e <=> dimp_p<br>h_p + din_p --> din_c + h_c<br>din_e <=> din_p<br>nadph_c + 25dkglcn_c + h_c --> nadp_c + 2dhguln_c<br>25dkglcn_c + h_c + nadh_c --> 5dglcn_c + nad_c<br>nadph_c + 25dkglcn_c + h_c --> nadp_c + 5dglcn_c<br>dmpp_c + ipdp_c --> grdp_c + ppi_c<br>h2mb4p_c + nadh_c + h_c --> dmpp_c + nad_c + h2o_c<br>amet_c + 2omhmbl_c --> q8h2_c + h_c + ahcys_c<br>dmso_c + mql8_c --> dms_c + h2o_c + mqn8_c<br>dmso_p + mql8_c --> h2o_p + mqn8_c + dms_p<br>2dmmql8_c + dmso_c --> dms_c + 2dmmq8_c + h2o_c<br>2dmmql8_c + dmso_p --> h2o_p + 2dmmq8_c + dms_p<br>dmso_e <=> dmso_p<br>dmso_p <=> dmso_c<br>dms_e <=> dms_p<br>dhpmp_c + h2o_c --> pi_c + dhnpt_c<br>ahdt_c + h2o_c --> ppi_c + dhpmp_c + h_c<br>23doguln_c + h_c + nadh_c --> nad_c + 3dhguln_c<br>dopa_e <=> dopa_p<br>doxrbcn_e <=> doxrbcn_p<br>h_p + doxrbcn_p --> doxrbcn_e + h_c<br>atp_c + dpcoa_c --> adp_c + coa_c + h_c<br>nadph_c + 2dhp_c + h_c --> nadp_c + pant__R_c<br>2dr5p_c --> g3p_c + acald_c<br>dsbard_p + q8_c --> q8h2_c + dsbaox_p<br>dsbard_p + mqn8_c --> dsbaox_p + mql8_c<br>dsbcox_p + 2.0 gthrd_p --> gthox_p + dsbcrd_p<br>trdrd_c + dsbdox_c --> trdox_c + dsbdrd_c<br>dsbgox_p + 2.0 gthrd_p --> gthox_p + dsbgrd_p<br>ser__D_c + nadp_c <=> nadph_c + h_c + 2amsa_c<br>ser__D_p + h_p --> ser__D_c + h_c<br>ser__D_e <=> ser__D_p<br>tartr__D_c --> oaa_c + h2o_c<br>atp_c + dtmp_c <=> dtdp_c + adp_c<br>dtmp_e <=> dtmp_p<br>dump_e <=> dump_p<br>56dura_c + nad_c <=> ura_c + nadh_c + h_c<br>atp_c + duri_c --> adp_c + h_c + dump_c<br>duri_c + pi_c <=> 2dr1p_c + ura_c<br>duri_p + h_p --> duri_c + h_c<br>duri_e <=> duri_p<br>h2o_c + dutp_c --> ppi_c + h_c + dump_c<br>nadph_c + h_c + dxyl5p_c --> 2me4p_c + nadp_c<br>g3p_c + pyr_c + h_c --> co2_c + dxyl5p_c<br>atp_c + dxyl_c --> h_c + adp_c + dxyl5p_c<br>e4p_c + nad_c + h2o_c <=> 4per_c + nadh_c + 2.0 h_c<br>tdec2eACP_c + h_c + nadh_c --> nad_c + dcaACP_c<br>nadph_c + h_c + tdec2eACP_c --> dcaACP_c + nadp_c<br>tddec2eACP_c + h_c + nadh_c --> nad_c + ddcaACP_c<br>nadph_c + tddec2eACP_c + h_c --> nadp_c + ddcaACP_c<br>t3c5ddeceACP_c + nadh_c + h_c --> nad_c + cddec5eACP_c<br>nadph_c + t3c5ddeceACP_c + h_c --> nadp_c + cddec5eACP_c<br>tmrs2eACP_c + h_c + nadh_c --> nad_c + myrsACP_c<br>nadph_c + h_c + tmrs2eACP_c --> myrsACP_c + nadp_c<br>nadh_c + t3c7mrseACP_c + h_c --> tdeACP_c + nad_c<br>nadph_c + h_c + t3c7mrseACP_c --> tdeACP_c + nadp_c<br>tpalm2eACP_c + nadh_c + h_c --> palmACP_c + nad_c<br>nadph_c + tpalm2eACP_c + h_c --> palmACP_c + nadp_c<br>t3c9palmeACP_c + h_c + nadh_c --> nad_c + hdeACP_c<br>nadph_c + t3c9palmeACP_c + h_c --> hdeACP_c + nadp_c<br>nadh_c + toctd2eACP_c + h_c --> nad_c + ocdcaACP_c<br>nadph_c + h_c + toctd2eACP_c --> ocdcaACP_c + nadp_c<br>t3c11vaceACP_c + nadh_c + h_c --> nad_c + octeACP_c<br>nadph_c + h_c + t3c11vaceACP_c --> nadp_c + octeACP_c<br>nadh_c + h_c + but2eACP_c --> nad_c + butACP_c<br>nadph_c + but2eACP_c + h_c --> butACP_c + nadp_c<br>h_c + nadh_c + thex2eACP_c --> nad_c + hexACP_c<br>nadph_c + h_c + thex2eACP_c --> hexACP_c + nadp_c<br>toct2eACP_c + h_c + nadh_c --> nad_c + ocACP_c<br>nadph_c + toct2eACP_c + h_c --> ocACP_c + nadp_c<br>atp_c + eca4colipa_p + h2o_c --> adp_c + eca4colipa_e + pi_c + h_c<br>eca4und_p + colipa_p --> h_p + eca4colipa_p + udcpdp_p<br>2.0 unagamuf_p --> h_p + eca2und_p + udcpdp_p<br>unagamuf_p + eca2und_p --> h_p + udcpdp_p + eca3und_p<br>unagamuf_p + eca3und_p --> h_p + udcpdp_p + eca4und_p<br>unagamuf_c --> unagamuf_p<br>3hbcoa_c <=> b2coa_c + h2o_c<br>3hhcoa_c <=> hx2coa_c + h2o_c<br>3hocoa_c <=> h2o_c + oc2coa_c<br>3hdcoa_c <=> h2o_c + dc2coa_c<br>3hddcoa_c <=> dd2coa_c + h2o_c<br>3htdcoa_c <=> h2o_c + td2coa_c<br>3hhdcoa_c <=> h2o_c + hdd2coa_c<br>3hodcoa_c <=> h2o_c + od2coa_c<br>2ddg6p_c --> g3p_c + pyr_c<br>6pgc_c --> 2ddg6p_c + h2o_c<br>kdo2lipid4_c + ddcaACP_c --> ACP_c + kdo2lipid4L_c<br>myrsACP_c + kdo2lipid4L_c --> ACP_c + lipa_c<br>kdo2lipid4_c + hdeACP_c --> kdo2lipid4p_c + ACP_c<br>kdo2lipid4p_c + myrsACP_c --> lipa_cold_c + ACP_c<br>nadph_c + egmeACP_c + h_c --> nadp_c + gmeACP_c<br>atp_c + enlipa_p + h2o_c --> adp_c + pi_c + enlipa_e + h_c<br>2pg_c <=> h2o_c + pep_c<br>3.0 seramp_c + 3.0 23dhba_c --> enter_c + 6.0 amp_c + 9.0 h_c<br>enter_c + 3.0 h2o_c --> 3.0 h_c + 3.0 23dhbzs_c<br>feenter_c + 3.0 h2o_c --> fe3_c + 3.0 h_c + 3.0 23dhbzs_c<br>nadph_c + h_c + epmeACP_c --> pmeACP_c + nadp_c<br>etha_c --> acald_c + nh4_c<br>h_p + etha_p --> etha_c + h_c<br>etha_e <=> etha_p<br>atp_c + ethso3_p + h2o_c --> adp_c + pi_c + ethso3_c + h_c<br>ethso3_e <=> ethso3_p<br>etoh_e <=> etoh_p<br>etoh_p <=> etoh_c<br>f6p_c <=> g3p_c + dha_c<br>f6p_c + h2o_c --> fru_c + pi_c<br>f6p_p + 2.0 pi_c --> f6p_c + 2.0 pi_p<br>f6p_e <=> f6p_p<br>h2o_c + dcaACP_c --> dca_c + ACP_c + h_c<br>h2o_c + ddcaACP_c --> ACP_c + h_c + ddca_c<br>h2o_c + myrsACP_c --> ACP_c + h_c + ttdca_c<br>tdeACP_c + h2o_c --> ACP_c + h_c + ttdcea_c<br>palmACP_c + h2o_c --> ACP_c + hdca_c + h_c<br>hdeACP_c + h2o_c --> hdcea_c + h_c + ACP_c<br>h2o_c + ocACP_c --> octa_c + ACP_c + h_c<br>h2o_c + dcacoa_c --> dca_c + coa_c + h_c<br>h2o_c + ddcacoa_c --> coa_c + h_c + ddca_c<br>tdcoa_c + h2o_c --> coa_c + h_c + ttdca_c<br>h2o_c + tdecoa_c --> coa_c + h_c + ttdcea_c<br>pmtcoa_c + h2o_c --> coa_c + hdca_c + h_c<br>h2o_c + hdcoa_c --> coa_c + hdcea_c + h_c<br>stcoa_c + h2o_c --> ocdca_c + coa_c + h_c<br>odecoa_c + h2o_c --> coa_c + ocdcea_c + h_c<br>h2o_c + hxcoa_c --> coa_c + hxa_c + h_c<br>h2o_c + occoa_c --> octa_c + coa_c + h_c<br>atp_c + h_p + coa_c + dca_p --> ppi_c + amp_c + h_c + dcacoa_c<br>atp_c + h_p + coa_c + ddca_p --> ppi_c + amp_c + ddcacoa_c + h_c<br>atp_c + coa_c + h_p + ttdca_p --> tdcoa_c + ppi_c + amp_c + h_c<br>atp_c + h_p + coa_c + ttdcea_p --> ppi_c + amp_c + h_c + tdecoa_c<br>h_p + coa_c + hdca_p + atp_c --> pmtcoa_c + ppi_c + amp_c + h_c<br>atp_c + h_p + coa_c + hdcea_p --> ppi_c + hdcoa_c + amp_c + h_c<br>atp_c + h_p + coa_c + ocdca_p --> ppi_c + stcoa_c + amp_c + h_c<br>atp_c + h_p + coa_c + ocdcea_p --> odecoa_c + ppi_c + amp_c + h_c<br>atp_c + h_p + coa_c + hxa_p --> ppi_c + amp_c + hxcoa_c + h_c<br>atp_c + h_p + coa_c + octa_p --> ppi_c + amp_c + occoa_c + h_c<br>fad_c + h_c + nadh_c --> nad_c + fadh2_c<br>nadph_c + fad_c + h_c --> nadp_c + fadh2_c<br>nad_c + hmgth_c <=> nadh_c + Sfglutth_c + h_c<br>fald_e <=> fald_p<br>fald_p <=> fald_c<br>fald_c + gthrd_c <=> hmgth_c<br>fdp_c <=> g3p_c + dhap_c<br>s17bp_c <=> e4p_c + dhap_c<br>h2o_c + fdp_c --> f6p_c + pi_c<br>fuc__L_c <=> fcl__L_c<br>atp_c + fcl__L_c --> h_c + fc1p_c + adp_c<br>fc1p_c <=> lald__L_c + dhap_c<br>fe2_c + ppp9_c --> pheme_c + 2.0 h_c<br>q8_c + for_p + 2.0 h_c --> h_p + co2_p + q8h2_c<br>mqn8_c + for_p + 2.0 h_c --> h_p + co2_p + mql8_c<br>o2_c + fmnh2_c + isetac_c --> fmn_c + h2o_c + gcald_c + so3_c + h_c<br>o2_c + fmnh2_c + mso3_c --> fmn_c + fald_c + h2o_c + so3_c + h_c<br>fmnh2_c + o2_c + ethso3_c --> fmn_c + acald_c + so3_c + h2o_c + h_c<br>o2_c + butso3_c + fmnh2_c --> fmn_c + h2o_c + btal_c + so3_c + h_c<br>sulfac_c + fmnh2_c + o2_c --> fmn_c + glx_c + h2o_c + so3_c + h_c<br>atp_c + fe2_p + h2o_c --> adp_c + pi_c + fe2_c + h_c<br>h_p + fe2_p --> fe2_c + h_c<br>fe2_c + h_p --> fe2_p + h_c<br>fe2_e <=> fe2_p<br>fe2_p --> fe2_c<br>atp_c + fe3dcit_p + h2o_c --> adp_c + fe3_c + pi_c + 2.0 cit_c + h_c<br>fe3dcit_e + h_p --> fe3dcit_p + h_c<br>fe3dhbzs_c --> 23dhbzs_c + fe3_c<br>atp_c + h2o_c + fe3dhbzs_p --> adp_c + pi_c + fe3dhbzs_c + h_c<br>h_p + fe3dhbzs_e --> h_c + fe3dhbzs_p<br>2.0 fe3hox_c + fadh2_c --> 2.0 fe2_c + fad_c + 2.0 h_c + 2.0 fe3hox_un_c<br>fmnh2_c + 2.0 fe3hox_c --> 2.0 fe2_c + fmn_c + 2.0 h_c + 2.0 fe3hox_un_c<br>rbflvrd_c + 2.0 fe3hox_c --> 2.0 fe2_c + ribflv_c + 2.0 h_c + 2.0 fe3hox_un_c<br>h_p + fe3hox_un_c --> fe3hox_un_p + h_c<br>h_p + fe3hox_un_p --> h_c + fe3hox_un_e<br>atp_c + fe3hox_p + h2o_c --> adp_c + pi_c + fe3hox_c + h_c<br>fe3hox_un_e + fe3_e --> fe3hox_e<br>fe3hox_e + h_p --> fe3hox_p + h_c<br>2.0 fe3_c + fadh2_c --> 2.0 fe2_c + fad_c + 2.0 h_c<br>atp_c + fe3_p + h2o_c --> adp_c + fe3_c + pi_c + h_c<br>fe3_e <=> fe3_p<br>2.0 fecrm_c + fadh2_c --> 2.0 fe2_c + 2.0 fecrm_un_c + fad_c + 2.0 h_c<br>2.0 fecrm_c + fmnh2_c --> 2.0 fe2_c + fmn_c + 2.0 fecrm_un_c + 2.0 h_c<br>rbflvrd_c + 2.0 fecrm_c --> 2.0 fe2_c + 2.0 fecrm_un_c + ribflv_c + 2.0 h_c<br>h_p + fecrm_un_p --> fecrm_un_e + h_c<br>h_p + fecrm_un_c --> fecrm_un_p + h_c<br>atp_c + fecrm_p + h2o_c --> adp_c + fecrm_c + pi_c + h_c<br>fecrm_un_e + fe3_e --> fecrm_e<br>h_p + fecrm_e --> h_c + fecrm_p<br>2.0 feenter_c + fadh2_c --> 2.0 fe2_c + 2.0 enter_c + fad_c + 2.0 h_c<br>fmnh2_c + 2.0 feenter_c --> 2.0 fe2_c + fmn_c + 2.0 enter_c + 2.0 h_c<br>rbflvrd_c + 2.0 feenter_c --> 2.0 fe2_c + 2.0 enter_c + ribflv_c + 2.0 h_c<br>atp_c + h2o_c + feenter_p --> adp_c + pi_c + feenter_c + h_c<br>enter_e + fe3_e --> feenter_e<br>h_p + enter_p --> h_c + enter_e<br>h_p + feenter_e --> feenter_p + h_c<br>h_p + enter_c --> enter_p + h_c<br>2.0 feoxam_c + fadh2_c --> 2.0 feoxam_un_c + 2.0 fe2_c + fad_c + 2.0 h_c<br>fmnh2_c + 2.0 feoxam_c --> 2.0 feoxam_un_c + 2.0 fe2_c + fmn_c + 2.0 h_c<br>rbflvrd_c + 2.0 feoxam_c --> 2.0 feoxam_un_c + 2.0 fe2_c + ribflv_c + 2.0 h_c<br>h_p + feoxam_un_p --> h_c + feoxam_un_e<br>feoxam_un_c + h_p --> feoxam_un_p + h_c<br>atp_c + feoxam_p + h2o_c --> adp_c + pi_c + feoxam_c + h_c<br>feoxam_un_e + fe3_e --> feoxam_e<br>h_p + feoxam_e --> h_c + feoxam_p<br>4.0 h_p + 4.0 fe2_p + o2_p --> 2.0 h2o_p + 4.0 fe3_p<br>2.0 h_c + 2.0 4fe4s_c + h2o2_c --> 2.0 fe3_c + 2.0 3fe4s_c + 2.0 h2o_c<br>2.0 no_c + 2.0 4fe4s_c + 2.0 h_c --> n2o_c + 2.0 fe3_c + 2.0 3fe4s_c + h2o_c<br>3fe4s_c + fe2_c --> 4fe4s_c<br>h2o_c + suc6p_c --> fru_c + g6p_c<br>for_c + h_c --> h2_c + co2_c<br>nadph_c + 2.0 flxso_c --> nadp_c + h_c + 2.0 flxr_c<br>nadph_c + ribflv_c + h_c --> rbflvrd_c + nadp_c<br>ribflv_c + h_c + nadh_c --> nad_c + rbflvrd_c<br>mettrna_c + 10fthf_c --> thf_c + fmettrna_c + h_c<br>atp_c + fmn_c + h_c --> ppi_c + fad_c<br>fmn_c + nadh_c + h_c --> nad_c + fmnh2_c<br>nadph_c + fmn_c + h_c --> nadp_c + fmnh2_c<br>h_c + 5fthf_c --> h2o_c + methf_c<br>forcoa_c + oxa_c <=> for_c + oxalcoa_c<br>h_p + for_p --> for_c + h_c<br>for_e <=> for_p<br>for_c --> for_p<br>mql8_c + fum_c --> mqn8_c + succ_c<br>2dmmql8_c + fum_c --> 2dmmq8_c + succ_c<br>atp_c + f1p_c --> adp_c + fdp_c + h_c<br>h2o_c + frulysp_c <=> g6p_c + lys__L_c<br>psclys_c <=> frulys_c<br>atp_c + frulys_c --> frulysp_c + h_c + adp_c<br>frulys_p + h_p --> frulys_c + h_c<br>frulys_e <=> frulys_p<br>h_p + fruur_p <=> fruur_c + h_c<br>fruur_e <=> fruur_p<br>fru_p + pep_c --> f6p_c + pyr_c<br>fru_p + pep_c --> pyr_c + f1p_c<br>fru_e <=> fru_p<br>h2o_c + 10fthf_c --> for_c + thf_c + h_c<br>atp_c + for_c + thf_c --> adp_c + pi_c + 10fthf_c<br>fuc__L_e <=> fuc__L_p<br>h_p + fuc__L_p <=> fuc__L_c + h_c<br>h2o_c + fum_c <=> mal__L_c<br>2.0 h_p + fum_p --> fum_c + 2.0 h_c<br>3.0 h_p + fum_p --> fum_c + 3.0 h_c<br>fum_e <=> fum_p<br>fusa_e <=> fusa_p<br>h_p + fusa_p --> fusa_e + h_c<br>gam1p_c + accoa_c --> acgam1p_c + coa_c + h_c<br>h2o_p + g1p_p --> glc__D_p + pi_p<br>dttp_c + g1p_c + h_c --> ppi_c + dtdpglu_c<br>g1p_e <=> g1p_p<br>glu1sa_c <=> 5aop_c<br>h2o_c + glyc2p_c --> glyc_c + pi_c<br>h2o_p + glyc2p_p --> pi_p + glyc_p<br>ddcaACP_c + glyc3p_c --> 1ddecg3p_c + ACP_c<br>myrsACP_c + glyc3p_c --> 1tdecg3p_c + ACP_c<br>tdeACP_c + glyc3p_c --> ACP_c + 1tdec7eg3p_c<br>palmACP_c + glyc3p_c --> 1hdecg3p_c + ACP_c<br>hdeACP_c + glyc3p_c --> 1hdec9eg3p_c + ACP_c<br>ocdcaACP_c + glyc3p_c --> 1odecg3p_c + ACP_c<br>octeACP_c + glyc3p_c --> 1odec11eg3p_c + ACP_c<br>atp_c + h2o_c + g3pc_p --> adp_c + pi_c + g3pc_c + h_c<br>g3pc_e <=> g3pc_p<br>nadp_c + glyc3p_c <=> nadph_c + dhap_c + h_c<br>q8_c + glyc3p_c --> q8h2_c + dhap_c<br>mqn8_c + glyc3p_c --> dhap_c + mql8_c<br>2dmmq8_c + glyc3p_c --> 2dmmql8_c + dhap_c<br>atp_c + g3pe_p + h2o_c --> adp_c + pi_c + g3pe_c + h_c<br>g3pe_e <=> g3pe_p<br>atp_c + h2o_c + g3pg_p --> adp_c + pi_c + g3pg_c + h_c<br>g3pg_e <=> g3pg_p<br>atp_c + h2o_c + g3pi_p --> adp_c + g3pi_c + pi_c + h_c<br>g3pi_e <=> g3pi_p<br>atp_c + g3ps_p + h2o_c --> adp_c + pi_c + g3ps_c + h_c<br>g3ps_e <=> g3ps_p<br>h2o_c + glyc3p_c --> glyc_c + pi_c<br>glu5sa_c --> 1pyr5c_c + h2o_c + h_c<br>glu5p_c + nadph_c + h_c --> glu5sa_c + pi_c + nadp_c<br>h2o_c + gam6p_c --> f6p_c + nh4_c<br>nadp_c + g6p_c <=> nadph_c + 6pgl_c + h_c<br>g6p_c + h2o_c --> glc__D_c + pi_c<br>2.0 pi_c + g6p_p --> 2.0 pi_p + g6p_c<br>g6p_e <=> g6p_p<br>h2o_p + gal1p_p --> gal_p + pi_p<br>gal1p_e <=> gal1p_p<br>gal_bD_e <=> gal_bD_p<br>galct__D_c --> 5dh4dglc_c + h2o_c<br>galctn__L_c + nad_c --> tagur_c + nadh_c + h_c<br>galctn__D_c --> h2o_c + 2dh3dgal_c<br>h_p + galctn__L_p --> galctn__L_c + h_c<br>galctn__L_e <=> galctn__L_p<br>galctn__D_p + h_p --> galctn__D_c + h_c<br>galctn__D_e <=> galctn__D_p<br>h_p + galct__D_p <=> h_c + galct__D_c<br>galct__D_e <=> galct__D_p<br>atp_c + gal_c <=> adp_c + h_c + gal1p_c<br>gal_bD_p --> gal_p<br>h2o_c + melib_c --> glc__D_c + gal_c<br>udpg_c + gicolipa_c --> udp_c + gagicolipa_c + h_c<br>pep_c + galt_p --> pyr_c + galt1p_c<br>galt_e <=> galt_p<br>h_p + galur_p <=> h_c + galur_c<br>galur_e <=> galur_p<br>utp_c + h_c + g1p_c --> ppi_c + udpg_c<br>atp_c + gal_p + h2o_c --> adp_c + gal_c + pi_c + h_c<br>h_p + gal_p --> gal_c + h_c<br>gal_e <=> gal_p<br>gam6p_p + 2.0 pi_c --> 2.0 pi_p + gam6p_c<br>gam6p_e <=> gam6p_p<br>pep_c + gam_p --> pyr_c + gam6p_c<br>gam_e <=> gam_p<br>g3p_c + pi_c + nad_c <=> 13dpg_c + nadh_c + h_c<br>gar_c + 10fthf_c <=> thf_c + fgam_c + h_c<br>atp_c + gar_c + for_c --> adp_c + fgam_c + pi_c + h_c<br>gbbtn_e <=> gbbtn_p<br>nad_c + h2o_c + gcald_c --> glyclt_c + 2.0 h_c + nadh_c<br>gdpddman_c --> gdpofuc_c<br>atp_c + gdp_c --> amp_c + ppgpp_c + h_c<br>gdpmann_c + h2o_c --> gdp_c + h_c + man_c<br>gdpmann_c + h2o_c --> gmp_c + 2.0 h_c + man1p_c<br>gdptp_c + h2o_c --> ppi_c + gtp_c<br>gdp_e <=> gdp_p<br>f6p_c + gln__L_c --> glu__L_c + gam6p_c<br>nadp_c + h2o_c + ggbutal_c <=> nadph_c + gg4abut_c + 2.0 h_c<br>gg4abut_c + h2o_c --> glu__L_c + 4abut_c<br>o2_c + ggptrc_c + h2o_c --> ggbutal_c + nh4_c + h2o2_c<br>atp_c + glu__L_c + ptrc_c --> adp_c + pi_c + ggptrc_c + h_c<br>sucsal_c + nadh_c + h_c <=> nad_c + ghb_c<br>thf_c + ser__L_c <=> gly_c + h2o_c + mlthf_c<br>atp_c + gmp_c <=> gdp_c + adp_c<br>glycogen_c --> bglycogen_c<br>glc__D_c + accoa_c <=> coa_c + acglc__D_c<br>q8_c + h2o_p + glc__D_p --> h_p + q8h2_c + glcn_p<br>h_p + glcn_p <=> glcn_c + h_c<br>glcn_e <=> glcn_p<br>pi_c + glycogen_c --> g1p_c<br>pi_c + bglycogen_c --> g1p_c<br>5dh4dglc_c --> 2h3oppan_c + pyr_c<br>glcr_c --> 5dh4dglc_c + h2o_c<br>h_p + glcr_p <=> glcr_c + h_c<br>glcr_e <=> glcr_p<br>adpglc_c --> adp_c + glycogen_c + h_c<br>icolipa_c + udpg_c --> udp_c + gicolipa_c + h_c<br>udpg_c + gagicolipa_c --> udp_c + ggagicolipa_c + h_c<br>ggagicolipa_c + udpg_c --> udp_c + gggagicolipa_c + h_c<br>glcur1p_e <=> glcur1p_p<br>h_p + glcur_p <=> glcur_c + h_c<br>glcur_e <=> glcur_p<br>atp_c + glc__D_p + h2o_c --> glc__D_c + pi_c + adp_c + h_c<br>glc__D_p + pep_c --> g6p_c + pyr_c<br>glc__D_p + h_p --> h_c + glc__D_c<br>glc__D_e <=> glc__D_p<br>glc__D_e --> glc__D_p<br>bglycogen_c --> glycogen_c<br>atp_c + g1p_c + h_c --> ppi_c + adpglc_c<br>atp_c + glu__L_c + nh4_c --> adp_c + pi_c + h_c + gln__L_c<br>atp_c + trnagln_c + gln__L_c --> glntrna_c + amp_c + ppi_c<br>atp_c + gln__L_p + h2o_c --> adp_c + pi_c + h_c + gln__L_c<br>gln__L_e <=> gln__L_p<br>nad_c + galt1p_c <=> nadh_c + h_c + tag6p__D_c<br>atp_c + glu__L_c --> glu5p_c + adp_c<br>4abut_c + glu__L_p <=> glu__L_c + 4abut_p<br>atp_c + glu__L_c + cys__L_c --> adp_c + pi_c + glucys_c + h_c<br>glu__L_c + h_c --> 4abut_c + co2_c<br>nadp_c + glu__L_c + h2o_c <=> akg_c + nadph_c + h_c + nh4_c<br>h2o_c + gln__L_c --> nh4_c + glu__L_c<br>h2o_p + gln__L_p --> nh4_p + glu__L_p<br>gln__L_c + h2o_c + prpp_c --> ppi_c + pram_c + glu__L_c<br>glu__D_c <=> glu__L_c<br>nadph_c + akg_c + h_c + gln__L_c --> nadp_c + 2.0 glu__L_c<br>nadph_c + glutrna_c + h_c --> nadp_c + trnaglu_c + glu1sa_c<br>atp_c + trnaglu_c + glu__L_c --> ppi_c + glutrna_c + amp_c<br>atp_c + h2o_c + glu__L_p --> adp_c + pi_c + glu__L_c + h_c<br>h_p + glu__L_p <=> glu__L_c + h_c<br>glu__L_p + na1_p --> glu__L_c + na1_c<br>glu__L_e <=> glu__L_p<br>h_c + 2.0 glx_c --> 2h3oppan_c + co2_c<br>glyald_e <=> glyald_p<br>glyald_p <=> glyald_c<br>accoa_c + gly_c <=> coa_c + 2aobut_c<br>atp_c + glyb_p + h2o_c --> adp_c + pi_c + glyb_c + h_c<br>h_p + glyb_p --> h_c + glyb_c<br>glyb_e <=> glyb_p<br>atp_c + glyc2p_p + h2o_c --> adp_c + glyc2p_c + pi_c + h_c<br>glyc2p_e <=> glyc2p_p<br>atp_c + h2o_c + glyc3p_p --> adp_c + pi_c + glyc3p_c + h_c<br>glyc3p_p + pi_c --> pi_p + glyc3p_c<br>glyc3p_e <=> glyc3p_p<br>h_p + glyc__R_p <=> glyc__R_c + h_c<br>glyc__R_e <=> glyc__R_p<br>nad_c + glyc_c --> dha_c + h_c + nadh_c<br>atp_c + glyc__R_c --> adp_c + 3pg_c + h_c<br>atp_c + glyc__R_c --> 2pg_c + h_c + adp_c<br>gly_c + thf_c + nad_c --> nh4_c + nadh_c + mlthf_c + co2_c<br>h_c + glx_c + nadh_c --> nad_c + glyclt_c<br>nadph_c + h_c + glx_c --> glyclt_c + nadp_c<br>h_p + glyclt_p <=> glyclt_c + h_c<br>glyclt_p + na1_p --> glyclt_c + na1_c<br>glyclt_e <=> glyclt_p<br>q8_c + glyclt_c --> q8h2_c + glx_c<br>glyclt_c + mqn8_c --> glx_c + mql8_c<br>glyclt_c + 2dmmq8_c --> 2dmmql8_c + glx_c<br>glyc_e <=> glyc_p<br>glyc_c <=> glyc_p<br>atp_c + glyc_c --> adp_c + h_c + glyc3p_c<br>h2o_c + lgt__S_c --> lac__D_c + h_c + gthrd_c<br>h2o_c + mthgxl_c --> lac__D_c + h_c<br>atp_c + gly_c + trnagly_c --> glytrna_c + ppi_c + amp_c<br>gly_p + h_p --> gly_c + h_c<br>gly_p + h_p <=> gly_c + h_c<br>gly_p + na1_p --> gly_c + na1_c<br>gly_e <=> gly_p<br>gdpmann_c --> h2o_c + gdpddman_c<br>atp_c + gmhep1p_c + h_c --> adphep_DD_c + ppi_c<br>atp_c + gmhep7p_c --> adp_c + gmhep17bp_c + h_c<br>gmhep17bp_c + h2o_c --> gmhep1p_c + pi_c<br>gmp_c + nadph_c + 2.0 h_c --> nadp_c + imp_c + nh4_c<br>atp_c + h2o_c + xmp_c + gln__L_c --> amp_c + gmp_c + glu__L_c + ppi_c + 2.0 h_c<br>gmp_e <=> gmp_p<br>6pgc_c + nadp_c --> nadph_c + ru5p__D_c + co2_c<br>glcn_c + atp_c --> 6pgc_c + adp_c + h_c<br>nadph_c + h_c + gdpofuc_c --> nadp_c + gdpfuc_c<br>h2o_c + gp4g_c --> 2.0 gdp_c + 2.0 h_c<br>h2o_c + g3pc_c --> chol_c + h_c + glyc3p_c<br>h2o_p + g3pc_p --> h_p + chol_p + glyc3p_p<br>g3pe_c + h2o_c --> etha_c + h_c + glyc3p_c<br>h2o_p + g3pe_p --> h_p + glyc3p_p + etha_p<br>g3ps_c + h2o_c --> ser__L_c + h_c + glyc3p_c<br>h2o_p + g3ps_p --> h_p + ser__L_p + glyc3p_p<br>h2o_c + g3pg_c --> glyc3p_c + glyc_c + h_c<br>h2o_p + g3pg_p --> h_p + glyc_p + glyc3p_p<br>g3pi_c + h2o_c --> inost_c + h_c + glyc3p_c<br>h2o_p + g3pi_p --> h_p + glyc3p_p + inost_p<br>grdp_c + ipdp_c --> ppi_c + frdp_c<br>grxox_c + 2.0 gthrd_c --> grxrd_c + gthox_c<br>gsn_c + atp_c --> gmp_c + h_c + adp_c<br>h_p + gsn_p --> gsn_c + h_c<br>gsn_e <=> gsn_p<br>h2o_c + gtspmd_c --> spmd_c + gthrd_c<br>atp_c + spmd_c + gthrd_c --> adp_c + pi_c + gtspmd_c + h_c<br>gthox_e <=> gthox_p<br>nadph_c + gthox_c + h_c <=> nadp_c + 2.0 gthrd_c<br>2.0 gthrd_c + h2o2_c --> 2.0 h2o_c + gthox_c<br>h2o_p + gthrd_p --> glu__L_p + cgly_p<br>atp_c + h2o_c + gthrd_c --> adp_c + pi_c + gthrd_p + h_c<br>atp_c + gthrd_p + h2o_c --> adp_c + pi_c + h_c + gthrd_c<br>gthrd_e <=> gthrd_p<br>atp_c + gly_c + glucys_c --> adp_c + pi_c + h_c + gthrd_c<br>gtp_c + h2o_c --> for_c + ahdt_c + h_c<br>gtp_c + 3.0 h2o_c --> for_c + ppi_c + 25drapp_c + 2.0 h_c<br>gdptp_c + h2o_c --> ppgpp_c + pi_c + h_c<br>atp_c + gtp_c --> gdptp_c + amp_c + h_c<br>gtp_c + h2o_c + h_c --> xtp_c + nh4_c<br>gtp_e <=> gtp_p<br>gtp_c --> 35cgmp_c + ppi_c<br>gua_c + h2o_c + h_c --> nh4_c + xan_c<br>gua_c + prpp_c --> ppi_c + gmp_c<br>h_p + gua_p --> h_c + gua_c<br>gua_e <=> gua_p<br>gua_p <=> gua_c<br>glcur_c <=> fruur_c<br>galur_c <=> tagur_c<br>h2o_p + glcur1p_p --> glcur_p + pi_p<br>h2o2_e <=> h2o2_p<br>h2o_e <=> h2o_p<br>h2o_p <=> h2o_c<br>h2s_c + 2.0 o2_c --> so4_c + 2.0 h_c<br>h2s_c --> h2s_p<br>h2s_e <=> h2s_p<br>h2_e <=> h2_p<br>h2_p <=> h2_c<br>aacoa_c + h_c + nadh_c <=> nad_c + 3hbcoa_c<br>nadh_c + h_c + 3ohcoa_c <=> nad_c + 3hhcoa_c<br>3oocoa_c + nadh_c + h_c <=> nad_c + 3hocoa_c<br>3odcoa_c + h_c + nadh_c <=> nad_c + 3hdcoa_c<br>3oddcoa_c + h_c + nadh_c <=> nad_c + 3hddcoa_c<br>3otdcoa_c + h_c + nadh_c <=> nad_c + 3htdcoa_c<br>nadh_c + h_c + 3ohdcoa_c <=> nad_c + 3hhdcoa_c<br>3oodcoa_c + nadh_c + h_c <=> nad_c + 3hodcoa_c<br>nad_c + 3hadpcoa_c <=> oxadpcoa_c + h_c + nadh_c<br>octdp_c + 4hbz_c --> 3ophb_c + ppi_c<br>h_p + 3hcinnm_p <=> 3hcinnm_c + h_c<br>3hcinnm_e <=> 3hcinnm_p<br>h2o_c + co2_c <=> h_c + hco3_c<br>hcys__L_c + amet_c --> met__L_c + h_c + ahcys_c<br>hcys__L_c + mmet_c --> 2.0 met__L_c + h_c<br>hdca_e --> hdca_p<br>hdcea_e --> hdcea_p<br>h2o_c + frdp_c + pheme_c --> ppi_c + hemeO_c<br>atp_c + hhlipa_c --> adp_c + phhlipa_c + h_c<br>atp_c + hphhlipa_c --> adp_c + phphhlipa_c + h_c<br>adphep_LD_c + lipa_c --> hlipa_c + adp_c + h_c<br>hlipa_c + adphep_LD_c --> adp_c + hhlipa_c + h_c<br>adphep_LD_c + phhlipa_c --> adp_c + hphhlipa_c + h_c<br>gggagicolipa_c + adphep_LD_c --> colipa_c + h_c + adp_c<br>atp_c + 4mhetz_c --> adp_c + 4mpetz_c + h_c<br>atp_c + glc__D_c --> g6p_c + h_c + adp_c<br>atp_c + man_c --> adp_c + h_c + man6p_c<br>atp_c + fru_c --> adp_c + h_c + f6p_c<br>h_p + hxa_p <=> hxa_c + h_c<br>atp_c + h2o_c + hg2_c --> adp_c + pi_c + hg2_p + h_c<br>h_p + hg2_c --> hg2_p + h_c<br>hg2_e <=> hg2_p<br>2.0 nad_c + h2o_c + histd_c --> his__L_c + 2.0 nadh_c + 3.0 h_c<br>hisp_c + h2o_c --> pi_c + histd_c<br>atp_c + his__L_c + trnahis_c --> histrna_c + ppi_c + amp_c<br>atp_c + h2o_c + his__L_p --> adp_c + pi_c + his__L_c + h_c<br>h_p + his__L_p <=> his__L_c + h_c<br>his__L_e <=> his__L_p<br>h2o_c + hkndd_c --> h_c + succ_c + op4en_c<br>hkntd_c + h2o_c --> fum_c + h_c + op4en_c<br>h2o_c + 4.0 ppbng_c --> hmbil_c + 4.0 nh4_c<br>4ahmmp_c + atp_c --> adp_c + 4ampm_c + h_c<br>h_p + hom__L_c --> h_c + hom__L_p<br>hom__L_e <=> hom__L_p<br>4h2opntn_c --> acald_c + pyr_c<br>atp_c + 6hmhpt_c --> amp_c + 6hmhptpp_c + h_c<br>dhpppn_c + o2_c --> h_c + hkndd_c<br>h_p + 3hpppn_p <=> 3hpppn_c + h_c<br>3hpppn_e <=> 3hpppn_p<br>hpyr_c <=> 2h3oppan_c<br>hpyr_c + h_c + nadh_c --> glyc__R_c + nad_c<br>nadph_c + hpyr_c + h_c --> glyc__R_c + nadp_c<br>nadp_c + hom__L_c <=> nadph_c + aspsa_c + h_c<br>atp_c + hom__L_c --> adp_c + h_c + phom_c<br>succoa_c + hom__L_c --> coa_c + suchms_c<br>glu__L_c + imacp_c --> hisp_c + akg_c<br>hxan_c + nad_c + h2o_c --> xan_c + h_c + nadh_c<br>hxa_e <=> hxa_p<br>accoa_c + hxa_c --> ac_c + hxcoa_c<br>prpp_c + hxan_c --> ppi_c + imp_c<br>h2_c + 2.0 h_c + q8_c --> 2.0 h_p + q8h2_c<br>mqn8_c + h2_c + 2.0 h_c --> 2.0 h_p + mql8_c<br>h2_c + 2.0 h_c + 2dmmq8_c --> 2.0 h_p + 2dmmql8_c<br>pyam5p_c + h2o_c --> pi_c + pydam_c<br>hxan_e <=> hxan_p<br>hxan_p <=> hxan_c<br>h_e <=> h_p<br>iscssh_c + 2fe1s_c + iscu_c --> iscu_2fe2s_c + iscs_c + 4.0 h_c<br>2.0 fe2_c + 2.0 iscssh_c + iscu_c + fadh2_c --> iscu_2fe2s_c + 2.0 iscs_c + fad_c + 6.0 h_c<br>2.0 fe2_c + iscu_2fe2s_c + 2.0 iscssh_c + fadh2_c --> iscu_2fe2s2_c + 2.0 iscs_c + fad_c + 6.0 h_c<br>iscu_2fe2s_c + 4.0 h_c --> iscu_c + 2fe2s_c<br>fadh2_c + 2.0 h_c + iscu_2fe2s2_c --> iscu_4fe4s_c + fad_c<br>iscu_4fe4s_c + 4.0 h_c --> iscu_c + 4fe4s_c<br>nadp_c + icit_c <=> nadph_c + akg_c + co2_c<br>chor_c <=> ichor_c<br>chor_c --> ichor_c<br>ichor_c + h2o_c --> pyr_c + 23ddhb_c<br>icit_c --> glx_c + succ_c<br>cys__L_c + iscs_c --> iscssh_c + ala__L_c<br>5dglcn_c + h_c + nadh_c <=> nad_c + idon__L_c<br>nadph_c + h_c + 5dglcn_c --> idon__L_c + nadp_c<br>h_p + idon__L_p <=> idon__L_c + h_c<br>idon__L_e <=> idon__L_p<br>gln__L_c + prlp_c --> eig3p_c + aicar_c + glu__L_c + h_c<br>eig3p_c --> imacp_c + h2o_c<br>h_c + 2cpr5p_c --> co2_c + h2o_c + 3ig3p_c<br>akg_c + ile__L_c <=> glu__L_c + 3mop_c<br>atp_c + ile__L_c + trnaile_c --> ppi_c + iletrna_c + amp_c<br>atp_c + h2o_c + ile__L_p --> adp_c + pi_c + ile__L_c + h_c<br>h_p + ile__L_p <=> ile__L_c + h_c<br>ile__L_e <=> ile__L_p<br>h2o_c + imp_c <=> fprica_c<br>imp_c + nad_c + h2o_c --> nadh_c + xmp_c + h_c<br>imp_e <=> imp_p<br>h_c + indole_c --> h_p + indole_p<br>h_p + indole_p <=> h_c + indole_c<br>indole_e <=> indole_p<br>inost_p + na1_p --> inost_c + na1_c<br>ins_c + h2o_c --> rib__D_c + hxan_c<br>atp_c + ins_c --> h_c + imp_c + adp_c<br>inost_e <=> inost_p<br>h_p + ins_p --> ins_c + h_c<br>h_p + ins_p <=> ins_c + h_c<br>ins_e <=> ins_p<br>ipdp_c <=> dmpp_c<br>h2mb4p_c + nadh_c + h_c --> ipdp_c + nad_c + h2o_c<br>nad_c + 3c2hmp_c --> nadh_c + h_c + 3c4mop_c<br>3c2hmp_c <=> h2o_c + 2ippm_c<br>h2o_c + 2ippm_c <=> 3c3hmp_c<br>h2o_c + accoa_c + 3mob_c --> coa_c + 3c3hmp_c + h_c<br>isetac_p + atp_c + h2o_c --> pi_c + adp_c + isetac_c + h_c<br>isetac_e <=> isetac_p<br>atp_c + kdo2lipid4_c + h2o_c --> adp_c + pi_c + kdo2lipid4_p + h_c<br>atp_c + kdo2lipid4_p + h2o_c --> adp_c + pi_c + h_c + kdo2lipid4_e<br>nadp_c + 23dhmb_c <=> nadph_c + alac__S_c + h_c<br>nadph_c + 2ahbut_c + h_c <=> 23dhmp_c + nadp_c<br>acACP_c + malACP_c + h_c --> actACP_c + ACP_c + co2_c<br>accoa_c + malACP_c + h_c --> actACP_c + co2_c + coa_c<br>ctp_c + kdo_c --> ppi_c + ckdo_c<br>h2o_c + kdo8p_c --> pi_c + kdo_c<br>ara5p_c + h2o_c + pep_c --> kdo8p_c + pi_c<br>h_c + 3dhgulnp_c --> co2_c + xu5p__L_c<br>atp_c + k_p + h2o_c --> adp_c + k_c + pi_c + h_c<br>h_p + k_p --> k_c + h_c<br>h_p + k_c --> k_p + h_c<br>k_e <=> k_p<br>lac__L_c + q8_c --> pyr_c + q8h2_c<br>lac__L_c + mqn8_c --> pyr_c + mql8_c<br>h_p + lac__L_p <=> lac__L_c + h_c<br>lac__L_e <=> lac__L_p<br>uLa4n_p + colipa_p --> acolipa_p + udcpp_p<br>h2o_c + lcts_c --> glc__D_c + gal_c<br>h2o_p + lcts_p --> glc__D_p + gal_p<br>LalaDgluMdap_c + h2o_c --> LalaDglu_c + 26dap__M_c<br>LalaDglu_e <=> LalaDglu_p<br>h_p + LalaDglu_p --> h_c + LalaDglu_c<br>LalaLglu_e <=> LalaLglu_p<br>h_p + LalaLglu_p --> LalaLglu_c + h_c<br>nadh_c + h_c + mthgxl_c --> nad_c + lald__D_c<br>nadph_c + mthgxl_c + h_c --> lald__L_c + nadp_c<br>LalaLglu_c + h2o_c --> glu__L_c + ala__L_c<br>lald__L_c + nad_c + h2o_c --> lac__L_c + nadh_c + 2.0 h_c<br>h_c + lald__D_c + nadh_c <=> nad_c + 12ppd__R_c<br>lald__L_c + h_c + nadh_c <=> 12ppd__S_c + nad_c<br>h_p + lcts_c --> lcts_p + h_c<br>lcts_e <=> lcts_p<br>h_p + lcts_p <=> lcts_c + h_c<br>nad_c + lac__D_c <=> pyr_c + h_c + nadh_c<br>q8_c + lac__D_c --> q8h2_c + pyr_c<br>glu__L_c + 4mop_c --> akg_c + leu__L_c<br>atp_c + leu__L_c + trnaleu_c --> ppi_c + leutrna_c + amp_c<br>atp_c + h2o_c + leu__L_p --> adp_c + leu__L_c + pi_c + h_c<br>h_p + leu__L_p <=> leu__L_c + h_c<br>leu__L_e <=> leu__L_p<br>mthgxl_c + gthrd_c --> lgt__S_c<br>atp_c + lipa_cold_c + h2o_c --> adp_c + pi_c + lipa_cold_p + h_c<br>colipa_e + h_e + hdca_e --> h2o_e + hacolipa_e<br>lipa_e + h_e + hdca_e --> halipa_e + h2o_e<br>lipoamp_c --> amp_c + lipopb_c<br>atp_c + lipoate_c --> lipoamp_c + ppi_c<br>atp_c + h2o_c + lipa_c --> adp_c + pi_c + lipa_p + h_c<br>lipa_p + h2o_c + atp_c --> lipa_e + adp_c + pi_c + h_c<br>ocACP_c + h_c --> ACP_c + octapb_c<br>2.0 amet_c + octapb_c + nad_c + h_c + 4fe4s_c --> 2.0 fe2_c + 2.0 met__L_c + lipopb_c + 2.0 dad_5_c + nadh_c + 2fe2s_c<br>lipoate_p + h_p --> lipoate_c + h_c<br>lipoate_e <=> lipoate_p<br>u23ga_c + lipidX_c --> udp_c + lipidAds_c + h_c<br>h2o_p + 1ddecg3p_p --> ddca_p + glyc3p_p + h_p<br>h2o_p + 1tdecg3p_p --> h_p + ttdca_p + glyc3p_p<br>h2o_p + 1tdec7eg3p_p --> h_p + ttdcea_p + glyc3p_p<br>h2o_p + 1hdecg3p_p --> h_p + hdca_p + glyc3p_p<br>h2o_p + 1hdec9eg3p_p --> h_p + glyc3p_p + hdcea_p<br>h2o_p + 1odecg3p_p --> h_p + ocdca_p + glyc3p_p<br>h2o_p + 1odec11eg3p_p --> glyc3p_p + ocdcea_p + h_p<br>h2o_p + 1agpe120_p --> ddca_p + g3pe_p + h_p<br>h2o_p + 1agpe140_p --> ttdca_p + g3pe_p + h_p<br>h2o_p + 1agpe141_p --> h_p + ttdcea_p + g3pe_p<br>h2o_p + 1agpe160_p --> h_p + hdca_p + g3pe_p<br>h2o_p + 1agpe161_p --> h_p + g3pe_p + hdcea_p<br>h2o_p + 1agpe180_p --> h_p + g3pe_p + ocdca_p<br>h2o_p + 1agpe181_p --> h_p + g3pe_p + ocdcea_p<br>h2o_p + 1agpg120_p --> ddca_p + g3pg_p + h_p<br>h2o_p + 1agpg140_p --> h_p + ttdca_p + g3pg_p<br>h2o_p + 1agpg141_p --> h_p + ttdcea_p + g3pg_p<br>h2o_p + 1agpg160_p --> h_p + hdca_p + g3pg_p<br>h2o_p + 1agpg161_p --> g3pg_p + hdcea_p + h_p<br>1agpg180_p + h2o_p --> h_p + g3pg_p + ocdca_p<br>h2o_p + 1agpg181_p --> h_p + g3pg_p + ocdcea_p<br>h2o_c + 2ddecg3p_c --> 2.0 h_c + ddca_c + glyc3p_c<br>h2o_c + 2tdecg3p_c --> ttdca_c + 2.0 h_c + glyc3p_c<br>2tdec7eg3p_c + h2o_c --> 2.0 h_c + ttdcea_c + glyc3p_c<br>2hdecg3p_c + h2o_c --> glyc3p_c + 2.0 h_c + hdca_c<br>h2o_c + 2hdec9eg3p_c --> hdcea_c + 2.0 h_c + glyc3p_c<br>h2o_c + 2odecg3p_c --> ocdca_c + glyc3p_c + 2.0 h_c<br>2odec11eg3p_c + h2o_c --> ocdcea_c + 2.0 h_c + glyc3p_c<br>pg120_c + 2agpe120_c --> apg120_c + g3pe_c<br>pg140_c + 2agpe140_c --> g3pe_c + apg140_c<br>pg141_c + 2agpe141_c --> apg141_c + g3pe_c<br>pg160_c + 2agpe160_c --> apg160_c + g3pe_c<br>2agpe161_c + pg161_c --> g3pe_c + apg161_c<br>pg180_c + 2agpe180_c --> apg180_c + g3pe_c<br>pg181_c + 2agpe181_c --> g3pe_c + apg181_c<br>pg120_c + 2agpg120_c --> g3pg_c + apg120_c<br>pg140_c + 2agpg140_c --> apg140_c + g3pg_c<br>pg141_c + 2agpg141_c --> apg141_c + g3pg_c<br>pg160_c + 2agpg160_c --> apg160_c + g3pg_c<br>2agpg161_c + pg161_c --> apg161_c + g3pg_c<br>pg180_c + 2agpg180_c --> apg180_c + g3pg_c<br>2agpg181_c + pg181_c --> apg181_c + g3pg_c<br>2agpe120_c + h2o_c --> g3pe_c + h_c + ddca_c<br>h2o_c + 2agpe140_c --> g3pe_c + h_c + ttdca_c<br>h2o_c + 2agpe141_c --> g3pe_c + h_c + ttdcea_c<br>2agpe160_c + h2o_c --> g3pe_c + h_c + hdca_c<br>2agpe161_c + h2o_c --> g3pe_c + hdcea_c + h_c<br>h2o_c + 2agpe180_c --> ocdca_c + g3pe_c + h_c<br>h2o_c + 2agpe181_c --> g3pe_c + ocdcea_c + h_c<br>h2o_c + 2agpg120_c --> ddca_c + g3pg_c + h_c<br>h2o_c + 2agpg140_c --> ttdca_c + h_c + g3pg_c<br>2agpg141_c + h2o_c --> h_c + ttdcea_c + g3pg_c<br>2agpg160_c + h2o_c --> hdca_c + g3pg_c + h_c<br>2agpg161_c + h2o_c --> hdcea_c + g3pg_c + h_c<br>h2o_c + 2agpg180_c --> ocdca_c + h_c + g3pg_c<br>2agpg181_c + h2o_c --> ocdcea_c + h_c + g3pg_c<br>ser__L_c + nadp_c <=> nadph_c + h_c + 2amsa_c<br>lys__L_c + h_c --> 15dap_c + co2_c<br>atp_c + trnalys_c + lys__L_c --> ppi_c + amp_c + lystrna_c<br>atp_c + lys__L_p + h2o_c --> adp_c + lys__L_c + pi_c + h_c<br>h_p + lys__L_p --> lys__L_c + h_c<br>h_p + lys__L_c --> h_c + lys__L_p<br>lys__L_e <=> lys__L_p<br>lyx__L_c --> xylu__L_c<br>h_p + lyx__L_p --> lyx__L_c + h_c<br>lyx__L_e <=> lyx__L_p<br>nad_c + mnl1p_c <=> f6p_c + nadh_c + h_c<br>malACP_c + h_c --> acACP_c + co2_c<br>malcoa_c + amet_c --> malcoame_c + ahcys_c<br>nad_c + mal__D_c --> pyr_c + nadh_c + co2_c<br>2.0 h_p + mal__D_p --> mal__D_c + 2.0 h_c<br>mal__D_e <=> mal__D_p<br>h2o_c + accoa_c + glx_c --> coa_c + mal__L_c + h_c<br>accoa_c + malt_c <=> coa_c + acmalt_c<br>atp_c + malthx_p + h2o_c --> adp_c + pi_c + malthx_c + h_c<br>malthx_e --> malthx_p<br>atp_c + h2o_c + maltpt_p --> adp_c + pi_c + h_c + maltpt_c<br>maltpt_e --> maltpt_p<br>atp_c + h2o_c + malttr_p --> adp_c + pi_c + malttr_c + h_c<br>malttr_e --> malttr_p<br>atp_c + maltttr_p + h2o_c --> adp_c + pi_c + maltttr_c + h_c<br>maltttr_e --> maltttr_p<br>atp_c + h2o_c + malt_p --> adp_c + pi_c + malt_c + h_c<br>pep_c + malt_p --> malt6p_c + pyr_c<br>malt_e --> malt_p<br>2.0 h_p + mal__L_p --> mal__L_c + 2.0 h_c<br>3.0 h_p + mal__L_p --> mal__L_c + 3.0 h_c<br>mal__L_c + h_p --> h_c + mal__L_p<br>mal__L_e <=> mal__L_p<br>gdp_c + man1p_c + h_c --> gdpmann_c + pi_c<br>man6p_c <=> f6p_c<br>man6p_p + 2.0 pi_c --> 2.0 pi_p + man6p_c<br>man6p_e <=> man6p_p<br>nad_c + mana_c <=> h_c + fruur_c + nadh_c<br>manglyc_p + pep_c --> man6pglyc_c + pyr_c<br>manglyc_e <=> manglyc_p<br>man6pglyc_c + h2o_c --> glyc__R_c + man6p_c<br>man_p + pep_c --> pyr_c + man6p_c<br>man_e <=> man_p<br>2mcit_c --> 2mcacn_c + h2o_c<br>micit_c <=> pyr_c + succ_c<br>h2o_c + ppcoa_c + oaa_c --> coa_c + 2mcit_c + h_c<br>malcoa_c + ACP_c <=> coa_c + malACP_c<br>cyan_c + mercppyr_c --> pyr_c + tcynt_c + h_c<br>murein5p5p_p --> murein5px4p_p + ala__D_p<br>murein5p5p_p --> murein5px3p_p + alaala_p<br>murein5p5p5p_p --> 2.0 ala__D_p + murein5px4px4p_p<br>h2o_p + murein5px4p_p --> ala__D_p + murein4px4p_p<br>h2o_p + murein5px4px4p_p --> ala__D_p + murein4px4px4p_p<br>h2o_p + murein5p5p_p --> ala__D_p + murein5p4p_p<br>h2o_p + murein5p4p_p --> ala__D_p + murein4p4p_p<br>h2o_p + murein5p3p_p --> murein4p3p_p + ala__D_p<br>h2o_p + murein4px4p_p --> murein4p4p_p<br>h2o_p + murein3px4p_p --> murein4p3p_p<br>h2o_p + murein5px4p_p --> murein5p4p_p<br>h2o_p + murein4px4px4p_p --> murein4px4p4p_p<br>mal__L_c + nad_c <=> oaa_c + h_c + nadh_c<br>mal__L_c + q8_c --> oaa_c + q8h2_c<br>mal__L_c + mqn8_c --> oaa_c + mql8_c<br>nad_c + mal__L_c --> pyr_c + nadh_c + co2_c<br>mal__L_c + nadp_c --> nadph_c + pyr_c + co2_c<br>2mecdp_c + 2.0 flxr_c + h_c --> h2mb4p_c + 2.0 flxso_c + h2o_c<br>2p4c2me_c --> 2mecdp_c + cmp_c<br>melib_p + h_p --> melib_c + h_c<br>melib_c + h_p --> melib_p + h_c<br>melib_e <=> melib_p<br>meoh_e <=> meoh_p<br>meoh_p <=> meoh_c<br>ctp_c + h_c + 2me4p_c --> 4c2me_c + ppi_c<br>atp_c + met__L_c + h2o_c --> ppi_c + pi_c + amet_c<br>atp_c + met__D_p + h2o_c --> adp_c + pi_c + h_c + met__D_c<br>met__D_e <=> met__D_p<br>met__L_c + h2o2_c --> h2o_c + metsox_S__L_c<br>met__L_c + h2o2_c --> metsox_R__L_c + h2o_c<br>5mthf_c + hcys__L_c --> thf_c + met__L_c + h_c<br>atp_c + metsox_S__L_p + h2o_c --> adp_c + pi_c + metsox_S__L_c + h_c<br>metsox_S__L_e <=> metsox_S__L_p<br>atp_c + h2o_c + metsox_R__L_p --> adp_c + pi_c + metsox_R__L_c + h_c<br>metsox_R__L_e <=> metsox_R__L_p<br>trdrd_c + metsox_S__L_c --> met__L_c + h2o_c + trdox_c<br>trdrd_c + metsox_R__L_c --> met__L_c + h2o_c + trdox_c<br>atp_c + met__L_c + trnamet_c --> mettrna_c + ppi_c + amp_c<br>atp_c + met__L_p + h2o_c --> adp_c + met__L_c + pi_c + h_c<br>met__L_e <=> met__L_p<br>2.0 h_c + mg2_p <=> 2.0 h_p + mg2_c<br>mg2_e <=> mg2_p<br>mg2_p --> mg2_c<br>atp_c + mg2_p + h2o_c --> adp_c + mg2_c + pi_c + h_c<br>dhap_c --> pi_c + mthgxl_c<br>h2o_c + mi1p__D_c --> inost_c + pi_c<br>2mcacn_c + h2o_c <=> micit_c<br>mincyc_e <=> mincyc_p<br>h_p + mincyc_p --> mincyc_e + h_c<br>minohp_e --> minohp_p<br>h2o_p + murein5px4p_p --> alaala_p + murein3px4p_p<br>h2o_p + murein4p4p_p --> murein4p3p_p + ala__D_p<br>h2o_p + murein5p5p_p --> murein5p3p_p + alaala_p<br>murein4p3p_p + h2o_p --> ala__D_p + murein3p3p_p<br>h2o_p + murein5px3p_p --> murein3px3p_p + alaala_p<br>h2o_p + murein3px3p_p --> murein3p3p_p<br>h2o_p + murein5px3p_p --> murein5p3p_p<br>malttr_c + h2o_c --> glc__D_c + malt_c<br>h2o_c + maltttr_c --> malttr_c + glc__D_c<br>h2o_c + maltpt_c --> glc__D_c + maltttr_c<br>malthx_c + h2o_c --> maltpt_c + glc__D_c<br>h2o_c + malthp_c --> malthx_c + glc__D_c<br>murein4p4p_p --> 2.0 anhgm4p_p<br>murein4p3p_p --> anhgm3p_p + anhgm4p_p<br>murein3p3p_p --> 2.0 anhgm3p_p<br>murein4px4p4p_p --> anhgm4p_p + murein4px4p_p<br>pi_c + maltpt_c <=> g1p_c + maltttr_c<br>malthx_c + pi_c <=> g1p_c + maltpt_c<br>pi_c + malthp_c <=> malthx_c + g1p_c<br>mmcoa__S_c + h_c --> ppcoa_c + co2_c<br>h_p + mmet_p --> mmet_c + h_c<br>mmet_e <=> mmet_p<br>succoa_c --> mmcoa__S_c<br>h_p + mn2_c --> mn2_p + h_c<br>mn2_p --> mn2_c<br>h2o_c + man6p_c --> pi_c + man_c<br>mnl_p + pep_c --> mnl1p_c + pyr_c<br>mnl_e <=> mnl_p<br>mana_c --> h2o_c + 2ddglcn_c<br>mn2_p + h_p --> h_c + mn2_c<br>mn2_e <=> mn2_p<br>moadamp_c + iscssh_c + nadh_c --> amp_c + moadcosh_c + nad_c + iscs_c<br>ckdo_c + lipidA_c --> kdolipid4_c + cmp_c + h_c<br>kdolipid4_c + ckdo_c --> kdo2lipid4_c + cmp_c + h_c<br>phphhlipa_c + ckdo_c --> kphphhlipa_c + h_c + cmp_c<br>atp_c + h2o_c + mobd_p --> adp_c + pi_c + mobd_c + h_c<br>mobd_e <=> mobd_p<br>ctp_c + h_c + moco_c --> ppi_c + mococdp_c<br>mobd_c + mptamp_c + 2.0 h_c --> cu2_c + amp_c + h2o_c + moco_c<br>gtp_c + moco_c + h_c --> ppi_c + mocogdp_c<br>h2o_c + mlthf_c + 3mob_c --> thf_c + 2dhp_c<br>mal__L_c + o2_c <=> oaa_c + h2o2_c<br>atp_c + mpt_c + h_c --> ppi_c + mptamp_c<br>2.0 uaagmda_c --> murein5p5p_p + 2.0 h_c + 2.0 udcpdp_c<br>uaagmda_c + murein5p5p_p --> murein5p5p5p_p + h_c + udcpdp_c<br>cu2_c + 2.0 moadcosh_c + cpmp_c --> 2.0 moadcoo_c + mpt_c + 5.0 h_c<br>atp_c + moadcoo_c + h_c --> moadamp_c + ppi_c<br>nadph_c + msa_c + h_c --> nadp_c + 3hpp_c<br>atp_c + mso3_p + h2o_c --> adp_c + pi_c + mso3_c + h_c<br>mso3_e <=> mso3_p<br>5mta_c + h2o_c --> ade_c + 5mtr_c<br>h2o_c + methf_c <=> 10fthf_c + h_c<br>nadp_c + mlthf_c <=> nadph_c + methf_c<br>2.0 h_c + nadh_c + mlthf_c --> nad_c + 5mthf_c<br>h2o_c + mdhdhf_c --> mththf_c<br>o2_c + h2o_c + Nmtrp_c --> fald_c + trp__L_c + h2o2_c<br>n2o_e <=> n2o_p<br>n2o_p <=> n2o_c<br>h2o_c + acg5sa_c --> ac_c + glu5sa_c<br>nac_e <=> nac_p<br>nac_p --> nac_c<br>h2o_c + nad_c --> nmn_c + amp_c + 2.0 h_c<br>mqn8_c + h_c + nadh_c --> nad_c + mql8_c<br>q8_c + nadh_c + 4.0 h_c --> 3.0 h_p + nad_c + q8h2_c<br>mqn8_c + 4.0 h_c + nadh_c --> 3.0 h_p + nad_c + mql8_c<br>nadh_c + 2dmmq8_c + 4.0 h_c --> 2dmmql8_c + 3.0 h_p + nad_c<br>nadh_c + q8_c + h_c --> nad_c + q8h2_c<br>nadh_c + 2dmmq8_c + h_c --> 2dmmql8_c + nad_c<br>atp_c + nad_c --> adp_c + nadp_c + h_c<br>nad_c + h2o_c --> ncam_c + h_c + adprib_c<br>nadph_c + h_c + q8_c --> q8h2_c + nadp_c<br>nadph_c + mqn8_c + h_c --> nadp_c + mql8_c<br>nadph_c + 2dmmq8_c + h_c --> 2dmmql8_c + nadp_c<br>h2o_c + nadp_c --> nad_c + pi_c<br>atp_c + dnad_c + nh4_c --> amp_c + ppi_c + nad_c + h_c<br>nadph_c + nad_c --> nadp_c + nadh_c<br>atp_c + nac_c + prpp_c + h2o_c --> adp_c + ppi_c + pi_c + nicrnt_c<br>3.0 h_p + 2.0 na1_c --> 3.0 h_c + 2.0 na1_p<br>2.0 h_p + na1_c --> 2.0 h_c + na1_p<br>h_p + na1_c --> h_c + na1_p<br>na1_e <=> na1_p<br>atp_c + gdp_c <=> adp_c + gtp_c<br>udp_c + atp_c <=> adp_c + utp_c<br>atp_c + cdp_c <=> ctp_c + adp_c<br>atp_c + dtdp_c <=> dttp_c + adp_c<br>atp_c + dgdp_c <=> adp_c + dgtp_c<br>atp_c + dudp_c <=> adp_c + dutp_c<br>atp_c + dcdp_c <=> adp_c + dctp_c<br>atp_c + dadp_c <=> adp_c + datp_c<br>nh4_e <=> nh4_p<br>nh4_p <=> nh4_c<br>2.0 no_c + h_c + nadh_c --> n2o_c + nad_c + h2o_c<br>atp_c + ni2_c + h2o_c --> adp_c + pi_c + ni2_p + h_c<br>h_p + ni2_c --> ni2_p + h_c<br>ni2_e <=> ni2_p<br>ni2_p --> ni2_c<br>atp_c + ni2_p + h2o_c --> adp_c + ni2_c + pi_c + h_c<br>atp_c + h_c + nmn_c --> nad_c + ppi_c<br>nmn_c + h2o_c --> nicrnt_c + nh4_c<br>nmn_c + h2o_c --> ncam_c + h_c + r5p_c<br>nmn_p --> nmn_c<br>h2o_c + nmn_p --> ncam_c + h_c + r5p_c<br>nmn_e <=> nmn_p<br>ncam_c + h2o_c --> nac_c + nh4_c<br>atp_c + h_c + nicrnt_c <=> ppi_c + dnad_c<br>dmbzid_c + nicrnt_c --> nac_c + h_c + 5prdmbz_c<br>quln_c + prpp_c + 2.0 h_c --> ppi_c + nicrnt_c + co2_c<br>h_p + no2_p <=> no2_c + h_c<br>no2_e <=> no2_p<br>q8h2_c + no3_p --> h2o_p + no2_p + q8_c<br>no3_c + q8h2_c + 2.0 h_c --> no2_c + 2.0 h_p + q8_c + h2o_c<br>no3_p + mql8_c --> h2o_p + no2_p + mqn8_c<br>no3_c + 2.0 h_c + mql8_c --> no2_c + 2.0 h_p + mqn8_c + h2o_c<br>no2_c + no3_p --> no2_p + no3_c<br>no3_e <=> no3_p<br>2.0 o2_c + 2.0 no_c + nadh_c --> 2.0 no3_c + nad_c + h_c<br>nadph_c + 2.0 o2_c + 2.0 no_c --> nadp_c + 2.0 no3_c + h_c<br>novbcn_e <=> novbcn_p<br>h_p + novbcn_p --> novbcn_e + h_c<br>no_e <=> no_p<br>no_p <=> no_c<br>h2o_c + dump_c --> duri_c + pi_c<br>h2o_c + xmp_c --> pi_c + xtsn_c<br>h2o_p + xmp_p --> pi_p + xtsn_p<br>h2o_c + imp_c --> ins_c + pi_c<br>h2o_p + imp_p --> pi_p + ins_p<br>h2o_c + dimp_c --> din_c + pi_c<br>h2o_p + dimp_p --> pi_p + din_p<br>h2o_p + dump_p --> duri_p + pi_p<br>ump_c + h2o_c --> uri_c + pi_c<br>h2o_p + ump_p --> uri_p + pi_p<br>dcmp_c + h2o_c --> dcyt_c + pi_c<br>h2o_p + dcmp_p --> pi_p + dcyt_p<br>h2o_c + cmp_c --> cytd_c + pi_c<br>h2o_p + cmp_p --> pi_p + cytd_p<br>h2o_c + dtmp_c --> pi_c + thymd_c<br>h2o_p + dtmp_p --> pi_p + thymd_p<br>h2o_c + damp_c --> dad_2_c + pi_c<br>h2o_p + damp_p --> dad_2_p + pi_p<br>h2o_c + amp_c --> adn_c + pi_c<br>h2o_p + amp_p --> pi_p + adn_p<br>dgmp_c + h2o_c --> dgsn_c + pi_c<br>h2o_p + dgmp_p --> dgsn_p + pi_p<br>gmp_c + h2o_c --> gsn_c + pi_c<br>h2o_p + gmp_p --> gsn_p + pi_p<br>atp_c + h2o_c --> adp_c + pi_c + h_c<br>h2o_c + itp_c --> h_c + pi_c + idp_c<br>h2o_c + ditp_c --> h_c + pi_c + didp_c<br>xtp_c + h2o_c --> xdp_c + pi_c + h_c<br>gtp_c + h2o_c --> gdp_c + pi_c + h_c<br>h2o_p + gtp_p --> h_p + pi_p + gdp_p<br>h2o_c + ctp_c --> cdp_c + h_c + pi_c<br>h2o_c + dgtp_c --> dgmp_c + ppi_c + h_c<br>h2o_c + ditp_c --> ppi_c + h_c + dimp_c<br>xtp_c + h2o_c --> ppi_c + h_c + xmp_c<br>h2o_c + gtp_c --> ppi_c + h_c + gmp_c<br>h2o_c + dctp_c --> dcmp_c + ppi_c + h_c<br>ctp_c + h2o_c --> ppi_c + h_c + cmp_c<br>h2o_c + datp_c --> ppi_c + damp_c + h_c<br>atp_c + h2o_c --> ppi_c + amp_c + h_c<br>dttp_c + h2o_c --> ppi_c + dtmp_c + h_c<br>h2o_c + utp_c --> ump_c + ppi_c + h_c<br>h2o_c + itp_c --> h_c + ppi_c + imp_c<br>h2o_c + dgtp_c --> dgsn_c + pppi_c<br>gtp_c + h2o_c --> gsn_c + pppi_c<br>no2_c + 5.0 h_c + 3.0 nadh_c --> nh4_c + 3.0 nad_c + 2.0 h2o_c<br>2.0 h_p + no2_p + 3.0 q8h2_c --> nh4_p + 3.0 q8_c + 2.0 h2o_p<br>2.0 h_p + no2_p + 3.0 mql8_c --> 3.0 mqn8_c + nh4_p + 2.0 h2o_p<br>atp_c + o16a4colipa_p + h2o_c --> adp_c + pi_c + o16a4colipa_e + h_c<br>o16a4und_p + colipa_p --> h_p + udcpdp_p + o16a4colipa_p<br>2.0 o16aund_p --> h_p + o16a2und_p + udcpdp_p<br>o16aund_p + o16a2und_p --> h_p + udcpdp_p + o16a3und_p<br>o16aund_p + o16a3und_p --> o16a4und_p + h_p + udcpdp_p<br>accoa_c + ragund_c --> coa_c + aragund_c<br>o16aund_c --> o16aund_p<br>garagund_c + udpgalfur_c --> udp_c + h_c + gfgaragund_c<br>aragund_c + udpg_c --> udp_c + h_c + garagund_c<br>udpg_c + gfgaragund_c --> udp_c + o16aund_c + h_c<br>o2s_e <=> o2s_p<br>o2_e <=> o2_p<br>o2_p <=> o2_c<br>oaa_c + h_c --> pyr_c + co2_c<br>coa_c + 2obut_c --> ppcoa_c + for_c<br>cbp_c + orn_c <=> citr__L_c + h_c + pi_c<br>ocdca_e --> ocdca_p<br>ocdcea_e --> ocdcea_p<br>octa_e <=> octa_p<br>5.0 ipdp_c + frdp_c --> octdp_c + 5.0 ppi_c<br>atp_c + octa_c + h_c --> amp_c + octapb_c + ppi_c<br>hgmeACP_c --> h2o_c + egmeACP_c<br>nadph_c + h_c + ogmeACP_c --> hgmeACP_c + nadp_c<br>malcoame_c + malACP_c + h_c --> coa_c + ogmeACP_c + co2_c<br>glu__L_c + ohpb_c <=> akg_c + phthr_c<br>2ohph_c + amet_c --> ahcys_c + 2omph_c + h_c<br>amet_c + 2ombzl_c --> 2ommbl_c + h_c + ahcys_c<br>h_c + 3c4mop_c --> 4mop_c + co2_c<br>0.5 o2_c + 2ommbl_c --> 2omhmbl_c<br>2.0 atp_c + 2ommbl_c + nad_c + 3.0 h2o_c --> 2.0 adp_c + 2.0 pi_c + 2omhmbl_c + 3.0 h_c + nadh_c<br>orot5p_c + h_c --> ump_c + co2_c<br>0.5 o2_c + 2omph_c --> 2ombzl_c<br>2.0 atp_c + nad_c + 3.0 h2o_c + 2omph_c --> 2.0 adp_c + 2ombzl_c + 2.0 pi_c + nadh_c + 3.0 h_c<br>h2o_c + op4en_c --> 4h2opntn_c<br>3ophb_c + h_c --> 2oph_c + co2_c<br>2oph_c + 0.5 o2_c --> 2ohph_c<br>2.0 atp_c + 2oph_c + nad_c + 3.0 h2o_c --> 2.0 adp_c + 2ohph_c + 2.0 pi_c + 3.0 h_c + nadh_c<br>hpmeACP_c --> h2o_c + epmeACP_c<br>nadph_c + opmeACP_c + h_c --> nadp_c + hpmeACP_c<br>gmeACP_c + malACP_c + h_c --> opmeACP_c + ACP_c + co2_c<br>orn_c + h_c --> ptrc_c + co2_c<br>atp_c + orn_p + h2o_c --> adp_c + orn_c + pi_c + h_c<br>orn_e <=> orn_p<br>2.0 h_p + orot_p --> orot_c + 2.0 h_c<br>orot_e <=> orot_p<br>ppi_c + orot5p_c <=> orot_c + prpp_c<br>oxur_c + pi_c --> oxam_c + cbp_c<br>h_c + oxalcoa_c --> forcoa_c + co2_c<br>nadp_c + 2oxpaccoa_c + 2.0 h2o_c --> 3oxdhscoa_c + nadph_c + 2.0 h_c<br>3oxdhscoa_c + coa_c --> accoa_c + 23dhacoa_c<br>1pyr5c_c + nad_c + 2.0 h2o_c --> glu__L_c + nadh_c + h_c<br>nadph_c + 1pyr5c_c + 2.0 h_c --> pro__L_c + nadp_c<br>atp_c + pa120_c + h2o_c --> adp_c + pi_c + pa120_p + h_c<br>atp_c + h2o_c + pa140_c --> adp_c + pi_c + pa140_p + h_c<br>atp_c + pa141_c + h2o_c --> adp_c + pi_c + h_c + pa141_p<br>atp_c + pa160_c + h2o_c --> adp_c + pa160_p + pi_c + h_c<br>atp_c + h2o_c + pa161_c --> adp_c + pi_c + h_c + pa161_p<br>atp_c + pa180_c + h2o_c --> adp_c + pi_c + pa180_p + h_c<br>atp_c + pa181_c + h2o_c --> adp_c + pa181_p + h_c + pi_c<br>pacald_p + h_p <=> pacald_c + h_c<br>pacald_e <=> pacald_p<br>phaccoa_c + nadph_c + o2_c + h_c --> h2o_c + nadp_c + rephaccoa_c<br>atp_c + coa_c + pac_c --> ppi_c + phaccoa_c + amp_c<br>atp_c + pant__R_c + ala_B_c --> amp_c + pnto__R_c + ppi_c + h_c<br>pa120_c + h2o_c --> 12dgr120_c + pi_c<br>h2o_p + pa120_p --> pi_p + 12dgr120_p<br>h2o_c + pa140_c --> 12dgr140_c + pi_c<br>h2o_p + pa140_p --> pi_p + 12dgr140_p<br>h2o_c + pa141_c --> 12dgr141_c + pi_c<br>h2o_p + pa141_p --> 12dgr141_p + pi_p<br>h2o_c + pa160_c --> pi_c + 12dgr160_c<br>h2o_p + pa160_p --> pi_p + 12dgr160_p<br>h2o_c + pa161_c --> 12dgr161_c + pi_c<br>h2o_p + pa161_p --> 12dgr161_p + pi_p<br>pa180_c + h2o_c --> pi_c + 12dgr180_c<br>h2o_p + pa180_p --> pi_p + 12dgr180_p<br>h2o_c + pa181_c --> 12dgr181_c + pi_c<br>h2o_p + pa181_p --> 12dgr181_p + pi_p<br>ugmda_c + udcpp_c --> ump_c + uagmda_c<br>trdrd_c + paps_c --> trdox_c + pap_c + 2.0 h_c + so3_c<br>grxrd_c + paps_c --> so3_c + pap_c + 2.0 h_c + grxox_c<br>h2o_c + camp_c --> amp_c + h_c<br>35cgmp_c + h2o_c --> gmp_c + h_c<br>coa_c + pyr_c + nad_c --> nadh_c + accoa_c + co2_c<br>pdx5p_c + nad_c --> pydx5p_c + nadh_c + h_c<br>pdx5p_c + o2_c --> pydx5p_c + h2o2_c<br>phthr_c + nad_c + dxyl5p_c --> pdx5p_c + pi_c + co2_c + 2.0 h2o_c + h_c + nadh_c<br>pdx5p_c + h2o_c --> pi_c + pydxn_c<br>atp_c + h2o_c + pe120_c --> adp_c + pe120_p + pi_c + h_c<br>atp_c + pe140_c + h2o_c --> adp_c + pi_c + pe140_p + h_c<br>atp_c + h2o_c + pe141_c --> adp_c + pi_c + pe141_p + h_c<br>atp_c + h2o_c + pe160_c --> adp_c + pi_c + h_c + pe160_p<br>atp_c + pe161_c + h2o_c --> adp_c + pi_c + pe161_p + h_c<br>atp_c + pe180_c + h2o_c --> adp_c + pe180_p + h_c + pi_c<br>atp_c + pe181_c + h2o_c --> adp_c + pi_c + pe181_p + h_c<br>h2o_p + peamn_p + o2_p --> pacald_p + nh4_p + h2o2_p<br>peamn_e <=> peamn_p<br>nad_c + 4per_c <=> ohpb_c + h_c + nadh_c<br>pe161_p + lipa_p --> enlipa_p + 12dgr161_p<br>pe181_p + lipa_p --> enlipa_p + 12dgr181_p<br>atp_c + f6p_c --> h_c + fdp_c + adp_c<br>atp_c + tag6p__D_c --> adp_c + tagdp__D_c + h_c<br>atp_c + s7p_c --> s17bp_c + h_c + adp_c<br>coa_c + pyr_c --> for_c + accoa_c<br>atp_c + pg120_c + h2o_c --> pi_c + adp_c + h_c + pg120_p<br>atp_c + pg140_c + h2o_c --> adp_c + pg140_p + pi_c + h_c<br>atp_c + pg141_c + h2o_c --> adp_c + pi_c + pg141_p + h_c<br>atp_c + h2o_c + pg160_c --> adp_c + pi_c + pg160_p + h_c<br>atp_c + pg161_c + h2o_c --> adp_c + pi_c + h_c + pg161_p<br>atp_c + pg180_c + h2o_c --> adp_c + pi_c + pg180_p + h_c<br>atp_c + pg181_c + h2o_c --> adp_c + pg181_p + pi_c + h_c<br>gam1p_c <=> gam6p_c<br>nad_c + 3pg_c --> 3php_c + h_c + nadh_c<br>g6p_c <=> f6p_c<br>atp_c + 3pg_c <=> 13dpg_c + adp_c<br>6pgl_c + h2o_c --> 6pgc_c + h_c<br>2pglyc_c + h2o_c --> glyclt_c + pi_c<br>2pg_c <=> 3pg_c<br>g1p_c <=> g6p_c<br>atp_c + pgp120_c + h2o_c --> adp_c + pi_c + pgp120_p + h_c<br>atp_c + h2o_c + pgp140_c --> adp_c + pi_c + pgp140_p + h_c<br>atp_c + pgp141_c + h2o_c --> adp_c + pi_c + h_c + pgp141_p<br>atp_c + pgp160_c + h2o_c --> adp_c + pgp160_p + pi_c + h_c<br>atp_c + pgp161_c + h2o_c --> adp_c + pi_c + pgp161_p + h_c<br>atp_c + h2o_c + pgp180_c --> adp_c + pi_c + pgp180_p + h_c<br>atp_c + pgp181_c + h2o_c --> adp_c + pi_c + pgp181_p + h_c<br>pgp120_c + h2o_c --> pi_c + pg120_c<br>h2o_p + pgp120_p --> pi_p + pg120_p<br>h2o_c + pgp140_c --> pg140_c + pi_c<br>h2o_p + pgp140_p --> pi_p + pg140_p<br>h2o_c + pgp141_c --> pg141_c + pi_c<br>h2o_p + pgp141_p --> pi_p + pg141_p<br>pgp160_c + h2o_c --> pg160_c + pi_c<br>h2o_p + pgp160_p --> pi_p + pg160_p<br>pgp161_c + h2o_c --> pi_c + pg161_c<br>h2o_p + pgp161_p --> pi_p + pg161_p<br>h2o_c + pgp180_c --> pg180_c + pi_c<br>h2o_p + pgp180_p --> pi_p + pg180_p<br>pgp181_c + h2o_c --> pi_c + pg181_c<br>h2o_p + pgp181_p --> pi_p + pg181_p<br>cdpdddecg_c + glyc3p_c --> pgp120_c + cmp_c + h_c<br>cdpdtdecg_c + glyc3p_c --> cmp_c + h_c + pgp140_c<br>cdpdtdec7eg_c + glyc3p_c --> cmp_c + pgp141_c + h_c<br>glyc3p_c + cdpdhdecg_c --> pgp160_c + cmp_c + h_c<br>cdpdhdec9eg_c + glyc3p_c --> cmp_c + h_c + pgp161_c<br>glyc3p_c + cdpdodecg_c --> cmp_c + h_c + pgp180_c<br>cdpdodec11eg_c + glyc3p_c --> pgp181_c + cmp_c + h_c<br>atp_c + h2o_c + pheme_c --> adp_c + pi_c + pheme_p + h_c<br>pheme_p --> pheme_e<br>akg_c + phe__L_c <=> glu__L_c + phpyr_c<br>trnaphe_c + phe__L_c + atp_c --> phetrna_c + ppi_c + amp_c<br>h_p + phe__L_p <=> h_c + phe__L_c<br>phe__L_e <=> phe__L_p<br>6.0 h2o_p + minohp_p --> inost_p + 6.0 pi_p<br>h_p + pi_p <=> pi_c + h_c<br>pi_e <=> pi_p<br>atp_c + pi_p + h2o_c --> adp_c + 2.0 pi_c + h_c<br>h2o_p + pa120_p --> ddca_p + 2ddecg3p_p<br>h2o_p + pa140_p --> ttdca_p + 2tdecg3p_p<br>h2o_p + pa141_p --> 2tdec7eg3p_p + ttdcea_p<br>h2o_p + pa160_p --> 2hdecg3p_p + hdca_p<br>h2o_p + pa161_p --> 2hdec9eg3p_p + hdcea_p<br>h2o_p + pa180_p --> ocdca_p + 2odecg3p_p<br>h2o_p + pa181_p --> ocdcea_p + 2odec11eg3p_p<br>h2o_p + pe120_p --> h_p + 2agpe120_p + ddca_p<br>h2o_p + pe140_p --> h_p + ttdca_p + 2agpe140_p<br>h2o_p + pe141_p --> h_p + ttdcea_p + 2agpe141_p<br>h2o_p + pe160_p --> h_p + hdca_p + 2agpe160_p<br>h2o_p + pe161_p --> h_p + 2agpe161_p + hdcea_p<br>h2o_p + pe180_p --> h_p + 2agpe180_p + ocdca_p<br>h2o_p + pe181_p --> ocdcea_p + 2agpe181_p + h_p<br>h2o_p + pg120_p --> ddca_p + 2agpg120_p + h_p<br>h2o_p + pg140_p --> h_p + 2agpg140_p + ttdca_p<br>h2o_p + pg141_p --> h_p + ttdcea_p + 2agpg141_p<br>h2o_p + pg160_p --> h_p + hdca_p + 2agpg160_p<br>h2o_p + pg161_p --> h_p + 2agpg161_p + hdcea_p<br>h2o_p + pg180_p --> h_p + 2agpg180_p + ocdca_p<br>h2o_p + pg181_p --> h_p + 2agpg181_p + ocdcea_p<br>h2o_p + pa120_p --> ddca_p + h_p + 1ddecg3p_p<br>h2o_p + pa140_p --> h_p + ttdca_p + 1tdecg3p_p<br>h2o_p + pa141_p --> h_p + 1tdec7eg3p_p + ttdcea_p<br>h2o_p + pa160_p --> h_p + hdca_p + 1hdecg3p_p<br>h2o_p + pa161_p --> h_p + 1hdec9eg3p_p + hdcea_p<br>h2o_p + pa180_p --> h_p + 1odecg3p_p + ocdca_p<br>h2o_p + pa181_p --> 1odec11eg3p_p + ocdcea_p + h_p<br>h2o_p + pe120_p --> h_p + 1agpe120_p + ddca_p<br>h2o_p + pe140_p --> 1agpe140_p + ttdca_p + h_p<br>h2o_p + pe141_p --> h_p + ttdcea_p + 1agpe141_p<br>h2o_p + pe160_p --> h_p + hdca_p + 1agpe160_p<br>h2o_p + pe161_p --> h_p + 1agpe161_p + hdcea_p<br>h2o_p + pe180_p --> h_p + ocdca_p + 1agpe180_p<br>h2o_p + pe181_p --> ocdcea_p + 1agpe181_p + h_p<br>h2o_p + pg120_p --> ddca_p + h_p + 1agpg120_p<br>h2o_p + pg140_p --> h_p + 1agpg140_p + ttdca_p<br>h2o_p + pg141_p --> h_p + ttdcea_p + 1agpg141_p<br>h2o_p + pg160_p --> h_p + hdca_p + 1agpg160_p<br>h2o_p + pg161_p --> h_p + 1agpg161_p + hdcea_p<br>h2o_p + pg180_p --> 1agpg180_p + h_p + ocdca_p<br>h2o_p + pg181_p --> ocdcea_p + 1agpg181_p + h_p<br>man1p_c <=> man6p_c<br>5aprbu_c + h2o_c --> 4r5au_c + pi_c<br>h2o_c + pmeACP_c --> meoh_c + pimACP_c<br>atp_c + 4ampm_c --> adp_c + 2mahmp_c<br>atp_c + pnto__R_c --> adp_c + 4ppan_c + h_c<br>pnto__R_p + na1_p --> pnto__R_c + na1_c<br>pnto__R_e <=> pnto__R_p<br>nadh_c + poaac_c --> nad_c + h2o_c + 3amac_c<br>coa_c + pyr_c + 2.0 flxso_c <=> 2.0 flxr_c + accoa_c + h_c + co2_c<br>pyr_c + q8_c + h2o_c --> ac_c + q8h2_c + co2_c<br>h2o_c + ppi_c --> 2.0 pi_c + h_c<br>h2o_c + pppi_c --> h_c + ppi_c + pi_c<br>ppap_c + adp_c <=> atp_c + ppa_c<br>ppal_e <=> ppal_p<br>ppal_c <=> ppal_p<br>ppa_p + na1_p --> ppa_c + na1_c<br>ppa_e <=> ppa_p<br>2.0 5aop_c --> 2.0 h2o_c + h_c + ppbng_c<br>co2_c + h2o_c + pep_c --> oaa_c + pi_c + h_c<br>4ppcys_c + h_c --> pan4p_c + co2_c<br>atp_c + oaa_c --> adp_c + pep_c + co2_c<br>ppcoa_c + succ_c --> succoa_c + ppa_c<br>ppgpp_c + h2o_c --> gdp_c + ppi_c<br>atp_c + ppi_c <=> adp_c + pppi_c<br>atp_c + pi_c <=> adp_c + ppi_c<br>r1p_c <=> r5p_c<br>2dr1p_c <=> 2dr5p_c<br>ctp_c + 4ppan_c + cys__L_c --> 4ppcys_c + ppi_c + cmp_c + h_c<br>nad_c + pphn_c --> co2_c + nadh_c + 34hpp_c<br>h_c + pphn_c --> co2_c + h2o_c + phpyr_c<br>1.5 o2_c + pppg9_c --> 3.0 h2o_c + ppp9_c<br>pppg9_c + 3.0 fum_c --> 3.0 succ_c + ppp9_c<br>pppn_c + o2_c + nadh_c + h_c --> cechddd_c + nad_c<br>h_p + pppn_p <=> h_c + pppn_c<br>pppn_e <=> pppn_p<br>atp_c + pyr_c + h2o_c --> amp_c + pi_c + pep_c + 2.0 h_c<br>h2o_p + ppt_p --> pi_p + h2_p<br>ppt_e <=> ppt_p<br>atp_c + gly_c + pram_c <=> adp_c + pi_c + gar_c + h_c<br>atp_c + fpram_c --> adp_c + pi_c + air_c + 2.0 h_c<br>pran_c --> 2cpr5p_c<br>h2o_c + prbamp_c --> prfp_c<br>atp_c + 5aizc_c + asp__L_c --> 25aics_c + pi_c + adp_c + h_c<br>h2o_c + prbatp_c --> prbamp_c + ppi_c + h_c<br>atp_c + fgam_c + gln__L_c + h2o_c --> adp_c + pi_c + glu__L_c + h_c + fpram_c<br>prfp_c <=> prlp_c<br>pro__L_c + fad_c --> 1pyr5c_c + h_c + fadh2_c<br>atp_c + progly_p + h2o_c --> adp_c + progly_c + pi_c + h_c<br>progly_e <=> progly_p<br>pro__L_c + atp_c + trnapro_c --> ppi_c + protrna_c + amp_c<br>pro__L_p + h2o_c + atp_c --> pro__L_c + pi_c + adp_c + h_c<br>pro__L_p + h_p <=> pro__L_c + h_c<br>pro__L_p + na1_p --> pro__L_c + na1_c<br>pro__L_e <=> pro__L_p<br>atp_c + r5p_c <=> amp_c + prpp_c + h_c<br>h_p + psclys_p --> h_c + psclys_c<br>psclys_e <=> psclys_p<br>skm5p_c + pep_c <=> 3psme_c + pi_c<br>h_c + ps120_c --> co2_c + pe120_c<br>ps140_c + h_c --> pe140_c + co2_c<br>h_c + ps141_c --> co2_c + pe141_c<br>ps160_c + h_c --> co2_c + pe160_c<br>ps161_c + h_c --> pe161_c + co2_c<br>ps180_c + h_c --> pe180_c + co2_c<br>ps181_c + h_c --> pe181_c + co2_c<br>glu__L_c + 3php_c --> pser__L_c + akg_c<br>pser__L_e <=> pser__L_p<br>pser__L_c + h2o_c --> pi_c + ser__L_c<br>h2o_p + pser__L_p --> pi_p + ser__L_p<br>cdpdddecg_c + ser__L_c --> cmp_c + h_c + ps120_c<br>cdpdtdecg_c + ser__L_c --> ps140_c + cmp_c + h_c<br>ser__L_c + cdpdtdec7eg_c --> cmp_c + h_c + ps141_c<br>ser__L_c + cdpdhdecg_c --> ps160_c + h_c + cmp_c<br>cdpdhdec9eg_c + ser__L_c --> cmp_c + ps161_c + h_c<br>cdpdodecg_c + ser__L_c --> ps180_c + cmp_c + h_c<br>cdpdodec11eg_c + ser__L_c --> h_c + cmp_c + ps181_c<br>ppcoa_c + pi_c --> ppap_c + coa_c<br>accoa_c + pi_c <=> actp_c + coa_c<br>h2o_p + thrp_p --> thr__L_p + pi_p<br>atp_c + pan4p_c + h_c --> dpcoa_c + ppi_c<br>orn_c + ptrc_p <=> orn_p + ptrc_c<br>akg_c + ptrc_c --> glu__L_c + 4abutn_c<br>atp_c + ptrc_p + h2o_c --> adp_c + h_c + ptrc_c + pi_c<br>h_p + ptrc_p --> ptrc_c + h_c<br>ptrc_e <=> ptrc_p<br>adn_c + pi_c <=> ade_c + r1p_c<br>dad_2_c + pi_c <=> ade_c + 2dr1p_c<br>gsn_c + pi_c <=> r1p_c + gua_c<br>dgsn_c + pi_c <=> gua_c + 2dr1p_c<br>ins_c + pi_c <=> r1p_c + hxan_c<br>din_c + pi_c <=> 2dr1p_c + hxan_c<br>pi_c + xtsn_c <=> r1p_c + xan_c<br>pyam5p_c + o2_c + h2o_c --> nh4_c + pydx5p_c + h2o2_c<br>atp_c + pydam_c --> adp_c + pyam5p_c + h_c<br>pydam_e <=> pydam_p<br>pydam_p --> pydam_c<br>atp_c + pydx_c --> adp_c + pydx5p_c + h_c<br>atp_c + pydxn_c --> pdx5p_c + h_c + adp_c<br>pydxn_e <=> pydxn_p<br>pydxn_p --> pydxn_c<br>h2o_c + pydx5p_c --> pydx_c + pi_c<br>pydx_e <=> pydx_p<br>pydx_p --> pydx_c<br>adp_c + pep_c + h_c --> atp_c + pyr_c<br>uri_c + pi_c <=> r1p_c + ura_c<br>ura_c + o2_c + nadh_c + h_c --> uracp_c + nad_c<br>h_p + pyr_p <=> pyr_c + h_c<br>pyr_e <=> pyr_p<br>q8h2_c + 2.0 o2_c --> 2.0 o2s_c + 2.0 h_c + q8_c<br>2.0 o2_c + mql8_c --> mqn8_c + 2.0 h_c + 2.0 o2s_c<br>quin_e <=> quin_p<br>quin_p <=> quin_c<br>nad_c + quin_c <=> 2.0 h_c + nadh_c + 3dhq_c<br>dhap_c + iasp_c --> quln_c + 2.0 h2o_c + pi_c<br>atp_c + r15bp_c --> adp_c + prpp_c<br>atp_c + r1p_c --> adp_c + r15bp_c + h_c<br>r5p_c + h2o_c --> rib__D_c + pi_c<br>h2o_p + r5p_p --> rib__D_p + pi_p<br>r5p_e <=> r5p_p<br>atp_c + ribflv_c --> adp_c + h_c + fmn_c<br>db4p_c + 4r5au_c --> 2.0 h2o_c + pi_c + dmlz_c<br>2.0 dmlz_c --> ribflv_c + 4r5au_c<br>atp_c + rib__D_c --> adp_c + r5p_c + h_c<br>atp_c + rbl__L_c --> adp_c + ru5p__L_c + h_c<br>ru5p__L_c <=> xu5p__D_c<br>rephaccoa_c <=> 2oxpaccoa_c<br>rfamp_e <=> rfamp_p<br>h_p + rfamp_p --> rfamp_e + h_c<br>dtdprmn_c + kphphhlipa_c --> icolipa_c + dtdp_c + h_c<br>rhcys_c --> dhptd_c + hcys__L_c<br>atp_c + rib__D_p + h2o_c --> adp_c + pi_c + rib__D_c + h_c<br>rib__D_e <=> rib__D_p<br>rmn_c <=> rml_c<br>atp_c + rml_c --> adp_c + rml1p_c + h_c<br>rmn_e <=> rmn_p<br>h_p + rmn_p --> rmn_c + h_c<br>rml1p_c <=> lald__L_c + dhap_c<br>trdrd_c + adp_c --> h2o_c + trdox_c + dadp_c<br>grxrd_c + adp_c --> h2o_c + grxox_c + dadp_c<br>gdp_c + trdrd_c --> h2o_c + dgdp_c + trdox_c<br>gdp_c + grxrd_c --> h2o_c + grxox_c + dgdp_c<br>trdrd_c + cdp_c --> h2o_c + trdox_c + dcdp_c<br>cdp_c + grxrd_c --> h2o_c + grxox_c + dcdp_c<br>udp_c + trdrd_c --> dudp_c + trdox_c + h2o_c<br>udp_c + grxrd_c --> dudp_c + grxox_c + h2o_c<br>atp_c + 2.0 flxr_c + 2.0 h_c --> 2.0 flxso_c + h2o_c + datp_c<br>gtp_c + 2.0 flxr_c + 2.0 h_c --> dgtp_c + 2.0 flxso_c + h2o_c<br>2.0 flxr_c + ctp_c + 2.0 h_c --> 2.0 flxso_c + h2o_c + dctp_c<br>2.0 flxr_c + utp_c + 2.0 h_c --> 2.0 flxso_c + h2o_c + dutp_c<br>ru5p__D_c <=> xu5p__D_c<br>r5p_c <=> ru5p__D_c<br>h2o_c + 5prdmbz_c --> rdmbzi_c + pi_c<br>sufsesh_c + 2fe1s_c + sufbcd_c + h2o_c + atp_c --> adp_c + sufbcd_2fe2s_c + pi_c + sufse_c + 5.0 h_c<br>atp_c + 2.0 fe2_c + h2o_c + 2.0 sufsesh_c + sufbcd_c + fadh2_c --> sufbcd_2fe2s_c + 2.0 sufse_c + adp_c + pi_c + fad_c + 7.0 h_c<br>atp_c + 2.0 fe2_c + sufbcd_2fe2s_c + h2o_c + 2.0 sufsesh_c + fadh2_c --> sufbcd_2fe2s2_c + 2.0 sufse_c + adp_c + pi_c + fad_c + 7.0 h_c<br>sufbcd_2fe2s_c + 4.0 h_c --> sufbcd_c + 2fe2s_c<br>2.0 h_c + fadh2_c + sufbcd_2fe2s2_c --> sufbcd_4fe4s_c + fad_c<br>sufbcd_4fe4s_c + 4.0 h_c --> sufbcd_c + 4fe4s_c<br>s7p_c --> gmhep7p_c<br>sucarg_c + 2.0 h2o_c + 2.0 h_c --> co2_c + sucorn_c + 2.0 nh4_c<br>atp_c + so4_c + gtp_c + h2o_c --> ppi_c + gdp_c + pi_c + aps_c<br>o2_c + sarcs_c + h2o_c --> gly_c + fald_c + h2o2_c<br>nad_c + sbt6p_c <=> f6p_c + nadh_c + h_c<br>sbt__D_p + pep_c --> sbt6p_c + pyr_c<br>sbt__D_e <=> sbt__D_p<br>sufse_c + cys__L_c --> sufsesh_c + ala__L_c<br>h2o_c + sl26da_c --> succ_c + 26dap_LL_c<br>akg_c + sl26da_c <=> glu__L_c + sl2a6o_c<br>selnp_c + sertrna_sec_c --> h_c + sectrna_c + pi_c<br>slnt_c + 2.0 h_c + 4.0 gthrd_c --> gthox_c + 3.0 h2o_c + dgslnt_c<br>nadph_c + dgslnt_c + h_c --> nadp_c + gslnt_c + gthrd_c<br>nadph_c + gslnt_c --> gthrd_c + nadp_c + seln_c<br>atp_c + h2o_c + seln_c --> selnp_c + amp_c + pi_c<br>sel_c + mql8_c --> h2o_c + slnt_c + mqn8_c<br>sel_e <=> sel_p<br>2.0 h_p + sel_p --> 2.0 h_c + sel_c<br>ichor_c + h_c + akg_c --> 2sephchc_c + co2_c<br>atp_c + h_c + ser__L_c <=> seramp_c + ppi_c<br>accoa_c + ser__L_c <=> coa_c + acser_c<br>ser__D_c --> pyr_c + nh4_c<br>ser__L_c --> nh4_c + pyr_c<br>atp_c + ser__L_c + trnaser_c --> ppi_c + amp_c + sertrna_c<br>atp_c + ser__L_c + trnasecys_c --> sertrna_sec_c + ppi_c + amp_c<br>h_p + ser__L_p <=> h_c + ser__L_c<br>ser__L_p + na1_p --> ser__L_c + na1_c<br>ser__L_e <=> ser__L_p<br>Sfglutth_c + h2o_c --> for_c + h_c + gthrd_c<br>h2o_c + sucglu_c --> glu__L_c + succ_c<br>sucgsa_c + nad_c + h2o_c --> sucglu_c + 2.0 h_c + nadh_c<br>2sephchc_c --> pyr_c + 2shchc_c<br>nad_c + dscl_c --> scl_c + nadh_c + h_c<br>fe2_c + scl_c --> 3.0 h_c + sheme_c<br>nadph_c + 3dhsk_c + h_c <=> skm_c + nadp_c<br>atp_c + skm_c --> skm5p_c + h_c + adp_c<br>suchms_c + cys__L_c --> cyst__L_c + succ_c + h_c<br>h_p + skm_p --> skm_c + h_c<br>skm_e <=> skm_p<br>slnt_e <=> slnt_p<br>h_p + slnt_p --> h_c + slnt_c<br>so2_e <=> so2_p<br>so2_p <=> so2_c<br>so3_e <=> so3_p<br>h_p + so4_p --> so4_c + h_c<br>so4_e <=> so4_p<br>akg_c + sucorn_c --> sucgsa_c + glu__L_c<br>accoa_c + spmd_c --> N1aspmd_c + coa_c + h_c<br>spmd_c + accoa_c --> coa_c + n8aspmd_c + h_c<br>atp_c + spmd_p + h2o_c --> adp_c + spmd_c + pi_c + h_c<br>spmd_c + h_p --> spmd_p + h_c<br>spmd_e <=> spmd_p<br>ametam_c + ptrc_c --> 5mta_c + spmd_c + h_c<br>2.0 h_c + 2.0 o2s_c --> o2_c + h2o2_c<br>2.0 h_p + 2.0 o2s_p --> h2o2_p + o2_p<br>sucsal_c + nad_c + h2o_c --> succ_c + nadh_c + 2.0 h_c<br>h2o_c + nadp_c + sucsal_c --> succ_c + nadph_c + 2.0 h_c<br>asp__L_p + succ_c --> succ_p + asp__L_c<br>atp_c + coa_c + sucbz_c --> ppi_c + amp_c + sbzcoa_c<br>2shchc_c --> sucbz_c + h2o_c<br>2.0 h_p + succ_p --> 2.0 h_c + succ_c<br>3.0 h_p + succ_p --> 3.0 h_c + succ_c<br>h_p + succ_c --> succ_p + h_c<br>succ_e <=> succ_p<br>q8_c + succ_c --> q8h2_c + fum_c<br>fum_p + succ_c --> fum_c + succ_p<br>succ_c + mal__L_p --> mal__L_c + succ_p<br>atp_c + coa_c + succ_c <=> adp_c + pi_c + succoa_c<br>sucr_e <=> sucr_p<br>succ_c + tartr__D_p --> tartr__D_c + succ_p<br>sucr_p + pep_c --> suc6p_c + pyr_c<br>atp_c + sulfac_p + h2o_c --> adp_c + sulfac_c + pi_c + h_c<br>sulfac_e <=> sulfac_p<br>3.0 nadph_c + 5.0 h_c + so3_c --> h2s_c + 3.0 nadp_c + 3.0 h2o_c<br>atp_c + h2o_c + so4_p --> so4_c + pi_c + adp_c + h_c<br>tdec2eACP_c <=> cdec3eACP_c<br>nad_c + altrn_c <=> tagur_c + h_c + nadh_c<br>g3p_c + s7p_c <=> e4p_c + f6p_c<br>tartr__L_c --> oaa_c + h2o_c<br>tartr__D_e <=> tartr__D_p<br>tartr__L_p + succ_c <=> tartr__L_c + succ_p<br>tartr__L_e <=> tartr__L_p<br>3.0 h_p + tartr__D_p --> 3.0 h_c + tartr__D_c<br>o2_c + taur_c + akg_c --> co2_c + aacald_c + succ_c + so3_c + h_c<br>atp_c + h2o_c + taur_p --> adp_c + pi_c + taur_c + h_c<br>taur_e <=> taur_p<br>tcynt_e <=> tcynt_p<br>thmpp_c + h2o_c --> thmmp_c + pi_c + h_c<br>accoa_c + dtdp4addg_c --> coa_c + h_c + dtdp4aaddg_c<br>glu__L_c + dtdp4d6dg_c --> akg_c + dtdp4addg_c<br>dtdp4d6dg_c --> dtdp4d6dm_c<br>nadph_c + dtdp4d6dm_c + h_c --> dtdprmn_c + nadp_c<br>dtdpglu_c --> h2o_c + dtdp4d6dg_c<br>atp_c + lipidAds_c --> adp_c + h_c + lipidA_c<br>dsbdrd_c + dsbcox_p --> dsbdox_c + dsbcrd_p<br>dsbdrd_c + dsbgox_p --> dsbgrd_p + dsbdox_c<br>tagdp__D_c <=> g3p_c + dhap_c<br>2.0 h_p + nadp_c + nadh_c --> nadph_c + nad_c + 2.0 h_c<br>succoa_c + h2o_c + thdp_c --> coa_c + sl2a6o_c<br>h2o_c + methf_c --> h_c + 5fthf_c<br>trdrd_c + h2o2_c --> 2.0 h2o_c + trdox_c<br>h_p + thymd_p --> h_c + thymd_c<br>h_p + thymd_p <=> h_c + thymd_c<br>thymd_e <=> thymd_p<br>atp_c + h2o_c + thm_p --> adp_c + thm_c + pi_c + h_c<br>thm_e <=> thm_p<br>athr__L_c --> gly_c + acald_c<br>thr__L_c --> acald_c + gly_c<br>nad_c + thr__L_c --> 2aobut_c + h_c + nadh_c<br>thr__L_c --> 2obut_c + nh4_c<br>thrp_e <=> thrp_p<br>h2o_c + phom_c --> thr__L_c + pi_c<br>atp_c + thr__L_c + trnathr_c --> ppi_c + thrtrna_c + amp_c<br>atp_c + thr__L_p + h2o_c --> adp_c + thr__L_c + pi_c + h_c<br>h_p + thr__L_c --> thr__L_p + h_c<br>thr__L_p + h_p <=> thr__L_c + h_c<br>thr__L_p + na1_p --> thr__L_c + na1_c<br>thr__L_e <=> thr__L_p<br>h_p + thym_c --> thym_p + h_c<br>thym_e <=> thym_p<br>atp_c + dhgly_c + nadph_c + iscssh_c + dxyl5p_c + h_c --> ppi_c + nadp_c + 2.0 h2o_c + co2_c + amp_c + iscs_c + 4mpetz_c<br>xu5p__D_c + r5p_c <=> g3p_c + s7p_c<br>e4p_c + xu5p__D_c <=> g3p_c + f6p_c<br>tmao_c + h_c + mql8_c --> tma_c + mqn8_c + h2o_c<br>h_p + tmao_p + mql8_c --> mqn8_c + h2o_p + tma_p<br>2dmmql8_c + tmao_c + h_c --> tma_c + h2o_c + 2dmmq8_c<br>2dmmql8_c + tmao_p + h_p --> h2o_p + 2dmmq8_c + tma_p<br>tmao_e <=> tmao_p<br>tma_e <=> tma_p<br>atp_c + thymd_c --> adp_c + dtmp_c + h_c<br>pi_c + thymd_c <=> thym_c + 2dr1p_c<br>mlthf_c + dump_c --> dhf_c + dtmp_c<br>atp_c + thm_c --> adp_c + thmmp_c + h_c<br>atp_c + thmmp_c --> adp_c + thmpp_c<br>2mahmp_c + 4mpetz_c + h_c --> ppi_c + thmmp_c<br>dhap_c <=> g3p_c<br>dpcoa_c + atp_c --> ade_c + 2tpr3dpcoa_c<br>nadph_c + h_c + trdox_c --> trdrd_c + nadp_c<br>tre6p_c + h2o_c --> glc__D_c + g6p_c<br>tre6p_c + h2o_c --> tre_c + pi_c<br>g6p_c + udpg_c --> udp_c + tre6p_c + h_c<br>h2o_c + tre_c --> 2.0 glc__D_c<br>h2o_p + tre_p --> 2.0 glc__D_p<br>pep_c + tre_p --> tre6p_c + pyr_c<br>tre_e <=> tre_p<br>h2o_c + trp__L_c <=> indole_c + pyr_c + nh4_c<br>3ig3p_c + ser__L_c --> g3p_c + h2o_c + trp__L_c<br>indole_c + ser__L_c --> h2o_c + trp__L_c<br>3ig3p_c --> indole_c + g3p_c<br>atp_c + trnatrp_c + trp__L_c --> trptrna_c + ppi_c + amp_c<br>trp__L_p + h_p <=> trp__L_c + h_c<br>trp__L_e <=> trp__L_p<br>2h3oppan_c + nadh_c + h_c <=> glyc__R_c + nad_c<br>atp_c + h2o_c + tsul_p --> adp_c + tsul_c + pi_c + h_c<br>tsul_e <=> tsul_p<br>ttdca_e --> ttdca_p<br>ttdcea_e --> ttdcea_p<br>ttrcyc_e <=> ttrcyc_p<br>h_p + ttrcyc_p --> ttrcyc_e + h_c<br>atp_c + tungs_p + h2o_c --> adp_c + tungs_c + h_c + pi_c<br>tungs_e <=> tungs_p<br>tym_e <=> tym_p<br>tyr__L_c + nadph_c + amet_c --> met__L_c + nadp_c + dhgly_c + 4crsol_c + h_c + dad_5_c<br>h2o_p + tym_p + o2_p --> 4hoxpacd_p + nh4_p + h2o2_p<br>h2o_p + tyrp_p --> tyr__L_p + pi_p<br>tyrp_e <=> tyrp_p<br>akg_c + tyr__L_c <=> glu__L_c + 34hpp_c<br>atp_c + trnatyr_c + tyr__L_c --> ppi_c + tyrtrna_c + amp_c<br>h_p + tyr__L_p <=> tyr__L_c + h_c<br>tyr__L_e <=> tyr__L_p<br>adp_c + thmpp_c + h_c --> athtp_c + pi_c<br>3hmrsACP_c + u3hga_c --> u23ga_c + ACP_c + h_c<br>atp_c + uamag_c + 26dap__M_c --> adp_c + pi_c + ugmd_c + h_c<br>h2o_p + udpacgal_p --> 2.0 h_p + acgal1p_p + ump_p<br>h2o_p + uacgam_p --> 2.0 h_p + acgam1p_p + ump_p<br>uacgam_e <=> uacgam_p<br>2.0 nad_c + h2o_c + uacmam_c --> 3.0 h_c + uacmamu_c + 2.0 nadh_c<br>uacgam_c <=> uacmam_c<br>uacgam_c + 3hmrsACP_c <=> ACP_c + u3aga_c<br>uacgam_c + pep_c --> pi_c + uaccg_c<br>acgam1p_c + h_c + utp_c --> uacgam_c + ppi_c<br>uacgam_c + uagmda_c --> udp_c + h_c + uaagmda_c<br>atp_c + glu__D_c + uama_c --> adp_c + pi_c + uamag_c + h_c<br>atp_c + ala__L_c + uamr_c --> adp_c + uama_c + pi_c + h_c<br>nadph_c + h_c + uaccg_c --> uamr_c + nadp_c<br>udcpdp_c + h2o_c --> udcpp_c + pi_c + h_c<br>8.0 ipdp_c + frdp_c --> udcpdp_c + 8.0 ppi_c<br>h2o_p + udcpdp_p --> h_p + pi_p + udcpp_p<br>udcpp_p --> udcpp_c<br>udpacgal_e <=> udpacgal_p<br>udpg_c <=> udpgal_c<br>udpgal_c --> udpgalfur_c<br>h2o_p + udpgal_p --> 2.0 h_p + ump_p + gal1p_p<br>udpgal_e <=> udpgal_p<br>udpg_c + 2.0 nad_c + h2o_c --> udpglcur_c + 2.0 nadh_c + 3.0 h_c<br>nad_c + udpglcur_c --> udpLa4o_c + nadh_c + co2_c<br>udpglcur_e <=> udpglcur_p<br>h2o_p + udpg_p --> 2.0 h_p + ump_p + g1p_p<br>udpg_e <=> udpg_p<br>glu__L_c + udpLa4o_c <=> akg_c + udpLa4n_c<br>h2o_p + udpglcur_p --> ump_p + glcur1p_p + 2.0 h_p<br>gal1p_c + udpg_c <=> udpgal_c + g1p_c<br>h2o_c + 2.0 h_c + urdglyc_c --> glx_c + co2_c + 2.0 nh4_c<br>atp_c + ugmd_c + alaala_c --> adp_c + pi_c + ugmda_c + h_c<br>h2o_c + u3aga_c --> ac_c + u3hga_c<br>udpLa4n_c + 10fthf_c --> thf_c + udpLa4fn_c + h_c<br>uLa4n_c --> uLa4n_p<br>atp_c + LalaDgluMdap_c + uamr_c --> adp_c + pi_c + ugmd_c + h_c<br>h2o_c + um4p_c --> ugmd_c + ala__D_c<br>atp_c + uamr_c + LalaDgluMdapDala_c --> adp_c + pi_c + h_c + um4p_c<br>atp_c + ump_c <=> udp_c + adp_c<br>ump_e <=> ump_p<br>uLa4fn_c + h2o_c --> for_c + uLa4n_c<br>udpLa4fn_c + udcpp_c --> udp_c + uLa4fn_c<br>uppg3_c + 2.0 amet_c --> h_c + dscl_c + 2.0 ahcys_c<br>hmbil_c --> uppg3_c + h2o_c<br>uppg3_c + 4.0 h_c --> 4.0 co2_c + cpppg3_c<br>ura_c + prpp_c --> ump_c + ppi_c<br>uracp_c + h2o_c --> cbm_c + h_c + poaac_c<br>h_p + ura_p --> ura_c + h_c<br>h_p + ura_p <=> ura_c + h_c<br>ura_e <=> ura_p<br>nad_c + urdglyc_c --> nadh_c + h_c + oxur_c<br>urea_e <=> urea_p<br>urea_p <=> urea_c<br>o2_c + urate_c + 2.0 h2o_c --> co2_c + alltn_c + h2o2_c<br>h2o_c + uri_c --> rib__D_c + ura_c<br>gtp_c + uri_c --> gdp_c + h_c + ump_c<br>h_p + uri_p --> uri_c + h_c<br>h_p + uri_p <=> uri_c + h_c<br>uri_e <=> uri_p<br>u23ga_c + h2o_c --> ump_c + lipidX_c + 2.0 h_c<br>akg_c + val__L_c <=> glu__L_c + 3mob_c<br>atp_c + trnaval_c + val__L_c --> valtrna_c + ppi_c + amp_c<br>atp_c + val__L_p + h2o_c --> adp_c + pi_c + val__L_c + h_c<br>h_p + val__L_p <=> val__L_c + h_c<br>val__L_e <=> val__L_p<br>ala__L_c + 3mob_c <=> pyr_c + val__L_c<br>tungs_c + mptamp_c + 2.0 h_c --> cu2_c + wco_c + h2o_c + amp_c<br>xu5p__L_c --> ru5p__L_c<br>xan_c + nad_c + h2o_c --> urate_c + nadh_c + h_c<br>h_p + xan_p --> h_c + xan_c<br>xan_e <=> xan_p<br>xan_p <=> xan_c<br>xmp_e <=> xmp_p<br>prpp_c + xan_c --> ppi_c + xmp_c<br>h2o_c + xtsn_c --> rib__D_c + xan_c<br>h_p + xtsn_p <=> xtsn_c + h_c<br>xtsn_e <=> xtsn_p<br>xyl__D_c <=> xylu__D_c<br>glc__D_c <=> fru_c<br>atp_c + xylu__D_c --> xu5p__D_c + h_c + adp_c<br>atp_c + xylu__L_c --> adp_c + xu5p__L_c + h_c<br>h_p + xylu__L_p --> h_c + xylu__L_c<br>xylu__L_e <=> xylu__L_p<br>atp_c + xyl__D_p + h2o_c --> adp_c + pi_c + h_c + xyl__D_c<br>xyl__D_p + h_p --> h_c + xyl__D_c<br>xyl__D_e <=> xyl__D_p<br>atp_c + zn2_c + h2o_c --> adp_c + pi_c + h_c + zn2_p<br>h_p + zn2_c --> h_c + zn2_p<br>zn2_p --> zn2_c<br>atp_c + h2o_c + zn2_p --> adp_c + pi_c + zn2_c + h_c<br>zn2_e <=> zn2_p</div></td>
    </tr>
    </table>



Step 2: Simulate a model
------------------------

The model can be simulated by executing
`~cameo.core.solver_based_model.SolverBasedModel.solve`.

.. code:: python

    solution = model.solve()

A quick overview of the solution can be obtained in form of a pandas
`DataFrame <http://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.html>`__
(all solution objects in cameo provide access to data frames through a
`data_frame` attribute).

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
          <td>0.000000</td>
        </tr>
        <tr>
          <th>DM_5drib_c</th>
          <td>0.000221</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>DM_aacald_c</th>
          <td>0.000000</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>DM_amob_c</th>
          <td>0.000002</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
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
          <td>0.000000</td>
        </tr>
      </tbody>
    </table>
    <p>2583 rows × 2 columns</p>
    </div>



The data frame is accessible through `solution.data_frame`.

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
          <td>0.000000</td>
        </tr>
        <tr>
          <th>DM_5drib_c</th>
          <td>0.000221</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>DM_aacald_c</th>
          <td>0.000000</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>DM_amob_c</th>
          <td>0.000002</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
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
          <td>0.000000</td>
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
          <td>0</td>
        </tr>
        <tr>
          <th>DM_5drib_c</th>
          <td>0.000221</td>
          <td>0</td>
        </tr>
        <tr>
          <th>DM_amob_c</th>
          <td>0.000002</td>
          <td>0</td>
        </tr>
        <tr>
          <th>DM_mththf_c</th>
          <td>0.000440</td>
          <td>0</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>USHD</th>
          <td>0.019113</td>
          <td>0</td>
        </tr>
        <tr>
          <th>VPAMTr</th>
          <td>0.415702</td>
          <td>0</td>
        </tr>
        <tr>
          <th>ZN2tpp</th>
          <td>0.000335</td>
          <td>0</td>
        </tr>
        <tr>
          <th>Zn2tex</th>
          <td>0.000335</td>
          <td>0</td>
        </tr>
      </tbody>
    </table>
    <p>439 rows × 2 columns</p>
    </div>



Step 3: Exploring a model
-------------------------

Objects—models, reactions, metabolites, genes—can easily be explored in
the Jupyter notebook, taking advantage of tab completion. For example,
place your cursor after the period in `model.reactions.` and press the
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
                    <td><strong>Stoichiometry</strong></td><td>atp_c + 3pg_c <=> 13dpg_c + adp_c</td>
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
                    <td><strong>Stoichiometry</strong></td><td>e4p_c + nad_c + h2o_c <=> 4per_c + nadh_c + 2.0 h_c</td>
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

In these cases you need to use the `model.reactions.get_by_id`.

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
            



Metabolites are accessible through `model.metabolites`. For example,
D-glucose in the cytosolic compartment.

.. code:: python

    model.metabolites.glc__D_c




.. parsed-literal::

    <Metabolite glc__D_c at 0x117fcd710>



And it is easy to find the associated reactions

.. code:: python

    model.metabolites.glc__D_c.reactions




.. parsed-literal::

    frozenset({<Reaction AMALT1 at 0x117d490b8>,
               <Reaction MLTG4 at 0x1149604e0>,
               <Reaction GLCabcpp at 0x114999940>,
               <Reaction LACZ at 0x114975d68>,
               <Reaction GALS3 at 0x1149a0978>,
               <Reaction AMALT4 at 0x117d49630>,
               <Reaction MLTG2 at 0x114960630>,
               <Reaction G6PP at 0x1149a0f98>,
               <Reaction AMALT3 at 0x117d49668>,
               <Reaction TREH at 0x1148faa90>,
               <Reaction GLCt2pp at 0x114999a90>,
               <Reaction MLTG1 at 0x1149606a0>,
               <Reaction XYLI2 at 0x1148ecac8>,
               <Reaction MLTG3 at 0x1149606d8>,
               <Reaction MLTG5 at 0x1149603c8>,
               <Reaction HEX1 at 0x114984f28>,
               <Reaction GLCATr at 0x114999f98>,
               <Reaction TRE6PH at 0x1148faba8>,
               <Reaction AMALT2 at 0x117d493c8>})



A list of the genes encoded in the model can be accessed via
`model.genes`.

.. code:: python

    model.genes[0:10]




.. parsed-literal::

    [<Gene b0241 at 0x118050128>,
     <Gene b1377 at 0x118050160>,
     <Gene b2215 at 0x118050198>,
     <Gene b0929 at 0x1180501d0>,
     <Gene b4034 at 0x118050208>,
     <Gene b4033 at 0x118050240>,
     <Gene b4035 at 0x118050278>,
     <Gene b4032 at 0x1180502b0>,
     <Gene b4036 at 0x1180502e8>,
     <Gene b4213 at 0x118050320>]



A few additional attributes have been added that are not available in a
`cobrapy <https://opencobra.github.io/cobrapy/>`__ model. For example,
exchange reactions that allow certain metabolites to enter or leave the
model can be accessed through `model.exchanges`.

.. code:: python

    model.exchanges[0:10]




.. parsed-literal::

    [<Reaction DM_4crsol_c at 0x116c7cd30>,
     <Reaction DM_5drib_c at 0x1124c2ac8>,
     <Reaction DM_aacald_c at 0x1124c2898>,
     <Reaction DM_amob_c at 0x1124c21d0>,
     <Reaction DM_mththf_c at 0x1124c2390>,
     <Reaction DM_oxam_c at 0x116cff588>,
     <Reaction EX_12ppd__R_e at 0x116cff978>,
     <Reaction EX_12ppd__S_e at 0x116cff7f0>,
     <Reaction EX_14glucan_e at 0x116cffa20>,
     <Reaction EX_15dap_e at 0x116cff278>]



Or, the current medium can be accessed through `model.medium`.

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
          <td>1000</td>
        </tr>
        <tr>
          <th>1</th>
          <td>EX_cbl1_e</td>
          <td>Cob(I)alamin exchange</td>
          <td>-0.01</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>2</th>
          <td>EX_cl_e</td>
          <td>Chloride exchange</td>
          <td>-1000.00</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>3</th>
          <td>EX_co2_e</td>
          <td>CO2 exchange</td>
          <td>-1000.00</td>
          <td>1000</td>
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
          <td>1000</td>
        </tr>
        <tr>
          <th>22</th>
          <td>EX_so4_e</td>
          <td>Sulfate exchange</td>
          <td>-1000.00</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>23</th>
          <td>EX_tungs_e</td>
          <td>Tungstate exchange</td>
          <td>-1000.00</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>24</th>
          <td>EX_zn2_e</td>
          <td>Zinc exchange</td>
          <td>-1000.00</td>
          <td>1000</td>
        </tr>
      </tbody>
    </table>
    <p>25 rows × 4 columns</p>
    </div>



It is also possible to get a list of essential reactions ...

.. code:: python

    model.essential_reactions()[0:10]




.. parsed-literal::

    [<Reaction DM_4crsol_c at 0x116c7cd30>,
     <Reaction DM_5drib_c at 0x1124c2ac8>,
     <Reaction DM_amob_c at 0x1124c21d0>,
     <Reaction DM_mththf_c at 0x1124c2390>,
     <Reaction BIOMASS_Ec_iJO1366_core_53p95M at 0x116cff2e8>,
     <Reaction EX_ca2_e at 0x116cfff60>,
     <Reaction EX_cl_e at 0x116d05198>,
     <Reaction EX_cobalt2_e at 0x116d05278>,
     <Reaction EX_cu2_e at 0x116d05470>,
     <Reaction EX_glc__D_e at 0x116d13128>]



... and essential genes.

.. code:: python

    model.essential_genes()[0:10]




.. parsed-literal::

    [<Gene b0784 at 0x1180c0048>,
     <Gene b3633 at 0x1180c0080>,
     <Gene b3959 at 0x1180640b8>,
     <Gene b0417 at 0x1180e8080>,
     <Gene b3993 at 0x1180e80b8>,
     <Gene b2818 at 0x1180640f0>,
     <Gene b3176 at 0x1180d0160>,
     <Gene b0182 at 0x1180b81d0>,
     <Gene b1261 at 0x1180e82e8>,
     <Gene b0134 at 0x1180c02e8>]


