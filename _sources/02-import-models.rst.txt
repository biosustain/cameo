
Import models
=============

Import models from files
------------------------

The function `~cameo.io.load_model` accepts a number of different
input formats.

1. `SBML <http://sbml.org/>`__ (Systems Biology Markup Language).
2. JSON
3. Pickle (pickled models)
4. Model identifiers (from the `BiGG Models <http://bigg.ucsd.edu>`__)

.. code:: ipython3

    less data/e_coli_core.xml

.. code:: ipython3

    from cameo import load_model
    model = load_model('data/e_coli_core.xml')

.. code:: ipython3

    model




.. raw:: html

    
            <table>
                <tr>
                    <td><strong>Name</strong></td>
                    <td>e_coli_core</td>
                </tr><tr>
                    <td><strong>Memory address</strong></td>
                    <td>0x01102bfdd8</td>
                </tr><tr>
                    <td><strong>Number of metabolites</strong></td>
                    <td>72</td>
                </tr><tr>
                    <td><strong>Number of reactions</strong></td>
                    <td>95</td>
                </tr><tr>
                    <td><strong>Objective expression</strong></td>
                    <td>-1.0*Biomass_Ecoli_core_w_GAM_reverse_1a29b + 1.0*Biomass_Ecoli_core_w_GAM</td>
                </tr><tr>
                    <td><strong>Compartments</strong></td>
                    <td>extracellular space, cytosol</td>
                </tr>
              </table>



Import models from the internet
-------------------------------

In the quick start `chapter <1-quick-start.ipynb>`__ we demonstrated how
to use `~cameo.io.load_model` to import a model by ID. But where did
the model come from? Cameo has currently access to two model
repositories on the internet, http://bigg.ucsd.edu and
http://darwin.di.uminho.pt/models.

.. code:: ipython3

    from cameo import models

.. code:: ipython3

    models.index_models_bigg()




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
          <th>bigg_id</th>
          <th>gene_count</th>
          <th>metabolite_count</th>
          <th>organism</th>
          <th>reaction_count</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>e_coli_core</td>
          <td>137</td>
          <td>72</td>
          <td>Escherichia coli str. K-12 substr. MG1655</td>
          <td>95</td>
        </tr>
        <tr>
          <th>1</th>
          <td>iAB_RBC_283</td>
          <td>346</td>
          <td>342</td>
          <td>Homo sapiens</td>
          <td>469</td>
        </tr>
        <tr>
          <th>2</th>
          <td>iAF1260</td>
          <td>1261</td>
          <td>1668</td>
          <td>Escherichia coli str. K-12 substr. MG1655</td>
          <td>2382</td>
        </tr>
        <tr>
          <th>3</th>
          <td>iAF1260b</td>
          <td>1261</td>
          <td>1668</td>
          <td>Escherichia coli str. K-12 substr. MG1655</td>
          <td>2388</td>
        </tr>
        <tr>
          <th>4</th>
          <td>iAF692</td>
          <td>692</td>
          <td>628</td>
          <td>Methanosarcina barkeri str. Fusaro</td>
          <td>690</td>
        </tr>
        <tr>
          <th>5</th>
          <td>iAF987</td>
          <td>987</td>
          <td>1109</td>
          <td>Geobacter metallireducens GS-15</td>
          <td>1285</td>
        </tr>
        <tr>
          <th>6</th>
          <td>iAPECO1_1312</td>
          <td>1313</td>
          <td>1942</td>
          <td>Escherichia coli APEC O1</td>
          <td>2735</td>
        </tr>
        <tr>
          <th>7</th>
          <td>iAT_PLT_636</td>
          <td>636</td>
          <td>738</td>
          <td>Homo sapiens</td>
          <td>1008</td>
        </tr>
        <tr>
          <th>8</th>
          <td>iB21_1397</td>
          <td>1337</td>
          <td>1943</td>
          <td>Escherichia coli BL21(DE3)</td>
          <td>2741</td>
        </tr>
        <tr>
          <th>9</th>
          <td>iBWG_1329</td>
          <td>1329</td>
          <td>1949</td>
          <td>Escherichia coli BW2952</td>
          <td>2741</td>
        </tr>
        <tr>
          <th>10</th>
          <td>ic_1306</td>
          <td>1307</td>
          <td>1936</td>
          <td>Escherichia coli CFT073</td>
          <td>2726</td>
        </tr>
        <tr>
          <th>11</th>
          <td>iCHOv1</td>
          <td>1766</td>
          <td>4456</td>
          <td>Cricetulus griseus</td>
          <td>6663</td>
        </tr>
        <tr>
          <th>12</th>
          <td>iE2348C_1286</td>
          <td>1287</td>
          <td>1919</td>
          <td>Escherichia coli O127:H6 str. E2348/69</td>
          <td>2703</td>
        </tr>
        <tr>
          <th>13</th>
          <td>iEC042_1314</td>
          <td>1314</td>
          <td>1926</td>
          <td>Escherichia coli 042</td>
          <td>2714</td>
        </tr>
        <tr>
          <th>14</th>
          <td>iEC55989_1330</td>
          <td>1330</td>
          <td>1953</td>
          <td>Escherichia coli 55989</td>
          <td>2756</td>
        </tr>
        <tr>
          <th>15</th>
          <td>iECABU_c1320</td>
          <td>1320</td>
          <td>1942</td>
          <td>Escherichia coli ABU 83972</td>
          <td>2731</td>
        </tr>
        <tr>
          <th>16</th>
          <td>iECB_1328</td>
          <td>1329</td>
          <td>1951</td>
          <td>Escherichia coli B str. REL606</td>
          <td>2748</td>
        </tr>
        <tr>
          <th>17</th>
          <td>iECBD_1354</td>
          <td>1354</td>
          <td>1952</td>
          <td>Escherichia coli 'BL21-Gold(DE3)pLysS AG'</td>
          <td>2748</td>
        </tr>
        <tr>
          <th>18</th>
          <td>iECD_1391</td>
          <td>1333</td>
          <td>1943</td>
          <td>Escherichia coli BL21(DE3)</td>
          <td>2741</td>
        </tr>
        <tr>
          <th>19</th>
          <td>iECDH10B_1368</td>
          <td>1327</td>
          <td>1947</td>
          <td>Escherichia coli str. K-12 substr. DH10B</td>
          <td>2742</td>
        </tr>
        <tr>
          <th>20</th>
          <td>iEcDH1_1363</td>
          <td>1363</td>
          <td>1949</td>
          <td>Escherichia coli DH1</td>
          <td>2750</td>
        </tr>
        <tr>
          <th>21</th>
          <td>iECDH1ME8569_1439</td>
          <td>1439</td>
          <td>1950</td>
          <td>Escherichia coli DH1</td>
          <td>2755</td>
        </tr>
        <tr>
          <th>22</th>
          <td>iEcE24377_1341</td>
          <td>1341</td>
          <td>1972</td>
          <td>Escherichia coli O139:H28 str. E24377A</td>
          <td>2763</td>
        </tr>
        <tr>
          <th>23</th>
          <td>iECED1_1282</td>
          <td>1279</td>
          <td>1929</td>
          <td>Escherichia coli ED1a</td>
          <td>2706</td>
        </tr>
        <tr>
          <th>24</th>
          <td>iECH74115_1262</td>
          <td>1262</td>
          <td>1918</td>
          <td>Escherichia coli O157:H7 str. EC4115</td>
          <td>2694</td>
        </tr>
        <tr>
          <th>25</th>
          <td>iEcHS_1320</td>
          <td>1321</td>
          <td>1963</td>
          <td>Escherichia coli HS</td>
          <td>2753</td>
        </tr>
        <tr>
          <th>26</th>
          <td>iECIAI1_1343</td>
          <td>1343</td>
          <td>1968</td>
          <td>Escherichia coli IAI1</td>
          <td>2765</td>
        </tr>
        <tr>
          <th>27</th>
          <td>iECIAI39_1322</td>
          <td>1321</td>
          <td>1953</td>
          <td>Escherichia coli IAI39</td>
          <td>2721</td>
        </tr>
        <tr>
          <th>28</th>
          <td>iECNA114_1301</td>
          <td>1301</td>
          <td>1927</td>
          <td>Escherichia coli NA114</td>
          <td>2718</td>
        </tr>
        <tr>
          <th>29</th>
          <td>iECO103_1326</td>
          <td>1327</td>
          <td>1958</td>
          <td>Escherichia coli O103:H2 str. 12009</td>
          <td>2758</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>54</th>
          <td>iLF82_1304</td>
          <td>1302</td>
          <td>1938</td>
          <td>Escherichia coli LF82</td>
          <td>2726</td>
        </tr>
        <tr>
          <th>55</th>
          <td>iLJ478</td>
          <td>482</td>
          <td>570</td>
          <td>Thermotoga maritima MSB8</td>
          <td>652</td>
        </tr>
        <tr>
          <th>56</th>
          <td>iML1515</td>
          <td>1516</td>
          <td>1877</td>
          <td>Escherichia coli str. K-12 substr. MG1655</td>
          <td>2712</td>
        </tr>
        <tr>
          <th>57</th>
          <td>iMM1415</td>
          <td>1375</td>
          <td>2775</td>
          <td>Mus musculus</td>
          <td>3726</td>
        </tr>
        <tr>
          <th>58</th>
          <td>iMM904</td>
          <td>905</td>
          <td>1226</td>
          <td>Saccharomyces cerevisiae S288C</td>
          <td>1577</td>
        </tr>
        <tr>
          <th>59</th>
          <td>iND750</td>
          <td>750</td>
          <td>1059</td>
          <td>Saccharomyces cerevisiae S288C</td>
          <td>1266</td>
        </tr>
        <tr>
          <th>60</th>
          <td>iNF517</td>
          <td>516</td>
          <td>650</td>
          <td>Lactococcus lactis subsp. cremoris MG1363</td>
          <td>754</td>
        </tr>
        <tr>
          <th>61</th>
          <td>iNJ661</td>
          <td>661</td>
          <td>825</td>
          <td>Mycobacterium tuberculosis H37Rv</td>
          <td>1025</td>
        </tr>
        <tr>
          <th>62</th>
          <td>iNRG857_1313</td>
          <td>1311</td>
          <td>1943</td>
          <td>Escherichia coli O83:H1 str. NRG 857C</td>
          <td>2735</td>
        </tr>
        <tr>
          <th>63</th>
          <td>iPC815</td>
          <td>815</td>
          <td>1552</td>
          <td>Yersinia pestis CO92</td>
          <td>1961</td>
        </tr>
        <tr>
          <th>64</th>
          <td>iRC1080</td>
          <td>1086</td>
          <td>1706</td>
          <td>Chlamydomonas reinhardtii</td>
          <td>2191</td>
        </tr>
        <tr>
          <th>65</th>
          <td>iS_1188</td>
          <td>1188</td>
          <td>1914</td>
          <td>Shigella flexneri 2a str. 2457T</td>
          <td>2619</td>
        </tr>
        <tr>
          <th>66</th>
          <td>iSB619</td>
          <td>619</td>
          <td>655</td>
          <td>Staphylococcus aureus subsp. aureus N315</td>
          <td>743</td>
        </tr>
        <tr>
          <th>67</th>
          <td>iSbBS512_1146</td>
          <td>1147</td>
          <td>1910</td>
          <td>Shigella boydii CDC 3083-94</td>
          <td>2591</td>
        </tr>
        <tr>
          <th>68</th>
          <td>iSBO_1134</td>
          <td>1134</td>
          <td>1908</td>
          <td>Shigella boydii Sb227</td>
          <td>2591</td>
        </tr>
        <tr>
          <th>69</th>
          <td>iSDY_1059</td>
          <td>1059</td>
          <td>1888</td>
          <td>Shigella dysenteriae Sd197</td>
          <td>2539</td>
        </tr>
        <tr>
          <th>70</th>
          <td>iSF_1195</td>
          <td>1195</td>
          <td>1917</td>
          <td>Shigella flexneri 2a str. 301</td>
          <td>2630</td>
        </tr>
        <tr>
          <th>71</th>
          <td>iSFV_1184</td>
          <td>1184</td>
          <td>1917</td>
          <td>Shigella flexneri 5 str. 8401</td>
          <td>2621</td>
        </tr>
        <tr>
          <th>72</th>
          <td>iSFxv_1172</td>
          <td>1169</td>
          <td>1918</td>
          <td>Shigella flexneri 2002017</td>
          <td>2638</td>
        </tr>
        <tr>
          <th>73</th>
          <td>iSSON_1240</td>
          <td>1240</td>
          <td>1936</td>
          <td>Shigella sonnei Ss046</td>
          <td>2693</td>
        </tr>
        <tr>
          <th>74</th>
          <td>iUMN146_1321</td>
          <td>1319</td>
          <td>1942</td>
          <td>Escherichia coli UM146</td>
          <td>2735</td>
        </tr>
        <tr>
          <th>75</th>
          <td>iUMNK88_1353</td>
          <td>1353</td>
          <td>1969</td>
          <td>Escherichia coli UMNK88</td>
          <td>2777</td>
        </tr>
        <tr>
          <th>76</th>
          <td>iUTI89_1310</td>
          <td>1310</td>
          <td>1940</td>
          <td>Escherichia coli UTI89</td>
          <td>2725</td>
        </tr>
        <tr>
          <th>77</th>
          <td>iWFL_1372</td>
          <td>1372</td>
          <td>1973</td>
          <td>Escherichia coli W</td>
          <td>2782</td>
        </tr>
        <tr>
          <th>78</th>
          <td>iY75_1357</td>
          <td>1358</td>
          <td>1953</td>
          <td>Escherichia coli str. K-12 substr. W3110</td>
          <td>2759</td>
        </tr>
        <tr>
          <th>79</th>
          <td>iYL1228</td>
          <td>1229</td>
          <td>1658</td>
          <td>Klebsiella pneumoniae subsp. pneumoniae MGH 78578</td>
          <td>2262</td>
        </tr>
        <tr>
          <th>80</th>
          <td>iYO844</td>
          <td>844</td>
          <td>990</td>
          <td>Bacillus subtilis subsp. subtilis str. 168</td>
          <td>1250</td>
        </tr>
        <tr>
          <th>81</th>
          <td>iZ_1308</td>
          <td>1308</td>
          <td>1923</td>
          <td>Escherichia coli O157:H7 str. EDL933</td>
          <td>2721</td>
        </tr>
        <tr>
          <th>82</th>
          <td>RECON1</td>
          <td>1905</td>
          <td>2766</td>
          <td>Homo sapiens</td>
          <td>3741</td>
        </tr>
        <tr>
          <th>83</th>
          <td>STM_v1_0</td>
          <td>1271</td>
          <td>1802</td>
          <td>Salmonella enterica subsp. enterica serovar Ty...</td>
          <td>2545</td>
        </tr>
      </tbody>
    </table>
    <p>84 rows × 5 columns</p>
    </div>



.. code:: ipython3

    models.index_models_minho()




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
          <th>id</th>
          <th>name</th>
          <th>doi</th>
          <th>author</th>
          <th>year</th>
          <th>formats</th>
          <th>organism</th>
          <th>taxonomy</th>
          <th>validated</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>1</td>
          <td>iJR904</td>
          <td>10.1186/gb-2003-4-9-r54</td>
          <td>Reed</td>
          <td>2003</td>
          <td>[sbml]</td>
          <td>Escherichia coli str. K12 substr. MG1655</td>
          <td>Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>1</th>
          <td>2</td>
          <td>iAF1260</td>
          <td>10.1038/msb4100155</td>
          <td>Feist</td>
          <td>2007</td>
          <td>[sbml]</td>
          <td>Escherichia coli str. K12 substr. MG1655</td>
          <td>Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>2</th>
          <td>3</td>
          <td>iMM904</td>
          <td>10.1186/1752-0509-3-37</td>
          <td>Mo</td>
          <td>2007</td>
          <td>[sbml]</td>
          <td>Saccharomyces cerevisiae</td>
          <td>Eukaryota; Opisthokonta; Fungi; Dikarya; Ascom...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>3</th>
          <td>4</td>
          <td>iJP815</td>
          <td>10.1371/journal.pcbi.1000210</td>
          <td>Puchalka</td>
          <td>2008</td>
          <td>[sbml]</td>
          <td>Pseudomonas putida str. KT2440</td>
          <td>Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>4</th>
          <td>5</td>
          <td>iMO1056</td>
          <td>10.1128/JB.01583-07</td>
          <td>Oberhardt</td>
          <td>2008</td>
          <td>[excel]</td>
          <td>Pseudomonas aeruginosa str. PAO1</td>
          <td>Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
          <td>False</td>
        </tr>
        <tr>
          <th>5</th>
          <td>6</td>
          <td>iIN800</td>
          <td>10.1186/1752-0509-2-71</td>
          <td>Nookaew</td>
          <td>2008</td>
          <td>[sbml]</td>
          <td>Saccharomyces cerevisiae</td>
          <td>Eukaryota; Opisthokonta; Fungi; Dikarya; Ascom...</td>
          <td>False</td>
        </tr>
        <tr>
          <th>6</th>
          <td>7</td>
          <td>iFF708</td>
          <td>10.1101/gr.234503</td>
          <td>Förster</td>
          <td>2003</td>
          <td>[sbml]</td>
          <td>Saccharomyces cerevisiae str. Sc288</td>
          <td>Eukaryota; Opisthokonta; Fungi; Dikarya; Ascom...</td>
          <td>False</td>
        </tr>
        <tr>
          <th>7</th>
          <td>8</td>
          <td>iCA1273</td>
          <td>10.1186/1471-2164-12-9</td>
          <td>Archer</td>
          <td>2011</td>
          <td>[sbml]</td>
          <td>Escherichia coli W</td>
          <td>Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>8</th>
          <td>9</td>
          <td>iJO1366</td>
          <td>10.1038/msb.2011.65</td>
          <td>Orth</td>
          <td>2011</td>
          <td>[sbml]</td>
          <td>Escherichia coli str. K12 substr. MG1655</td>
          <td>Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>9</th>
          <td>10</td>
          <td>yeast 4.0</td>
          <td>10.1186/1752-0509-4-145</td>
          <td>Dobson</td>
          <td>2010</td>
          <td>[]</td>
          <td>Yeast</td>
          <td>Eukaryota; Opisthokonta; Fungi;</td>
          <td>False</td>
        </tr>
        <tr>
          <th>10</th>
          <td>11</td>
          <td>iJN746</td>
          <td>10.1186/1752-0509-2-79</td>
          <td>Nogales</td>
          <td>2008</td>
          <td>[sbml]</td>
          <td>Pseudomonas putida str. KT2440</td>
          <td>Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>11</th>
          <td>12</td>
          <td>AbyMBEL891</td>
          <td>10.1039/B916446D</td>
          <td>Kim</td>
          <td>2010</td>
          <td>[sbml]</td>
          <td>Acinetobacter baumannii AYE</td>
          <td>Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
          <td>False</td>
        </tr>
        <tr>
          <th>12</th>
          <td>13</td>
          <td>iJP962</td>
          <td>10.1371/journal.pcbi.1001116</td>
          <td>Oberhardt</td>
          <td>2011</td>
          <td>[sbml]</td>
          <td>Pseudomonas putida str. KT2440</td>
          <td>Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>13</th>
          <td>14</td>
          <td>iYL1228</td>
          <td>10.1128/JB.01218-10</td>
          <td>Liao</td>
          <td>2011</td>
          <td>[sbml]</td>
          <td>Klebsiella pneumoniae str. MGH 78578</td>
          <td>Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>14</th>
          <td>15</td>
          <td>iSR432</td>
          <td>10.1186/1752-0509-4-31</td>
          <td>Roberts</td>
          <td>2010</td>
          <td>[sbml]</td>
          <td>Clostridium thermocellum str. ATCC 27405</td>
          <td>Bacteria; Firmicutes; Clostridia; Clostridial...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>15</th>
          <td>16</td>
          <td>iNJ661m</td>
          <td>10.1186/1752-0509-4-160</td>
          <td>Fang</td>
          <td>2010</td>
          <td>[sbml]</td>
          <td>Mycobacterium tuberculosis str. H37Rv</td>
          <td>Bacteria; Actinobacteria; Actinobacteria; Acti...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>16</th>
          <td>17</td>
          <td>iCM925</td>
          <td>10.1186/1752-0509-5-130</td>
          <td>Milne</td>
          <td>2011</td>
          <td>[sbml]</td>
          <td>Clostridium beijerinckii str. NCIMB 8052</td>
          <td>Bacteria; Firmicutes; Clostridia; Clostridiale...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>17</th>
          <td>18</td>
          <td>iBsu1103</td>
          <td>10.1186/gb-2009-10-6-r69</td>
          <td>Henry</td>
          <td>2009</td>
          <td>[sbml]</td>
          <td>Bacillus subtilis 168</td>
          <td>Bacteria; Firmicutes; Bacilli; Bacillales; Bac...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>18</th>
          <td>19</td>
          <td>iAI549</td>
          <td>10.1371/journal.pcbi.1000887</td>
          <td>Islam</td>
          <td>2010</td>
          <td>[sbml]</td>
          <td>Dehalococcoides ethenogenes str. 2061</td>
          <td>Bacteria; Chloroflexi; Dehalococcoidetes; Deha...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>19</th>
          <td>20</td>
          <td>iAF692</td>
          <td>10.1038/msb4100046</td>
          <td>Feist</td>
          <td>2006</td>
          <td>[sbml]</td>
          <td>Methanosarcina barkeri</td>
          <td>Archaea; Euryarchaeota; Methanomicrobia; Metha...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>20</th>
          <td>21</td>
          <td>AraGEM</td>
          <td>10.1186/1471-2164-12-S4-S5</td>
          <td>de Oliveira Dal'Molin</td>
          <td>2010</td>
          <td>[sbml]</td>
          <td>Arabidopsis thaliana</td>
          <td>Eukaryota; Viridiplantae; Streptophyta; Strept...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>21</th>
          <td>22</td>
          <td>Ecoli core Model</td>
          <td>doi:10.1128/ecosalplus.10.2.1</td>
          <td>Orth</td>
          <td>2010</td>
          <td>[sbml]</td>
          <td>Escherichia coli str. K12 substr. MG1655</td>
          <td>Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>22</th>
          <td>23</td>
          <td>iIT341</td>
          <td>10.1128/JB.187.16.5818-5830.2005</td>
          <td>Thiele</td>
          <td>2005</td>
          <td>[sbml]</td>
          <td>Helicobacter pylori str. 26695</td>
          <td>Bacteria; Proteobacteria; delta/epsilon subdi...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>23</th>
          <td>24</td>
          <td>iMH805/775</td>
          <td>10.1038/nbt1492</td>
          <td>Herrgård</td>
          <td>2008</td>
          <td>[sbml]</td>
          <td>Saccharomyces cerevisiae str. Sc288</td>
          <td>Eukaryota; Opisthokonta; Fungi; Dikarya; Ascom...</td>
          <td>False</td>
        </tr>
        <tr>
          <th>24</th>
          <td>25</td>
          <td>iND750</td>
          <td>10.1101/gr.2250904</td>
          <td>Duarte</td>
          <td>2004</td>
          <td>[sbml]</td>
          <td>Saccharomyces cerevisiae str. Sc288</td>
          <td>Eukaryota; Opisthokonta; Fungi; Dikarya; Ascom...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>25</th>
          <td>26</td>
          <td>iRC1080</td>
          <td>10.1038/msb.2011.52</td>
          <td>Chang</td>
          <td>2011</td>
          <td>[sbml]</td>
          <td>Chlamydomonas reinhardtii</td>
          <td>Eukaryota; Viridiplantae; Chlorophyta; Chlorop...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>26</th>
          <td>27</td>
          <td>iSB619</td>
          <td>10.1186/1471-2180-5-8</td>
          <td>Becker</td>
          <td>2005</td>
          <td>[sbml]</td>
          <td>Staphylococcus aureus</td>
          <td>Bacteria; Firmicutes; Bacilli; Bacillales; St...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>27</th>
          <td>28</td>
          <td>iTH366</td>
          <td>10.1038/msb.2010.60</td>
          <td>Plata</td>
          <td>2010</td>
          <td>[sbml]</td>
          <td>Plasmodium falciparum</td>
          <td>Eukaryota; Alveolata; Apicomplexa; Aconoidasid...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>28</th>
          <td>29</td>
          <td>iTZ479</td>
          <td>10.1126/science.1174671</td>
          <td>Zhang</td>
          <td>2009</td>
          <td>[sbml]</td>
          <td>Thermotoga maritima str. MSB8</td>
          <td>Bacteria; Thermotogae; Thermotogae; Thermotoga...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>29</th>
          <td>30</td>
          <td>recon2</td>
          <td>10.1038/nbt.2488</td>
          <td>Thiele</td>
          <td>2013</td>
          <td>[sbml]</td>
          <td>Homo sapiens</td>
          <td>Eukaryota; Opisthokonta; Metazoa; Eumetazoa; B...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>118</th>
          <td>123</td>
          <td>iCyc792</td>
          <td>10.1186/1752-0509-7-142</td>
          <td>Mueller</td>
          <td>2013</td>
          <td>[sbml, excel]</td>
          <td>Cyanothece sp. PCC 7424</td>
          <td>Bacteria; Cyanobacteria; Oscillatoriophycideae...</td>
          <td>False</td>
        </tr>
        <tr>
          <th>119</th>
          <td>124</td>
          <td>iCyn731</td>
          <td>10.1186/1752-0509-7-142</td>
          <td>Mueller</td>
          <td>2013</td>
          <td>[sbml, excel]</td>
          <td>Cyanothece sp. PCC 7425</td>
          <td>Bacteria; Cyanobacteria; Oscillatoriophycideae...</td>
          <td>False</td>
        </tr>
        <tr>
          <th>120</th>
          <td>125</td>
          <td>iCyj826</td>
          <td>10.1186/1752-0509-7-142</td>
          <td>Mueller</td>
          <td>2013</td>
          <td>[sbml, excel]</td>
          <td>Cyanothece sp. PCC 7822</td>
          <td>Bacteria; Cyanobacteria; Oscillatoriophycideae...</td>
          <td>False</td>
        </tr>
        <tr>
          <th>121</th>
          <td>126</td>
          <td>iCyp752</td>
          <td>10.1186/1752-0509-7-142</td>
          <td>Mueller</td>
          <td>2013</td>
          <td>[sbml, excel]</td>
          <td>Cyanothece sp. PCC 8801</td>
          <td>Bacteria; Cyanobacteria; Oscillatoriophycideae...</td>
          <td>False</td>
        </tr>
        <tr>
          <th>122</th>
          <td>127</td>
          <td>iCyh755</td>
          <td>10.1186/1752-0509-7-142</td>
          <td>Mueller</td>
          <td>2013</td>
          <td>[sbml, excel]</td>
          <td>Cyanothece sp. PCC 8802</td>
          <td>Bacteria; Cyanobacteria; Oscillatoriophycideae...</td>
          <td>False</td>
        </tr>
        <tr>
          <th>123</th>
          <td>128</td>
          <td>iNF518</td>
          <td>10.1007/s00253-013-5140-2</td>
          <td>Flahaut</td>
          <td>2013</td>
          <td>[sbml]</td>
          <td>Lactococcus lactis subsp. cremoris MG1363</td>
          <td>Bacteria; Firmicutes; Bacilli; Lactobacillales...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>124</th>
          <td>129</td>
          <td>iJL1454</td>
          <td>10.1039/C3MB70090A</td>
          <td>Jie Liu</td>
          <td>2013</td>
          <td>[excel]</td>
          <td>Aspergillus terreus NIH2624</td>
          <td>Eukaryota; Opisthokonta; Fungi; Dikarya; Ascom...</td>
          <td>False</td>
        </tr>
        <tr>
          <th>125</th>
          <td>130</td>
          <td>iBif452</td>
          <td>10.1186/1752-0509-8-41</td>
          <td>El-Semman</td>
          <td>2014</td>
          <td>[sbml, excel]</td>
          <td>Bifidobacterium adolescentis L2-32</td>
          <td>Bacteria; Actinobacteria; Actinobacteria; Acti...</td>
          <td>False</td>
        </tr>
        <tr>
          <th>126</th>
          <td>131</td>
          <td>iFap484</td>
          <td>10.1186/1752-0509-8-41</td>
          <td>El-Semman</td>
          <td>2014</td>
          <td>[sbml, excel]</td>
          <td>Faecalibacterium prausnitzii A2-165</td>
          <td>Bacteria; Firmicutes; Clostridia; Clostridiale...</td>
          <td>False</td>
        </tr>
        <tr>
          <th>127</th>
          <td>132</td>
          <td>iAM388</td>
          <td>10.1186/1471-2164-12-535</td>
          <td>Aline Metris</td>
          <td>2011</td>
          <td>[excel]</td>
          <td>Campylobacter jejuni subsp. jejuni NCTC 11168</td>
          <td>Bacteria; Proteobacteria; delta/epsilon subdiv...</td>
          <td>False</td>
        </tr>
        <tr>
          <th>128</th>
          <td>133</td>
          <td>iTT548</td>
          <td>10.1186/1475-2859-13-61</td>
          <td>Na-Rae Lee</td>
          <td>2014</td>
          <td>[sbml, excel]</td>
          <td>Thermus thermophilus</td>
          <td>Bacteria; Deinococcus-Thermus; Deinococci; The...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>129</th>
          <td>134</td>
          <td>EctoGEM-1.0</td>
          <td>10.1111/tpj.12627</td>
          <td>Prigent</td>
          <td>2014</td>
          <td>[sbml]</td>
          <td>Ectocarpus siliculosus</td>
          <td>Eukaryota; Stramenopiles; PX clade; Phaeophyce...</td>
          <td>False</td>
        </tr>
        <tr>
          <th>130</th>
          <td>135</td>
          <td>iMF721</td>
          <td>10.1111/1462-2920.12513</td>
          <td>Fondi</td>
          <td>2014</td>
          <td>[sbml]</td>
          <td>Pseudoalteromonas haloplanktis TAC125</td>
          <td>Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
          <td>False</td>
        </tr>
        <tr>
          <th>131</th>
          <td>136</td>
          <td>Arabidopsis core model</td>
          <td>10.1104/pp.114.235358</td>
          <td>Arnold</td>
          <td>2014</td>
          <td>[sbml]</td>
          <td>Arabidopsis thaliana</td>
          <td>Eukaryota; Viridiplantae; Streptophyta; Strept...</td>
          <td>False</td>
        </tr>
        <tr>
          <th>132</th>
          <td>137</td>
          <td>iHN637</td>
          <td>10.1186/1475-2859-12-118</td>
          <td>Harish Nagarajan</td>
          <td>2013</td>
          <td>[excel]</td>
          <td>Clostridium ljungdahlii</td>
          <td>Bacteria; Firmicutes; Clostridia; Clostridiale...</td>
          <td>False</td>
        </tr>
        <tr>
          <th>133</th>
          <td>138</td>
          <td>iCac802</td>
          <td>10.1186/s13068-014-0144-4</td>
          <td>Satyakam Dash</td>
          <td>2014</td>
          <td>[sbml]</td>
          <td>Clostridium acetobutylicum ATCC 824</td>
          <td>Bacteria; Firmicutes; Clostridia; Clostridiale...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>134</th>
          <td>139</td>
          <td>iMLTC806cdf</td>
          <td>10.1186/s12918-014-0117-z</td>
          <td>M. Larocque</td>
          <td>2014</td>
          <td>[sbml, excel]</td>
          <td>Clostridium difficile 630</td>
          <td>Bacteria; Firmicutes; Clostridia; Clostridiale...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>135</th>
          <td>140</td>
          <td>iCY1106</td>
          <td>10.1186/s12918-014-0137-8</td>
          <td>Chao Ye</td>
          <td>2015</td>
          <td>[sbml]</td>
          <td>Mortierella alpina</td>
          <td>Eukaryota; Opisthokonta; Fungi; Fungi incertae...</td>
          <td>False</td>
        </tr>
        <tr>
          <th>136</th>
          <td>141</td>
          <td>iMM518</td>
          <td>10.1039/c3mb70421a</td>
          <td>N. Goyal</td>
          <td>2014</td>
          <td>[excel]</td>
          <td>Methanococcus maripaludis S2</td>
          <td>Archaea; Euryarchaeota; Methanococci; Methanoc...</td>
          <td>False</td>
        </tr>
        <tr>
          <th>137</th>
          <td>142</td>
          <td>iPC1209</td>
          <td>10.1016/j.febslet.2014.12.010</td>
          <td>Cheng Wang</td>
          <td>2015</td>
          <td>[]</td>
          <td>Pectobacterium carotovorum subsp. carotovorum PC1</td>
          <td>Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
          <td>False</td>
        </tr>
        <tr>
          <th>138</th>
          <td>143</td>
          <td>iNV706</td>
          <td>10.1128/AEM.03279-14</td>
          <td>N. Veith</td>
          <td>2014</td>
          <td>[]</td>
          <td>Enterococcus faecalis V583</td>
          <td>Bacteria; Firmicutes; Bacilli; Lactobacillales...</td>
          <td>False</td>
        </tr>
        <tr>
          <th>139</th>
          <td>144</td>
          <td>KoxGSC1457</td>
          <td>10.1186/1475-2859-12-20</td>
          <td>J. Park</td>
          <td>2013</td>
          <td>[]</td>
          <td>Klebsiella oxytoca</td>
          <td>Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
          <td>False</td>
        </tr>
        <tr>
          <th>140</th>
          <td>145</td>
          <td>iJL846</td>
          <td>10.1016/j.gene.2014.10.034</td>
          <td>Jie Liu</td>
          <td>2015</td>
          <td>[]</td>
          <td>Lactobacillus casei LC2W</td>
          <td>Bacteria; Firmicutes; Bacilli; Lactobacillales...</td>
          <td>False</td>
        </tr>
        <tr>
          <th>141</th>
          <td>146</td>
          <td>yeast 7.6</td>
          <td>10.1089/ind.2013.0013</td>
          <td>Aung (updated)</td>
          <td>2015</td>
          <td>[sbml]</td>
          <td>Yeast</td>
          <td>Eukaryota; Opisthokonta; Fungi;</td>
          <td>False</td>
        </tr>
        <tr>
          <th>142</th>
          <td>147</td>
          <td>iJK849</td>
          <td>10.1111/tpj.13081</td>
          <td>Joomi Kim</td>
          <td>2015</td>
          <td>[sbml, excel]</td>
          <td>Phaeodactylum tricornutum</td>
          <td>Eukaryota; Stramenopiles; Bacillariophyta; Bac...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>143</th>
          <td>148</td>
          <td>iNL895</td>
          <td>10.1186/1752-0509-6-35</td>
          <td>Nicolas Loira</td>
          <td>2012</td>
          <td>[sbml]</td>
          <td>Yarrowia lipolytica</td>
          <td>Eukaryota; Fungi; Dikarya; Ascomycota; Sacchar...</td>
          <td>False</td>
        </tr>
        <tr>
          <th>144</th>
          <td>149</td>
          <td>iYali4</td>
          <td>10.1038/npjsba.2016.5</td>
          <td>Kerkhoven</td>
          <td>2016</td>
          <td>[sbml]</td>
          <td>Yarrowia lipolytica</td>
          <td>Eukaryota; Fungi; Dikarya; Ascomycota; Sacchar...</td>
          <td>False</td>
        </tr>
        <tr>
          <th>145</th>
          <td>150</td>
          <td>iLB1027_lipid</td>
          <td>10.1371/journal.pone.0155038</td>
          <td>Jennifer Levering</td>
          <td>2016</td>
          <td>[sbml]</td>
          <td>Phaeodactylum tricornutum</td>
          <td>Eukaryota; Stramenopiles; Bacillariophyta; Bac...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>146</th>
          <td>151</td>
          <td>iLB1027</td>
          <td>10.1371/journal.pone.0155038</td>
          <td>Jennifer Levering</td>
          <td>2016</td>
          <td>[sbml]</td>
          <td>Phaeodactylum tricornutum</td>
          <td>Eukaryota; Stramenopiles; Bacillariophyta; Bac...</td>
          <td>True</td>
        </tr>
        <tr>
          <th>147</th>
          <td>152</td>
          <td>iMT1174</td>
          <td>10.1186/s12918-015-0190-y</td>
          <td>Mohammad Tajparast</td>
          <td>2015</td>
          <td>[excel]</td>
          <td>Rhodococcus jostii RHA1</td>
          <td>Bacteria; Terrabacteria group; Actinobacteria;...</td>
          <td>False</td>
        </tr>
      </tbody>
    </table>
    <p>148 rows × 9 columns</p>
    </div>



Models from `BiGG <http://bigg.ucsd.edu>`__ and the `University of
Minho <http://darwin.di.uminho.pt/models>`__ can conveniently be accessd
via `~cameo.models.bigg` and `~cameo.models.minho` respectively.

.. code:: ipython3

    models.bigg.iJN746




.. raw:: html

    
            <table>
                <tr>
                    <td><strong>Name</strong></td>
                    <td>iJN746</td>
                </tr><tr>
                    <td><strong>Memory address</strong></td>
                    <td>0x01106d16d8</td>
                </tr><tr>
                    <td><strong>Number of metabolites</strong></td>
                    <td>907</td>
                </tr><tr>
                    <td><strong>Number of reactions</strong></td>
                    <td>1054</td>
                </tr><tr>
                    <td><strong>Objective expression</strong></td>
                    <td>-1.0*BIOMASS_KT_TEMP_reverse_d18f7 + 1.0*BIOMASS_KT_TEMP</td>
                </tr><tr>
                    <td><strong>Compartments</strong></td>
                    <td>extracellular space, cytosol, periplasm</td>
                </tr>
              </table>



.. code:: ipython3

    models.minho.iMM904




.. raw:: html

    
            <table>
                <tr>
                    <td><strong>Name</strong></td>
                    <td>iMM904</td>
                </tr><tr>
                    <td><strong>Memory address</strong></td>
                    <td>0x0115e79a58</td>
                </tr><tr>
                    <td><strong>Number of metabolites</strong></td>
                    <td>1228</td>
                </tr><tr>
                    <td><strong>Number of reactions</strong></td>
                    <td>1577</td>
                </tr><tr>
                    <td><strong>Objective expression</strong></td>
                    <td>-1.0*biomass_SC5_notrace_reverse_e32ff + 1.0*biomass_SC5_notrace</td>
                </tr><tr>
                    <td><strong>Compartments</strong></td>
                    <td>Golgi_Apparatus, Extra_organism, Nucleus, Endoplasmic_Reticulum, Cytosol, Peroxisome, Mitochondria, Vacuole</td>
                </tr>
              </table>



Models in the Minho database have been manually verified. The subset of
models shown bellow can be used to run simulations as described in the
publications.

.. code:: ipython3

    models.minho.validated.VvuMBEL943 # use TAB completion to see the other models




.. raw:: html

    
            <table>
                <tr>
                    <td><strong>Name</strong></td>
                    <td>HyunUkKim2010_VvuMBEL943_MetabolicModeling</td>
                </tr><tr>
                    <td><strong>Memory address</strong></td>
                    <td>0x010c1676a0</td>
                </tr><tr>
                    <td><strong>Number of metabolites</strong></td>
                    <td>912</td>
                </tr><tr>
                    <td><strong>Number of reactions</strong></td>
                    <td>1019</td>
                </tr><tr>
                    <td><strong>Objective expression</strong></td>
                    <td>0</td>
                </tr><tr>
                    <td><strong>Compartments</strong></td>
                    <td>cell</td>
                </tr>
              </table>


