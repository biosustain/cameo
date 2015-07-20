
.. code:: python

    from cameo import webmodels

.. code:: python

    models = webmodels.index_models()

.. code:: python

    models




.. raw:: html

    <div style="max-height:1000px;max-width:1500px;overflow:auto;">
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
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0  </th>
          <td>   1</td>
          <td>                 iJR904</td>
          <td>          10.1186/gb-2003-4-9-r54</td>
          <td>                  Reed</td>
          <td> 2003</td>
          <td>        [sbml]</td>
          <td>          Escherichia coli str. K12 substr. MG1655</td>
          <td> Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
        </tr>
        <tr>
          <th>1  </th>
          <td>   2</td>
          <td>                iAF1260</td>
          <td>               10.1038/msb4100155</td>
          <td>                 Feist</td>
          <td> 2007</td>
          <td>        [sbml]</td>
          <td>          Escherichia coli str. K12 substr. MG1655</td>
          <td> Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
        </tr>
        <tr>
          <th>2  </th>
          <td>   3</td>
          <td>                 iMM904</td>
          <td>           10.1186/1752-0509-3-37</td>
          <td>                    Mo</td>
          <td> 2007</td>
          <td>        [sbml]</td>
          <td>                          Saccharomyces cerevisiae</td>
          <td> Eukaryota; Opisthokonta; Fungi; Dikarya; Ascom...</td>
        </tr>
        <tr>
          <th>3  </th>
          <td>   4</td>
          <td>                 iJP815</td>
          <td>     10.1371/journal.pcbi.1000210</td>
          <td>              Puchalka</td>
          <td> 2008</td>
          <td>        [sbml]</td>
          <td>                    Pseudomonas putida str. KT2440</td>
          <td> Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
        </tr>
        <tr>
          <th>4  </th>
          <td>   5</td>
          <td>                iMO1056</td>
          <td>              10.1128/JB.01583-07</td>
          <td>             Oberhardt</td>
          <td> 2008</td>
          <td>            []</td>
          <td>                  Pseudomonas aeruginosa str. PAO1</td>
          <td> Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
        </tr>
        <tr>
          <th>5  </th>
          <td>   6</td>
          <td>                 iIN800</td>
          <td>           10.1186/1752-0509-2-71</td>
          <td>               Nookaew</td>
          <td> 2008</td>
          <td>        [sbml]</td>
          <td>                          Saccharomyces cerevisiae</td>
          <td> Eukaryota; Opisthokonta; Fungi; Dikarya; Ascom...</td>
        </tr>
        <tr>
          <th>6  </th>
          <td>   7</td>
          <td>                 iFF708</td>
          <td>                10.1101/gr.234503</td>
          <td>               Förster</td>
          <td> 2003</td>
          <td>        [sbml]</td>
          <td>              Saccharomyces cerevisiae str. Sc288 </td>
          <td> Eukaryota; Opisthokonta; Fungi; Dikarya; Ascom...</td>
        </tr>
        <tr>
          <th>7  </th>
          <td>   8</td>
          <td>                iCA1273</td>
          <td>           10.1186/1471-2164-12-9</td>
          <td>                Archer</td>
          <td> 2011</td>
          <td>        [sbml]</td>
          <td>                                Escherichia coli W</td>
          <td> Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
        </tr>
        <tr>
          <th>8  </th>
          <td>   9</td>
          <td>                iJO1366</td>
          <td>              10.1038/msb.2011.65</td>
          <td>                  Orth</td>
          <td> 2011</td>
          <td>        [sbml]</td>
          <td>          Escherichia coli str. K12 substr. MG1655</td>
          <td> Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
        </tr>
        <tr>
          <th>9  </th>
          <td>  10</td>
          <td>           Consensus v4</td>
          <td>                  10.1038/nbt1492</td>
          <td>              Herrgård</td>
          <td> 2008</td>
          <td>            []</td>
          <td>                                             Yeast</td>
          <td>                   Eukaryota; Opisthokonta; Fungi;</td>
        </tr>
        <tr>
          <th>10 </th>
          <td>  11</td>
          <td>                 iJN746</td>
          <td>           10.1186/1752-0509-2-79</td>
          <td>               Nogales</td>
          <td> 2008</td>
          <td>        [sbml]</td>
          <td>                    Pseudomonas putida str. KT2440</td>
          <td> Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
        </tr>
        <tr>
          <th>11 </th>
          <td>  12</td>
          <td>             AbyMBEL891</td>
          <td>                 10.1039/B916446D</td>
          <td>                   Kim</td>
          <td> 2010</td>
          <td>        [sbml]</td>
          <td>                       Acinetobacter baumannii AYE</td>
          <td> Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
        </tr>
        <tr>
          <th>12 </th>
          <td>  13</td>
          <td>                 iJP962</td>
          <td>     10.1371/journal.pcbi.1001116</td>
          <td>             Oberhardt</td>
          <td> 2011</td>
          <td>        [sbml]</td>
          <td>                    Pseudomonas putida str. KT2440</td>
          <td> Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
        </tr>
        <tr>
          <th>13 </th>
          <td>  14</td>
          <td>                iYL1228</td>
          <td>              10.1128/JB.01218-10</td>
          <td>                  Liao</td>
          <td> 2011</td>
          <td>        [sbml]</td>
          <td>              Klebsiella pneumoniae str. MGH 78578</td>
          <td> Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
        </tr>
        <tr>
          <th>14 </th>
          <td>  15</td>
          <td>                 iSR432</td>
          <td>           10.1186/1752-0509-4-31</td>
          <td>               Roberts</td>
          <td> 2010</td>
          <td>        [sbml]</td>
          <td>          Clostridium thermocellum str. ATCC 27405</td>
          <td>  Bacteria; Firmicutes; Clostridia; Clostridial...</td>
        </tr>
        <tr>
          <th>15 </th>
          <td>  16</td>
          <td>                iNJ661m</td>
          <td>          10.1186/1752-0509-4-160</td>
          <td>                  Fang</td>
          <td> 2010</td>
          <td>        [sbml]</td>
          <td>             Mycobacterium tuberculosis str. H37Rv</td>
          <td> Bacteria; Actinobacteria; Actinobacteria; Acti...</td>
        </tr>
        <tr>
          <th>16 </th>
          <td>  17</td>
          <td>                 iCM925</td>
          <td>          10.1186/1752-0509-5-130</td>
          <td>                 Milne</td>
          <td> 2011</td>
          <td>        [sbml]</td>
          <td>          Clostridium beijerinckii str. NCIMB 8052</td>
          <td> Bacteria; Firmicutes; Clostridia; Clostridiale...</td>
        </tr>
        <tr>
          <th>17 </th>
          <td>  18</td>
          <td>               iBsu1103</td>
          <td>         10.1186/gb-2009-10-6-r69</td>
          <td>                 Henry</td>
          <td> 2009</td>
          <td>        [sbml]</td>
          <td>                             Bacillus subtilis 168</td>
          <td> Bacteria; Firmicutes; Bacilli; Bacillales; Bac...</td>
        </tr>
        <tr>
          <th>18 </th>
          <td>  19</td>
          <td>                 iAI549</td>
          <td>     10.1371/journal.pcbi.1000887</td>
          <td>                 Islam</td>
          <td> 2010</td>
          <td>        [sbml]</td>
          <td>             Dehalococcoides ethenogenes str. 2061</td>
          <td> Bacteria; Chloroflexi; Dehalococcoidetes; Deha...</td>
        </tr>
        <tr>
          <th>19 </th>
          <td>  20</td>
          <td>                 iAF692</td>
          <td>               10.1038/msb4100046</td>
          <td>                 Feist</td>
          <td> 2006</td>
          <td>        [sbml]</td>
          <td>                            Methanosarcina barkeri</td>
          <td> Archaea; Euryarchaeota; Methanomicrobia; Metha...</td>
        </tr>
        <tr>
          <th>20 </th>
          <td>  21</td>
          <td>                 AraGEM</td>
          <td>       10.1186/1471-2164-12-S4-S5</td>
          <td> de Oliveira Dal'Molin</td>
          <td> 2010</td>
          <td>        [sbml]</td>
          <td>                              Arabidopsis thaliana</td>
          <td> Eukaryota; Viridiplantae; Streptophyta; Strept...</td>
        </tr>
        <tr>
          <th>21 </th>
          <td>  22</td>
          <td>       Ecoli core Model</td>
          <td>            10.1128/ecosal.10.2.1</td>
          <td>                  Orth</td>
          <td> 2010</td>
          <td>        [sbml]</td>
          <td>          Escherichia coli str. K12 substr. MG1655</td>
          <td> Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
        </tr>
        <tr>
          <th>22 </th>
          <td>  23</td>
          <td>                 iIT341</td>
          <td> 10.1128/JB.187.16.5818-5830.2005</td>
          <td>                Thiele</td>
          <td> 2005</td>
          <td>        [sbml]</td>
          <td>                    Helicobacter pylori str. 26695</td>
          <td>  Bacteria; Proteobacteria; delta/epsilon subdi...</td>
        </tr>
        <tr>
          <th>23 </th>
          <td>  24</td>
          <td>             iMH805/775</td>
          <td>                  10.1038/nbt1492</td>
          <td>              Herrgård</td>
          <td> 2008</td>
          <td>        [sbml]</td>
          <td>              Saccharomyces cerevisiae str. Sc288 </td>
          <td> Eukaryota; Opisthokonta; Fungi; Dikarya; Ascom...</td>
        </tr>
        <tr>
          <th>24 </th>
          <td>  25</td>
          <td>                 iND750</td>
          <td>               10.1101/gr.2250904</td>
          <td>                Duarte</td>
          <td> 2004</td>
          <td>        [sbml]</td>
          <td>              Saccharomyces cerevisiae str. Sc288 </td>
          <td> Eukaryota; Opisthokonta; Fungi; Dikarya; Ascom...</td>
        </tr>
        <tr>
          <th>25 </th>
          <td>  26</td>
          <td>                iRC1080</td>
          <td>              10.1038/msb.2011.52</td>
          <td>                 Chang</td>
          <td> 2011</td>
          <td>        [sbml]</td>
          <td>                         Chlamydomonas reinhardtii</td>
          <td> Eukaryota; Viridiplantae; Chlorophyta; Chlorop...</td>
        </tr>
        <tr>
          <th>26 </th>
          <td>  27</td>
          <td>                 iSB619</td>
          <td>            10.1186/1471-2180-5-8</td>
          <td>                Becker</td>
          <td> 2005</td>
          <td>        [sbml]</td>
          <td>                             Staphylococcus aureus</td>
          <td>  Bacteria; Firmicutes; Bacilli; Bacillales; St...</td>
        </tr>
        <tr>
          <th>27 </th>
          <td>  28</td>
          <td>                 iTH366</td>
          <td>              10.1038/msb.2010.60</td>
          <td>                 Plata</td>
          <td> 2010</td>
          <td>        [sbml]</td>
          <td>                             Plasmodium falciparum</td>
          <td> Eukaryota; Alveolata; Apicomplexa; Aconoidasid...</td>
        </tr>
        <tr>
          <th>28 </th>
          <td>  29</td>
          <td>                 iTZ479</td>
          <td>          10.1126/science.1174671</td>
          <td>                 Zhang</td>
          <td> 2009</td>
          <td>        [sbml]</td>
          <td>                     Thermotoga maritima str. MSB8</td>
          <td> Bacteria; Thermotogae; Thermotogae; Thermotoga...</td>
        </tr>
        <tr>
          <th>29 </th>
          <td>  30</td>
          <td>                 recon2</td>
          <td>                 10.1038/nbt.2488</td>
          <td>                Thiele</td>
          <td> 2013</td>
          <td>        [sbml]</td>
          <td>                                      Homo sapiens</td>
          <td> Eukaryota; Opisthokonta; Metazoa; Eumetazoa; B...</td>
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
        </tr>
        <tr>
          <th>114</th>
          <td> 115</td>
          <td>                 iAK692</td>
          <td>           10.1186/1752-0509-6-71</td>
          <td>              Klanchui</td>
          <td> 2012</td>
          <td>        [sbml]</td>
          <td>                            Spirulina platensis C1</td>
          <td> Bacteria; Cyanobacteria; Oscillatoriophycideae...</td>
        </tr>
        <tr>
          <th>115</th>
          <td> 116</td>
          <td>                 iSS352</td>
          <td>    10.1016/j.jbiotec.2013.01.023</td>
          <td>        Schatschneider</td>
          <td> 2013</td>
          <td>        [sbml]</td>
          <td>   Xanthomonas campestris pv. campestris str. B100</td>
          <td> Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
        </tr>
        <tr>
          <th>116</th>
          <td> 117</td>
          <td>                 iSS352</td>
          <td>    10.1016/j.jbiotec.2013.01.023</td>
          <td>        Schatschneider</td>
          <td> 2013</td>
          <td>        [sbml]</td>
          <td>   Xanthomonas campestris pv. campestris str. B100</td>
          <td> Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
        </tr>
        <tr>
          <th>117</th>
          <td> 118</td>
          <td>                 iOD907</td>
          <td>                                 </td>
          <td>            Oscar Dias</td>
          <td> 2013</td>
          <td>        [sbml]</td>
          <td>                  Kluyveromyces lactis NRRL Y-1140</td>
          <td> Eukaryota; Opisthokonta; Fungi; Dikarya; Ascom...</td>
        </tr>
        <tr>
          <th>118</th>
          <td> 119</td>
          <td>                iCAC490</td>
          <td>           10.1186/1752-0509-6-42</td>
          <td>              McAnulty</td>
          <td> 2012</td>
          <td>        [sbml]</td>
          <td>               Clostridium acetobutylicum ATCC 824</td>
          <td> Bacteria; Firmicutes; Clostridia; Clostridiale...</td>
        </tr>
        <tr>
          <th>119</th>
          <td> 120</td>
          <td>              iTZ479_v2</td>
          <td>   10.1016/j.ijhydene.2012.06.032</td>
          <td>               Nogales</td>
          <td> 2012</td>
          <td>            []</td>
          <td>                               Thermotoga maritima</td>
          <td> Bacteria; Thermotogae; Thermotogae; Thermotoga...</td>
        </tr>
        <tr>
          <th>120</th>
          <td> 121</td>
          <td>                 iAH991</td>
          <td>               10.4161/gmic.22370</td>
          <td>               Heinken</td>
          <td> 2012</td>
          <td>            []</td>
          <td>                      Bacteroides thetaiotaomicron</td>
          <td> Bacteria; Bacteroidetes/Chlorobi group; Bacter...</td>
        </tr>
        <tr>
          <th>121</th>
          <td> 122</td>
          <td>                SPNV1.0</td>
          <td>          10.1186/1475-2859-13-41</td>
          <td>                  Wang</td>
          <td> 2014</td>
          <td>       [excel]</td>
          <td>              Saccharopolyspora spinosa NRRL 18395</td>
          <td> Bacteria; Actinobacteria; Actinobacteria; Acti...</td>
        </tr>
        <tr>
          <th>122</th>
          <td> 123</td>
          <td>                iCyc792</td>
          <td>          10.1186/1752-0509-7-142</td>
          <td>               Mueller</td>
          <td> 2013</td>
          <td> [sbml, excel]</td>
          <td>                           Cyanothece sp. PCC 7424</td>
          <td> Bacteria; Cyanobacteria; Oscillatoriophycideae...</td>
        </tr>
        <tr>
          <th>123</th>
          <td> 124</td>
          <td>                iCyn731</td>
          <td>          10.1186/1752-0509-7-142</td>
          <td>               Mueller</td>
          <td> 2013</td>
          <td> [sbml, excel]</td>
          <td>                           Cyanothece sp. PCC 7425</td>
          <td> Bacteria; Cyanobacteria; Oscillatoriophycideae...</td>
        </tr>
        <tr>
          <th>124</th>
          <td> 125</td>
          <td>                iCyj826</td>
          <td>          10.1186/1752-0509-7-142</td>
          <td>               Mueller</td>
          <td> 2013</td>
          <td> [sbml, excel]</td>
          <td>                           Cyanothece sp. PCC 7822</td>
          <td> Bacteria; Cyanobacteria; Oscillatoriophycideae...</td>
        </tr>
        <tr>
          <th>125</th>
          <td> 126</td>
          <td>                iCyp752</td>
          <td>          10.1186/1752-0509-7-142</td>
          <td>               Mueller</td>
          <td> 2013</td>
          <td> [sbml, excel]</td>
          <td>                           Cyanothece sp. PCC 8801</td>
          <td> Bacteria; Cyanobacteria; Oscillatoriophycideae...</td>
        </tr>
        <tr>
          <th>126</th>
          <td> 127</td>
          <td>                iCyh755</td>
          <td>          10.1186/1752-0509-7-142</td>
          <td>               Mueller</td>
          <td> 2013</td>
          <td> [sbml, excel]</td>
          <td>                           Cyanothece sp. PCC 8802</td>
          <td> Bacteria; Cyanobacteria; Oscillatoriophycideae...</td>
        </tr>
        <tr>
          <th>127</th>
          <td> 128</td>
          <td>                 iNF518</td>
          <td>        10.1007/s00253-013-5140-2</td>
          <td>               Flahaut</td>
          <td> 2013</td>
          <td>        [sbml]</td>
          <td>         Lactococcus lactis subsp. cremoris MG1363</td>
          <td> Bacteria; Firmicutes; Bacilli; Lactobacillales...</td>
        </tr>
        <tr>
          <th>128</th>
          <td> 129</td>
          <td>                iJL1454</td>
          <td>               10.1039/C3MB70090A</td>
          <td>               Jie Liu</td>
          <td> 2013</td>
          <td>       [excel]</td>
          <td>                       Aspergillus terreus NIH2624</td>
          <td> Eukaryota; Opisthokonta; Fungi; Dikarya; Ascom...</td>
        </tr>
        <tr>
          <th>129</th>
          <td> 130</td>
          <td>                iBif452</td>
          <td>           10.1186/1752-0509-8-41</td>
          <td>             El-Semman</td>
          <td> 2014</td>
          <td> [sbml, excel]</td>
          <td>                Bifidobacterium adolescentis L2-32</td>
          <td> Bacteria; Actinobacteria; Actinobacteria; Acti...</td>
        </tr>
        <tr>
          <th>130</th>
          <td> 131</td>
          <td>                iFap484</td>
          <td>           10.1186/1752-0509-8-41</td>
          <td>             El-Semman</td>
          <td> 2014</td>
          <td> [sbml, excel]</td>
          <td>               Faecalibacterium prausnitzii A2-165</td>
          <td> Bacteria; Firmicutes; Clostridia; Clostridiale...</td>
        </tr>
        <tr>
          <th>131</th>
          <td> 132</td>
          <td>                 iAM388</td>
          <td>         10.1186/1471-2164-12-535</td>
          <td>          Aline Metris</td>
          <td> 2011</td>
          <td>       [excel]</td>
          <td>     Campylobacter jejuni subsp. jejuni NCTC 11168</td>
          <td> Bacteria; Proteobacteria; delta/epsilon subdiv...</td>
        </tr>
        <tr>
          <th>132</th>
          <td> 133</td>
          <td>                 iTT548</td>
          <td>          10.1186/1475-2859-13-61</td>
          <td>            Na-Rae Lee</td>
          <td> 2014</td>
          <td> [sbml, excel]</td>
          <td>                              Thermus thermophilus</td>
          <td> Bacteria; Deinococcus-Thermus; Deinococci; The...</td>
        </tr>
        <tr>
          <th>133</th>
          <td> 134</td>
          <td>            EctoGEM-1.0</td>
          <td>                10.1111/tpj.12627</td>
          <td>               Prigent</td>
          <td> 2014</td>
          <td>        [sbml]</td>
          <td>                            Ectocarpus siliculosus</td>
          <td> Eukaryota; Stramenopiles; PX clade; Phaeophyce...</td>
        </tr>
        <tr>
          <th>134</th>
          <td> 135</td>
          <td>                 iMF721</td>
          <td>          10.1111/1462-2920.12513</td>
          <td>                 Fondi</td>
          <td> 2014</td>
          <td>        [sbml]</td>
          <td>             Pseudoalteromonas haloplanktis TAC125</td>
          <td> Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
        </tr>
        <tr>
          <th>135</th>
          <td> 136</td>
          <td> Arabidopsis core model</td>
          <td>            10.1104/pp.114.235358</td>
          <td>                Arnold</td>
          <td> 2014</td>
          <td>        [sbml]</td>
          <td>                              Arabidopsis thaliana</td>
          <td> Eukaryota; Viridiplantae; Streptophyta; Strept...</td>
        </tr>
        <tr>
          <th>136</th>
          <td> 137</td>
          <td>                 iHN637</td>
          <td>         10.1186/1475-2859-12-118</td>
          <td>      Harish Nagarajan</td>
          <td> 2013</td>
          <td>            []</td>
          <td>                           Clostridium ljungdahlii</td>
          <td> Bacteria; Firmicutes; Clostridia; Clostridiale...</td>
        </tr>
        <tr>
          <th>137</th>
          <td> 138</td>
          <td>                iCac802</td>
          <td>        10.1186/s13068-014-0144-4</td>
          <td>         Satyakam Dash</td>
          <td> 2014</td>
          <td>        [sbml]</td>
          <td>               Clostridium acetobutylicum ATCC 824</td>
          <td> Bacteria; Firmicutes; Clostridia; Clostridiale...</td>
        </tr>
        <tr>
          <th>138</th>
          <td> 139</td>
          <td>            iMLTC806cdf</td>
          <td>        10.1186/s12918-014-0117-z</td>
          <td>           M. Larocque</td>
          <td> 2014</td>
          <td> [sbml, excel]</td>
          <td>                         Clostridium difficile 630</td>
          <td> Bacteria; Firmicutes; Clostridia; Clostridiale...</td>
        </tr>
        <tr>
          <th>139</th>
          <td> 140</td>
          <td>                iCY1106</td>
          <td>        10.1186/s12918-014-0137-8</td>
          <td>               Chao Ye</td>
          <td> 2015</td>
          <td>        [sbml]</td>
          <td>                                Mortierella alpina</td>
          <td> Eukaryota; Opisthokonta; Fungi; Fungi incertae...</td>
        </tr>
        <tr>
          <th>140</th>
          <td> 141</td>
          <td>                 iMM518</td>
          <td>               10.1039/c3mb70421a</td>
          <td>              N. Goyal</td>
          <td> 2014</td>
          <td>       [excel]</td>
          <td>                      Methanococcus maripaludis S2</td>
          <td> Archaea; Euryarchaeota; Methanococci; Methanoc...</td>
        </tr>
        <tr>
          <th>141</th>
          <td> 142</td>
          <td>                iPC1209</td>
          <td>    10.1016/j.febslet.2014.12.010</td>
          <td>            Cheng Wang</td>
          <td> 2015</td>
          <td>            []</td>
          <td> Pectobacterium carotovorum subsp. carotovorum PC1</td>
          <td> Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
        </tr>
        <tr>
          <th>142</th>
          <td> 143</td>
          <td>                 iNV706</td>
          <td>             10.1128/AEM.03279-14</td>
          <td>              N. Veith</td>
          <td> 2014</td>
          <td>            []</td>
          <td>                        Enterococcus faecalis V583</td>
          <td> Bacteria; Firmicutes; Bacilli; Lactobacillales...</td>
        </tr>
        <tr>
          <th>143</th>
          <td> 144</td>
          <td>             KoxGSC1457</td>
          <td>          10.1186/1475-2859-12-20</td>
          <td>               J. Park</td>
          <td> 2013</td>
          <td>            []</td>
          <td>                                Klebsiella oxytoca</td>
          <td> Bacteria; Proteobacteria; Gammaproteobacteria;...</td>
        </tr>
      </tbody>
    </table>
    <p>144 rows × 8 columns</p>
    </div>




