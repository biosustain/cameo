
.. code:: python

    from IPython.display import display
    import re

Predict heterologous pathways
=============================

Predicting heterologous pathways is an important strategy to generate
new viable strains. Because portfolio of available reactions is very
large, computer assisted pathway design becomes essential. **Cameo**
implements a pathway search algorithm using an universal biochemical
reaction database that enumerates the shortest pathways.

.. raw:: html

   <div class="alert alert-warning">

If you're running this notebook on
`try.cameo.bio <http://try.cameo.bio>`__, things might run very slow due
to our inability to provide access to the proprietary
`CPLEX <https://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/>`__
solver on a public webserver. Furthermore, Jupyter kernels might crash
and restart due to memory limitations on the server.

.. raw:: html

   </div>

.. code:: python

    from cameo import models
    from cameo.strain_design import pathway_prediction



.. raw:: html

    




.. code:: python

    model = models.bigg.iMM904

.. code:: python

    predictor = pathway_prediction.PathwayPredictor(model=model, compartment_regexp=re.compile(".*_c$"))

.. code:: python

    pathways = predictor.run(product="vanillin", max_predictions=4)



.. raw:: html

    <span>Pathway 1</span>



.. raw:: html

    <div>
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
          <th>MNXR5336</th>
          <td>vanillin + H2O + NAD(+) &lt;=&gt; NADH(2-) + 2.0 H(+...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR5340</th>
          <td>formaldehyde + H2O + NAD(+) + 3,4-dihydroxyben...</td>
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



.. raw:: html

    <span>Pathway 2</span>



.. raw:: html

    <div>
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
          <th>MNXR230</th>
          <td>NADP(+) + H2O + 3,4-dihydroxybenzoate &lt;=&gt; O2 +...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR5336</th>
          <td>vanillin + H2O + NAD(+) &lt;=&gt; NADH(2-) + 2.0 H(+...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR5340</th>
          <td>formaldehyde + H2O + NAD(+) + 3,4-dihydroxyben...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
      </tbody>
    </table>
    </div>



.. raw:: html

    <span>Pathway 3</span>



.. raw:: html

    <div>
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
          <th>MNXR5336</th>
          <td>vanillin + H2O + NAD(+) &lt;=&gt; NADH(2-) + 2.0 H(+...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR5340</th>
          <td>formaldehyde + H2O + NAD(+) + 3,4-dihydroxyben...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR14769</th>
          <td>H2O + 3,4-dihydroxybenzoate + NAD(+) &lt;=&gt; NADH(...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
      </tbody>
    </table>
    </div>



.. raw:: html

    <span>Pathway 4</span>



.. raw:: html

    <div>
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
          <th>MNXR5336</th>
          <td>vanillin + H2O + NAD(+) &lt;=&gt; NADH(2-) + 2.0 H(+...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR5340</th>
          <td>formaldehyde + H2O + NAD(+) + 3,4-dihydroxyben...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR6101</th>
          <td>anthranilate + NADH(2-) + O2 + 3.0 H(+) &lt;=&gt; NH...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR7067</th>
          <td>3,4-dihydroxybenzoate + H(+) &lt;=&gt; catechol + CO(2)</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
      </tbody>
    </table>
    </div>


.. code:: python

    pathways.pathways[0].reactions[0]




.. raw:: html

    
            <table>
                <tr>
                    <td><strong>Id</strong></td><td>MNXR5336</td>
                </tr>
                <tr>
                    <td><strong>Name</strong></td><td>rhea:13309</td>
                </tr>
                <tr>
                    <td><strong>Stoichiometry</strong></td><td>MNXM754 + MNXM2 + MNXM8 <=> MNXM10 + 2.0 MNXM1 + MNXM982</td>
                </tr>
                <tr>
                    <td><strong>Lower bound</strong></td><td>-1000.000000</td>
                </tr>
                <tr>
                    <td><strong>Upper bound</strong></td><td>1000.000000</td>
                </tr>
            </table>
            



.. code:: python

    pathways.plot_production_envelopes(model, objective=model.reactions.BIOMASS_SC5_notrace)



.. raw:: html

    
    
        <div class="plotdiv" id="7e6285e3-e28e-48dd-8bb1-84277fe73b86"></div>
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
        };var element = document.getElementById("7e6285e3-e28e-48dd-8bb1-84277fe73b86");
        if (element == null) {
          console.log("Bokeh: ERROR: autoload.js configured with elementid '7e6285e3-e28e-48dd-8bb1-84277fe73b86' but no matching script tag was found. ")
          return false;
        }
      
        var js_urls = [];
      
        var inline_js = [
          function(Bokeh) {
            Bokeh.$(function() {
                var docs_json = {"44f7e98d-a65e-4a02-be5b-0131c21d0d70":{"roots":{"references":[{"attributes":{"callback":null},"id":"ecef4944-e951-47d6-9c6b-8dd0c286ce99","type":"DataRange1d"},{"attributes":{"plot":{"id":"f9334575-dce2-4399-a5b1-f8add3827aae","subtype":"Figure","type":"Plot"}},"id":"c69af69c-cf6a-4597-88e9-794eca854224","type":"WheelZoomTool"},{"attributes":{"callback":null},"id":"354ec457-b55d-452c-872c-dbb98a448482","type":"DataRange1d"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"x":{"field":"x"},"y":{"field":"y"}},"id":"ee7e88ac-28a4-4ccc-b070-c58085f9135c","type":"Line"},{"attributes":{"callback":null,"column_names":["y","x"],"data":{"x":[0.2878657037040178,0.27271487719328,0.25756405068254223,0.24241322417180447,0.2272623976610667,0.21211157115032891,0.19696074463959112,0.18180991812885336,0.1666590916181156,0.1515082651073778,0.13635743859664,0.12120661208590225,0.10605578557516446,0.0909049590644267,0.0757541325536889,0.06060330604295114,0.04545247953221335,0.030301653021475583,0.015150826510737792,0.0],"y":[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]}},"id":"e6fa0350-7819-43f4-8b7b-e6df0535dbc4","type":"ColumnDataSource"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"x":{"field":"x"},"y":{"field":"y"}},"id":"9a070b4b-eba8-4b93-90db-e0da56c1113d","type":"Line"},{"attributes":{"plot":{"id":"2e7b3ef2-cd9e-4a57-a834-7d0fd45f429d","subtype":"Figure","type":"Plot"}},"id":"389d7cd6-c4a1-4b49-8cdc-8e708c3ad1a7","type":"ResetTool"},{"attributes":{"axis_label":"BIOMASS_SC5_notrace","formatter":{"id":"6c53d353-3c29-4386-9c1b-f07fc6ccc32b","type":"BasicTickFormatter"},"plot":{"id":"d3054252-9bf3-4b4a-946d-ec21f1b67e14","subtype":"Figure","type":"Plot"},"ticker":{"id":"a0269c26-dc16-404b-b08e-6ee4d9b830c3","type":"BasicTicker"}},"id":"9476800b-d926-4273-89d1-4b9c22a70313","type":"LinearAxis"},{"attributes":{"data_source":{"id":"5f01d43f-eacb-41e7-b9eb-f0210823f29f","type":"ColumnDataSource"},"glyph":{"id":"550fac0e-96d8-4312-9e3e-fa4de2c0d6ea","type":"Patch"},"hover_glyph":null,"nonselection_glyph":{"id":"82eee6dc-d668-4203-b195-e70ac7be2856","type":"Patch"},"selection_glyph":null},"id":"2166f63b-071f-4c3f-90e1-257720611083","type":"GlyphRenderer"},{"attributes":{"callback":null,"column_names":["y","x"],"data":{"x":[0.0,0.0],"y":[1.9053254437869822,0.0]}},"id":"0b02c023-bca9-4fe4-8874-2c3251a04e1f","type":"ColumnDataSource"},{"attributes":{"plot":{"id":"d3054252-9bf3-4b4a-946d-ec21f1b67e14","subtype":"Figure","type":"Plot"}},"id":"88ba76b9-2a6b-4fdb-860f-3f4346590c26","type":"PreviewSaveTool"},{"attributes":{},"id":"31cde394-d34e-4b99-83c6-cc4e20e23791","type":"ToolEvents"},{"attributes":{"fill_alpha":{"value":0.1},"fill_color":{"value":"#1f77b4"},"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"x":{"field":"x"},"y":{"field":"y"}},"id":"2d7597e4-5b7b-442d-841d-2651845836a3","type":"Patch"},{"attributes":{"plot":{"id":"d3054252-9bf3-4b4a-946d-ec21f1b67e14","subtype":"Figure","type":"Plot"}},"id":"7ea42d9a-6a1d-4d3b-a948-5f80d8a14608","type":"HelpTool"},{"attributes":{"data_source":{"id":"834237a1-58c7-46c0-804c-64a6f6fd205f","type":"ColumnDataSource"},"glyph":{"id":"e10d1364-e882-429d-94b3-642a5c152bb5","type":"Line"},"hover_glyph":null,"nonselection_glyph":{"id":"f7395237-34d4-41ea-a7aa-e1df7150e431","type":"Line"},"selection_glyph":null},"id":"e518ce08-cb0c-456f-be5f-cf207f098c2c","type":"GlyphRenderer"},{"attributes":{"plot":{"id":"2e7b3ef2-cd9e-4a57-a834-7d0fd45f429d","subtype":"Figure","type":"Plot"}},"id":"a91df55d-7ffc-4a04-9e09-5b807a4bc5e7","type":"PreviewSaveTool"},{"attributes":{"line_color":{"value":"#B3E2CD"},"x":{"field":"x"},"y":{"field":"y"}},"id":"4f7a9bf0-09ee-4da4-92a2-e204096af5c0","type":"Line"},{"attributes":{"callback":null},"id":"d33aa947-4804-4d32-8db3-eccfc6f4eb6c","type":"DataRange1d"},{"attributes":{"callback":null},"id":"1bcea24a-f780-46a2-975b-28daf6bf0d28","type":"DataRange1d"},{"attributes":{"fill_alpha":{"value":0.1},"fill_color":{"value":"#1f77b4"},"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"x":{"field":"x"},"y":{"field":"y"}},"id":"bb308c6c-a0d7-4956-9972-292d75b05765","type":"Patch"},{"attributes":{"callback":null,"column_names":["y","x"],"data":{"x":[0.0,0.0],"y":[5.05751391465677,0.0]}},"id":"f37780b3-d40b-4716-995a-fc8c2d8b8489","type":"ColumnDataSource"},{"attributes":{"below":[{"id":"9476800b-d926-4273-89d1-4b9c22a70313","type":"LinearAxis"}],"left":[{"id":"8dc8eb70-230e-4790-9ec3-763290941f98","type":"LinearAxis"}],"plot_height":278,"plot_width":450,"renderers":[{"id":"9476800b-d926-4273-89d1-4b9c22a70313","type":"LinearAxis"},{"id":"d0d0a844-36e8-4d20-a34e-2b043f77d4b6","type":"Grid"},{"id":"8dc8eb70-230e-4790-9ec3-763290941f98","type":"LinearAxis"},{"id":"4f01b310-2e1f-462b-8f32-9e9f96c613a5","type":"Grid"},{"id":"e834602c-32e4-4f38-8ccc-9ee6a4f185f9","type":"BoxAnnotation"},{"id":"da51939b-a9be-4795-86d8-cff5fd8b3a81","type":"Legend"},{"id":"2166f63b-071f-4c3f-90e1-257720611083","type":"GlyphRenderer"},{"id":"868f0c29-b3b0-4780-be91-71bc68ec8c07","type":"GlyphRenderer"},{"id":"2d73811b-4506-4047-8ef5-33df7bf720d2","type":"GlyphRenderer"},{"id":"b15aeabc-8a1f-4fb3-b1ce-8b62f406c792","type":"GlyphRenderer"}],"title":"Pathway 3","tool_events":{"id":"31cde394-d34e-4b99-83c6-cc4e20e23791","type":"ToolEvents"},"tools":[{"id":"4b7fe162-b870-4385-ab6e-6c03b5b0a941","type":"PanTool"},{"id":"def23848-bd59-4215-831c-3216c3638295","type":"WheelZoomTool"},{"id":"a323b27b-d527-4e63-a6ca-0d82e3700efe","type":"BoxZoomTool"},{"id":"88ba76b9-2a6b-4fdb-860f-3f4346590c26","type":"PreviewSaveTool"},{"id":"0bd41cd9-3ec1-4b26-96a2-4a385db9185d","type":"ResizeTool"},{"id":"bcd4f175-eb5c-4873-bc66-3976362abbd4","type":"ResetTool"},{"id":"7ea42d9a-6a1d-4d3b-a948-5f80d8a14608","type":"HelpTool"}],"x_range":{"id":"1bcea24a-f780-46a2-975b-28daf6bf0d28","type":"DataRange1d"},"y_range":{"id":"ecef4944-e951-47d6-9c6b-8dd0c286ce99","type":"DataRange1d"}},"id":"d3054252-9bf3-4b4a-946d-ec21f1b67e14","subtype":"Figure","type":"Plot"},{"attributes":{"callback":null,"column_names":["y","x"],"data":{"x":[0.28786894488710274,0.2727179477877815,0.25756695068846036,0.24241595358913914,0.22726495648981795,0.21211395939049676,0.19696296229117555,0.18181196519185436,0.16666096809253317,0.15150997099321198,0.1363589738938908,0.12120797679456957,0.10605697969524838,0.09090598259592719,0.07575498549660598,0.060603988397284786,0.045452991297963596,0.03030199419864238,0.015150997099321217,0.0],"y":[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]}},"id":"18cb1b3a-24ea-4aa6-805c-707f7908a3b7","type":"ColumnDataSource"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"x":{"field":"x"},"y":{"field":"y"}},"id":"e7d46aee-548b-49d6-bcb6-75718d0531f1","type":"Line"},{"attributes":{"dimension":1,"plot":{"id":"25762579-4c02-4e7a-ba6c-eaceea6d3a7a","subtype":"Figure","type":"Plot"},"ticker":{"id":"7e71eadf-4e0c-49aa-bb38-4401dd36c094","type":"BasicTicker"}},"id":"36fb3f0a-5d41-411e-8067-5776dd6fdfd3","type":"Grid"},{"attributes":{"line_color":{"value":"#B3E2CD"},"x":{"field":"x"},"y":{"field":"y"}},"id":"4b107b7c-d69e-41b8-a424-afd922409d37","type":"Line"},{"attributes":{"data_source":{"id":"47606863-2b03-45d9-9c46-0cb0af003d1e","type":"ColumnDataSource"},"glyph":{"id":"4b107b7c-d69e-41b8-a424-afd922409d37","type":"Line"},"hover_glyph":null,"nonselection_glyph":{"id":"5a2d41ba-f33a-4c7c-8e9f-7686181b297d","type":"Line"},"selection_glyph":null},"id":"868f0c29-b3b0-4780-be91-71bc68ec8c07","type":"GlyphRenderer"},{"attributes":{},"id":"9fc8f2ea-3ce2-42d4-9850-2434e471c2bd","type":"BasicTicker"},{"attributes":{"legends":[["WT",[{"id":"cdcf2bae-4e98-4622-b235-f0bed5e8b828","type":"GlyphRenderer"}]]],"plot":{"id":"25762579-4c02-4e7a-ba6c-eaceea6d3a7a","subtype":"Figure","type":"Plot"}},"id":"0a0d9308-899d-4070-8b8a-17eb3d754b46","type":"Legend"},{"attributes":{"data_source":{"id":"c51cdf46-cfe3-40ab-8886-0c817ed48c9c","type":"ColumnDataSource"},"glyph":{"id":"905b4f3f-7acf-4f70-8002-d328d6cb235e","type":"Line"},"hover_glyph":null,"nonselection_glyph":{"id":"42744764-61ac-4a3d-aec5-a1d2114b6ad3","type":"Line"},"selection_glyph":null},"id":"67eac2d3-d959-454d-9343-077dcd8cc4cb","type":"GlyphRenderer"},{"attributes":{"overlay":{"id":"c6293838-84ee-4e6f-91d3-6dd92e8e80dc","type":"BoxAnnotation"},"plot":{"id":"2e7b3ef2-cd9e-4a57-a834-7d0fd45f429d","subtype":"Figure","type":"Plot"}},"id":"91cfcf31-a124-45d6-acdd-63f0b7f02ed0","type":"BoxZoomTool"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"x":{"field":"x"},"y":{"field":"y"}},"id":"f7395237-34d4-41ea-a7aa-e1df7150e431","type":"Line"},{"attributes":{},"id":"4d28dbd9-f11b-4adc-98e5-e2664da60062","type":"BasicTickFormatter"},{"attributes":{"axis_label":"DM_MNXM754","formatter":{"id":"dfce2671-b2f8-4829-9f3a-1215547e04bd","type":"BasicTickFormatter"},"plot":{"id":"d3054252-9bf3-4b4a-946d-ec21f1b67e14","subtype":"Figure","type":"Plot"},"ticker":{"id":"46b88b7a-4b01-4a5d-b0c7-ab6a8b177b9d","type":"BasicTicker"}},"id":"8dc8eb70-230e-4790-9ec3-763290941f98","type":"LinearAxis"},{"attributes":{"bottom_units":"screen","fill_alpha":{"value":0.5},"fill_color":{"value":"lightgrey"},"left_units":"screen","level":"overlay","line_alpha":{"value":1.0},"line_color":{"value":"black"},"line_dash":[4,4],"line_width":{"value":2},"plot":null,"render_mode":"css","right_units":"screen","top_units":"screen"},"id":"ff6e0877-dee3-4a08-a193-576e63e69290","type":"BoxAnnotation"},{"attributes":{"axis_label":"DM_MNXM754","formatter":{"id":"aa771144-c208-4dbb-8bba-40dbc10260a0","type":"BasicTickFormatter"},"plot":{"id":"2e7b3ef2-cd9e-4a57-a834-7d0fd45f429d","subtype":"Figure","type":"Plot"},"ticker":{"id":"e0764d81-f5e9-48ef-9015-9fb45681d606","type":"BasicTicker"}},"id":"68829b21-2f64-4776-a0f9-9406c3d0cf67","type":"LinearAxis"},{"attributes":{"plot":{"id":"f9334575-dce2-4399-a5b1-f8add3827aae","subtype":"Figure","type":"Plot"}},"id":"3df98c73-4108-4c21-a533-c3331dd3b4d3","type":"ResizeTool"},{"attributes":{"line_color":{"value":"#B3E2CD"},"x":{"field":"x"},"y":{"field":"y"}},"id":"b8741388-c46e-42ae-b878-ded478136ba8","type":"Line"},{"attributes":{},"id":"50431d28-03ab-4d1d-8d82-c8f8acdf35fa","type":"ToolEvents"},{"attributes":{"plot":{"id":"f9334575-dce2-4399-a5b1-f8add3827aae","subtype":"Figure","type":"Plot"},"ticker":{"id":"b8e042dd-7420-449b-aec5-29b9dadf52d9","type":"BasicTicker"}},"id":"13a65df7-29e3-4f7b-baba-1cf09009e5a1","type":"Grid"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"x":{"field":"x"},"y":{"field":"y"}},"id":"5a2d41ba-f33a-4c7c-8e9f-7686181b297d","type":"Line"},{"attributes":{"callback":null,"column_names":["y","x"],"data":{"x":[0.0,0.0],"y":[2.697588126159554,0.0]}},"id":"d9e3490f-2df1-4add-af42-581002819eca","type":"ColumnDataSource"},{"attributes":{"plot":{"id":"2e7b3ef2-cd9e-4a57-a834-7d0fd45f429d","subtype":"Figure","type":"Plot"}},"id":"d8703d02-f497-46ad-802f-0bfbd09a5339","type":"HelpTool"},{"attributes":{"overlay":{"id":"18dd9651-0a3b-4464-81e0-4803d3853f34","type":"BoxAnnotation"},"plot":{"id":"f9334575-dce2-4399-a5b1-f8add3827aae","subtype":"Figure","type":"Plot"}},"id":"1e0ad84c-69ab-410c-b4fa-ba0926d73f31","type":"BoxZoomTool"},{"attributes":{"bottom_units":"screen","fill_alpha":{"value":0.5},"fill_color":{"value":"lightgrey"},"left_units":"screen","level":"overlay","line_alpha":{"value":1.0},"line_color":{"value":"black"},"line_dash":[4,4],"line_width":{"value":2},"plot":null,"render_mode":"css","right_units":"screen","top_units":"screen"},"id":"e834602c-32e4-4f38-8ccc-9ee6a4f185f9","type":"BoxAnnotation"},{"attributes":{"data_source":{"id":"18cb1b3a-24ea-4aa6-805c-707f7908a3b7","type":"ColumnDataSource"},"glyph":{"id":"e1eb2083-a83f-45ab-a718-8f21722250de","type":"Line"},"hover_glyph":null,"nonselection_glyph":{"id":"ee7e88ac-28a4-4ccc-b070-c58085f9135c","type":"Line"},"selection_glyph":null},"id":"2d73811b-4506-4047-8ef5-33df7bf720d2","type":"GlyphRenderer"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"x":{"field":"x"},"y":{"field":"y"}},"id":"bbe07337-d1b5-41ae-b4b0-054ed337081c","type":"Line"},{"attributes":{"axis_label":"DM_MNXM754","formatter":{"id":"89d19c4f-affa-4084-990b-245765d15ee5","type":"BasicTickFormatter"},"plot":{"id":"f9334575-dce2-4399-a5b1-f8add3827aae","subtype":"Figure","type":"Plot"},"ticker":{"id":"ce3197d2-2e84-4271-8740-c1a6cdc50cbb","type":"BasicTicker"}},"id":"4f682599-da23-42de-be91-bde5fa43c724","type":"LinearAxis"},{"attributes":{},"id":"ce3197d2-2e84-4271-8740-c1a6cdc50cbb","type":"BasicTicker"},{"attributes":{"callback":null,"column_names":["y","x"],"data":{"x":[0.2878657037040163,0.27271487719327864,0.2575640506825409,0.24241322417180322,0.2272623976610655,0.2121115711503278,0.19696074463959012,0.18180991812885242,0.1666590916181147,0.151508265107377,0.1363574385966393,0.12120661208590161,0.1060557855751639,0.0909049590644262,0.07575413255368851,0.060603306042950805,0.0454524795322131,0.030301653021475417,0.01515082651073768,0.0],"y":[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]}},"id":"c51cdf46-cfe3-40ab-8886-0c817ed48c9c","type":"ColumnDataSource"},{"attributes":{},"id":"7e71eadf-4e0c-49aa-bb38-4401dd36c094","type":"BasicTicker"},{"attributes":{"axis_label":"BIOMASS_SC5_notrace","formatter":{"id":"f3ff280e-8a98-4320-9797-87f66db9d076","type":"BasicTickFormatter"},"plot":{"id":"25762579-4c02-4e7a-ba6c-eaceea6d3a7a","subtype":"Figure","type":"Plot"},"ticker":{"id":"9fc8f2ea-3ce2-42d4-9850-2434e471c2bd","type":"BasicTicker"}},"id":"8a2fb091-1ae8-4977-b855-b55fd0f83cee","type":"LinearAxis"},{"attributes":{"below":[{"id":"9535c367-8fcf-4c46-b90f-df6058734517","type":"LinearAxis"}],"left":[{"id":"68829b21-2f64-4776-a0f9-9406c3d0cf67","type":"LinearAxis"}],"plot_height":278,"plot_width":450,"renderers":[{"id":"9535c367-8fcf-4c46-b90f-df6058734517","type":"LinearAxis"},{"id":"317417ba-dd06-4070-bb81-b50a69e853d8","type":"Grid"},{"id":"68829b21-2f64-4776-a0f9-9406c3d0cf67","type":"LinearAxis"},{"id":"572d704a-4121-47bd-9359-91537244347a","type":"Grid"},{"id":"c6293838-84ee-4e6f-91d3-6dd92e8e80dc","type":"BoxAnnotation"},{"id":"bd9f14ec-9396-4ee4-81a6-fcf48437a905","type":"Legend"},{"id":"5e50fc34-66cf-4e9e-94d3-ed458deafc54","type":"GlyphRenderer"},{"id":"18eee7a6-a67d-4825-b407-802e4cb365cf","type":"GlyphRenderer"},{"id":"6d61e318-a476-4643-8389-cec4948df1b5","type":"GlyphRenderer"},{"id":"071556e7-d41e-408d-95a2-ed1ca693e72e","type":"GlyphRenderer"}],"title":"Pathway 0","tool_events":{"id":"50431d28-03ab-4d1d-8d82-c8f8acdf35fa","type":"ToolEvents"},"tools":[{"id":"2e95323a-3f8b-4751-8246-774696b64549","type":"PanTool"},{"id":"df61542c-f583-47f6-a171-a056f6fd0f17","type":"WheelZoomTool"},{"id":"91cfcf31-a124-45d6-acdd-63f0b7f02ed0","type":"BoxZoomTool"},{"id":"a91df55d-7ffc-4a04-9e09-5b807a4bc5e7","type":"PreviewSaveTool"},{"id":"6d6f68b8-aba7-49da-bf6d-622f77c5d295","type":"ResizeTool"},{"id":"389d7cd6-c4a1-4b49-8cdc-8e708c3ad1a7","type":"ResetTool"},{"id":"d8703d02-f497-46ad-802f-0bfbd09a5339","type":"HelpTool"}],"x_range":{"id":"2238db27-74f1-4f8b-98ae-b66e55770e88","type":"DataRange1d"},"y_range":{"id":"354ec457-b55d-452c-872c-dbb98a448482","type":"DataRange1d"}},"id":"2e7b3ef2-cd9e-4a57-a834-7d0fd45f429d","subtype":"Figure","type":"Plot"},{"attributes":{"axis_label":"BIOMASS_SC5_notrace","formatter":{"id":"a2668286-a7c3-446b-8031-35edb7c3b65f","type":"BasicTickFormatter"},"plot":{"id":"2e7b3ef2-cd9e-4a57-a834-7d0fd45f429d","subtype":"Figure","type":"Plot"},"ticker":{"id":"6a69fad3-7600-4f97-bb5b-375827f6ba5c","type":"BasicTicker"}},"id":"9535c367-8fcf-4c46-b90f-df6058734517","type":"LinearAxis"},{"attributes":{"data_source":{"id":"dc1e2d82-40b8-4dbe-aa1f-30e3bffda4b1","type":"ColumnDataSource"},"glyph":{"id":"e4ce0ca9-6af3-401c-a237-2eb95d6ec131","type":"Patch"},"hover_glyph":null,"nonselection_glyph":{"id":"bb308c6c-a0d7-4956-9972-292d75b05765","type":"Patch"},"selection_glyph":null},"id":"cdcf2bae-4e98-4622-b235-f0bed5e8b828","type":"GlyphRenderer"},{"attributes":{},"id":"6a69fad3-7600-4f97-bb5b-375827f6ba5c","type":"BasicTicker"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"x":{"field":"x"},"y":{"field":"y"}},"id":"24ef094f-3125-41b9-8260-e17117bdde54","type":"Line"},{"attributes":{"callback":null},"id":"7185c3f2-9054-47fa-9337-6898ebc3466b","type":"DataRange1d"},{"attributes":{"children":[[{"id":"2e7b3ef2-cd9e-4a57-a834-7d0fd45f429d","subtype":"Figure","type":"Plot"},{"id":"f9334575-dce2-4399-a5b1-f8add3827aae","subtype":"Figure","type":"Plot"}],[{"id":"25762579-4c02-4e7a-ba6c-eaceea6d3a7a","subtype":"Figure","type":"Plot"},{"id":"d3054252-9bf3-4b4a-946d-ec21f1b67e14","subtype":"Figure","type":"Plot"}]],"name":"Production envelops for Demand vanillin"},"id":"04f9bfaa-d2b6-4cb3-ba36-939ef9ab0715","type":"GridPlot"},{"attributes":{"plot":{"id":"d3054252-9bf3-4b4a-946d-ec21f1b67e14","subtype":"Figure","type":"Plot"}},"id":"bcd4f175-eb5c-4873-bc66-3976362abbd4","type":"ResetTool"},{"attributes":{"dimension":1,"plot":{"id":"2e7b3ef2-cd9e-4a57-a834-7d0fd45f429d","subtype":"Figure","type":"Plot"},"ticker":{"id":"e0764d81-f5e9-48ef-9015-9fb45681d606","type":"BasicTicker"}},"id":"572d704a-4121-47bd-9359-91537244347a","type":"Grid"},{"attributes":{"data_source":{"id":"e6fa0350-7819-43f4-8b7b-e6df0535dbc4","type":"ColumnDataSource"},"glyph":{"id":"f6173ad4-d77a-46b5-9456-2fccc789e618","type":"Line"},"hover_glyph":null,"nonselection_glyph":{"id":"24ef094f-3125-41b9-8260-e17117bdde54","type":"Line"},"selection_glyph":null},"id":"cf3aeb0f-6b81-4716-98f6-190cfb6e1bba","type":"GlyphRenderer"},{"attributes":{"plot":{"id":"d3054252-9bf3-4b4a-946d-ec21f1b67e14","subtype":"Figure","type":"Plot"}},"id":"4b7fe162-b870-4385-ab6e-6c03b5b0a941","type":"PanTool"},{"attributes":{},"id":"a2668286-a7c3-446b-8031-35edb7c3b65f","type":"BasicTickFormatter"},{"attributes":{"line_color":{"value":"#B3E2CD"},"x":{"field":"x"},"y":{"field":"y"}},"id":"9ba7116c-79c7-4087-b09f-44ee38523a46","type":"Line"},{"attributes":{"data_source":{"id":"a3125ab1-dca0-47be-a84a-3a2b02f28ddf","type":"ColumnDataSource"},"glyph":{"id":"532b3024-5ab3-4449-8c7e-eaeef4969ff6","type":"Line"},"hover_glyph":null,"nonselection_glyph":{"id":"2bf7ed1d-581f-4c3f-ae41-b10e2e8e86b6","type":"Line"},"selection_glyph":null},"id":"6b55213d-52c7-45b3-8139-3692b97c4b30","type":"GlyphRenderer"},{"attributes":{},"id":"2bec14f8-77c4-45c4-8768-b610dffcfdbe","type":"BasicTickFormatter"},{"attributes":{"plot":{"id":"25762579-4c02-4e7a-ba6c-eaceea6d3a7a","subtype":"Figure","type":"Plot"}},"id":"2d3fcec0-6418-454e-b883-6579bc7957fe","type":"WheelZoomTool"},{"attributes":{"line_color":{"value":"#B3E2CD"},"x":{"field":"x"},"y":{"field":"y"}},"id":"f6173ad4-d77a-46b5-9456-2fccc789e618","type":"Line"},{"attributes":{"data_source":{"id":"aefd52e4-8476-40b0-bc58-80701af7bb30","type":"ColumnDataSource"},"glyph":{"id":"6a7ae084-5cbc-42c9-a60a-a58f3dcbf73f","type":"Line"},"hover_glyph":null,"nonselection_glyph":{"id":"bbe07337-d1b5-41ae-b4b0-054ed337081c","type":"Line"},"selection_glyph":null},"id":"63deb866-5517-4f3f-a443-93bb46088240","type":"GlyphRenderer"},{"attributes":{"callback":null,"column_names":["y","x"],"data":{"x":[0.2878657037040178,0.27271487719328,0.25756405068254223,0.24241322417180447,0.2272623976610667,0.21211157115032891,0.19696074463959112,0.18180991812885336,0.1666590916181156,0.1515082651073778,0.13635743859664,0.12120661208590225,0.10605578557516446,0.0909049590644267,0.0757541325536889,0.06060330604295114,0.04545247953221335,0.030301653021475583,0.015150826510737792,0.0],"y":[-3.817527732558974e-14,0.16023362324817922,0.319970122258177,0.4797066212681532,0.6394431202781523,0.7991796192881254,0.9589161182981297,1.1186526173081195,1.2783891163181074,1.438125615328088,1.5978621143380762,1.7575986133480708,1.9173351123580575,2.077071611368046,2.2368081103780315,2.396544609388022,2.556281108398011,2.7160176074079976,2.8757541064179857,3.0354906054279764]}},"id":"834237a1-58c7-46c0-804c-64a6f6fd205f","type":"ColumnDataSource"},{"attributes":{"legends":[["WT",[{"id":"1c471e33-5a32-450a-984a-f34d53e7c8e4","type":"GlyphRenderer"}]]],"plot":{"id":"f9334575-dce2-4399-a5b1-f8add3827aae","subtype":"Figure","type":"Plot"}},"id":"0cecdc58-cb67-45fb-b69a-27e447a5021d","type":"Legend"},{"attributes":{"line_color":{"value":"#B3E2CD"},"x":{"field":"x"},"y":{"field":"y"}},"id":"e1eb2083-a83f-45ab-a718-8f21722250de","type":"Line"},{"attributes":{},"id":"e0764d81-f5e9-48ef-9015-9fb45681d606","type":"BasicTicker"},{"attributes":{"bottom_units":"screen","fill_alpha":{"value":0.5},"fill_color":{"value":"lightgrey"},"left_units":"screen","level":"overlay","line_alpha":{"value":1.0},"line_color":{"value":"black"},"line_dash":[4,4],"line_width":{"value":2},"plot":null,"render_mode":"css","right_units":"screen","top_units":"screen"},"id":"c6293838-84ee-4e6f-91d3-6dd92e8e80dc","type":"BoxAnnotation"},{"attributes":{"fill_alpha":{"value":0.3},"fill_color":{"value":"#B3E2CD"},"line_color":{"value":"#1f77b4"},"x":{"field":"x"},"y":{"field":"y"}},"id":"3d63e360-3183-41e1-bed6-676e6f04c7ae","type":"Patch"},{"attributes":{"fill_alpha":{"value":0.3},"fill_color":{"value":"#B3E2CD"},"line_color":{"value":"#1f77b4"},"x":{"field":"x"},"y":{"field":"y"}},"id":"550fac0e-96d8-4312-9e3e-fa4de2c0d6ea","type":"Patch"},{"attributes":{"plot":{"id":"f9334575-dce2-4399-a5b1-f8add3827aae","subtype":"Figure","type":"Plot"}},"id":"b6b1577f-a475-4b1e-b11d-6b9113f475fe","type":"ResetTool"},{"attributes":{"plot":{"id":"25762579-4c02-4e7a-ba6c-eaceea6d3a7a","subtype":"Figure","type":"Plot"}},"id":"bf6d30a3-2b12-498c-9513-0d46073b150e","type":"HelpTool"},{"attributes":{},"id":"f536e755-41d3-4eef-94b8-f6a39df42f9e","type":"ToolEvents"},{"attributes":{},"id":"89d19c4f-affa-4084-990b-245765d15ee5","type":"BasicTickFormatter"},{"attributes":{"callback":null,"column_names":["y","x"],"data":{"x":[0.2878657037040163,0.27271487719327864,0.2575640506825409,0.24241322417180322,0.2272623976610655,0.2121115711503278,0.19696074463959012,0.18180991812885242,0.1666590916181147,0.151508265107377,0.1363574385966393,0.12120661208590161,0.1060557855751639,0.0909049590644262,0.07575413255368851,0.060603306042950805,0.0454524795322131,0.030301653021475417,0.01515082651073768,0.0,0.0,0.01515082651073768,0.030301653021475417,0.0454524795322131,0.060603306042950805,0.07575413255368851,0.0909049590644262,0.1060557855751639,0.12120661208590161,0.1363574385966393,0.151508265107377,0.1666590916181147,0.18180991812885242,0.19696074463959012,0.2121115711503278,0.2272623976610655,0.24241322417180322,0.2575640506825409,0.27271487719327864,0.2878657037040163],"y":[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.9053254437869822,1.808092191580202,1.710858939373418,1.6136256871666381,1.5163924349598485,1.4191591827530734,1.3219259305462905,1.2246926783395076,1.1274594261327182,1.0302261739259404,0.9329929217191623,0.8354770791227152,0.7342299103776272,0.6299031297951706,0.5253460487659517,0.420788967736725,0.3160236519224981,0.210946833464358,0.10587001500621629,2.9290963817141727e-15]}},"id":"f98e07d4-ee87-4b6d-ad52-9615d12539c9","type":"ColumnDataSource"},{"attributes":{},"id":"aa771144-c208-4dbb-8bba-40dbc10260a0","type":"BasicTickFormatter"},{"attributes":{},"id":"6c53d353-3c29-4386-9c1b-f07fc6ccc32b","type":"BasicTickFormatter"},{"attributes":{"bottom_units":"screen","fill_alpha":{"value":0.5},"fill_color":{"value":"lightgrey"},"left_units":"screen","level":"overlay","line_alpha":{"value":1.0},"line_color":{"value":"black"},"line_dash":[4,4],"line_width":{"value":2},"plot":null,"render_mode":"css","right_units":"screen","top_units":"screen"},"id":"18dd9651-0a3b-4464-81e0-4803d3853f34","type":"BoxAnnotation"},{"attributes":{"callback":null,"column_names":["y","x"],"data":{"x":[0.28786894488710274,0.2727179477877815,0.25756695068846036,0.24241595358913914,0.22726495648981795,0.21211395939049676,0.19696296229117555,0.18181196519185436,0.16666096809253317,0.15150997099321198,0.1363589738938908,0.12120797679456957,0.10605697969524838,0.09090598259592719,0.07575498549660598,0.060603988397284786,0.045452991297963596,0.03030199419864238,0.015150997099321217,0.0,0.0,0.015150997099321217,0.03030199419864238,0.045452991297963596,0.060603988397284786,0.07575498549660598,0.09090598259592719,0.10605697969524838,0.12120797679456957,0.1363589738938908,0.15150997099321198,0.16666096809253317,0.18181196519185436,0.19696296229117555,0.21211395939049676,0.22726495648981795,0.24241595358913914,0.25756695068846036,0.2727179477877815,0.28786894488710274],"y":[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,2.697588126159554,2.5556314572877974,2.4136747884160377,2.271718119544287,2.129761450672519,1.9878047818007563,1.8458481129290167,1.7038914440572328,1.561934775185485,1.4199781063137196,1.2780214374419605,1.13606476857021,0.9941080996984444,0.8521514308266925,0.7101947619549311,0.5682380930831765,0.4262814242114241,0.28432475533965035,0.14236808646788648,4.6379102849828174e-15]}},"id":"5f01d43f-eacb-41e7-b9eb-f0210823f29f","type":"ColumnDataSource"},{"attributes":{"plot":{"id":"25762579-4c02-4e7a-ba6c-eaceea6d3a7a","subtype":"Figure","type":"Plot"}},"id":"b78302ca-9e8f-41ef-84e2-7940499738c0","type":"PreviewSaveTool"},{"attributes":{"below":[{"id":"edd685a3-6c8b-42d1-9256-44b4f925f375","type":"LinearAxis"}],"left":[{"id":"4f682599-da23-42de-be91-bde5fa43c724","type":"LinearAxis"}],"plot_height":278,"plot_width":450,"renderers":[{"id":"edd685a3-6c8b-42d1-9256-44b4f925f375","type":"LinearAxis"},{"id":"13a65df7-29e3-4f7b-baba-1cf09009e5a1","type":"Grid"},{"id":"4f682599-da23-42de-be91-bde5fa43c724","type":"LinearAxis"},{"id":"75d398a1-1343-48f3-9dff-5f0276c19758","type":"Grid"},{"id":"18dd9651-0a3b-4464-81e0-4803d3853f34","type":"BoxAnnotation"},{"id":"0cecdc58-cb67-45fb-b69a-27e447a5021d","type":"Legend"},{"id":"1c471e33-5a32-450a-984a-f34d53e7c8e4","type":"GlyphRenderer"},{"id":"6b55213d-52c7-45b3-8139-3692b97c4b30","type":"GlyphRenderer"},{"id":"67eac2d3-d959-454d-9343-077dcd8cc4cb","type":"GlyphRenderer"},{"id":"b8696d98-e479-4edf-b150-e52ca7f838c3","type":"GlyphRenderer"}],"title":"Pathway 1","tool_events":{"id":"f536e755-41d3-4eef-94b8-f6a39df42f9e","type":"ToolEvents"},"tools":[{"id":"fd474d26-da3c-4ce3-93bc-130e4b3b2380","type":"PanTool"},{"id":"c69af69c-cf6a-4597-88e9-794eca854224","type":"WheelZoomTool"},{"id":"1e0ad84c-69ab-410c-b4fa-ba0926d73f31","type":"BoxZoomTool"},{"id":"815aef33-d194-4873-a26c-3ae0f61fa288","type":"PreviewSaveTool"},{"id":"3df98c73-4108-4c21-a533-c3331dd3b4d3","type":"ResizeTool"},{"id":"b6b1577f-a475-4b1e-b11d-6b9113f475fe","type":"ResetTool"},{"id":"786d0d39-5327-4ea4-b4cd-df9bf243c910","type":"HelpTool"}],"x_range":{"id":"d33aa947-4804-4d32-8db3-eccfc6f4eb6c","type":"DataRange1d"},"y_range":{"id":"7185c3f2-9054-47fa-9337-6898ebc3466b","type":"DataRange1d"}},"id":"f9334575-dce2-4399-a5b1-f8add3827aae","subtype":"Figure","type":"Plot"},{"attributes":{"data_source":{"id":"f37780b3-d40b-4716-995a-fc8c2d8b8489","type":"ColumnDataSource"},"glyph":{"id":"4f7a9bf0-09ee-4da4-92a2-e204096af5c0","type":"Line"},"hover_glyph":null,"nonselection_glyph":{"id":"995e6a45-6dae-4e14-9d4b-718213625e61","type":"Line"},"selection_glyph":null},"id":"071556e7-d41e-408d-95a2-ed1ca693e72e","type":"GlyphRenderer"},{"attributes":{"legends":[["WT",[{"id":"5e50fc34-66cf-4e9e-94d3-ed458deafc54","type":"GlyphRenderer"}]]],"plot":{"id":"2e7b3ef2-cd9e-4a57-a834-7d0fd45f429d","subtype":"Figure","type":"Plot"}},"id":"bd9f14ec-9396-4ee4-81a6-fcf48437a905","type":"Legend"},{"attributes":{"line_color":{"value":"#B3E2CD"},"x":{"field":"x"},"y":{"field":"y"}},"id":"5b386884-5a8b-407a-8c59-32474910c21d","type":"Line"},{"attributes":{"fill_alpha":{"value":0.3},"fill_color":{"value":"#B3E2CD"},"line_color":{"value":"#1f77b4"},"x":{"field":"x"},"y":{"field":"y"}},"id":"e4ce0ca9-6af3-401c-a237-2eb95d6ec131","type":"Patch"},{"attributes":{"callback":null,"column_names":["y","x"],"data":{"x":[0.0,0.0],"y":[3.0354906054279764,0.0]}},"id":"aefd52e4-8476-40b0-bc58-80701af7bb30","type":"ColumnDataSource"},{"attributes":{"overlay":{"id":"ff6e0877-dee3-4a08-a193-576e63e69290","type":"BoxAnnotation"},"plot":{"id":"25762579-4c02-4e7a-ba6c-eaceea6d3a7a","subtype":"Figure","type":"Plot"}},"id":"5a1ebe83-a671-4a57-9295-2eec5a5ab58d","type":"BoxZoomTool"},{"attributes":{"line_color":{"value":"#B3E2CD"},"x":{"field":"x"},"y":{"field":"y"}},"id":"532b3024-5ab3-4449-8c7e-eaeef4969ff6","type":"Line"},{"attributes":{"callback":null},"id":"f0e30ec5-0283-41dd-8dda-12cafc521bb2","type":"DataRange1d"},{"attributes":{},"id":"b8e042dd-7420-449b-aec5-29b9dadf52d9","type":"BasicTicker"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"x":{"field":"x"},"y":{"field":"y"}},"id":"995e6a45-6dae-4e14-9d4b-718213625e61","type":"Line"},{"attributes":{"callback":null,"column_names":["y","x"],"data":{"x":[0.28786894488710274,0.2727179477877815,0.25756695068846036,0.24241595358913914,0.22726495648981795,0.21211395939049676,0.19696296229117555,0.18181196519185436,0.16666096809253317,0.15150997099321198,0.1363589738938908,0.12120797679456957,0.10605697969524838,0.09090598259592719,0.07575498549660598,0.060603988397284786,0.045452991297963596,0.03030199419864238,0.015150997099321217,0.0],"y":[4.6379102849828174e-15,0.14236808646788648,0.28432475533965035,0.4262814242114241,0.5682380930831765,0.7101947619549311,0.8521514308266925,0.9941080996984444,1.13606476857021,1.2780214374419605,1.4199781063137196,1.561934775185485,1.7038914440572328,1.8458481129290167,1.9878047818007563,2.129761450672519,2.271718119544287,2.4136747884160377,2.5556314572877974,2.697588126159554]}},"id":"47606863-2b03-45d9-9c46-0cb0af003d1e","type":"ColumnDataSource"},{"attributes":{"data_source":{"id":"f98e07d4-ee87-4b6d-ad52-9615d12539c9","type":"ColumnDataSource"},"glyph":{"id":"af44b093-56c9-4ab7-8ddf-a8d628cd3c97","type":"Patch"},"hover_glyph":null,"nonselection_glyph":{"id":"59b9c2e1-ee04-4a3d-ad64-093b08af9ffe","type":"Patch"},"selection_glyph":null},"id":"1c471e33-5a32-450a-984a-f34d53e7c8e4","type":"GlyphRenderer"},{"attributes":{"plot":{"id":"d3054252-9bf3-4b4a-946d-ec21f1b67e14","subtype":"Figure","type":"Plot"}},"id":"def23848-bd59-4215-831c-3216c3638295","type":"WheelZoomTool"},{"attributes":{},"id":"a0269c26-dc16-404b-b08e-6ee4d9b830c3","type":"BasicTicker"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"x":{"field":"x"},"y":{"field":"y"}},"id":"42744764-61ac-4a3d-aec5-a1d2114b6ad3","type":"Line"},{"attributes":{"below":[{"id":"8a2fb091-1ae8-4977-b855-b55fd0f83cee","type":"LinearAxis"}],"left":[{"id":"d86293f4-c790-4007-a9a9-fda055ec8e54","type":"LinearAxis"}],"plot_height":278,"plot_width":450,"renderers":[{"id":"8a2fb091-1ae8-4977-b855-b55fd0f83cee","type":"LinearAxis"},{"id":"8eddf051-7cb0-42d5-9563-be9630981c64","type":"Grid"},{"id":"d86293f4-c790-4007-a9a9-fda055ec8e54","type":"LinearAxis"},{"id":"36fb3f0a-5d41-411e-8067-5776dd6fdfd3","type":"Grid"},{"id":"ff6e0877-dee3-4a08-a193-576e63e69290","type":"BoxAnnotation"},{"id":"0a0d9308-899d-4070-8b8a-17eb3d754b46","type":"Legend"},{"id":"cdcf2bae-4e98-4622-b235-f0bed5e8b828","type":"GlyphRenderer"},{"id":"e518ce08-cb0c-456f-be5f-cf207f098c2c","type":"GlyphRenderer"},{"id":"cf3aeb0f-6b81-4716-98f6-190cfb6e1bba","type":"GlyphRenderer"},{"id":"63deb866-5517-4f3f-a443-93bb46088240","type":"GlyphRenderer"}],"title":"Pathway 2","tool_events":{"id":"19149192-c4d5-49c5-ae6a-59bf776330cd","type":"ToolEvents"},"tools":[{"id":"4a41d67a-5ab5-4e81-898f-6e0fdd7377cb","type":"PanTool"},{"id":"2d3fcec0-6418-454e-b883-6579bc7957fe","type":"WheelZoomTool"},{"id":"5a1ebe83-a671-4a57-9295-2eec5a5ab58d","type":"BoxZoomTool"},{"id":"b78302ca-9e8f-41ef-84e2-7940499738c0","type":"PreviewSaveTool"},{"id":"44de6c7e-278e-40ef-9b18-48952f45768b","type":"ResizeTool"},{"id":"dcdd92c2-71de-4e77-9203-63a25d1d8979","type":"ResetTool"},{"id":"bf6d30a3-2b12-498c-9513-0d46073b150e","type":"HelpTool"}],"x_range":{"id":"c1d9b749-8afb-4e84-890b-43d073bebe35","type":"DataRange1d"},"y_range":{"id":"f0e30ec5-0283-41dd-8dda-12cafc521bb2","type":"DataRange1d"}},"id":"25762579-4c02-4e7a-ba6c-eaceea6d3a7a","subtype":"Figure","type":"Plot"},{"attributes":{"callback":null,"column_names":["y","x"],"data":{"x":[0.2878657037040172,0.27271487719327947,0.25756405068254173,0.24241322417180397,0.2272623976610662,0.21211157115032847,0.19696074463959073,0.18180991812885297,0.1666590916181152,0.15150826510737747,0.13635743859663974,0.12120661208590197,0.10605578557516424,0.09090495906442647,0.07575413255368874,0.06060330604295097,0.045452479532213236,0.030301653021475472,0.015150826510737736,0.0,0.0,0.015150826510737736,0.030301653021475472,0.045452479532213236,0.06060330604295097,0.07575413255368874,0.09090495906442647,0.10605578557516424,0.12120661208590197,0.13635743859663974,0.15150826510737747,0.1666590916181152,0.18180991812885297,0.19696074463959073,0.21211157115032847,0.2272623976610662,0.24241322417180397,0.25756405068254173,0.27271487719327947,0.2878657037040172],"y":[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,5.05751391465677,4.941044891092806,4.8204595634977405,4.699874235902684,4.579288908307622,4.458703580712558,4.338118253117486,4.2175329255224225,4.02808326275171,3.674450813452812,3.3208183641539106,2.967185914855023,2.613553465556123,2.2599210162572385,1.8857489538868117,1.5088337665676597,1.1319185792485487,0.7550033919293992,0.37808820461026316,-4.574357385936312e-17]}},"id":"05793d2d-f9ef-4b1c-b768-113ed1651947","type":"ColumnDataSource"},{"attributes":{"plot":{"id":"d3054252-9bf3-4b4a-946d-ec21f1b67e14","subtype":"Figure","type":"Plot"}},"id":"0bd41cd9-3ec1-4b26-96a2-4a385db9185d","type":"ResizeTool"},{"attributes":{"plot":{"id":"25762579-4c02-4e7a-ba6c-eaceea6d3a7a","subtype":"Figure","type":"Plot"}},"id":"44de6c7e-278e-40ef-9b18-48952f45768b","type":"ResizeTool"},{"attributes":{"axis_label":"DM_MNXM754","formatter":{"id":"4d28dbd9-f11b-4adc-98e5-e2664da60062","type":"BasicTickFormatter"},"plot":{"id":"25762579-4c02-4e7a-ba6c-eaceea6d3a7a","subtype":"Figure","type":"Plot"},"ticker":{"id":"7e71eadf-4e0c-49aa-bb38-4401dd36c094","type":"BasicTicker"}},"id":"d86293f4-c790-4007-a9a9-fda055ec8e54","type":"LinearAxis"},{"attributes":{"data_source":{"id":"d9e3490f-2df1-4add-af42-581002819eca","type":"ColumnDataSource"},"glyph":{"id":"b8741388-c46e-42ae-b878-ded478136ba8","type":"Line"},"hover_glyph":null,"nonselection_glyph":{"id":"9a070b4b-eba8-4b93-90db-e0da56c1113d","type":"Line"},"selection_glyph":null},"id":"b15aeabc-8a1f-4fb3-b1ce-8b62f406c792","type":"GlyphRenderer"},{"attributes":{"axis_label":"BIOMASS_SC5_notrace","formatter":{"id":"2bec14f8-77c4-45c4-8768-b610dffcfdbe","type":"BasicTickFormatter"},"plot":{"id":"f9334575-dce2-4399-a5b1-f8add3827aae","subtype":"Figure","type":"Plot"},"ticker":{"id":"b8e042dd-7420-449b-aec5-29b9dadf52d9","type":"BasicTicker"}},"id":"edd685a3-6c8b-42d1-9256-44b4f925f375","type":"LinearAxis"},{"attributes":{"callback":null},"id":"2238db27-74f1-4f8b-98ae-b66e55770e88","type":"DataRange1d"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"x":{"field":"x"},"y":{"field":"y"}},"id":"d735bcf0-423b-45f4-b18f-20bbd609efe8","type":"Line"},{"attributes":{"callback":null,"column_names":["y","x"],"data":{"x":[0.2878657037040178,0.27271487719328,0.25756405068254223,0.24241322417180447,0.2272623976610667,0.21211157115032891,0.19696074463959112,0.18180991812885336,0.1666590916181156,0.1515082651073778,0.13635743859664,0.12120661208590225,0.10605578557516446,0.0909049590644267,0.0757541325536889,0.06060330604295114,0.04545247953221335,0.030301653021475583,0.015150826510737792,0.0,0.0,0.015150826510737792,0.030301653021475583,0.04545247953221335,0.06060330604295114,0.0757541325536889,0.0909049590644267,0.10605578557516446,0.12120661208590225,0.13635743859664,0.1515082651073778,0.1666590916181156,0.18180991812885336,0.19696074463959112,0.21211157115032891,0.2272623976610667,0.24241322417180447,0.25756405068254223,0.27271487719328,0.2878657037040178],"y":[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,3.0354906054279764,2.8757541064179857,2.7160176074079976,2.556281108398011,2.396544609388022,2.2368081103780315,2.077071611368046,1.9173351123580575,1.7575986133480708,1.5978621143380762,1.438125615328088,1.2783891163181074,1.1186526173081195,0.9589161182981297,0.7991796192881254,0.6394431202781523,0.4797066212681532,0.319970122258177,0.16023362324817922,-3.817527732558974e-14]}},"id":"dc1e2d82-40b8-4dbe-aa1f-30e3bffda4b1","type":"ColumnDataSource"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"x":{"field":"x"},"y":{"field":"y"}},"id":"9c181553-e192-4d26-924d-c30b2656ea26","type":"Line"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"x":{"field":"x"},"y":{"field":"y"}},"id":"2bf7ed1d-581f-4c3f-ae41-b10e2e8e86b6","type":"Line"},{"attributes":{"dimension":1,"plot":{"id":"f9334575-dce2-4399-a5b1-f8add3827aae","subtype":"Figure","type":"Plot"},"ticker":{"id":"ce3197d2-2e84-4271-8740-c1a6cdc50cbb","type":"BasicTicker"}},"id":"75d398a1-1343-48f3-9dff-5f0276c19758","type":"Grid"},{"attributes":{"callback":null,"column_names":["y","x"],"data":{"x":[0.2878657037040172,0.27271487719327947,0.25756405068254173,0.24241322417180397,0.2272623976610662,0.21211157115032847,0.19696074463959073,0.18180991812885297,0.1666590916181152,0.15150826510737747,0.13635743859663974,0.12120661208590197,0.10605578557516424,0.09090495906442647,0.07575413255368874,0.06060330604295097,0.045452479532213236,0.030301653021475472,0.015150826510737736,0.0],"y":[-4.574357385936312e-17,0.37808820461026316,0.7550033919293992,1.1319185792485487,1.5088337665676597,1.8857489538868117,2.2599210162572385,2.613553465556123,2.967185914855023,3.3208183641539106,3.674450813452812,4.02808326275171,4.2175329255224225,4.338118253117486,4.458703580712558,4.579288908307622,4.699874235902684,4.8204595634977405,4.941044891092806,5.05751391465677]}},"id":"07d75331-5104-45c3-8ddf-03a470e14612","type":"ColumnDataSource"},{"attributes":{"fill_alpha":{"value":0.1},"fill_color":{"value":"#1f77b4"},"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"x":{"field":"x"},"y":{"field":"y"}},"id":"82eee6dc-d668-4203-b195-e70ac7be2856","type":"Patch"},{"attributes":{},"id":"dfce2671-b2f8-4829-9f3a-1215547e04bd","type":"BasicTickFormatter"},{"attributes":{"plot":{"id":"2e7b3ef2-cd9e-4a57-a834-7d0fd45f429d","subtype":"Figure","type":"Plot"}},"id":"df61542c-f583-47f6-a171-a056f6fd0f17","type":"WheelZoomTool"},{"attributes":{"line_color":{"value":"#B3E2CD"},"x":{"field":"x"},"y":{"field":"y"}},"id":"9171c190-82b5-4fc3-ad35-02d1a183c4fa","type":"Line"},{"attributes":{"plot":{"id":"25762579-4c02-4e7a-ba6c-eaceea6d3a7a","subtype":"Figure","type":"Plot"},"ticker":{"id":"9fc8f2ea-3ce2-42d4-9850-2434e471c2bd","type":"BasicTicker"}},"id":"8eddf051-7cb0-42d5-9563-be9630981c64","type":"Grid"},{"attributes":{"line_color":{"value":"#B3E2CD"},"x":{"field":"x"},"y":{"field":"y"}},"id":"905b4f3f-7acf-4f70-8002-d328d6cb235e","type":"Line"},{"attributes":{"plot":{"id":"25762579-4c02-4e7a-ba6c-eaceea6d3a7a","subtype":"Figure","type":"Plot"}},"id":"dcdd92c2-71de-4e77-9203-63a25d1d8979","type":"ResetTool"},{"attributes":{"plot":{"id":"d3054252-9bf3-4b4a-946d-ec21f1b67e14","subtype":"Figure","type":"Plot"},"ticker":{"id":"a0269c26-dc16-404b-b08e-6ee4d9b830c3","type":"BasicTicker"}},"id":"d0d0a844-36e8-4d20-a34e-2b043f77d4b6","type":"Grid"},{"attributes":{"dimension":1,"plot":{"id":"d3054252-9bf3-4b4a-946d-ec21f1b67e14","subtype":"Figure","type":"Plot"},"ticker":{"id":"46b88b7a-4b01-4a5d-b0c7-ab6a8b177b9d","type":"BasicTicker"}},"id":"4f01b310-2e1f-462b-8f32-9e9f96c613a5","type":"Grid"},{"attributes":{"data_source":{"id":"0b02c023-bca9-4fe4-8874-2c3251a04e1f","type":"ColumnDataSource"},"glyph":{"id":"9ba7116c-79c7-4087-b09f-44ee38523a46","type":"Line"},"hover_glyph":null,"nonselection_glyph":{"id":"9c181553-e192-4d26-924d-c30b2656ea26","type":"Line"},"selection_glyph":null},"id":"b8696d98-e479-4edf-b150-e52ca7f838c3","type":"GlyphRenderer"},{"attributes":{"plot":{"id":"f9334575-dce2-4399-a5b1-f8add3827aae","subtype":"Figure","type":"Plot"}},"id":"815aef33-d194-4873-a26c-3ae0f61fa288","type":"PreviewSaveTool"},{"attributes":{"plot":{"id":"2e7b3ef2-cd9e-4a57-a834-7d0fd45f429d","subtype":"Figure","type":"Plot"},"ticker":{"id":"6a69fad3-7600-4f97-bb5b-375827f6ba5c","type":"BasicTicker"}},"id":"317417ba-dd06-4070-bb81-b50a69e853d8","type":"Grid"},{"attributes":{"callback":null},"id":"c1d9b749-8afb-4e84-890b-43d073bebe35","type":"DataRange1d"},{"attributes":{"fill_alpha":{"value":0.1},"fill_color":{"value":"#1f77b4"},"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"x":{"field":"x"},"y":{"field":"y"}},"id":"59b9c2e1-ee04-4a3d-ad64-093b08af9ffe","type":"Patch"},{"attributes":{"plot":{"id":"2e7b3ef2-cd9e-4a57-a834-7d0fd45f429d","subtype":"Figure","type":"Plot"}},"id":"6d6f68b8-aba7-49da-bf6d-622f77c5d295","type":"ResizeTool"},{"attributes":{"callback":null,"column_names":["y","x"],"data":{"x":[0.2878657037040172,0.27271487719327947,0.25756405068254173,0.24241322417180397,0.2272623976610662,0.21211157115032847,0.19696074463959073,0.18180991812885297,0.1666590916181152,0.15150826510737747,0.13635743859663974,0.12120661208590197,0.10605578557516424,0.09090495906442647,0.07575413255368874,0.06060330604295097,0.045452479532213236,0.030301653021475472,0.015150826510737736,0.0],"y":[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]}},"id":"008bc650-0ddb-4f90-866f-2a1b326cf40b","type":"ColumnDataSource"},{"attributes":{"data_source":{"id":"07d75331-5104-45c3-8ddf-03a470e14612","type":"ColumnDataSource"},"glyph":{"id":"5b386884-5a8b-407a-8c59-32474910c21d","type":"Line"},"hover_glyph":null,"nonselection_glyph":{"id":"e7d46aee-548b-49d6-bcb6-75718d0531f1","type":"Line"},"selection_glyph":null},"id":"18eee7a6-a67d-4825-b407-802e4cb365cf","type":"GlyphRenderer"},{"attributes":{"plot":{"id":"f9334575-dce2-4399-a5b1-f8add3827aae","subtype":"Figure","type":"Plot"}},"id":"786d0d39-5327-4ea4-b4cd-df9bf243c910","type":"HelpTool"},{"attributes":{"fill_alpha":{"value":0.3},"fill_color":{"value":"#B3E2CD"},"line_color":{"value":"#1f77b4"},"x":{"field":"x"},"y":{"field":"y"}},"id":"af44b093-56c9-4ab7-8ddf-a8d628cd3c97","type":"Patch"},{"attributes":{"callback":null,"column_names":["y","x"],"data":{"x":[0.2878657037040163,0.27271487719327864,0.2575640506825409,0.24241322417180322,0.2272623976610655,0.2121115711503278,0.19696074463959012,0.18180991812885242,0.1666590916181147,0.151508265107377,0.1363574385966393,0.12120661208590161,0.1060557855751639,0.0909049590644262,0.07575413255368851,0.060603306042950805,0.0454524795322131,0.030301653021475417,0.01515082651073768,0.0],"y":[2.9290963817141727e-15,0.10587001500621629,0.210946833464358,0.3160236519224981,0.420788967736725,0.5253460487659517,0.6299031297951706,0.7342299103776272,0.8354770791227152,0.9329929217191623,1.0302261739259404,1.1274594261327182,1.2246926783395076,1.3219259305462905,1.4191591827530734,1.5163924349598485,1.6136256871666381,1.710858939373418,1.808092191580202,1.9053254437869822]}},"id":"a3125ab1-dca0-47be-a84a-3a2b02f28ddf","type":"ColumnDataSource"},{"attributes":{},"id":"f3ff280e-8a98-4320-9797-87f66db9d076","type":"BasicTickFormatter"},{"attributes":{"data_source":{"id":"008bc650-0ddb-4f90-866f-2a1b326cf40b","type":"ColumnDataSource"},"glyph":{"id":"9171c190-82b5-4fc3-ad35-02d1a183c4fa","type":"Line"},"hover_glyph":null,"nonselection_glyph":{"id":"d735bcf0-423b-45f4-b18f-20bbd609efe8","type":"Line"},"selection_glyph":null},"id":"6d61e318-a476-4643-8389-cec4948df1b5","type":"GlyphRenderer"},{"attributes":{"overlay":{"id":"e834602c-32e4-4f38-8ccc-9ee6a4f185f9","type":"BoxAnnotation"},"plot":{"id":"d3054252-9bf3-4b4a-946d-ec21f1b67e14","subtype":"Figure","type":"Plot"}},"id":"a323b27b-d527-4e63-a6ca-0d82e3700efe","type":"BoxZoomTool"},{"attributes":{"plot":{"id":"2e7b3ef2-cd9e-4a57-a834-7d0fd45f429d","subtype":"Figure","type":"Plot"}},"id":"2e95323a-3f8b-4751-8246-774696b64549","type":"PanTool"},{"attributes":{},"id":"46b88b7a-4b01-4a5d-b0c7-ab6a8b177b9d","type":"BasicTicker"},{"attributes":{"line_color":{"value":"#B3E2CD"},"x":{"field":"x"},"y":{"field":"y"}},"id":"e10d1364-e882-429d-94b3-642a5c152bb5","type":"Line"},{"attributes":{"legends":[["WT",[{"id":"2166f63b-071f-4c3f-90e1-257720611083","type":"GlyphRenderer"}]]],"plot":{"id":"d3054252-9bf3-4b4a-946d-ec21f1b67e14","subtype":"Figure","type":"Plot"}},"id":"da51939b-a9be-4795-86d8-cff5fd8b3a81","type":"Legend"},{"attributes":{"plot":{"id":"25762579-4c02-4e7a-ba6c-eaceea6d3a7a","subtype":"Figure","type":"Plot"}},"id":"4a41d67a-5ab5-4e81-898f-6e0fdd7377cb","type":"PanTool"},{"attributes":{"data_source":{"id":"05793d2d-f9ef-4b1c-b768-113ed1651947","type":"ColumnDataSource"},"glyph":{"id":"3d63e360-3183-41e1-bed6-676e6f04c7ae","type":"Patch"},"hover_glyph":null,"nonselection_glyph":{"id":"2d7597e4-5b7b-442d-841d-2651845836a3","type":"Patch"},"selection_glyph":null},"id":"5e50fc34-66cf-4e9e-94d3-ed458deafc54","type":"GlyphRenderer"},{"attributes":{},"id":"19149192-c4d5-49c5-ae6a-59bf776330cd","type":"ToolEvents"},{"attributes":{"line_color":{"value":"#B3E2CD"},"x":{"field":"x"},"y":{"field":"y"}},"id":"6a7ae084-5cbc-42c9-a60a-a58f3dcbf73f","type":"Line"},{"attributes":{"plot":{"id":"f9334575-dce2-4399-a5b1-f8add3827aae","subtype":"Figure","type":"Plot"}},"id":"fd474d26-da3c-4ce3-93bc-130e4b3b2380","type":"PanTool"}],"root_ids":["2e7b3ef2-cd9e-4a57-a834-7d0fd45f429d","f9334575-dce2-4399-a5b1-f8add3827aae","25762579-4c02-4e7a-ba6c-eaceea6d3a7a","d3054252-9bf3-4b4a-946d-ec21f1b67e14","04f9bfaa-d2b6-4cb3-ba36-939ef9ab0715"]},"title":"Bokeh Application","version":"0.11.1"}};
                var render_items = [{"docid":"44f7e98d-a65e-4a02-be5b-0131c21d0d70","elementid":"7e6285e3-e28e-48dd-8bb1-84277fe73b86","modelid":"04f9bfaa-d2b6-4cb3-ba36-939ef9ab0715","notebook_comms_target":"2bb10ac4-0689-440b-b198-698d6bfc4809"}];
                
                Bokeh.embed.embed_items(docs_json, render_items);
            });
          },
          function(Bokeh) {
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


