
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

.. code:: python

    from cameo import models
    from cameo.strain_design import pathway_prediction



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
          <td>NAD(+) + vanillin + H2O &lt;=&gt; vanillate + 2.0 H(...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR5340</th>
          <td>formaldehyde + NAD(+) + 3,4-dihydroxybenzoate ...</td>
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
          <td>H2O + 3,4-dihydroxybenzoate + NADP(+) &lt;=&gt; O2 +...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR5336</th>
          <td>NAD(+) + vanillin + H2O &lt;=&gt; vanillate + 2.0 H(...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR5340</th>
          <td>formaldehyde + NAD(+) + 3,4-dihydroxybenzoate ...</td>
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
          <td>NAD(+) + vanillin + H2O &lt;=&gt; vanillate + 2.0 H(...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR5340</th>
          <td>formaldehyde + NAD(+) + 3,4-dihydroxybenzoate ...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR14769</th>
          <td>NAD(+) + H2O + 3,4-dihydroxybenzoate &lt;=&gt; O2 + ...</td>
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
          <td>NAD(+) + vanillin + H2O &lt;=&gt; vanillate + 2.0 H(...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR5340</th>
          <td>formaldehyde + NAD(+) + 3,4-dihydroxybenzoate ...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR5836</th>
          <td>NADPH + O2 + anthranilate + 3.0 H(+) &lt;=&gt; CO(2)...</td>
          <td>-1000</td>
          <td>1000</td>
        </tr>
        <tr>
          <th>MNXR7067</th>
          <td>H(+) + 3,4-dihydroxybenzoate &lt;=&gt; CO(2) + catechol</td>
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
                    <td><strong>Stoichiometry</strong></td><td>MNXM8 + MNXM754 + MNXM2 <=> MNXM982 + 2.0 MNXM1 + MNXM10</td>
                </tr>
                <tr>
                    <td><strong>Lower bound</strong></td><td>-1000.000000</td>
                </tr>
                <tr>
                    <td><strong>Upper bound</strong></td><td>1000.000000</td>
                </tr>
            </table>
            



.. code:: python

    pathways.plot_production_envelopes(model, variables=[model.reactions.biomass_SC5_notrace])



.. raw:: html

    <script type="text/javascript">
            Bokeh.$(function() {
            var all_models = [{"attributes": {"line_color": {"value": "#1f77b4"}, "line_alpha": {"value": 0.1}, "tags": [], "doc": null, "y": {"field": "y"}, "x": {"field": "x"}, "id": "7e8b5c77-2705-4580-a13f-6369cf82cbdd"}, "type": "Line", "id": "7e8b5c77-2705-4580-a13f-6369cf82cbdd"}, {"subtype": "Figure", "type": "Plot", "id": "72310695-e051-442c-915a-4f17c796aad0", "attributes": {"x_range": {"type": "DataRange1d", "id": "549169ff-3a58-47c9-9b1a-2f96a6b720e7"}, "right": [], "tags": [], "tools": [{"type": "PreviewSaveTool", "id": "c8bf4d67-8079-422c-95b4-de888378e7af"}], "title": "Pathway 1", "extra_y_ranges": {}, "plot_width": 450, "renderers": [{"type": "LinearAxis", "id": "9073553f-b946-484f-905a-99ac0c1a5821"}, {"type": "Grid", "id": "5c57fb05-0955-4595-a816-18479fb7a656"}, {"type": "LinearAxis", "id": "0b826ae5-8513-4a04-a2b2-500395960315"}, {"type": "Grid", "id": "c272cb81-262a-4e57-aba8-fc6051c48bd2"}, {"type": "GlyphRenderer", "id": "c5ecf777-6732-4deb-861e-d3818c5b6ff8"}, {"type": "GlyphRenderer", "id": "3f4fec63-0ed5-4135-a6d6-94273ed0a4cb"}, {"type": "GlyphRenderer", "id": "d0af3a9f-7567-4289-81ce-019f646003a5"}, {"type": "GlyphRenderer", "id": "d48e2b8c-7c99-4dac-9553-1267f03ccd3d"}], "extra_x_ranges": {}, "plot_height": 278, "tool_events": {"type": "ToolEvents", "id": "7526f5b8-6928-4fbf-8228-14aaa406cc2b"}, "above": [], "doc": null, "id": "72310695-e051-442c-915a-4f17c796aad0", "y_range": {"type": "DataRange1d", "id": "a7da8d9e-c326-4d65-a697-c481e845cc0d"}, "below": [{"type": "LinearAxis", "id": "9073553f-b946-484f-905a-99ac0c1a5821"}], "left": [{"type": "LinearAxis", "id": "0b826ae5-8513-4a04-a2b2-500395960315"}]}}, {"attributes": {"column_names": ["y", "x"], "tags": [], "doc": null, "selected": {"2d": {"indices": []}, "1d": {"indices": []}, "0d": {"indices": [], "flag": false}}, "callback": null, "data": {"y": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], "x": [0.0, 0.015150826510737745, 0.03030165302147549, 0.045452479532213236, 0.06060330604295098, 0.07575413255368872, 0.09090495906442647, 0.10605578557516421, 0.12120661208590196, 0.1363574385966397, 0.15150826510737744, 0.16665909161811518, 0.18180991812885294, 0.19696074463959068, 0.21211157115032842, 0.22726239766106618, 0.24241322417180391, 0.2575640506825417, 0.2727148771932794, 0.28786570370401715]}, "id": "79da9bd1-a853-4ee1-8d38-e7ea693fb87a"}, "type": "ColumnDataSource", "id": "79da9bd1-a853-4ee1-8d38-e7ea693fb87a"}, {"attributes": {"line_color": {"value": "blue"}, "line_alpha": {"value": 1.0}, "tags": [], "doc": null, "y": {"field": "y"}, "x": {"field": "x"}, "id": "6837bb38-ace6-4acd-a460-fb1a695849b1"}, "type": "Line", "id": "6837bb38-ace6-4acd-a460-fb1a695849b1"}, {"attributes": {"line_color": {"value": "#1f77b4"}, "line_alpha": {"value": 0.1}, "fill_color": {"value": "#1f77b4"}, "tags": [], "doc": null, "fill_alpha": {"value": 0.1}, "y": {"field": "y"}, "x": {"field": "x"}, "id": "009df02f-4fbd-4704-acec-66be116a7ec4"}, "type": "Patch", "id": "009df02f-4fbd-4704-acec-66be116a7ec4"}, {"attributes": {"nonselection_glyph": {"type": "Patch", "id": "009df02f-4fbd-4704-acec-66be116a7ec4"}, "data_source": {"type": "ColumnDataSource", "id": "9cda4cf2-2c98-4e7c-b17e-01f3fc19f73b"}, "tags": [], "doc": null, "selection_glyph": null, "id": "300ab335-337b-46e1-909d-76e7f86f317a", "glyph": {"type": "Patch", "id": "864b1b9e-eb55-4143-aca3-db8f6ff2c279"}}, "type": "GlyphRenderer", "id": "300ab335-337b-46e1-909d-76e7f86f317a"}, {"attributes": {"line_color": {"value": "#1f77b4"}, "line_alpha": {"value": 0.1}, "tags": [], "doc": null, "y": {"field": "y"}, "x": {"field": "x"}, "id": "26f072a9-126f-4bc3-8ffe-f4603a07b524"}, "type": "Line", "id": "26f072a9-126f-4bc3-8ffe-f4603a07b524"}, {"attributes": {"nonselection_glyph": {"type": "Line", "id": "95ab062e-7df7-4873-8d83-ff0569d07d24"}, "data_source": {"type": "ColumnDataSource", "id": "055c53ea-376a-454e-a6bd-00e688406ae6"}, "tags": [], "doc": null, "selection_glyph": null, "id": "3f4fec63-0ed5-4135-a6d6-94273ed0a4cb", "glyph": {"type": "Line", "id": "28fd3577-0f25-4a8e-8677-d70e14368fa5"}}, "type": "GlyphRenderer", "id": "3f4fec63-0ed5-4135-a6d6-94273ed0a4cb"}, {"attributes": {"nonselection_glyph": {"type": "Line", "id": "26f072a9-126f-4bc3-8ffe-f4603a07b524"}, "data_source": {"type": "ColumnDataSource", "id": "bd45635a-0a4d-47d1-8b79-3505e7c8191b"}, "tags": [], "doc": null, "selection_glyph": null, "id": "731f4d07-55b2-41b8-a1e2-e382bba5f33f", "glyph": {"type": "Line", "id": "d2a57cd9-aafb-4c70-81d2-9f74b1d9058c"}}, "type": "GlyphRenderer", "id": "731f4d07-55b2-41b8-a1e2-e382bba5f33f"}, {"attributes": {"nonselection_glyph": {"type": "Patch", "id": "6417a7a6-e78c-4d9d-9e9a-7c89603c1066"}, "data_source": {"type": "ColumnDataSource", "id": "3e58a8cd-c65a-4593-8e2e-00fcee6ba964"}, "tags": [], "doc": null, "selection_glyph": null, "id": "c5ecf777-6732-4deb-861e-d3818c5b6ff8", "glyph": {"type": "Patch", "id": "06652ede-7ffe-455e-b628-5cd34e4d69e6"}}, "type": "GlyphRenderer", "id": "c5ecf777-6732-4deb-861e-d3818c5b6ff8"}, {"attributes": {"column_names": ["y", "x"], "tags": [], "doc": null, "selected": {"2d": {"indices": []}, "1d": {"indices": []}, "0d": {"indices": [], "flag": false}}, "callback": null, "data": {"y": [2.9080000000000004, 2.7549724339484305, 2.601944867896865, 2.4489173018452974, 2.295889735793735, 2.1428621697421546, 1.9898346036906003, 1.8368070376390269, 1.6837794715874503, 1.5307519055358954, 1.377724339484304, 1.2246967734327623, 1.0716692073811938, 0.9186416413296291, 0.765614075278056, 0.6125865092264888, 0.45955894317491236, 0.3065313771233533, 0.1535038110717895, 1.194849404806364e-14], "x": [0.0, 0.015150826510737668, 0.030301653021475337, 0.04545247953221301, 0.06060330604295067, 0.07575413255368835, 0.09090495906442601, 0.10605578557516368, 0.12120661208590135, 0.136357438596639, 0.1515082651073767, 0.16665909161811435, 0.18180991812885203, 0.19696074463958968, 0.21211157115032736, 0.227262397661065, 0.2424132241718027, 0.25756405068254035, 0.272714877193278, 0.2878657037040157]}, "id": "055c53ea-376a-454e-a6bd-00e688406ae6"}, "type": "ColumnDataSource", "id": "055c53ea-376a-454e-a6bd-00e688406ae6"}, {"attributes": {"line_color": {"value": "blue"}, "line_alpha": {"value": 1.0}, "tags": [], "doc": null, "y": {"field": "y"}, "x": {"field": "x"}, "id": "28fd3577-0f25-4a8e-8677-d70e14368fa5"}, "type": "Line", "id": "28fd3577-0f25-4a8e-8677-d70e14368fa5"}, {"attributes": {"line_color": {"value": "#1f77b4"}, "line_alpha": {"value": 0.1}, "tags": [], "doc": null, "y": {"field": "y"}, "x": {"field": "x"}, "id": "95ab062e-7df7-4873-8d83-ff0569d07d24"}, "type": "Line", "id": "95ab062e-7df7-4873-8d83-ff0569d07d24"}, {"attributes": {"tags": [], "doc": null, "renderers": [], "callback": null, "names": [], "id": "549169ff-3a58-47c9-9b1a-2f96a6b720e7"}, "type": "DataRange1d", "id": "549169ff-3a58-47c9-9b1a-2f96a6b720e7"}, {"attributes": {"tags": [], "doc": null, "renderers": [], "callback": null, "names": [], "id": "a7da8d9e-c326-4d65-a697-c481e845cc0d"}, "type": "DataRange1d", "id": "a7da8d9e-c326-4d65-a697-c481e845cc0d"}, {"attributes": {"column_names": ["y", "x"], "tags": [], "doc": null, "selected": {"2d": {"indices": []}, "1d": {"indices": []}, "0d": {"indices": [], "flag": false}}, "callback": null, "data": {"y": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.3545890025884279e-14, 0.3780882046102813, 0.7550033919294131, 1.1319185792485353, 1.5088337665676044, 1.8857489538868188, 2.259921016257238, 2.6135534655561408, 2.9671859148550452, 3.320818364153951, 3.674450813452809, 4.0280832627517125, 4.217532925522428, 4.338118253117486, 4.458703580712547, 4.5792889083076105, 4.699874235902686, 4.820459563497744, 4.9410448910928055, 5.057513914656775], "x": [0.0, 0.0, 0.015150826510737745, 0.03030165302147549, 0.045452479532213236, 0.06060330604295098, 0.07575413255368872, 0.09090495906442647, 0.10605578557516421, 0.12120661208590196, 0.1363574385966397, 0.15150826510737744, 0.16665909161811518, 0.18180991812885294, 0.19696074463959068, 0.21211157115032842, 0.22726239766106618, 0.24241322417180391, 0.2575640506825417, 0.2727148771932794, 0.28786570370401715, 0.28786570370401715, 0.2727148771932794, 0.2575640506825417, 0.24241322417180391, 0.22726239766106618, 0.21211157115032842, 0.19696074463959068, 0.18180991812885294, 0.16665909161811518, 0.15150826510737744, 0.1363574385966397, 0.12120661208590196, 0.10605578557516421, 0.09090495906442647, 0.07575413255368872, 0.06060330604295098, 0.045452479532213236, 0.03030165302147549, 0.015150826510737745, 0.0]}, "id": "9cda4cf2-2c98-4e7c-b17e-01f3fc19f73b"}, "type": "ColumnDataSource", "id": "9cda4cf2-2c98-4e7c-b17e-01f3fc19f73b"}, {"attributes": {"x_range": null, "right": [], "tags": [], "y_range": null, "title": "Production envelops for DM_MNXM754", "extra_y_ranges": {}, "renderers": [], "below": [], "extra_x_ranges": {}, "tool_events": {"type": "ToolEvents", "id": "941cc2ff-b1b7-4b7d-bb19-4448ae0a2eba"}, "above": [], "doc": null, "id": "39b39971-2492-4168-ab62-a0016c553553", "tools": [], "children": [[{"subtype": "Figure", "type": "Plot", "id": "9aefa5de-c6ca-4522-a00a-e8b5b2a8389f"}, {"subtype": "Figure", "type": "Plot", "id": "72310695-e051-442c-915a-4f17c796aad0"}], [{"subtype": "Figure", "type": "Plot", "id": "4027feab-709a-443f-a8f7-abb57c0b22aa"}, {"subtype": "Figure", "type": "Plot", "id": "6a26f3c4-8ada-4486-a678-87e0774840e5"}]], "left": []}, "type": "GridPlot", "id": "39b39971-2492-4168-ab62-a0016c553553"}, {"attributes": {"line_color": {"value": "#99d8c9"}, "line_alpha": {"value": 0.3}, "fill_color": {"value": "#99d8c9"}, "tags": [], "doc": null, "fill_alpha": {"value": 0.3}, "y": {"field": "y"}, "x": {"field": "x"}, "id": "864b1b9e-eb55-4143-aca3-db8f6ff2c279"}, "type": "Patch", "id": "864b1b9e-eb55-4143-aca3-db8f6ff2c279"}, {"attributes": {"nonselection_glyph": {"type": "Line", "id": "41bba77a-1b5e-42de-8af8-c262184a48a1"}, "data_source": {"type": "ColumnDataSource", "id": "6fbe677f-ac77-4fd4-bb4c-4fa3a3bb731d"}, "tags": [], "doc": null, "selection_glyph": null, "id": "8539d51c-8c93-4245-b529-b72bf28be0dd", "glyph": {"type": "Line", "id": "45b29d5d-e76d-4668-af4f-06f298b9f6ec"}}, "type": "GlyphRenderer", "id": "8539d51c-8c93-4245-b529-b72bf28be0dd"}, {"attributes": {"column_names": ["y", "x"], "tags": [], "doc": null, "selected": {"2d": {"indices": []}, "1d": {"indices": []}, "0d": {"indices": [], "flag": false}}, "callback": null, "data": {"y": [-5.339390491992528e-15, 0.0], "x": [0.28786570370401654, 0.28786570370401654]}, "id": "18de7214-b09f-4206-b154-cf7ae9c13ce4"}, "type": "ColumnDataSource", "id": "18de7214-b09f-4206-b154-cf7ae9c13ce4"}, {"attributes": {"doc": null, "id": "08f4d570-3ccd-4596-98cd-70423f8eb9b5", "tags": []}, "type": "BasicTickFormatter", "id": "08f4d570-3ccd-4596-98cd-70423f8eb9b5"}, {"attributes": {"tags": [], "doc": null, "renderers": [], "callback": null, "names": [], "id": "1e083eda-e9a0-4b77-9c51-d1e18674dab0"}, "type": "DataRange1d", "id": "1e083eda-e9a0-4b77-9c51-d1e18674dab0"}, {"attributes": {"plot": {"subtype": "Figure", "type": "Plot", "id": "4027feab-709a-443f-a8f7-abb57c0b22aa"}, "tags": [], "doc": null, "dimension": 0, "ticker": {"type": "BasicTicker", "id": "1fb316db-e2fc-48dc-89b4-76fc181a8d28"}, "id": "5e1bc548-8bd9-4ed0-9295-6ebd70fbda26"}, "type": "Grid", "id": "5e1bc548-8bd9-4ed0-9295-6ebd70fbda26"}, {"attributes": {"line_color": {"value": "#1f77b4"}, "line_alpha": {"value": 0.1}, "tags": [], "doc": null, "y": {"field": "y"}, "x": {"field": "x"}, "id": "41bba77a-1b5e-42de-8af8-c262184a48a1"}, "type": "Line", "id": "41bba77a-1b5e-42de-8af8-c262184a48a1"}, {"attributes": {"tags": [], "doc": null, "renderers": [], "callback": null, "names": [], "id": "61d9418b-12eb-4927-b98d-978c105840d1"}, "type": "DataRange1d", "id": "61d9418b-12eb-4927-b98d-978c105840d1"}, {"attributes": {"line_color": {"value": "blue"}, "line_alpha": {"value": 1.0}, "tags": [], "doc": null, "y": {"field": "y"}, "x": {"field": "x"}, "id": "20076ffd-4d90-4411-a60c-9ef2db84f9a3"}, "type": "Line", "id": "20076ffd-4d90-4411-a60c-9ef2db84f9a3"}, {"attributes": {"line_color": {"value": "#1f77b4"}, "line_alpha": {"value": 0.1}, "tags": [], "doc": null, "y": {"field": "y"}, "x": {"field": "x"}, "id": "ffe2a69a-6c28-4c40-ac28-5b18fdaddd08"}, "type": "Line", "id": "ffe2a69a-6c28-4c40-ac28-5b18fdaddd08"}, {"attributes": {"doc": null, "id": "b06ec8ce-cc24-4d50-b18e-89d91566f062", "tags": []}, "type": "BasicTickFormatter", "id": "b06ec8ce-cc24-4d50-b18e-89d91566f062"}, {"attributes": {"doc": null, "id": "f52999f8-ea6a-43f9-b96a-6104696c4fc6", "tags": []}, "type": "BasicTickFormatter", "id": "f52999f8-ea6a-43f9-b96a-6104696c4fc6"}, {"attributes": {"line_color": {"value": "#99d8c9"}, "line_alpha": {"value": 0.3}, "fill_color": {"value": "#99d8c9"}, "tags": [], "doc": null, "fill_alpha": {"value": 0.3}, "y": {"field": "y"}, "x": {"field": "x"}, "id": "7bfb434f-2938-40d9-bc25-90eec58a5989"}, "type": "Patch", "id": "7bfb434f-2938-40d9-bc25-90eec58a5989"}, {"attributes": {"line_color": {"value": "#1f77b4"}, "line_alpha": {"value": 0.1}, "fill_color": {"value": "#1f77b4"}, "tags": [], "doc": null, "fill_alpha": {"value": 0.1}, "y": {"field": "y"}, "x": {"field": "x"}, "id": "b4855a40-4878-4fc8-b73d-b81bc2db403c"}, "type": "Patch", "id": "b4855a40-4878-4fc8-b73d-b81bc2db403c"}, {"attributes": {"nonselection_glyph": {"type": "Patch", "id": "b4855a40-4878-4fc8-b73d-b81bc2db403c"}, "data_source": {"type": "ColumnDataSource", "id": "06432304-04aa-4c55-9757-258468234482"}, "tags": [], "doc": null, "selection_glyph": null, "id": "6f014756-f630-4279-9dc0-07d2797369e5", "glyph": {"type": "Patch", "id": "7bfb434f-2938-40d9-bc25-90eec58a5989"}}, "type": "GlyphRenderer", "id": "6f014756-f630-4279-9dc0-07d2797369e5"}, {"attributes": {"column_names": ["y", "x"], "tags": [], "doc": null, "selected": {"2d": {"indices": []}, "1d": {"indices": []}, "0d": {"indices": [], "flag": false}}, "callback": null, "data": {"y": [2.0506329113924022, 1.9470191529452257, 1.8399479138329933, 1.7318713798659124, 1.6237854379772216, 1.5155875259015064, 1.4073896138258037, 1.2991917017500922, 1.1909937896743874, 1.0827958775986766, 0.9745979655229724, 0.8664000534472495, 0.7582021413715513, 0.6500042292958477, 0.5418063172201237, 0.4336084051444136, 0.3254104930687175, 0.21721258099299529, 0.10901466891729979, -5.339390491992528e-15], "x": [0.0, 0.015150826510737713, 0.030301653021475427, 0.04545247953221314, 0.060603306042950854, 0.07575413255368857, 0.09090495906442628, 0.106055785575164, 0.12120661208590171, 0.13635743859663943, 0.15150826510737714, 0.16665909161811485, 0.18180991812885255, 0.19696074463959026, 0.212111571150328, 0.2272623976610657, 0.24241322417180342, 0.2575640506825411, 0.27271487719327886, 0.28786570370401654]}, "id": "9ee16a14-b994-43b4-b866-5f2cc88e65e5"}, "type": "ColumnDataSource", "id": "9ee16a14-b994-43b4-b866-5f2cc88e65e5"}, {"attributes": {"doc": null, "id": "0c161790-a97f-405d-baee-7503951b4ec0", "tags": []}, "type": "BasicTickFormatter", "id": "0c161790-a97f-405d-baee-7503951b4ec0"}, {"subtype": "Figure", "type": "Plot", "id": "4027feab-709a-443f-a8f7-abb57c0b22aa", "attributes": {"x_range": {"type": "DataRange1d", "id": "43e0ce86-d462-4000-aba5-3122aa251edf"}, "right": [], "tags": [], "tools": [{"type": "PreviewSaveTool", "id": "f9f80872-06ea-4a46-a9ae-f688a381dcea"}], "title": "Pathway 2", "extra_y_ranges": {}, "plot_width": 450, "renderers": [{"type": "LinearAxis", "id": "279fb185-2e30-4309-b5e8-7f1af4677324"}, {"type": "Grid", "id": "5e1bc548-8bd9-4ed0-9295-6ebd70fbda26"}, {"type": "LinearAxis", "id": "7bfce0db-7531-479f-88da-5cb8c65533bf"}, {"type": "Grid", "id": "fa39704d-9378-4f53-b526-73ab5ae6ecf8"}, {"type": "GlyphRenderer", "id": "6f014756-f630-4279-9dc0-07d2797369e5"}, {"type": "GlyphRenderer", "id": "b885521a-024f-4890-bf9e-eba91eb8ce31"}, {"type": "GlyphRenderer", "id": "8539d51c-8c93-4245-b529-b72bf28be0dd"}, {"type": "GlyphRenderer", "id": "10c2dff1-2966-4614-b86d-d5b828da04be"}], "extra_x_ranges": {}, "plot_height": 278, "tool_events": {"type": "ToolEvents", "id": "72239015-4685-49a3-87e0-ffefb3a59d32"}, "above": [], "doc": null, "id": "4027feab-709a-443f-a8f7-abb57c0b22aa", "y_range": {"type": "DataRange1d", "id": "1e083eda-e9a0-4b77-9c51-d1e18674dab0"}, "below": [{"type": "LinearAxis", "id": "279fb185-2e30-4309-b5e8-7f1af4677324"}], "left": [{"type": "LinearAxis", "id": "7bfce0db-7531-479f-88da-5cb8c65533bf"}]}}, {"attributes": {"geometries": [], "tags": [], "doc": null, "id": "7526f5b8-6928-4fbf-8228-14aaa406cc2b"}, "type": "ToolEvents", "id": "7526f5b8-6928-4fbf-8228-14aaa406cc2b"}, {"attributes": {"plot": {"subtype": "Figure", "type": "Plot", "id": "4027feab-709a-443f-a8f7-abb57c0b22aa"}, "tags": [], "doc": null, "id": "f9f80872-06ea-4a46-a9ae-f688a381dcea"}, "type": "PreviewSaveTool", "id": "f9f80872-06ea-4a46-a9ae-f688a381dcea"}, {"attributes": {"plot": {"subtype": "Figure", "type": "Plot", "id": "6a26f3c4-8ada-4486-a678-87e0774840e5"}, "tags": [], "doc": null, "axis_label": "biomass_SC5_notrace", "formatter": {"type": "BasicTickFormatter", "id": "b06ec8ce-cc24-4d50-b18e-89d91566f062"}, "ticker": {"type": "BasicTicker", "id": "109bf3ef-e8ca-4e66-8165-4d93af32d234"}, "id": "e05ddef5-1167-4427-9f39-eadf1dfd3d45"}, "type": "LinearAxis", "id": "e05ddef5-1167-4427-9f39-eadf1dfd3d45"}, {"attributes": {"nonselection_glyph": {"type": "Line", "id": "fb87b82e-43c3-431c-916c-3044a7d2f3c7"}, "data_source": {"type": "ColumnDataSource", "id": "18de7214-b09f-4206-b154-cf7ae9c13ce4"}, "tags": [], "doc": null, "selection_glyph": null, "id": "10c2dff1-2966-4614-b86d-d5b828da04be", "glyph": {"type": "Line", "id": "a933c8ec-7fb5-407b-97d1-ae6c2f711012"}}, "type": "GlyphRenderer", "id": "10c2dff1-2966-4614-b86d-d5b828da04be"}, {"attributes": {"plot": {"subtype": "Figure", "type": "Plot", "id": "4027feab-709a-443f-a8f7-abb57c0b22aa"}, "tags": [], "doc": null, "dimension": 1, "ticker": {"type": "BasicTicker", "id": "3a89e554-83f2-4981-b30e-96a0f7a2bc36"}, "id": "fa39704d-9378-4f53-b526-73ab5ae6ecf8"}, "type": "Grid", "id": "fa39704d-9378-4f53-b526-73ab5ae6ecf8"}, {"attributes": {"tags": [], "doc": null, "mantissas": [2, 5, 10], "id": "109bf3ef-e8ca-4e66-8165-4d93af32d234", "num_minor_ticks": 5}, "type": "BasicTicker", "id": "109bf3ef-e8ca-4e66-8165-4d93af32d234"}, {"attributes": {"line_color": {"value": "blue"}, "line_alpha": {"value": 1.0}, "tags": [], "doc": null, "y": {"field": "y"}, "x": {"field": "x"}, "id": "9992fa1e-a6d0-4c5f-9b39-57f1b76db7f1"}, "type": "Line", "id": "9992fa1e-a6d0-4c5f-9b39-57f1b76db7f1"}, {"attributes": {"line_color": {"value": "#1f77b4"}, "line_alpha": {"value": 0.1}, "tags": [], "doc": null, "y": {"field": "y"}, "x": {"field": "x"}, "id": "969136cd-7a21-42a7-92d5-022cfcbd4464"}, "type": "Line", "id": "969136cd-7a21-42a7-92d5-022cfcbd4464"}, {"attributes": {"line_color": {"value": "blue"}, "line_alpha": {"value": 1.0}, "tags": [], "doc": null, "y": {"field": "y"}, "x": {"field": "x"}, "id": "d18785f0-e395-40af-9c70-c601e3c63bcd"}, "type": "Line", "id": "d18785f0-e395-40af-9c70-c601e3c63bcd"}, {"attributes": {"plot": {"subtype": "Figure", "type": "Plot", "id": "6a26f3c4-8ada-4486-a678-87e0774840e5"}, "tags": [], "doc": null, "id": "08e4649f-a434-4b09-980c-de742e358409"}, "type": "PreviewSaveTool", "id": "08e4649f-a434-4b09-980c-de742e358409"}, {"attributes": {"line_color": {"value": "#1f77b4"}, "line_alpha": {"value": 0.1}, "tags": [], "doc": null, "y": {"field": "y"}, "x": {"field": "x"}, "id": "f7dc8e8c-5e50-4337-bef0-7391603ac02e"}, "type": "Line", "id": "f7dc8e8c-5e50-4337-bef0-7391603ac02e"}, {"attributes": {"plot": {"subtype": "Figure", "type": "Plot", "id": "4027feab-709a-443f-a8f7-abb57c0b22aa"}, "tags": [], "doc": null, "axis_label": "DM_MNXM754", "formatter": {"type": "BasicTickFormatter", "id": "0c161790-a97f-405d-baee-7503951b4ec0"}, "ticker": {"type": "BasicTicker", "id": "3a89e554-83f2-4981-b30e-96a0f7a2bc36"}, "id": "7bfce0db-7531-479f-88da-5cb8c65533bf"}, "type": "LinearAxis", "id": "7bfce0db-7531-479f-88da-5cb8c65533bf"}, {"attributes": {"nonselection_glyph": {"type": "Line", "id": "f7dc8e8c-5e50-4337-bef0-7391603ac02e"}, "data_source": {"type": "ColumnDataSource", "id": "fec62e3f-cfe0-497c-ae16-4f067940cc4c"}, "tags": [], "doc": null, "selection_glyph": null, "id": "dca95100-936a-4b3c-9229-723affe3dd01", "glyph": {"type": "Line", "id": "d18785f0-e395-40af-9c70-c601e3c63bcd"}}, "type": "GlyphRenderer", "id": "dca95100-936a-4b3c-9229-723affe3dd01"}, {"attributes": {"column_names": ["y", "x"], "tags": [], "doc": null, "selected": {"2d": {"indices": []}, "1d": {"indices": []}, "0d": {"indices": [], "flag": false}}, "callback": null, "data": {"y": [-2.7187337821916038e-14, 0.0], "x": [0.28786570370401754, 0.28786570370401754]}, "id": "be00c777-867b-4424-8fc2-570869322786"}, "type": "ColumnDataSource", "id": "be00c777-867b-4424-8fc2-570869322786"}, {"attributes": {"nonselection_glyph": {"type": "Line", "id": "969136cd-7a21-42a7-92d5-022cfcbd4464"}, "data_source": {"type": "ColumnDataSource", "id": "be00c777-867b-4424-8fc2-570869322786"}, "tags": [], "doc": null, "selection_glyph": null, "id": "f0573ec2-dc8a-48e3-92cc-ae5ec23d60c9", "glyph": {"type": "Line", "id": "9992fa1e-a6d0-4c5f-9b39-57f1b76db7f1"}}, "type": "GlyphRenderer", "id": "f0573ec2-dc8a-48e3-92cc-ae5ec23d60c9"}, {"attributes": {"geometries": [], "tags": [], "doc": null, "id": "941cc2ff-b1b7-4b7d-bb19-4448ae0a2eba"}, "type": "ToolEvents", "id": "941cc2ff-b1b7-4b7d-bb19-4448ae0a2eba"}, {"attributes": {"tags": [], "doc": null, "mantissas": [2, 5, 10], "id": "df25aef3-8a6d-4fad-ada9-ac402a234c46", "num_minor_ticks": 5}, "type": "BasicTicker", "id": "df25aef3-8a6d-4fad-ada9-ac402a234c46"}, {"attributes": {"nonselection_glyph": {"type": "Patch", "id": "b2b461b7-85f6-473f-963d-1fa49346688b"}, "data_source": {"type": "ColumnDataSource", "id": "f5e6f6e1-b5ad-4529-80f4-886ffc20c739"}, "tags": [], "doc": null, "selection_glyph": null, "id": "46da8a7f-d25f-4eeb-8ec3-4d781db012be", "glyph": {"type": "Patch", "id": "a8539196-7923-403b-8f5a-c09ba728fcf1"}}, "type": "GlyphRenderer", "id": "46da8a7f-d25f-4eeb-8ec3-4d781db012be"}, {"attributes": {"column_names": ["y", "x"], "tags": [], "doc": null, "selected": {"2d": {"indices": []}, "1d": {"indices": []}, "0d": {"indices": [], "flag": false}}, "callback": null, "data": {"y": [2.5964285714285706, 2.4597968160253862, 2.3231650606222014, 2.1865333052190072, 2.049901549815823, 1.9132697944126387, 1.7766380390094503, 1.6400062836062712, 1.5033745282030804, 1.366742772799897, 1.2301110173967018, 1.0934792619935207, 0.9568475065903349, 0.8202157511871553, 0.6835839957839618, 0.5469522403807715, 0.41032048497757967, 0.2736887295743981, 0.13705697417122378, -2.7187337821916038e-14], "x": [0.0, 0.015150826510737765, 0.03030165302147553, 0.0454524795322133, 0.06060330604295106, 0.07575413255368883, 0.0909049590644266, 0.10605578557516436, 0.12120661208590212, 0.1363574385966399, 0.15150826510737767, 0.16665909161811543, 0.1818099181288532, 0.19696074463959096, 0.21211157115032872, 0.22726239766106648, 0.24241322417180425, 0.257564050682542, 0.2727148771932798, 0.28786570370401754]}, "id": "e135db10-280f-4aca-92dd-9f7a14cb8643"}, "type": "ColumnDataSource", "id": "e135db10-280f-4aca-92dd-9f7a14cb8643"}, {"attributes": {"line_color": {"value": "#99d8c9"}, "line_alpha": {"value": 0.3}, "fill_color": {"value": "#99d8c9"}, "tags": [], "doc": null, "fill_alpha": {"value": 0.3}, "y": {"field": "y"}, "x": {"field": "x"}, "id": "a8539196-7923-403b-8f5a-c09ba728fcf1"}, "type": "Patch", "id": "a8539196-7923-403b-8f5a-c09ba728fcf1"}, {"attributes": {"line_color": {"value": "#1f77b4"}, "line_alpha": {"value": 0.1}, "fill_color": {"value": "#1f77b4"}, "tags": [], "doc": null, "fill_alpha": {"value": 0.1}, "y": {"field": "y"}, "x": {"field": "x"}, "id": "b2b461b7-85f6-473f-963d-1fa49346688b"}, "type": "Patch", "id": "b2b461b7-85f6-473f-963d-1fa49346688b"}, {"subtype": "Figure", "type": "Plot", "id": "9aefa5de-c6ca-4522-a00a-e8b5b2a8389f", "attributes": {"x_range": {"type": "DataRange1d", "id": "3ebaa444-f82a-4d2d-80f1-4ad6f37ac490"}, "right": [], "tags": [], "tools": [{"type": "PreviewSaveTool", "id": "36c418ac-dac9-4961-ad2f-1afb60cd8572"}], "title": "Pathway 0", "extra_y_ranges": {}, "plot_width": 450, "renderers": [{"type": "LinearAxis", "id": "9e912a7d-1979-4c9a-80ee-6257d53aa655"}, {"type": "Grid", "id": "f1eb2f4d-8d0e-4fe6-b47f-cfc4401dde46"}, {"type": "LinearAxis", "id": "419ec5e1-4f25-4aee-bb5b-2f9b10852444"}, {"type": "Grid", "id": "cc589ba5-d5a2-4024-a5d4-14f725c9c9c2"}, {"type": "GlyphRenderer", "id": "300ab335-337b-46e1-909d-76e7f86f317a"}, {"type": "GlyphRenderer", "id": "731f4d07-55b2-41b8-a1e2-e382bba5f33f"}, {"type": "GlyphRenderer", "id": "8702af68-3b1a-49f2-86bf-cc350b639db4"}, {"type": "GlyphRenderer", "id": "54d6f530-5a0e-4bc2-b093-47dea3baa908"}], "extra_x_ranges": {}, "plot_height": 278, "tool_events": {"type": "ToolEvents", "id": "8a9a5d1a-2511-4a4c-a53f-0df39580732d"}, "above": [], "doc": null, "id": "9aefa5de-c6ca-4522-a00a-e8b5b2a8389f", "y_range": {"type": "DataRange1d", "id": "46f7f3d0-1092-47aa-b007-b74542d63bee"}, "below": [{"type": "LinearAxis", "id": "9e912a7d-1979-4c9a-80ee-6257d53aa655"}], "left": [{"type": "LinearAxis", "id": "419ec5e1-4f25-4aee-bb5b-2f9b10852444"}]}}, {"attributes": {"line_color": {"value": "#1f77b4"}, "line_alpha": {"value": 0.1}, "tags": [], "doc": null, "y": {"field": "y"}, "x": {"field": "x"}, "id": "09a9b505-0415-4c29-810a-5d8e02456993"}, "type": "Line", "id": "09a9b505-0415-4c29-810a-5d8e02456993"}, {"attributes": {"line_color": {"value": "blue"}, "line_alpha": {"value": 1.0}, "tags": [], "doc": null, "y": {"field": "y"}, "x": {"field": "x"}, "id": "d2a57cd9-aafb-4c70-81d2-9f74b1d9058c"}, "type": "Line", "id": "d2a57cd9-aafb-4c70-81d2-9f74b1d9058c"}, {"attributes": {"plot": {"subtype": "Figure", "type": "Plot", "id": "6a26f3c4-8ada-4486-a678-87e0774840e5"}, "tags": [], "doc": null, "axis_label": "DM_MNXM754", "formatter": {"type": "BasicTickFormatter", "id": "08f4d570-3ccd-4596-98cd-70423f8eb9b5"}, "ticker": {"type": "BasicTicker", "id": "d14b880e-3e28-4291-b4ef-bf8e1b562aa6"}, "id": "3617ec1f-7f81-444f-b0de-b67325badf1f"}, "type": "LinearAxis", "id": "3617ec1f-7f81-444f-b0de-b67325badf1f"}, {"attributes": {"column_names": ["y", "x"], "tags": [], "doc": null, "selected": {"2d": {"indices": []}, "1d": {"indices": []}, "0d": {"indices": [], "flag": false}}, "callback": null, "data": {"y": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], "x": [0.0, 0.015150826510737765, 0.03030165302147553, 0.0454524795322133, 0.06060330604295106, 0.07575413255368883, 0.0909049590644266, 0.10605578557516436, 0.12120661208590212, 0.1363574385966399, 0.15150826510737767, 0.16665909161811543, 0.1818099181288532, 0.19696074463959096, 0.21211157115032872, 0.22726239766106648, 0.24241322417180425, 0.257564050682542, 0.2727148771932798, 0.28786570370401754]}, "id": "fec62e3f-cfe0-497c-ae16-4f067940cc4c"}, "type": "ColumnDataSource", "id": "fec62e3f-cfe0-497c-ae16-4f067940cc4c"}, {"attributes": {"plot": {"subtype": "Figure", "type": "Plot", "id": "9aefa5de-c6ca-4522-a00a-e8b5b2a8389f"}, "tags": [], "doc": null, "axis_label": "DM_MNXM754", "formatter": {"type": "BasicTickFormatter", "id": "8b1843c9-023f-4cd7-aaa4-62f036d67eb9"}, "ticker": {"type": "BasicTicker", "id": "df5b2924-fd0f-40aa-86ea-441831cb4bd9"}, "id": "419ec5e1-4f25-4aee-bb5b-2f9b10852444"}, "type": "LinearAxis", "id": "419ec5e1-4f25-4aee-bb5b-2f9b10852444"}, {"attributes": {"tags": [], "doc": null, "mantissas": [2, 5, 10], "id": "df5b2924-fd0f-40aa-86ea-441831cb4bd9", "num_minor_ticks": 5}, "type": "BasicTicker", "id": "df5b2924-fd0f-40aa-86ea-441831cb4bd9"}, {"attributes": {"line_color": {"value": "blue"}, "line_alpha": {"value": 1.0}, "tags": [], "doc": null, "y": {"field": "y"}, "x": {"field": "x"}, "id": "a933c8ec-7fb5-407b-97d1-ae6c2f711012"}, "type": "Line", "id": "a933c8ec-7fb5-407b-97d1-ae6c2f711012"}, {"attributes": {"doc": null, "id": "d0bf4b88-e703-4e93-927e-a7c808263215", "tags": []}, "type": "BasicTickFormatter", "id": "d0bf4b88-e703-4e93-927e-a7c808263215"}, {"attributes": {"plot": {"subtype": "Figure", "type": "Plot", "id": "9aefa5de-c6ca-4522-a00a-e8b5b2a8389f"}, "tags": [], "doc": null, "dimension": 0, "ticker": {"type": "BasicTicker", "id": "df25aef3-8a6d-4fad-ada9-ac402a234c46"}, "id": "f1eb2f4d-8d0e-4fe6-b47f-cfc4401dde46"}, "type": "Grid", "id": "f1eb2f4d-8d0e-4fe6-b47f-cfc4401dde46"}, {"subtype": "Figure", "type": "Plot", "id": "6a26f3c4-8ada-4486-a678-87e0774840e5", "attributes": {"x_range": {"type": "DataRange1d", "id": "61d9418b-12eb-4927-b98d-978c105840d1"}, "right": [], "tags": [], "tools": [{"type": "PreviewSaveTool", "id": "08e4649f-a434-4b09-980c-de742e358409"}], "title": "Pathway 3", "extra_y_ranges": {}, "plot_width": 450, "renderers": [{"type": "LinearAxis", "id": "e05ddef5-1167-4427-9f39-eadf1dfd3d45"}, {"type": "Grid", "id": "2339bc1b-2c44-4978-8a52-e3b135a29cf4"}, {"type": "LinearAxis", "id": "3617ec1f-7f81-444f-b0de-b67325badf1f"}, {"type": "Grid", "id": "de5ccbe4-9171-48b3-9bd1-9d0d3c28b66f"}, {"type": "GlyphRenderer", "id": "46da8a7f-d25f-4eeb-8ec3-4d781db012be"}, {"type": "GlyphRenderer", "id": "e2cb29da-4389-47c6-a293-15a3361cced3"}, {"type": "GlyphRenderer", "id": "dca95100-936a-4b3c-9229-723affe3dd01"}, {"type": "GlyphRenderer", "id": "f0573ec2-dc8a-48e3-92cc-ae5ec23d60c9"}], "extra_x_ranges": {}, "plot_height": 278, "tool_events": {"type": "ToolEvents", "id": "54037db1-5492-4b36-9dd9-c06bd3971ec9"}, "above": [], "doc": null, "id": "6a26f3c4-8ada-4486-a678-87e0774840e5", "y_range": {"type": "DataRange1d", "id": "cfee35b8-7e4b-4de2-83c6-99f62070221f"}, "below": [{"type": "LinearAxis", "id": "e05ddef5-1167-4427-9f39-eadf1dfd3d45"}], "left": [{"type": "LinearAxis", "id": "3617ec1f-7f81-444f-b0de-b67325badf1f"}]}}, {"attributes": {"plot": {"subtype": "Figure", "type": "Plot", "id": "6a26f3c4-8ada-4486-a678-87e0774840e5"}, "tags": [], "doc": null, "dimension": 1, "ticker": {"type": "BasicTicker", "id": "d14b880e-3e28-4291-b4ef-bf8e1b562aa6"}, "id": "de5ccbe4-9171-48b3-9bd1-9d0d3c28b66f"}, "type": "Grid", "id": "de5ccbe4-9171-48b3-9bd1-9d0d3c28b66f"}, {"attributes": {"plot": {"subtype": "Figure", "type": "Plot", "id": "9aefa5de-c6ca-4522-a00a-e8b5b2a8389f"}, "tags": [], "doc": null, "axis_label": "biomass_SC5_notrace", "formatter": {"type": "BasicTickFormatter", "id": "d0bf4b88-e703-4e93-927e-a7c808263215"}, "ticker": {"type": "BasicTicker", "id": "df25aef3-8a6d-4fad-ada9-ac402a234c46"}, "id": "9e912a7d-1979-4c9a-80ee-6257d53aa655"}, "type": "LinearAxis", "id": "9e912a7d-1979-4c9a-80ee-6257d53aa655"}, {"attributes": {"tags": [], "doc": null, "renderers": [], "callback": null, "names": [], "id": "43e0ce86-d462-4000-aba5-3122aa251edf"}, "type": "DataRange1d", "id": "43e0ce86-d462-4000-aba5-3122aa251edf"}, {"attributes": {"line_color": {"value": "#1f77b4"}, "line_alpha": {"value": 0.1}, "tags": [], "doc": null, "y": {"field": "y"}, "x": {"field": "x"}, "id": "fb87b82e-43c3-431c-916c-3044a7d2f3c7"}, "type": "Line", "id": "fb87b82e-43c3-431c-916c-3044a7d2f3c7"}, {"attributes": {"plot": {"subtype": "Figure", "type": "Plot", "id": "9aefa5de-c6ca-4522-a00a-e8b5b2a8389f"}, "tags": [], "doc": null, "id": "36c418ac-dac9-4961-ad2f-1afb60cd8572"}, "type": "PreviewSaveTool", "id": "36c418ac-dac9-4961-ad2f-1afb60cd8572"}, {"attributes": {"doc": null, "id": "8b1843c9-023f-4cd7-aaa4-62f036d67eb9", "tags": []}, "type": "BasicTickFormatter", "id": "8b1843c9-023f-4cd7-aaa4-62f036d67eb9"}, {"attributes": {"plot": {"subtype": "Figure", "type": "Plot", "id": "9aefa5de-c6ca-4522-a00a-e8b5b2a8389f"}, "tags": [], "doc": null, "dimension": 1, "ticker": {"type": "BasicTicker", "id": "df5b2924-fd0f-40aa-86ea-441831cb4bd9"}, "id": "cc589ba5-d5a2-4024-a5d4-14f725c9c9c2"}, "type": "Grid", "id": "cc589ba5-d5a2-4024-a5d4-14f725c9c9c2"}, {"attributes": {"tags": [], "doc": null, "mantissas": [2, 5, 10], "id": "e273d631-4014-46f9-ad2f-c8e504ad2d7f", "num_minor_ticks": 5}, "type": "BasicTicker", "id": "e273d631-4014-46f9-ad2f-c8e504ad2d7f"}, {"attributes": {"nonselection_glyph": {"type": "Line", "id": "ffe2a69a-6c28-4c40-ac28-5b18fdaddd08"}, "data_source": {"type": "ColumnDataSource", "id": "9ee16a14-b994-43b4-b866-5f2cc88e65e5"}, "tags": [], "doc": null, "selection_glyph": null, "id": "b885521a-024f-4890-bf9e-eba91eb8ce31", "glyph": {"type": "Line", "id": "20076ffd-4d90-4411-a60c-9ef2db84f9a3"}}, "type": "GlyphRenderer", "id": "b885521a-024f-4890-bf9e-eba91eb8ce31"}, {"attributes": {"column_names": ["y", "x"], "tags": [], "doc": null, "selected": {"2d": {"indices": []}, "1d": {"indices": []}, "0d": {"indices": [], "flag": false}}, "callback": null, "data": {"y": [5.057513914656775, 4.9410448910928055, 4.820459563497744, 4.699874235902686, 4.5792889083076105, 4.458703580712547, 4.338118253117486, 4.217532925522428, 4.0280832627517125, 3.674450813452809, 3.320818364153951, 2.9671859148550452, 2.6135534655561408, 2.259921016257238, 1.8857489538868188, 1.5088337665676044, 1.1319185792485353, 0.7550033919294131, 0.3780882046102813, -1.3545890025884279e-14], "x": [0.0, 0.015150826510737745, 0.03030165302147549, 0.045452479532213236, 0.06060330604295098, 0.07575413255368872, 0.09090495906442647, 0.10605578557516421, 0.12120661208590196, 0.1363574385966397, 0.15150826510737744, 0.16665909161811518, 0.18180991812885294, 0.19696074463959068, 0.21211157115032842, 0.22726239766106618, 0.24241322417180391, 0.2575640506825417, 0.2727148771932794, 0.28786570370401715]}, "id": "bd45635a-0a4d-47d1-8b79-3505e7c8191b"}, "type": "ColumnDataSource", "id": "bd45635a-0a4d-47d1-8b79-3505e7c8191b"}, {"attributes": {"plot": {"subtype": "Figure", "type": "Plot", "id": "72310695-e051-442c-915a-4f17c796aad0"}, "tags": [], "doc": null, "id": "c8bf4d67-8079-422c-95b4-de888378e7af"}, "type": "PreviewSaveTool", "id": "c8bf4d67-8079-422c-95b4-de888378e7af"}, {"attributes": {"tags": [], "doc": null, "mantissas": [2, 5, 10], "id": "d14b880e-3e28-4291-b4ef-bf8e1b562aa6", "num_minor_ticks": 5}, "type": "BasicTicker", "id": "d14b880e-3e28-4291-b4ef-bf8e1b562aa6"}, {"attributes": {"tags": [], "doc": null, "renderers": [], "callback": null, "names": [], "id": "3ebaa444-f82a-4d2d-80f1-4ad6f37ac490"}, "type": "DataRange1d", "id": "3ebaa444-f82a-4d2d-80f1-4ad6f37ac490"}, {"attributes": {"plot": {"subtype": "Figure", "type": "Plot", "id": "72310695-e051-442c-915a-4f17c796aad0"}, "tags": [], "doc": null, "axis_label": "DM_MNXM754", "formatter": {"type": "BasicTickFormatter", "id": "f52999f8-ea6a-43f9-b96a-6104696c4fc6"}, "ticker": {"type": "BasicTicker", "id": "e273d631-4014-46f9-ad2f-c8e504ad2d7f"}, "id": "0b826ae5-8513-4a04-a2b2-500395960315"}, "type": "LinearAxis", "id": "0b826ae5-8513-4a04-a2b2-500395960315"}, {"attributes": {"geometries": [], "tags": [], "doc": null, "id": "54037db1-5492-4b36-9dd9-c06bd3971ec9"}, "type": "ToolEvents", "id": "54037db1-5492-4b36-9dd9-c06bd3971ec9"}, {"attributes": {"tags": [], "doc": null, "mantissas": [2, 5, 10], "id": "1fb316db-e2fc-48dc-89b4-76fc181a8d28", "num_minor_ticks": 5}, "type": "BasicTicker", "id": "1fb316db-e2fc-48dc-89b4-76fc181a8d28"}, {"attributes": {"column_names": ["y", "x"], "tags": [], "doc": null, "selected": {"2d": {"indices": []}, "1d": {"indices": []}, "0d": {"indices": [], "flag": false}}, "callback": null, "data": {"y": [-1.3545890025884279e-14, 0.0], "x": [0.28786570370401715, 0.28786570370401715]}, "id": "5a4f496b-ff6b-4175-b306-f1edc504c0fc"}, "type": "ColumnDataSource", "id": "5a4f496b-ff6b-4175-b306-f1edc504c0fc"}, {"attributes": {"plot": {"subtype": "Figure", "type": "Plot", "id": "72310695-e051-442c-915a-4f17c796aad0"}, "tags": [], "doc": null, "dimension": 1, "ticker": {"type": "BasicTicker", "id": "e273d631-4014-46f9-ad2f-c8e504ad2d7f"}, "id": "c272cb81-262a-4e57-aba8-fc6051c48bd2"}, "type": "Grid", "id": "c272cb81-262a-4e57-aba8-fc6051c48bd2"}, {"attributes": {"tags": [], "doc": null, "renderers": [], "callback": null, "names": [], "id": "cfee35b8-7e4b-4de2-83c6-99f62070221f"}, "type": "DataRange1d", "id": "cfee35b8-7e4b-4de2-83c6-99f62070221f"}, {"attributes": {"column_names": ["y", "x"], "tags": [], "doc": null, "selected": {"2d": {"indices": []}, "1d": {"indices": []}, "0d": {"indices": [], "flag": false}}, "callback": null, "data": {"y": [1.194849404806364e-14, 0.0], "x": [0.2878657037040157, 0.2878657037040157]}, "id": "5dacbb9c-5790-4ee3-b1c8-18269f9fc2f7"}, "type": "ColumnDataSource", "id": "5dacbb9c-5790-4ee3-b1c8-18269f9fc2f7"}, {"attributes": {"tags": [], "doc": null, "renderers": [], "callback": null, "names": [], "id": "46f7f3d0-1092-47aa-b007-b74542d63bee"}, "type": "DataRange1d", "id": "46f7f3d0-1092-47aa-b007-b74542d63bee"}, {"attributes": {"geometries": [], "tags": [], "doc": null, "id": "72239015-4685-49a3-87e0-ffefb3a59d32"}, "type": "ToolEvents", "id": "72239015-4685-49a3-87e0-ffefb3a59d32"}, {"attributes": {"nonselection_glyph": {"type": "Line", "id": "09a9b505-0415-4c29-810a-5d8e02456993"}, "data_source": {"type": "ColumnDataSource", "id": "e135db10-280f-4aca-92dd-9f7a14cb8643"}, "tags": [], "doc": null, "selection_glyph": null, "id": "e2cb29da-4389-47c6-a293-15a3361cced3", "glyph": {"type": "Line", "id": "76a6fdf1-e63e-46b4-ac16-d220b75d1883"}}, "type": "GlyphRenderer", "id": "e2cb29da-4389-47c6-a293-15a3361cced3"}, {"attributes": {"tags": [], "doc": null, "mantissas": [2, 5, 10], "id": "3a89e554-83f2-4981-b30e-96a0f7a2bc36", "num_minor_ticks": 5}, "type": "BasicTicker", "id": "3a89e554-83f2-4981-b30e-96a0f7a2bc36"}, {"attributes": {"doc": null, "id": "254a8ec2-0a36-4d01-bc4a-f6f23a925e60", "tags": []}, "type": "BasicTickFormatter", "id": "254a8ec2-0a36-4d01-bc4a-f6f23a925e60"}, {"attributes": {"line_color": {"value": "blue"}, "line_alpha": {"value": 1.0}, "tags": [], "doc": null, "y": {"field": "y"}, "x": {"field": "x"}, "id": "76a6fdf1-e63e-46b4-ac16-d220b75d1883"}, "type": "Line", "id": "76a6fdf1-e63e-46b4-ac16-d220b75d1883"}, {"attributes": {"line_color": {"value": "#99d8c9"}, "line_alpha": {"value": 0.3}, "fill_color": {"value": "#99d8c9"}, "tags": [], "doc": null, "fill_alpha": {"value": 0.3}, "y": {"field": "y"}, "x": {"field": "x"}, "id": "06652ede-7ffe-455e-b628-5cd34e4d69e6"}, "type": "Patch", "id": "06652ede-7ffe-455e-b628-5cd34e4d69e6"}, {"attributes": {"plot": {"subtype": "Figure", "type": "Plot", "id": "72310695-e051-442c-915a-4f17c796aad0"}, "tags": [], "doc": null, "axis_label": "biomass_SC5_notrace", "formatter": {"type": "BasicTickFormatter", "id": "c46d3c2e-6fba-43dc-95c5-43cd21df8458"}, "ticker": {"type": "BasicTicker", "id": "7bc1751b-a5be-4d7b-8f3e-64b30ef1ca39"}, "id": "9073553f-b946-484f-905a-99ac0c1a5821"}, "type": "LinearAxis", "id": "9073553f-b946-484f-905a-99ac0c1a5821"}, {"attributes": {"column_names": ["y", "x"], "tags": [], "doc": null, "selected": {"2d": {"indices": []}, "1d": {"indices": []}, "0d": {"indices": [], "flag": false}}, "callback": null, "data": {"y": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], "x": [0.0, 0.015150826510737713, 0.030301653021475427, 0.04545247953221314, 0.060603306042950854, 0.07575413255368857, 0.09090495906442628, 0.106055785575164, 0.12120661208590171, 0.13635743859663943, 0.15150826510737714, 0.16665909161811485, 0.18180991812885255, 0.19696074463959026, 0.212111571150328, 0.2272623976610657, 0.24241322417180342, 0.2575640506825411, 0.27271487719327886, 0.28786570370401654]}, "id": "6fbe677f-ac77-4fd4-bb4c-4fa3a3bb731d"}, "type": "ColumnDataSource", "id": "6fbe677f-ac77-4fd4-bb4c-4fa3a3bb731d"}, {"attributes": {"line_color": {"value": "blue"}, "line_alpha": {"value": 1.0}, "tags": [], "doc": null, "y": {"field": "y"}, "x": {"field": "x"}, "id": "3a7db4c2-6317-44b2-9784-34b3dbea9eec"}, "type": "Line", "id": "3a7db4c2-6317-44b2-9784-34b3dbea9eec"}, {"attributes": {"line_color": {"value": "#1f77b4"}, "line_alpha": {"value": 0.1}, "tags": [], "doc": null, "y": {"field": "y"}, "x": {"field": "x"}, "id": "79044af9-01b2-4a7a-9df1-327228f9139d"}, "type": "Line", "id": "79044af9-01b2-4a7a-9df1-327228f9139d"}, {"attributes": {"geometries": [], "tags": [], "doc": null, "id": "8a9a5d1a-2511-4a4c-a53f-0df39580732d"}, "type": "ToolEvents", "id": "8a9a5d1a-2511-4a4c-a53f-0df39580732d"}, {"attributes": {"nonselection_glyph": {"type": "Line", "id": "79044af9-01b2-4a7a-9df1-327228f9139d"}, "data_source": {"type": "ColumnDataSource", "id": "5dacbb9c-5790-4ee3-b1c8-18269f9fc2f7"}, "tags": [], "doc": null, "selection_glyph": null, "id": "d48e2b8c-7c99-4dac-9553-1267f03ccd3d", "glyph": {"type": "Line", "id": "3a7db4c2-6317-44b2-9784-34b3dbea9eec"}}, "type": "GlyphRenderer", "id": "d48e2b8c-7c99-4dac-9553-1267f03ccd3d"}, {"attributes": {"plot": {"subtype": "Figure", "type": "Plot", "id": "4027feab-709a-443f-a8f7-abb57c0b22aa"}, "tags": [], "doc": null, "axis_label": "biomass_SC5_notrace", "formatter": {"type": "BasicTickFormatter", "id": "254a8ec2-0a36-4d01-bc4a-f6f23a925e60"}, "ticker": {"type": "BasicTicker", "id": "1fb316db-e2fc-48dc-89b4-76fc181a8d28"}, "id": "279fb185-2e30-4309-b5e8-7f1af4677324"}, "type": "LinearAxis", "id": "279fb185-2e30-4309-b5e8-7f1af4677324"}, {"attributes": {"line_color": {"value": "#1f77b4"}, "line_alpha": {"value": 0.1}, "fill_color": {"value": "#1f77b4"}, "tags": [], "doc": null, "fill_alpha": {"value": 0.1}, "y": {"field": "y"}, "x": {"field": "x"}, "id": "6417a7a6-e78c-4d9d-9e9a-7c89603c1066"}, "type": "Patch", "id": "6417a7a6-e78c-4d9d-9e9a-7c89603c1066"}, {"attributes": {"column_names": ["y", "x"], "tags": [], "doc": null, "selected": {"2d": {"indices": []}, "1d": {"indices": []}, "0d": {"indices": [], "flag": false}}, "callback": null, "data": {"y": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -5.339390491992528e-15, 0.10901466891729979, 0.21721258099299529, 0.3254104930687175, 0.4336084051444136, 0.5418063172201237, 0.6500042292958477, 0.7582021413715513, 0.8664000534472495, 0.9745979655229724, 1.0827958775986766, 1.1909937896743874, 1.2991917017500922, 1.4073896138258037, 1.5155875259015064, 1.6237854379772216, 1.7318713798659124, 1.8399479138329933, 1.9470191529452257, 2.0506329113924022], "x": [0.0, 0.0, 0.015150826510737713, 0.030301653021475427, 0.04545247953221314, 0.060603306042950854, 0.07575413255368857, 0.09090495906442628, 0.106055785575164, 0.12120661208590171, 0.13635743859663943, 0.15150826510737714, 0.16665909161811485, 0.18180991812885255, 0.19696074463959026, 0.212111571150328, 0.2272623976610657, 0.24241322417180342, 0.2575640506825411, 0.27271487719327886, 0.28786570370401654, 0.28786570370401654, 0.27271487719327886, 0.2575640506825411, 0.24241322417180342, 0.2272623976610657, 0.212111571150328, 0.19696074463959026, 0.18180991812885255, 0.16665909161811485, 0.15150826510737714, 0.13635743859663943, 0.12120661208590171, 0.106055785575164, 0.09090495906442628, 0.07575413255368857, 0.060603306042950854, 0.04545247953221314, 0.030301653021475427, 0.015150826510737713, 0.0]}, "id": "06432304-04aa-4c55-9757-258468234482"}, "type": "ColumnDataSource", "id": "06432304-04aa-4c55-9757-258468234482"}, {"attributes": {"doc": null, "id": "c46d3c2e-6fba-43dc-95c5-43cd21df8458", "tags": []}, "type": "BasicTickFormatter", "id": "c46d3c2e-6fba-43dc-95c5-43cd21df8458"}, {"attributes": {"line_color": {"value": "#1f77b4"}, "line_alpha": {"value": 0.1}, "tags": [], "doc": null, "y": {"field": "y"}, "x": {"field": "x"}, "id": "4b7efa88-4a92-46a1-983d-618670c5748b"}, "type": "Line", "id": "4b7efa88-4a92-46a1-983d-618670c5748b"}, {"attributes": {"plot": {"subtype": "Figure", "type": "Plot", "id": "6a26f3c4-8ada-4486-a678-87e0774840e5"}, "tags": [], "doc": null, "dimension": 0, "ticker": {"type": "BasicTicker", "id": "109bf3ef-e8ca-4e66-8165-4d93af32d234"}, "id": "2339bc1b-2c44-4978-8a52-e3b135a29cf4"}, "type": "Grid", "id": "2339bc1b-2c44-4978-8a52-e3b135a29cf4"}, {"attributes": {"nonselection_glyph": {"type": "Line", "id": "4b7efa88-4a92-46a1-983d-618670c5748b"}, "data_source": {"type": "ColumnDataSource", "id": "5a4f496b-ff6b-4175-b306-f1edc504c0fc"}, "tags": [], "doc": null, "selection_glyph": null, "id": "54d6f530-5a0e-4bc2-b093-47dea3baa908", "glyph": {"type": "Line", "id": "7d33b2b8-3948-4a0a-9f70-ca8db23f2ee9"}}, "type": "GlyphRenderer", "id": "54d6f530-5a0e-4bc2-b093-47dea3baa908"}, {"attributes": {"nonselection_glyph": {"type": "Line", "id": "7e8b5c77-2705-4580-a13f-6369cf82cbdd"}, "data_source": {"type": "ColumnDataSource", "id": "79da9bd1-a853-4ee1-8d38-e7ea693fb87a"}, "tags": [], "doc": null, "selection_glyph": null, "id": "8702af68-3b1a-49f2-86bf-cc350b639db4", "glyph": {"type": "Line", "id": "6837bb38-ace6-4acd-a460-fb1a695849b1"}}, "type": "GlyphRenderer", "id": "8702af68-3b1a-49f2-86bf-cc350b639db4"}, {"attributes": {"line_color": {"value": "blue"}, "line_alpha": {"value": 1.0}, "tags": [], "doc": null, "y": {"field": "y"}, "x": {"field": "x"}, "id": "45b29d5d-e76d-4668-af4f-06f298b9f6ec"}, "type": "Line", "id": "45b29d5d-e76d-4668-af4f-06f298b9f6ec"}, {"attributes": {"nonselection_glyph": {"type": "Line", "id": "a8b95395-5a8f-42c7-98b9-62c4471964ac"}, "data_source": {"type": "ColumnDataSource", "id": "aa70dace-37c2-4683-81dc-380f19b0cd71"}, "tags": [], "doc": null, "selection_glyph": null, "id": "d0af3a9f-7567-4289-81ce-019f646003a5", "glyph": {"type": "Line", "id": "76e3720f-7808-4778-8d02-8ff9a35505c7"}}, "type": "GlyphRenderer", "id": "d0af3a9f-7567-4289-81ce-019f646003a5"}, {"attributes": {"plot": {"subtype": "Figure", "type": "Plot", "id": "72310695-e051-442c-915a-4f17c796aad0"}, "tags": [], "doc": null, "dimension": 0, "ticker": {"type": "BasicTicker", "id": "7bc1751b-a5be-4d7b-8f3e-64b30ef1ca39"}, "id": "5c57fb05-0955-4595-a816-18479fb7a656"}, "type": "Grid", "id": "5c57fb05-0955-4595-a816-18479fb7a656"}, {"attributes": {"column_names": ["y", "x"], "tags": [], "doc": null, "selected": {"2d": {"indices": []}, "1d": {"indices": []}, "0d": {"indices": [], "flag": false}}, "callback": null, "data": {"y": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.7187337821916038e-14, 0.13705697417122378, 0.2736887295743981, 0.41032048497757967, 0.5469522403807715, 0.6835839957839618, 0.8202157511871553, 0.9568475065903349, 1.0934792619935207, 1.2301110173967018, 1.366742772799897, 1.5033745282030804, 1.6400062836062712, 1.7766380390094503, 1.9132697944126387, 2.049901549815823, 2.1865333052190072, 2.3231650606222014, 2.4597968160253862, 2.5964285714285706], "x": [0.0, 0.0, 0.015150826510737765, 0.03030165302147553, 0.0454524795322133, 0.06060330604295106, 0.07575413255368883, 0.0909049590644266, 0.10605578557516436, 0.12120661208590212, 0.1363574385966399, 0.15150826510737767, 0.16665909161811543, 0.1818099181288532, 0.19696074463959096, 0.21211157115032872, 0.22726239766106648, 0.24241322417180425, 0.257564050682542, 0.2727148771932798, 0.28786570370401754, 0.28786570370401754, 0.2727148771932798, 0.257564050682542, 0.24241322417180425, 0.22726239766106648, 0.21211157115032872, 0.19696074463959096, 0.1818099181288532, 0.16665909161811543, 0.15150826510737767, 0.1363574385966399, 0.12120661208590212, 0.10605578557516436, 0.0909049590644266, 0.07575413255368883, 0.06060330604295106, 0.0454524795322133, 0.03030165302147553, 0.015150826510737765, 0.0]}, "id": "f5e6f6e1-b5ad-4529-80f4-886ffc20c739"}, "type": "ColumnDataSource", "id": "f5e6f6e1-b5ad-4529-80f4-886ffc20c739"}, {"attributes": {"line_color": {"value": "blue"}, "line_alpha": {"value": 1.0}, "tags": [], "doc": null, "y": {"field": "y"}, "x": {"field": "x"}, "id": "7d33b2b8-3948-4a0a-9f70-ca8db23f2ee9"}, "type": "Line", "id": "7d33b2b8-3948-4a0a-9f70-ca8db23f2ee9"}, {"attributes": {"tags": [], "doc": null, "mantissas": [2, 5, 10], "id": "7bc1751b-a5be-4d7b-8f3e-64b30ef1ca39", "num_minor_ticks": 5}, "type": "BasicTicker", "id": "7bc1751b-a5be-4d7b-8f3e-64b30ef1ca39"}, {"attributes": {"line_color": {"value": "blue"}, "line_alpha": {"value": 1.0}, "tags": [], "doc": null, "y": {"field": "y"}, "x": {"field": "x"}, "id": "76e3720f-7808-4778-8d02-8ff9a35505c7"}, "type": "Line", "id": "76e3720f-7808-4778-8d02-8ff9a35505c7"}, {"attributes": {"line_color": {"value": "#1f77b4"}, "line_alpha": {"value": 0.1}, "tags": [], "doc": null, "y": {"field": "y"}, "x": {"field": "x"}, "id": "a8b95395-5a8f-42c7-98b9-62c4471964ac"}, "type": "Line", "id": "a8b95395-5a8f-42c7-98b9-62c4471964ac"}, {"attributes": {"column_names": ["y", "x"], "tags": [], "doc": null, "selected": {"2d": {"indices": []}, "1d": {"indices": []}, "0d": {"indices": [], "flag": false}}, "callback": null, "data": {"y": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], "x": [0.0, 0.015150826510737668, 0.030301653021475337, 0.04545247953221301, 0.06060330604295067, 0.07575413255368835, 0.09090495906442601, 0.10605578557516368, 0.12120661208590135, 0.136357438596639, 0.1515082651073767, 0.16665909161811435, 0.18180991812885203, 0.19696074463958968, 0.21211157115032736, 0.227262397661065, 0.2424132241718027, 0.25756405068254035, 0.272714877193278, 0.2878657037040157]}, "id": "aa70dace-37c2-4683-81dc-380f19b0cd71"}, "type": "ColumnDataSource", "id": "aa70dace-37c2-4683-81dc-380f19b0cd71"}, {"attributes": {"column_names": ["y", "x"], "tags": [], "doc": null, "selected": {"2d": {"indices": []}, "1d": {"indices": []}, "0d": {"indices": [], "flag": false}}, "callback": null, "data": {"y": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.194849404806364e-14, 0.1535038110717895, 0.3065313771233533, 0.45955894317491236, 0.6125865092264888, 0.765614075278056, 0.9186416413296291, 1.0716692073811938, 1.2246967734327623, 1.377724339484304, 1.5307519055358954, 1.6837794715874503, 1.8368070376390269, 1.9898346036906003, 2.1428621697421546, 2.295889735793735, 2.4489173018452974, 2.601944867896865, 2.7549724339484305, 2.9080000000000004], "x": [0.0, 0.0, 0.015150826510737668, 0.030301653021475337, 0.04545247953221301, 0.06060330604295067, 0.07575413255368835, 0.09090495906442601, 0.10605578557516368, 0.12120661208590135, 0.136357438596639, 0.1515082651073767, 0.16665909161811435, 0.18180991812885203, 0.19696074463958968, 0.21211157115032736, 0.227262397661065, 0.2424132241718027, 0.25756405068254035, 0.272714877193278, 0.2878657037040157, 0.2878657037040157, 0.272714877193278, 0.25756405068254035, 0.2424132241718027, 0.227262397661065, 0.21211157115032736, 0.19696074463958968, 0.18180991812885203, 0.16665909161811435, 0.1515082651073767, 0.136357438596639, 0.12120661208590135, 0.10605578557516368, 0.09090495906442601, 0.07575413255368835, 0.06060330604295067, 0.04545247953221301, 0.030301653021475337, 0.015150826510737668, 0.0]}, "id": "3e58a8cd-c65a-4593-8e2e-00fcee6ba964"}, "type": "ColumnDataSource", "id": "3e58a8cd-c65a-4593-8e2e-00fcee6ba964"}];
            Bokeh.load_models(all_models);
            var plots = [{'modeltype': 'GridPlot', 'elementid': '7d6a5139-4e04-4026-bfbb-571aa549a7ce', 'modelid': '39b39971-2492-4168-ab62-a0016c553553'}];
            for (idx in plots) {
            	var plot = plots[idx];
            	var model = Bokeh.Collections(plot.modeltype).get(plot.modelid);
            	Bokeh.logger.info('Realizing plot:')
            	Bokeh.logger.info(' - modeltype: ' + plot.modeltype);
            	Bokeh.logger.info(' - modelid: ' + plot.modelid);
            	Bokeh.logger.info(' - elementid: ' + plot.elementid);
            	var view = new model.default_view({
            		model: model,
            		el: '#' + plot.elementid
            	});
            	Bokeh.index[plot.modelid] = view;
            }
        });
        </script>
    <div class="plotdiv" id="7d6a5139-4e04-4026-bfbb-571aa549a7ce"></div>
    


