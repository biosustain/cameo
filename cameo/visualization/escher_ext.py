# Copyright 2015 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import json
import escher

from uuid import uuid4

try:
    from IPython.display import Javascript, display
except ImportError:
    pass


def _builders_var():
    js = Javascript("if(window.escher_builders == undefined) {\n"
                    "\twindow.escher_builders = {};\n"
                    "}")
    display(js)


class NotebookBuilder(escher.Builder):
    def __init__(self, *args, **kwargs):
        kwargs['id'] = "escher_" + str(uuid4()).replace("-", "_")
        super(NotebookBuilder, self).__init__(*args, **kwargs)

    def update(self, reaction_data=None, metabolite_data=None, gene_data=None, reaction_scale=None):
        self.gene_data = gene_data if gene_data is not None else self.gene_data
        self.reaction_data = reaction_data if reaction_data is not None else self.reaction_data
        self.metabolite_data = metabolite_data if metabolite_data is not None else self.metabolite_data

        js = Javascript('var reaction_data_{the_id} = {reaction_data};\n'
                        'var metabolite_data_{the_id} = {metabolite_data};\n'
                        'var gene_data_{the_id} = {gene_data};\n'
                        'var builder = window.escher_builders["{the_id}"];\n'
                        'if(builder) {{\n'
                        '\tbuilder.options.reaction_scale = {reaction_scale}\n'
                        '\tbuilder.set_reaction_data(reaction_data_{the_id});\n'
                        '\tbuilder.set_metabolite_data(metabolite_data_{the_id});\n'
                        '\tbuilder.set_gene_data(gene_data_{the_id});\n'
                        '}}'.format(
            the_id=self.the_id,
            reaction_data=(json.dumps(self.reaction_data) if self.reaction_data else 'null'),
            metabolite_data=(json.dumps(self.metabolite_data) if self.metabolite_data else 'null'),
            gene_data=(json.dumps(self.gene_data) if self.gene_data else 'null'),
            reaction_scale=json.dumps(reaction_scale) if reaction_scale else 'builder.options.reaction_scale'))
        display(js)

    def _draw_js(self, the_id, enable_editing, menu, enable_keys, dev,
                 fill_screen, scroll_behavior, never_ask_before_quit,
                 static_site_index_json):
        _builders_var()
        draw = ("options = {{\n"
                "\tenable_editing: {enable_editing},\n"
                "\tmenu: {menu},\n"
                "\tenable_keys: {enable_keys},\n"
                "\tscroll_behavior: {scroll_behavior},\n"
                "\tfill_screen: {fill_screen},\n"
                "\treaction_data: reaction_data_{the_id},\n"
                "\tmetabolite_data: metabolite_data_{the_id},\n"
                "\tgene_data: gene_data_{the_id},\n"
                "\tnever_ask_before_quit: {never_ask_before_quit}\n").format(
            the_id=the_id,
            enable_editing=json.dumps(enable_editing),
            menu=json.dumps(menu),
            enable_keys=json.dumps(enable_keys),
            scroll_behavior=json.dumps(scroll_behavior),
            fill_screen=json.dumps(fill_screen),
            never_ask_before_quit=json.dumps(never_ask_before_quit))
        # Add the specified options
        for option in self.options:
            val = getattr(self, option)
            if val is None:
                continue
            draw += ",\t{option}: {value}\n".format(option=option, value=json.dumps(val))

        draw += "};\n\n"

        # dev needs escher.
        dev_str = '' if dev else 'escher.'
        # parse the url in javascript, for building a static site
        draw += ('window.escher_builders["{the_id}"] = {dev_str}Builder(\n'
                 '\tmap_data_{the_id},\n'
                 '\tmodel_data_{the_id},\n'
                 '\tembedded_css_{the_id},\n'
                 '\td3.select("#{the_id}"),\n'
                 '\toptions);\n'.format(dev_str=dev_str, the_id=the_id))

        return draw

    def _repr_html_(self):
        return self.display_in_notebook()
