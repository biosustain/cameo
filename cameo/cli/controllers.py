# Copyright 2017 The Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
from cement.core.controller import CementBaseController, expose
from six.moves import input

from cameo.api.designer import design
from cameo.api.hosts import hosts
from cameo.api.products import products
from cameo.parallel import MultiprocessingView, SequentialView

VALID_OUTPUT_FORMATS = ["xlsx", "csv", "tsv"]
OUTPUT_WRITER = {
    "xlsx": lambda df, path: df.to_excel(path),
    "csv": lambda df, path: df.to_csv(path, sep=","),
    "tsv": lambda df, path: df.to_csv(path, sep="\t"),
}
OUTPUT_FORMATS_EXPLANATION = {
    "xlsx": "Excel file",
    "csv": "Comma separated values",
    "tsv": "Tab separated values"
}


class BaseController(CementBaseController):
    class Meta:
        label = 'base'
        description = "cameo high-level API"
        arguments = [
            (['--host'], dict(help='Host to test', action='append')),
            (['-p', '--product'], dict(help='The product to optimize')),
            (['-a', '--anaerobic'], dict(help="Run anaerobic", action='store_true')),
            (['-m', '--multiprocess'], dict(help='Run multiprocess mode', action='store_true')),
            (['-c', '--cores'], dict(help="Number of cores (if multiprocess)")),
            (['-o', '--output'], dict(help="Output file")),
            (['-of', '--output-format'], dict(help="Output file format (default xlsx)\nOptions:%s" %
                                                   VALID_OUTPUT_FORMATS)),
            (['-y', '--yes-all'], dict(help="Auto select suggested options", action="store_true")),
            (['-t', '--test'], dict(help="Test mode", action="store_true"))]

    @expose(hide=True)
    def default(self):
        product = None
        output_format = None
        auto_select = self.app.pargs.yes_all

        if self.app.pargs.product is None:
            print("Argument --product must be provided")
            exit(1)
        else:
            product = self.app.pargs.product

        if self.app.pargs.output_format is None:
            if auto_select:
                output_format = "xlsx"
                print("Output format will be Excel")
            else:
                formats = ["%i - %s (.%s)\n" % (i + 1, OUTPUT_FORMATS_EXPLANATION[ext], ext)
                           for i, ext in enumerate(VALID_OUTPUT_FORMATS)]

                choice = input("Choose your output format:\n" + "".join(formats) + "Enter [1]:")
                if choice == "":
                    choice = "1"
                try:
                    choice = int(choice)
                    output_format = VALID_OUTPUT_FORMATS[choice - 1]
                except TypeError:
                    print("Invalid choice %s" % choice)
                    exit(1)

        else:
            output_format = self.app.pargs.output_format

        if output_format not in VALID_OUTPUT_FORMATS:
            print("Invalid input format")

        if self.app.pargs.output is None:
            output = self.app.pargs.product + "." + output_format
            if auto_select:
                print("Output will be written to %s" % output)
            else:
                choice = input("Enter your output path (%s):" % output)
                if choice != "":
                    output = choice
        else:
            output = self.app.pargs.output

        if self.app.pargs.host is None:
            _hosts = ['ecoli', 'yeast']
        else:
            _hosts = self.app.pargs.host

        _hosts = [getattr(hosts, host) for host in _hosts]

        if self.app.pargs.multiprocess:
            view = MultiprocessingView(processes=self.app.pargs.cores)
        else:
            view = SequentialView()

        design.debug = self.app.pargs.test

        results = design(product=product, hosts=_hosts,
                         view=view, aerobic=not self.app.pargs.anaerobic)

        results['heterologous_pathway'] = results.heterologous_pathway.apply(str)
        results['manipulations'] = results.manipulations.apply(str)

        OUTPUT_WRITER[output_format](results, output)

    @expose(help="Search for products in our internal database")
    def search(self):
        if self.app.pargs.product is None:
            print("Argument --product must be provided")
            exit(1)

        result = products.search(self.app.pargs.product)
        if len(result) > 0:
            print(result[['name', 'formula', 'mass', 'search_rank']])
        else:
            print("No results found")
