#!/usr/bin/env python

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

__all__ = ['main']

import warnings
warnings.simplefilter("ignore", UserWarning)

import multiprocessing
import pickle
import sys

import click
import pandas

pandas.options.display.width = None
pandas.options.display.max_rows = 30

import logging

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel('WARNING')

cameo_logger = logging.getLogger('cameo')
cameo_logger.setLevel('WARNING')

OUTPUT_WRITER = {
    "xlsx": lambda df, path: df.to_excel(path),
    "csv": lambda df, path: df.to_csv(path, sep=","),
    "tsv": lambda df, path: df.to_csv(path, sep="\t"),
}
VALID_OUTPUT_FORMATS = list(OUTPUT_WRITER.keys())
VALID_OUTPUT_FORMATS.append('pickle')


def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    from cameo import __version__
    click.echo(__version__)
    ctx.exit()


@click.group()
@click.option('--version', is_eager=True, callback=print_version, default=False, is_flag=True,
              help='Display version information.')
def main(version):
    """cameo command line interface"""
    pass


@main.command()
@click.option('-o', '--output', default=['-'], multiple=True, type=click.Path(),
              help='Output filename. Multiple output files can be provided (pair with respective format options).')
@click.option('-f', '--format', multiple=True, default=['csv'],
              type=click.Choice(VALID_OUTPUT_FORMATS), help='Output file format (default csv).')
@click.option('-h', '--host', default=['ecoli'], multiple=True,
              type=click.Choice(['ecoli', 'scerevisiae']),
              help='The host organisms to consider (default: all). '
                   'Multiple hosts can be specified by repeating --host HOST.')
@click.option('--aerobic/--anaerobic', default=True, is_flag=True,
              help='Make oxygen available to the host organism (default).')
@click.option('--cores', default=1, type=click.IntRange(1, multiprocessing.cpu_count()),
              help='Number of CPU cores to use (default 1).')
@click.option('--differential-fva/--no-differential-fva', default=True,
              help='Perform differential flux variability analysis to determine flux modulation targets (default).')
@click.option('--heuristic-optimization/--no-heuristic-optimization', default=True,
              help='Find gene knockout targets through heuristic optimization (default).')
@click.option('--max-pathway-predictions', default=1,
              help='Maximum number of heterologous pathways to predict (default 1).')
@click.option('--differential-fva-points', default=10, help='Grid points for differential FVA (default 10).')
@click.option('--pathway-prediction-timeout', default=10 * 60,
              help='Time limit (min) for individual pathway predictions (default 10 min).')
@click.option('--heuristic-optimization-timeout', default=45,
              help='Time limit (min) on individual heuristic optimizations (default 45 min).')
@click.option('--logging', default='ERROR',
              type=click.Choice(['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG']),
              help='Logging level (default WARNING).')
@click.argument('product')
def design(product, host, output, format, cores, aerobic, differential_fva, heuristic_optimization,
           max_pathway_predictions, differential_fva_points, pathway_prediction_timeout,
           heuristic_optimization_timeout, logging):
    """Compute strain designs for desired product."""
    from cameo.api.designer import design
    from cameo.api.hosts import hosts
    from cameo.parallel import MultiprocessingView, SequentialView

    logger.setLevel(logging)
    cameo_logger.setLevel(logging)

    hosts = [getattr(hosts, host) for host in host]

    if cores > 1:
        if sys.platform == 'darwin':  # check if numpy is compatible with multiprocessing
            logger.debug('On OS X, testing if numpy is compatible with multiprocessing')
            import numpy
            logger.debug('numpy.config: {}'.format(numpy.__config__.blas_opt_info))
            if 'Accelerate' in str(numpy.__config__.blas_opt_info):
                click.echo(click.style('WARNING! It looks like the present installation of numpy '
                                       'is linked against the accelerate framework, which '
                                       'unfortunately is incompatible with multiprocessing'
                                       'and might cause the program to hang at some point. '
                                       'Either set --cores 1 to run sequentially or consider '
                                       'installing numpy via conda (https://conda.io/docs/installation.html), '
                                       'which should fix this issue. \n'
                                       'Proceed with caution ...',
                                       fg='red'),
                           err=True)
                click.confirm('Do you want to proceed?', abort=True)

        view = MultiprocessingView(processes=cores)
    elif cores == 1:
        view = SequentialView()

    design.options.pathway_prediction_timeout = pathway_prediction_timeout * 60
    design.options.heuristic_optimization_timeout = heuristic_optimization_timeout
    design.options.max_pathway_predictions = max_pathway_predictions
    design.options.differential_fva = differential_fva
    design.options.heuristic_optimization = heuristic_optimization
    design.options.differential_fva_points = differential_fva_points

    results = design(product=product, hosts=hosts,
                     view=view, aerobic=aerobic)
    click.echo(results)

    for output, format in zip(output, format):
        if format == 'pickle':
            click.echo('Pickling results into {}'.format(output))
            with open(output, 'wb') as f:
                pickle.dump(results, f)
        else:
            click.echo('Writing results into {} using format "{}"'.format(output, format))
            results['heterologous_pathway'] = results.heterologous_pathway.apply(str)
            results['manipulations'] = results.manipulations.apply(str)
            OUTPUT_WRITER[format](results, output)

    click.echo(results)


@main.command()
@click.argument('product')
def search(product):
    """Search for available products.

    PRODUCT: The target product. You can search by name, InChI, and metanetx ID.

    \b
    Examples
    --------
    $ cameo search chebi:30838  # search for itaconate
    """
    click.echo('Loading product database ...')
    from cameo.api.products import products
    click.echo('Searching product database ...')
    result = products.search(product)
    if len(result) > 0:
        click.echo('The following products were found:')
        click.echo(result[['name', 'formula', 'mass', 'InChI', 'search_rank']])
    else:
        click.echo("No products could be found that match the search query.")


if __name__ == '__main__':
    main()
