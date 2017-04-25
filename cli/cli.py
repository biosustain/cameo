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

import pickle
import multiprocessing

import click
import pandas

pandas.options.display.width = None
pandas.options.display.max_rows = 30

import logging

logging.basicConfig()
logger = logging.getLogger('cameo')
logger.setLevel('WARNING')

inspyred_logger = logging.getLogger('inspyred.ec')
inspyred_logger.setLevel('DEBUG')

OUTPUT_WRITER = {
    "xlsx": lambda df, path: df.to_excel(path),
    "csv": lambda df, path: df.to_csv(path, sep=","),
    "tsv": lambda df, path: df.to_csv(path, sep="\t"),
}
VALID_OUTPUT_FORMATS = list(OUTPUT_WRITER.keys())
VALID_OUTPUT_FORMATS.append('pickle')


@click.group()
def main():
    """Cameo command line interface"""
    pass


@main.command()
@click.option('-o', '--output', default='-', help='Number of greetings.', multiple=True, type=click.Path())
@click.option('-f', '--format', default='pickle', prompt='Format', multiple=True,
              type=click.Choice(VALID_OUTPUT_FORMATS), help='Output file format (default xlsx).')
@click.option('-h', '--host', default='ecoli', multiple=True,
              type=click.Choice(['ecoli', 'scerevisiae']),
              help='The host organisms to consider (default: all). '
                   'Multiple hosts can be specified by repeating --host HOST')
@click.option('--cores', default=1, type=click.IntRange(1, multiprocessing.cpu_count()),
              help='Number of CPU cores to use.')
@click.option('--aerobic', default=True, is_flag=True, help='Anaerobic condicitons.')
@click.option('--differential-fva/--no-differential-fva', default=True, help='...')
@click.option('--heuristic-optimization/--no-heuristic-optimization', default=True, help='...')
@click.option('--max-pathway-predictions', default=1, help='Maximum number of predicted pathways.')
@click.option('--differential-fva-points', default=10, help='Grid pints for differential flux variability analysis.')
@click.option('--pathway-prediction-timeout', default=10 * 60,
              help='Time limit (min) for for individual pathway predictions.')
@click.option('--heuristic-optimization-timeout', default=45,
              help='Time limit (min) on individual heuristic optimizations.')
@click.option('-v', '--verbose', default=True, is_flag=True, help='Extra verbose')
@click.option('--logging', default='ERROR',
              type=click.Choice(['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG']),
              help='Choose between ERROR, INFO, DEBUG (default is WARNING).')
@click.argument('product')
def design(product, host, output, format, cores, aerobic, differential_fva, heuristic_optimization,
           max_pathway_predictions, differential_fva_points, pathway_prediction_timeout,
           heuristic_optimization_timeout, verbose, logging):
    """Compute strain designs for desired products and host organisms."""
    from cameo.api.designer import design
    from cameo.api.hosts import hosts
    from cameo.parallel import MultiprocessingView, SequentialView

    logger.setLevel(logging)

    hosts = [getattr(hosts, host) for host in host]

    if cores > 1:
        # view = MultiprocessingView(processes=cores)
        from multiprocessing import Pool
        view = Pool(processes=cores)
    elif cores == 1:
        view = SequentialView()

    design.options.pathway_prediction_timeout = pathway_prediction_timeout * 60
    click.echo(heuristic_optimization_timeout)
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

    PRODUCT: The target product. You can search by name, InChI, MNX ID, or source and ...

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
