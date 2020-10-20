# -*- coding: utf-8 -*-
# Copyright 2013 Novo Nordisk Foundation Center for Biosustainability,
# Technical University of Denmark.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import absolute_import, print_function

from itertools import chain

from setuptools import setup, find_packages

import versioneer


requirements = ['numpy>=1.9.1',
                'scipy>=0.14.0',
                'blessings>=1.5.1',
                'pandas>=1.1.3',
                'ordered-set>=1.2',
                'cobra>=0.11.1',
                'future>=0.15.2',
                'optlang==1.4.2',
                'numexpr>=2.4',
                'requests>=2.10.0',
                'networkx>=2.4',
                'escher>=1.1.2',
                'IProgress>=0.4',
                'inspyred>=1.0',
                'lazy-object-proxy>=1.2.0',
                'palettable>=2.1.1',
                'gnomic==1.0.1',
                'openpyxl>=2.4.5',
                'click>=6.7']

extra_requirements = {
    'docs': ['Sphinx>=1.3.5', 'numpydoc>=0.5'],
    'plotly': ['plotly>=3.0.0'],
    'bokeh': ['bokeh<=0.12.1'],
    'jupyter': ['jupyter>=1.0.0', 'ipywidgets>=4.1.1'],
    'test': ['pytest', 'pytest-cov', 'pytest-benchmark'],
    'parallel': ['redis>=2.10.5', 'ipyparallel>=5.0.1'],
    'sbml': ['python-libsbml>=5.13.0', 'lxml>=3.6.0']
}
extra_requirements['all'] = list(set(chain(*extra_requirements.values())))

# Run
# pandoc --from=markdown --to=rst README.md -o README.rst
# from time to time, to keep README.rst updated
try:
    with open('README.rst', 'r') as f:
        description = f.read()
except:
    description = ''

setup(
    name='cameo',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    install_requires=requirements,
    extras_require=extra_requirements,
    include_package_data=True,
    entry_points='''
            [console_scripts]
            cameo=cli:main
        ''',
    author='Nikolaus Sonnenschein, Joao Cardoso, Emre Özdemir, Kristian Jensen',
    author_email='niko.sonnenschein@gmail.com',
    description='cameo - computer aided metabolic engineering & optimization',
    license='Apache License Version 2.0',
    keywords='biology metabolism bioinformatics',
    url='http://cameo.bio',
    long_description=description,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Utilities',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'License :: OSI Approved :: Apache Software License'
    ],
)
