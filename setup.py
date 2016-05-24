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

versioneer.VCS = 'git'
versioneer.versionfile_source = 'cameo/_version.py'
versioneer.versionfile_build = 'cameo/_version.py'
versioneer.tag_prefix = ''  # tags are like 1.2.0
versioneer.parentdir_prefix = 'myproject-'  # dirname like 'myproject-1.2.0'

requirements = ['numpy>=1.9.1',
                'scipy>=0.14.0',
                'blessings>=1.5.1',
                'pandas>=0.15.2',
                'ordered-set>=1.2',
                'cobra==0.4.0b6',
                'optlang>=0.3.0',
                'requests>=2.5.0',
                'numexpr>=2.4',
                'networkx>=1.9.1',
                'six>=1.9.0',
                'escher>=1.1.2',
                'IProgress>=0.2',
                'inspyred>=1.0',
                'lazy-object-proxy>=1.2.0',
                'palettable>=2.1.1']

extra_requirements = {
    'docs': ['Sphinx>=1.3.5', 'numpydoc>=0.5'],
    'swiglpk': ['swiglpk>=1.2.14'],
    'plotly': ['plotly>=1.9.6'],
    'bokeh': ['bokeh>=0.11.1'],
    'jupyter': ['jupyter>=1.0.0', 'ipywidgets>=4.1.1', 'runipy>=0.1.5'],
    'test': ['nose>=1.3.7', 'rednose>=0.4.3', 'coverage>=4.0.3'],
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
    author='Nikolaus Sonnenschein, Joao Cardoso, Emre Özdemir, Kristian Jensen',
    author_email='niko.sonnenschein@gmail.com',
    description='cameo - computer aided metabolic engineering & optimziation',
    license='Apache License Version 2.0',
    keywords='biology metabolism bioinformatics',
    url='http://cameo.bio',
    long_description=description,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Utilities',
        'Programming Language :: Python :: 2.5',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'License :: OSI Approved :: Apache Software License'
    ],
)
