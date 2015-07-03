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

import os
import sys
from setuptools import setup, find_packages
import versioneer
versioneer.VCS = 'git'
versioneer.versionfile_source = 'cameo/_version.py'
versioneer.versionfile_build = 'cameo/_version.py'
versioneer.tag_prefix = '' # tags are like 1.2.0
versioneer.parentdir_prefix = 'myproject-' # dirname like 'myproject-1.2.0'

on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
if on_rtd:
    requirements = []
else:
    requirements = ['numpy>=1.9.1',
                    'pyzmq>=14.3.1',
                    'ipython>=2.1.0',
                    'scipy>=0.14.0',
                    'blessings>=1.5.1',
                    'Jinja2>=2.7.3',
                    'pandas>=0.15.2',
                    'ordered-set>=1.2',
                    'inspyred>=1.0',
                    'cobra>=0.3.2',
                    'optlang>=0.2.6',
                    'requests>=2.5.0',
                    'numexpr>=2.4',
                    'networkx>=1.9.1',
                    'six>=1.9.0',
                    'escher>=1.1.2'
    ]
    if sys.version_info[0] < 3:
        requirements.extend(['inspyred>=1.0', 'progressbar-ipython>=2.3.1', 'bashplotlib>=0.6.1'])
    else:
        requirements.extend(['progressbar'])

dependency_links = [
    'https://github.com/coagulant/progressbar-python3.git'
]

# Run
# pandoc --from=markdown --to=rst README.md -o README.rst
# from time to time, to keep README.rst updated
with open('README.rst', 'r') as f:
    description = f.read()

setup(
    name='cameo',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    install_requires=requirements,
    dependency_links=dependency_links,
    include_package_data = True,
    author='Nikolaus Sonnenschein, Joao Cardoso, Emre Özdemir',
    author_email='niko.sonnenschein@gmail.com',
    description='cameo - computer aided metabolic engineering & optimziation',
    license='Apache License Version 2.0',
    keywords='biology metabolism bioinformatics',
    url='TBD',
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
