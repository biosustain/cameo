# Copyright 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from vbench.api import Benchmark
from datetime import datetime

import os

modules = ['model', 'flux_analysis']

by_module = {}
benchmarks = []

for modname in modules:
    ref = __import__(modname)
    by_module[modname] = [v for v in ref.__dict__.values()
                          if isinstance(v, Benchmark)]
    benchmarks.extend(by_module[modname])

for bm in benchmarks:
    assert(bm.name is not None)

import getpass
import sys

USERNAME = getpass.getuser()

if sys.platform == 'darwin':
    HOME = '/Users/%s' % USERNAME
else:
    HOME = '/home/%s' % USERNAME

try:
    import ConfigParser

    config = ConfigParser.ConfigParser()
    config.readfp(open(os.path.expanduser('~/.vbenchcfg')))

    REPO_PATH = config.get('setup', 'repo_path')
    REPO_URL = config.get('setup', 'repo_url')
    DB_PATH = config.get('setup', 'db_path')
    TMP_DIR = config.get('setup', 'tmp_dir')
except:
    REPO_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "../"))
    REPO_URL = 'git@github.com:biosustain/cameo.git'
    DB_PATH = os.path.join(REPO_PATH, 'vbench/benchmarks.db')
    TMP_DIR = os.path.join('/tmp/vbench_cameo')

PREPARE = """
python setup.py clean
"""
BUILD = """
python setup.py build_ext --inplace
"""
dependencies = ['common.py']

START_DATE = datetime.now()
