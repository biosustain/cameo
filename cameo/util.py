# Copyright 2013 Novo Nordisk Foundation Center for Biosustainability, DTU.
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

import six
from six.moves import range

from collections import OrderedDict
from uuid import uuid1
from time import time
from datetime import datetime
import colorsys
import pip
import platform

import progressbar
from numpy.random import RandomState

import logging
logger = logging.getLogger(__name__)


class RandomGenerator():
    def __init__(self, seed=None):
        self._random = RandomState(seed=seed)

    def random(self):
        return self._random.rand()

    def randint(self, a, b=None):
        if b is None:
            b = a
            a = 0
        return self._random.randint(a, b)

    def gauss(self, mean, stdev):
        return self._random.normal(mean, stdev)

    def sample(self, population, k):
        if k == 0:
            return []
        return self._random.choice(population, size=k, replace=True)

    def __getattr__(self, attr):
        return getattr(self._random, attr)

    def __getstate__(self):
        return {'_random': self._random}

    def __setstate__(self, d):
        self._random = d['_random']


class Singleton(object):
    """
    Singleton class to be extended
    """
    _instance = None

    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            cls._instance = super(Singleton, cls).__new__(cls, *args, **kwargs)
        return cls._instance


class AutoVivification(dict):

    """Implementation of perl's autovivification feature. Checkout http://stackoverflow.com/a/652284/280182"""

    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


class TimeMachine(object):

    """Travel back and forth in time."""

    def __init__(self):
        super(TimeMachine, self).__init__()
        self.history = OrderedDict()

    def __call__(self, do=None, undo=None, bookmark=None):
        output = do()
        current_time = time()
        if bookmark is None:
            entry_id = uuid1()
        else:
            entry_id = bookmark
        # make sure that entry is added to the end of history
        self.history.pop(entry_id, None)
        self.history[entry_id] = {'unix_epoch':
                                  current_time, 'undo': undo, 'redo': do}
        return entry_id, output

    def __str__(self):
        info = '\n'
        for item in six.iteritems(self.history):
            info += self._history_item_to_str(item)
        return info

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.reset()

    @staticmethod
    def _history_item_to_str(item):
        info = ''
        uuid, entry = item
        info += datetime.fromtimestamp(entry['unix_epoch']).strftime('%Y-%m-%d %H:%M:%S') + '\n'
        undo_entry = entry['undo']
        try:
            elements = undo_entry.func, undo_entry.args, undo_entry.keywords  # partial
            info += 'undo: ' + ' '.join([str(elem) for elem in elements]) + '\n'
        except AttributeError:  # normal python function
            info += 'undo: ' + undo_entry.__name__ + '\n'

        redo_entry = entry['redo']
        try:
            elements = redo_entry.func, redo_entry.args, redo_entry.keywords  # partial
            info += 'redo: ' + ' '.join([str(elem) for elem in elements]) + '\n'
        except AttributeError:
            info += 'redo: ' + redo_entry.__name__ + '\n'
        return info

    def undo(self, bookmark=None):
        if bookmark is None:
            try:
                (uuid, entry) = self.history.popitem()
                entry['undo']()
            except KeyError:  # history is empty
                pass
        elif bookmark in list(self.history.keys()):
            uuid = False
            while uuid is not bookmark:
                (uuid, entry) = self.history.popitem()
                entry['undo']()
        else:
            raise Exception(
                'Provided bookmark %s cannot be found in the time machine.')

    def redo(self):
        raise NotImplementedError

    def reset(self):
        if self.history:  # history is not empty
            self.undo(bookmark=list(self.history.keys())[0])


def partition(lst, n):
    """Partition a list into n bite size chunks."""
    division = len(lst) / float(n)
    return [lst[int(round(division * i)): int(round(division * (i + 1)))] for i in range(n)]


def generate_colors(n):
    hsv_tuples = [(v*1.0/n, 0.5, 0.5) for v in range(n)]
    color_map = {}
    for i in range(n):
        rgb = colorsys.hsv_to_rgb(*hsv_tuples[i])
        color = tuple([rgb[0]*256, rgb[1]*256, rgb[2]*256])
        color_map[i] = '#%02x%02x%02x' % color
    return color_map


class Timer(object):
    """Taken from http://stackoverflow.com/a/5849861/280182"""

    def __init__(self, name=None):
        self.name = name

    def __enter__(self):
        self.tstart = time()

    def __exit__(self, type, value, traceback):
        if self.name:
            print('[%s]' % self.name, end=' ')
        print('Elapsed: %s' % (time() - self.tstart))


def memoize(function, memo={}):
    def wrapper(*args):
        logger.debug("key = %s" % str(args))
        if args in memo:
            logger.debug("Key found")
            return memo[args]
        else:
            logger.debug("Key not found")
            rv = function(*args)
            memo[args] = rv
            return rv
    return wrapper


class IntelliContainer(object):
    def __init__(self, **kwargs):
        self._dict = dict(**kwargs)

    def __getattr__(self, value):
        return self._dict.get(value)

    def __setitem__(self, key, value):
        self._dict[key] = value

    def __iter__(self):
        return six.itervalues(self._dict)

    def __dir__(self):
        return list(self._dict.keys())


class DisplayItemsWidget(progressbar.widgets.Widget):
    """Display an items[pbar.currval]

    Examples
    --------
    import time
    from progressbar import Progressbar, widges
    pbar = ProgressBar(widgets=[DisplayItemsWidget(["asdf"+str(i) for i in range(10)]), widgets.Bar()])
    pbar.maxval = 10
    pbar.start()
    for i in range(10):
        time.sleep(.2)
        pbar.update(i)
    pbar.finish()
    """

    def __init__(self, items):
        self.items = items

    def update(self, pbar):
        try:
            return "%s" % self.items[pbar.currval]
        except IndexError:
            return ""

def get_system_info():
    # pip freeze (adapted from http://stackoverflow.com/a/24322465/280182)
    package_info = list()
    for dist in pip.get_installed_distributions():
        req = str(dist.as_requirement())
        package_info.append(req)
    return dict(package_info=package_info,
                platform=platform.platform(),
                machine=platform.machine(),
                system=platform.system())
