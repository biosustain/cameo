# Copyright 2013 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from collections import OrderedDict
from itertools import product, combinations
from uuid import uuid1
from time import time
from datetime import datetime
import colorsys
from pandas.core.common import in_ipnb
import progressbar
import ipython_notebook_utils
import logging
from numpy.random import RandomState

logger = logging.getLogger('cameo')


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
        do()
        current_time = time()
        if bookmark is None:
            entry_id = uuid1()
        else:
            entry_id = bookmark
        # make sure that entry is added to the end of history
        self.history.pop(entry_id, None)
        self.history[entry_id] = {'unix_epoch':
                                  current_time, 'undo': undo, 'redo': do}
        return entry_id

    def __str__(self):
        info = '\n'
        for item in self.history.iteritems():
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
            info += 'undo: ' + undo_entry.func_name + '\n'

        redo_entry = entry['redo']
        try:
            elements = redo_entry.func, redo_entry.args, redo_entry.keywords  # partial
            info += 'redo: ' + ' '.join([str(elem) for elem in elements]) + '\n'
        except AttributeError:
            info += 'redo: ' + redo_entry.func_name + '\n'
        return info

    def undo(self, bookmark=None):
        if bookmark is None:
            try:
                (uuid, entry) = self.history.popitem()
                entry['undo']()
            except KeyError:  # history is empty
                pass
        elif bookmark in self.history.keys():
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
            self.undo(bookmark=self.history.keys()[0])


def partition(lst, n):
    """Partition a list into n bite size chunks."""
    division = len(lst) / float(n)
    return [lst[int(round(division * i)): int(round(division * (i + 1)))] for i in xrange(n)]


def generate_colors(n):
    hsv_tuples = [(v*1.0/n, 0.5, 0.5) for v in xrange(n)]
    color_map = {}
    for i in xrange(n):
        rgb = colorsys.hsv_to_rgb(*hsv_tuples[i])
        color = tuple([rgb[0]*256, rgb[1]*256, rgb[2]*256])
        color_map[i] = '#%02x%02x%02x' % color
    return color_map


class ProgressBar(object):
    def __init__(self, *args, **kwargs):
        if in_ipnb():
            self.progress_bar = ipython_notebook_utils.ProgressBar(*args, **kwargs)
        else:
            self.progress_bar = progressbar.ProgressBar(*args, **kwargs)

    def start(self):
        self.progress_bar.start()

    def update(self, value):
        if isinstance(self.progress_bar, ipython_notebook_utils.ProgressBar):
            self.progress_bar.set(value)
        else:
            self.progress_bar.update(value)

    def __call__(self, iterable):
        if isinstance(self.progress_bar, ipython_notebook_utils.ProgressBar):
            self.progress_bar.size = len(iterable)

            def _(iterable):
                count = 0
                self.progress_bar.set(0)
                for item in iterable:
                    count += 1
                    self.progress_bar.set(count)
                    yield item


        else:
            return self.progress_bar(iterable)


class Timer(object):
    """Taken from http://stackoverflow.com/a/5849861/280182"""

    def __init__(self, name=None):
        self.name = name

    def __enter__(self):
        self.tstart = time()

    def __exit__(self, type, value, traceback):
        if self.name:
            print '[%s]' % self.name,
        print 'Elapsed: %s' % (time() - self.tstart)


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

ID_SANITIZE_RULES_SIMPHENY = [('_DASH_', '__'), ('_FSLASH_', '/'),('_BSLASH_', "\\"), ('_LPAREN_', '('), ('_LSQBKT_', '['),
                     ('_RSQBKT_', ']'), ('_RPAREN_', ')'), ('_COMMA_', ','), ('_PERIOD_', '.'), ('_APOS_', "'"),
                     ('&amp;', '&'), ('&lt;', '<'), ('&gt;', '>'), ('&quot;', '"'), ('__', '-')]

ID_SANITIZE_RULES_TAB_COMPLETION = [('_DASH_', '_dsh_'), ('_FSLASH_', '_fsh_'),('_BSLASH_', "_bsh_"), ('_LPAREN_', '_lp_'), ('_LSQBKT_', '_lb_'),
                     ('_RSQBKT_', '_rb_'), ('_RPAREN_', '_rp_'), ('_COMMA_', '_cm_'), ('_PERIOD_', '_prd_'), ('_APOS_', "_apo_"),
                     ('&amp;', '_amp_'), ('&lt;', '_lt_'), ('&gt;', '_gt_'), ('&quot;', '_qot_'), ('__', '_dsh_')]

def sanitize_ids(model):
    """Makes IDs crippled by the XML specification less annoying.

    For example, 'EX_glc_LPAREN_e_RPAREN_' will be converted to 'EX_glc_lp_e_rp_'.
    Furthermore, reactions and metabolites will be equipped with a `nice_id` attribute
    that provides the original ID, i.e., 'EX_glc(d)'.

    Parameters
    ----------
    model : model

    Notes
    -----

    Will add a nice_id attribute

    """

    def _apply_sanitize_rules(id, rules):
        for rule in rules:
            id = id.replace(*rule)
        return id

    for metabolite in model.metabolites:
        met_id = metabolite.id
        metabolite.id = _apply_sanitize_rules(met_id, ID_SANITIZE_RULES_TAB_COMPLETION)
        metabolite.nice_id = _apply_sanitize_rules(met_id, ID_SANITIZE_RULES_SIMPHENY)
        metabolite.name = _apply_sanitize_rules(metabolite.name, ID_SANITIZE_RULES_SIMPHENY)

    for reaction in model.reactions:
        rxn_id = reaction.id
        reaction.id = _apply_sanitize_rules(rxn_id, ID_SANITIZE_RULES_TAB_COMPLETION)
        reaction.nice_id = _apply_sanitize_rules(rxn_id, ID_SANITIZE_RULES_SIMPHENY)
        reaction.name = _apply_sanitize_rules(reaction.name, ID_SANITIZE_RULES_SIMPHENY)

    for variable in model.solver.variables:
        variable.name = _apply_sanitize_rules(variable.name, ID_SANITIZE_RULES_TAB_COMPLETION)

    for constraint in model.solver.constraints:
        constraint.name = _apply_sanitize_rules(constraint.name, ID_SANITIZE_RULES_TAB_COMPLETION)

    model.repair()
