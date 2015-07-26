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

import re
from collections import OrderedDict
from uuid import uuid1
from time import time
from datetime import datetime
import colorsys

import pip
import platform
from itertools import islice
from functools import partial

from numpy.random import RandomState

import logging

logger = logging.getLogger(__name__)


class ProblemCache(object):
    def __init__(self, model):
        self.time_machine = None
        self._model = model
        self.variables = {}
        self.constraints = {}
        self.original_objective = model.objective.expression
        self.time_machine = TimeMachine()
        self.transaction_id = None

    def begin_transaction(self):
        self.transaction_id = uuid1()
        self.time_machine(do=int, undo=int, bookmark=self.transaction_id)

    @property
    def model(self):
        return self._model

    def _rebuild_constraint(self, constraint):
        (expression, lb, ub, name) = constraint.expression, constraint.lb, constraint.ub, constraint.name

        def rebuild():
            self.model.solver._remove_constraint(constraint)
            new_constraint = self.model.solver.interface.Constraint(expression, lb=lb, ub=ub, name=name)
            self.constraints[name] = new_constraint
            self.model.solver._add_constraint(new_constraint, sloppy=True)

        return rebuild

    def _append_constraint(self, constraint_id, create, *args, **kwargs):
        self.constraints[constraint_id] = create(self.model, constraint_id, *args, **kwargs)
        self.model.solver._add_constraint(self.constraints[constraint_id])

    def _remove_constraint(self, constraint_id):
        constraint = self.constraints.pop(constraint_id)
        self.model.solver._remove_constraint(constraint)

    def _append_variable(self, variable_id, create, *args, **kwargs):
        self.variables[variable_id] = create(self.model, variable_id, *args, **kwargs)
        self.model.solver._add_variable(self.variables[variable_id])

    def _remove_variable(self, variable_id):
        variable = self.variables.pop(variable_id)
        self.model.solver._remove_variable(variable)

    def _rebuild_variable(self, variable):
        (type, lb, ub, name) = variable.type, variable.lb, variable.ub, variable.name

        def rebuild():
            self.model.solver._remove_variable(variable)
            new_variable = self.model.solver.interface.Variable(name, lb=lb, ub=ub, type=type)
            self.variables[name] = variable
            self.model.solver._add_variable(new_variable, sloppy=True)

        return rebuild

    def add_constraint(self, constraint_id, create, update, *args, **kwargs):
        if constraint_id in self.constraints:
            if update is not None:
                self.time_machine(
                    do=partial(update, self.model, self.constraints[constraint_id], *args, **kwargs),
                    undo=self._rebuild_constraint(self.constraints[constraint_id]))
        else:
            self.time_machine(
                do=partial(self._append_constraint, constraint_id, create, *args, **kwargs),
                undo=partial(self._remove_constraint, constraint_id)
            )

    def add_variable(self, variable_id, create, update, *args, **kwargs):
        if variable_id in self.variables:
            if update is not None:
                self.time_machine(
                    do=partial(update, self.model, self.variables[variable_id], *args, **kwargs),
                    undo=self._rebuild_variable(self.variables[variable_id]))
        else:
            self.time_machine(
                do=partial(self._append_variable, variable_id, create, *args, **kwargs),
                undo=partial(self._remove_variable, variable_id)
            )

    def reset(self):
        self.model.solver._remove_constraints(self.constraints.values())
        self.model.solver._remove_variables(self.variables.values())
        self.model.objective = self.original_objective
        self.variables = {}
        self.constraints = {}
        self.transaction_id = None
        self.time_machine.history.clear()

    def rollback(self):
        if self.transaction_id is None:
            raise RuntimeError("Start transaction must be called before rollback")
        self.time_machine.undo(self.transaction_id)
        self.transaction_id = None


class RandomGenerator(object):
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
        output = do()
        current_time = time()
        if bookmark is None:
            entry_id = uuid1()
        else:
            entry_id = bookmark
        # make sure that entry is added to the end of history
        self.history.pop(entry_id, None)
        self.history[entry_id] = {'unix_epoch': current_time, 'undo': undo, 'redo': do}
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


# class DisplayItemsWidget(progressbar.widgets.Widget):
#     """Display an items[pbar.currval]
#
#     Examples
#     --------
#     import time
#     from progressbar import Progressbar, widges
#     pbar = ProgressBar(widgets=[DisplayItemsWidget(["asdf"+str(i) for i in range(10)]), widgets.Bar()])
#     pbar.maxval = 10
#     pbar.start()
#     for i in range(10):
#         time.sleep(.2)
#         pbar.update(i)
#     pbar.finish()
#     """
#
#     def __init__(self, items):
#         self.items = items
#
#     def update(self, pbar):
#         try:
#             return "%s" % self.items[pbar.currval]
#         except IndexError:
#             return ""


def partition_(lst, n):
    """Partition a list into n bite size chunks."""
    division = len(lst) / float(n)
    return [lst[int(round(division * i)): int(round(division * (i + 1)))] for i in range(n)]


def partition(ite, n):
    """Partition an iterable into n bite size chunks."""
    try:
        length = len(ite)
    except TypeError:
        ite = list(ite)
        length = len(ite)
    division = length / float(n)
    iterator = iter(ite)
    return [list(islice(iterator, 0, round(division * (i + 1)) - round(division * i))) for i in range(n)]


def flatten(l):
    return [item for sublist in l for item in sublist]


def generate_colors(n):
    hsv_tuples = [(v * 1.0 / n, 0.5, 0.5) for v in range(n)]
    color_map = {}
    for i in range(n):
        rgb = colorsys.hsv_to_rgb(*hsv_tuples[i])
        color = tuple([rgb[0] * 256, rgb[1] * 256, rgb[2] * 256])
        color_map[i] = '#%02x%02x%02x' % color
    return color_map


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


def in_ipnb():
    """
    Check if it is running inside an IPython Notebook (updated for new notebooks)
    """
    try:
        import IPython

        ip = IPython.get_ipython()

        front_end = None
        if "IPKernelApp" in ip.config:
            front_end = ip.config.get('IPKernelApp').get("parent_appname")
        elif "KernelApp" in ip.config:
            front_end = ip.config.get('KernelApp').get("parent_appname")

        if isinstance(front_end, IPython.config.loader.LazyConfigValue) or front_end is None:
            if isinstance(ip, IPython.kernel.zmq.zmqshell.ZMQInteractiveShell):
                return True
            else:
                return False
        elif isinstance(front_end, six.string_types):
            if 'ipython-notebook' in front_end.lower():
                return True
            elif 'notebook' in front_end.lower():
                return True
    except Exception as e:
        logger.debug("Cannot determine if running a notebook because of %s" % e)
        return False
    return False


def str_to_valid_variable_name(s):
    """Adapted from http://stackoverflow.com/a/3303361/280182"""

    # Remove invalid characters
    s = re.sub('[^0-9a-zA-Z_]', '_', s)

    # Remove leading characters until we find a letter or underscore
    s = re.sub('^[^a-zA-Z_]+', '', s)

    return s
