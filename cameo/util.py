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

import colorsys
import inspect
import itertools
import logging
import math
import platform
import re
from collections import OrderedDict
from datetime import datetime
from functools import partial
from itertools import islice
from time import time
from uuid import uuid1

import numpy
import pandas
import pkg_resources
from cobra.util.context import HistoryManager
from numpy.random import RandomState

logger = logging.getLogger(__name__)

_BIOMASS_RE_ = re.compile("biomass", re.IGNORECASE)


class frozendict(dict):
    def __init__(self, iterable, **kwargs):
        super(frozendict, self).__init__(iterable, **kwargs)

    def popitem(self):
        raise AttributeError("'frozendict' object has no attribute 'popitem")

    def pop(self, k, d=None):
        raise AttributeError("'frozendict' object has no attribute 'pop")

    def __setitem__(self, key, value):
        raise AttributeError("'frozendict' object has no attribute '__setitem__")

    def setdefault(self, k, d=None):
        raise AttributeError("'frozendict' object has no attribute 'setdefault")

    def __delitem__(self, key):
        raise AttributeError("'frozendict' object has no attribute '__delitem__")

    def __hash__(self):
        return hash(tuple(sorted(self.items())))

    def update(self, E=None, **F):
        raise AttributeError("'frozendict' object has no attribute 'update")


def float_ceil(val, decimals=0):
    """
    Like ceil but the number of decimals can be set.

    Equivalent of $$round(val + 1e^{-decimals}/2, decimals)$$

    val: float, numpy.array
        The initial value.
    decimals: int
        The number of decimal places.

    Returns
    -------
    float, numpy.array

    """
    aux = math.pow(10, -decimals) / 2

    return numpy.round(val + aux, decimals)


def float_floor(val, decimals=0):
    """
    Like floor but the number of decimals can be set.

    Equivalent of $$round(val - 1e^{-decimals}/2, decimals)$$

    val: float, numpy.array
        The initial value.
    decimals: int
        The number of decimal places.

    Returns
    -------
    float, numpy.array

    """
    aux = math.pow(10, -decimals) / 2

    return numpy.round(val - aux, decimals)


class ProblemCache(object):
    """
    Variable and constraint cache for models.

    To be used in complex methods that require many extra variables and constraints when one must run
    simulations with the same method many times.

    It allows rollback to the previous state in case one iteration fails to build the problem or
    generates an invalid state.
    """

    def __init__(self, model):
        self.history_manager = None
        self._model = model
        self.variables = {}
        self.constraints = {}
        self.objective = None
        self.original_objective = model.solver.objective
        self._contexts = [HistoryManager()]
        self.transaction_id = None

    def begin_transaction(self):
        """
        Creates a time point. If rollback is called, the variables and constrains will be reverted to this point.
        """
        self._contexts.append(HistoryManager())

    @property
    def model(self):
        return self._model

    def _append_constraint(self, constraint_id, create, *args, **kwargs):
        constraint = self.constraints[constraint_id] = create(self._model, constraint_id, *args, **kwargs)
        assert constraint_id in self.constraints
        self._model.solver.add(constraint)

    def _remove_constraint(self, constraint_id):
        constraint = self.constraints.pop(constraint_id)
        self._model.solver.remove(constraint)

    def _append_variable(self, variable_id, create, *args, **kwargs):
        variable = self.variables[variable_id] = create(self._model, variable_id, *args, **kwargs)
        assert variable_id in self.variables
        self._model.solver.add(variable)

    def _remove_variable(self, variable_id):
        variable = self.variables.pop(variable_id)
        self._model.solver.remove(variable)

    def _rebuild_variable(self, variable):
        (type, lb, ub, name) = variable.type, variable.lb, variable.ub, variable.name

        def rebuild():
            self._model.solver.remove(variable)
            new_variable = self.model.solver.interface.Variable(name, lb=lb, ub=ub, type=type)
            self.variables[name] = variable
            self._model.solver.add(new_variable, sloppy=True)

        return rebuild

    def add_constraint(self, constraint_id, create, update, *args, **kwargs):
        """
        Adds a new cached constraint.

        The create and update functions must have the following signatures:
        >>> create(model, constraint_id, *args)
        >>> update(model, constraint, *args)

        "args" in the first example must match args on the second example.

        Parameters
        ----------
        constraint_id : str
            The identifier of the constraint
        create : function
            A function that creates an optlang.interface.Constraint
        update : function
            a function that updates an optlang.interface.Constraint

        """
        context = self._contexts[-1]

        if constraint_id not in self.constraints:
            self._append_constraint(constraint_id, create, *args, **kwargs)
            context(partial(self._remove_constraint, constraint_id))
        elif update is not None:
            update(self._model, self.constraints[constraint_id], *args, **kwargs)

        assert constraint_id in self.constraints

    def add_variable(self, variable_id, create, update, *args, **kwargs):
        """
        Adds a new cached variable.

        The create and update functions must have the following signatures:
        >>> create(model, variable_id, *args)
        >>> update(model, variable, *args)

        "args" in the first example must match args on the second example.

        Parameters
        ----------
        variable_id : str
            The identifier of the constraint
        create : function
            A function that creates an optlang.interface.Variable
        update : function
            a function that updates an optlang.interface.Variable

        """
        context = self._contexts[-1]
        if variable_id not in self.variables:
            self._append_variable(variable_id, create, *args, **kwargs)
            context(partial(self._remove_variable, variable_id))
        elif update is not None:
            # rebuild_function = self._rebuild_variable(self.variables[variable_id])
            update(self._model, self.variables[variable_id], *args, **kwargs)
            # context(rebuild_function)

        assert variable_id in self.variables

    def add_objective(self, create, update, *args):
        context = self._contexts[-1]
        if self.objective is None:
            previous_objective = self._model.solver.objective
            self.model.solver.objective = self.objective = create(self._model, *args)
            context(partial(setattr, self._model.solver, 'objective', previous_objective))

        elif update:
            previous_objective = self._model.solver.objective
            self.model.solver.objective = self.objective = update(self._model, *args)
            context(partial(setattr, self._model.solver, 'objective', previous_objective))

    def reset(self):
        """
        Removes all constraints and variables from the cache.
        """
        variables = list(self.variables.keys())
        constraints = list(self.constraints.keys())
        while len(self._contexts) > 0:
            manager = self._contexts.pop()
            manager.reset()
        self._contexts.append(HistoryManager())
        assert all(var_id not in self._model.solver.variables for var_id in variables)
        assert all(const_id not in self._model.solver.constraints for const_id in constraints)
        self.variables = {}
        self.constraints = {}
        self._model.objective = self.original_objective
        self.objective = None

    def rollback(self):
        """
        Returns to the previous transaction start point.
        """
        if len(self._contexts) < 2:
            raise RuntimeError("Start transaction must be called before rollback")
        self._contexts.pop().reset()

    def __enter__(self):
        """
        Allows problem cache to be used with a _with_ statement.

        Examples
        --------
        You want to run room/lmoma for every single knockout.
        >>> with ProblemCache(model) as cache:
        >>>     for reaction in reactions:
        >>>         reaction.knock_out()
        >>>         result = lmoma(model, reference=reference, cache=cache)

        Returns
        -------
        ProblemCache
            returns itself
        """
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.reset()


class RandomGenerator(object):
    def __init__(self, seed=None):
        self._random = RandomState(seed=seed)

    def seed(self, seed):
        self._random = RandomState(seed=seed)

    def random(self):
        return self._random.rand()

    def randint(self, a, b=None):
        if b is None:
            b = a
            a = 0
        r = self._random.randint(a, high=b, size=1)
        return r[0]

    def sample(self, population, k):
        if k == 0:
            return []
        return list(self._random.choice(population, size=k, replace=False))

    def __getattr__(self, attr):
        return getattr(self._random, attr)

    def __getstate__(self):
        return {'_random': self._random}

    def __setstate__(self, d):
        self._random = d['_random']

    def uniform(self, low=0.0, high=1.0, size=None):
        return self._random.uniform(low, high, size)


class Singleton(object):
    """
    Singleton class to be extended
    """
    _instance = None

    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            cls._instance = super(Singleton, cls).__new__(cls)
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
        for item in self.history.items():
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
            # partial  (if .keywords is None print {} instead)
            elements = undo_entry.func, undo_entry.args, undo_entry.keywords or {}
            info += 'undo: ' + ' '.join([str(elem) for elem in elements]) + '\n'
        except AttributeError:  # normal python function
            info += 'undo: ' + undo_entry.__name__ + '\n'

        redo_entry = entry['redo']
        try:
            elements = redo_entry.func, redo_entry.args, redo_entry.keywords or {}  # partial
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
        return iter(self._dict.values())

    def __dir__(self):
        return list(self._dict.keys())


def inheritdocstring(name, bases, attrs):
    """Use as metaclass to inherit class and method docstrings from parent.
    Adapted from http://stackoverflow.com/questions/13937500/inherit-a-parent-class-docstring-as-doc-attribute"""
    temp = type('temporaryclass', bases, {})
    if '__doc__' not in attrs or not attrs["__doc__"]:
        # create a temporary 'parent' to (greatly) simplify the MRO search
        for cls in inspect.getmro(temp):
            if cls.__doc__ is not None:
                attrs['__doc__'] = cls.__doc__
                break

    for attr_name, attr in attrs.items():
        if not attr.__doc__:
            for cls in inspect.getmro(temp):
                try:
                    if getattr(cls, attr_name).__doc__ is not None:
                        attr.__doc__ = getattr(cls, attr_name).__doc__
                        break
                except (AttributeError, TypeError):
                    continue

    return type(name, bases, attrs)


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


def flatten(input_list):
    return [item for sublist in input_list for item in sublist]


def generate_colors(n):
    hsv_tuples = [(v * 1.0 / n, 0.5, 0.5) for v in range(n)]
    color_map = {}
    for i in range(n):
        rgb = colorsys.hsv_to_rgb(*hsv_tuples[i])
        color = tuple(int(channel * 256) for channel in rgb)
        color_map[i] = '#%02x%02x%02x' % color
    return color_map


def memoize(function, memo={}):
    def wrapper(*args):
        if args in memo:
            return memo[args]
        else:
            rv = function(*args)
            memo[args] = rv
            return rv

    return wrapper


def get_system_info():
    package_info = list()
    for dist in pkg_resources.working_set:
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
    return pandas.io.formats.console.in_ipython_frontend()


def str_to_valid_variable_name(s):
    """Adapted from http://stackoverflow.com/a/3303361/280182"""

    # Remove invalid characters
    s = re.sub('[^0-9a-zA-Z_]', '_', s)

    # Remove leading characters until we find a letter or underscore
    s = re.sub('^[^a-zA-Z_]+', '', s)

    return s


def zip_repeat(long_iter, short_iter):
    """
    Zips two iterable objects but repeats the second one if it is shorter than the first one.

    Parameters
    ----------
    long_iter: iterable

    short_iter: iterable

    Returns
    -------
    generator
    """
    for i, j in zip(long_iter, itertools.cycle(short_iter)):
        yield i, j


def pick_one(iterable):
    """
    Helper function that returns an element of an iterable (it the iterable is ordered this will be the first
    element).
    """
    it = iter(iterable)
    return next(it)


def reduce_reaction_set(reaction_set, groups):
    """
    Reduces a set of reactions according to a number of groups of reactions.
    The reduction will be performed so that the resulting set will contain no more than 1 reaction from each group.
    Reactions that are not in any of the groups will remain in the set.

    Parameters
    ----------
    reaction_set: Set
    groups: Iterable of sets

    Returns
    -------
    Set
    """
    reaction_set = set(reaction_set)  # Make a shallow copy
    result = []
    for group in groups:
        intersection = group & reaction_set
        if intersection:  # If any elements of group are in reaction_set, add one of these to result
            result.append(pick_one(intersection))
            reaction_set = reaction_set - intersection
    result = set(result) | reaction_set  # Add the remaining reactions to result
    return result


def decompose_reaction_groups(reaction_groups, reactions):
    """

    reaction_groups : list
        A list with dictionaries (element: relative_coefficient)
    reactions : list, set, tuple
        A collection of reactions.

    Returns
    -------
    list
        A list of all possible group substitutions.

    """
    to_keep = []
    to_replace = {}
    for element in reactions:
        for g in reaction_groups:
            if element in g:
                to_replace[element] = g.keys()
                break
        if element not in to_replace:
            to_keep.append(element)

    return (list(combo) + to_keep for combo in itertools.product(*to_replace.values()))


def current_solver_name(model):
    """Give a string representation for an optlang interface.

    Parameters
    ----------
    model : cobra.Model
        A model

    Returns
    -------
    string
       The name of the interface as a string
    """
    interface = model.solver.interface.__name__
    return re.sub(r"optlang.|.interface", "", interface)
