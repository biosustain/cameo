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
from uuid import uuid1
from time import time
from datetime import datetime
import colorsys
import numpy as np
from numpy.linalg import svd


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


# Taken from http://wiki.scipy.org/Cookbook/RankNullspace
def nullspace(A, atol=1e-13, rtol=0):
    """Compute an approximate basis for the nullspace of A.

    The algorithm used by this function is based on the singular value
    decomposition of `A`.

    Parameters
    ----------
    A : ndarray
        A should be at most 2-D.  A 1-D array with length k will be treated
        as a 2-D with shape (1, k)
    atol : float
        The absolute tolerance for a zero singular value.  Singular values
        smaller than `atol` are considered to be zero.
    rtol : float
        The relative tolerance.  Singular values less than rtol*smax are
        considered to be zero, where smax is the largest singular value.

    If both `atol` and `rtol` are positive, the combined tolerance is the
    maximum of the two; that is::
        tol = max(atol, rtol * smax)
    Singular values smaller than `tol` are considered to be zero.

    Return value
    ------------
    ns : ndarray
        If `A` is an array with shape (m, k), then `ns` will be an array
        with shape (k, n), where n is the estimated dimension of the
        nullspace of `A`.  The columns of `ns` are a basis for the
        nullspace; each element in numpy.dot(A, ns) will be approximately
        zero.
    """

    A = np.atleast_2d(A)
    u, s, vh = svd(A)
    tol = max(atol, rtol * s[0])
    nnz = (s >= tol).sum()
    ns = vh[nnz:].conj().T
    return ns