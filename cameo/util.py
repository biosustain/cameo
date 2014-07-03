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
        for uuid, entry in self.history.iteritems():
            info += datetime.fromtimestamp(entry['unix_epoch']
                                           ).strftime('%Y-%m-%d %H:%M:%S') + '\n'
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

    def __repr__(self):
        return self.__str__()

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
        color = colorsys.hsv_to_rgb(*hsv_tuples[i])
        color = tuple(map(lambda v: v*256, color))
        color_map[i] = '#%02x%02x%02x' % color
    return color_map