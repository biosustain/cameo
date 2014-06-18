# Copyright (c) 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.
# See LICENSE for details.

import unittest
from functools import partial

from cameo.util import TimeMachine


class TimeMachineTestCase(unittest.TestCase):
    def setUp(self):
        self.tm = TimeMachine()

    def test_one_change_list(self):
        l = [1, 2, 3, 4]
        self.tm(do=partial(l.append, 5), undo=l.pop)
        self.assertEqual(l, [1, 2, 3, 4, 5])
        self.tm.reset()
        self.assertEqual(l, [1, 2, 3, 4, 5])
