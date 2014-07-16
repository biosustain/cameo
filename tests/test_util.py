# Copyright (c) 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.
# See LICENSE for details.

import unittest
from functools import partial

from cameo.util import TimeMachine, generate_colors, Singleton


class TimeMachineTestCase(unittest.TestCase):
    def setUp(self):
        self.tm = TimeMachine()

    def test_one_change_list(self):
        l = [1, 2, 3, 4]
        self.tm(do=partial(l.append, 5), undo=l.pop)
        self.assertEqual(l, [1, 2, 3, 4, 5])
        self.tm.reset()
        self.assertEqual(l, [1, 2, 3, 4])

    def test_str_handles_different_types_of_stored_operations(self):
        def normal_function():
            pass

        partial_function = partial(str, 1)
        self.tm(do=normal_function, undo=partial_function)
        self.assertEqual(self.tm.__str__().split('\n')[2:-1], ["undo: <type 'str'> (1,) None", 'redo: normal_function'])


class TestUtils(unittest.TestCase):
    def test_color_generation(self):
        for i in xrange(1, 100):
            color_map = generate_colors(i)
            self.assertEqual(len(color_map), i)
            self.assertEqual(len(color_map), len(set(color_map.values())))


class TestSingleton(unittest.TestCase):
    def test_singleton(self):
        s1 = Singleton()
        s2 = Singleton()
        self.assertEqual(s1, s2)