# Copyright 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.
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

import unittest
from functools import partial
from itertools import chain

from cobra import Metabolite
from six.moves import range

from cameo.network_analysis.util import distance_based_on_molecular_formula
from cameo.util import TimeMachine, generate_colors, Singleton, partition, RandomGenerator


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
        self.assertEqual(self.tm.__str__().split('\n')[2:-1],
                         ["undo: " + str(str) + " (1,) {}", 'redo: normal_function'])

    def test_with_statement(self):
        l = [1, 2, 3, 4]
        with TimeMachine() as tm:
            tm(do=partial(l.append, 33), undo=partial(l.pop))
            tm(do=partial(l.append, 66), undo=partial(l.pop))
            tm(do=partial(l.append, 99), undo=partial(l.pop))
        self.assertEqual(l, [1, 2, 3, 4])


class TestRandomGenerator(unittest.TestCase):
    def setUp(self):
        self.seed = 1234

    def test_random(self):
        random = RandomGenerator()
        for _ in range(1000):
            self.assertGreaterEqual(random.random(), 0)
            self.assertLessEqual(random.random(), 1)

    def test_randint(self):
        random = RandomGenerator()
        lower = 0
        upper = 10
        for _ in range(10000):
            self.assertGreaterEqual(random.randint(lower, upper), lower)
            self.assertLessEqual(random.randint(lower, upper), upper)

        lower = -10
        upper = 100
        for _ in range(10000):
            self.assertGreaterEqual(random.randint(lower, upper), lower)
            self.assertLessEqual(random.randint(lower, upper), upper)

        lower = 5
        upper = 21
        for _ in range(10000):
            self.assertGreaterEqual(random.randint(lower, upper), lower)
            self.assertLessEqual(random.randint(lower, upper), upper)

        lower = -5
        upper = 5
        for _ in range(10000):
            self.assertGreaterEqual(random.randint(lower, upper), lower)
            self.assertLessEqual(random.randint(lower, upper), upper)

    def test_seeded_methods(self):
        random = RandomGenerator()

        random.seed(self.seed)
        value = random.random()
        random.seed(self.seed)
        self.assertEqual(value, random.random())

        random.seed(self.seed)
        value = random.randint(1, 10)
        random.seed(self.seed)
        self.assertEqual(value, random.randint(1, 10))

        random.seed(self.seed)
        population = [1, 2, 3, 4, 5]
        value = random.sample(population, 2)
        random.seed(self.seed)
        self.assertEqual(value, random.sample(population, 2))

        random.seed(self.seed)
        value = random.uniform()
        random.seed(self.seed)
        self.assertEqual(value, random.uniform())


class TestUtils(unittest.TestCase):
    def test_color_generation(self):
        for i in range(1, 100):
            color_map = generate_colors(i)
            self.assertEqual(len(color_map), i)
            self.assertEqual(len(color_map), len(set(color_map.values())))

    def test_partition(self):
        chunks = 3
        iterables = [
            [1, 2, 3, 4, 5, 6, 7, 8, 9],
            {5, 3, 8, 3, 8, 5, 8, 0, 10, 11, 15},
            range(29)
        ]
        for fixture in iterables:
            test_output = partition(fixture, chunks)
            self.assertEqual(len(fixture), sum(map(len, test_output)))
            self.assertEqual(len(test_output), chunks)
            self.assertEqual(list(fixture), list(chain(*test_output)))
            for out_chunk in test_output:
                self.assertTrue(set(out_chunk).issubset(set(fixture)))

        bad_input = 5
        self.assertRaises(TypeError, partition, bad_input, chunks)

    @unittest.skip  # Development API changes to cobra.core.Metabolite
    def test_distance_based_on_molecular_formula(self):  # from network_analysis.util
        met1 = Metabolite("H2O", formula="H2O")
        met2 = Metabolite("H2O2", formula="H2O2")
        met3 = Metabolite("C6H12O6", formula="C6H12O6")

        self.assertEqual(distance_based_on_molecular_formula(met1, met2, normalize=False), 1)
        self.assertEqual(distance_based_on_molecular_formula(met1, met2, normalize=True), 1. / 7)

        self.assertEqual(distance_based_on_molecular_formula(met2, met3, normalize=False), 20)
        self.assertEqual(distance_based_on_molecular_formula(met2, met3, normalize=True), 20. / 28)

        self.assertEqual(distance_based_on_molecular_formula(met1, met3, normalize=False), 21)
        self.assertEqual(distance_based_on_molecular_formula(met1, met3, normalize=True), 21. / 27)


class TestSingleton(unittest.TestCase):
    def test_singleton(self):
        s1 = Singleton()
        s2 = Singleton()
        self.assertIs(s1, s2)


if __name__ == "__main__":
    import nose

    nose.runmodule()
