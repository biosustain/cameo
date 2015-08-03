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

import six.moves.queue

import unittest
from cameo.parallel import SequentialView
import subprocess
from time import sleep
from multiprocessing import cpu_count

from IPython.parallel import Client, interactive
from six.moves import range

SOLUTION = [x ** 2 for x in range(100)]


@interactive
def to_the_power_of_2_interactive(arg):
    return arg ** 2


def to_the_power_of_2(arg):
    return arg ** 2


class TestSequentialView(unittest.TestCase):
    def setUp(self):
        self.view = SequentialView()

    def test_map(self):
        self.assertEqual(self.view.map(to_the_power_of_2, list(range(100))), SOLUTION)

    def test_apply(self):
        for i in range(100):
            self.assertEqual(self.view.apply(to_the_power_of_2, i), SOLUTION[i])


try:
    from cameo.parallel import MultiprocessingView

    class TestMultiprocessingView(unittest.TestCase):
        def setUp(self):
            self.view = MultiprocessingView()

        def test_map(self):
            self.assertEqual(self.view.map(to_the_power_of_2, list(range(100))), SOLUTION)

        def test_apply(self):
            for i in range(100):
                self.assertEqual(self.view.apply(to_the_power_of_2, i), SOLUTION[i])

        def test_length(self):
            self.assertEqual(len(self.view), cpu_count())

            view = MultiprocessingView(4)
            self.assertEqual(len(view), 4)

            view = MultiprocessingView(processes=3)
            self.assertEqual(len(view), 3)

except ImportError:
    print("Skipping MultiprocessingView tests ...")

try:
    from cameo.parallel import RedisQueue

    class TestRedisQueue(unittest.TestCase):
        def test_queue_size(self):
            queue = RedisQueue("test-queue-size-1", maxsize=1)
            queue.put(1)
            self.assertRaises(six.moves.queue.Full, queue.put, 1)

            queue = RedisQueue("test-queue-size-2", maxsize=2)
            queue.put(1)
            queue.put(1)
            self.assertRaises(six.moves.queue.Full, queue.put, 1)
            queue.get()
            queue.get()
            self.assertRaises(six.moves.queue.Empty, queue.get_nowait)

        def test_queue_objects(self):
            queue = RedisQueue("test-queue", maxsize=100)
            # put int
            queue.put(1)
            v = queue.get_nowait()
            self.assertEqual(v, 1)
            self.assertIsInstance(v, int)

            # put str
            queue.put("a")
            v = queue.get_nowait()
            self.assertEqual(v, "a")
            self.assertIsInstance(v, str)

            # put float
            queue.put(1.)

            v = queue.get_nowait()
            self.assertEqual(v, 1.)
            self.assertIsInstance(v, float)

            # put list
            queue.put([1, 3, 4, 5, "a", "b", "c", 1., 2., 3.])
            v = queue.get_nowait()
            self.assertEqual(v, [1, 3, 4, 5, "a", "b", "c", 1., 2., 3.])
            self.assertIsInstance(v, list)

            # put dict
            queue.put({"x": "y"})
            v = queue.get_nowait()
            self.assertEqual(v, {"x": "y"})
            self.assertIsInstance(v, dict)

        def test_queue_len(self):
            queue = RedisQueue("test-queue-len", maxsize=100)
            self.assertEqual(queue.length, 0)
            queue.put(1)
            self.assertEqual(queue.length, 1)
            queue.put(1)
            self.assertEqual(queue.length, 2)
            queue.put(1)
            self.assertEqual(queue.length, 3)
            queue.get_nowait()
            self.assertEqual(queue.length, 2)
            queue.get_nowait()
            self.assertEqual(queue.length, 1)
            queue.get_nowait()
            self.assertEqual(queue.length, 0)

except ImportError:
    print("Skipping MultiprocessingView tests ...")


# class TestIPythonParallelView(unittest.TestCase):
# def setUp(self):
# try:
#             subprocess.check_output(["ipcluster", "start", "--daemonize"])
#             sleep(6)
#         except subprocess.CalledProcessError:
#             subprocess.check_output(["ipcluster", "stop"])
#             subprocess.check_output(["ipcluster", "start", "--daemonize"])
#             sleep(6)
#
#     def test_map_load_balanced_view(self):
#         client = Client()
#         self.view = client.load_balanced_view()
#         self.view.block = False
#         self.assertEqual(self.view.map_sync(to_the_power_of_2_interactive, range(100)),
#                          map(lambda x: x ** 2, range(100)))
#
#     def test_map_direct_view(self):
#         client = Client()
#         self.view = client.direct_view()
#         self.view.block = False
#         self.assertEqual(self.view.map_sync(to_the_power_of_2_interactive, range(100)),
#                          map(lambda x: x ** 2, range(100))
#
#     def tearDown(self):
#         try:
#             subprocess.check_output(["ipcluster", "stop"])
#             sleep(4)
#         except subprocess.CalledProcessError:
#             pass


if __name__ == '__main__':
    import nose

    nose.runmodule()
