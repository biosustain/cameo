# Copyright (c) 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.
# See LICENSE for details.
import Queue

import unittest
from cameo.parallel import SequentialView
import subprocess
from time import sleep

from IPython.parallel import Client, interactive

SOLUTION = map(lambda x: x ** 2, range(100))


@interactive
def to_the_power_of_2_interactive(arg):
    return arg ** 2


def to_the_power_of_2(arg):
    return arg ** 2


class TestSequentialView(unittest.TestCase):
    def setUp(self):
        self.view = SequentialView()

    def test_map(self):
        self.assertEqual(self.view.map(to_the_power_of_2, range(100)), SOLUTION)

    def test_apply(self):
        for i in xrange(100):
            self.assertEqual(self.view.apply(to_the_power_of_2, i), SOLUTION[i])


try:
    from cameo.parallel import MultiprocessingView

    class TestMultiprocessingView(unittest.TestCase):
        def setUp(self):
            self.view = MultiprocessingView()

        def test_map(self):
            self.assertEqual(self.view.map(to_the_power_of_2, range(100)), SOLUTION)

        def test_apply(self):
            for i in xrange(100):
                self.assertEqual(self.view.apply(to_the_power_of_2, i), SOLUTION[i])

except ImportError:
    print "Skipping MultiprocessingView tests ..."

try:
    from cameo.parallel import RedisQueue

    class TestRedisQueue(unittest.TestCase):
        def test_queue_size(self):
            queue = RedisQueue("test-queue-size-1", maxsize=1)
            queue.put(1)
            self.assertRaises(Queue.Full, queue.put, 1)

            queue = RedisQueue("test-queue-size-2", maxsize=2)
            queue.put(1)
            queue.put(1)
            self.assertRaises(Queue.Full, queue.put, 1)
            queue.get()
            queue.get()
            self.assertRaises(Queue.Empty, queue.get_nowait)

        def test_queue_len(self):
            queue = RedisQueue("test-queue-len", maxsize=100)
            self.assertEqual(queue.length(), 0)
            queue.put(1)
            self.assertEqual(queue.length(), 1)
            queue.put(1)
            self.assertEqual(queue.length(), 2)
            queue.put(1)
            self.assertEqual(queue.length(), 3)
            queue.get_nowait()
            self.assertEqual(queue.length(), 2)
            queue.get_nowait()
            self.assertEqual(queue.length(), 1)
            queue.get_nowait()
            self.assertEqual(queue.length(), 0)

except ImportError:
    print "Skipping MultiprocessingView tests ..."


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