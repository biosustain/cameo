# Copyright (c) 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.
# See LICENSE for details.

import unittest
import subprocess
from time import sleep

from IPython.parallel import Client, interactive

try:
    from cameo.parallel import MultiprocessingView


    @interactive
    def to_the_power_of_2_interactive(arg):
        return arg ** 2


    def to_the_power_of_2(arg):
        return arg ** 2


    class TestMultiprocessingView(unittest.TestCase):
        def setUp(self):
            self.view = MultiprocessingView()

        def test_map(self):
            self.assertEqual(self.view.map(to_the_power_of_2, range(100)), map(lambda x: x ** 2, range(100)))
except ImportError:
    print "Skipping MultiprocessingView tests ..."


# class TestIPythonParallelView(unittest.TestCase):
#     def setUp(self):
#         try:
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

    # def func(arg):
    #     return arg**2

    # client = Client()
    # view = client.direct_view()
    # view.block = False
    # print view.map_sync(to_the_power_of_2, range(100))