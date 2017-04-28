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

import os
import warnings
from multiprocessing import cpu_count

import pytest

import six.moves.queue
from cameo.parallel import SequentialView
from six.moves import range

views = [SequentialView()]

try:
    from cameo.parallel import MultiprocessingView
    views.append(MultiprocessingView())
except ImportError:
    MultiprocessingView = None

try:
    from cameo.parallel import RedisQueue
except ImportError:
    RedisQueue = None

try:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from IPython.parallel import Client, interactive
except ImportError:
    try:
        from ipyparallel import Client, interactive
    except ImportError:
        def interactive(f):
            return f

SOLUTION = [x ** 2 for x in range(100)]


if os.getenv('REDIS_PORT_6379_TCP_ADDR'):
    REDIS_HOST = os.getenv('REDIS_PORT_6379_TCP_ADDR')  # wercker
else:
    REDIS_HOST = 'localhost'


@interactive
def to_the_power_of_2_interactive(arg):
    return arg ** 2


def to_the_power_of_2(arg):
    return arg ** 2


class TestView:
    @pytest.mark.parametrize('view', views)
    def test_map(self, view):
        assert view.map(to_the_power_of_2, list(range(100))) == SOLUTION

    @pytest.mark.parametrize('view', views)
    def test_apply(self, view):
        for i in range(100):
            assert view.apply(to_the_power_of_2, i) == SOLUTION[i]

    @pytest.mark.skipif(not MultiprocessingView, reason="no multiprocessing available")
    def test_length(self):
        view = MultiprocessingView()
        assert len(view) == cpu_count()

        view = MultiprocessingView(4)
        assert len(view) == 4

        view = MultiprocessingView(processes=3)
        assert len(view) == 3


@pytest.mark.skipif(not RedisQueue, reason='no redis queue available')
class TestRedisQueue:

    def test_queue_size(self):
        print(REDIS_HOST)
        print(os.getenv('REDIS_PORT_6379_TCP_ADDR'))
        queue = RedisQueue("test-queue-size-1", maxsize=1, host=REDIS_HOST)
        queue.put(1)
        with pytest.raises(six.moves.queue.Full):
            queue.put(1)

        queue = RedisQueue("test-queue-size-2", maxsize=2, host=REDIS_HOST)
        queue.put(1)
        queue.put(1)
        with pytest.raises(six.moves.queue.Full):
            queue.put(1)
        queue.get()
        queue.get()
        with pytest.raises(six.moves.queue.Empty):
            queue.get_nowait()

    def test_queue_objects(self):
        queue = RedisQueue("test-queue", maxsize=100, host=REDIS_HOST)
        # put int
        queue.put(1)
        v = queue.get_nowait()
        assert v == 1
        assert isinstance(v, int)

        # put str
        queue.put("a")
        v = queue.get_nowait()
        assert v == "a"
        assert isinstance(v, six.string_types)

        # put float
        queue.put(1.)

        v = queue.get_nowait()
        assert v == 1.
        assert isinstance(v, float)

        # put list
        queue.put([1, 3, 4, 5, "a", "b", "c", 1., 2., 3.])
        v = queue.get_nowait()
        assert v == [1, 3, 4, 5, "a", "b", "c", 1., 2., 3.]
        assert isinstance(v, list)

        # put dict
        queue.put({"x": "y"})
        v = queue.get_nowait()
        assert v == {"x": "y"}
        assert isinstance(v, dict)

    def test_queue_len(self):
        queue = RedisQueue("test-queue-len", maxsize=100, host=REDIS_HOST)
        assert queue.length == 0
        queue.put(1)
        assert queue.length == 1
        queue.put(1)
        assert queue.length == 2
        queue.put(1)
        assert queue.length == 3
        queue.get_nowait()
        assert queue.length == 2
        queue.get_nowait()
        assert queue.length == 1
        queue.get_nowait()
        assert queue.length == 0
