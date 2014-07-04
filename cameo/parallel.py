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
import Queue

from pathos.multiprocessing import ProcessingPool
import redis


class MultiprocessingView(object):

    """docstring for Parallel"""

    def __init__(self, *args, **kwargs):
        # super(MultiprocessingView, self).__init__(*args, **kwargs)
        self.pool = ProcessingPool(*args, **kwargs)

    def map(self, *args, **kwargs):
        return self.pool.map(*args, **kwargs)

    def apply(self, *args, **kwargs):
        raise NotImplementedError

    def __len__(self):
        return self.pool.nodes

    def shutdown(self):
        self.pool.terminate()


class SequentialView(object):

    def map(self, *args, **kwargs):
        return map(*args, **kwargs)

    def apply(self, *args, **kwargs):
        return apply(*args, **kwargs)

    def __len__(self):
        return 1

    def shutdown(self):
        pass


class RedisQueue(object):
    """
    Queue with Redis Backend
    http://peter-hoffmann.com/2012/python-simple-queue-redis-queue.html
    """

    MAX_REDIS_LIST_SIZE = 4294967295L

    def __init__(self, name, maxsize=0, namespace='queue', **kwargs):
        """The default connection parameters are: host='localhost', port=6379, db=0"""
        if maxsize <= 0:
            maxsize = self.MAX_REDIS_LIST_SIZE
        self._maxsize = maxsize
        self._db = redis.Redis(**kwargs)
        self._key = '%s:%s' % (namespace, name)

    def __len__(self):
        return self._db.llen(self._key)

    def empty(self):
        return self.qsize() == 0

    def put(self, item):
        if len(self) > self._maxsize:
            raise Queue.Full

        self._db.rpush(self._key, item)

    def get(self, block=True, timeout=None):
        if block:
            item = self._db.blpop(self._key, timeout=timeout)
        else:
            item = self._db.lpop(self._key)

        if item:
            item = item[1]
        else:
            raise Queue.Empty
        return item

    def get_nowait(self):
        """Equivalent to get(False)."""
        return self.get(False)


if __name__ == '__main__':

    from time import sleep

    class FunctionObject(object):

        """docstring for FunctionObject"""

        def __init__(self):
            super(FunctionObject, self).__init__()

        def __call__(self, arg):
            sleep(.1)
            return arg ** 2

    view = MultiprocessingView(processes=4)
    print view.map(FunctionObject(), range(0, 100))

