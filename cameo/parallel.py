# Copyright 2013 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import Queue
from multiprocessing import Pool, cpu_count
from cameo.util import Singleton


class MultiprocessingView(Singleton):

    """docstring for Parallel"""

    def __init__(self, *args, **kwargs):
        self._args = args
        self._kwargs = kwargs
        self.pool = None

    def map(self, *args, **kwargs):
        if self.pool is None:
            self.pool = Pool(*self._args, **self._kwargs)
        return self.pool.map(*args, **kwargs)

    def apply(self, func, *args, **kwargs):
        if self.pool is None:
            self.pool = Pool(*self._args, **self._kwargs)
        return self.pool.apply(func, args=args, **kwargs)

    def apply_async(self, func, *args, **kwargs):
        if self.pool is None:
            self.pool = Pool(*self._args, **self._kwargs)
        self.pool.apply_async(func, args=args, **kwargs)

    def __len__(self):
        return cpu_count()

    def shutdown(self):
        if not self.pool is None:
            self.pool.terminate()


try:
    import redis
    import cPickle as pickle

    class RedisQueue(object):
        """
        Queue with Redis Backend
        http://peter-hoffmann.com/2012/python-simple-queue-redis-queue.html
        """

        MAX_REDIS_LIST_SIZE = 4294967295L

        default_host = "localhost"
        default_port = "6379"
        default_db = 0

        def __init__(self, name, maxsize=0, namespace='queue', **connection_args):
            """The default connection parameters are: host='localhost', port=6379, db=0"""
            if maxsize <= 0:
                maxsize = self.MAX_REDIS_LIST_SIZE
            self._maxsize = maxsize
            self._connection_args = {
                'host': self.default_host,
                'port': self.default_port,
                'db': self.default_db
            }
            for key, val in connection_args.iteritems():
                self._connection_args[key] = val
            self._db = redis.Redis(**self._connection_args)
            self._key = '%s:%s' % (namespace, name)

        def __getstate__(self):
            return {
                '_maxsize': self._maxsize,
                '_connection_args': self._connection_args,
                '_key': self._key
            }

        def __setstate__(self, d):
            self._maxsize = d['_maxsize']
            self._connection_args = d['_connection_args']
            self._key = d['_key']
            self._db = redis.Redis(**self._connection_args)

        def __len__(self):
            return self.length()

        def length(self):
            return self._db.llen(self._key)

        def empty(self):
            return self.length() == 0

        def put(self, item):
            if self.length() >= self._maxsize:
                raise Queue.Full

            item = pickle.dumps(item)

            self._db.rpush(self._key, item)

        def put_nowait(self, item):
            """
            Same as put, for backward compatibility with Python Queue.

            See also
            --------
            put

            :param item: object to put in the queue

            """
            self.put(item)

        def get(self, block=True, timeout=None):
            if block:
                item = self._db.blpop(self._key, timeout=timeout)
            else:
                item = self._db.lpop(self._key)

            if item:
                if isinstance(item, str):
                    return pickle.loads(item)
                else:
                    return pickle.loads(item[1])
            else:
                raise Queue.Empty

        def get_nowait(self):
            """Equivalent to get(False)."""
            return self.get(False)

        def __del__(self):
            self._db.delete(self._key)
            self._db.connection_pool.disconnect()

except ImportError:
    pass


class SequentialView(object):
    def map(self, *args, **kwargs):
        return map(*args, **kwargs)

    def apply(self, func, *args, **kwargs):
        return func(*args, **kwargs)

    def apply_async(self, func, *args, **kwargs):
        return func(*args, **kwargs)

    def __len__(self):
        return 1

    def shutdown(self):
        pass


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