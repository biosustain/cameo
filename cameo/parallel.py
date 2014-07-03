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

from pathos.multiprocessing import ProcessingPool


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
        return self.pool._processes

    def shutdown(self):
        pass


class SequentialView(object):

    def map(self, *args, **kwargs):
        return map(*args, **kwargs)

    def apply(self, *args, **kwargs):
        return apply(*args, **kwargs)

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
