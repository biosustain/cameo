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
from cameo.parallel import RedisQueue
from uuid import uuid4

import logging

__all__ = ['MultiprocessingMigrator']


logger = logging.getLogger(__name__)


class MultiprocessingMigrator(object):
    """Migrate among processes on multiple machines. This is possible by
    having a Queue implemented on redis.

    This callable class allows individuals to migrate from one process
    to another. It maintains a queue of migrants whose maximum length
    can be fixed via the ``max_migrants`` parameter in the constructor.
    If the number of migrants in the queue reaches this value, new migrants
    are not added until earlier ones are consumed. The unreliability of a
    multiprocessing environment makes it difficult to provide guarantees.
    However, migrants are theoretically added and consumed at the same rate,
    so this value should determine the "freshness" of individuals, where
    smaller queue sizes provide more recency.

    An optional keyword argument in ``args`` requires the migrant to be
    evaluated by the current evolutionary computation before being inserted
    into the population. This can be important when different populations
    use different evaluation functions and you need to be able to compare
    "apples with apples," so to speak.

    The migration takes the current individual *I* out of the queue, if
    one exists. It then randomly chooses an individual *E* from the population
    to insert into the queue. Finally, if *I* exists, it replaces *E* in the
    population (re-evaluating fitness if necessary). Otherwise, *E* remains in
    the population and also exists in the queue as a migrant.

    Optional keyword arguments in args:

    - *evaluate_migrant* -- should new migrants be evaluated before
      adding them to the population (default False)

    Arguments
    ---------
    max_migrants: int
        Number of migrants in the queue at the same time.
    connection_kwargs: keyword arguments:
        see parallel.RedisQueue

    """

    def __init__(self, max_migrants=1, **connection_kwargs):
        self.max_migrants = max_migrants
        self.migrants = RedisQueue(uuid4(), **connection_kwargs)
        self.__name__ = self.__class__.__name__

    def __call__(self, random, population, args):
        evaluate_migrant = args.setdefault('evaluate_migrant', False)
        migrant_index = random.randint(0, len(population) - 1)
        old_migrant = population[migrant_index]
        try:
            migrant = self.migrants.get(block=False)
            logger.debug("Robinson Crusoe arrives on an island")
            if evaluate_migrant:
                fit = args["_ec"].evaluator([migrant.candidate], args)
                migrant.fitness = fit[0]
                args["_ec"].num_evaluations += 1
            population[migrant_index] = migrant
        except six.moves.queue.Empty:
            logger.debug("Empty queue")
        try:
            logger.debug("Robinson Crusoe leaves an island")
            self.migrants.put_nowait(old_migrant)
        except six.moves.queue.Full:
            logger.debug("Full queue")
        return population
