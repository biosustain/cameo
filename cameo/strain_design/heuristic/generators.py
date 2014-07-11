# Copyright 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
from inspyred.ec.generators import diversify


def set_generator(random, args):
    representation = args.get('representation')
    max_size = args.get('max_candidate_size', 9)
    variable_size = args.get('variable_candidate_size', True)
    if variable_size:
        size = random.randint(1, max_size)
    else:
        size = max_size
    candidate = random.sample(xrange(len(representation)), size)
    return candidate


@diversify
def unique_set_generator(random, args):
    representation = args.get('representation')
    max_size = args.get('max_candidate_size', 9)
    variable_size = args.get('variable_candidate_size', True)
    if variable_size:
        size = random.randint(1, max_size)
    else:
        size = max_size
    candidate = random.sample(xrange(len(representation)), size)
    return candidate


