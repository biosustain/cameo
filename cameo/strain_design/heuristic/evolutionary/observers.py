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

from IProgress.progressbar import ProgressBar
from IProgress.widgets import Percentage, Bar


class ProgressObserver(object):
    """
    Progress bar to in command line. It keeps track of the progress during heuristic optimization.
    """
    __name__ = "Progress Observer"

    def __init__(self):
        self.progress = None

    def __call__(self, population, num_generations, num_evaluations, args):
        if self.progress is None:
            self.max_evaluations = args.get('max_evaluations', 50000)
            self.progress = ProgressBar(maxval=self.max_evaluations,
                                        widgets=["Running Optimization", Bar(), Percentage()])
            self.progress.start()

        if num_evaluations % args.get('n', 1) == 0:
            if num_evaluations > self.max_evaluations:
                self.progress.update(self.max_evaluations)
            else:
                self.progress.update(num_evaluations)

    def reset(self):
        self.progress = None

    def end(self):
        self.progress.finish()
