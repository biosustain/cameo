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

from ipython_notebook_utils import ProgressBar as IPythonProgressBar
from progressbar import ProgressBar as CLIProgressBar


class CLIProgressObserver():
    """
    Progress bar to in command line
    """
    __name__ = "CLI progress"

    def __init__(self):
        self.progress = None

    def __call__(self, population, num_generations, num_evaluations, args):
        if self.progress is None:
            self.progress = CLIProgressBar(args.get('max_evaluations', 50000))
            self.progress.start()

        if num_evaluations % args.get('n', 1) == 0:
            self.progress.update(num_evaluations)

    def reset(self):
        self.progress = None


class IPythonNotebookProgressObserver():
    """
    Progress bar to use within IPython Notebook
    """
    __name__ = "IPython Notebook progress"

    def __init__(self):
        self.progress = IPythonProgressBar()

    def __call__(self, population, num_generations, num_evaluations, args):
        if num_evaluations % args.get('n', 1) == 0:
            p = (float(num_evaluations) / float(args.get('max_evaluations', 50000))) * 100.0
            self.progress.set(p)

    def reset(self):
        self.progress.start()