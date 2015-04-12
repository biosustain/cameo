# Copyright 2015 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
from progressbar import ProgressBar, Bar, RotatingMarker, Percentage
from cameo.reporters import AbstractReporter


class ProgressReporter(AbstractReporter):
    widgets = [' ', Bar(marker=RotatingMarker()), ' ', Percentage()]

    def __init__(self, label, size, initial=0):
        self.size = size
        self.initial = initial
        self.pb = ProgressBar(widgets=[label] + self.widgets, maxval=self.size - self.initial)
        self.started = False

    def __call__(self, current_value):
        if not self.started:
            self.pb.start()
            self.started = True

        value = current_value - self.initial
        self.pb.update(int(value))

    def close(self):
        self.pb.finish()
        self.started = False
