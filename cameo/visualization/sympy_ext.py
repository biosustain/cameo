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

from sympy import And, Or


def _str_or(self):
    return "(" + " OR ".join([str(v) for v in self.args]) + ")"


Or.__str__ = _str_or


def _str_and(self):
    return "(" + " AND ".join([str(v) for v in self.args]) + ")"


And.__str__ = _str_and
