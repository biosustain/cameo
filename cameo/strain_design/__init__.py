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
from cameo.core.result import Result


class StrainDesignMethod(object):
    def __init__(self, *args, **kwargs):
        super(StrainDesignMethod, self).__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):
        self.run(*args, **kwargs)

    def run(self, *args, **kwargs):
        raise NotImplementedError


class StrainDesign(object):
    def __init__(self, knockouts=None, knock_ins=None, over_expression=None, down_regulation=None):
        self.knockouts = knockouts
        self.knock_ins = knock_ins
        self.over_expression = over_expression
        self.down_regulation = down_regulation


class StrainDesignResult(Result):
    def __init__(self, *args, **kwargs):
        super(StrainDesignResult, self).__init__(*args, **kwargs)

    def __iter__(self):
        """

        Returns an iterator that yields StrainDesign objects.
        """
        raise NotImplementedError