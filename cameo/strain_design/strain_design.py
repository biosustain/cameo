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
import six

from cameo import ui
from pandas import DataFrame

from cameo.core.result import Result


class StrainDesignMethod(object):
    def __init__(self, *args, **kwargs):
        super(StrainDesignMethod, self).__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):
        self.run(*args, **kwargs)

    def run(self, *args, **kwargs):
        raise NotImplementedError


class StrainDesign(object):
    def __init__(self, knockouts=None, knock_ins=None, over_expression=None, down_regulation=None,
                 manipulation_type="gene"):
        self.knockouts = knockouts
        self.knock_ins = knock_ins
        self.over_expression = over_expression
        self.down_regulation = down_regulation
        self.manipulation_type = manipulation_type

    def __repr__(self):
        return "".join([ui.delta() + ko for ko in self.knockouts] +
                       [ui.knockin() + ki for ki in self.knock_ins] +
                       [ui.upreg(coeff) + oe for oe, coeff in six.iteritems(self.over_expression)] +
                       [ui.downreg(coeff) + dr for dr, coeff in six.iteritems(self.down_regulation)])


class StrainDesignResult(Result):

    __method_name__ = ""

    def __init__(self, *args, **kwargs):
        super(StrainDesignResult, self).__init__(*args, **kwargs)

    def __len__(self):
        """
        Return the number of designs found.
        """
        raise NotImplementedError

    def __iter__(self):
        """
        Returns an iterator that yields StrainDesign objects.
        """
        raise NotImplementedError

    def data_frame(self):
        return DataFrame([str(design) for design in self])

    def display_on_map(self, map_name):
        raise NotImplementedError
