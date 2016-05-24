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
from cameo.util import frozendict


class StrainDesignMethod(object):
    def __init__(self, *args, **kwargs):
        super(StrainDesignMethod, self).__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):
        self.run(*args, **kwargs)

    def run(self, *args, **kwargs):
        raise NotImplementedError


class StrainDesign(object):
    def __init__(self, knockouts=None, knock_ins=None, over_expression=None, down_regulation=None,
                 manipulation_type="genes"):
        self.knockouts = knockouts or []
        self.knock_ins = knock_ins or []
        self.over_expression = frozendict(over_expression or {})
        self.down_regulation = frozendict(down_regulation or {})
        self.manipulation_type = manipulation_type

    def __repr__(self):
        return "".join([ui.delta() + ko for ko in self.knockouts] +
                       [ui.knockin() + ki for ki in self.knock_ins] +
                       [ui.upreg(coeff) + oe for oe, coeff in six.iteritems(self.over_expression)] +
                       [ui.downreg(coeff) + dr for dr, coeff in six.iteritems(self.down_regulation)])

    def __iter__(self):
        yield tuple(set(self.knockouts))
        yield tuple(set(self.knock_ins))
        yield self.over_expression
        yield self.down_regulation


class StrainDesignResult(Result):
    __method_name__ = None

    _aggreate_functions_ = {
        "method": lambda series: tuple(sum(series.values.tolist(), []))
    }

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

    @property
    def data_frame(self):
        return DataFrame([design for design in self],
                         columns=["knockouts", "knock_ins", "over_expression", "down_regulation"])

    def display_on_map(self, map_name):
        raise NotImplementedError

    def __add__(self, other):
        df = DataFrame(columns=["knockouts", "knock_ins", "over_expression", "down_regulation", "type", "method"])
        i = 0
        for i, design in enumerate(self):
            df.loc[i] = list(design) + [design.manipulation_type, [self.__method_name__]]

        for j, design in enumerate(other):
            df.loc[i + j] = list(design) + [design.manipulation_type, [self.__method_name__]]

        df = df.groupby(["knockouts", "knock_ins",
                         "over_expression", "down_regulation", "type"]).aggregate(self._aggreate_functions_)

        return StrainDesignEnsemble(df.index.tolist(), df['method'].tolist())

    def _repr_html_(self):
        return self.data_frame._repr_html_()


class StrainDesignEnsemble(StrainDesignResult):
    def __init__(self, designs, methods, *args, **kwargs):
        super(StrainDesignEnsemble, self).__init__(*args, **kwargs)
        self._designs = designs
        self._methods = methods

    def __iter__(self):
        for design in self._designs:
            yield design

    def __len__(self):
        return len(self._designs)

    def data_frame(self):
        data = []
        for design, method in zip(self._designs, self._methods):
            row = list(design)
            row.append(method)

        return DataFrame(data, columns=["knockouts", "knock_ins", "over_expression", "down_regulation", "method"])

    def __add__(self, other):
        df = DataFrame(columns=["knockouts", "knock_ins", "over_expression", "down_regulation", "type", "method"])
        i = 0
        for i, design in enumerate(self):
            df.loc[i] = list(design) + [design.manipulation_type, self._methods]

        for j, design in enumerate(other):
            if isinstance(other, StrainDesignEnsemble):
                df.loc[i + j] = list(design) + [design.manipulation_type, self._methods]
            else:
                df.loc[i + j] = list(design) + [design.manipulation_type, [self.__method_name__]]

        df = df.groupby(["knockouts", "knock_ins",
                         "over_expression", "down_regulation", "type"]).aggregate(self._aggreate_functions_)

        designs = [StrainDesign(row.values[:-1]) for _, row in df.iterrows()]

        return StrainDesignEnsemble(designs, df['method'].tolist())
