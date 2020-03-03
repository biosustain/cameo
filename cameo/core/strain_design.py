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
"""
Core implementation of strain design. It contains core structures.
Targets (and subclasses) are identified by strain design methods.
"""

import sys

from cobra import DictList
from pandas import DataFrame

from cameo.core.result import Result
from cameo.core.target import EnsembleTarget, Target
from cameo.visualization.plotting import plotter

from gnomic import Genotype

__all__ = ['StrainDesign', 'StrainDesignMethod', 'StrainDesignMethodResult']


class StrainDesignMethod(object):
    def __init__(self, *args, **kwargs):
        super(StrainDesignMethod, self).__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):
        self.run(*args, **kwargs)

    def run(self, *args, **kwargs):
        raise NotImplementedError


class StrainDesign(object):
    """
    A StrainDesign is a collection of targets in a COBRA model. The targets, identified by a StrainDesignMethod,
    map elements in the model that need to be modified to achieve the objective of the design method.
    """
    def __init__(self, targets):
        self.targets = DictList(targets)

    def __str__(self):
        return ", ".join(str(t) for t in self.targets)

    def __repr__(self):
        return "<StrainDesign [" + ";".join(repr(t) for t in self.targets) + "]>"

    def __iter__(self):
        return iter(self.targets)

    def __len__(self):
        return len(self.targets)

    def __contains__(self, item):
        if isinstance(item, Target):
            if item in self.targets:
                return item == self.targets.get_by_id(item.id)
            else:
                return False
        elif isinstance(item, str):
            return item in self.targets
        else:
            return False

    def __eq__(self, other):
        if isinstance(other, StrainDesign):
            if len(self) != len(other):
                return False
            else:
                return all(t in other for t in self.targets) and all(t in self for t in other.targets)
        else:
            return False

    def apply(self, model):
        for target in self.targets:
            target.apply(model)

    def __add__(self, other):
        if not isinstance(other, StrainDesign):
            raise AssertionError("Only instances of StrainDesign can be added together")

        targets = {}
        for target in self.targets:
            if target.id not in targets:
                targets[target.id] = set()
            targets[target.id].add(target)

        for target in other.targets:
            if target.id not in targets:
                targets[target.id] = set()
            targets[target.id].add(target)

        targets = [next(iter(t)) if len(t) == 1 else EnsembleTarget(id, t) for id, t in targets.items()]

        return StrainDesign(targets)

    def __iadd__(self, other):
        if not isinstance(other, StrainDesign):
            raise AssertionError("Only instances of StrainDesign can be added together")

        targets = {}
        for target in self.targets:
            if target.id not in targets:
                targets[target.id] = set()
            targets[target.id].add(target)

        for target in other.targets:
            if target.id not in targets:
                targets[target.id] = set()
            targets[target.id].add(target)

        targets = [next(iter(t)) if len(t) == 1 else EnsembleTarget(id, t) for id, t in targets.items()]

        self.targets = DictList(targets)

        return self

    def _repr_html_(self):
        return " ".join(t._repr_html_() for t in self.targets)

    def to_gnomic(self):
        return Genotype([target.to_gnomic() for target in self.targets])


class StrainDesignMethodResult(Result):
    __method_name__ = None

    _aggreate_functions_ = {
        "method": lambda series: tuple(sum(series.values.tolist(), []))
    }

    def __init__(self, designs, *args, **kwargs):
        super(StrainDesignMethodResult, self).__init__(*args, **kwargs)
        self._designs = designs

    def __len__(self):
        """
        Return the number of designs found.
        """
        return len(self._designs)

    def __iter__(self):
        """
        Returns an iterator that yields StrainDesign objects.
        """
        return iter(self._designs)

    @property
    def data_frame(self):
        df = DataFrame(columns=['targets'])
        for i, design in enumerate(self._designs):
            df.loc[i] = [design.targets]
        return df

    def display_on_map(self, index, map_name):
        raise NotImplementedError

    def __add__(self, other):
        df = DataFrame([design.targets for design in self._designs], columns=["targets"])
        df['method'] = self.__method_name__

        i = len(df)
        for j, design in enumerate(other):
            df.loc[i + j] = [design, [self.__method_name__]]

        df = df.groupby(["design"]).aggregate(self._aggreate_functions_)

        return StrainDesignMethodEnsemble([StrainDesign(targets) for targets in df.targets], df.method.tolist())

    def _repr_html_(self):
        return self.data_frame._repr_html_()

    def plot(self, grid=None, width=None, height=None, title=None, *args, **kwargs):
        if title is None:
            title = "Target frequency plot for %s result" % self.__method_name__

        counts = DataFrame(columns=['count'])
        for design in self._designs:
            for target in design:
                if target.id not in counts.index:
                    counts.loc[target.id, 'count'] = 0
                counts.loc[target.id, 'count'] += 1
            counts['frequency'] = counts['count'].apply(lambda c: c / len(self._designs))

        plotter.frequency(counts, title=title, width=width, height=height, **kwargs)


class StrainDesignMethodEnsemble(StrainDesignMethodResult):
    def __init__(self, designs, methods, *args, **kwargs):
        super(StrainDesignMethodEnsemble, self).__init__(designs, *args, **kwargs)
        self._methods = methods

    @property
    def data_frame(self):
        data = []
        for design, method in zip(self._designs, self._methods):
            data.append([design.targets, method])

        return DataFrame(data, columns=["design", "method"])

    def __add__(self, other):
        df = self.data_frame

        i = len(df)
        for j, design in enumerate(other):
            if isinstance(other, StrainDesignMethodEnsemble):
                df.loc[i + j] = [design, self._methods]
            else:
                df.loc[i + j] = [design, [self.__method_name__]]

        df = df.groupby("design").aggregate(self._aggreate_functions_)

        return StrainDesignMethodEnsemble([StrainDesign(targets) for targets in df.targets], df['method'].tolist())
