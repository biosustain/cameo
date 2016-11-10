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
from functools import partial

from cobra import DictList
from pandas import DataFrame

from cameo import ui
from cameo.core.result import Result


class Target(object):
    def __init__(self, id):
        self.id = id

    def apply(self, model, time_machine=None):
        raise NotImplemented

    def __eq__(self, other):
        if isinstance(other, Target):
            return other.id == self.id
        else:
            return False


class FluxModulationTarget(Target):
    def __init__(self, id, value, reference_value):
        super(FluxModulationTarget, self).__init__(id)
        self._reference_value = reference_value
        self._value = value

    def get_target(self, model):
        raise NotImplemented

    def apply(self, model, time_machine=None):
        target = self.get_target(model)

        if self._value > 0:
            target.overexpress(self._value, reference_value=self._reference_value, time_machine=time_machine)
        elif self._value < 0:
            target.downregulate(self._value, reference_value=self._reference_value, time_machine=time_machine)
        else:
            target.knock_out(time_machine=time_machine)

    def __eq__(self, other):
        if isinstance(other, FluxModulationTarget):
            return self.id == other.id and self._value == other._value
        else:
            return False

    def __str__(self):
        fold_change = (self._value - self._reference_value) / self._reference_value

        if self._value == 0:
            return ui.delta() + self.id
        elif self._value > self._reference_value:
            return ui.upreg(fold_change) + self.id
        elif self._value < self._reference_value:
            return ui.downreg(fold_change) + self.id


class SwapTarget(Target):
    def __init__(self, id):
        super(SwapTarget, self).__init__(id)


class ReactionCofactorSwapTarget(Target):
    def __init__(self, id, swap_pairs):
        super(ReactionCofactorSwapTarget, self).__init__(id)
        self.swap_pairs = swap_pairs

    def apply(self, model, time_machine=None):
        reaction = model.reactions.get_by_id(self.id)
        reaction.swap_cofactors(self.swap_pairs, time_machine=time_machine)


class KnockinTarget(Target):
    def __init__(self, id, value):
        super(KnockinTarget, self).__init__(id)
        self._value = value


class ReactionKnockinTarget(KnockinTarget):
    def __init__(self, id, value):
        super(ReactionKnockinTarget, self).__init__(id, value)

    def apply(self, model, time_machine=None):
        if time_machine is None:
            model.add_reaction(self._value)
        else:
            time_machine(do=partial(model.add_reaction, self._value),
                         undo=partial(model.remove_reactions, self._value, delete=False, remove_orphans=True))


class GeneModulationTarget(FluxModulationTarget):
    def __init__(self, id, value, reference_value):
        super(GeneModulationTarget, self).__init__(id, value, reference_value)

    def get_target(self, model):
        return model.genes.get_by_id(self.id)


class GeneKnockoutTarget(GeneModulationTarget):
    def __init__(self, id):
        super(GeneKnockoutTarget, self).__init__(id, 0, None)

    def apply(self, model, time_machine=None):
        target = self.get_target(model)
        target.knock_out(time_machine=time_machine)


class ReactionModulationTarget(FluxModulationTarget):
    def __init__(self, id, value, reference_value):
        super(ReactionModulationTarget, self).__init__(id, value, reference_value)

    def get_target(self, model):
        return model.reactions.get_by_id(self.id)


class ReactionKnockoutTarget(ReactionModulationTarget):
    def __init__(self, id):
        super(ReactionKnockoutTarget, self).__init__(id, 0, None)

    def apply(self, model, time_machine=None):
        target = self.get_target(model)
        target.knock_out(time_machine=time_machine)


class StrainDesignMethod(object):
    def __init__(self, *args, **kwargs):
        super(StrainDesignMethod, self).__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):
        self.run(*args, **kwargs)

    def run(self, *args, **kwargs):
        raise NotImplementedError


class StrainDesign(object):
    def __init__(self, targets):
        self.targets = DictList(targets)

    def __repr__(self):
        return "".join(str(t) for t in self.targets)

    def __iter__(self):
        for t in self.targets:
            yield t

    def __len__(self):
        return len(self.targets)

    def __eq__(self, other):
        if isinstance(other, StrainDesign):
            if len(self) != len(other):
                return False
            else:
                return all(self.targets[i] == other.targets[i] for i in range(self))
        else:
            return False

    def apply(self, model, time_machine=None):
        for target in self.targets:
            target.apply(model, time_machine)

    def __add__(self, other):
        return StrainDesign(self.targets + other.targets)


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
        for design in self._designs:
            yield design

    @property
    def data_frame(self):
        return DataFrame([design.targets for design in self._designs], columns=["targets"])

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
