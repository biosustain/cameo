# Copyright 2016 The Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import sys
from functools import partial

import numpy
import six

from cameo.core.manipulation import swap_cofactors, increase_flux, decrease_flux, reverse_flux

try:
    from gnomic import Accession, Feature, Del, Mutation, Sub, Ins

    _gnomic_available_ = True
except (ImportError, SyntaxError):  # SyntaxError from py2 incompatible syntax
    _gnomic_available_ = False

from cameo import ui
from cameo.exceptions import IncompatibleTargets

__all__ = ["GeneModulationTarget", "GeneKnockoutTarget", "ReactionCofactorSwapTarget", "ReactionKnockinTarget",
           "ReactionKnockoutTarget", "ReactionModulationTarget"]


# TODO: This should be ABC.
class Target(object):
    """
    A Target is element of a COBRA model that will change to achieve a given phenotype.
    They are identified by StrainDesignMethod.run when identifying which manipulations
    are necessary to improve some fitness function.

    Attributes
    ----------
    id: str
        The identifier of the target. The id must be present in the COBRA model.

    """

    def __init__(self, id):
        self.id = id

    def apply(self, model, time_machine=None):
        """
        Applies the modification on the target, depending on the target type.

        See Also
        --------
        Subclass implementations

        """
        raise NotImplementedError

    def __eq__(self, other):
        if isinstance(other, Target):
            return other.id == self.id
        else:
            return False

    def __gt__(self, other):
        return self.id > other.id

    def to_gnomic(self):
        """
        If gnomic is available, return a Gnomic representation of the Target.

        """
        if _gnomic_available_:
            return Accession(identifier=self.id)
        else:
            raise SystemError("Gnomic is only compatible with python >= 3 (%i.%i)" %
                              (sys.version_info.major, sys.version_info.minor))

    def __repr__(self):
        return "<Target %s>" % self.id

    def __str__(self):
        return self.id

    def __hash__(self):
        return hash(str(self))


class FluxModulationTarget(Target):
    """
    A Target that changes the flux constraints (knockouts, over-expression or down-regulation).

    See Also
    --------
    ReactionModulationTarget, ReactionKnockoutTarget, GeneModulationTarget and GeneKnockoutTarget for implementation.

    """
    __gnomic_feature_type__ = 'flux'

    def __init__(self, id, value, reference_value):
        super(FluxModulationTarget, self).__init__(id)
        self._value = value
        self._reference_value = reference_value

    def get_model_target(self, model):
        raise NotImplementedError

    def apply(self, model, time_machine=None):
        """
        Applies a change to the flux. If the fold change is higher than 0 it increases the flux.
        If the target flux is zero it applies a knockout. If the fold change is lower than 0, then it
        decreases the flux.

        Parameters
        ----------
        model: cameo.core.solver_based_model.SolverBasedModel
            A model.
        time_machine: cameo.util.TimeMachine
            A TimeMachine for action stack.
        """
        target = self.get_model_target(model)

        if self._value == 0:
            target.knock_out(time_machine=time_machine)
        elif self.fold_change > 0:
            increase_flux(target, self._reference_value, self._value, time_machine=time_machine)
        elif self.fold_change < 0:
            decrease_flux(target, self._reference_value, self._value, time_machine=time_machine)

    def __eq__(self, other):
        if isinstance(other, FluxModulationTarget):
            return self.id == other.id and self._value == other._value
        else:
            return False

    @property
    def fold_change(self):
        """
        (B - A)/A

        Return
        ------
            float
        """
        try:
            ref = abs(self._reference_value)
            val = abs(self._value)
            return (val - ref) / ref
        except ZeroDivisionError:
            return numpy.inf

    def __str__(self):
        if self._value == 0:
            s = ui.delta() + self.id
        elif self.fold_change > 0:
            s = ui.upreg(self.fold_change) + self.id
        elif self.fold_change < 0:
            s = ui.downreg(self.fold_change) + self.id
        else:
            s = RuntimeError("fold_change shouldn't be 0")
        if six.PY2:
            return s.encode('utf-8')

        return s

    def __repr__(self):
        if self._value == 0:
            return "<FluxModulationTarget KO-%s>" % self.id
        elif self.fold_change > 0:
            return "<FluxModulationTarget UP(%.3f)-%s>" % (self.fold_change, self.id)
        elif self.fold_change < 0:
            return "<FluxModulationTarget DOWN(%.3f)-%s>" % (self.fold_change, self.id)
        else:
            raise RuntimeError("fold_change shouldn't be 0")

    def __hash__(self):
        return hash(str(self))

    def _repr_html_(self):
        if self._value == 0:
            return "&Delta;%s" % self.id
        elif self.fold_change > 0:
            return "&uarr;(%.3f)%s" % (self.fold_change, self.id)
        elif self.fold_change < 0:
            return "&darr;(%.3f)%s" % (self.fold_change, self.id)
        else:
            raise RuntimeError("fold_change shouldn't be 0")

    def to_gnomic(self):
        accession = Target.to_gnomic(self)
        feature = Feature(accession=accession, type=self.__gnomic_feature_type__)
        if self._value == 0:
            return Del(feature)
        elif self.fold_change > 0:
            over_expression = Feature(accession=accession, type=self.__gnomic_feature_type__,
                                      variant="over-expression(%.3f)" % self.fold_change)
            return Mutation(feature, over_expression)
        elif self.fold_change < 0:
            under_expression = Feature(accession=accession, type=self.__gnomic_feature_type__,
                                       variant="down-regulation(%.3f)" % self.fold_change)
            return Mutation(feature, under_expression)


class ReactionCofactorSwapTarget(Target):
    """
    Swap cofactors of a given reaction.
    """

    def __init__(self, id, swap_pairs):
        super(ReactionCofactorSwapTarget, self).__init__(id)
        self.swap_pairs = swap_pairs

    def apply(self, model, time_machine=None):
        reaction = model.reactions.get_by_id(self.id)
        swap_cofactors(reaction, model, self.swap_pairs, time_machine=time_machine)

    @property
    def swap_str(self):
        return "+".join(str(s) for s in self.swap_pairs[0]) + "<->" + "+".join(str(s) for s in self.swap_pairs[1])

    def to_gnomic(self):
        accession = Target.to_gnomic(self)
        new_accession = Accession(self.id + self.swap_str)
        original_feature = Feature(accession=accession, type='reaction')
        new_feature = Feature(accession=new_accession, type='reaction')
        return Sub(original_feature, new_feature)

    def __repr__(self):
        return "<ReactionCofactorSwap %s swap=%s>" % (self.id, self.swap_str)

    def __str__(self):
        return "%s[%s]" % (self.id, self.swap_str)

    def __gt__(self, other):
        if self.id == other.id:
            if isinstance(other, ReactionKnockoutTarget):
                raise IncompatibleTargets(self, other)
            if isinstance(other, ReactionModulationTarget):
                return False
            if isinstance(other, ReactionKnockinTarget):
                return True
            elif isinstance(other, EnsembleTarget):
                return not other > self
            else:
                raise IncompatibleTargets(self, other)
        else:
            return self.id > other.id

    def __eq__(self, other):
        if isinstance(other, ReactionCofactorSwapTarget):
            return self.id == other.id and self.swap_pairs == other.swap_pairs
        else:
            return False

    def _repr_html_(self):
        return self.id + "|" + self.swap_str.replace("<->", "&rlarr;")

    def __hash__(self):
        return hash(str(self))


class KnockinTarget(Target):
    def __init__(self, id, value):
        super(KnockinTarget, self).__init__(id)
        self._value = value

    def to_gnomic(self):
        raise NotImplementedError

    def __repr__(self):
        return "<KnockinTarget %s>" % self.id

    def __str__(self):
        return "::%s" % self.id

    def __hash__(self):
        return hash(str(self))


class ReactionKnockinTarget(KnockinTarget):
    def __init__(self, id, value):
        super(ReactionKnockinTarget, self).__init__(id, value)

    def apply(self, model, time_machine=None):
        # TODO: this is an hack. The objective_coefficient becomes 1 w/o explicitly being changed.
        self._value.objective_coefficient = 0
        if time_machine is None:
            model.add_reaction(self._value)
        else:
            time_machine(do=partial(model.add_reaction, self._value),
                         undo=partial(model.remove_reactions, [self._value], delete=False, remove_orphans=True))

    def to_gnomic(self):
        accession = Target.to_gnomic(self)
        feature = Feature(accession=accession, type='reaction')
        return Ins(feature)

    def __gt__(self, other):
        if self.id == other.id:
            if isinstance(other, ReactionKnockoutTarget):
                raise IncompatibleTargets(self, other)
            elif isinstance(other, ReactionModulationTarget):
                if other._value == 0:
                    raise IncompatibleTargets(self, other)
                else:
                    return False
            elif isinstance(other, ReactionKnockinTarget):
                return False
            elif isinstance(other, ReactionCofactorSwapTarget):
                return False
            elif isinstance(other, EnsembleTarget):
                return not other > self
            else:
                raise IncompatibleTargets(self, other)
        else:
            return self.id > other.id

    def __eq__(self, other):
        if isinstance(other, ReactionKnockinTarget):
            return self.id == other.id
        else:
            return False

    def __str__(self):
        return "::%s" % self.id

    def __repr__(self):
        return "<ReactionKnockin %s>" % self.id

    def __hash__(self):
        return hash(str(self))

    def _repr_html_(self):
        return "::%s" % self.id


class GeneModulationTarget(FluxModulationTarget):
    __gnomic_feature_type__ = "gene"

    def __init__(self, id, value, reference_value):
        super(GeneModulationTarget, self).__init__(id, value, reference_value)

    def get_model_target(self, model):
        return model.genes.get_by_id(self.id)

    def __gt__(self, other):
        if self.id == other.id:
            if isinstance(other, GeneKnockoutTarget):
                return False
            elif isinstance(other, GeneModulationTarget) and not isinstance(other, GeneKnockoutTarget):
                return self.fold_change > other.fold_change
            elif isinstance(other, EnsembleTarget):
                return not other > self
            else:
                raise IncompatibleTargets(self, other)
        else:
            return self.id > other.id

    def __eq__(self, other):
        if isinstance(other, GeneModulationTarget):
            return (self.id == other.id and self._value == other._value and
                    self._reference_value == other._reference_value)
        else:
            return False

    def __hash__(self):
        return hash(str(self))


class GeneKnockoutTarget(GeneModulationTarget):
    """
    Gene Knockout Target. Knockout a gene present in a COBRA model.
    """

    def __init__(self, id):
        super(GeneKnockoutTarget, self).__init__(id, 0, None)

    def apply(self, model, time_machine=None):
        target = self.get_model_target(model)
        target.knock_out(time_machine=time_machine)

    def __gt__(self, other):
        if self.id == other.id:
            if isinstance(other, GeneModulationTarget) and not isinstance(other, GeneKnockoutTarget):
                return self.fold_change > other.fold_change
            elif isinstance(other, GeneKnockoutTarget):
                return False
            elif isinstance(other, EnsembleTarget):
                return not other > self
            else:
                raise IncompatibleTargets(self, other)
        else:
            return self.id > other.id

    def __eq__(self, other):
        if isinstance(other, GeneKnockoutTarget):
            return self.id == other.id
        elif isinstance(other, GeneModulationTarget):
            return self.id == other.id and other._value == 0
        else:
            return False

    def __repr__(self):
        return "<GeneKnockout %s>" % self.id

    def __hash__(self):
        return hash(str(self))


class ReactionModulationTarget(FluxModulationTarget):
    __gnomic_feature_type__ = "reaction"

    def __init__(self, id, value, reference_value):
        super(ReactionModulationTarget, self).__init__(id, value, reference_value)

    def get_model_target(self, model):
        return model.reactions.get_by_id(self.id)

    def __gt__(self, other):
        if self.id == other.id:
            if isinstance(other, ReactionKnockinTarget):
                if self._value == 0:
                    raise IncompatibleTargets(self, other)
                else:
                    return True
            elif isinstance(other, ReactionCofactorSwapTarget):
                return True
            elif isinstance(other, ReactionModulationTarget) and not isinstance(other, ReactionKnockoutTarget):
                return self.fold_change > other.fold_change
            elif isinstance(other, EnsembleTarget):
                return not other > self
            else:
                raise IncompatibleTargets(self, other)
        else:
            return self.id > other.id

    def __eq__(self, other):
        if isinstance(other, ReactionModulationTarget):
            return (self.id == other.id and self._value == other._value and
                    self._reference_value == other._reference_value)
        else:
            return False

    def __hash__(self):
        return hash(str(self))


class ReactionKnockoutTarget(ReactionModulationTarget):
    """
    Reaction Knockout Target. Knockout a reaction present in a COBRA model.
    """

    def __init__(self, id):
        super(ReactionKnockoutTarget, self).__init__(id, 0, None)

    def apply(self, model, time_machine=None):
        target = self.get_model_target(model)
        target.knock_out(time_machine=time_machine)

    def __gt__(self, other):
        if self.id == other.id:
            if isinstance(other, ReactionModulationTarget):
                if other._value == 0:
                    return True
                else:
                    raise IncompatibleTargets(self, other)
            elif isinstance(other, ReactionCofactorSwapTarget):
                raise IncompatibleTargets(self, other)
            elif isinstance(other, EnsembleTarget):
                return not other > self
            else:
                raise IncompatibleTargets(self, other)
        else:
            return self.id > other.id

    def __eq__(self, other):
        if isinstance(other, ReactionKnockoutTarget):
            return self.id == other.id
        elif isinstance(other, ReactionModulationTarget):
            return self.id == other.id and other._value == 0
        else:
            return False

    def __repr__(self):
        return "<ReactionKnockout %s>" % self.id

    def __hash__(self):
        return hash(str(self))


class ReactionInversionTarget(ReactionModulationTarget):
    def __str__(self):
        return "INV(%.3f -> %.3f)-%s" % (self._reference_value, self._value, self.id)

    def __repr__(self):
        return "<ReactionInversion %s (%.5f -> %.5f)>" % (self.id, self._reference_value, self._value)

    def _repr_html_(self):
        return "&#8645;-%s" % self.id

    def __gt__(self, other):
        if self.id == other.id:
            if isinstance(other, ReactionModulationTarget):
                if other._value == self._value:
                    return True
                else:
                    raise IncompatibleTargets(self, other)
            elif isinstance(other, ReactionKnockoutTarget):
                raise IncompatibleTargets(self, other)
            elif isinstance(other, ReactionKnockinTarget):
                return True
            else:
                raise IncompatibleTargets(self, other)
        else:
            return self.id > other.id

    def __eq__(self, other):
        if isinstance(other, ReactionInversionTarget):
            return self.id == other.id and other._value == self._value
        else:
            return False

    def apply(self, model, time_machine=None):
        reaction = self.get_model_target(model)
        reverse_flux(reaction, self._reference_value, self._value, time_machine=time_machine)

    def __hash__(self):
        return hash(str(self))


class EnsembleTarget(Target):
    """
    A Target resulting of merging multiple Targets.
    For example:
        1. ReactionKnockinTarget
        2. ReactionModulationTarget

    The Targets are prioritized by what should happen first in order to don't break.

    """

    def __init__(self, id, targets):
        super(EnsembleTarget, self).__init__(id)
        assert all(t.id == id for t in targets)
        self.targets = list(sorted(targets))

    def apply(self, model, time_machine=None):
        for target in self.targets:
            target.apply(model, time_machine=time_machine)

    def __repr__(self):
        head = "<EnsembleTarget %s" % self.id
        body = ", ".join("%i - %s" % (i, str(target)) for i, target in enumerate(self.targets))
        end = ">"

        return head + "|" + body + "|" + end

    # TODO implement gnomic compatibility
    def to_gnomic(self):
        raise NotImplementedError

    def _repr_html_(self):
        return "|" + ";".join(target._repr_html_() for target in self.targets) + '|'

    def __gt__(self, other):
        if isinstance(other, EnsembleTarget):
            return self.id > other.id

        else:
            is_greater = False
            for t in self.targets:
                is_greater = is_greater or t >= other

    def __eq__(self, other):
        if isinstance(other, EnsembleTarget):
            if len(self.targets) == len(other.targets):
                return all(target == other.targets[i] for i, target in enumerate(self.targets))

        return False

    def __str__(self):
        return "%s[%s]" % (self.id, ", ".join(str(t) for t in self.targets))

    def __hash__(self):
        return hash(str(self))
