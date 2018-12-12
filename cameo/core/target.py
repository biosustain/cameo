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

import numpy
from gnomic.genotype import Genotype
from gnomic.types import Accession, Feature, Change
from gnomic.utils import genotype_to_string, genotype_to_text

from cameo.core.manipulation import swap_cofactors, increase_flux, decrease_flux, reverse_flux
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
    id : str
        The identifier of the target. The id must be present in the COBRA model.
    accession_id : str (optional)
        An accession identifier in a database.
    accession_db : str (optional)
        The database corresponding to `accession_id`.
    """

    def __init__(self, id, accession_id=None, accession_db=None):
        self.id = id
        self.accession_id = accession_id
        self.accession_db = accession_db

    def apply(self, model):
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
        if self.accession_id:
            if self.accession_db:
                return Accession(identifier=self.accession_id, database=self.accession_db)
            else:
                return Accession(identifier=self.accession_id)
        else:
            return None

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

    def __init__(self, id, value, reference_value, *args, **kwargs):
        super(FluxModulationTarget, self).__init__(id, *args, **kwargs)
        # TODO: refactor targets to make the following work
        # if value == 0:
        #     raise ValueError("Please use ReactionKnockoutTarget if flux = 0.")
        self._value = value
        self._reference_value = reference_value

    def get_model_target(self, model):
        raise NotImplementedError

    def apply(self, model):
        """
        Applies a change to the flux. If the fold change is higher than 0 it increases the flux.
        If the target flux is zero it applies a knockout. If the fold change is lower than 0, then it
        decreases the flux.

        Parameters
        ----------
        model: cobra.Model
            A model.
        """
        target = self.get_model_target(model)

        if self._value == 0:
            target.knock_out()
        elif self.fold_change > 0:
            increase_flux(target, self._reference_value, self._value)
        elif self.fold_change < 0:
            decrease_flux(target, self._reference_value, self._value)

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
        return genotype_to_string(Genotype([self.to_gnomic()]))

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
        return genotype_to_text(Genotype([self.to_gnomic()]))

    def to_gnomic(self):
        accession = Target.to_gnomic(self)
        if self._value == 0:
            old_feature = Feature(name=self.id, accession=accession, type=self.__gnomic_feature_type__)
            new_feature = None
        elif self.fold_change != 0:
            old_feature = Feature(name=self.id, type=self.__gnomic_feature_type__,
                                  variant=tuple(["value={}".format(self._reference_value)]))
            new_feature = Feature(name=self.id, accession=accession, type=self.__gnomic_feature_type__,
                                  variant=tuple(["value={}".format(self._value)]))
        else:
            raise RuntimeError("fold_change shouldn't be 0")

        return Change(before=old_feature, after=new_feature)


class ReactionCofactorSwapTarget(Target):
    """
    Swap cofactors of a given reaction.
    """

    def __init__(self, id, swap_pairs, *args, **kwargs):
        super(ReactionCofactorSwapTarget, self).__init__(id, *args, **kwargs)
        self.swap_pairs = swap_pairs

    def apply(self, model):
        reaction = model.reactions.get_by_id(self.id)
        swap_cofactors(reaction, model, self.swap_pairs)

    @property
    def swap_str(self):
        return "+".join(str(s) for s in self.swap_pairs[0]) + "<->" + "+".join(str(s) for s in self.swap_pairs[1])

    def to_gnomic(self):
        accession = Target.to_gnomic(self)
        pairs = [(m.id for m in pair) for pair in self.swap_pairs]
        wt_feature = Feature(name=self.id, accession=accession, type='reaction',
                             variant=["cofactors=%s" % ",".join(pairs[0])])

        feature = Feature(name=self.id, accession=accession, type='reaction',
                          variant=["cofactors=%s" % ",".join(pairs[1])])
        return Change(before=wt_feature, after=feature)

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
    def __init__(self, id, value, *args, **kwargs):
        super(KnockinTarget, self).__init__(id, *args, **kwargs)
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
    def __init__(self, id, value, *args, **kwargs):
        super(ReactionKnockinTarget, self).__init__(id, value, *args, **kwargs)

    def apply(self, model):
        model.add_reactions([self._value])

    def to_gnomic(self):
        accession = Target.to_gnomic(self)
        feature = Feature(accession=accession, type='reaction', name=self.id)
        return Change(after=feature)

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
        return genotype_to_string(Genotype([self.to_gnomic()]))

    def __repr__(self):
        return "<ReactionKnockin %s>" % self.id

    def __hash__(self):
        return hash(str(self))

    def _repr_html_(self):
        return genotype_to_text(Genotype([self.to_gnomic()]))


class GeneModulationTarget(FluxModulationTarget):
    __gnomic_feature_type__ = "gene"

    def __init__(self, id, value, reference_value, *args, **kwargs):
        super(GeneModulationTarget, self).__init__(id, value, reference_value, *args, **kwargs)

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
            return (
                (self.id == other.id) and (
                    self._value == other._value) and (
                        self._reference_value == other._reference_value)
            )
        else:
            return False

    def __hash__(self):
        return hash(str(self))


class GeneKnockoutTarget(GeneModulationTarget):
    """
    Gene Knockout Target. Knockout a gene present in a COBRA model.
    """

    def __init__(self, id, *args, **kwargs):
        super(GeneKnockoutTarget, self).__init__(id, 0, None, *args, **kwargs)

    def apply(self, model):
        target = self.get_model_target(model)
        target.knock_out()

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

    def to_gnomic(self):
        accession = Target.to_gnomic(self)
        feature = Feature(name=self.id, accession=accession)
        return Change(before=feature)


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
            return (
                (self.id == other.id) and (
                    self._value == other._value) and (
                        self._reference_value == other._reference_value)
            )
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

    def apply(self, model):
        target = self.get_model_target(model)
        target.knock_out()

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

    def apply(self, model):
        reaction = self.get_model_target(model)
        reverse_flux(reaction, self._reference_value, self._value)

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

    def apply(self, model):
        for target in self.targets:
            target.apply(model)

    def __repr__(self):
        head = "<EnsembleTarget %s" % self.id
        body = ", ".join("%i - %s" % (i, str(target)) for i, target in enumerate(self.targets))
        end = ">"

        return head + "|" + body + "|" + end

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
