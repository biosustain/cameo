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

from __future__ import absolute_import, print_function

try:
    import cPickle as pickle
except ImportError:
    import pickle

import os
import six
import optlang
from cobra.io import read_sbml_model, load_json_model

from cameo.core.solver_based_model import SolverBasedModel, to_solver_based_model
from cameo.config import solvers

import logging

logger = logging.getLogger(__name__)


def load_model(path_or_handle, solver_interface=optlang, sanitize=True):
    """Read a metabolic model .

    Parameters
    ----------
    path_or_handle : path, fhandle or name.
        One of:
            * file path of a model file;
            * file handle to a SBML or pickled model; or
            * the identifier of a model in a web database (optflux.org/models)
    solver_interface : solver_interface, optional
        E.g. optlang.glpk_interface or any other optlang interface.
    sanitize : boolean, optional
        If reaction and metabolite IDs should be sanitized (works only for SBML models).
    """
    solver_interface = solvers.get(solver_interface, solver_interface)

    if isinstance(path_or_handle, six.string_types) and not os.path.isfile(path_or_handle):
        from cameo.models.webmodels import load_webmodel
        logger.debug("Given path is not a file. Trying to load from webmodels")
        model = load_webmodel(path_or_handle, solver_interface)
    else:
        if isinstance(path_or_handle, six.string_types):
            # Open the given file
            path = path_or_handle
            handle = open(path_or_handle, 'rb')
        elif hasattr(path_or_handle, 'read'):
            # Argument is already an open file
            path = path_or_handle.name
            handle = path_or_handle
        else:
            raise ValueError('Provided argument %s has to be either a string or a file handle' % path_or_handle)
        model = _load_model_from_file(path, handle)  # Parse model from the file

    if sanitize:
        sanitize_ids(model)

    if not isinstance(model, SolverBasedModel):
        if solver_interface is not None:
            logger.debug("Changing solver interface to %s" % solver_interface)
            model = to_solver_based_model(model, solver_interface=solver_interface)
    else:
        if solver_interface is not None and not isinstance(model.solver, solver_interface.Model):
            logger.debug("Changing solver interface to %s" % solver_interface)
            model.solver = solver_interface

    return model


def _load_model_from_file(path, handle):
    """Try to parse a model from a file handle using different encodings."""
    logger.debug('Reading file from %s assuming pickled model.' % path)
    try:
        model = pickle.load(handle)
    except (TypeError, pickle.UnpicklingError):
        logger.debug('Cannot unpickle %s. Assuming json model next.' % path)
        try:
            model = load_json_model(path)
        except ValueError:
            logger.debug("Cannot import %s as json model. Assuming sbml model next." % path)
            try:
                model = read_sbml_model(path)
            except AttributeError as e:
                logger.error("cobrapy doesn't raise a proper exception if a file does not contain an SBML model")
                raise e
            except Exception as e:
                logger.error(
                    "Looks like something blow up while trying to import {} as a SBML model."
                    "Try validating the model at http://sbml.org/Facilities/Validator/ to get more information.".format(
                        path))
                raise e
    return model


ID_SANITIZE_RULES_SIMPHENY = [('_DASH_', '-'), ('_FSLASH_', '/'), ('_BSLASH_', "\\"), ('_LPAREN_', '('),
                              ('_LSQBKT_', '['),
                              ('_RSQBKT_', ']'), ('_RPAREN_', ')'), ('_COMMA_', ','), ('_PERIOD_', '.'),
                              ('_APOS_', "'"),
                              ('&amp;', '&'), ('&lt;', '<'), ('&gt;', '>'), ('&quot;', '"')]

ID_SANITIZE_RULES_TAB_COMPLETION = [('_DASH_', '_dsh_'), ('_FSLASH_', '_fsh_'), ('_BSLASH_', "_bsh_"),
                                    ('_LPAREN_', '_lp_'), ('_LSQBKT_', '_lb_'),
                                    ('_RSQBKT_', '_rb_'), ('_RPAREN_', '_rp_'), ('_COMMA_', '_cm_'),
                                    ('_PERIOD_', '_prd_'), ('_APOS_', "_apo_"),
                                    ('&amp;', '_amp_'), ('&lt;', '_lt_'), ('&gt;', '_gt_'), ('&quot;', '_qot_')]


def _apply_sanitize_rules(id, rules):
    for rule in rules:
        id = id.replace(*rule)
    return id


def sanitize_ids(model):
    """Makes IDs crippled by the XML specification less annoying.

    For example, EX_glc_LPAREN_e_RPAREN_ will be converted to EX_glc_lp_e_rp_.
    Furthermore, reactions and metabolites will be equipped with a nice_id attribute
    that provides the original ID, i.e., EX_glc(d).

    Parameters
    ----------
    model : model

    Notes
    -----

    Will add a nice_id attribute.

    """

    for metabolite in model.metabolites:
        met_id = metabolite.id
        metabolite.id = _apply_sanitize_rules(met_id, ID_SANITIZE_RULES_TAB_COMPLETION)
        metabolite.nice_id = _apply_sanitize_rules(met_id, ID_SANITIZE_RULES_SIMPHENY)
        metabolite.name = _apply_sanitize_rules(metabolite.name, ID_SANITIZE_RULES_SIMPHENY)

    for reaction in model.reactions:
        if isinstance(model, SolverBasedModel):
            forward_variable = reaction.forward_variable
            reverse_variable = reaction.reverse_variable
        rxn_id = reaction.id
        reaction.id = _apply_sanitize_rules(rxn_id, ID_SANITIZE_RULES_TAB_COMPLETION)
        reaction.nice_id = _apply_sanitize_rules(rxn_id, ID_SANITIZE_RULES_SIMPHENY)
        reaction.name = _apply_sanitize_rules(reaction.name, ID_SANITIZE_RULES_SIMPHENY)
        if isinstance(model, SolverBasedModel):
            forward_variable.name = reaction.id
            if reverse_variable is not None:
                reverse_variable.name = reaction._get_reverse_id()

    if isinstance(model, SolverBasedModel):
        for constraint in model.solver.constraints:
            constraint.name = _apply_sanitize_rules(constraint.name, ID_SANITIZE_RULES_TAB_COMPLETION)

    model.repair()
