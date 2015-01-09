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

import logging
import types
import pickle
from multiprocessing import Process, Queue
import optlang
from cobra.io import read_sbml_model
import time
import progressbar
from cameo.solver_based_model import SolverBasedModel, to_solver_based_model
from cameo import webmodels

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def load_model(path_or_handle, solver_interface=optlang.glpk_interface, sanitize=True):
    """Read a model from a file.

    Parameters
    ----------
    path_or_handle : path, fhandle
        Either a file path or a file handle to a SBML or pickled model
    solver_interface : solver_interface, optional
        e.g. optlang.glpk_interface or any other optlang interface
    sanitize : boolean, optional
        if reaction and metabolite IDs should be sanitized (works only for SBML models)
    """

    if isinstance(path_or_handle, types.StringType):
        path = path_or_handle
        try:
            handle = open(path_or_handle)
        except IOError:
            logger.debug('%s not a file path. Querying webmodels ...' % path)
            try:
                df = webmodels.index_models()
                index = df.query('name == "%s"' % path_or_handle).id.values[0]
                handle = webmodels.get_sbml_file(index)
                path = handle.name
            except IndexError:
                raise ValueError("%s is neither a file nor a model ID." % path)
    elif hasattr(path_or_handle, 'read'):
        path = path_or_handle.name
        handle = path_or_handle
    else:
        raise ValueError('Provided argument %s has to be either a file path or handle' % path_or_handle)
    logger.debug('Reading file from %s assuming pickled model.' % path)
    try:
        model = pickle.load(handle)
    except Exception:
        logger.debug('Cannot unpickle %s. Assuming sbml model next.' % path)
        try:
            interactive = False  # TODO: make this figure out if run inside a shell or IPython notebook
            if interactive:
                def read_sbml_model_subprocess(path, queue):
                    queue.put(read_sbml_model(path))
                pbar = progressbar.ProgressBar()
                queue = Queue()
                proc = Process(target=read_sbml_model_subprocess, args=(path, queue))
                proc.start()
                while proc.is_alive():
                    time.sleep(.5)
                model = queue.get()
            else:
                model = read_sbml_model(path)
            if sanitize:
                sanitize_ids(model)
        except AttributeError:  # TODO: cobrapy doesn't raise a proper exception if a file does not contain an SBML model
            raise ValueError('FIXME')

    if not isinstance(model, SolverBasedModel):
        if solver_interface is not None:
            logger.debug("Changing solver interface to %s" % solver_interface)
            model = to_solver_based_model(model, solver_interface=solver_interface)
    else:
        if model.interface is not solver_interface and solver_interface is not None:
            logger.debug("Changing solver interface to %s" % solver_interface)
            model.solver = solver_interface

    return model


ID_SANITIZE_RULES_SIMPHENY = [('_DASH_', '-'), ('_FSLASH_', '/'),('_BSLASH_', "\\"), ('_LPAREN_', '('), ('_LSQBKT_', '['),
                     ('_RSQBKT_', ']'), ('_RPAREN_', ')'), ('_COMMA_', ','), ('_PERIOD_', '.'), ('_APOS_', "'"),
                     ('&amp;', '&'), ('&lt;', '<'), ('&gt;', '>'), ('&quot;', '"')]

ID_SANITIZE_RULES_TAB_COMPLETION = [('_DASH_', '_dsh_'), ('_FSLASH_', '_fsh_'),('_BSLASH_', "_bsh_"), ('_LPAREN_', '_lp_'), ('_LSQBKT_', '_lb_'),
                     ('_RSQBKT_', '_rb_'), ('_RPAREN_', '_rp_'), ('_COMMA_', '_cm_'), ('_PERIOD_', '_prd_'), ('_APOS_', "_apo_"),
                     ('&amp;', '_amp_'), ('&lt;', '_lt_'), ('&gt;', '_gt_'), ('&quot;', '_qot_')]

def _apply_sanitize_rules(id, rules):
        for rule in rules:
            id = id.replace(*rule)
        return id

def sanitize_ids(model):
    """Makes IDs crippled by the XML specification less annoying.

    For example, 'EX_glc_LPAREN_e_RPAREN_' will be converted to 'EX_glc_lp_e_rp_'.
    Furthermore, reactions and metabolites will be equipped with a `nice_id` attribute
    that provides the original ID, i.e., 'EX_glc(d)'.

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
            variable = reaction.variable
            reverse_variable = reaction.reverse_variable
        rxn_id = reaction.id
        reaction.id = _apply_sanitize_rules(rxn_id, ID_SANITIZE_RULES_TAB_COMPLETION)
        reaction.nice_id = _apply_sanitize_rules(rxn_id, ID_SANITIZE_RULES_SIMPHENY)
        reaction.name = _apply_sanitize_rules(reaction.name, ID_SANITIZE_RULES_SIMPHENY)
        if isinstance(model, SolverBasedModel):
            variable.name = reaction.id
            if reverse_variable is not None:
                reverse_variable.name = reaction._get_reverse_id()

    if isinstance(model, SolverBasedModel):
        for constraint in model.solver.constraints:
            constraint.name = _apply_sanitize_rules(constraint.name, ID_SANITIZE_RULES_TAB_COMPLETION)

    model.repair()