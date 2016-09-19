import re
from cameo.data import metanetx
from cameo.core.reaction import Reaction
from cameo.core.metabolite import Metabolite
import logging
logger = logging.getLogger(__name__)


def clean_bigg_id(string):
    return re.sub(r"bigg:|dsh", "", string)


def get_existing_metabolite(mnx_id, model, compartment):
    """Find compartment in the model by Metanetx id.

    Parameters
    ----------
    mnx_id : string
        Metanetx id
    model
        cobra model
    compartment : string
        f.e "_c"

    Returns
    -------
    Metabolite or None

    """
    if not mnx_id:
        return
    try:
        clean_id = clean_bigg_id(metanetx.mnx2bigg[mnx_id])
        return model.metabolites.get_by_id(clean_id + compartment)
    except KeyError:
        try:
            return model.metabolites.get_by_id(mnx_id + compartment)
        except KeyError:
            pass


def contains_carbon(metabolite):  # TODO: use method from Metabolite class when this change is merged
    if not metabolite.formula:
        raise ValueError("No formula for metabolite {}, it's unknown if there is carbon in it")
    return 'C' in metabolite.elements


def find_metabolite_info(met_id):
    """Find chemical formula of metabolite in metanetx.chem_prop dictionary

    Parameters
    ----------
    met_id : string
        string of format "<metabolite_id>_<compartment_id>", where <metabolite_id> is a BIGG id or a Metanetx id

    Returns
    -------
    pandas row or None

    """
    met_id = met_id[:-2]
    try:
        if met_id in metanetx.chem_prop.index:
            return metanetx.chem_prop.loc[met_id]
        return metanetx.chem_prop.loc[metanetx.all2mnx['bigg:' + met_id]]
    except KeyError:
        return None


class ModelModification(object):
    """
    Base model modification class, providing methods for adding adapter and exchange reactions for new metabolites
    """
    model = None
    added_reactions = None

    def create_exchange(self, gene_name, metabolite):
        """For given metabolite A_c from c compartment, create:
        a) corresponding metabolite A_e from e compartment;
        b) adapter reaction A_c <--> A_e
        c) exchange reaction A_e -->

        Parameters
        ----------
        gene_name : string
            metabolite id in format <bigg_id>_c, f.e. Nacsertn_c
        metabolite : Metabolite
            gene associated with the metabolite

        Returns
        -------

        """
        exchange_metabolite = Metabolite(metabolite.id.replace('_c', '_e'), formula=metabolite.formula, compartment='e')
        self.add_adapter_reaction(metabolite, exchange_metabolite, gene_name)
        self.add_exchange_reaction(exchange_metabolite)

    def add_exchange_reaction(self, metabolite):
        """Add exchange reaction "A --> " for given metabolite A
         If reaction exists, log and pass

        Parameters
        ----------
        metabolite : basestring
            metabolite id in format <bigg_id>_<compartment_id>, f.e. Nacsertn_c

        Returns
        -------

        """
        try:
            logger.debug('Add exchange reaction for metabolite: {}'.format(metabolite.id))
            exchange_reaction = self.model.add_exchange(metabolite, prefix='EX_')
            self.added_reactions.add(exchange_reaction.id)
        except ValueError:
            logger.debug('Exchange reaction exists for metabolite {}'.format(metabolite.id))

    def add_adapter_reaction(self, metabolite, existing_metabolite, gene_name):
        """Add adapter reaction A <--> B for metabolites A and B

        Parameters
        ----------
        metabolite : Metabolite
            metabolite A
        existing_metabolite : Metabolite
            metabolite B
        gene_name : basestring
            gene associated with the metabolites

        Returns
        -------

        """
        try:
            adapter_reaction = Reaction(str('adapter_' + metabolite.id + '_' + existing_metabolite.id))
            adapter_reaction.lower_bound = -1000
            adapter_reaction.add_metabolites({metabolite: -1, existing_metabolite: 1})
            adapter_reaction.gene_reaction_rule = gene_name
            self.model.add_reactions([adapter_reaction])
            self.added_reactions.add(adapter_reaction.id)
            logger.debug('Adapter reaction added: {} <--> {}'.format(metabolite.id, existing_metabolite.id))
        except Exception:  # TODO: raise a reasonable exception on cameo side if the reaction exists
            logger.debug('Adapter reaction exists: {} <--> {}'.format(metabolite.id, existing_metabolite.id))

    def add_demand_reaction(self, metabolite):
        """For metabolite in e compartment with existing exchange reaction, make it possible to consume metabolite
        by decreasing the lower bound of exchange reaction

        Parameters
        ----------
        metabolite : Metabolite
            metabolite from e compartment, f.e. melatn_e

        Returns
        -------

        """
        exchange_reaction = list(set(metabolite.reactions).intersection(self.model.exchanges))[0]
        if exchange_reaction.lower_bound >= 0:
            exchange_reaction.lower_bound = -1 if contains_carbon(metabolite) else -1000

    @staticmethod
    def annotate_new_metabolite(metabolite):
        """Find information about new metabolite in chem_prop dictionary and add it to model

        Parameters
        ----------
        metabolite : Metabolite
            new metabolite

        Returns
        -------

        """
        info = find_metabolite_info(metabolite.id)
        if info is not None:
            metabolite.formula = info['formula']
            metabolite.name = info['name']
            metabolite.annotation = info.to_dict()
        else:
            logger.debug('No formula for {}'.format(metabolite.id))
