from cobra.util import add_exchange


def get_reaction_for(model, value, add=True):
    """Get or create a reaction for a metabolite or a reaction.

    If value is a Metabolite or a Metabolite id, return any already existing demand or exchange reaction.
    If *add* is true, add a demand reaction if it does not already exist.

    Parameters
    ----------
    model : cameo.core.SolverBasedModel
        The model to for which to get / create a reaction
    value: str, Reaction or Metabolite
        A reaction identifier, a Reaction or a Metabolite for which an exchange reaction is to be created.
    add: bool
        Adds a demand reaction for a metabolite if a metabolite is found for *value*

    Returns
    -------
    Reaction

    Raises
    ------
    KeyError
        If *value* does not match any Reaction or Metabolite

    """
    try:
        reactions = model.reactions.get_by_any(value)
    except (ValueError, KeyError, TypeError):
        metabolite = model.metabolites.get_by_any(value)[0]
        reactions = model.reactions.query("^(EX|DM)_{}$".format(metabolite.id))
        if len(reactions) == 0:
            if add:
                reactions = [add_exchange(model, metabolite)]
            else:
                raise KeyError('Invalid target %s' % value)
    return reactions[0]
