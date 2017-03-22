import cameo


def add_exchange(model, metabolite, demand=True, prefix='DM_', bound=1000.0):
    """Add an exchange reaction for a metabolite (demand=TRUE: metabolite --> Ø or demand=False: 0 --> metabolite )

    Parameters
    ----------
    model : cameo.core.SolverBasedModel
        The model to add the exchange reaction to.
    metabolite : cameo.core.Metabolite
    demand : bool, optional
        True for sink type exchange, False for uptake type exchange
    prefix : str, optional
        A prefix that will be added to the metabolite ID to be used as the demand reaction's ID (defaults to 'DM_').
    bound : float, optional
        Upper bound for sink reaction / lower bound for uptake (multiplied by -1)

    Returns
    -------
    Reaction
        The created exchange reaction.
    """
    reaction_id = str(prefix + metabolite.id)
    name = "Exchange %s" % metabolite.name if prefix != "DM_" else "Demand %s" % metabolite.name
    if reaction_id in model.reactions:
        raise ValueError("The metabolite already has a demand reaction.")

    reaction = cameo.core.Reaction(id=reaction_id, name=name)
    reaction.add_metabolites({metabolite: -1})
    reaction.bounds = (0, bound) if demand else (-bound, 0)
    model.add_reactions([reaction])
    return reaction


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
