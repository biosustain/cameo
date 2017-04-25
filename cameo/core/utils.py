import six
import csv

from pandas import DataFrame


def get_reaction_for(model, value, add=True):
    """Get or create a reaction for a metabolite or a reaction.

    If value is a Metabolite or a Metabolite id, return any already existing demand or exchange reaction.
    If *add* is true, add a demand reaction if it does not already exist.

    Parameters
    ----------
    model : cobra.Model
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
                reactions = [model.add_boundary(metabolite, type='demand')]
            else:
                raise KeyError('Invalid target %s' % value)
    return reactions[0]


def medium(model):
    """Current medium for this model."""
    reaction_ids = []
    reaction_names = []
    lower_bounds = []
    upper_bounds = []
    for ex in model.exchanges:
        metabolite = list(ex.metabolites.keys())[0]
        coeff = ex.metabolites[metabolite]
        if coeff * ex.lower_bound > 0:
            reaction_ids.append(ex.id)
            reaction_names.append(ex.name)
            lower_bounds.append(ex.lower_bound)
            upper_bounds.append(ex.upper_bound)

    return DataFrame({'reaction_id': reaction_ids,
                      'reaction_name': reaction_names,
                      'lower_bound': lower_bounds,
                      'upper_bound': upper_bounds},
                     index=None, columns=['reaction_id', 'reaction_name', 'lower_bound', 'upper_bound'])


def load_medium(model, medium_def, copy=False, delimiter="\t"):
    """
    Loads a medium into the model. If copy is true it will return
    a copy of the model. Otherwise it applies the medium to itself.
    Supported formats
    TODO

    Parameters
    ----------
    model : cobra.Model
        The model to load medium for
    medium_def: str, pandas.DataFrame, dict.
        The medium to load
    copy: boolean, optional
        If True copies the model, otherwise the changes will happen inplace.
    delimiter: str
        Only if loading the medium from a file.

    Returns
    -------
    cobra.Model
        If copy=True, returns a copy of the model.

    """

    if copy:
        model = model.copy()
    else:
        model = model
    if isinstance(medium_def, dict):
        _load_medium_from_dict(model, medium_def)
    elif isinstance(medium_def, DataFrame):
        _load_medium_from_dataframe(model, medium_def)
    elif isinstance(medium_def, six.string_types):
        _load_medium_from_file(model, medium_def, delimiter=delimiter)
    else:
        raise AssertionError("input type (%s) is not valid" % type(medium))

    return model


def _load_medium_from_dict(model, medium_def):
    assert isinstance(medium_def, dict)
    for ex_reaction in model.exchanges:
        ex_reaction.lower_bound = medium_def.get(ex_reaction.id, 0)


def _load_medium_from_file(model, file_path, delimiter="\t"):
    this_medium = {}

    with open(file_path, "rb") as csv_file:
        reader = csv.reader(csv_file, delimiter=delimiter)

        for row in reader:
            model.reactions.get_by_id(row[0])
            this_medium[row[0]] = row[1]

    _load_medium_from_dict(model, this_medium)


def _load_medium_from_dataframe(model, medium_df):
    assert isinstance(medium_df, DataFrame)
    for ex_reaction in model.exchanges:
        if ex_reaction.id in medium_df.reaction_id.values:
            medium_row = medium_df[medium_df.reaction_id == ex_reaction.id]
            ex_reaction.lower_bound = medium_row.lower_bound.values[0]
        else:
            ex_reaction.lower_bound = 0
