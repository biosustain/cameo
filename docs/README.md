Cameo's source is documented using numpy style docstrings (http://sphinxcontrib-napoleon.readthedocs.org/en/latest/example_numpy.html#example-numpy).

Cameo tutorials come in the form of jupyter notebooks (./notebooks) that can be converted to .rst files using the following command.

```make ipy2rst```

The notebooks folder is a git submodule pointing to [http://github.com/biosustain/cameo-notebooks](http://github.com/biosustain/cameo-notebooks).

To regenerate the API docs run

```make apidoc```

