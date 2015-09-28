# run the following command in the root directory of cameo to (re-)generate the API .rst files using sphinx-apidoc
# -f will overwrite any existing files
sphinx-apidoc --no-toc -f -o docs/apidoc_output/ cameo cameo/strain_design/heuristic/multiprocess
