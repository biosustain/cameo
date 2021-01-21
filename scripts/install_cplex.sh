#!/usr/bin/env bash

set -eu

# Build on master and tags.
if [[ "${PYTHON_VERSION}" == "3.6" ]];then
    # this should not be logged by GH actions by make it -s silent just in case
    curl -s -o cplex.zip -L $CPLEX_SECRET 
    unzip -q cplex.zip
    ROOT_DIR=`pwd`
    cd "python/${PYTHON_VERSION}/x86-64_linux"
    python setup.py install
    cd $ROOT_DIR
fi
