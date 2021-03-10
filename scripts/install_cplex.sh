#!/usr/bin/env bash

set -eu

# Build on master and tags.
if [[ "${PYTHON_VERSION}" == "3.6" ]];then
    # this should not be logged by GH actions by make it -s silent just in case
    curl -s -o cplex36.tar.gz -L $CPLEX_SECRET 
    tar xzf cplex36.tar.gz
    ROOT_DIR=`pwd`
    cd "cplex_128/python/${PYTHON_VERSION}/x86-64_linux"
    python setup.py install
    cd $ROOT_DIR
fi
