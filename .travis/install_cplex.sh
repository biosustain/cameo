#!/usr/bin/env bash

set -exu

# Build on master and tags.
# $CPLEX_URL is defined in the Travis repository settings.
if [[ ("${TRAVIS_BRANCH}" == "master" || -n "${TRAVIS_TAG}") \
&& ("${TRAVIS_PYTHON_VERSION}" == "2.7" || "${TRAVIS_PYTHON_VERSION}" == "3.5") ]];then
    curl -L "${CPLEX_URL}" -o cplex.tar.gz
    tar xzf cplex.tar.gz
    cd "cplex/python/${TRAVIS_PYTHON_VERSION}/x86-64_linux"
    python setup.py install
    cd "${TRAVIS_BUILD_DIR}"
fi
