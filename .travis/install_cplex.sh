#!/usr/bin/env bash

set -eu

# Build on master and tags.
if [[ ("${TRAVIS_BRANCH}" == "master" || -n "${TRAVIS_TAG}") \
&& ("${TRAVIS_PYTHON_VERSION}" == "2.7" || "${TRAVIS_PYTHON_VERSION}" == "3.5") ]];then
    curl -O $CPLEX_SECRET  # check lastpass
    tar xzf cplex-python3.6.tar.gz
    PYTHON_VERSION=${TRAVIS_PYTHON_VERSION}
    if [ "${PYTHON_VERSION}" == "3.5" ]
    then
    PYTHON_VERSION="3.4"
    else
    PYTHON_VERSION="2.6"
    fi
    cd "cplex/python/${PYTHON_VERSION}/x86-64_linux"
    python setup.py install
    cd "${TRAVIS_BUILD_DIR}"
fi