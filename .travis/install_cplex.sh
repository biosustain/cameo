#!/usr/bin/env bash

set -eu

# Build on master and tags.
# $CPLEX_URL is defined in the Travis repository settings.
if [[ ("${TRAVIS_BRANCH}" == "master" || -n "${TRAVIS_TAG}") \
&& ("${TRAVIS_PYTHON_VERSION}" == "2.7" || "${TRAVIS_PYTHON_VERSION}" == "3.5") ]];then
    python ./.travis/load_dependency.py "${GH_TOKEN}" "cplex-python${TRAVIS_PYTHON_VERSION}.tar.gz"
    tar xzf cplex-python${TRAVIS_PYTHON_VERSION}.tar.gz
    cd "cplex/python/3.4/x86-64_linux"
    python setup.py install
    cd "${TRAVIS_BUILD_DIR}"
fi
