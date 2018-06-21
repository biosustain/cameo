#!/usr/bin/env bash

set -eu

# Build on master and tags.
# $CPLEX_URL and $GH_TOKEN is defined in the Travis repository settings.
if [[ ("${TRAVIS_BRANCH}" == "master" || -n "${TRAVIS_TAG}") \
&& ("${TRAVIS_PYTHON_VERSION}" == "2.7" || "${TRAVIS_PYTHON_VERSION}" == "3.5") ]];then
    python ./.travis/load_dependency.py "${GH_TOKEN}" "cplex-python${TRAVIS_PYTHON_VERSION}.tar.gz"
    tar xzf cplex-python${TRAVIS_PYTHON_VERSION}.tar.gz
    PYTHON_VERSION=${TRAVIS_PYTHON_VERSION}
    if [ "${PYTHON_VERSION}" == "3.5" ]
    then
    PYTHON_VERSION="3.4"
    fi
    cd "cplex/python/${PYTHON_VERSION}/x86-64_linux"
    python setup.py install
    cd "${TRAVIS_BUILD_DIR}"
fi
