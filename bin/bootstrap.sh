#!/bin/bash

# set up a virtualenv

# `DEENURP_BRANCH` optionally defines the git branch to check out from
# the deenurp repository (default master)

# `WHEELHOUSE` optionally provides a path or url containing python
# wheels for `pip install --find-links=...`

set -e

if [[ -z $1 ]]; then
    venv=$(basename $(pwd))-env
else
    venv=$1
fi

MKREFPKG_DIR=$(cd $(dirname $BASH_SOURCE) && cd .. && pwd)
DEENURP=src/deenurp
mkdir -p src

if [[ -z "$DEENURP_BRANCH" ]]; then
    DEENURP_BRANCH=master
fi

if [ ! -d $DEENURP ]; then
    git clone -b "$DEENURP_BRANCH" git@github.com:fhcrc/deenurp.git $DEENURP
fi

${DEENURP}/bin/bootstrap.sh $venv

source $venv/bin/activate

if [[ -z "$WHEELHOUSE" ]]; then
    pip install -U -r ${MKREFPKG_DIR}/requirements.txt
else
    pip install -U --use-wheel --find-links="$WHEELHOUSE" -r ${MKREFPKG_DIR}/requirements.txt
fi



