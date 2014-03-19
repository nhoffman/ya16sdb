#!/bin/bash

set -e

if [[ -z $1 ]]; then
    venv=$(basename $(pwd))-env
else
    venv=$1
fi

MKREFPKG_DIR=$(cd $(dirname $BASH_SOURCE) && cd .. && pwd)
DEENURP_DIR=src/deenurp

if [ ! -d $DEENURP_DIR ]; then
  git clone git@github.com:fhcrc/deenurp.git $DEENURP_DIR
fi

${DEENURP_DIR}/bin/bootstrap.sh $venv

source $venv/bin/activate

pip install -U -r ${MKREFPKG_DIR}/requirements.txt
