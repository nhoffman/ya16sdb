#!/bin/bash

set -e

MKREFPKG_DIR=$(cd $(dirname $BASH_SOURCE) && cd .. && pwd)
DEENURP_DIR=${MKREFPKG_DIR}/deenurp

if [ ! -d $DEENURP_DIR ]; then
  git clone git@github.com:fhcrc/deenurp.git $DEENURP_DIR
fi

${DEENURP_DIR}/bin/bootstrap.sh

pip install -U -r ${MKREFPKG_DIR}/requirements.txt
