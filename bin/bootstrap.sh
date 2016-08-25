#!/bin/bash

# set up a virtualenv

# `DEENURP_BRANCH` optionally defines the git branch to check out from
# the deenurp repository (default master)

# `PIP_FIND_LINKS` optionally provides a path or url containing python
# wheels https://pip.pypa.io/en/latest/user_guide.html#environment-variables

set -e

if [[ -n "$1" ]]; then
    venv="$1"
elif [[ -n $VIRTUAL_ENV ]]; then
    venv=$VIRTUAL_ENV
else
  venv=$(pwd)/$(basename $(pwd))-env
fi

MKREFPKG_DIR=$(cd $(dirname $BASH_SOURCE) && cd .. && pwd)
DEENURP=src/deenurp
mkdir -p src

if [[ -z "$DEENURP_BRANCH" ]]; then
    DEENURP_BRANCH=v0.1.8
fi

if [[ ! -d $DEENURP ]]; then
    git clone -b "$DEENURP_BRANCH" git@github.com:fhcrc/deenurp.git $DEENURP
fi

${DEENURP}/bin/bootstrap.sh $venv

if [ ! -f $venv/bin/makeblastdb ]; then
  BLAST_GZ=ncbi-blast-*-x64-linux.tar.gz
  (cd src &&
   wget --quiet -nc ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/$BLAST_GZ &&
   tar tzf $BLAST_GZ |
   grep makeblastdb |
   xargs tar xzf $BLAST_GZ --strip-components 2 --directory $venv/bin)
else
    echo "makeblastdb is already installed: $(makeblastdb -version)"
fi

source $venv/bin/activate

# set PIP_FIND_LINKS to use wheels https://pip.pypa.io/en/latest/user_guide.html#environment-variables
pip install --upgrade --requirement ${MKREFPKG_DIR}/requirements.txt

virtualenv --relocatable $venv
