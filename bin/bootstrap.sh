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
    DEENURP_BRANCH=master
fi

if [[ ! -d $DEENURP ]]; then
    git clone -b "$DEENURP_BRANCH" git@github.com:fhcrc/deenurp.git $DEENURP
fi

# install python envs
python3 -m venv $venv
virtualenv --quiet --python python2 $venv
${DEENURP}/bin/bootstrap.sh $venv
source $venv/bin/activate
$venv/bin/pip2 install --requirement ${MKREFPKG_DIR}/requirements2.txt
$venv/bin/pip3 install --requirement ${MKREFPKG_DIR}/requirements3.txt

if [[ ! -f $venv/bin/makeblastdb ]]; then
  BLAST_GZ=ncbi-blast-*-x64-linux.tar.gz
  (cd src &&
   wget -nc --user anonymous --password $(git config user.email) \
      ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/$BLAST_GZ &&
   tar tzf $BLAST_GZ | grep makeblastdb | xargs tar xzf $BLAST_GZ --strip-components 2 --directory $venv/bin)
else
  echo "makeblastdb is already installed: $(makeblastdb -version)"
fi

if [[ ! -f $venv/bin/esearch ]]; then
  (cd src &&
   wget ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz &&
   tar tzf edirect.tar.gz |
   grep -E "esearch|edirect.pl" |
   xargs tar xzf edirect.tar.gz --strip-components 1 --directory $venv/bin)
else
  echo "esearch already installed: $(esearch -version)"
fi

echo 'All done!'
