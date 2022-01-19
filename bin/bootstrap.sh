#!/bin/bash

set -e

if [[ -n "$1" ]]; then
    venv="$1"
elif [[ -n $VIRTUAL_ENV ]]; then
    venv=$VIRTUAL_ENV
else
  venv=$(pwd)/$(basename $(pwd))-env
fi

bin_dir=$(dirname "$(readlink -f "$0")")
pipeline_dir=$(dirname $bin_dir)

$venv/bin/pip3 install --requirement $pipeline_dir/requirements.txt

mkdir src

(cd src &&
 git clone --branch 070-python3 https://github.com/fhcrc/deenurp.git  &&
 cd deenurp &&
 PYTHON=$venv/bin/python3 \
 DEENURP=$(pwd) \
 bin/bootstrap.sh $venv)

INFERNAL_GZ=infernal-1.1.4-linux-intel-gcc.tar.gz
(cd src &&
 wget -nc --quiet http://eddylab.org/infernal/$INFERNAL_GZ &&
 tar xvf $INFERNAL_GZ --strip-components 2 --directory $venv/bin infernal-1.1.4-linux-intel-gcc/binaries/cmsearch)

BLAST_GZ=ncbi-blast-*-x64-linux.tar.gz
(cd src && \
 wget --quiet -nc --user anonymous --password anonymous \
   ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/$BLAST_GZ &&
 tar tzf $BLAST_GZ | grep makeblastdb | xargs tar xzf $BLAST_GZ --strip-components 2 --directory $venv/bin)

VSEARCH_GZ=vsearch-2.13.0-linux-x86_64.tar.gz
(cd src &&
 wget --quiet https://github.com/torognes/vsearch/releases/download/v2.13.0/$VSEARCH_GZ &&
 tar tzf $VSEARCH_GZ | \
 grep bin/vsearch | xargs tar xzf $VSEARCH_GZ --strip-components 2 --directory $venv/bin)

(cd src &&
 wget --quiet https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz &&
 tar tzf edirect.tar.gz | grep -E 'esearch$|ecommon.sh|nquire$|xtract$' | xargs tar xzf edirect.tar.gz --strip-components 1 --directory $venv/bin)


rm -rf src

echo 'All done!'
