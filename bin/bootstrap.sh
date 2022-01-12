#!/bin/bash

set -e

if [[ -n "$1" ]]; then
    venv="$1"
elif [[ -n $VIRTUAL_ENV ]]; then
    venv=$VIRTUAL_ENV
else
  venv=$(pwd)/$(basename $(pwd))-env
fi

$venv/bin/pip3 install --requirement requirements.txt

mkdir src

if [[ ! -f $venv/bin/cmsearch ]]; then
  (cd src &&
   wget -nc --quiet http://eddylab.org/infernal/infernal-1.1.4-linux-intel-gcc.tar.gz && \
   tar xvf infernal-1.1.4-linux-intel-gcc.tar.gz --strip-components 2 --directory $venv/bin infernal-1.1.4-linux-intel-gcc/binaries/cmsearch)
else
  echo "cmsearch is already installed"
fi


if [[ ! -f $venv/bin/makeblastdb ]]; then
  BLAST_GZ=ncbi-blast-*-x64-linux.tar.gz
  (cd src &&
   wget -nc --user anonymous --password anonymous \
      ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/$BLAST_GZ &&
   tar tzf $BLAST_GZ | grep makeblastdb | xargs tar xzf $BLAST_GZ --strip-components 2 --directory $venv/bin)
else
  echo "makeblastdb is already installed: $(makeblastdb -version)"
fi

if [[ ! -f $venv/bin/vsearch ]]; then
  (cd src &&
   wget https://github.com/torognes/vsearch/releases/download/v2.13.0/vsearch-2.13.0-linux-x86_64.tar.gz &&
   tar tzf vsearch-2.13.0-linux-x86_64.tar.gz |
   grep bin/vsearch | xargs tar xzf vsearch-2.13.0-linux-x86_64.tar.gz --strip-components 2 --directory $venv/bin)
else
  echo "vsearch already installed"
fi

rm -rf src

echo 'All done!'
