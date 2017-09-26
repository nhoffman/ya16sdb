========================
 NCBI curation pipeline
========================

https://bitbucket.org/uwlabmed/mkrefpkg

Pipeline executed using ``scons -f SConstruct-ncbi``

output directories
==================

``output-ncbi`` contains the pipleine output, with intermediate files
in the top level, and a dated subdirectory for each run of the
pipeline. The most recent output is typically symlinked to
``output-ncbi/LATEST``. The layout of each of these subdirectories is
as follows::

  % ls -1 output-ncbi/LATEST/
  1200bp/
  appended_info.csv
  git_version.txt
  is_type_info.csv
  ncbi.log
  ncbi.txt
  new/
  pubmed_ids.csv
  references.csv
  requirements.txt
  seq_info.csv
  seqs.fasta
  types.txt
  versions.txt

* in the top level, seq_info.csv and seqs.fasta contain the entire
  corpus of downloaded 16 sequences.
* ``new/`` - new records downloaded relative to the last run
  (typically not used for any downstream applications)
* ``1200bp/`` - records limited to those >= 1200bp
* ``1200bp/valid/`` - records in ``1200bp/`` with species-level names
  not matching the regular expression defined here:
  https://github.com/fhcrc/taxtastic/blob/master/taxtastic/ncbi.py#L79
* ``1200bp/valid/dedup`` - "valid" records deduplicated by first
  grouping by accession, then choosing a representative of each set of
  identical sequences (ie, same md5 hash),
* ``1200bp/valid/filtered`` - deduplicated sequences from which
  outliers have been discarded (using ``deenurp filter_outliers``).
* ``1200bp/valid/types`` - type strains, as identified by NCBI
  (correspond to accessions downloaded with the genbank query filter
  ``sequence_from_type[Filter]``)

As an (extreme) example of the incremental filtering process::

  LATEST % find . -name seq_info.csv | xargs grep -c 'Streptococcus pneumoniae' seq_info.csv
  seq_info.csv:12446
  ./new/seq_info.csv:0
  ./new/valid/seq_info.csv:0
  ./new/invalid/seq_info.csv:0
  ./seq_info.csv:12446
  ./1200bp/seq_info.csv:11393
  ./1200bp/valid/seq_info.csv:11393
  ./1200bp/valid/dedup/seq_info.csv:11208
  ./1200bp/valid/filtered/seq_info.csv:11045
  ./1200bp/valid/types/seq_info.csv:58
