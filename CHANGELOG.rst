===================
Changes for ya16sdb
===================

0.1-dev
=======
* imported mkrefpkg SConstruct-ncbi
* deduplicated refseq and original sequences, preferring the refseq sequence (GH: 20)
* sequences determined obsolete by ncbi via esearch query are removed from database (GH: 13)
* type strain sequences in filter_outliers plot will not align to themselves as nearest type strain 
  (but might if another allele of itself is closest match) (GH: 29)
* using JUST ncbi sequence_from_type filter to mark type strain sequences (GH: 14)
