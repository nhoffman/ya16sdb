===================
Changes for ya16sdb
===================

0.7.2
====
* Preferring same species in type strain match_species report (GH: 47)
* Added SUCCESS or FAILED "at exit" of pipeline

0.7.1
=====
* Added "No Data" to ANI filters in Dash app
* Added link to matching ANI assembly in Dash app table

0.7
===
* Added ANI tax check annotations
* Dependency updates
* Cleaner pipeline
* Cache files are consolidated into a hidden directory
* All type strains are considered trusted unless part of do_not_trust list GL: 79
* Bug fix to correct merge with vsearch nearest type strain GH: 39

0.6.1
=====
* Bug fix to include only records with updated taxids to be re-downloaded (GL: 72)

0.6
===
* Sort and cap reference sequences by confidence (and other seq_info attributes) in blast db(s) (GL: 38)
* Redownload records that have an updated modified_date to pick up changes in pubmed_id (GL: 17)
* Do not re-download vsearch filtered sequences if training set has not been updated (GL: 44)
* Move data/ignore.txt to /molmicro/lists dir and settings.conf and change file name to do_not_download.txt (GL: 47)
* Include copy of /molmicro/common/uw/taxonomy/taxdmp.zip in dated output dir (GL: 48)
* Update record tax_ids using accession2taxid from ncbi ftp site (GL: 54)
* Create a "Trusted" list to compliment "Do Not Trust List" (GL: 57)

0.5
===
* New dash based filter outlier plots (GL: 37)
* Inclusion of tm7 Candidatus Saccharibacteria records (GL: 36)
* Updated medirect to fix "Too Many Requests" issue with NCBI (GL: 42)
* New api-key feature to increase the number of NCBI reqs/sec to 10 (GL: 43)

0.4
=======
* all paths (binaries, image locations, database credentials, etc) will be defined in a user defined file paths.conf (GH: 8, 9)
* species only records to be paritioned (GH: 13)
* new filtered/trusted/types directory (GH: 14)
* new filtered/details_out.feather, [named, named/trusted, types, types/trusted]/lineages.[csv, txt] outputs (GH: 5)

0.3
===
* new folder structure dedup -> 1200bp -> named -> filtered -> trusted

0.2
===
* deduplicated refseq and original sequences, preferring the refseq sequence (GH: 20)
* sequences determined obsolete by ncbi via esearch query are removed from database (GH: 13)
* type strain sequences in filter_outliers plot will not align to themselves as nearest type strain 
  (but might if another allele of itself is closest match) (GH: 29)
* using JUST ncbi sequence_from_type filter to mark type strain sequences (GH: 14)

0.1
=======
* imported mkrefpkg SConstruct-ncbi
