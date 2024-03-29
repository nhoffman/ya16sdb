===============================
 Yet Another 16S rRNA database
===============================

ya16sdb is a pipeline for downloading, curating, and annotating a
database of bacterial 16S rRNA sequences. This repository also
implements a web application (https://ya16sdb.labmed.uw.edu/) that can
be used to visualize the distance-based relationships among sequences
for a given species.

The purpose of the project is to provide a high quality source of
bacterial 16S rRNA sequences that is up to date with NCBI, in a format
that is useful as an input for various bioinformatics pipleines such
as blast searching, phylogenetic reference set creation,
sequence-based taxonomic assignment, etc.

Project information
===================

This project is a product of ongoing research interests of Noah
Hoffman (https://faculty.washington.edu/ngh2/home/pages/software.html)
at the University of Washington in the Department of Laboratory
Medicine.

Christopher Rosenthal is the primary author of the pipeline.

The pipeline heavily relies on *taxtastic*
(https://github.com/fhcrc/taxtastic) and *deenurp*
(https://github.com/fhcrc/deenurp), both of which began as
collaborations with Erick Matsen at The Fred Hutchinson Cancer
Research Center in Seattle, WA.

Please cite this project as

Rosenthal C and Hoffman NG. 2019. ya16sdb: a pipeline for creating a
collection of high-quality bacterial 16S rRNA sequences from
NCBI. Version 0.6.1. University of Washington. https://github.com/nhoffman/ya16sdb

Overview
========

At a high level, this pipeline does the following:

* Downloads annotation for all available sequence records from the
  NCBI matching search terms for 16S rRNA.
* Retrieves sequence records for corresponding full length (or near
  full-length) 16S rRNA genes; this involves extracting subsequences
  from genome sequences or contigs.
* Ensures that all records are 16S rRNA genes
* Ensures that sequences are in a consistent orientation.
* Identifies the taxonomic lineage of each record.
* Annotates records as a "type strain" (according to NCBI's definition
  of type strain), "published" (annotation has an accompanying PubMed
  ID), "refseq" (belonging to the Genbank refseq collection), or
  "direct" (direct submissions).
* Discards records likely to be mis-annotated using ``deenurp filter-outliers``.
* Provides various subsets of annotated sequences. Each record subset
  provides sequence metadata, sequences, taxonomic lineages, and a
  blast database. For example:

  * only records with taxonomic name consistent with species-level classifications
  * type strains only
  * outliers removed
  * downsampled to a subset of sequences for each species, prioritizing type strains and "published" records.

Docker
======

Docker image can be built with the following:
::

  docker build --tag ya16sdb:latest .

Once a Docker image has been built a Singularity image can be built using the docker daemon:
::

  singularity build ya16sdb.img docker-daemon://ya16sdb:latest

A Singularity image can also be built using a Singularity Docker container:
::

  docker run --volume /var/run/:/var/run/ --volume $(pwd):$(pwd) --workdir $(pwd) singularity:latest build ya16sdb.img docker-daemon://ya16sdb:latest

Pipeline execution
------------------

The virtual containers have a predefined entry point to the SConstruct pipeline file.

To execute using Docker just a settings.conf file is required and can be run as follows:
::

  docker run --volume $(pwd):$(pwd) --workdir $(pwd) ya16sdb:latest

And with Singularity
::

  singularity run --bind $(pwd) --pwd $(pwd) ya16sdb.img
