"""
Download and curate the 16S rRNA sequences from NCBI

NOTE: run with `scons -j 1` only
"""
import atexit
import configparser
import csv
import errno
import os
import sys
import time
import warnings

from SCons.Script import ARGUMENTS, Environment, GetBuildFailures, Depends

this_dir = os.path.dirname((lambda x: x).__code__.co_filename)
venv = os.environ.get('VIRTUAL_ENV')
if not venv:
    warnings.warn('No active virtualenv detected, using system environment')
settings_file = ARGUMENTS.get('settings', 'settings.conf')
if not os.path.isfile(settings_file):
    sys.exit("Can't find settings.conf")
conf = configparser.ConfigParser()
conf.optionxform = str  # preserve key case-sensitivity
conf.read(settings_file)
settings = conf['DEFAULT']
if 'ncbi_conf' in settings:
    conf.read(settings['ncbi_conf'])
else:
    conf.read(os.path.join(this_dir, 'ncbi.conf'))
true_vals = ['t', 'y', '1']
release = ARGUMENTS.get('release', 'no').lower()[0] in true_vals
out = os.path.join(settings['outdir'], time.strftime('%Y%m%d'))
log = os.path.join(out, 'log.txt')
cachedir = settings['cachedir']


def blast_db(env, sequence_file, output_base):
    '''
    Create a blast database and file md5sum
    '''
    extensions = ['.nhr', '.nin', '.nsq']
    blast_out = env.Command(
        target=[output_base + ext for ext in extensions],
        source=sequence_file,
        action='makeblastdb -dbtype nucl '
               '-in $SOURCE -out ' + output_base)
    env.Command(
        target=output_base,
        source=blast_out,
        action='md5sum $SOURCES > $TARGET')
    return blast_out


def taxonomy(fa, info, path):
    """
    Make filtered taxtable with ranked columns and no_rank rows

    METHODS:
    """
    taxtable = env.Command(
        target=os.path.join(path, 'taxonomy.csv'),
        source=info,
        action='taxit taxtable '
               '--seq-info $SOURCE '
               '--out $TARGET '
               '$tax_url')

    """
    Taxtable output replacing tax_ids with taxnames

    2.
    """
    lineages = env.Command(
        target=os.path.join(path, 'lineages.csv'),
        source=[taxtable, info],
        action='taxit lineage_table --csv-table $TARGET $SOURCES')

    """
    Mothur output - https://mothur.org/wiki/Taxonomy_File

    3.
    """
    mothur = env.Command(
        target=os.path.join(path, 'lineages.txt'),
        source=[taxtable, info],
        action='taxit lineage_table --taxonomy-table $TARGET $SOURCES')

    blast = blast_db(env, fa, os.path.join(path, 'blast'))

    return taxtable, lineages, mothur, blast


cfiles = {
    'genbank_cache': os.path.join(cachedir, 'records.gb'),
    'outliers_cache': os.path.join(cachedir, 'filter_outliers.csv'),
    'pubmed_info_cache': os.path.join(cachedir, 'pubmed_info.csv'),
    'seqs_cache': os.path.join(cachedir, 'seqs.fasta'),
    'seq_info_cache': os.path.join(cachedir, 'seq_info.csv'),
    'records_cache': os.path.join(cachedir, 'records.txt'),
    'references_cache': os.path.join(cachedir, 'references.csv'),
    'refseq_info_cache': os.path.join(cachedir, 'refseq_info.csv')
    }

# create cache files if they do not exist
for k, v in cfiles.items():
    if not os.path.exists(v):
        d = os.path.dirname(v)
        if not os.path.exists(d):
            os.makedirs(d)
        open(v, 'w').close()

# use default ENV and add/append values after
env = Environment(
    log=log,
    out=out,
    pipeline=this_dir,
    tax_url=os.path.abspath(settings['taxonomy']),
    **settings
)
env.PrependENVPath('PATH', os.path.join(this_dir, 'bin'))

for k, v in os.environ.items():
    if k.startswith('MEFETCH_'):
        env['ENV'][k] = v
env['ENV'].update(conf['ENV'])
env['ENV']['MEFETCH_DB'] = 'nucleotide'
env['ENV']['MEFETCH_LOG'] = log
env['ENV']['MEFETCH_MODE'] = 'text'

env.Decider('MD5-timestamp')

"""
METHODS: NCBI sequence accessions indexed as Bacteria and Candidatus
Saccharibacteria 16S and rRNA are queried using the NCBI esearch Entrez
Programming tool and downloaded using mefetch. Mefetch is a command
line, multithreaded Python application that breaks up and downloads large
efetch requests in parallel.
"""
classified = env.Command(
    source=None,
    target='$out/ncbi/classified.txt',
    action='esearch -db nucleotide -query "$classified" | '
           'mefetch -vv -format acc -out $TARGET -reqs 3'
    )

"""
Candidatus Saccharibacteria
https://gitlab.labmed.uw.edu/molmicro/mkrefpkg/issues/36
"""
tm7 = env.Command(
    source=None,
    target='$out/ncbi/tm7.txt',
    action='esearch -db nucleotide -query "$tm7" | '
           'mefetch -vv -format acc -out $TARGET')

"""
Check the cache for last download_date and download list of modified
records in order to re-download modified records that we have previously
downloaded.

METHODS: Records with a modified date after the last download date are
re-downloaded and updated in the record cache.
"""
si = csv.DictReader(open(env.subst(cfiles['seq_info_cache'])))
si = [time.strptime(r['download_date'], '%d-%b-%Y') for r in si]
if si:
    env['download_date'] = time.strftime('%Y/%m/%d', sorted(si)[-1])
else:
    env['download_date'] = time.strftime('%Y/%m/%d')
modified = env.Command(
    source=None,
    target='$out/ncbi/modified.txt',
    action='esearch -db nucleotide -query "($classified OR $tm7) AND '
           '$download_date[Modification Date] : 3000[Modification Date]" | '
           'mefetch -vv -format acc -out $TARGET -reqs 3')

"""
type strains records
NCBI web blast uses `sequence_from_type[Filter]` so we will use that
http://www.ncbi.nlm.nih.gov/news/01-21-2014-sequence-by-type/

METHODS: Type strains are queried using the NCBI filter sequence_from_type and
downloaded using mefetch.
"""
types = env.Command(
    source=None,
    target='$out/ncbi/types.txt',
    action='esearch -db nucleotide -query "$types" | '
           'mefetch -vv -format acc -out $TARGET -reqs 3')

"""
Trim accession2taxid with 16s records and update taxids

METHODS: For each potential 16s accession taxonomy ids are downloaded from
NCBI's ftp accession2taxid site.  Accessions without a taxonomy id are dropped.
"""
accession2taxid = env.Command(
    target='$out/ncbi/accession2taxid.csv',
    source=[settings['accession2taxid'], classified, tm7],
    action='accession2taxid.py --out $TARGET $SOURCES')

"""
METHODS: Taxonomy ids are updated using the Taxtastic update_taxids subcommand
Taxtastic is a command line Python application with a suite of subcommands
that help update and manage NCBI taxonomy ids.
"""
accession2taxid = env.Command(
    target='$out/ncbi/accession2update_taxids.csv',
    source=accession2taxid,
    action='taxit update_taxids '
           '--outfile $TARGET '
           '--unknown-action drop '
           '$SOURCE $tax_url')

"""
Create a list of cached records removing:

METHODS: Modified records
2. Records not included in esearch query
3. Records with new tax_ids

Modified and records with new tax_ids will be redownloaded

FIXME: call this cache_in.py (cache_out.py will replace refresh.py)

METHODS: A list of potential 16s accessions is compiled for download if they are not
present in or have a different taxonomy id in the record cache or have been
modified.  A list of accessions can be automatically excluded from download
can be supplied.
"""
cache = env.Command(
    target='$out/ncbi/cache.txt',
    source=[accession2taxid, cfiles['seq_info_cache'],
            modified, cfiles['records_cache']],
    action='cache.py --out $TARGET $SOURCES')

"""
Find records to download and filter ncbi master list
with tax_ids that exist in our taxonomy
"""
download = env.Command(
    target='$out/ncbi/download.txt',
    source=[accession2taxid, cache, settings['do_not_download']],
    action='download.py --out $TARGET $SOURCES')

"""
METHODS: Potential 16s record feature tables are downloaded using mefetch.  Record
coordinates cooresponding to feature rrna product=16S ribosomal RNA are
selected using medirect ftract.  Ftract is part of the medirect application
suite that parses NCBI feature tables for string patterns.  It parses feature
tables in raw text format and is a faster alternative to NCBI's xml xtract
utility.
"""
coordinates = env.Command(
    target='$out/ncbi/genbank.csv',
    source=download,
    action='mefetch -vv '
           '-failed $out/ncbi/ft_failed.txt '
           '-format ft '
           '-id $SOURCE '
           '-retmax 10 '
           '| ftract -feature "rrna:product:16S ribosomal RNA" -out $TARGET'
           )
env.Precious(coordinates)

"""
download genbank records

-retmax 1 when reading from -csv file

METHODS: 16s records are downloaded at the feature coordinates using mefetch.
"""
gbs = env.Command(
    target='$out/ncbi/download.gb',
    source=coordinates,
    action='mefetch -vv '
           '-csv '
           '-failed $out/ncbi/failed.gb '
           '-format gb '
           '-id $SOURCE '
           '-out $TARGET')

"""
METHODS: Downloaded Genbank records are parsed into fasta and csv formats. For a full
list of parsed annotations see Table [].
"""
today = time.strftime('%d-%b-%Y')
fa, seq_info, pubmed_info, references, refseq_info = env.Command(
    target=['$out/ncbi/extract_genbank/seqs.fasta',
            '$out/ncbi/extract_genbank/seq_info.csv',
            '$out/ncbi/extract_genbank/pubmed_info.csv',
            '$out/ncbi/extract_genbank/references.csv',
            '$out/ncbi/extract_genbank/refseq_info.csv'],
    source=gbs,
    action='extract_genbank.py $SOURCE ' + today + ' $TARGETS')

"""
Record versions returned from esearch that had no actual 16s features or were
below `ftract -min-length`

FIXME: write script that adds the annotations to the info with no seqnames

METHODS: Record accession returned from the esearch query but do not contain
rRNA product "16S ribosomal RNA" feautres are saved in a no_features.txt
file to avoid redownloaded in future pipeline runs.
"""
no_features = env.Command(
    target='$out/ncbi/no_features.txt',
    source=[seq_info, download],
    action=('csvcut.py --columns version ${SOURCES[0]} | '
            'tail -n +2 | '
            'grep '
            '--file - '
            '--fixed-strings '
            '--invert-match '
            '--line-regexp '
            '${SOURCES[1]} > $TARGET '
            '|| true'))

"""
Remove new records with tax_ids not in our taxonomy database

METHODS: Taxonomy ids are again updated this time from the downloaded Genbank
extracted taxonomy ids.  Records where are taxonomy id is unknown are
discarded.
"""
fa, seq_info = env.Command(
    target=['$out/ncbi/extract_genbank/update_taxids/seqs.fasta',
            '$out/ncbi/extract_genbank/update_taxids/seq_info.csv'],
    source=[fa, seq_info],
    action='taxit update_taxids '
           '--unknown-action drop '
           '${SOURCES[1]} $tax_url | '
           'partition_refs.py ${SOURCES[0]} - $TARGETS')

"""
cmsearch new sequences against rfam model

METHODS: Record sequence orientation and quality are assessed using Infernal cmsearch
against a curated 16S coveriance model available from the publicly available
rfam database.  Sequences in 3'-5' orientation are put into their
reverse complement 5'-3' orientation.  Sequences with an E value
below 0.01 are considered of low sequence quality or missannotated
and are discarded.
"""
cmsearch = env.Command(
    target='$out/ncbi/extract_genbank/update_taxids/cmsearch/table.tsv',
    source=['$pipeline/data/SSU_rRNA_bacteria.cm', fa],
    action='cmsearch --cpu 14 -E 0.01 --hmmonly -o /dev/null '
           '--tblout $TARGET $SOURCES || true')

"""
Fix record orientations
"""
fa, seq_info = env.Command(
    target=['$out/ncbi/extract_genbank/update_taxids/cmsearch/seqs.fasta',
            '$out/ncbi/extract_genbank/update_taxids/cmsearch/seq_info.csv'],
    source=[cmsearch, fa, seq_info],
    action='cmsearch.py $SOURCES $TARGETS')

"""
Refresh/append with older records

METHODS: Drop seqnames missing either a sequence or row in seq_info
2. Append seqs, seq_info, pubmed_info and references to previous data set
3. Drop records not in the ncbi records.txt file
4. Drop sequences that have a refseq equivalent
5. Deduplicate pubmeds and references
6. Copy and cache

FIXME: Call this cache_out.py ... with --cachedir argument
FIXME: Move refseq dedup to partition_refs.py --drop-duplicates-sequences

METHODS: Newly downloaded records are combined and added to the record cache.
"""
fa, seq_info, pubmed_info, _, refseq_info, _ = env.Command(
    target=['$out/ncbi/seqs.fasta',
            '$out/ncbi/seq_info.csv',
            '$out/pubmed_info.csv',
            '$out/references.csv',
            '$out/refseq_info.csv',
            '$out/records.txt'],
    source=[accession2taxid,
            fa, cfiles['seqs_cache'],
            seq_info, cfiles['seq_info_cache'],
            pubmed_info, cfiles['pubmed_info_cache'],
            references, cfiles['references_cache'],
            refseq_info, cfiles['refseq_info_cache'],
            no_features, cfiles['records_cache']],
    action=['refresh.py $SOURCES $TARGETS',
            # cache
            'cp ${TARGETS[0]} ' + cfiles['seqs_cache'],
            'cp ${TARGETS[1]} ' + cfiles['seq_info_cache'],
            'cp ${TARGETS[2]} ' + cfiles['pubmed_info_cache'],
            'cp ${TARGETS[3]} ' + cfiles['references_cache'],
            'cp ${TARGETS[4]} ' + cfiles['refseq_info_cache'],
            'cp ${TARGETS[5]} ' + cfiles['records_cache']])

"""
append new records to global list

TODO: create bin/dedup_gb.py for genbank record cache maintenance
"""
env.Command(
    target=None,
    source=gbs,
    action='cat $SOURCE >> ' + cfiles['genbank_cache'])

"""
METHODS: A csv taxtable is created for classification using pplacer or Moose
classifiers.
"""
taxtable = env.Command(
    target='$out/taxonomy.csv',
    source=seq_info,
    action='taxit -v taxtable --seq-info $SOURCE --out $TARGET $tax_url')

"""
Map for WGS records without a refseq assembly accession

METHODS:
"""
asm = env.Command(
    target='$out/ani/assembly_summary_genbank.txt',
    source=None,
    action='wget '
           '--output-document $TARGET '
           '--quiet '
           '--retry-on-http-error 403 '
           '--tries 100 '
           'https://ftp.ncbi.nlm.nih.gov/'
           'genomes/genbank/assembly_summary_genbank.txt')

"""
The ANI tax check report

METHODS: ANI Tax Check reports are downloaded from NCBI and added to the
seq_info.csv file using an genbank to assembly accession map file if needed.
"""
ani = env.Command(
    target='$out/ani/ANI_report_prokaryotes.txt',
    source=None,
    action='wget '
           '--output-document $TARGET '
           '--quiet '
           '--retry-on-http-error 403 '
           '--tries 100 '
           'https://ftp.ncbi.nlm.nih.gov/'
           'genomes/ASSEMBLY_REPORTS/ANI_report_prokaryotes.txt')

"""
Create feather file and initial columns
FIXME: Switch $SOURCE and $TARGET around for each script

METHODS: The seq_info.csv file is converted to a feather file (make this sqlite file)
Non Genbank annotations are added including for each record species and genus
taxonomy ids and taxonomy names.  A boolean column is_type is added and record
set to true if record accession appears in the downloaded sequence_from _type
list of accessions. A boolean column is_published is added and set to true if
record has an associated pubmed_id.  A boolean column is_refseq is added and
set to true if record accession contains an underscore indicating the record
is a refseq.  A boolean is_valid is added and set to true if a record species
name passes the Taxtastic classified regex expression.  Unclassified species
names, in general, are typically environmental or non-specific taxonomy names
[]. A string column confidence is added and each record is designated as either
type, published, refseq or direct depending on which true attribute comes first
to best describes each record. A column seqhash is added and for each record
sequence a sha1 hash is generated and added to this column and used for
deduplication purposes.  Finally, if a cap on species representatives is
specified then an optional sort_by value can be set to prioritize records by
column attributes.  By default, the number of species representatives is capped
at 5,000 records and sorted by columns is_type, is_published, is_refseq,
ambig_count, modified_date and download_date.
"""
feather = env.Command(
    target='$out/seq_info.feather',
    source=[
        seq_info,
        taxtable,
        types,
        pubmed_info,
        refseq_info,
        ani,
        asm,
        fa],
    action=['to_feather.py ${SOURCES[0]} $TARGET',
            # Add 'species' tax_id, species_name, genus taxid and genus name
            'taxonomy.py $TARGET ${SOURCES[1]}',
            # https://github.com/nhoffman/ya16sdb/issues/11
            'is_type.py $TARGET ${SOURCES[2]}',
            'is_published.py $TARGET ${SOURCES[3]}',
            # add is_refseq and original column with refseq accession
            'is_refseq.py $TARGET ${SOURCES[4]}',
            # taxtastic GH: 172
            'is_valid.py $TARGET $tax_url',
            # https://gitlab.labmed.uw.edu/uwlabmed/mkrefpkg/issues/40
            'confidence.py $TARGET',
            'ani.py $TARGET ${SOURCES[5:7]}',
            'seqhash.py $TARGET ${SOURCES[7]}',
            'sort_values.py $TARGET %(sort_by)s' % settings]
            )

"""
METHODS: Records are paritioned into Blast databases with the following sequence and
annotaiton filtering.  Record sequences must have a minimum sequence length of
1200 base pairs, a 1% maximimum proportion of ambiguous bases and an
alignment to the 16S coveriance model with an E-value of 0.01 or higher.
Record annotations must bee of species rank or below with a Taxtastic
classified species level name.  Sequences are grouped by accession and
deduplicated.  Following initial filtering a partion of type strain
records and a parition of all records group by species and
capped at a maxium of 5,000 representatives is created.
"""
fa, seq_info = env.Command(
    target=['$out/seqs.fasta', '$out/seq_info.csv'],
    source=[fa, feather],
    action='partition_refs.py --drop-noaligns $SOURCES $TARGETS')

raw = fa

"""
A type strain database is partitioned
"""
type_fa, type_info = env.Command(
    target=['$out/dedup/1200bp/types/seqs.fasta',
            '$out/dedup/1200bp/types/seq_info.csv'],
    source=[fa, feather],
    action='partition_refs.py '
           '--drop-duplicate-sequences '
           '--is_species '
           '--is_type '
           '--is_valid '
           '--min-length 1200 '
           '--prop-ambig-cutoff 0.01 '
           '$SOURCES $TARGETS')

type_tax, type_lineages, types_mothur, type_blast = taxonomy(
    type_fa, type_info, '$out/dedup/1200bp/types/')

"""
filter into named set and other criteria
"""
fa, seq_info = env.Command(
    target=['$out/dedup/1200bp/named/seqs.fasta',
            '$out/dedup/1200bp/named/seq_info.csv'],
    source=[fa, feather],
    action='partition_refs.py '
           '--drop-duplicate-sequences '
           '--is_species '
           '--is_valid '
           '--min-length 1200 '
           '--prop-ambig-cutoff 0.01 '
           '--species-cap %(species_cap)s '
           '${SOURCES[:2]} $TARGETS' % settings)

named = fa


"""

"""
taxtable, lineages, mothur, blast = taxonomy(
    fa, seq_info, '$out/dedup/1200bp/named/')

"""
Trimmed seqname,tax_id map file that can be easily cached by filter_outliers

METHODS:
"""
taxid_map = env.Command(
    target='$out/dedup/1200bp/named/tax_id_map.csv',
    source=seq_info,
    action='csvcut.py --columns seqname,tax_id --out $TARGET $SOURCE')

"""
Filter sequences. Use --threads if you need to to limit the number
of processes - otherwise deenurp will use all of them!

METHODS: Sequence outliers are determined using the deenurp filter_outliers
outlier detection algorithm.
"""
fa, details = env.Command(
    source=[fa, taxid_map, taxtable, cfiles['outliers_cache']],
    target=['$out/dedup/1200bp/named/filtered/unsorted.fasta',
            '$out/dedup/1200bp/named/filtered/outliers.csv'],
    action=['deenurp -vvv filter_outliers '
            '--cluster-type single '
            '--detailed-seqinfo ${TARGETS[1]} '
            '--distance-percentile 90.0 '
            '--filter-rank species '
            '--jobs 1 '
            '--log $log '
            '--max-distance 0.02 '
            '--min-distance 0.01 '
            '--min-seqs-for-filtering 5 '
            '--output-seqs ${TARGETS[0]}  '
            '--previous-details ${SOURCES[3]} '
            '--strategy cluster '
            '--threads-per-job 14 '
            '${SOURCES[:3]}',
            'cp ${TARGETS[1]} ' + cfiles['outliers_cache']])

"""
add distance metrics to feather file
"""
filter_outliers = env.Command(
    target='$out/.feather/filter_outliers.md5',
    source=[feather, details],
    action=['filter_outliers.py $SOURCES', 'md5sum ${SOURCES[0]} > $TARGET'])

"""
METHODS: After species outliers are determined a partition of records is
created with species outlier records discarded.
"""
taxtable, filtered_lineages, mothur, blast = taxonomy(
    fa, seq_info, '$out/dedup/1200bp/named/filtered/')

"""
METHODS: From the species outlier filtered partition a partition of type
strain records is created.
"""
filtered_type_fa, filtered_type_info = env.Command(
    target=['$out/dedup/1200bp/named/filtered/types/seqs.fasta',
            '$out/dedup/1200bp/named/filtered/types/seq_info.csv'],
    source=[fa, feather],
    action='partition_refs.py --is_type $SOURCES $TARGETS')

filtered_type_tax, filtered_type_lineages, filtered_type_mothur, _ = taxonomy(
    filtered_type_fa,
    filtered_type_info,
    '$out/dedup/1200bp/named/filtered/types/')

"""
METHODS:  All named records are parwise aligned with the outlier filtered
type strain database using vsearch to determine the top five closest type
strain hits.
"""
filtered_type_hits = env.Command(
    target='$out/dedup/1200bp/named/types_vsearch.tsv',
    source=[named, filtered_type_fa],
    action='vsearch --usearch_global ${SOURCES[0]} '
           '--blast6out $TARGET '
           '--db ${SOURCES[1]} '
           '--id 0.75 '
           '--maxaccepts 5 '
           '--self '  # reject same sequence hits
           '--strand plus '
           '--threads 14 '
           '--top_hits_only')

"""
This output will be used in the filter plots

METHODS: Using the top five type strain hits a classification is generated
and used to determine taxonomy annotation accuracy.
"""
filtered_type_classifications = env.Command(
    target='$out/dedup/1200bp/named/classifications.csv',
    source=[filtered_type_hits,
            filtered_type_info,
            filtered_type_tax,
            os.path.join('$pipeline', 'data/classifier_thresholds.csv')],
    action='classify -vv '
           '--lineages ${SOURCES[2]} '
           '--rank-thresholds ${SOURCES[3]} '
           '--seq-info ${SOURCES[1]} '
           '--starred 101 '
           '--out $TARGET '
           '${SOURCES[0]}')

"""
Adds type_classification to feather file for Dash application

METHODS:
"""
type_classifications = env.Command(
    target='$out/.feather/type_classifications.md5',
    source=[feather, filtered_type_classifications],
    action=['type_classifications.py $SOURCES',
            'md5sum ${SOURCES[0]} > $TARGET'])

"""
expand taxids into descendants

METHODS:
"""
trusted_taxids = env.Command(
    target='$out/dedup/1200bp/named/filtered/trusted/taxids.txt',
    source=settings['trust'],
    action='taxit get_descendants --out $TARGET $tax_url $SOURCE')

"""
expand taxids into descendants
"""
dnt_ids = env.Command(
    target='$out/dedup/1200bp/named/filtered/trusted/dnt_ids.txt',
    source=settings['do_not_trust'],
    action='taxit get_descendants --out $TARGET $tax_url $SOURCE')

"""
"""
trusted = env.Command(
    target='$out/dedup/1200bp/named/filtered/trusted/trust_ids.txt',
    source=[settings['trust'], trusted_taxids],
    action='cat $SOURCES > $TARGET')

"""
"""
dnt = env.Command(
    target='$out/dedup/1200bp/named/filtered/trusted/dnt.txt',
    source=[settings['do_not_trust'], dnt_ids],
    action='cat $SOURCES > $TARGET')

"""
Same as named set with inliers and trust/no_trust records
"""
fa, seq_info = env.Command(
    target=['$out/dedup/1200bp/named/filtered/trusted/seqs.fasta',
            '$out/dedup/1200bp/named/filtered/trusted/seq_info.csv'],
    source=[named, feather, trusted, dnt],
    action='partition_refs.py '
           '--do_not_trust ${SOURCES[3]} '
           '--drop-duplicate-sequences '
           '--inliers '  # filter_outliers = True & is_out = False
           '--is_species '
           '--is_valid '
           '--min-length 1200 '
           '--prop-ambig-cutoff 0.01 '
           '--species-cap %(species_cap)s '
           '--trusted ${SOURCES[2]} '
           '${SOURCES[:2]} $TARGETS' % settings)
Depends([fa, seq_info], filter_outliers)

taxtable, lineages, mothur, blast = taxonomy(
    fa, seq_info, '$out/dedup/1200bp/named/filtered/trusted/')

'''
Pull type strains from trusted db
'''
fa, seq_info = env.Command(
    target=['$out/dedup/1200bp/named/filtered/trusted/types/seqs.fasta',
            '$out/dedup/1200bp/named/filtered/trusted/types/seq_info.csv'],
    source=[fa, feather],
    action='partition_refs.py --is_type $SOURCES $TARGETS')


taxtable, lineages, mothur, blast = taxonomy(
    fa, seq_info, '$out/dedup/1200bp/named/filtered/trusted/types/')

'''
find top hit for each sequence among type strains

NOTE: alleles will never align with themselves (--self) BUT
can align with other alleles in the same genome accession

FIXME: use named types not trusted types

METHODS:
'''
named_type_hits = env.Command(
    target='$out/dedup/1200bp/named/trusted_vsearch.tsv',
    source=[named, fa],
    action='vsearch --usearch_global ${SOURCES[0]} '
           '--blast6out $TARGET '
           '--db ${SOURCES[1]} '
           '--id 0.75 '
           '--maxaccepts 5 '
           '--self '  # reject same sequence hits
           '--strand plus '
           '--threads 14 '
           '--top_hits_only')

"""
Creates match_seqname, match_pct, match_version, match_species and
match_species_id columns for best type strain hits from the trusted
dataset

METHODS:
"""
match_hits = env.Command(
    target='$out/.feather/match_hits.md5',
    source=[feather, named_type_hits],
    action=['match_hits.py $SOURCES', 'md5sum ${SOURCES[0]} > $TARGET'])

"""
gzip feather file as last step for the filter_outlier plots

METHODS:
"""
gzip = env.Command(
    target='${SOURCE}.gz',
    source=feather,
    action='gzip --to-stdout $SOURCE > $TARGET')
Depends([gzip], match_hits)

"""
TODO: create a RUNNING and SUCCCESS file
https://scons.org/doc/2.0.1/HTML/scons-user/x2045.html
"""

"""
copy taxdmp file into output dir so a ``taxit new_database``
can be built again in the future if needed

METHODS:
"""
taxdmp = env.Command(
    source=settings['taxdmp'],
    target='$out/' + os.path.basename(settings['taxdmp']),
    action='cp $SOURCE $TARGET')

"""
git version used to generate output
"""
commit = env.Command(
    target='$out/git_version.txt',
    source=None,
    action='(echo $$(hostname):$pipeline;version.py) > $TARGET')


def write_build_status():
    """
    Create a file SUCCESS or FAILED depending on how the pipeline finishes
    """
    success = os.path.join(out, 'SUCCESS')
    failed = os.path.join(out, 'FAILED')
    if os.path.isfile(success):
        os.remove(success)
    if os.path.isfile(failed):
        os.remove(failed)
    if GetBuildFailures():
        status = failed
    else:
        status = success
    try:
        open(status, 'w').close()
    except FileNotFoundError:
        # scons --dry-run
        pass


atexit.register(write_build_status)

"""
release steps
"""
if release:
    def SymLink(directory='.'):
        '''
        scons does not manage files outside its base directory so we work
        around that by passing the symlink directory as a separate argument not
        managed by scons
        '''
        def SymLinkAction(target, source, env):
            src = os.path.abspath(str(source[0]))
            targ = os.path.abspath(os.path.join(directory, str(target[0])))
            try:
                os.symlink(src, targ)
            except OSError as e:
                if e.errno == errno.EEXIST:
                    os.remove(targ)
                    os.symlink(src, targ)
                else:
                    raise e

        return SymLinkAction

    '''
    create symbolic link LATEST on directory up to point to $out dir
    '''
    latest = env.Command(
        target=os.path.join('$outdir', 'LATEST'),
        source='$out',
        action=SymLink())

    if venv:
        freeze = env.Command(
            target='$out/requirements.txt',
            source=venv,
            action='pip freeze > $TARGET')
