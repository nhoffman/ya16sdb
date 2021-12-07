"""
Download and curate the 16S rRNA sequences from NCBI

NOTE: run with `scons -j 1` only
"""
import atexit
import configparser
import csv
import errno
import os
import SCons
import sys
import time

from SCons.Script import (
        ARGUMENTS,
        GetBuildFailures,
        Depends,
        Environment,
        Help,
        PathVariable,
        Variables,
        )

venv = os.environ.get('VIRTUAL_ENV')
if not venv:
    sys.exit('--> an active virtualenv is required')
if not os.path.exists('settings.conf'):
    sys.exit("Can't find settings.conf")


def PathIsFileCreate(key, val, env):
    """check if Path is a cache file and
    creating it if it does not exist."""
    if os.path.isdir(val):
        m = 'Path for option %s is a directory, not a file: %s'
        raise SCons.Errors.UserError(m % (key, val))
    if not os.path.exists(val):
        d = os.path.dirname(val)
        if not os.path.exists(d):
            os.makedirs(d)
        open(val, 'w').close()


def blast_db(env, sequence_file, output_base):
    '''
    Create a blast database and file md5sum
    '''
    extensions = ['.nhr', '.nin', '.nsq']
    blast_out = env.Command(
        target=[output_base + ext for ext in extensions],
        source=sequence_file,
        action=('$blast makeblastdb -dbtype nucl '
                '-in $SOURCE -out ' + output_base))
    env.Command(
        target=output_base,
        source=blast_out,
        action='md5sum $SOURCES > $TARGET')
    return blast_out


true_vals = ['t', 'y', '1']
release = ARGUMENTS.get('release', 'no').lower()[0] in true_vals
test = ARGUMENTS.get('test', 'no').lower()[0] in true_vals
conf = configparser.SafeConfigParser()
conf.read('settings.conf')
settings = conf['TEST'] if test else conf['DEFAULT']
out = os.path.join(settings['outdir'], time.strftime('%Y%m%d'))
vrs = Variables(None, ARGUMENTS)
vrs.Add('base', help='Path to output directory', default=settings['outdir'])
vrs.Add('cache', default=os.path.join('$base', '.cache'))
vrs.Add('out', default=out)
vrs.Add('tax_url', default=settings['taxonomy'], help='database url')
# cache vars
vrs.Add(PathVariable(
    'genbank_cache', '', '$cache/records.gb', PathIsFileCreate))
vrs.Add(PathVariable(
    'outliers_cache', '', '$cache/filter_outliers.csv', PathIsFileCreate))
vrs.Add(PathVariable(
    'pubmed_info_cache', '', '$cache/pubmed_info.csv', PathIsFileCreate))
vrs.Add(PathVariable(
    'seqs_cache', '', '$cache/seqs.fasta', PathIsFileCreate))
vrs.Add(PathVariable(
    'seq_info_cache', '', '$cache/seq_info.csv', PathIsFileCreate))
vrs.Add(PathVariable(
    'records_cache', '', '$cache/records.txt', PathIsFileCreate))
vrs.Add(PathVariable(
    'references_cache', '', '$cache/references.csv', PathIsFileCreate))
vrs.Add(PathVariable(
    'refseq_info_cache', '', '$cache/refseq_info.csv', PathIsFileCreate))
vrs.Add(PathVariable(
    'unknown_cache', '', '$cache/unknown.txt', PathIsFileCreate))

environment_variables = dict(
    os.environ,
    PATH=':'.join([
        'bin',
        os.path.join(venv, 'bin'),
        '/usr/local/bin',
        '/usr/bin',
        '/bin']),
)

env = Environment(
    ENV=environment_variables,
    variables=vrs,
    shell='bash',
    blast=settings['blast'],
    eutils=settings['eutils'],
    taxit=settings['taxit'],
    deenurp=settings['deenurp'],
)

env.EnsureSConsVersion(3, 0, 5)
env.Decider('MD5')

Help(vrs.GenerateHelpText(env))

classified = env.Command(
    source=None,
    target='$out/ncbi/classified.txt',
    action=['$eutils %(esearch)s "%(classified)s" | '
            '%(acc)s -out $TARGET' % settings])

"""
Candidatus Saccharibacteria
https://gitlab.labmed.uw.edu/molmicro/mkrefpkg/issues/36
"""
tm7 = env.Command(
    source=None,
    target='$out/ncbi/tm7.txt',
    action='$eutils %(esearch)s "%(tm7)s" | %(acc)s -out $TARGET' % settings)

"""
Check the cache for last download_date and download list of modified
records in order to re-download modified records that we have previously
downloaded.
"""
si = csv.DictReader(open(env.subst('$seq_info_cache')))
si = [time.strptime(r['download_date'], '%d-%b-%Y') for r in si]
if si:
    download_date = time.strftime('%Y/%m/%d', sorted(si)[-1])
else:
    download_date = time.strftime('%Y/%m/%d')
modified = env.Command(
    source=None,
    target='$out/ncbi/modified.txt',
    action=['$eutils %(esearch)s "(%(classified)s OR %(tm7)s) '
            'AND %(download_date)s[Modification Date] : '
            '3000[Modification Date]" | %(acc)s -out $TARGET' %
            {'download_date': download_date, **settings}])

"""
type strains records
NCBI web blast uses `sequence_from_type[Filter]` so we will use that
http://www.ncbi.nlm.nih.gov/news/01-21-2014-sequence-by-type/
"""
types = env.Command(
    source=None,
    target='$out/ncbi/types.txt',
    action='$eutils %(esearch)s "%(types)s" | %(acc)s -out $TARGET' % settings)

"""
Trim accession2taxid with 16s records and update taxids
"""
accession2taxid = env.Command(
    target='$out/ncbi/accession2taxid.csv',
    source=[settings['accession2taxid'], classified, tm7],
    action=['accession2taxid.py $SOURCES | '
            '$taxit update_taxids '
            '--outfile $TARGET '
            '--unknown-action drop '
            '- $tax_url'])

"""
FIXME: call this cache_in.py (cache_out.py will replace refresh.py)
"""
cache = env.Command(
    target='$out/ncbi/cache.txt',
    source=[accession2taxid, '$seq_info_cache', modified,
            '$records_cache', '$unknown_cache'],
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
download genbank records
"""
gbs = env.Command(
    target='$out/ncbi/download.gb',
    source=download,
    action=['%(fts)s -id $SOURCE | '
            '%(ftract)s | '
            '%(gbs)s -out $TARGET' % settings])

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

FIXME: get rid of the ugly grep
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
"""
fa, seq_info = env.Command(
    target=['$out/ncbi/extract_genbank/update_taxids/seqs.fasta',
            '$out/ncbi/extract_genbank/update_taxids/seq_info.csv'],
    source=[fa, seq_info],
    action=['$taxit update_taxids '
            '--unknown-action drop '
            '${SOURCES[1]} $tax_url | '
            'partition_refs.py ${SOURCES[0]} - $TARGETS'])

"""
vsearch new sequences with training set to test sequence orientation
and 16s region

FIXME: consider using cmsearch or generating a
single sequence from the rfam covariance model
"""
vsearch = env.Command(
    target='$out/ncbi/vsearch/vsearch.tsv',
    source=[fa, 'data/rdp_16s_type_strains.fasta.bz2'],
    action=('vsearch '
            '--db ${SOURCES[1]} '
            '--id 0.70 '
            '--iddef 2 '  # matching cols / alignment len excluding term gaps
            '--maxaccepts 1 '  # default is 1
            '--maxrejects 32 '  # default is 32
            '--mincols 350 '  # 500 (min len downloaded) * 0.70 (--query_col)
            '--output_no_hits '
            '--query_cov 0.70 '
            '--strand both '
            '--threads 14 '
            '--top_hits_only '
            '--usearch_global ${SOURCES[0]} '
            '--userfields query+target+qstrand+id+tilo+tihi '
            '--userout $TARGET'
            ' || true'))

"""
Fix record orientation and ignore sequences with no vsearch alignments.

NOTE: unknown.txt will contain records (accession.version) with ANY
filtered 16s allele.

FIXME: I think this unknown_cache should be appended not copied
"""
fa, seq_info, _, _ = env.Command(
    target=['$out/ncbi/vsearch/seqs.fasta',
            '$out/ncbi/vsearch/seq_info.csv',
            '$out/ncbi/vsearch/unknown.fasta',
            '$out/ncbi/vsearch/unknown.txt'],
    source=[vsearch, fa, seq_info, '$unknown_cache'],
    action=['vsearch.py $SOURCES $TARGETS',
            'cp ${TARGETS[3]} $unknown_cache'])

"""
Refresh/append with older records

1. Drop seqnames missing either a sequence or row in seq_info
2. Append seqs, seq_info, pubmed_info and references to previous data set
3. Drop records not in the ncbi records.txt file
4. Drop sequences that have a refseq equivalent
5. Deduplicate pubmeds and references
6. Copy and cache

FIXME: Call this cache_out.py ... with --cachedir argument
FIXME: Move refseq dedup to partition_refs.py --drop-duplicates-sequences
"""
fa, seq_info, pubmed_info, _, refseq_info, _ = env.Command(
    target=['$out/ncbi/seqs.fasta',
            '$out/ncbi/seq_info.csv',
            '$out/pubmed_info.csv',
            '$out/references.csv',
            '$out/refseq_info.csv',
            '$out/records.txt'],
    source=[accession2taxid,
            fa, '$seqs_cache',
            seq_info, '$seq_info_cache',
            pubmed_info, '$pubmed_info_cache',
            references, '$references_cache',
            refseq_info, '$refseq_info_cache',
            no_features, '$records_cache'],
    action=['refresh.py $SOURCES $TARGETS',
            # cache
            'cp ${TARGETS[0]} $seqs_cache',
            'cp ${TARGETS[1]} $seq_info_cache',
            'cp ${TARGETS[2]} $pubmed_info_cache',
            'cp ${TARGETS[3]} $references_cache',
            'cp ${TARGETS[4]} $refseq_info_cache',
            'cp ${TARGETS[5]} $records_cache'])

"""
append new records to global list

TODO: create bin/dedup_gb.py for genbank record cache maintenance
"""
env.Command(
    target=None,
    source=gbs,
    action='cat $SOURCE >> $genbank_cache')

"""
taxonomy
"""
taxonomy = env.Command(
    target='$out/taxonomy.csv',
    source=seq_info,
    action='$taxit -v taxtable --seq-info $SOURCE --out $TARGET $tax_url')

"""
Map for WGS records without a refseq assembly accession
"""
asm = env.Command(
    target='$out/ani/assembly_summary_genbank.txt',
    source=None,
    action=('wget '
            '--output-document $TARGET '
            '--retry-connrefused '
            'https://ftp.ncbi.nlm.nih.gov/'
            'genomes/genbank/assembly_summary_genbank.txt'))

"""
The ANI tax check report
"""
ani = env.Command(
    target='$out/ani/ANI_report_prokaryotes.txt',
    source=None,
    action=('wget '
            '--output-document $TARGET '
            '--retry-connrefused '
            'https://ftp.ncbi.nlm.nih.gov/'
            'genomes/ASSEMBLY_REPORTS/ANI_report_prokaryotes.txt'))

"""
Create feather file and initial columns
FIXME: Switch $SOURCE and $TARGET around for each script
"""
feather = env.Command(
    target='$out/seq_info.feather',
    source=[
        seq_info,
        taxonomy,
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

fa, seq_info = env.Command(
    target=['$out/seqs.fasta', '$out/seq_info.csv'],
    source=[fa, feather],
    action='partition_refs.py $SOURCES $TARGETS')

raw = fa

"""
pull type sequences
"""
type_fa, type_info = env.Command(
    target=['$out/dedup/1200bp/types/seqs.fasta',
            '$out/dedup/1200bp/types/seq_info.csv'],
    source=[fa, feather],
    action=['partition_refs.py '
            '--drop-duplicate-sequences '
            '--is_species '
            '--is_type '
            '--is_valid '
            '--min-length 1200 '
            '--prop-ambig-cutoff 0.01 '
            '$SOURCES $TARGETS'])

"""
Make sequence_from_type taxtable with all ranks included
"""
type_tax = env.Command(
    target='$out/dedup/1200bp/types/taxonomy.csv',
    source=type_info,
    action='$taxit -v taxtable --seq-info $SOURCE --out $TARGET $tax_url')

"""
Create taxtable output with replacing tax_ids with taxnames
"""
type_lineages = env.Command(
    target='$out/dedup/1200bp/types/lineages.csv',
    source=[type_tax, type_info],
    action='$taxit lineage_table --csv-table $TARGET $SOURCES')

"""
Create mothur output

https://mothur.org/wiki/Taxonomy_File
"""
types_mothur = env.Command(
    target='$out/dedup/1200bp/types/lineages.txt',
    source=[type_tax, type_info],
    action='$taxit lineage_table --taxonomy-table $TARGET $SOURCES')

blast_db(env, type_fa, '$out/dedup/1200bp/types/blast')

"""
filter into named set and other criteria
"""
fa, seq_info = env.Command(
    target=['$out/dedup/1200bp/named/seqs.fasta',
            '$out/dedup/1200bp/named/seq_info.csv'],
    source=[fa, feather],
    action=('partition_refs.py '
            '--drop-duplicate-sequences '
            '--is_species '
            '--is_valid '
            '--min-length 1200 '
            '--prop-ambig-cutoff 0.01 '
            '--species-cap %(species_cap)s '
            '${SOURCES[:2]} $TARGETS' % settings))

named = fa

"""
Make named taxtable for filter_outliers
"""
taxonomy = env.Command(
    target='$out/dedup/1200bp/named/taxonomy.csv',
    source=seq_info,
    action='$taxit -v taxtable --seq-info $SOURCE --out $TARGET $tax_url')

"""
Trimmed seqname,tax_id map file that can be easily cached by filter_outliers
"""
taxid_map = env.Command(
    target='$out/dedup/1200bp/named/tax_id_map.csv',
    source=seq_info,
    action='csvcut.py --columns seqname,tax_id --out $TARGET $SOURCE')

"""
Create taxtable output with replacing tax_ids with taxnames
"""
lineages = env.Command(
    target='$out/dedup/1200bp/named/lineages.csv',
    source=[taxonomy, seq_info],
    action='$taxit lineage_table --csv-table $TARGET $SOURCES')

"""
Create taxtable output with replacing tax_ids with taxnames

https://mothur.org/wiki/Taxonomy_File
"""
mothur = env.Command(
    target='$out/dedup/1200bp/named/lineages.txt',
    source=[taxonomy, seq_info],
    action='$taxit lineage_table --taxonomy-table $TARGET $SOURCES')

"""
Filter sequences. Use --threads if you need to to limit the number
of processes - otherwise deenurp will use all of them!
"""
fa, details = env.Command(
    source=[fa, taxid_map, taxonomy, '$outliers_cache'],
    target=['$out/dedup/1200bp/named/filtered/unsorted.fasta',
            '$out/dedup/1200bp/named/filtered/outliers.csv'],
    action=['$deenurp -vvv filter_outliers '
            '--cluster-type single '
            '--detailed-seqinfo ${TARGETS[1]} '
            '--distance-percentile 90.0 '
            '--filter-rank species '
            '--jobs 1 '
            '--log $out/dedup/1200bp/named/filtered/deenurp.log '
            '--max-distance 0.02 '
            '--min-distance 0.01 '
            '--min-seqs-for-filtering 5 '
            '--output-seqs ${TARGETS[0]}  '
            '--previous-details ${SOURCES[3]} '
            '--strategy cluster '
            '--threads-per-job 14 '
            '${SOURCES[:3]}',
            'cp ${TARGETS[1]} $outliers_cache'])

"""
add distance metrics to feather file
"""
filter_outliers = env.Command(
    target='$out/.feather/filter_outliers.md5',
    source=[feather, details],
    action=['filter_outliers.py $SOURCES', 'md5sum ${SOURCES[0]} > $TARGET'])

"""
Create a filtered seq_info.csv file of filter_outliers fasta
"""
fa, seq_info = env.Command(
    target=['$out/dedup/1200bp/named/filtered/seqs.fasta',
            '$out/dedup/1200bp/named/filtered/seq_info.csv'],
    source=[fa, feather],
    action='partition_refs.py $SOURCES $TARGETS')

"""
Make general named taxtable with all ranks included for filter_outliers
"""
taxonomy = env.Command(
    target='$out/dedup/1200bp/named/filtered/taxonomy.csv',
    source=seq_info,
    action='$taxit -v taxtable --seq-info $SOURCE --out $TARGET $tax_url')

"""
Create taxtable output with replacing tax_ids with taxnames
"""
filtered_lineages = env.Command(
    target='$out/dedup/1200bp/named/filtered/lineages.csv',
    source=[taxonomy, seq_info],
    action='$taxit lineage_table --csv-table $TARGET $SOURCES')

"""
Create mothur output

https://mothur.org/wiki/Taxonomy_File
"""
mothur = env.Command(
    target='$out/dedup/1200bp/named/filtered/lineages.txt',
    source=[taxonomy, seq_info],
    action='$taxit lineage_table --taxonomy-table $TARGET $SOURCES')

"""
expand taxids into descendants
"""
trusted_taxids = env.Command(
    target='$out/dedup/1200bp/named/filtered/trusted/taxids.txt',
    source=settings['trust'],
    action='$taxit get_descendants --out $TARGET $tax_url $SOURCE')

"""
expand taxids into descendants
"""
dnt_ids = env.Command(
    target='$out/dedup/1200bp/named/filtered/trusted/dnt_ids.txt',
    source=settings['do_not_trust'],
    action='$taxit get_descendants --out $TARGET $tax_url $SOURCE')

trusted = env.Command(
    target='$out/dedup/1200bp/named/filtered/trusted/trust_ids.txt',
    source=[settings['trust'], trusted_taxids],
    action='cat $SOURCES > $TARGET')

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
    source=[raw, feather, trusted, dnt],
    action=['partition_refs.py '
            '--do_not_trust ${SOURCES[3]} '
            '--drop-duplicate-sequences '
            '--inliers '  # filter_outliers = True & is_out = False
            '--is_species '
            '--is_valid '
            '--min-length 1200 '
            '--prop-ambig-cutoff 0.01 '
            '--species-cap %(species_cap)s '
            '--trusted ${SOURCES[2]} '
            '${SOURCES[:2]} $TARGETS' % settings])
Depends([fa, seq_info], filter_outliers)

"""
Make filtered taxtable with ranked columns and no_rank rows
"""
taxonomy = env.Command(
    target='$out/dedup/1200bp/named/filtered/trusted/taxonomy.csv',
    source=seq_info,
    action='$taxit taxtable --seq-info $SOURCE --out $TARGET $tax_url')

"""
Taxtable output replacing tax_ids with taxnames
"""
lineages = env.Command(
    target='$out/dedup/1200bp/named/filtered/trusted/lineages.csv',
    source=[taxonomy, seq_info],
    action='$taxit lineage_table --csv-table $TARGET $SOURCES')

"""
Mothur output - https://mothur.org/wiki/Taxonomy_File
"""
mothur = env.Command(
    target='$out/dedup/1200bp/named/filtered/trusted/lineages.txt',
    source=[taxonomy, seq_info],
    action='$taxit lineage_table --taxonomy-table $TARGET $SOURCES')

blast_db(env, fa, '$out/dedup/1200bp/named/filtered/trusted/blast')

'''
Pull type strains from trusted db
'''
fa, seq_info = env.Command(
    target=['$out/dedup/1200bp/named/filtered/trusted/types/seqs.fasta',
            '$out/dedup/1200bp/named/filtered/trusted/types/seq_info.csv'],
    source=[fa, feather],
    action='partition_refs.py --is_type $SOURCES $TARGETS')

"""
Make trusted type taxtable with ranked columns and no_rank rows
"""
taxonomy = env.Command(
    target='$out/dedup/1200bp/named/filtered/trusted/types/taxonomy.csv',
    source=seq_info,
    action=('$taxit taxtable '
            '--seq-info $SOURCE '
            '--out $TARGET '
            '$tax_url'))

"""
Create taxtable output with replacing tax_ids with taxnames
"""
lineages = env.Command(
    target='$out/dedup/1200bp/named/filtered/trusted/types/lineages.csv',
    source=[taxonomy, seq_info],
    action='$taxit lineage_table --csv-table $TARGET $SOURCES')

"""
Mothur output

https://mothur.org/wiki/Taxonomy_File
"""
mothur = env.Command(
    target='$out/dedup/1200bp/named/filtered/trusted/types/lineages.txt',
    source=[taxonomy, seq_info],
    action='$taxit lineage_table --taxonomy-table $TARGET $SOURCES')

blast_db(env, fa, '$out/dedup/1200bp/named/filtered/trusted/types/blast')

'''
find top hit for each sequence among type strains

NOTE: alleles will never align with themselves (--self) BUT
can align with other alleles in the same genome accession
'''
named_type_hits = env.Command(
    target='$out/dedup/1200bp/named/vsearch.tsv',
    source=[named, fa],
    action=('vsearch --usearch_global ${SOURCES[0]} '
            '--blast6out $TARGET '
            '--db ${SOURCES[1]} '
            '--id 0.75 '
            '--maxaccepts 5 '
            '--self '  # reject same sequence hits
            '--strand plus '
            '--threads 14 '
            '--top_hits_only'))

"""
add named_type_hits match columns to feather file
"""
match_hits = env.Command(
    target='$out/.feather/match_hits.md5',
    source=[feather, named_type_hits],
    action=['match_hits.py $SOURCES', 'md5sum ${SOURCES[0]} > $TARGET'])

"""
gzip feather file as last step for the filter_outlier plots
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
    source='.git/objects',
    action=['(echo $$(hostname):$$(pwd); '
            'git describe --tags --dirty) > $TARGET'])


def write_build_status():
    """
    Create a file SUCCESS or FAILED depending on how the pipeline finishes
    """
    if GetBuildFailures():
        status = 'FAILED'
    else:
        status = 'SUCCESS'
    try:
        open(os.path.join(out, status), 'w').close()
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
        target=os.path.join('$base', 'LATEST'),
        source='$out',
        action=SymLink())

    freeze = env.Command(
        target='$out/requirements.txt',
        source=venv,
        action='pip freeze > $TARGET')
