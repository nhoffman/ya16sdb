"""
Download and curate the NCBI 16S rRNA sequences
"""
import configparser
import csv
import errno
import os
import sys
import SCons
import time

from SCons.Script import (
        ARGUMENTS, Depends, Environment, Help, PathVariable, Variables)

venv = os.environ.get('VIRTUAL_ENV')
if not venv:
    sys.exit('--> an active virtualenv is required'.format(venv))
if not os.path.exists('settings.conf'):
    sys.exit("Can't find settings.conf")
conf = configparser.SafeConfigParser()
conf.read('settings.conf')
settings = conf['DEFAULT']


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
        action='makeblastdb -dbtype nucl -in $SOURCE -out ' + output_base)
    env.Command(
        target=output_base,
        source=blast_out,
        action='md5sum $SOURCES > $TARGET')
    return blast_out


true_vals = ['t', 'y', '1']
release = ARGUMENTS.get('release', 'no').lower()[0] in true_vals
test = ARGUMENTS.get('test', 'no').lower()[0] in true_vals
vrs = Variables(None, ARGUMENTS)
out = 'test_output' if test else 'output'
vrs.Add('base', help='Path to output directory', default=out)
vrs.Add('out', default=os.path.join('$base', time.strftime('%Y%m%d')))
vrs.Add('tax_url', default=settings['taxonomy'], help='database url')
# cache vars
vrs.Add(PathVariable(
    'genbank_cache', '', '$base/records.gb', PathIsFileCreate))
vrs.Add(PathVariable(
    'outliers_cache', '', '$base/filter_outliers.csv', PathIsFileCreate))
vrs.Add(PathVariable(
    'pubmed_info_cache', '', '$base/pubmed_info.csv', PathIsFileCreate))
vrs.Add(PathVariable(
    'seqs_cache', '', '$base/seqs.fasta', PathIsFileCreate))
vrs.Add(PathVariable(
    'seq_info_cache', '', '$base/seq_info.csv', PathIsFileCreate))
vrs.Add(PathVariable(
    'records_cache', '', '$base/records.txt', PathIsFileCreate))
vrs.Add(PathVariable(
    'references_cache', '', '$base/references.csv', PathIsFileCreate))
vrs.Add(PathVariable(
    'refseq_info_cache', '', '$base/refseq_info.csv', PathIsFileCreate))
vrs.Add(PathVariable(
    'unknown_cache', '', '$base/unknown.txt', PathIsFileCreate))

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
    eutils=settings['eutils'],
    taxit=settings['taxit'],
    deenurp=settings['deenurp'],
)

env.Decider('MD5')

Help(vrs.GenerateHelpText(env))

"""
get accessions (versions) of records considered type strains
NCBI web blast uses `sequence_from_type[Filter]` so we will use that
http://www.ncbi.nlm.nih.gov/news/01-21-2014-sequence-by-type/
"""
types = env.Command(
    source=None,
    target='$out/ncbi/types.txt',
    action='$eutils %(esearch)s "%(types)s" | %(acc)s' % settings)

if test:
    taxids = (i.strip() for i in open('testfiles/tax_ids.txt'))
    taxids = ('txid' + i + '[Organism]' for i in taxids if i)
    taxids = ' OR '.join(taxids)
    ncbi = env.Command(
        target='$out/ncbi/records.txt',
        source='testfiles/tax_ids.txt',
        action=('$eutils %(esearch)s "%(16s)s AND (%(taxids)s)" | %(acc)s' %
                {'taxids': taxids, **settings}))
else:
    tm7 = env.Command(
        source=None,
        target='$out/ncbi/tm7/versions.txt',
        action='$eutils %(esearch)s "%(tm7)s" | %(acc)s' % settings)

    classified = env.Command(
        source=None,
        target='$out/ncbi/classified.txt',
        action='$eutils %(esearch)s "%(classified)s" | %(acc)s' % settings)

    '''
    We will use this records.txt list to remove any records that we had
    downloaded previously
    '''
    ncbi = env.Command(
        source=[classified, tm7],
        target='$out/ncbi/records.txt',
        action='cat $SOURCES > $TARGET')

"""
Do not download record accessions in the ignore list or that have been
previously downloaded in the records_cache or in unknown_cache.
Exit script if no new records exist.

TODO: Do not exit, just continue
"""
new = env.Command(
    target='$out/ncbi/new.txt',
    source=[ncbi, settings['ignore'], '$records_cache', '$unknown_cache'],
    action=['cat ${SOURCES[1:]} | '
            'grep '
            '--invert-match '
            '--fixed-strings '
            '--file /dev/stdin '
            '${SOURCES[0]} > $TARGET '
            '|| (echo "No new records" && exit 1)'])

"""
Check the cache for last download_date and download list of modified records
in order to redownload modified records that we have previously downloaded.
"""
seq_info = csv.DictReader(open(env.subst('$seq_info_cache')))
seq_info = [time.strptime(r['download_date'], '%d-%b-%Y') for r in seq_info]
if seq_info:
    download_date = time.strftime('%Y/%m/%d', sorted(seq_info)[-1])
else:
    download_date = time.strftime('%Y/%m/%d')
modified = env.Command(
    source=None,
    target='$out/ncbi/modified.txt',
    action=['$eutils %(esearch)s "%(classified)s '
            'AND %(download_date)s[Modification Date] : '
            '3000[Modification Date]" | %(acc)s' %
            {'download_date': download_date, **settings}])

"""
`sort --unique` so we are not downloading records twice
"""
download = env.Command(
    target='$out/ncbi/download.txt',
    source=[new, modified],
    action='cat $SOURCES | sort --unique > $TARGET')

"""
download genbank records
"""
gbs = env.Command(
    target='$out/ncbi/records.gb',
    source=download,
    action='%(fts)s | %(ftract)s | %(gbs)s' % settings)

"""
extract
"""
today = time.strftime('%d-%b-%Y')
fa, seq_info, pubmed_info, references, refseq_info = env.Command(
    target=['$out/ncbi/seqs.fasta',
            '$out/ncbi/seq_info.csv',
            '$out/ncbi/pubmed_info.csv',
            '$out/ncbi/references.csv',
            '$out/ncbi/refseq_info.csv'],
    source=gbs,
    action='extract_genbank.py $SOURCE ' + today + ' $TARGETS')

"""
Record versions returned from esearch that had no actual 16s features
"""
no_features = env.Command(
    target='$out/ncbi/no_features.txt',
    source=[seq_info, download],
    action=('csvcut.py --columns version ${SOURCES[0]} | '
            'tail -n +2 | '
            'grep '
            '--invert-match '
            '--fixed-strings '
            '--file /dev/stdin '
            '${SOURCES[1]} > $TARGET '
            # force continue if grep finds no matches
            '|| true'))

"""
filter unnamed tax_ids

Do nothing with the unknown records because it might simply mean
the ncbi taxonomy pipeline is out of sync with the latest 16s records

TODO: replace this with new accession2taxid script
"""
seq_info, _ = env.Command(
    target=['$out/ncbi/taxit/seq_info.csv',
            '$out/ncbi/taxit/unknown.csv'],
    source=seq_info,
    action=['$taxit update_taxids '
            '--unknown-action drop '
            '--unknowns ${TARGETS[1]} '
            '--outfile ${TARGETS[0]} '
            '$SOURCE $tax_url'])

"""
vsearch new sequences with training set to test sequence orientation
and 16s region
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
            '--userout $TARGET'))

"""
Fix record orientation and ignore sequences with no vsearch alignments.

NOTE: unknown.txt will contain records (accession.version) ids with ANY
16s coordinates so we can potentially re-download these
coordinates later when data/rdp_16s_type_strains.fasta.bz2 is updated.

TODO: concat previous $unknown_cache file
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
Append with older records

1. Drop seqnames missing either a sequence or row in seq_info, sequences
   filtered out of the vsearch 16s alignment or sequences with unknown tax_ids
2. Append seqs, seq_info, pubmed_info and references to previous data set
3. Drop records not in records.txt file
4. Drop sequences that have a refseq equivalent
5. Deduplicate pubmeds and references
6. Copy and cache full dataset

NOTE: seqs that failed either the vsearch or taxit update_taxids will
not be appended to the records.txt file and will therefore be re-downloaded
the next time this pipeline is run
"""
fa, seq_info, pubmed_info, _, refseq_info, _ = env.Command(
    target=['$out/ncbi/combined/seqs.fasta',
            '$out/ncbi/combined/seq_info.csv',
            '$out/pubmed_info.csv',
            '$out/references.csv',
            '$out/refseq_info.csv',
            '$out/records.txt'],
    source=[ncbi,
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
update tax_ids
"""
seq_info = env.Command(
    target='$out/ncbi/combined/update_taxids/seq_info.csv',
    source=seq_info,
    action='$taxit -v update_taxids --out $TARGET $SOURCE $tax_url')

"""
taxonomy
"""
taxonomy = env.Command(
    target='$out/taxonomy.csv',
    source=seq_info,
    action='$taxit -v taxtable --seq-info $SOURCE --out $TARGET $tax_url')

"""
feather file
"""
feather = env.Command(
    target='$out/seq_info.feather',
    source=seq_info,
    action='to_feather.py $SOURCE $TARGET')

"""Start feather file columns"""
"""
'species' taxonomy id column
"""
species = env.Command(
    target='$out/.feather/species.md5',
    source=[feather, taxonomy],
    action=['species.py $SOURCES', 'md5sum ${SOURCES[0]} > $TARGET'])

"""
https://github.com/nhoffman/ya16sdb/issues/11
"""
is_type = env.Command(
    target='$out/.feather/is_type.md5',
    source=[feather, types],
    action=['is_type.py $SOURCES', 'md5sum ${SOURCES[0]} > $TARGET'])

is_published = env.Command(
    target='$out/.feather/is_publshed.md5',
    source=[feather, pubmed_info],
    action=['is_published.py $SOURCES', 'md5sum ${SOURCES[0]} > $TARGET'])

"""
add is_refseq and original column with refseq accession
"""
is_refseq = env.Command(
    target='$out/.feather/is_refseq.md5',
    source=[feather, refseq_info],
    action=['is_refseq.py $SOURCES', 'md5sum ${SOURCES[0]} > $TARGET'])

'''
is_valid attribute from taxtastic

https://github.com/fhcrc/taxtastic/blob/master/taxtastic/ncbi.py#L172
'''
is_valid = env.Command(
    target='$out/.feather/is_valid.md5',
    source=feather,
    action=['is_valid.py $SOURCE $tax_url', 'md5sum $SOURCE > $TARGET'])

"""
https://gitlab.labmed.uw.edu/uwlabmed/mkrefpkg/issues/40
"""
confidence = env.Command(
    target='$out/.feather/confidence.md5',
    source=feather,
    action=['confidence.py $SOURCE', 'md5sum $SOURCE > $TARGET'])
Depends(confidence, [is_type, is_refseq, is_published])

"""
calculate md5hash of sequences
"""
seqhash = env.Command(
    target='$out/.feather/seqhash.md5',
    source=[feather, fa],
    action=['seqhash.py $SOURCES', 'md5sum ${SOURCES[0]} > $TARGET'])

"""End feather columns"""

sort_values = env.Command(
    target='$out/.feather/sorted.md5',
    source=feather,
    action=['sort_values.py $SOURCE %(sort_by)s' % settings,
            'md5sum $SOURCE > $TARGET'])
Depends(sort_values, [is_type, is_refseq, is_published, seqhash])

fa, seq_info = env.Command(
    target=['$out/seqs.fasta', '$out/seq_info.csv'],
    source=[fa, feather],
    action='partition_refs.py $SOURCES $TARGETS')
Depends([fa, seq_info], sort_values)

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
filter using various criteria
"""
fa, info = env.Command(
    target=['$out/dedup/1200bp/named/seqs.fasta',
            '$out/dedup/1200bp/named/seq_info.csv'],
    source=[fa, feather],
    action=('partition_refs.py '
            '--drop-duplicate-sequences '
            '--is_species '
            '--is_valid '
            '--min-length 1200 '
            '--prop-ambig-cutoff 0.01 '
            '--species-cap 5000 '
            '${SOURCES[:2]} $TARGETS'))
Depends([fa, seq_info], [is_valid, species])

named = fa

"""
Make general named taxtable with all ranks included for filter_outliers
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
update tax_ids in details_in cache
"""
details_cache = env.Command(
    source='$outliers_cache',
    target='$out/dedup/1200bp/named/filtered/details_in.csv',
    action=['$taxit update_taxids '
            '--outfile $TARGET '
            '$SOURCE $tax_url'
            # continue if filter_outliers cache is empty
            ' || echo "Continuing without filter_outliers cache"'])

"""
Filter sequences. Use --threads if you need to to limit the number
of processes - otherwise deenurp will use all of them!
"""
fa, details = env.Command(
    source=[fa, taxid_map, taxonomy, details_cache],
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
Create a filtered seq_info.csv file
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
get taxid descendants from any taxids in trust.txt
"""
trusted_taxids = env.Command(
    target='$out/dedup/1200bp/named/filtered/trusted/taxids.txt',
    source=settings['trust'],
    action='$taxit get_descendants --out $TARGET $tax_url $SOURCE')

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
labmed trust/do_not_trust seqs
"""
fa, seq_info = env.Command(
    target=['$out/dedup/1200bp/named/filtered/trusted/seqs.fasta',
            '$out/dedup/1200bp/named/filtered/trusted/seq_info.csv'],
    source=[raw, feather, trusted, dnt],
    action=['partition_refs.py '
            '--do_not_trust ${SOURCES[3]} '
            '--drop-duplicate-sequences '
            '--inliers '  # column is_out = False
            '--is_species '
            '--is_valid '
            '--min-length 1200 '
            '--prop-ambig-cutoff 0.01 '
            '--trusted ${SOURCES[2]} '
            '--trusted-taxids ${SOURCES[2]} '
            '--species-cap 5000 '
            '${SOURCES[:2]} $TARGETS'])
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
            '--db ${SOURCES[1]} '
            '--blast6out $TARGET '
            '--id 0.75 '
            '--threads 14 '
            '--self '  # reject same sequence hits
            '--threads 12 '
            '--maxaccepts 1 '
            '--strand plus'))

"""
add named_type_hits match columns to feather file
"""
match_hits = env.Command(
    target='$out/.feather/match_hits.md5',
    source=[feather, named_type_hits],
    action=['match_hits.py $SOURCES', 'md5sum ${SOURCES[0]} > $TARGET'])

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
