"""
Download and curate the NCBI 16S rRNA sequences

TODO: download records updated before last download
"""

import errno
import os
import sys
import time

# requirements installed in the virtualenv
import SCons
from SCons.Script import (
    Variables, ARGUMENTS, Help, Copy, Environment, PathVariable)


def PathIsFileCreate(key, val, env):
    """Validator to check if Path is a cache file,
    creating it if it does not exist."""
    if os.path.isdir(val):
        m = 'Path for option %s is a directory, not a file: %s'
        raise SCons.Errors.UserError(m % (key, val))
    if not os.path.exists(val):
        d = os.path.dirname(val)
        if not os.path.exists(d):
            os.makedirs(d)
        open(val, 'w').close()


def blast_db(env, sequence_file, output_base, dbtype='nucl'):
    '''
    Create a blast database
    '''
    prefix = dbtype[0]
    extensions = ['.{0}{1}'.format(prefix, suffix)
                  for suffix in ('hr', 'in', 'sq')]

    blast_out = env.Command(
        target=[output_base + ext for ext in extensions],
        source=sequence_file,
        action=('makeblastdb -dbtype {0} '
                '-in $SOURCE '
                '-out {1}'.format(dbtype, output_base)))

    env.Command(
        target=output_base,
        source=blast_out,
        action=('md5sum $SOURCES > $TARGET'))

    return blast_out


venv = os.environ.get('VIRTUAL_ENV')
if not venv:
    sys.exit('--> an active virtualenv is required'.format(venv))

true_vals = ['t', 'y', '1']
release = ARGUMENTS.get('release', 'no').lower()[0] in true_vals
test = ARGUMENTS.get('test', 'no').lower()[0] in true_vals

vrs = Variables(None, ARGUMENTS)
if test:
    vrs.Add('base', help='Path to output directory', default='test_output')
else:
    vrs.Add('base', help='Path to output directory', default='output')
vrs.Add(
    'out',
    help='Path to dated output sub directory',
    default=os.path.join('$base', time.strftime('%Y%m%d')))
vrs.Add('email', 'email address for ncbi', 'crosenth@uw.edu')
vrs.Add('retry', 'ncbi retry milliseconds', '60000')
# nreq should be set to 3 during weekdays
vrs.Add('nreq', ('Number of concurrent http requests to ncbi'), 12)
vrs.Add(
    'tax_url',
    default='"postgresql://crosenth:password@db3.labmed.uw.edu/molmicro"',
    help='database url')
vrs.Add(
    'schema',
    default='ncbi_taxonomy',
    help='postgres database schema')

# cache vars
vrs.Add(PathVariable(
    'genbank_cache', '', '$base/records.gb', PathIsFileCreate))
vrs.Add(PathVariable(
    'seqs_cache', '', '$base/seqs.fasta', PathIsFileCreate))
vrs.Add(PathVariable(
    'annotations_cache', '', '$base/annotations.csv', PathIsFileCreate))
vrs.Add(PathVariable(
    'pubmed_info_cache', '', '$base/pubmed_info.csv', PathIsFileCreate))
vrs.Add(PathVariable(
    'references_cache', '', '$base/references.csv', PathIsFileCreate))
vrs.Add(PathVariable(
    'records_cache', '', '$base/records.txt', PathIsFileCreate))
vrs.Add(PathVariable(
    'refseq_info_cache', '', '$base/refseq_info.csv', PathIsFileCreate))
vrs.Add(PathVariable(
    'outliers_cache', '', '$base/filter_outliers.csv', PathIsFileCreate))

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
    taxit=(
        'singularity exec '
        '--bind $$(pwd$)):$$(pwd$)) '
        '--workdir $$(pwd$)) '
        '/molmicro/common/singularity/taxtastic-0.8.2.img taxit'),
    deenurp=(
        'singularity exec '
        '--bind $$(pwd$)):$$(pwd$)) '
        '--workdir $$(pwd$)) '
        '/molmicro/common/singularity/deenurp-v0.2.0.img deenurp')

)

env.Decider('MD5-timestamp')

Help(vrs.GenerateHelpText(env))

rrna_16s = ('16s[All Fields] '
            'AND rRNA[Feature Key] '
            'AND Bacteria[Organism] '
            'AND 500 : 99999999999[Sequence Length]')

mefetch_acc = ('mefetch -vv '
               '-email $email '
               '-mode text '
               '-format acc '
               '-retry $retry '
               '-proc $nreq '
               '-log $out/ncbi.log '
               '-out $TARGET')

"""
Download accessions for bacteria with 16s rrna annotations
NOT environmental or unclassified

TODO: Limit to a recent date frame as initial search.  If latest records
are not found then run the whole thing in full.  This will save time in
execution of this step.
"""
classified = env.Command(
    source=None,
    target='$out/esearch/classified.txt',
    action=('esearch -db nucleotide -query "' + rrna_16s +
            'NOT(environmental samples[Organism] '
            'OR unclassified Bacteria[Organism])" | ' + mefetch_acc))

"""
Download TM7 accessions
"""
tm7 = env.Command(
    source=None,
    target='$out/esearch/tm7.txt',
    action=('esearch -db nucleotide -query "' + rrna_16s +
            ' AND Candidatus Saccharibacteria[Organism]" | ' + mefetch_acc))

"""
get accessions (versions) of records considered type strains
NCBI web blast uses `sequence_from_type[Filter]` so we will use that
http://www.ncbi.nlm.nih.gov/news/01-21-2014-sequence-by-type/
"""
types = env.Command(
    source=None,
    target='$out/1200bp/valid/types/esearch.txt',
    action=('esearch -db nucleotide -query "' + rrna_16s +
            ' AND sequence_from_type[Filter]" | ' + mefetch_acc))

if test:
    tax_ids = (i.strip() for i in open('testfiles/tax_ids.txt') if i)
    tax_ids = ('txid' + i + '[Organism]' for i in tax_ids)
    records = env.Command(
        target='$out/esearch.txt',
        source='testfiles/tax_ids.txt',
        action=('esearch -db nucleotide -query "' + rrna_16s +
                ' AND (' + ' OR '.join(tax_ids) + ')" | ' + mefetch_acc))
else:
    """
    concat our download set
    """
    records = env.Command(
        source=[classified, tm7],
        target='$out/esearch.txt',
        action='cat $SOURCES > $TARGET')

"""
Do not download record accessions in the ignore list or
that have been previously downloaded in the records_cache.
Exit script if no new records exist.
"""
new = env.Command(
    target='$out/new/records.txt',
    source=[records, '$records_cache', 'data/ignore.txt'],
    action=['cat ${SOURCES[-2:]} | '
            'grep '
            '--invert-match '
            '--fixed-strings '
            '--file /dev/stdin '
            '${SOURCES[0]} > $TARGET '
            '|| (echo "No new records" && exit 1)'])

gbs = env.Command(
    target='$out/new/records.gb',
    source=new,
    action=['mefetch -vv '  # download feature tables
            '-email $email '
            '-retry $retry '
            '-id $SOURCE '
            '-db nucleotide '
            '-format ft '
            '-retmax 1 '
            '-log $out/ncbi.log '
            '-proc $nreq | '
            'ftract -feature rrna::16s | '   # extract 16s features
            'mefetch '  # download genbank records
            '-vv '
            '-email $email '
            '-retmax 1 '
            '-csv '
            '-db nucleotide '
            '-format gbwithparts '
            '-log $out/ncbi.log '
            '-proc $nreq '
            '-retry $retry > $TARGET'])

"""
extract
"""
today = time.strftime('%d-%b-%Y')
new_fa, new_annotations, new_pub_info, new_refs, new_refseq_info = env.Command(
    target=['$out/new/seqs.fasta',
            '$out/new/annotations.csv',
            '$out/new/pubmed_info.csv',
            '$out/new/references.csv',
            '$out/new/refseq_info.csv'],
    source=gbs,
    action='extract_genbank.py $SOURCE ' + today + ' $TARGETS')

"""
Record versions returned from esearch that had no actual 16s features
"""
no_features = env.Command(
    target='$out/new/no_features.txt',
    source=[new_annotations, new],
    action=('csvcut.py --columns version ${SOURCES[0]} | '
            'tail -n +2 | '
            'grep '
            '--invert-match '
            '--fixed-strings '
            '--file /dev/stdin '
            '${SOURCES[1]} > $TARGET '
            # avoid the grep no match exit code 1 that scons hates
            '|| true'))

"""
filter invalid tax_ids

Do nothing with the unknown records for now because it might simply mean
the ncbi taxonomy pipeline is out of sync with the latest 16s records
"""
known_info, _ = env.Command(
    target=['$out/new/taxit/annotations.csv', '$out/new/taxit/unknown.csv'],
    source=new_annotations,
    action=['$taxit update_taxids '
            '--unknown-action drop '
            '--unknowns ${TARGETS[1]} '
            '--outfile ${TARGETS[0]} '
            '--schema $schema '
            '$SOURCE $tax_url'])

"""
vsearch new sequences with training set to test sequence orientation
and 16s region
"""
vsearch = env.Command(
    target='$out/new/vsearch.csv',
    source=[new_fa, 'data/rdp_16s_type_strains.fasta.bz2'],
    action=('vsearch '
            '--usearch_global ${SOURCES[0]} '
            '--db ${SOURCES[1]} '
            '--id 0.70 '
            '--threads 14 '
            '--userfields query+target+qstrand+id+tilo+tihi '
            '--strand both '
            '--top_hits_only '
            '--output_no_hits '
            '--maxaccepts 1 '  # default is 1
            '--maxrejects 32 '  # default is 32
            '--userout $TARGET'))

"""
Fix record orientation.  Drop sequences with no alignments.
These records are rare and we want to be able to re-download later
them when the 16s training set updates.
"""
vsearch_fa, _ = env.Command(
    target=['$out/new/vsearch/seqs.fa',
            '$out/new/vsearch/unknown.fa'],
    source=[vsearch, new_fa],
    action='vsearch.py --unknowns ${TARGETS[1]} --out ${TARGETS[0]} $SOURCES')

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
fa, refresh_annotations, pubmed_info, references, refseq_info, _ = env.Command(
    target=['$out/seqs.fasta',
            '$out/refresh/annotations.csv',
            '$out/pubmed_info.csv',
            '$out/references.csv',
            '$out/refseq_info.csv',
            '$out/records.txt'],
    source=[records,
            vsearch_fa, '$seqs_cache',
            known_info, '$annotations_cache',
            new_pub_info, '$pubmed_info_cache',
            new_refs, '$references_cache',
            new_refseq_info, '$refseq_info_cache',
            no_features, '$records_cache'],
    action=['refresh.py $SOURCES $TARGETS',
            # cache
            Copy('$seqs_cache', '${TARGETS[0]}'),
            Copy('$pubmed_info_cache', '${TARGETS[2]}'),
            Copy('$references_cache', '${TARGETS[3]}'),
            Copy('$refseq_info_cache', '${TARGETS[4]}'),
            Copy('$records_cache', '${TARGETS[5]}')])

"""
append new records to global list

NOTE: see bin/dedup_gb.py for genbank record cache maintenance
TODO: create bin/dedup_gb.py
"""
env.Command(
    target=None,
    source=gbs,
    action='cat $SOURCE >> $genbank_cache')

"""
update all tax_ids
"""
annotations = env.Command(
    target='$out/annotations.csv',
    source=refresh_annotations,
    action=['$taxit -v update_taxids '
            '--out $TARGET '
            '--schema $schema '
            '$SOURCE $tax_url',
            Copy('$annotations_cache', '$TARGET')])

"""
pull sequences at least 1200 bp and less than 1% ambiguous
"""
full_annotations = env.Command(
    target='$out/1200bp/annotations.csv',
    source=annotations,
    action=('partition_refs.py '
            # filtering
            '--min-length 1200 '
            '--prop-ambig-cutoff 0.01 '
            '--out-annotations $TARGET '
            '$SOURCE'))

"""
filter for valid annotations
"""
valid_annotations = env.Command(
    target='$out/1200bp/valid/annotations.csv',
    source=full_annotations,
    action='is_valid.py --out $TARGET  --schema $schema $SOURCE $tax_url')

"""
pull valid sequences
"""
valid_fa = env.Command(
    target='$out/1200bp/valid/seqs.fasta',
    source=[fa, valid_annotations],
    action=('partition_refs.py '
            '--fasta ${SOURCES[0]} '
            '--out-fa $TARGET '
            '${SOURCES[1]}'))

"""
Make general valid taxtable with all ranks included
"""
valid_tax = env.Command(
    target='$out/1200bp/valid/taxonomy.csv',
    source=valid_annotations,
    action=('$taxit -v taxtable '
            '--seq-info $SOURCE '
            '--out $TARGET '
            '--schema $schema '
            '$tax_url'))

blast_db(env, valid_fa, '$out/1200bp/valid/blast')

"""
Remove and re-append is_type column with sequences per discussion:

https://gitlab.labmed.uw.edu/uwlabmed/mkrefpkg/issues/14
"""
type_annotations = env.Command(
    target='$out/1200bp/valid/types/annotations.csv',
    source=[valid_annotations, types],
    action='is_type.py $SOURCES $TARGET')

"""
pull type sequences
"""
type_fa = env.Command(
    target='$out/1200bp/valid/types/seqs.fasta',
    source=[valid_fa, type_annotations],
    action=('partition_refs.py '
            '--out-fa $TARGET '
            '--fasta ${SOURCES[0]} '
            '${SOURCES[1]}'))

"""
Deduplicate sequences by isolate (accession)

Appends a weight column showing number of sequences being represented
"""
dedup_fa, dedup_annotations = env.Command(
    target=['$out/1200bp/valid/dedup/seqs.fasta',
            '$out/1200bp/valid/dedup/annotations.csv'],
    source=[valid_fa, valid_annotations],
    action=('$deenurp -v deduplicate_sequences '
            '--group-by accession '
            '$SOURCES $TARGETS'))

"""
dedup info file for filter_outliers
"""
dedup_info = env.Command(
    target='$out/1200bp/valid/dedup/seq_info.csv',
    source=dedup_annotations,
    action='csvcut.py --columns seqname,tax_id --out $TARGET $SOURCE')

"""
update tax_ids in details_in cache
"""
filtered_details_in = env.Command(
    source='$outliers_cache',
    target='$out/1200bp/valid/dedup/filtered/details_in.csv',
    action=('$taxit update_taxids '
            '--schema $schema '
            '--outfile $TARGET '
            '$SOURCE $tax_url'
            # continue if filter_outliers cache is empty
            ' || true'))

"""
Filter sequences. Use --threads if you need to to limit the number
of processes - otherwise deenurp will use all of them!
"""
filtered_fa, filtered_info, filtered_details, deenurp_log = env.Command(
    source=[dedup_fa, dedup_info, valid_tax, filtered_details_in],
    target=['$out/1200bp/valid/dedup/filtered/seqs.fasta',
            '$out/1200bp/valid/dedup/filtered/seq_info.csv',
            '$out/1200bp/valid/dedup/filtered/details_out.csv',
            '$out/1200bp/valid/dedup/filtered/deenurp.log'],
    action=['$deenurp -vvv filter_outliers '
            '--log ${TARGETS[3]} '
            '--filter-rank species '
            '--threads-per-job 14 '
            '--jobs 1 '
            '--output-seqs ${TARGETS[0]}  '
            '--filtered-seqinfo ${TARGETS[1]} '
            '--detailed-seqinfo ${TARGETS[2]} '
            '--previous-details ${SOURCES[3]} '
            '--strategy cluster '
            '--cluster-type single '
            '--distance-percentile 90.0 '
            '--min-distance 0.01 '
            '--max-distance 0.02 '
            '--min-seqs-for-filtering 5 '
            '${SOURCES[:3]}',
            # cache it
            Copy('$outliers_cache', '${TARGETS[2]}')])

"""
Make filtered taxtable with ranked columns and no_rank rows
"""
filtered_tax = env.Command(
    target='$out/1200bp/valid/dedup/filtered/taxonomy.csv',
    source=filtered_info,
    action=('$taxit taxtable --seq-info $SOURCE --schema $schema $tax_url '
            '| ranked.py --columns --out $TARGET'))

blast_db(env, filtered_fa, '$out/1200bp/valid/dedup/filtered/blast')

'''
find top hit for each sequence among type strains
'''
type_hits = env.Command(
    target='$out/1200bp/valid/dedup/filtered/vsearch.blast6out',
    source=[dedup_fa, type_fa],
    action=('vsearch --usearch_global ${SOURCES[0]} '
            '--db ${SOURCES[1]} '
            '--blast6out $TARGET '
            '--id 0.75 '
            '--threads 14 '
            '--self '  # reject same sequence hits
            '--threads 12 '
            '--maxaccepts 1 '
            '--strand plus'))

'''
bokeh plot filtered sequences

hard coded: sort column 2 (records) desc
'''
env.Command(
    target=['$out/1200bp/valid/dedup/filtered/index.html',
            '$out/1200bp/valid/dedup/filtered/plots/map.csv'],
    source=[filtered_details,
            type_hits,
            annotations,
            valid_tax,
            types,
            deenurp_log],
    action=('plot_details.py ${SOURCES[:5]} '
            '--param strategy:cluster '
            '--param cluster_type:single '
            '--param distance_percentile:90.0 '
            '--param min_distance:0.01 '
            '--param max_distance:0.02 '
            '--log-in ${SOURCES[5]} '
            '--plot-dir $out/1200bp/valid/dedup/filtered/plots '
            '--plot-map ${TARGETS[1]} '
            '--plot-index ${TARGETS[0]}'))

"""
Append contributers
"""
contributors = env.Command(
    source='.git/logs/HEAD',
    target='contributors.txt',
    action='git log --all --format="%cN <%cE>" | sort | uniq > $TARGET')

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

    commit = env.Command(
        target='$out/git_version.txt',
        source='.git/objects',
        action='git describe --tags > $TARGET')

    freeze = env.Command(
        target='$out/requirements.txt',
        source=venv,
        action='pip freeze > $TARGET')
