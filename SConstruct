"""
Download and curate the NCBI 16S rRNA sequences

TODO: download records updated before last download
"""

import errno
import os
import sys
import time

# requirements installed in the virtualenv
from SCons.Script import (
    Variables, BoolVariable, ARGUMENTS, Help, Copy, Environment)

venv = os.environ.get('VIRTUAL_ENV')
if not venv:
    sys.exit('--> an active virtualenv is required'.format(venv))

true_vals = ['t', 'y', '1']
release = ARGUMENTS.get('release', 'no').lower()[0] in true_vals

vrs = Variables(None, ARGUMENTS)
vrs.Add('email', 'email address for ncbi', 'crosenth@uw.edu')
vrs.Add('retry', 'ncbi retry milliseconds', '60000')
vrs.Add('nproc', ('Number of concurrent processes '), 14)
# nreq should be set to 3 during weekdays
vrs.Add('nreq', ('Number of concurrent http requests to ncbi'), 12)
vrs.Add(BoolVariable('use_cluster', 'Use slurm', False))
vrs.Add('base', help='Path to base output directory', default='output-ncbi')
vrs.Add(
    'tax_url',
    default='"postgresql://crosenth:password@db3.labmed.uw.edu/molmicro"',
    help='database url')
vrs.Add(
    'schema',
    default='ncbi_taxonomy',
    help='postgres database schema')
vrs.Add(
    'out',
    help='Path to dated output directory',
    default=os.path.join('$base', time.strftime('%Y%m%d')))

# cache
vrs.Add('genbank_cache', default='$base/records.gb')
vrs.Add('seqs_cache', default='$base/seqs.fasta')
vrs.Add('annotations_cache', default='$base/seq_info.csv')
vrs.Add('tax_ids_cache', default='$base/tax_ids.csv')
vrs.Add('pubmed_info', default='$base/pubmed_info.csv')
vrs.Add('references_cache', default='$base/references.csv')
vrs.Add('versions_cache', default='$base/versions.txt')
vrs.Add('refseq_info_cache', default='$base/refseq_info.csv')

# Provides access to options prior to instantiation of env object
# below; it's better to access variables through the env object.
varargs = dict({opt.key: opt.default for opt in vrs.options}, **vrs.args)
truevals = {True, 'yes', 'y', 'True', 'true', 't'}
nproc = varargs['nproc']

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
    time=False
)

env.Decider('MD5')

Help(vrs.GenerateHelpText(env))


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

    env.Local(
        target=output_base,
        source=blast_out,
        action=('md5sum $SOURCES > $TARGET'))

    return blast_out


rrna_16s = ('16s[All Fields] '
            'AND rRNA[Feature Key] '
            'AND Bacteria[Organism] '
            'AND 500 : 99999999999[Sequence Length]')

mefetch_acc = ('mefetch '
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
    target='$out/classified.txt',
    action=('esearch -db nucleotide -query "' + rrna_16s +
            'NOT(environmental samples[Organism] '
            'OR unclassified Bacteria[Organism])" | ' + mefetch_acc))

"""
Download pubmed record accessions
"""
pubmed = env.Command(
    source=None,
    target='$out/pubmed.txt',
    action=('esearch -db nucleotide -query "' + rrna_16s +
            ' AND nuccore pubmed[Filter]" | ' + mefetch_acc))

"""
Download TM7 accessions
"""
tm7 = env.Command(
    source=None,
    target='$out/tm7.txt',
    action=('esearch -db nucleotide -query "' + rrna_16s +
            ' AND Candidatus Saccharibacteria[Organism]" | ' + mefetch_acc))

"""
get accessions (versions) of records considered type strains
NCBI web blast uses `sequence_from_type[Filter]` so we will use that
http://www.ncbi.nlm.nih.gov/news/01-21-2014-sequence-by-type/
"""
types = env.Command(
    source=None,
    target='$out/types.txt',
    action=('esearch -db nucleotide -query "' + rrna_16s +
            ' AND sequence_from_type[Filter]" | ' + mefetch_acc))

"""
concat our download set
"""
records = env.Command(
    source=[classified, tm7],
    target='$out/records.txt',
    action='cat $SOURCES > $TARGET')

"""
Do not download record accessions in the ignore list or
that have been previously downloaded in the versions_cache
"""
new = env.Command(
    target='$out/new/versions.txt',
    source=['data/text.txt', '$versions_cache', 'data/ignore.txt'],
    action=['cat ${SOURCES[-2]} | '
            'grep '
            '--invert-match '
            '--fixed-strings '
            '--file /dev/stdin '
            '${SOURCES[0]} > $TARGET '
            # avoid the grep no match exit code 1 that scons hates
            '|| true'])

gbs = env.Command(
    target='$out/new/records.gb',
    source=new,
    action=['mefetch '  # download feature tables
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
            '-retry $retry'])

"""
extract
"""
today = time.strftime('%d-%b-%Y')
new_fa, new_annotations, new_pub_info, new_refs, new_refseq_info = env.Command(
    target=['$out/new/seqs.fasta'
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
    target=['$out/new/taxit/annotations.csv',
            '$out/new/taxit/invalid.csv'],
    source=new_annotations,
    action=['taxit update_taxids '
            '--unknowns ${TARGETS[1]} '
            '--out ${TARGETS[0]} '
            '--schema $schema '
            '$SOURCES $tax_url'])

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
            '$out/new/vsearch/invalid.fa'],
    source=[vsearch, new_fa],
    action='vsearch.py --unknowns ${TARGETS[1]} --out ${TARGETS[0]} $SOURCES')

"""
1. Deduplicate pubmeds and references
2. Append seqs, seq_info, pubmed_ids and references to previous data set
3. Drop records not in records.txt file
4. Drop seqnames missing eigher a sequence or row in seq_info, sequences
   filtered out of the vsearch 16s alignment or sequences with unknown tax_ids
5. Drop sequences that have a refseq equivalent
6. Copy full dataset back to base dir for next download

NOTE: seqs that failed either the vsearch or taxit update_taxids will
not be appended to the versions.txt file and will therefore be re-downloaded
the next time this pipeline is run
"""
fa, refresh_annotations, pubmed_info, references, refseq_info, _ = env.Command(
    target=['$out/seqs.fasta',
            '$out/refresh/annotations.csv',
            '$out/pubmed_ids.csv',
            '$out/references.csv',
            '$out/refseq_info.csv',
            '$out/versions.csv'],
    source=[records,
            vsearch_fa, '$seqs_cache',
            known_info, '$annotations_cache',
            new_pub_info, '$pubmed_info',
            new_refs, '$references_cache',
            new_refseq_info, '$refseq_map_cache',
            no_features, '$versions_cache'],
    action=['refresh.py $SOURCES $TARGETS',
            # cache
            Copy('$seqs_cache', '${TARGETS[0]}'),
            Copy('$pubmed_info', '${TARGETS[2]}'),
            Copy('$references_cache', '${TARGETS[3]}'),
            Copy('$refseq_info', '${TARGETS[4]}'),
            Copy('$versions_cache', '$TARGETS[5]')])

"""
update all tax_ids
"""
annotations = env.Command(
    target='$out/annotations.csv',
    source=refresh_annotations,
    action=['taxit -v update_taxids '
            '--name-column organism '
            '--out $TARGET '
            '--schema $schema '
            '$SOURCE $tax_url',
            Copy('$annotations_cache', '$TARGET')])

'''
append new records to global list

NOTE: see bin/dedup_gb.py for global genbank record maintenance
'''
env.Command(
    target=None,
    source=gbs,
    action='cat $SOURCE >> $genbank_cache')

"""
Remove and re-append is_type column with sequences per discussion:

https://gitlab.labmed.uw.edu/uwlabmed/mkrefpkg/issues/14
"""
is_type = env.Command(
    target='$out/is_type.csv',
    source=[annotations, types],
    action=('is_type.py $SOURCES $TARGET'))

"""
pull sequences at least 1200 bp and less than 1% ambiguous
"""
full_fa, full_annotations = env.Command(
    target=['$out/1200bp/seqs.fasta', '$out/1200bp/seq_info.csv'],
    source=[fa, is_type],
    action=('partition_refs.py '
            # filtering
            '--min-length 1200 '
            '--prop-ambig-cutoff 0.01 '
            '--out-fa ${TARGETS[0]} '
            '--out-annotations ${TARGETS[1]} '
            '$SOURCES'))

"""
filter for valid annotations
"""
valid_annotations = env.Command(
    target='$out/1200bp/valid/annotations.csv',
    source=full_annotations,
    action='is_valid.py --out $TARGET  --schema $schema $SOURCE $tax_url')

"""
Make general valid taxtable with all ranks included
"""
valid_tax = env.Command(
    target='$out/1200bp/valid/taxonomy.csv',
    source=valid_annotations,
    action=('taxit -v taxtable '
            '--seq-info $SOURCE '
            '--out $TARGET '
            '--schema $schema '
            '$tax_url'))

"""
pull valid sequences based on valid and ranked tax_ids
"""
valid_fa = env.Command(
    target='$out/1200bp/valid/seqs.fasta',
    source=[full_fa, valid_annotations],
    action=('partition_refs.py '
            '--tax-ids ${SOURCES[1]} '
            '--out-fa $TARGET '
            '${SOURCES[0]}'))

"""
Count reference sequences that can be used for classification filtering
"""
env.Command(
    target='$out/1200bp/valid/tax_counts.csv',
    source=valid_tax,
    action='taxit -v count_taxids --out $TARGET $SOURCE')

"""
Create blast database valid seqs
"""
blast_db(env, valid_fa, '$out/1200bp/valid/blast')

"""
Partition type strains
"""
types_fasta, types_annotations = env.Command(
    target=['$out/1200bp/valid/types/seqs.fasta',
            '$out/1200bp/valid/types/seq_info.csv'],
    source=[valid_fa, valid_annotations],
    action=('partition_refs.py '
            '--types '
            '--out-fa ${TARGETS[0]} '
            '--out-annotations ${TARGETS[1]} '
            '$SOURCES'))

"""
Make type strain taxtable
"""
types_tax = env.Command(
    target='$out/1200bp/valid/types/taxonomy.csv',
    source=[valid_tax, types_annotations],
    action=('taxit -v taxtable '
            '--ranked columns '
            '--taxtable ${SOURCES[0]} '
            '--seq-info ${SOURCES[1]} '
            '--out $TARGET '
            '--schema $schema '
            '$tax_url'))

"""
Count reference sequences that can be used for classification filtering
"""
env.Command(
    target='$out/1200bp/valid/types/tax_counts.csv',
    source=types_tax,
    action='taxit -v count_taxids --out $TARGET $SOURCE')

"""
Make type strain Blast database
"""
blast_db(env, types_fasta, '$out/1200bp/valid/types/blast')

"""
Deduplicate sequences by isolate (accession)

Appends a weight column showing number of sequences being represented
"""
dedup_fa, dedup_info = env.Command(
    target=['$out/1200bp/valid/dedup/seqs.fasta',
            '$out/1200bp/valid/dedup/seq_info.csv'],
    source=[valid_fa, valid_annotations],
    action=('deenurp -v deduplicate_sequences '
            '--group-by accession '
            '--prefer-columns is_type '
            '$SOURCES $TARGETS'))

filtered_details_in = env.Command(
    source='$base/filtered_details.csv',
    target='$out/1200bp/valid/filtered/details_in.csv',
    action=('taxit update_taxids '
            '--schema $schema '
            '--out $TARGET '
            '$SOURCE $tax_url'))

"""
Filter sequences. Use --threads if you need to to limit the number
of processes - otherwise deenurp will use all of them!
"""
filtered_fa, filtered_info, filtered_details, deenurp_log = env.Command(
    source=[dedup_fa, dedup_info, valid_tax, filtered_details_in],
    target=['$out/1200bp/valid/filtered/seqs.fasta',
            '$out/1200bp/valid/filtered/seq_info.csv',
            '$out/1200bp/valid/filtered/details_out.csv',
            '$out/1200bp/valid/filtered/deenurp.log'],
    action=['deenurp -vvv filter_outliers '
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
            # cache this
            Copy('$base/filtered_details.csv', '${TARGETS[2]}')])

"""
Make filtered taxtable
"""
filtered_tax = env.Command(
    target='$out/1200bp/valid/filtered/taxonomy.csv',
    source=[valid_tax, filtered_info],
    action=('taxit -v taxtable '
            '--ranked columns '
            '--taxtable ${SOURCES[0]} '
            '--seq-info ${SOURCES[1]} '
            '--out $TARGET '
            '--schema $schema '
            '$tax_url'))

'''
find top hit for each sequence among type strains

take --maxaccepts 2 so we can filter out hits where the query and target
sequences are the same
'''
type_hits = env.Command(
    target='$out/1200bp/valid/filtered/vsearch.blast6out',
    source=[dedup_fa, types_fasta],
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

FIXME: BokehDeprecationWarning 2> /dev/null

TODO: include derep_map.csv to input for
"dereplicated" sequence counts per accession
'''
env.Command(
    target=['$out/1200bp/valid/filtered/index.html',
            '$out/1200bp/valid/filtered/plots/map.csv'],
    source=[filtered_details,
            valid_tax,
            type_hits,
            types_annotations,
            deenurp_log],
    action=('plot_details.py ${SOURCES[:4]} '
            '--param strategy:cluster '
            '--param cluster_type:single '
            '--param distance_percentile:90.0 '
            '--param min_distance:0.01 '
            '--param max_distance:0.02 '
            '--log-in ${SOURCES[4]} '
            '--plot-dir $out/1200bp/valid/filtered/plots '
            '--plot-map ${TARGETS[1]} '
            '--plot-index ${TARGETS[0]}'
            ))

blast_db(env, filtered_fa, '$out/1200bp/valid/filtered/blast')

"""
Count reference sequences that can be used for classification filtering
"""
env.Command(
    target='$out/1200bp/valid/filtered/tax_counts.csv',
    source=filtered_tax,
    action='taxit -v count_taxids --out $TARGET $SOURCE')

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
