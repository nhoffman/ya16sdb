"""
Download and curate the NCBI 16S rRNA sequences

TODO: download records updated before last download
"""
import configparser
import errno
import os
import sys
import time

# requirements installed in the virtualenv
import SCons
from SCons.Script import Variables, ARGUMENTS, Help, Environment, PathVariable

venv = os.environ.get('VIRTUAL_ENV')
if not venv:
    sys.exit('--> an active virtualenv is required'.format(venv))
if not os.path.exists('settings.conf'):
    sys.exit("Can't find settings.conf")
conf = configparser.SafeConfigParser()
conf.read('settings.conf')
paths = conf['paths']


def PathIsFileCreate(key, val, env):
    """namedator to check if Path is a cache file,
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


true_vals = ['t', 'y', '1']
release = ARGUMENTS.get('release', 'no').lower()[0] in true_vals
test = ARGUMENTS.get('test', 'no').lower()[0] in true_vals
vrs = Variables(None, ARGUMENTS)
if test:
    vrs.Add('base', help='Path to output directory', default='test_output')
else:
    vrs.Add('base', help='Path to output directory', default='output')
vrs.Add('out', default=os.path.join('$base', time.strftime('%Y%m%d')))
vrs.Add('api_key', 'ncbi api_key for downloading data', paths['api_key'])
vrs.Add('email', 'email address for ncbi', 'crosenth@uw.edu')
vrs.Add('retry', 'ncbi retry milliseconds', '60000')
vrs.Add('tax_url', default=paths['taxonomy'], help='database url')
# cache vars
vrs.Add(PathVariable(
    'genbank_cache', '', '$base/records.gb', PathIsFileCreate))
vrs.Add(PathVariable(
    'seqs_cache', '', '$base/seqs.fasta', PathIsFileCreate))
vrs.Add(PathVariable(
    'seq_info_cache', '', '$base/seq_info.csv', PathIsFileCreate))
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
        '{singularity} exec '
        '--bind {taxonomy} '
        '--bind $$(readlink -f $$(pwd)) '
        '--pwd $$(readlink -f $$(pwd)) '
        '{taxtastic} taxit'.format(**paths)),
    deenurp=(
        '{singularity} exec '
        '--bind $$(readlink -f $$(pwd)) '
        '--pwd $$(readlink -f $$(pwd)) '
        '{deenurp} deenurp'.format(**paths))
)

env.Decider('MD5-timestamp')

Help(vrs.GenerateHelpText(env))

rrna_16s = ('16s[All Fields] '
            'AND rRNA[Feature Key] '
            'AND Bacteria[Organism] '
            'AND 500 : 99999999999[Sequence Length]')

mefetch_acc = ('mefetch -vv '
               '-api-key $api_key '
               '-email $email '
               '-format acc '
               '-log $out/ncbi.log '
               '-max-retry -1 '  # continuous retries
               '-mode text '
               '-out $TARGET '
               '-retry $retry')

"""
get accessions (versions) of records considered type strains
NCBI web blast uses `sequence_from_type[Filter]` so we will use that
http://www.ncbi.nlm.nih.gov/news/01-21-2014-sequence-by-type/
"""
types = env.Command(
    source=None,
    target='$out/dedup/1200bp/types/esearch.txt',
    action=('esearch -db nucleotide -query "' + rrna_16s +
            ' AND sequence_from_type[Filter]" | ' + mefetch_acc))

"""
Download TM7 accessions
"""
tm7 = env.Command(
    source=None,
    target='$out/esearch/tm7.txt',
    action=('esearch -db nucleotide -query "' + rrna_16s +
            ' AND Candidatus Saccharibacteria[Organism]" | '
            '' + mefetch_acc))

if test:
    tax_ids = (i.strip() for i in open('testfiles/tax_ids.txt') if i)
    tax_ids = ('txid' + i + '[Organism]' for i in tax_ids)
    esearch = env.Command(
        target='$out/esearch/records.txt',
        source='testfiles/tax_ids.txt',
        action=('esearch -db nucleotide -query "' + rrna_16s +
                ' AND (' + ' OR '.join(tax_ids) + ')" | ' + mefetch_acc))
else:
    """
    Download accessions for bacteria with 16s rrna seq_info
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
    concat our download set
    """
    esearch = env.Command(
        source=[classified, tm7],
        target='$out/esearch/records.txt',
        action='cat $SOURCES > $TARGET')

'''
filter out ignored accessions
'''
esearch = env.Command(
    source=['data/ignore.txt', esearch],
    target='$out/esearch.txt',
    action='grep --invert-match --file $SOURCES > $TARGET')

"""
Do not download record accessions in the ignore list or
that have been previously downloaded in the records_cache.
Exit script if no new records exist.
"""
new = env.Command(
    target='$out/new/records.txt',
    source=[esearch, '$records_cache'],
    action=['cat ${SOURCES[1]} | '
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
            '-api-key $api_key '
            '-db nucleotide '
            '-email $email '
            '-format ft '
            '-id $SOURCE '
            '-log $out/ncbi.log '
            '-max-retry -1 '  # continuous retry
            '-mode text '
            '-retmax 1 '
            '-retry $retry '
            '| '
            # extract 16s features
            'ftract -feature "rrna:product:16S ribosomal RNA" '
            '-log $out/ncbi.log -on-error continue | '
            'accession_version.py | '  # parse accession.version from id column
            'mefetch '  # download genbank records
            '-vv '
            '-api-key $api_key '
            '-csv '
            '-db nucleotide '
            '-email $email '
            '-format gbwithparts '
            '-log $out/ncbi.log '
            '-max-retry -1 '  # continuous retry
            '-mode text '
            '-retmax 1 '
            '-retry $retry > $TARGET'])

"""
extract
"""
today = time.strftime('%d-%b-%Y')
new_fa, new_info, new_pub_info, new_refs, new_refseq_info = env.Command(
    target=['$out/new/seqs.fasta',
            '$out/new/seq_info.csv',
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
    source=[new_info, new],
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
filter innamed tax_ids

Do nothing with the unknown records for now because it might simply mean
the ncbi taxonomy pipeline is out of sync with the latest 16s records
"""
known_info, _ = env.Command(
    target=['$out/new/taxit/seq_info.csv',
            '$out/new/taxit/unknown.csv'],
    source=new_info,
    action=['$taxit update_taxids '
            '--unknown-action drop '
            '--unknowns ${TARGETS[1]} '
            '--outfile ${TARGETS[0]} '
            '$SOURCE $tax_url'])

"""
vsearch new sequences with training set to test sequence orientation
and 16s region

TODO: move to cmsearch
"""
vsearch = env.Command(
    target='$out/new/vsearch.csv',
    source=[new_fa, 'data/rdp_16s_type_strains.fasta.bz2'],
    action=('vsearch '
            '--usearch_global ${SOURCES[0]} '
            '--db ${SOURCES[1]} '
            '--iddef 2 '  # matching cols / alignment len excluding term gaps
            '--id 0.70 '
            '--query_cov 0.70 '
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
fa, refresh_info, pubmed_info, references, refseq_info, _ = env.Command(
    target=['$out/seqs.fasta',
            '$out/refresh/seq_info.csv',
            '$out/pubmed_info.csv',
            '$out/references.csv',
            '$out/refseq_info.csv',
            '$out/records.txt'],
    source=[esearch,
            vsearch_fa, '$seqs_cache',
            known_info, '$seq_info_cache',
            new_pub_info, '$pubmed_info_cache',
            new_refs, '$references_cache',
            new_refseq_info, '$refseq_info_cache',
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
update_seq_info = env.Command(
    target='$out/refresh/taxit/seq_info.csv',
    source=refresh_info,
    action='$taxit -v update_taxids --out $TARGET $SOURCE $tax_url')

"""
is_type column

https://github.com/nhoffman/ya16sdb/issues/11 about is_type column
https://gitlab.labmed.uw.edu/uwlabmed/mkrefpkg/issues/40 about confidence col
"""
is_type_seq_info = env.Command(
    target='$out/refresh/taxit/is_type/seq_info.csv',
    source=[update_seq_info, types],
    action='is_type.py --out $TARGET $SOURCES')

"""
confidence column

https://gitlab.labmed.uw.edu/uwlabmed/mkrefpkg/issues/40 about confidence col
"""
confidence_seq_info = env.Command(
    target='$out/refresh/taxit/is_type/confidence/seq_info.csv',
    source=[is_type_seq_info, pubmed_info],
    action='confidence.py --out $TARGET $SOURCES')

"""
taxonomy
"""
taxonomy = env.Command(
    target='$out/taxonomy.csv',
    source=confidence_seq_info,
    action='$taxit -v taxtable --seq-info $SOURCE --out $TARGET $tax_url')

"""
seq_info with species column
"""
seq_info = env.Command(
    target='$out/seq_info.csv',
    source=[confidence_seq_info, taxonomy],
    action='merge.py --right-columns tax_id,species --out $TARGET $SOURCES')

"""
Deduplicate sequences by isolate (accession)

Weight column is meaningless and removed to avoid confusion
"""
dedup_fa, dedup_info = env.Command(
    target=['$out/dedup/seqs.fasta', '$out/dedup/seq_info.csv'],
    source=[fa, seq_info],
    action=('$deenurp deduplicate_sequences '
            '--group-by accession '
            '$SOURCES ${TARGETS[0]} /dev/stdout | '
            # remove weight column because it is meaningless
            'csvcut.py --not-columns weight --out ${TARGETS[1]}'))

"""
pull sequences at least 1200 bp, less than 1% ambiguous and with species tax id
"""
full_fa, full_seq_info = env.Command(
    target=['$out/dedup/1200bp/seqs.fasta',
            '$out/dedup/1200bp/seq_info.csv'],
    source=[dedup_fa, dedup_info],
    action=('partition_refs.py '
            '--min-length 1200 '
            '--prop-ambig-cutoff 0.01 '
            '--species '
            '$SOURCES $TARGETS'))

"""
generate full taxonomy table. Needed (at least) when running to_feather.py
"""
full_tax = env.Command(
    target='$out/dedup/1200bp/taxonomy.csv',
    source=full_seq_info,
    action='$taxit -v taxtable --seq-info $SOURCE --out $TARGET $tax_url')

"""
pull type sequences
"""
type_fa, type_info = env.Command(
    target=['$out/dedup/1200bp/types/seqs.fasta',
            '$out/dedup/1200bp/types/seq_info.csv'],
    source=[fa, full_seq_info],
    action='partition_refs.py --is_type $SOURCES $TARGETS')

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
filter for named and tm7 sequence records
"""
named_fa, named_info = env.Command(
    target=['$out/dedup/1200bp/named/seqs.fasta',
            '$out/dedup/1200bp/named/seq_info.csv'],
    source=[fa, full_seq_info, tm7],
    action=('named.py $tax_url | '
            'partition_refs.py '
            '--versions ${SOURCES[2]} '
            '--tax-ids /dev/stdin '
            '${SOURCES[:2]} $TARGETS'))

"""
Make general named taxtable with all ranks included for filter_outliers
"""
named_tax = env.Command(
    target='$out/dedup/1200bp/named/taxonomy.csv',
    source=named_info,
    action='$taxit -v taxtable --seq-info $SOURCE --out $TARGET $tax_url')

"""
Trimmed seqname,tax_id map file that can be easily cached by filter_outliers
"""
named_taxid_map = env.Command(
    target='$out/dedup/1200bp/named/tax_id_map.csv',
    source=named_info,
    action='csvcut.py --columns seqname,tax_id --out $TARGET $SOURCE')

"""
Create taxtable output with replacing tax_ids with taxnames
"""
named_lineages = env.Command(
    target='$out/dedup/1200bp/named/lineages.csv',
    source=[named_tax, named_info],
    action='$taxit lineage_table --csv-table $TARGET $SOURCES')

"""
Create taxtable output with replacing tax_ids with taxnames

https://mothur.org/wiki/Taxonomy_File
"""
named_mothur = env.Command(
    target='$out/dedup/1200bp/named/lineages.txt',
    source=[named_tax, named_info],
    action='$taxit lineage_table --taxonomy-table $TARGET $SOURCES')

"""
get all taxid descendants from trusted_taxids.txt file
"""
trusted_taxids = env.Command(
    target='$out/dedup/1200bp/named/filtered/trusted_taxids.txt',
    source=paths['trusted_taxids'],
    action=[
        # copy trusted_taxids.txt to current dir bind point
        'cp $SOURCE $$(dirname $TARGET)/$$(basename $SOURCE)',
        '$taxit get_descendants --out $TARGET $tax_url $SOURCE'])

"""
update tax_ids in details_in cache
"""
filtered_details_in = env.Command(
    source='$outliers_cache',
    target='$out/dedup/1200bp/named/filtered/details_in.csv',
    action=('$taxit update_taxids '
            '--outfile $TARGET '
            '$SOURCE $tax_url'
            # continue if filter_outliers cache is empty
            ' || echo "Continuing without filter_outliers cache"'))

"""
Filter sequences. Use --threads if you need to to limit the number
of processes - otherwise deenurp will use all of them!
"""
filtered_fa, filtered_taxid_map, filtered_details, deenurp_log = env.Command(
    source=[named_fa, named_taxid_map, named_tax,
            filtered_details_in, trusted_taxids],
    target=['$out/dedup/1200bp/named/filtered/seqs.fasta',
            '$out/dedup/1200bp/named/filtered/tax_id_map.csv',
            '$out/dedup/1200bp/named/filtered/details_out.csv',
            '$out/dedup/1200bp/named/filtered/deenurp.log'],
    action=['$deenurp -vvv filter_outliers '
            '--cluster-type single '
            '--detailed-seqinfo ${TARGETS[2]} '
            '--distance-percentile 90.0 '
            '--filter-rank species '
            '--filtered-seqinfo ${TARGETS[1]} '
            '--jobs 1 '
            '--log ${TARGETS[3]} '
            '--max-distance 0.02 '
            '--min-distance 0.01 '
            '--min-seqs-for-filtering 5 '
            '--no-filter ${SOURCES[4]} '
            '--output-seqs ${TARGETS[0]}  '
            '--previous-details ${SOURCES[3]} '
            '--strategy cluster '
            '--threads-per-job 14 '
            '${SOURCES[:3]}',
            # cache it
            'cp ${TARGETS[2]} $outliers_cache'])

"""
Create a filtered seq_info.csv file
"""
filtered_info = env.Command(
    target='$out/dedup/1200bp/named/filtered/seq_info.csv',
    source=[filtered_taxid_map, named_info],
    action='merge.py --out $TARGET $SOURCES')

"""
Make general named taxtable with all ranks included for filter_outliers
"""
filtered_tax = env.Command(
    target='$out/dedup/1200bp/named/filtered/taxonomy.csv',
    source=filtered_info,
    action='$taxit -v taxtable --seq-info $SOURCE --out $TARGET $tax_url')

"""
Create taxtable output with replacing tax_ids with taxnames
"""
filtered_lineages = env.Command(
    target='$out/dedup/1200bp/named/filtered/lineages.csv',
    source=[filtered_tax, filtered_info],
    action='$taxit lineage_table --csv-table $TARGET $SOURCES')

"""
Create mothur output

https://mothur.org/wiki/Taxonomy_File
"""
filtered_mothur = env.Command(
    target='$out/dedup/1200bp/named/filtered/lineages.txt',
    source=[filtered_tax, filtered_info],
    action='$taxit lineage_table --taxonomy-table $TARGET $SOURCES')

"""
labmed do-not-trust filtering

FIXME:  This will process will eventually go in its own repo but consider
doing a reverse grep into partition_refs.py
"""
trusted_fa, trusted_info = env.Command(
    target=['$out/dedup/1200bp/named/filtered/trusted/seqs.fasta',
            '$out/dedup/1200bp/named/filtered/trusted/seq_info.csv'],
    source=[paths['do_not_trust'], filtered_fa, named_info],
    action='do_not_trust.py $SOURCES $TARGETS')

"""
Make filtered taxtable with ranked columns and no_rank rows
"""
trusted_tax = env.Command(
    target='$out/dedup/1200bp/named/filtered/trusted/taxonomy.csv',
    source=filtered_info,
    action='$taxit taxtable --seq-info $SOURCE --out $TARGET $tax_url')

"""
Taxtable output replacing tax_ids with taxnames
"""
trusted_lineages = env.Command(
    target='$out/dedup/1200bp/named/filtered/trusted/lineages.csv',
    source=[trusted_tax, trusted_info],
    action='$taxit lineage_table --csv-table $TARGET $SOURCES')

"""
Mothur output

https://mothur.org/wiki/Taxonomy_File
"""
trusted_mothur = env.Command(
    target='$out/dedup/1200bp/named/filtered/trusted/lineages.txt',
    source=[trusted_tax, trusted_info],
    action='$taxit lineage_table --taxonomy-table $TARGET $SOURCES')

blast_db(env, trusted_fa, '$out/dedup/1200bp/named/filtered/trusted/blast')

'''
Pull type strains from trusted db
'''
trusted_type_fa, trusted_type_info = env.Command(
    target=['$out/dedup/1200bp/named/filtered/trusted/types/seqs.fasta',
            '$out/dedup/1200bp/named/filtered/trusted/types/seq_info.csv'],
    source=[trusted_fa, trusted_info],
    action='partition_refs.py --is_type $SOURCES $TARGETS')

"""
Make trusted type taxtable with ranked columns and no_rank rows
"""
trusted_type_tax = env.Command(
    target='$out/dedup/1200bp/named/filtered/trusted/types/taxonomy.csv',
    source=trusted_type_info,
    action=('$taxit taxtable '
            '--seq-info $SOURCE '
            '--out $TARGET '
            '$tax_url'))

"""
Create taxtable output with replacing tax_ids with taxnames
"""
trusted_type_lineages = env.Command(
    target='$out/dedup/1200bp/named/filtered/trusted/types/lineages.csv',
    source=[trusted_type_tax, trusted_type_info],
    action='$taxit lineage_table --csv-table $TARGET $SOURCES')

"""
Mothur output

https://mothur.org/wiki/Taxonomy_File
"""
trusted_type_mothur = env.Command(
    target='$out/dedup/1200bp/named/filtered/trusted/types/lineages.txt',
    source=[trusted_type_tax, trusted_type_info],
    action='$taxit lineage_table --taxonomy-table $TARGET $SOURCES')

blast_db(
    env,
    trusted_type_fa,
    '$out/dedup/1200bp/named/filtered/trusted/types/blast')

'''
find top hit for each sequence among type strains

NOTE: alleles will never align with themselves (--self) BUT
can align with other alleles in the same genome accession
'''
named_type_hits = env.Command(
    target='$out/dedup/1200bp/named/vsearch.tsv',
    source=[named_fa, trusted_type_fa],
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
    target=['$out/dedup/1200bp/named/filtered/index.html',
            '$out/dedup/1200bp/named/filtered/plots/map.csv'],
    source=[filtered_details,
            named_type_hits,
            seq_info,
            named_tax,
            types,
            deenurp_log],
    action=('plot_details.py ${SOURCES[:5]} '
            '--param strategy:cluster '
            '--param cluster_type:single '
            '--param distance_percentile:90.0 '
            '--param min_distance:0.01 '
            '--param max_distance:0.02 '
            '--log-in ${SOURCES[5]} '
            '--plot-dir $out/dedup/1200bp/named/filtered/plots '
            '--plot-map ${TARGETS[1]} '
            '--plot-index ${TARGETS[0]}'))

"""
feather output - https://github.com/wesm/feather
"""
filtered_feather = env.Command(
    target='$out/dedup/1200bp/named/filtered/filter_details.feather.gz',
    source=[filtered_details, full_seq_info, full_tax,
            named_type_hits, pubmed_info],
    action='to_feather.py --out $TARGET $SOURCES')

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
