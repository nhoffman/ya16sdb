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
import warnings

from SCons.Script import ARGUMENTS, GetBuildFailures, Depends

venv = os.environ.get('VIRTUAL_ENV')
if not venv:
    warnings.warn('No active virtualenv detected, using system environment')
if not os.path.exists('settings.conf'):
    sys.exit("Can't find settings.conf")


class Environment(SCons.Environment.Environment):
    def __init__(self, singularity=None, verbosity=0, **kws):
        self.singularity = singularity
        self.verbosity = verbosity
        SCons.Environment.Environment.__init__(self, **kws)

    def Command(self,
                target,
                source,
                action,
                singularity=None,
                options=None,
                **kws):
        # if docker and self.docker:  # TODO: implement this
        if singularity and self.singularity:
            self.Depends(target, singularity)
            sactions = []
            for a in self.Flatten(action):
                sa = self.singularity
                sa = '{} --bind $pipeline'.format(sa)
                if options:
                    sa = '{} {}'.format(sa, options)
                sa = '{} {} {}'.format(sa, singularity, a)
                sactions.append(VirtualAction(a, sa, self.verbosity))
            action = sactions
        return SCons.Environment.Environment.Command(
            self, target, source, action, **kws)


class VirtualAction(SCons.Action.CommandAction):
    def __init__(self, command, singularity, verbosity, **kw):
        self.command = singularity if verbosity else command
        SCons.Action.CommandAction.__init__(self, singularity, **kw)

    def print_cmd_line(self, _, target, source, env):
        c = env.subst(self.command, SCons.Subst.SUBST_RAW, target, source)
        SCons.Action.CommandAction.print_cmd_line(self, c, target, source, env)


def blast_db(env, sequence_file, output_base):
    '''
    Create a blast database and file md5sum
    '''
    extensions = ['.nhr', '.nin', '.nsq']
    blast_out = env.Command(
        target=[output_base + ext for ext in extensions],
        source=sequence_file,
        action=('makeblastdb -dbtype nucl '
                '-in $SOURCE -out ' + output_base),
        singularity=settings['blast'])
    env.Command(
        target=output_base,
        source=blast_out,
        action='md5sum $SOURCES > $TARGET')
    return blast_out


def taxonomy(fa, info, path):
    """
    Make filtered taxtable with ranked columns and no_rank rows
    """
    taxtable = env.Command(
        target=os.path.join(path, 'taxonomy.csv'),
        source=info,
        action=['taxit taxtable '
                '--seq-info $SOURCE '
                '--out $TARGET '
                '$tax_url'],
        singularity=settings['taxit'],
        options='--bind $tax_url')

    """
    Taxtable output replacing tax_ids with taxnames
    """
    lineages = env.Command(
        target=os.path.join(path, 'lineages.csv'),
        source=[taxtable, info],
        action='taxit lineage_table --csv-table $TARGET $SOURCES',
        singularity=settings['taxit'])

    """
    Mothur output - https://mothur.org/wiki/Taxonomy_File
    """
    mothur = env.Command(
        target=os.path.join(path, 'lineages.txt'),
        source=[taxtable, info],
        action='taxit lineage_table --taxonomy-table $TARGET $SOURCES',
        singularity=settings['taxit'])

    blast = blast_db(env, fa, os.path.join(path, 'blast'))

    return taxtable, lineages, mothur, blast


true_vals = ['t', 'y', '1']
release = ARGUMENTS.get('release', 'no').lower()[0] in true_vals
test = ARGUMENTS.get('test', 'no').lower()[0] in true_vals
absolute_dir = os.path.dirname((lambda x: x).__code__.co_filename)
conf = configparser.SafeConfigParser()
conf.read(['settings.conf', os.path.join(absolute_dir, 'ncbi.conf')])
settings = conf['TEST'] if test else conf['DEFAULT']
out = os.path.join(settings['outdir'], time.strftime('%Y%m%d'))
cachedir = settings['cachedir']

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

environment_variables = dict(
    os.environ,
    PATH=':'.join([
        os.path.join(absolute_dir, 'bin'),
        (os.path.join(venv, 'bin') if venv else ''),
        '/usr/local/bin',
        '/usr/bin',
        '/bin']),
)

env = Environment(
    ENV=environment_variables,
    log=os.path.join(out, 'log.txt'),
    out=out,
    pipeline=absolute_dir,
    shell='bash',
    tax_url=os.path.abspath(settings['taxonomy']),
    verbosity=1,
    **settings
)

env.EnsureSConsVersion(3, 0, 5)
env.Decider('MD5')

mefetch = ('mefetch -vv '
           '-api-key $api_key '  # FIXME: what happens if no $api_key?
           '-db nucleotide '
           '-email $email '
           '-log $log '
           '-max-retry -1 '
           '-mode text '
           '-proc $nreq '
           '-retry $retry')

classified = env.Command(
    source=None,
    target='$out/ncbi/classified.txt',
    action=['esearch -db nucleotide -query "$classified" | ' +
            mefetch + ' -format acc -out $TARGET'],
    singularity=settings['eutils'])

"""
Candidatus Saccharibacteria
https://gitlab.labmed.uw.edu/molmicro/mkrefpkg/issues/36
"""
tm7 = env.Command(
    source=None,
    target='$out/ncbi/tm7.txt',
    action=['esearch -db nucleotide -query "$tm7" | ' +
            mefetch + ' -format acc -out $TARGET'],
    singularity=settings['eutils'])

"""
Check the cache for last download_date and download list of modified
records in order to re-download modified records that we have previously
downloaded.
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
    action=[
        'esearch -db nucleotide -query "($classified OR $tm7) '
        'AND $download_date[Modification Date] : 3000[Modification Date]" | ' +
        mefetch + ' -format acc -out $TARGET'],
    singularity=settings['eutils'])

"""
type strains records
NCBI web blast uses `sequence_from_type[Filter]` so we will use that
http://www.ncbi.nlm.nih.gov/news/01-21-2014-sequence-by-type/
"""
types = env.Command(
    source=None,
    target='$out/ncbi/types.txt',
    action=['esearch -db nucleotide -query "$types" | ' +
            mefetch + ' -format acc -out $TARGET'],
    singularity=settings['eutils'])

"""
Trim accession2taxid with 16s records and update taxids
"""
accession2taxid = env.Command(
    target='$out/ncbi/accession2taxid.csv',
    source=[settings['accession2taxid'], classified, tm7],
    action='accession2taxid.py --out $TARGET $SOURCES')

accession2taxid = env.Command(
    target='$out/ncbi/accession2update_taxids.csv',
    source=accession2taxid,
    action=['taxit update_taxids '
            '--outfile $TARGET '
            '--unknown-action drop '
            '$SOURCE $tax_url'],
    singularity=settings['taxit'],
    options='--bind $tax_url')

"""
Create a list of cached records removing:

1. Modified records
2. Records not included in esearch query
3. Records with new tax_ids

Modified and records with new tax_ids will be redownloaded

FIXME: call this cache_in.py (cache_out.py will replace refresh.py)
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
download genbank records
"""
gbs = env.Command(
    target='$out/ncbi/download.gb',
    source=download,
    action=[mefetch + ' -format ft -id $SOURCE -retmax 1 | '
            'ftract -feature "rrna:product:16S ribosomal RNA" '
            '-log $log '
            '-on-error continue '
            '-min-length 1200 | ' +
            mefetch + ' -csv -format gbwithparts -out $TARGET -retmax 1'])

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
    action=['taxit update_taxids '
            '--unknown-action drop '
            '${SOURCES[1]} $tax_url | '
            'partition_refs.py ${SOURCES[0]} - $TARGETS'],
    singularity=settings['taxit'],
    options='--bind $tax_url')

"""
cmsearch new sequences against rfam model
"""
cmsearch = env.Command(
    target='$out/ncbi/extract_genbank/update_taxids/cmsearch/table.tsv',
    source=['$pipeline/data/SSU_rRNA_bacteria.cm', fa],
    action=('cmsearch --cpu 14 -E 0.01 --hmmonly -o /dev/null '
            '--tblout $TARGET $SOURCES || true'),
    singularity=settings['infernal'])

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

taxtable = env.Command(
    target='$out/taxonomy.csv',
    source=seq_info,
    action='taxit -v taxtable --seq-info $SOURCE --out $TARGET $tax_url',
    singularity=settings['taxit'],
    options='--bind $tax_url')

"""
Map for WGS records without a refseq assembly accession
"""
asm = env.Command(
    target='$out/ani/assembly_summary_genbank.txt',
    source=None,
    action=('wget '
            '--output-document $TARGET '
            '--quiet '
            '--retry-on-http-error 403 '
            '--tries 100 '
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
            '--quiet '
            '--retry-on-http-error 403 '
            '--tries 100 '
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

fa, seq_info = env.Command(
    target=['$out/seqs.fasta', '$out/seq_info.csv'],
    source=[fa, feather],
    action='partition_refs.py --drop-noaligns $SOURCES $TARGETS')

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

type_tax, type_lineages, types_mothur, type_blast = taxonomy(
    type_fa, type_info, '$out/dedup/1200bp/types/')

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

taxtable, lineages, mothur, blast = taxonomy(
    fa, seq_info, '$out/dedup/1200bp/named/')

"""
Trimmed seqname,tax_id map file that can be easily cached by filter_outliers
"""
taxid_map = env.Command(
    target='$out/dedup/1200bp/named/tax_id_map.csv',
    source=seq_info,
    action='csvcut.py --columns seqname,tax_id --out $TARGET $SOURCE')

"""
Filter sequences. Use --threads if you need to to limit the number
of processes - otherwise deenurp will use all of them!
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
            'cp ${TARGETS[1]} ' + cfiles['outliers_cache']],
    singularity=settings['deenurp'])

"""
add distance metrics to feather file
"""
filter_outliers = env.Command(
    target='$out/.feather/filter_outliers.md5',
    source=[feather, details],
    action=['filter_outliers.py $SOURCES', 'md5sum ${SOURCES[0]} > $TARGET'])

taxtable, filtered_lineages, mothur, blast = taxonomy(
    fa, seq_info, '$out/dedup/1200bp/named/filtered/')

filtered_type_fa, filtered_type_info = env.Command(
    target=['$out/dedup/1200bp/named/filtered/types/seqs.fasta',
            '$out/dedup/1200bp/named/filtered/types/seq_info.csv'],
    source=[fa, feather],
    action='partition_refs.py --is_type $SOURCES $TARGETS')

filtered_type_tax, filtered_type_lineages, filtered_type_mothur, _ = taxonomy(
    filtered_type_fa,
    filtered_type_info,
    '$out/dedup/1200bp/named/filtered/types/')

filtered_type_hits = env.Command(
    target='$out/dedup/1200bp/named/types_vsearch.tsv',
    source=[named, filtered_type_fa],
    action=('vsearch --usearch_global ${SOURCES[0]} '
            '--blast6out $TARGET '
            '--db ${SOURCES[1]} '
            '--id 0.75 '
            '--maxaccepts 5 '
            '--self '  # reject same sequence hits
            '--strand plus '
            '--threads 14 '
            '--top_hits_only'),
    singularity=settings['vsearch'])

"""
This output will be used in the filter plots
"""
filtered_type_classifications = env.Command(
    target='$out/dedup/1200bp/named/classifications.csv',
    source=[filtered_type_hits,
            filtered_type_info,
            filtered_type_tax,
            os.path.join('$pipeline', 'data/classifier_thresholds.csv')],
    action=['classify -vv '
            '--lineages ${SOURCES[2]} '
            '--rank-thresholds ${SOURCES[3]} '
            '--seq-info ${SOURCES[1]} '
            '--starred 101 '
            '--out $TARGET '
            '${SOURCES[0]}'])

"""
Adds type_classification to feather file for Dash application
"""
type_classifications = env.Command(
    target='$out/.feather/type_classifications.md5',
    source=[feather, filtered_type_classifications],
    action=['type_classifications.py $SOURCES',
            'md5sum ${SOURCES[0]} > $TARGET'])

"""
expand taxids into descendants
"""
trusted_taxids = env.Command(
    target='$out/dedup/1200bp/named/filtered/trusted/taxids.txt',
    source=settings['trust'],
    action='taxit get_descendants --out $TARGET $tax_url $SOURCE',
    singularity=settings['taxit'],
    options='--bind $tax_url')

"""
expand taxids into descendants
"""
dnt_ids = env.Command(
    target='$out/dedup/1200bp/named/filtered/trusted/dnt_ids.txt',
    source=settings['do_not_trust'],
    action='taxit get_descendants --out $TARGET $tax_url $SOURCE',
    singularity=settings['taxit'],
    options='--bind $tax_url')

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
    source=[named, feather, trusted, dnt],
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
'''
named_type_hits = env.Command(
    target='$out/dedup/1200bp/named/trusted_vsearch.tsv',
    source=[named, fa],
    action=('vsearch --usearch_global ${SOURCES[0]} '
            '--blast6out $TARGET '
            '--db ${SOURCES[1]} '
            '--id 0.75 '
            '--maxaccepts 5 '
            '--self '  # reject same sequence hits
            '--strand plus '
            '--threads 14 '
            '--top_hits_only'),
    singularity=settings['deenurp'])

"""
Creates match_seqname, match_pct, match_version, match_species and
match_species_id columns for best type strain hits from the trusted
dataset
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
    source=os.path.join('$pipeline', '.git/objects'),
    action=['(echo $$(hostname):$pipeline;'
            'git --git-dir $pipeline/.git describe --tags --dirty) > $TARGET'])


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
