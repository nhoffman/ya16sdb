import re

# assumes 4 or 6 A-Z prefix WGS record prefix
WGS_PREFIX = re.compile(r"^([A-Z]{6}|[A-Z]{4})")

SEQ_INFO_COLS = [
    'seqname',
    'version',
    'accession',
    'name',
    'description',
    'tax_id',
    'modified_date',
    'download_date',
    'version_num',
    'source',
    'keywords',
    'organism',
    'length',
    'ambig_count',
    'strain',
    'mol_type',
    'isolate',
    'isolation_source',
    'seq_start',
    'seq_stop',
    '16s_start',
    '16s_stop',
    'assembly_refseq',
    'master',
    'locus_tag',
    'old_locus_tag'
]

DTYPES = {
    '16s_start': int,
    '16s_stop': int,
    'accession': str,
    'ambig_count': int,
    'assembly_refseq': str,
    'description': str,
    'download_date': str,
    'isolate': str,
    'isolation_source': str,
    'keywords': str,
    'length': int,
    'locus_tag': str,
    'master': str,
    'modified_date': str,
    'mol_type': str,
    'name': str,
    'old_locus_tag': str,
    'organism': str,
    'seq_start': int,
    'seq_stop': int,
    'seqname': str,
    'source': str,
    'strain': str,
    'tax_id': str,
    'version': str,
    'version_num': str,
}


def open_clean(fl):
    if isinstance(fl, str):
        fl = [fl]
    for f in fl:
        f = (row.strip() for row in open(f))
        f = (row for row in f if row)
        for row in f:
            yield row
