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
    'assembly_refseq'
]

DTYPES = {
    'accession': str,
    'ambig_count': int,
    'description': str,
    'download_date': str,
    'isolate': str,
    'isolation_source': str,
    'keywords': str,
    'length': int,
    'modified_date': str,
    'mol_type': str,
    'name': str,
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
