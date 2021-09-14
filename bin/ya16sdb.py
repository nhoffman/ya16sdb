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
