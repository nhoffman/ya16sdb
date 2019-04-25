#!/usr/bin/env python
"""
Create a seq_info.csv feather file optionally gzipped
"""
import argparse
import pandas
import sys

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


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('seq_info')
    parser.add_argument('out')
    args = parser.parse_args(arguments)
    seq_info = pandas.read_csv(
        args.seq_info,
        dtype=DTYPES,
        parse_dates=['download_date', 'modified_date'],
        date_parser=lambda x: pandas.datetime.strptime(x, '%d-%b-%Y'))
    seq_info.to_feather(args.out)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
