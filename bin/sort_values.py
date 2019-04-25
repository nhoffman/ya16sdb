#!/usr/bin/env python3
"""
sort sequences in REVERSE preferred order for makeblastdb
"""
import argparse
import hashlib
import pandas
import sys

ascending = {
    'ambig_count': False,
    'download_date': True,
    'is_published': True,
    'is_refseq': True,
    'is_type': True,
    'length': True,
    'modified_date': True,
    'seqhash': False,
    'version_num': True
}


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        'feather', help='full seq_info file with description column')
    p.add_argument(
        'sort_values',
        help='Database string URI or filename.')
    p.add_argument(
        '--out',
        default=sys.stdout,
        type=argparse.FileType('w'),
        help='feather file md5sum [stdout')
    args = p.parse_args()
    info = pandas.read_feather(args.feather)
    by = args.sort_values.split(',')
    info = info.sort_values(by=by, ascending=[ascending[i] for i in by])
    info.reset_index(drop=True).to_feather(args.feather)
    args.out.write(hashlib.md5(open(args.feather, 'rb').read()).hexdigest())


if __name__ == '__main__':
    main()
