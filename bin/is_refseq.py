#!/usr/bin/env python3
"""
Add is_type column based on presence in type strain file
"""
import argparse
import hashlib
import pandas
import sys


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('feather')
    p.add_argument('refseqs')
    p.add_argument(
        '--out',
        default=sys.stdout,
        type=argparse.FileType('w'),
        help='feather file md5sum [stdout]')
    args = p.parse_args()
    refseqs = pandas.read_csv(
        args.refseqs, squeeze=True, dtype=str, usecols=['seqname'])
    info = pandas.read_feather(args.feather)
    info['is_refseq'] = False
    info.loc[info['seqname'].isin(refseqs), 'is_refseq'] = True
    info.to_feather(args.feather)
    args.out.write(hashlib.md5(open(args.feather, 'rb').read()).hexdigest())


if __name__ == '__main__':
    main()
