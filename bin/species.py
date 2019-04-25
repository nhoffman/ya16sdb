#!/usr/bin/env python3
"""
Add is_published column
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
    p.add_argument('taxonomy')
    p.add_argument(
        '--out',
        default=sys.stdout,
        type=argparse.FileType('w'),
        help='feather file md5sum [stdout]')
    args = p.parse_args()
    taxonomy = pandas.read_csv(
        args.taxonomy, dtype=str, usecols=['tax_id', 'species'])
    info = pandas.read_feather(args.feather)
    info = info.merge(taxonomy, how='left')
    info.to_feather(args.feather)
    args.out.write(hashlib.md5(open(args.feather, 'rb').read()).hexdigest())


if __name__ == '__main__':
    main()
