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
    p.add_argument('pubmed_ids')
    p.add_argument(
        '--out',
        default=sys.stdout,
        type=argparse.FileType('w'),
        help='feather file md5sum [stdout]')
    args = p.parse_args()
    published = pandas.read_csv(
        args.pubmed_ids, dtype=str, squeeze=True, usecols=['version'])
    info = pandas.read_feather(args.feather)
    info['is_published'] = False
    info.loc[info['version'].isin(published), 'is_published'] = True
    info.to_feather(args.feather)
    args.out.write(hashlib.md5(open(args.feather, 'rb').read()).hexdigest())


if __name__ == '__main__':
    main()
