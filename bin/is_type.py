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
    p.add_argument(
        'feather', help='full seq_info file with description column')
    p.add_argument(
        'types',
        type=argparse.FileType('r'),
        help='txt file of version numbers that are type strains')
    p.add_argument(
        '--out',
        default=sys.stdout,
        type=argparse.FileType('w'),
        help='feather file md5sum [stdout]')
    args = p.parse_args()
    types = (t.strip() for t in args.types)
    types = set(t for t in types if t)
    info = pandas.read_feather(args.feather)
    info['is_type'] = False
    info.loc[info['version'].isin(types), 'is_type'] = True
    info.to_feather(args.feather)
    args.out.write(hashlib.md5(open(args.feather, 'rb').read()).hexdigest())


if __name__ == '__main__':
    main()
