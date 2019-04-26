#!/usr/bin/env python3
"""
Add is_type column based on presence in type strain file
"""
import argparse
import pandas


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
    args = p.parse_args()
    types = (t.strip() for t in args.types)
    types = set(t for t in types if t)
    info = pandas.read_feather(args.feather)
    info['is_type'] = False
    info.loc[info['version'].isin(types), 'is_type'] = True
    info.to_feather(args.feather)


if __name__ == '__main__':
    main()
