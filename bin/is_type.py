#!/usr/bin/env python3
"""
Filter csv file with version column
"""

import argparse
import csv
import sys


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        'seq_info',
        type=argparse.FileType('r'),
        help='full seq_info file with description column')
    p.add_argument(
        'types',
        type=argparse.FileType('r'),
        help='txt file of version numbers that are type strains')
    p.add_argument(
        '--in-description',
        help=('comma delimited list of annotation(s) that may be '
              'in the description field that would that would '
              'designate a type strain'))
    p.add_argument(
        'out',
        type=argparse.FileType('w'),
        default=sys.stdout,
        help='list of all version downloaded')

    args = p.parse_args()
    types = set(t.strip() for t in args.types if t)
    info = csv.reader(args.seq_info, quotechar='"')
    header = next(info)
    version = header.index('version')
    writer = csv.writer(args.out)
    writer.writerows(i for i in info if i[version] in types)


if __name__ == '__main__':
    main()
