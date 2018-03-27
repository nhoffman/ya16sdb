#!/usr/bin/env python3
"""
Add is_type column based on presence in type strain file
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
        nargs='?',
        default=sys.stdin,
        help='full seq_info file with description column')
    p.add_argument(
        'types',
        type=argparse.FileType('r'),
        help='txt file of version numbers that are type strains')
    p.add_argument(
        '--out',
        type=argparse.FileType('w'),
        default=sys.stdout,
        help='list of all version downloaded')
    args = p.parse_args()
    types = (t.strip() for t in args.types)
    types = set(t for t in types if t)
    info = csv.DictReader(args.seq_info)
    writer = csv.DictWriter(args.out, fieldnames=info.fieldnames + ['is_type'])
    writer.writeheader()
    for row in info:
        row['is_type'] = row['version'] in types
        writer.writerow(row)


if __name__ == '__main__':
    main()
