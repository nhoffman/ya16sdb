#!/usr/bin/env python3
"""
Return taxtable with ranked only taxons
"""

import argparse
import csv
import sys


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        'taxtable',
        nargs='?',
        default=sys.stdin,
        type=argparse.FileType('r'),
        help='taxit taxtable')
    p.add_argument(
        '--columns',
        action='store_false',
        dest='rows',
        help='ranked columns with no_rank rows')
    p.add_argument(
        '--endswith',
        default='_',
        help='not ranked')
    p.add_argument(
        '--out',
        default=sys.stdout,
        type=argparse.FileType('w'),
        help='list of all version downloaded')

    args = p.parse_args()
    tax = csv.DictReader(args.taxtable)
    columns = tax.fieldnames
    columns = [c for c in columns if not c.endswith(args.endswith)]
    if args.rows:
        tax = (r for r in tax if not r['rank'].endswith(args.endswith))
    writer = csv.DictWriter(args.out, columns, extrasaction='ignore')
    writer.writeheader()
    writer.writerows(tax)


if __name__ == '__main__':
    main()
