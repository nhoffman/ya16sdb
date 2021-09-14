#!/usr/bin/env python3
import argparse
import csv
import sys


def open_clean(fl):
    if isinstance(fl, str):
        fl = [fl]
    for f in fl:
        f = (row.strip() for row in open(f))
        f = (row for row in f if row)
        for row in f:
            yield row


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('a2t')
    p.add_argument('seq_info')
    p.add_argument('modified')
    p.add_argument('cache', nargs=argparse.REMAINDER)
    p.add_argument('--out', type=argparse.FileType('w'), default=sys.stdout)
    args = p.parse_args()
    taxids = []
    for row in csv.DictReader(open(args.a2t)):
        taxids.append(row['tax_id'])
    taxids = set(taxids)
    seqinfo = csv.DictReader(open(args.seq_info))
    seqinfo = (row for row in seqinfo if row['tax_id'] not in taxids)
    removed = set(row['version'] for row in seqinfo)
    modified = set(open_clean(args.modified))
    cache = set(open_clean(args.cache))
    cache -= removed - modified
    args.out.write('\n'.join(cache))


if __name__ == '__main__':
    main()
