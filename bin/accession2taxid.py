#!/usr/bin/env python3
"""
create a filtered csv version of NCBI's accession2taxid
"""
import argparse
import csv
import gzip
import sys


def open_clean(fl):
    for f in list(fl):
        f = (row.strip() for row in open(f))
        f = (row for row in f if row)
        for row in f:
            yield row


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('accession2taxid')
    p.add_argument('accessions', nargs=argparse.REMAINDER)
    p.add_argument(
        '--out',
        default=sys.stdout,
        type=argparse.FileType('w'))
    args = p.parse_args()
    accessions = set(open_clean(args.accessions))
    a2t = csv.reader(gzip.open(args.accession2taxid, 'rt'), delimiter='\t')
    next(a2t)  # header
    acc, tax_id = 1, 2
    a2t = ((a[acc], a[tax_id]) for a in a2t if a[acc] in accessions)
    out = csv.writer(args.out)
    out.writerow(['version', 'tax_id'])
    out.writerows(a2t)


if __name__ == '__main__':
    main()
