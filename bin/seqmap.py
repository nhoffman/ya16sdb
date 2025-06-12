#!/usr/bin/env python3
import argparse
import csv
import sys


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('seq_info', type=argparse.FileType('r')),
    parser.add_argument(
        '--out', type=argparse.FileType('w'), default=sys.stdout)
    return parser.parse_args()


def main():
    args = get_args()
    out = csv.DictWriter(
        args.out,
        delimiter=' ',
        fieldnames=['seqname', 'tax_id'],
        extrasaction='ignore')
    out.writerows(csv.DictReader(args.seq_info))


if __name__ == '__main__':
    sys.exit(main())
