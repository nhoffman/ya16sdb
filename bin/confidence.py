#!/usr/bin/env python3
"""
Add confidence column to seq_info file

Confidence values:

type - is_type is TRUE (may also have pubmed_id)
published - has pubmed_id and is_type is FALSE
direct - is_type is FALSE and no pubmed_id (the rest)
"""
import argparse
import csv
import sys


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('seq_info', type=argparse.FileType('r'))
    p.add_argument('pubmed_ids', type=argparse.FileType('r'))
    p.add_argument('--out', type=argparse.FileType('w'), default=sys.stdout)
    args = p.parse_args()
    pids = set(r['version'] for r in csv.DictReader(args.pubmed_ids))
    info = csv.DictReader(args.seq_info)
    writer = csv.DictWriter(
        args.out, fieldnames=info.fieldnames + ['confidence'])
    writer.writeheader()
    for row in info:
        if row['is_type'] == 'True':
            row['confidence'] = 'type'
        elif row['version'] in pids and row['is_type'] == 'False':
            row['confidence'] = 'published'
        else:
            row['confidence'] = 'direct'
        writer.writerow(row)


if __name__ == '__main__':
    main()
