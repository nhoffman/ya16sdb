#!/usr/bin/env python3
"""
Filters and aligns fasta records based on vsearch alignment and
presence of invalid sequence characters.
"""
import argparse
import csv

from Bio import SeqIO


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument('do_not_trust', type=argparse.FileType('r'))
    p.add_argument('annotations', type=argparse.FileType('r'))

    p.add_argument('fasta')
    p.add_argument('seq_info', type=argparse.FileType('r'))
    p.add_argument('details', type=argparse.FileType('r'))

    p.add_argument('out_fa', type=argparse.FileType('w'))
    p.add_argument('out_info', type=argparse.FileType('w'))
    p.add_argument('out_details', type=argparse.FileType('w'))

    args = p.parse_args()

    do_not_trust = (a.strip() for a in args.do_not_trust)
    do_not_trust = (a for a in do_not_trust if a)
    do_not_trust = (a for a in do_not_trust if not a.startswith('#'))
    do_not_trust = set(do_not_trust)

    keep = []
    for row in csv.DictReader(args.annotations):
        if row['accession'] not in do_not_trust:
            keep.append(row['seqname'])
    keep = set(keep)

    for seq in SeqIO.parse(args.fasta, 'fasta'):
        if seq.id in keep:
            args.out_fa.write('>{}\n{}\n'.format(seq.description, seq.seq))

    seq_info = csv.DictReader(args.seq_info)
    out_info = csv.DictWriter(args.out_info, seq_info.fieldnames)
    out_info.writeheader()
    out_info.writerows(i for i in seq_info if i['seqname'] in keep)

    details = csv.DictReader(args.details)
    out_details = csv.DictWriter(args.out_details, details.fieldnames)
    out_details.writeheader()
    out_details.writerows(d for d in details if d['seqname'] in keep)


if __name__ == '__main__':
    main()
