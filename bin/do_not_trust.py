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

    p.add_argument('fasta')
    p.add_argument('annotations', type=argparse.FileType('r'))
    p.add_argument('details', type=argparse.FileType('r'))

    p.add_argument('out_fa', type=argparse.FileType('w'))
    p.add_argument('out_info', type=argparse.FileType('w'))
    p.add_argument('out_details', type=argparse.FileType('w'))

    args = p.parse_args()

    do_not_trust = set(args.do_not_trust)

    keep = []

    annos = csv.DictReader(args.annotations)
    out_annos = csv.DictWriter(args.out_info, annos.fieldnames)
    out_annos.writeheader()
    for row in annos:
        if row['accession'] not in do_not_trust:
            out_annos.writerow(row)
            keep.append(row['seqname'])

    keep = set(keep)

    for seq in SeqIO.parse(args.fasta, 'fasta'):
        if seq.id in keep:
            args.out_fa.write('>{}\n{}\n'.format(seq.description, seq.seq))

    details = csv.DictReader(args.details)
    out_details = csv.DictWriter(args.out_details, details.fieldnames)
    out_details.writeheader()
    out_details.writerows(d for d in details if d['seqname'] in keep)


if __name__ == '__main__':
    main()
