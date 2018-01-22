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
    p.add_argument('do_not_trust')
    p.add_argument('fasta')
    p.add_argument('seq_info')
    p.add_argument('out_fa')
    p.add_argument('out_info')
    args = p.parse_args()

    do_not_trust = (a.strip() for a in open(args.do_not_trust))
    do_not_trust = (a for a in do_not_trust if a)
    do_not_trust = (a for a in do_not_trust if not a.startswith('#'))
    do_not_trust = set(do_not_trust)

    keep = []
    with open(args.seq_info) as in_info:
        for row in csv.DictReader(in_info):
            if row['accession'] not in do_not_trust:
                keep.append(row['seqname'])
    keep = set(keep)

    new_keep = []
    with open(args.fasta) as in_fasta, open(args.out_fa, 'w') as out_fa:
        for seq in SeqIO.parse(in_fasta, 'fasta'):
            if seq.id in keep:
                out_fa.write('>{}\n{}\n'.format(seq.description, seq.seq))
                new_keep.append(seq.id)
    new_keep = set(new_keep)

    with open(args.seq_info) as in_info, open(args.out_info, 'w') as out_info:
        in_info = csv.DictReader(in_info)
        out_info = csv.DictWriter(out_info, in_info.fieldnames)
        out_info.writeheader()
        out_info.writerows(i for i in in_info if i['seqname'] in new_keep)


if __name__ == '__main__':
    main()
