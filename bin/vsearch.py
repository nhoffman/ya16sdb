#!/usr/bin/env python3
"""
Filters and aligns fasta records based on vsearch alignment and
presence of invalid sequence characters.
"""
import argparse
import pandas
from Bio import SeqIO


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        'vsearch',
        help='vsearch alignments')
    p.add_argument(
        'fasta',
        help='vsearch fasta file')
    p.add_argument('seq_info')
    p.add_argument('prev_unknowns')
    p.add_argument(
        'out_fa',
        help='fasta output of sequences in forward orientation')
    p.add_argument(
        'out_info',
        help='fasta output of sequences in forward orientation')
    p.add_argument(
        'unknowns',
        help=('fasta format output of sequences not '
              'aligned or with invalid sequence characters'))
    p.add_argument(
        'unknown_versions',
        help=('accession.versions of records that '
              'contain at least one non-16s alignment'))
    args = p.parse_args()
    vsearch = (row.strip().split('\t') for row in open(args.vsearch))
    vsearch = (row for row in vsearch if row[1] != '*')
    vsearch = sorted(vsearch, key=lambda x: x[2], reverse=True)
    # we want the plus aligns at the end of x[2]
    vsearch = {row[0]: row[2] for row in vsearch}
    with open(args.out_fa, 'w') as out, open(args.unknowns, 'w') as unknowns:
        parsed = SeqIO.parse(args.fasta, 'fasta')
        for s in parsed:
            if s.id in vsearch:
                if vsearch[s.id] == '-':
                    s.seq = s.seq.reverse_complement()
                out.write('>{}\n{}\n'.format(s.description, s.seq))
            else:
                unknowns.write('>{}\n{}\n'.format(s.description, s.seq))
    info = pandas.read_csv(args.seq_info, dtype=str)
    unknowns = info[~info['seqname'].isin(vsearch)]['version'].tolist()
    info = info[info['seqname'].isin(vsearch)]
    info.to_csv(args.out_info, index=False)
    with open(args.prev_unknowns, 'r') as unknowns_in:
        unknowns_in = (u.strip() for u in unknowns_in)
        unknowns_in = (u for u in unknowns_in if u)
        unknowns.extend(unknowns_in)
    with open(args.unknown_versions, 'w') as unknowns_out:
        for u in set(unknowns):
            unknowns_out.write(u + '\n')


if __name__ == '__main__':
    main()
