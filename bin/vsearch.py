#!/usr/bin/env python3
"""
Filters and aligns fasta records based on vsearch alignment and
presence of invalid sequence characters.
"""
import argparse
from Bio import SeqIO, Alphabet


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
    p.add_argument(
        'out',
        help='fasta output of sequences in forward orientation')
    p.add_argument(
        'unknowns',
        help=('fasta format output of sequences not '
              'aligned or with invalid sequence characters'))
    args = p.parse_args()
    vsearch = (row.strip().split('\t') for row in open(args.vsearch))
    vsearch = (row for row in vsearch if row[1] != '*')
    vsearch = sorted(vsearch, key=lambda x: x[2], reverse=True)
    # we want the plus (x[2]) aligns and they will be at the end
    vsearch = {row[0]: row[2] for row in vsearch}
    with open(args.out, 'w') as out, open(args.unknowns, 'w') as unknowns:
        seqs = SeqIO.parse(args.fasta, 'fasta', Alphabet.IUPAC.ambiguous_dna)
        for s in seqs:
            if s.id in vsearch and Alphabet._verify_alphabet(s.seq):
                if vsearch[s.id] == '-':
                    s.seq = s.seq.reverse_complement()
                out.write('>{}\n{}\n'.format(s.description, s.seq))
            else:
                unknowns.write('>{}\n{}\n'.format(s.description, s.seq))


if __name__ == '__main__':
    main()
