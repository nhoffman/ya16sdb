#!/usr/bin/env python3
"""
Filters and aligns fasta records based on vsearch alignment and
presence of invalid sequence characters.
"""
import argparse
import sys

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
        '--unknowns',
        required=True,  # mimics ``taxit update_taxids``
        metavar='fasta',
        help=('fasta format output of sequences not '
              'aligned or with invalid sequence characters'))
    p.add_argument(
        '--out',
        default=sys.stdout,
        metavar='fasta',
        help='fasta output of sequences in forward orientation')

    args = p.parse_args()

    vsearch = (row.split('\t') for row in open(args.vsearch))
    vsearch = {row[0]: row[2] for row in vsearch if row[1] != '*'}

    for r in SeqIO.parse(args.fasta, 'fasta', Alphabet.IUPAC.ambiguous_dna):
        if r.id in vsearch and Alphabet._verify_alphabet(r.seq):
            if vsearch[r.id] == '-':
                r.seq = r.seq.reverse_complement()
            SeqIO.write(r, args.out, 'fasta')
        else:
            SeqIO.write(r, args.unknowns, 'fasta')


if __name__ == '__main__':
    main()
