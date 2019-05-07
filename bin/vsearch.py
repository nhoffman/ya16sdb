#!/usr/bin/env python3
"""
Filters and aligns fasta records based on vsearch alignment and
presence of invalid sequence characters.
"""
import argparse
import pandas
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
        'seq_info',
        help='seq_info file')
    p.add_argument(
        'out_fa',
        help='fasta output of sequences in forward orientation')
    p.add_argument(
        'out_info',
        help='fasta output of sequences in forward orientation')
    p.add_argument(
        'unknown',
        help=('fasta format output of sequences not '
              'aligned or with invalid sequence characters'))
    p.add_argument(
        'unknown_accessions',
        help=('accession.versions of records that '
              'contain at least one non-16s alignment'))
    args = p.parse_args()
    vsearch = (row.strip().split('\t') for row in open(args.vsearch))
    vsearch = (row for row in vsearch if row[1] != '*')
    vsearch = sorted(vsearch, key=lambda x: x[2], reverse=True)
    # we want the plus (x[2]) aligns and they will be at the end
    vsearch = {row[0]: row[2] for row in vsearch}
    with open(args.out_fa, 'w') as out, open(args.unknown, 'w') as unknown:
        parsed = SeqIO.parse(args.fasta, 'fasta', Alphabet.IUPAC.ambiguous_dna)
        for s in parsed:
            if s.id in vsearch and Alphabet._verify_alphabet(s.seq):
                if vsearch[s.id] == '-':
                    s.seq = s.seq.reverse_complement()
                out.write('>{}\n{}\n'.format(s.description, s.seq))
            else:
                unknown.write('>{}\n{}\n'.format(s.description, s.seq))
    info = pandas.read_csv(args.seq_info, dtype=str)
    info = info[info['seqname'].isin(vsearch)]
    info.to_csv(args.out_info, index=False)
    unknown = info[~info['seqname'].isin(vsearch)]
    unknown['version'].to_csv(
        args.unknown_accessions, header=False, index=False)


if __name__ == '__main__':
    main()
