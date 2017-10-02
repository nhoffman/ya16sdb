#!/usr/bin/env python3
"""Splits sequence files into partions and optionally
filters by length and percent ambiguity.
"""
import argparse
import pandas

from Bio import SeqIO


def build_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument(
        'annotations',
        metavar='CSV',
        help="""Sequence metadata""")
    # inputs
    p.add_argument(
        '--fasta',
        metavar='FASTA',
        help="""sequence file""",
        type=argparse.FileType('r'))

    # outputs
    p.add_argument(
        '--out-annotations',
        metavar='CSV',
        type=argparse.FileType('w'),
        help='seq info out')
    p.add_argument(
        '--out-fa',
        type=argparse.FileType('w'),
        help='fasta out')
    # filtering switches
    flt = p.add_argument_group('filtering options')
    flt.add_argument(
        '-a', '--prop-ambig-cutoff',
        type=float,
        help=('Maximum proportion of characters in '
              'sequence which may be ambiguous'))
    flt.add_argument(
        '--min-length',
        type=int,
        help='Minimum sequence length')
    flt.add_argument(
        '--tax-ids',
        help='column tax_id')

    return p


def main():
    args = build_parser().parse_args()

    annotations = pandas.read_csv(args.annotations, dtype=str)
    annotations = annotations.set_index('seqname')

    # raw min_length filtering
    if args.min_length:
        annotations['length'] = annotations['length'].astype(int)
        annotations = annotations[annotations['length'] > args.min_length]

    # raw prop_ambig filtering
    if args.prop_ambig_cutoff:
        annotations['ambig_count'] = annotations['ambig_count'].astype(int)
        annotations['length'] = annotations['length'].astype(int)
        annotations['prop_ambig'] = (
            annotations['ambig_count'] / annotations['length'])
        annotations = (
            annotations[annotations['prop_ambig'] < args.prop_ambig_cutoff])
        annotations = annotations.drop('prop_ambig', axis=1)

    if args.types:
        annotations = annotations[annotations['is_type'].str.lower() == 'true']

    if args.tax_ids:
        tax_ids = pandas.read_csv(
            args.tax_ids, usecols=['tax_id'], squeeze=True, dtype=str)
        annotations = annotations[annotations['tax_id'].isin(tax_ids)]

    if args.out_fa:
        seqs = SeqIO.parse(args.fasta, 'fasta')
        seqs = (s for s in seqs if s.id in annotations.index)
        SeqIO.write(seqs, args.out_fa, 'fasta')

    if args.out_annotations:
        annotations.to_csv(args.out_annotations)


if __name__ == '__main__':
    main()
