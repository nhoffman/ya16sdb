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
    # inputs
    p.add_argument(
        'fasta',
        metavar='FASTA',
        help="""sequence file""",
        type=argparse.FileType('r'))
    p.add_argument(
        'annotations',
        metavar='CSV',
        help="""Sequence metadata""")
    # outputs
    p.add_argument(
        'out_fa',
        type=argparse.FileType('w'),
        help='fasta out')
    p.add_argument(
        'out_annotations',
        metavar='CSV',
        type=argparse.FileType('w'),
        help='seq info out')
    # filtering switches
    flt = p.add_argument_group('filtering options')
    flt.add_argument(
        '--types',
        type=argparse.FileType('r'),
        help=('text file of accessions to be marked '
              'is_type=True in seq_info file'))
    flt.add_argument(
        '--is_type',
        action='store_true',
        help='filter for records is_type=True')
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

    if args.tax_ids:
        tax_ids = (i.strip() for i in open(args.tax_ids))
        tax_ids = set(i for i in tax_ids if i)
        annotations = annotations[annotations['tax_id'].isin(tax_ids)]

    if args.types:
        types = (t.strip() for t in args.types)
        types = set(t for t in types if t)
        annotations['is_type'] = annotations['version'].isin(types)

    if args.is_type:
        annotations = annotations[annotations['is_type'] == 'True']

    if args.out_fa:
        for s in SeqIO.parse(args.fasta, 'fasta'):
            if s.id in annotations.index:
                args.out_fa.write('>{}\n{}\n'.format(s.description, s.seq))

    if args.out_annotations:
        annotations.to_csv(args.out_annotations)


if __name__ == '__main__':
    main()
