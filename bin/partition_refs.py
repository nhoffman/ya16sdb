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
        'seqinfo',
        metavar='CSV',
        help="""Sequence metadata""")

    # outputs
    p.add_argument(
        '--out-info',
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
        '--types',
        action='store_true',
        help='select is_type == TRUE')
    flt.add_argument(
        '--tax-ids',
        help='column tax_id')

    return p


def main():
    args = build_parser().parse_args()

    seqs = SeqIO.parse(args.fasta, 'fasta')
    dtype = {'tax_id': str, 'gi': str, 'species': str,
             'ambig_count': float, 'is_type': bool,
             'length': int, 'description': str, 'accession': str,
             'seq_start': int, 'seq_stop': int, 'gi': str, 'name': str,
             'organism': str, 'date': str, 'keywords': str, 'seqname': str,
             'source': str, 'version': str}
    seq_info = pandas.read_csv(args.seqinfo, dtype=dtype).set_index('seqname')

    # raw min_length filtering
    if args.min_length:
        seq_info = seq_info[seq_info['length'] > args.min_length]

    # raw prop_ambig filtering
    if args.prop_ambig_cutoff:
        seq_info['prop_ambig'] = seq_info['ambig_count'] / seq_info['length']
        seq_info = seq_info[seq_info['prop_ambig'] < args.prop_ambig_cutoff]
        seq_info = seq_info.drop('prop_ambig', axis=1)

    if args.types:
        seq_info = seq_info[seq_info['is_type']]

    if args.tax_ids:
        tax_ids = pandas.read_csv(
            args.tax_ids, usecols=['tax_id'], squeeze=True, dtype=str)
        seq_info = seq_info[seq_info['tax_id'].isin(tax_ids)]

    if args.out_fa:
        seqs = (s for s in seqs if s.id in seq_info.index)
        SeqIO.write(seqs, args.out_fa, 'fasta')

    if args.out_info:
        seq_info.to_csv(args.out_info)


if __name__ == '__main__':
    main()
