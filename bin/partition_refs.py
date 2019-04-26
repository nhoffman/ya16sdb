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
    p.add_argument('fasta')
    p.add_argument('feather')
    # outputs
    p.add_argument('out_fa')
    p.add_argument('out_info')
    p.add_argument(
        '--drop-duplicate-sequences',
        action='store_true',
        help='group by accession and drop rows with duplicate seqhashes')
    # filtering switches
    o = p.add_argument_group('OR options')
    o.add_argument(
        '--is_valid',
        action='store_true',
        help='filter for named (is_valid=true) records')
    o.add_argument(
        '--tax-ids',
        metavar='FILE',
        help='if tax_id is in FILE')
    a = p.add_argument_group('AND options')
    a.add_argument(
        '--is_species',
        action='store_true',
        help='filter for records with tax id in species column')
    a.add_argument(
        '--is_type',
        action='store_true',
        help='filter for type straing records')
    a.add_argument(
        '--min-length',
        metavar='',
        type=int,
        help='Minimum sequence length')
    a.add_argument(
        '--prop-ambig-cutoff',
        metavar='',
        type=float,
        help=('Maximum proportion of characters in '
              'sequence which may be ambiguous'))
    a.add_argument(
        '--species-cap',
        metavar='INT',
        type=int,
        help='group records by species taxid and accept ony top nth')
    return p


def main():
    args = build_parser().parse_args()

    info = pandas.read_feather(args.feather)
    info['keep'] = True

    if args.is_valid:
        info.loc[~info['is_valid'], 'keep'] = False

    if args.tax_ids:
        tax_ids = (i.strip() for i in open(args.tax_ids))
        tax_ids = set(i for i in tax_ids if i)
        info.loc[info['tax_id'].isin(tax_ids), 'keep'] = True

    # drop OR rows
    info = info[info['keep']].drop('keep', axis='columns')

    if args.min_length:
        info = info[info['length'] >= args.min_length]

    # raw prop_ambig filtering
    if args.prop_ambig_cutoff:
        info['prop_ambig'] = (
            info['ambig_count'] / info['length'])
        info = info[info['prop_ambig'] <= args.prop_ambig_cutoff]
        info = info.drop('prop_ambig', axis='columns')

    if args.is_type:
        info = info[info['is_type']]

    if args.is_species:
        info = info[~info['species'].isna()]

    if args.drop_duplicate_sequences:
        info = info[~info[['accession', 'seqhash']].duplicated(keep='first')]

    if args.species_cap:
        info = info.groupby(by='species', as_index=False)
        info = info.nth(list(range(args.species_cap)))

    seqs = {s.id: s.seq for s in SeqIO.parse(args.fasta, 'fasta')}

    drop = []
    with open(args.out_fa, 'w') as out_fa:
        for s in info['seqname'].values:
            if s in seqs:
                out_fa.write('>{}\n{}\n'.format(s, seqs[s]))
            else:
                drop.append(s)
    drop = set(drop)
    info[~info['seqname'].isin(drop)].to_csv(args.out_info, index=False)


if __name__ == '__main__':
    main()
