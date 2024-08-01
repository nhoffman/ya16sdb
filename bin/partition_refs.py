#!/usr/bin/env python3
"""Splits sequence files into partions and optionally
filters by length and percent ambiguity.
"""
import argparse
import pandas
import sys
import ya16sdb

from Bio import SeqIO


def build_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    # inputs
    p.add_argument('fasta')
    p.add_argument('info')
    # outputs
    p.add_argument('out_fa')
    p.add_argument('out_info')
    p.add_argument(
        '--drop-duplicate-sequences',
        action='store_true',
        help='group by accession and drop rows with duplicate seqhashes')
    f = p.add_argument_group('filtering options')
    f.add_argument(
        '--is_valid',
        action='store_true',
        help='filter for named (is_valid=true) records')
    f.add_argument(
        '--do_not_trust',
        help='drop these sequences or tax_ids')
    f.add_argument(
        '--drop-noaligns',
        action='store_true',
        help=('drop sequences that did not align '
              'to the cmsearch covariance model'))
    f.add_argument(
        '--inliers',
        action='store_true',
        help='choose sequences that were filtered and designated inliers')
    f.add_argument(
        '--is_species',
        action='store_true',
        help='filter for records with tax id in species column')
    f.add_argument(
        '--is_type',
        action='store_true',
        help='filter for type straing records')
    f.add_argument(
        '--min-length',
        metavar='',
        type=int,
        help='Minimum sequence length')
    f.add_argument(
        '--prop-ambig-cutoff',
        metavar='',
        type=float,
        help=('Maximum proportion of characters in '
              'sequence which may be ambiguous'))
    f.add_argument(
        '--species-cap',
        metavar='INT',
        type=int,
        help='group records by species taxid and accept ony top nth')
    f.add_argument(
        '--trusted',
        help='trusted record accessions and versions')
    return p


def main():
    args = build_parser().parse_args()

    if args.info == '-':
        info = pandas.read_csv(sys.stdin, dtype=ya16sdb.DTYPES)
    elif args.info.endswith('.csv'):
        info = pandas.read_csv(args.info, dtype=ya16sdb.DTYPES)
    elif args.info.endswith('.feather'):
        info = pandas.read_feather(args.info)
    else:
        raise TypeError('file extension not yet supported ' + args.info)

    if args.trusted:
        # all type strains are trusted - GL #79
        trusted = info[info['is_type']]
        recs = (i.strip() for i in open(args.trusted) if not i.startswith('#'))
        recs = set(i for i in recs if i)
        # add special trusted recs after type strains
        # because they are preferred
        trusted = pandas.concat([trusted, info[
            (info['version'].isin(recs)) |
            (info['accession'].isin(recs)) |
            (info['tax_id'].isin(recs)) |
            (info['seqname'].isin(recs))]])

    if args.is_valid:
        info = info[info['is_valid']]

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

    if args.inliers:
        info = info[info['filter_outliers'] & ~info['is_out']]

    if args.drop_noaligns:
        info = info[info['16s_stop'] != 0]

    if args.drop_duplicate_sequences:
        old_loci = info[
            ~info['locus_tag'].isna() &
            ~info['assembly_genbank'].isna() &
            info['locus_tag'].isin(info['old_locus_tag'])]
        old_loci = old_loci[['seqname', 'locus_tag', 'assembly_genbank']]
        old_loci = old_loci.merge(
            info[['old_locus_tag', 'assembly_genbank']],
            left_on=['locus_tag', 'assembly_genbank'],
            right_on=['old_locus_tag', 'assembly_genbank'])
        # remove old locus tags
        info = info[~info['seqname'].isin(old_loci['seqname'])]
        info = info.drop_duplicates(
            subset=['accession', 'seqhash'], keep='first')

    if args.trusted:
        # trusted sequences are preferred and appended
        # last to ensure they pass species_cap
        trusted = trusted[~trusted['seqname'].isin(info['seqname'])]
        info = pandas.concat([info, trusted])

    # apply args.do_not_trust after adding args.trusted in case a
    # type strain needs to be removed
    if args.do_not_trust:
        dnt = (i for i in open(args.do_not_trust) if not i.startswith('#'))
        dnt = (i.strip() for i in dnt)
        dnt = set(i for i in dnt if i)
        info = info[~info['tax_id'].isin(dnt)]
        info = info[~info['accession'].isin(dnt)]
        info = info[~info['version'].isin(dnt)]
        info = info[~info['seqname'].isin(dnt)]

    if args.species_cap:
        # remember seqs are sorted in reverse preferred order
        info = info.groupby(by='species').tail(args.species_cap)

    seqs = {s.id: s.seq for s in SeqIO.parse(args.fasta, 'fasta')}

    drop = []
    with open(args.out_fa, 'w') as out_fa:
        for s in info['seqname'].values:
            if s in seqs:
                out_fa.write('>{}\n{}\n'.format(s, seqs[s]))
            else:
                drop.append(s)

    info = info[~info['seqname'].isin(set(drop))]

    info.to_csv(args.out_info, index=False, date_format='%d-%b-%Y')


if __name__ == '__main__':
    main()
