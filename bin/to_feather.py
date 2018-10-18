#!/usr/bin/env python
"""
Create a feather file combining sequence details,
seq_info, taxonomy and publication annotations
"""
import argparse
import pandas
import sys

DTYPES = {
    'accession': str,
    'ambig_count': int,
    'authors': str,
    'centroid': str,
    'cluster': float,
    'comment': str,
    'consrtm': str,
    'description': str,
    'dist': float,
    'download_date': str,
    'genus': str,
    'is_out': bool,
    'is_type': bool,
    'isolate': str,
    'isolation_source': str,
    'journal': str,
    'keywords': str,
    'length': int,
    'match_pct': float,
    'match_seqname': str,
    'modified_date': str,
    'mol_type': str,
    'name': str,
    'organism': str,
    'pubmed_id': str,
    'seq_start': int,
    'seq_stop': int,
    'seqname': str,
    'source': str,
    'species': str,
    'strain': str,
    'tax_id': str,
    'title': str,
    'version': str,
    'version_num': str,
    'x': float,
    'y': float,
    }


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('details')
    parser.add_argument('seq_info')
    parser.add_argument('taxonomy')
    parser.add_argument('hits')
    parser.add_argument('pubmed_ids')
    parser.add_argument('--out', default=sys.stdout, help='[%(default)s]')
    args = parser.parse_args(arguments)

    details = pandas.read_csv(args.details, dtype=DTYPES)
    seq_info = pandas.read_csv(
        args.seq_info,
        dtype=DTYPES,
        parse_dates=['download_date', 'modified_date'],
        date_parser=lambda x: pandas.datetime.strptime(x, '%d-%b-%Y'))
    taxonomy = pandas.read_csv(
        args.taxonomy,
        dtype=DTYPES,
        usecols=['tax_id', 'tax_name', 'species', 'genus'])
    hits = pandas.read_table(
        args.hits,
        dtype=DTYPES,
        header=None,
        names=['seqname', 'match_seqname', 'match_pct'],
        usecols=['seqname', 'match_seqname', 'match_pct'])
    pubmed_ids = pandas.read_csv(
        args.pubmed_ids, dtype=DTYPES, usecols=['pubmed_id', 'version'])

    hits = hits.merge(
        seq_info[['seqname', 'tax_id', 'version']].rename(
            columns={'version': 'match_version'}),
        left_on='match_seqname',
        suffixes=['', '_'],
        right_on='seqname').drop(columns='seqname_', axis='columns')
    hits = hits.merge(taxonomy[['tax_id', 'species']])
    hits = hits.merge(taxonomy[['tax_id', 'tax_name']].rename(
        columns={'tax_name': 'match_species'}),
        left_on='species',
        right_on='tax_id',
        suffixes=['', '_'])
    hits = hits.drop(columns=['tax_id_', 'species'], axis='columns')

    details = details.merge(seq_info)
    details = details.merge(taxonomy)
    details = details.merge(hits, how='left')
    details = details.merge(pubmed_ids, how='left')

    taxonomy = taxonomy[['tax_id', 'tax_name']]

    details = details.merge(
        taxonomy.rename(columns={'tax_name': 'genus_name'}),
        how='left',
        left_on='genus',
        right_on='tax_id',
        suffixes=('', '_')).drop(columns='tax_id_', axis='columns')
    details = details.merge(
        taxonomy.rename(columns={'tax_name': 'species_name'}),
        left_on='species',
        right_on='tax_id',
        suffixes=('', '_')).drop(columns='tax_id_', axis='columns')

    # setup data for dist_pct v rank_order plot
    details = details.sort_values(by='dist')
    species_groups = details.groupby(by='species')['species']
    details['rank_order'] = species_groups.transform(lambda x: range(len(x)))
    details['dist_pct'] = details['dist'].apply(
        lambda x: '{:.2f}'.format(x * 100))

    details.reset_index(drop=True).to_feather(args.out)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
