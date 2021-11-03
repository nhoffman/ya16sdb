#!/usr/bin/env python3
"""
"""
import argparse
import pandas


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        'feather', help='full seq_info file with description column')
    p.add_argument(
        'vsearch',
        help='txt file of version numbers that are type strains')
    args = p.parse_args()
    vsearch = pandas.read_csv(
        args.vsearch,
        dtype={'seqname': str, 'match_seqname': str, 'match_pct': float},
        header=None,
        names=['seqname', 'match_seqname', 'match_pct'],
        sep='\t',
        usecols=['seqname', 'match_seqname', 'match_pct'])
    info = pandas.read_feather(args.feather)
    vsearch = vsearch.merge(
        info[['seqname', 'version', 'species_name', 'species']],
        left_on='match_seqname',
        right_on='seqname',
        suffixes=['', '_'])
    vsearch = vsearch.drop('seqname_', axis='columns')
    vsearch = vsearch.rename(
            columns={
                'version': 'match_version',
                'species_name': 'match_species',
                'species': 'match_species_id'})
    info = info.merge(vsearch, how='left', on='seqname')
    matching = info[info['species'] == info['match_species_id']]
    matching = matching.drop_duplicates(subset=['seqname'], keep='first')
    not_matching = info[~info['seqname'].isin(matching['seqname'])]
    not_matching = not_matching.drop_duplicates(
        subset=['seqname'], keep='first')
    info = pandas.concat([matching, not_matching]).reset_index(drop=True)
    info.to_feather(args.feather)


if __name__ == '__main__':
    main()
