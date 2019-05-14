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
    vsearch = vsearch.merge(info[['seqname', 'version', 'species_name']])
    vsearch = vsearch.rename(
        columns={'version': 'match_version', 'species_name': 'match_species'})
    info = info.merge(vsearch)
    info.to_feather(args.feather)


if __name__ == '__main__':
    main()
