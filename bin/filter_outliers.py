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
    p.add_argument('details')
    args = p.parse_args()
    info = pandas.read_feather(args.feather)
    for c in ['centroid', 'cluster', 'dist', 'is_out',
              'x', 'y', 'filter_outliers']:
        if c in info.columns:
            info = info.drop(c, axis='columns')
    details = pandas.read_csv(
        args.details,
        dtype={
            'seqname': str,
            'tax_id': str,
            'centroid': str,
            'cluster': float,
            'dist': float,
            'is_out': bool,
            'species': str,
            'x': float,
            'y': float
            }
        )
    details['filter_outliers'] = True
    info = info.merge(details, how='left')
    info.loc[info['filter_outliers'].isna(), 'filter_outliers'] = False
    info.loc[info['is_out'].isna(), 'is_out'] = False
    info.to_feather(args.feather)


if __name__ == '__main__':
    main()
