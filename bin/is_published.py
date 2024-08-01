#!/usr/bin/env python3
"""
Add is_published column
"""
import argparse
import pandas


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('feather')
    p.add_argument('pubmed_ids')
    args = p.parse_args()
    published = pandas.read_csv(
        args.pubmed_ids, dtype=str, usecols=['version']).squeeze()
    info = pandas.read_feather(args.feather)
    info['is_published'] = False
    info.loc[info['version'].isin(published), 'is_published'] = True
    info.to_feather(args.feather)


if __name__ == '__main__':
    main()
