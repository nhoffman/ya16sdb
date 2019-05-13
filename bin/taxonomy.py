#!/usr/bin/env python3
"""
"""
import argparse
import pandas


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('feather')
    p.add_argument('taxonomy')
    args = p.parse_args()
    taxonomy = pandas.read_csv(
        args.taxonomy, dtype=str, usecols=['tax_id', 'species'])
    info = pandas.read_feather(args.feather)
    info = info.merge(taxonomy, how='left')
    info.to_feather(args.feather)


if __name__ == '__main__':
    main()
