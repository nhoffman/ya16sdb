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
        args.taxonomy,
        dtype=str,
        usecols=['tax_id', 'tax_name', 'species', 'genus'])
    taxonomy = taxonomy.set_index('tax_id')
    info = pandas.read_feather(args.feather)
    for c in ['species', 'genus', 'species_name', 'genus_name']:
        if c in info.columns:
            info = info.drop(c, axis='columns')
    info = info.join(taxonomy[['species', 'genus']], on='tax_id')
    info = info.join(taxonomy['tax_name'], on='species')
    info = info.rename(columns={'tax_name': 'species_name'})
    info = info.join(taxonomy['tax_name'], on='genus')
    info = info.rename(columns={'tax_name': 'genus_name'})
    info.to_feather(args.feather)


if __name__ == '__main__':
    main()
