#!/usr/bin/env python3
"""
Add is_type column based on presence in type strain file
"""
import argparse
import pandas


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('feather')
    p.add_argument('refseqs')
    args = p.parse_args()
    refseqs = pandas.read_csv(
        args.refseqs, squeeze=True, dtype=str, usecols=['seqname'])
    info = pandas.read_feather(args.feather)
    info['is_refseq'] = False
    info.loc[info['seqname'].isin(refseqs), 'is_refseq'] = True
    info.to_feather(args.feather)


if __name__ == '__main__':
    main()
