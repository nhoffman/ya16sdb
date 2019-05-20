#!/usr/bin/env python3
"""
Add confidence column to feather file

Confidence values:

type - is_type is TRUE (may also have pubmed_id)
published - has pubmed_id and is_type is FALSE
direct - is_type is FALSE and no pubmed_id (the rest)
"""
import argparse
import pandas


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('feather')
    args = p.parse_args()
    info = pandas.read_feather(args.feather)
    if 'confidence' in info.columns:
        info = info.drop('confidence', axis='columns')
    info.loc[info['is_type'], 'confidence'] = 'type'
    is_refseq = info['confidence'].isna() & info['is_refseq']
    info.loc[is_refseq, 'confidence'] = 'refseq'
    is_published = info['confidence'].isna() & info['is_published']
    info.loc[is_published, 'confidence'] = 'published'
    info.loc[info['confidence'].isna(), 'confidence'] = 'direct'
    info.to_feather(args.feather)


if __name__ == '__main__':
    main()
