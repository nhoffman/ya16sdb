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

    def confidence(s):
        if s['is_type']:
            s['confidence'] = 'type'
        elif s['is_refseq']:
            s['confidence'] = 'refseq'
        elif s['is_published']:
            s['confidence'] = 'published'
        else:
            s['confidence'] = 'direct'
        return s

    info = info.apply(confidence, axis='columns')
    info.to_feather(args.feather)


if __name__ == '__main__':
    main()
