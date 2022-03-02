#!/usr/bin/env python3
import argparse
import pandas


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        'feather', help='full seq_info file with description column')
    p.add_argument(
        'classifications',
        help='csv output from Moose classifier')
    args = p.parse_args()
    classifications = pandas.read_csv(
        args.classifications,
        dtype={'specimen': str, 'assignment': str},
        usecols=['specimen', 'assignment'])
    classifications = classifications.rename(
        columns={
            'specimen': 'seqname',
            'assignment': 'type_classification'})
    info = pandas.read_feather(args.feather)
    if 'type_classification' in info.columns:
        info = info.drop('type_classification', axis='columns')
    info = info.merge(classifications, how='left')
    info.to_feather(args.feather)


if __name__ == '__main__':
    main()
