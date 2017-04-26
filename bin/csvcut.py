#!/usr/bin/env python3

import pandas
import sys
import argparse


def get_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('file',
                        help="The dir corresponding to the run")
    columns_parser = parser.add_mutually_exclusive_group(required=True)
    columns_parser.add_argument('--columns',
                                help=('Comma delimited list of column '
                                      'names or indices if --no-header'))
    columns_parser.add_argument('--not-columns',
                                help=('Comma delimited list of column '
                                      'names or indices if --no-header'))
    parser.add_argument('--out',
                        type=argparse.FileType('w'),
                        default=sys.stdout)
    return parser.parse_args()


def main():
    args = get_args()
    if args.columns:
        columns = args.columns.split(',')
        df = pandas.read_csv(
            args.file, dtype=str, na_filter=False, usecols=columns)
    else:
        df = pandas.read_csv(args.file, dtype=str, na_filter=False)
        df = df.drop(args.not_columns.tolist(), axis=1)
    df.to_csv(args.out, index=False)


if __name__ == '__main__':
    sys.exit(main())
