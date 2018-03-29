#!/usr/bin/env python
'''
merge two csv files using Pandas

https://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.merge.html
'''
import argparse
import logging
import pandas
import sys

log = logging.getLogger(__name__)


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        'left',
        help="csv file to append column")
    parser.add_argument(
        'right',
        help="csv file to append file from")
    parser.add_argument(
        '--left-columns',
        type=lambda x: x.split(','),
        metavar='',
        help="column delimited list of columns to include")
    parser.add_argument(
        '--right-columns',
        type=lambda x: x.split(','),
        metavar='',
        help="column delimited list of columns to include")
    parser.add_argument(
        '--on',
        metavar='',
        type=lambda x: x.split(','),
        help="comma delimited list of columns to merge on")
    parser.add_argument(
        '--left-on',
        metavar='',
        type=lambda x: x.split(','),
        help="comma delimited list of columns to merge on")
    parser.add_argument(
        '--right-on',
        metavar='',
        type=lambda x: x.split(','),
        help="comma delimited list of columns to merge on")
    parser.add_argument(
        '--how',
        metavar='',
        choices=['left', 'right', 'outer', 'inner'],
        default='inner',
        help="")
    parser.add_argument(
        '--out',
        metavar='',
        default=sys.stdout,
        help="csv output")
    args = parser.parse_args(arguments)

    left = pandas.read_csv(
        args.left,
        usecols=args.left_columns,
        dtype=str,
        compression='infer')
    right = pandas.read_csv(
        args.right,
        dtype=str,
        usecols=args.right_columns,
        compression='infer')

    merged = left.merge(
        right,
        how=args.how,
        on=args.on,
        left_on=args.left_on,
        right_on=args.right_on)

    ext = None if args.out is sys.stdout else args.out.split('.')[-1]
    if ext in {'gzip', 'bz2', 'xz'}:
        merged.to_csv(args.out, index=False, compression=ext)
    else:
        merged.to_csv(args.out, index=False)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
