#!/usr/bin/env python
"""
Create a seq_info.csv feather file optionally gzipped
"""
import argparse
import pandas
import sys
import ya16sdb


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('seq_info')
    parser.add_argument('out')
    args = parser.parse_args(arguments)
    seq_info = pandas.read_csv(args.seq_info, dtype=ya16sdb.DTYPES)
    seq_info.to_feather(args.out)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
