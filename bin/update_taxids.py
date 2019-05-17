#!/usr/bin/env python3
"""
"""
import argparse
import os
import pandas


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('seq_info')
    p.add_argument('a2t')
    p.add_argument('out')
    p.add_argument('unknowns', nargs='?')
    args = p.parse_args()
    if os.stat(args.seq_info).st_size > 0:
        info = pandas.read_csv(args.seq_info, dtype=str)
        columns = info.columns  # column order
        a2t = pandas.read_csv(args.a2t, dtype=str)
        if args.unknowns:
            unknowns = info[~info['version'].isin(a2t['version'])]
            unknowns.to_csv(args.unknowns, index=False)
        info = info[info['version'].isin(a2t['version'])]
        info = info.drop('tax_id', axis='columns')
        info = info.merge(a2t)
        info.to_csv(args.out, columns=columns, index=False)
    else:
        print('seq_info file is empty')


if __name__ == '__main__':
    main()
