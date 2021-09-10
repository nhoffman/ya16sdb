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
    p.add_argument('a2t')
    p.add_argument('ncbi')
    p.add_argument('accessions')
    p.add_argument('seq_info')
    p.add_argument('out_ncbi')
    p.add_argument('out_download')
    args = p.parse_args()
    a2t = pandas.read_csv(args.a2t, dtype=str)
    ncbi = pandas.read_csv(
        args.ncbi,
        dtype=str,
        header=None,
        names=['version'])
    ncbi = ncbi[ncbi['version'].isin(a2t['version'])]
    ncbi.to_csv(args.out_ncbi, header=None, index=False)
    accessions = pandas.read_csv(
        args.accessions,
        dtype=str,
        header=None,
        names=['version'])
    accessions = accessions[accessions['version'].isin(a2t['version'])]
    if os.stat(args.seq_info).st_size > 0:
        seq_info = pandas.read_csv(
            args.seq_info, dtype=str, usecols=['version', 'tax_id'])
        same = seq_info.merge(a2t)
        # only add back sequences with updated taxids
        seq_info = seq_info[(~seq_info['version'].isin(same['version'])) &
                            (seq_info['version'].isin(a2t['version']))]
        accessions = accessions.append(seq_info[['version']])
    accessions = accessions.drop_duplicates()
    accessions.to_csv(args.out_download, header=None, index=False)


if __name__ == '__main__':
    main()
