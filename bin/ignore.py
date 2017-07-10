#!/usr/bin/env python3
"""
Remove records in accession list or past as a RefSeq under the original column
"""
import pandas
import sys
import argparse


def get_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'ncbi', help='list of ncbi record versions')
    parser.add_argument(
        'versions', help='list of ncbi records by versions')
    parser.add_argument(
        'accessions',
        help='list of records by accession')
    parser.add_argument(
        '--refseq-originals',
        help='list of accessions replaced by RefSeqs with column original')
    parser.add_argument(
        '--out',
        default=sys.stdout,
        help='output list of sequences')
    return parser.parse_args()


def main():
    args = get_args()
    ncbi = pandas.read_csv(
        args.ncbi,
        header=None,
        names=['version'],
        squeeze=True,
        dtype=str)
    versions = pandas.read_csv(
        args.versions,
        header=None,
        names=['version'],
        squeeze=True,
        dtype=str)
    accessions = pandas.read_csv(
        args.accessions,
        header=None,
        comment='#',
        names=['accession'],
        squeeze=True,
        dtype=str)
    original_accessions = pandas.read_csv(
        args.refseq_originals,
        usecols=['original'],
        squeeze=True,
        dtype=str)

    # drop records
    ncbi = ncbi[~ncbi.isin(versions)]
    ncbi = ncbi[~ncbi.str[:-2].isin(accessions)]
    ncbi = ncbi[~ncbi.str[:-2].isin(original_accessions)]

    # output
    ncbi.to_csv(args.out, index=False)


if __name__ == '__main__':
    sys.exit(main())
