#!/usr/bin/env python

"""
TODO: use taxit to calculae lineages here or do it in pipeline

eg:

datadir=/mnt/disk2/molmicro/common/ncbi/16s/output/20180214
singularity exec -B $datadir --pwd $datadir/dedup/1200bp/named /molmicro/common/singularity/taxtastic-0.8.5-singularity2.4.img taxit lineage_table taxonomy.csv seq_info.csv -c lineages.csv

Then:

export named=$datadir/dedup/1200bp/named
bin/prepare_plot_data.py \
    --details $named/filtered/details_out.csv \
    --seq-info $named/seq_info.csv \
    --lineages $named/lineages.csv \
    --hits $named/vsearch.tsv \
    -o $named/filtered_details.feather
"""

import sys
from os import path
import argparse

import pandas as pd

def main(arguments):

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--details', help='[default %(default)s]', default='details_out.csv')
    parser.add_argument('--seq-info', help='[default %(default)s]', default='seq_info.csv')
    parser.add_argument('--lineages', help='[default %(default)s]', default='lineages.csv')
    parser.add_argument('--hits', help='[default %(default)s]', default='vsearch.tsv')
    parser.add_argument('-o', '--outfile', help='[default %(default)s]', default='filtered_details.feather')

    args = parser.parse_args(arguments)

    outfile = args.outfile

    details = pd.read_csv(
        args.details,
        dtype={
            'cluster': str,
            'tax_id': str,
            'species': str,
            'is_out': str,
        })

    seq_info = pd.read_csv(
        args.seq_info,
        usecols=['seqname', 'accession', 'description', 'is_type'],
    )

    details = pd.merge(details, seq_info, how='left', on=['seqname'])

    lineages = pd.read_csv(
        args.lineages,
        usecols=['seqname', 'species'],
    )
    lineages = lineages.rename(columns={'species': 'tax_name'})

    hits = pd.read_table(
        args.hits,
        usecols=[0, 1, 2],
        header=None,
        names=['seqname', 'match_seqname', 'match_pct'],
    )
    hits = pd.merge(
        hits, lineages.rename(
            columns={'seqname': 'match_seqname', 'tax_name': 'match_species'}
        ),
        how='left', on=['match_seqname']
    )

    df = pd.merge(details, lineages, how='left', on=['seqname'])
    df = pd.merge(df, hits, how='left', on=['seqname'])

    print(df.head(1))
    print(pd.isnull(df).describe())
    df.to_feather(outfile)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

# read it back in:
# time python -c "import pandas as pd; df = pd.read_feather('output/filtered_details.feather'); print(df.head())"
