import sys
from os import path

import pandas as pd

# TODO: implement a proper CLI and replace hard-coded input files with
# command line arguments.

try:
    topdir = sys.argv[1]
except IndexError:
    sys.exit('usage {} <path> (output of ya16sdb)'.format(sys.argv[0]))

named = path.join(topdir, 'dedup', '1200bp', 'named')

details = pd.read_csv(
    path.join(named, 'filtered', 'details_out.csv'),
    dtype={
        'cluster': str,
        'tax_id': str,
        'species': str,
        'is_out': str,
    })

seq_info = pd.read_csv(
    path.join(named, 'seq_info.csv'),
    usecols=['seqname', 'accession', 'description', 'is_type'],
)

details = pd.merge(details, seq_info, how='left', on=['seqname'])

# TODO: use taxit to calculae lineages here or do it in pipeline
# eg:

# cd ~/src/taxtastic
# datadir=~/working/2018-01-11-bokeh-outlier-plots/20180112
# named=$datadir/dedup/1200bp/named
# ./taxit.py lineage_table $named/taxonomy.csv $named/seq_info.csv -c $named/lineages.csv

lineages = pd.read_csv(
    path.join(named, 'lineages.csv'),
    usecols=['seqname', 'species'],
)
lineages = lineages.rename(columns={'species': 'tax_name'})

hits = pd.read_table(
    path.join(named, 'vsearch.tsv'),
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
df.to_feather('output/filtered_details.feather')

# read it back in:
# time python -c "import pandas as pd; df = pd.read_feather('output/filtered_details.feather'); print(df.head())"
