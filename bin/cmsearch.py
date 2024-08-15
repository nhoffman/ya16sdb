#!/usr/bin/env python3
"""
Filters and aligns fasta records based on vsearch alignment and
presence of invalid sequence characters.
"""
import argparse
import os
import pandas
from Bio import SeqIO

CMSEARCH_COLS = [
    'seqname', 'accession', 'query name', 'query accession', 'mdl',
    '16s_start', '16s_stop', 'seq from', 'seq to', 'strand', 'trunc',
    'pass', 'gc', 'bias', 'score', 'E-value', 'inc', 'description of target']


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        'cmsearch',
        help='vsearch alignments')
    p.add_argument(
        'fasta',
        help='vsearch fasta file')
    p.add_argument('seq_info')
    p.add_argument(
        'out_fa',
        help='fasta output of sequences in forward orientation')
    p.add_argument(
        'out_info',
        help='fasta output of sequences in forward orientation')
    args = p.parse_args()
    dtypes = {
        'seqname': str, '16s_start': int, '16s_stop': int, 'strand': str}
    if not os.path.isfile(args.cmsearch):
        cmsearch = pandas.DataFrame(columns=dtypes.keys())
    else:
        cmsearch = pandas.read_csv(
            args.cmsearch,
            comment='#',
            dtype=dtypes,
            header=None,
            names=CMSEARCH_COLS,
            sep='\\s+',
            usecols=dtypes.keys())
    # alignments are already sorted by bitscore quality
    cmsearch = cmsearch.drop_duplicates(subset='seqname', keep='first')
    cmsearch = cmsearch.set_index('seqname')
    with open(args.out_fa, 'w') as out:
        for s in SeqIO.parse(args.fasta, 'fasta'):
            if s.id in cmsearch.index and cmsearch.loc[s.id]['strand'] == '-':
                s.seq = s.seq.reverse_complement()
            out.write('>{}\n{}\n'.format(s.description, s.seq))
    cmsearch = cmsearch.drop('strand', axis='columns')
    info = pandas.read_csv(args.seq_info, dtype=str, index_col='seqname')
    info = info.join(cmsearch)
    info['16s_start'] = info['16s_start'].fillna(0).astype(int)
    info['16s_stop'] = info['16s_stop'].fillna(0).astype(int)
    info.to_csv(args.out_info)


if __name__ == '__main__':
    main()
