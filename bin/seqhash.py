#!/usr/bin/env python3
"""
"""
import argparse
import hashlib
import pandas

from Bio import SeqIO


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('feather')
    p.add_argument('fasta')
    args = p.parse_args()
    info = pandas.read_feather(args.feather)
    hashes = {}
    with open(args.fasta) as fasta:
        for s in SeqIO.parse(fasta, 'fasta'):
            seq = str(s.seq).replace('\n', '').lower().encode('utf-8')
            hashes[s.id] = hashlib.sha1(seq).hexdigest()
    info['seqhash'] = info['seqname'].map(hashes)
    info.to_feather(args.feather)


if __name__ == '__main__':
    main()
