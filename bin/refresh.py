#!/usr/bin/env python3
"""
Append new records with old records

Sequences that have no_features or passed all filtering into the
seq_info_out are appended to the records_out file.

Sequences that did not pass the vsearch or tax_id update_taxids steps will be
not included in the records_out file and consequently will be re-downloaded

Sequence ids in the fasta file are compared with the seq_seq_info. If any
seq.id or seqname not present in either source all records by accession
associated will be dropped and re-downloaded at a later time.
"""

import argparse
import itertools
import pandas
import os

from Bio import SeqIO


def is_empty(fpath):
    return os.path.getsize(fpath) == 0


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument(
        'ncbi',
        help='list of all the latest ncbi versions')

    p.add_argument(
        'new_fasta',
        help='records file in fasta format')
    p.add_argument(
        'previous_fasta',
        help='records file in fasta format')

    p.add_argument(
        'new_seq_info',
        help='records file in csv format')
    p.add_argument(
        'previous_seq_info',
        help='records file in csv format')

    p.add_argument(
        'new_pubmed_seq_info',
        help='pubmed_seq_info file')
    p.add_argument(
        'previous_pubmed_seq_info',
        help='pubmed_seq_info file')

    p.add_argument(
        'new_references',
        help='references file')
    p.add_argument(
        'previous_references',
        help='references file')

    p.add_argument(
        'new_refseq_info',
        help='new refseq file')
    p.add_argument(
        'previous_refseq_info',
        help='previous refseq file')

    p.add_argument(
        'no_features',
        type=argparse.FileType('r'),
        help='list of versions downloaded with no 16s regions')
    p.add_argument(
        'previous_records',
        type=argparse.FileType('r'),
        help='list of previously downloaded record versions')

    p.add_argument(
        'fasta_out',
        type=argparse.FileType('w'),
        help='records file in same Bio.SeqIO format as input')
    p.add_argument(
        'seq_info_out',
        help='seq_info output')
    p.add_argument(
        'pubmed_seq_info_out',
        help='pubmed_info output')
    p.add_argument(
        'references_out',
        help='references output')
    p.add_argument(
        'refseq_info_out',
        help='refseqs output')
    p.add_argument(
        'records_out',
        type=argparse.FileType('w'),
        help='list of all version downloaded')

    args = p.parse_args()

    new_seq_info = pandas.read_csv(args.new_seq_info, dtype=str)

    if os.path.getsize(args.previous_seq_info) == 0:
        prev_seq_info = pandas.DataFrame(columns=new_seq_info.columns)
    else:
        prev_seq_info = pandas.read_csv(args.previous_seq_info, dtype=str)

    # remove previous seq_info in new seq_info with same accessions
    prev_seq_info = prev_seq_info[
        ~prev_seq_info['accession'].isin(new_seq_info['accession'])]

    # combine seq_info files and retain column order
    columns = prev_seq_info.columns.tolist()
    for c in new_seq_info.columns:
        if c not in columns:
            columns.append(c)
    seq_info = prev_seq_info.append(new_seq_info)[columns]

    # remove anything dropped from ncbi
    ncbi = pandas.read_csv(
        args.ncbi, dtype=str, squeeze=True, usecols=['version'])
    seq_info = seq_info[seq_info['version'].isin(ncbi)]

    # include with records_out
    current_records = set(v for v in seq_info['version'].tolist())

    # read refseqs and append to old refseq list
    new_refseq_info = pandas.read_csv(args.new_refseq_info, dtype=str)
    if os.path.getsize(args.previous_refseq_info) == 0:
        prev_refseqs = pandas.DataFrame(columns=new_refseq_info.columns)
    else:
        prev_refseqs = pandas.read_csv(args.previous_refseq_info, dtype=str)
    prev_refseqs = prev_refseqs[
        ~prev_refseqs['seqname'].isin(new_refseq_info['seqname'])]
    refseqs = new_refseq_info.append(prev_refseqs)
    refseqs = refseqs[refseqs['seqname'].isin(seq_info['seqname'])]
    assert(len(refseqs) == len(refseqs['seqname'].drop_duplicates()))
    refseqs.to_csv(args.refseq_info_out, index=False)

    # remove seq_info with refseq duplicates
    refseqs = refseqs[~refseqs['accession'].isnull()]  # avoid Null TypeError
    seq_info = seq_info[~seq_info['accession'].isin(refseqs['accession'])]

    # check for and raise exception if duplicate seqname rows
    if len(seq_info['seqname']) != len(seq_info['seqname'].drop_duplicates()):
        dups = seq_info.groupby(by='seqname').filter(lambda x: len(x) > 1)
        dups = dups['seqname'].drop_duplicates()
        err = ', '.join(dups.values)
        raise DuplicateSeqnameError(err)

    """
    deduplicate and write fasta
    """
    to_write = set(seq_info['seqname'].values)
    new_fa = SeqIO.parse(args.new_fasta, 'fasta')
    prev_fa = SeqIO.parse(args.previous_fasta, 'fasta')
    for r in itertools.chain(new_fa, prev_fa):
        if r.id in to_write:
            args.fasta_out.write('>{}\n{}\n'.format(r.description, r.seq))
            to_write.remove(r.id)

    # write seq_info
    seq_info.to_csv(args.seq_info_out, index=False, date_format='%d-%b-%Y')

    '''
    just like with the seq_seq_info, replace old pubmed_info with newly
    downloaded pubmed_info, by accession
    '''
    new_pubmed = pandas.read_csv(args.new_pubmed_seq_info, dtype=str)
    if os.path.getsize(args.previous_pubmed_seq_info) == 0:
        prev_pubmed = pandas.DataFrame(columns=new_pubmed.columns)
    else:
        prev_pubmed = pandas.read_csv(args.previous_pubmed_seq_info, dtype=str)
    prev_pubmed = prev_pubmed[~prev_pubmed.isin(prev_pubmed['accession'])]
    pubmed_info = prev_pubmed.append(new_pubmed).drop_duplicates()
    pubmed_info = pubmed_info[pubmed_info['version'].isin(seq_info['version'])]
    pubmed_info.to_csv(args.pubmed_seq_info_out, index=False)

    '''
    append new references. Use latest pubmed_info and
    only if present in pubmed_info dataframe
    '''
    new_refs = pandas.read_csv(args.new_references, dtype=str)
    if os.path.getsize(args.previous_references) == 0:
        prev_refs = pandas.DataFrame(columns=new_refs.columns)
    else:
        prev_refs = pandas.read_csv(args.previous_references, dtype=str)
    prev_refs = prev_refs[~prev_refs['pubmed_id'].isin(new_refs['pubmed_id'])]
    refs = prev_refs.append(new_refs)
    refs = refs[refs['pubmed_id'].isin(pubmed_info['pubmed_id'])]
    refs = refs.drop_duplicates()
    refs.to_csv(args.references_out, index=False)

    '''
    Filter the previous_records with only existing records in ncbi. That way
    if any previous records appear again in a future esearch
    we can re-download.
    '''
    previous_records = set(v.strip() for v in args.previous_records)
    previous_records = set(v for v in previous_records if v in ncbi)

    '''
    output downloaded versions we do not want to download again. This
    includes records with no features and records that passed filtering
    up to this point
    '''
    versions = set(seq_info['version'].tolist())
    versions |= set(n.strip() for n in args.no_features)
    versions |= current_records
    versions |= previous_records

    for v in sorted(versions):
        if v:
            args.records_out.write(v + '\n')


class DuplicateSeqnameError(Exception):
    """Raised when a seqname has more than one row of seq_info"""
    pass


if __name__ == '__main__':
    main()
