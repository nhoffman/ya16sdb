#!/usr/bin/env python3
"""
Append new records with old records

Sequences that have no_features or passed all filtering into the
annotations_out are appended to the versions_out file.

Sequences that did not pass the vsearch or tax_id update_taxids steps will be
not included in the versions_out file and consequently will be re-downloaded

Sequence ids in the fasta file are compared with the seq_info. If any seq.id
or seqname not present in either source all records by accession associated
will be dropped and re-downloaded at a later time.
"""

import argparse
import pandas

from Bio import SeqIO


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
        'new_annos',
        help='records file in csv format')
    p.add_argument(
        'previous_annotations',
        help='records file in csv format')

    p.add_argument(
        'new_pubmed_info',
        help='pubmed_info file')
    p.add_argument(
        'previous_pubmed_info',
        help='pubmed_info file')

    p.add_argument(
        'new_references',
        help='references file')
    p.add_argument(
        'previous_references',
        help='references file')

    p.add_argument(
        'no_features',
        type=argparse.FileType('r'),
        help='list of versions downloaded with no 16s regions')
    p.add_argument(
        'previous_versions',
        type=argparse.FileType('r'),
        help='list of previously downloaded record versions')

    p.add_argument(
        'fasta_out',
        type=argparse.FileType('w'),
        help='records file in same Bio.SeqIO format as input')
    p.add_argument(
        'annotations_out',
        help='seq_info file to output')
    p.add_argument(
        'pubmed_info_out',
        help='pubmed_ids file to output')
    p.add_argument(
        'references_out',
        help='references file to output')
    p.add_argument(
        'versions_out',
        type=argparse.FileType('w'),
        help='list of all version downloaded')

    args = p.parse_args()

    new_annos = pandas.read_csv(args.new_annos, dtype=str)

    # sync fasta and annotations
    new_seqs = (s.id for s in SeqIO.parse(args.new_fasta, 'fasta'))
    new_seqs = (i for i in new_seqs if i in new_annos['seqnames'].values)
    new_annos = new_annos[new_annos['seqname'].isin(new_seqs)]

    # deduplicate records by version and modified date
    new_annos = new_annos.sort_values(
        by=['version_num', 'modified_date', 'download_date'],
        ascending=True)
    new_annos = new_annos.drop_duplicates(subset='seqname', keep='last')

    # read in prev_info ignoring any accessions in the new data set
    prev_info = pandas.read_csv(args.previous_annotations, dtype=str)
    prev_info = prev_info[~prev_info['accession'].isin(new_annos['accession'])]

    # combine info files and retain column order
    columns = prev_info.columns.tolist()
    for c in new_annos.columns:
        if c not in columns:
            columns.append(c)
    info = prev_info.append(new_annos)[columns]

    assert(len(info['seqname']) == len(info['seqname'].drop_duplicates()))

    # remove anything dropped in ncbi
    ncbi = set(v.strip() for v in open(args.ncbi))
    info = info[info['version'].isin(ncbi)]

    # remove original records with refseq equivalent
    info = info[~info['accession'].isin(info['original'])]

    """
    filter accession records if any record that
    is inconsistent or has illegal chars
    """
    info_names = set(info['seqname'].tolist())
    keep = set()
    drop = set()

    for r in SeqIO.parse(args.new_fasta, 'fasta'):
        if r.id in info_names:
            keep.add(r.id)
        else:
            drop.add(r.id)
    for r in SeqIO.parse(args.previous_fasta, 'fasta'):
        if r.id in info_names:
            keep.add(r.id)
        else:
            drop.add(r.id)

    # drop and keep seqnames by accession
    invalid_acc = info[info['seqname'].isin(drop)]['accession']
    info = info[~info['accession'].isin(invalid_acc)]
    info = info[info['seqname'].isin(keep)]

    '''
    Second pass:
    1. Write sequences in info starting with new sequences
    2. Drop records that can exist between versions favoring newer recs
    '''
    seqnames = set(info['seqname'].tolist())
    for r in SeqIO.parse(args.new_fasta, 'fasta'):
        if r.id in seqnames:
            SeqIO.write(r, args.fasta_out, 'fasta')
            seqnames.remove(r.id)

    info.to_csv(args.annotations_out, index=False)

    '''
    just like with the seq_info, replace old pubmed_ids with newly
    downloaded pubmed_ids, by accession
    '''
    new_pubmed = pandas.read_csv(args.new_pubmed_info, dtype=str)
    prev_pubmed = pandas.read_csv(args.previous_pubmed_info, dtype=str)
    prev_pubmed = prev_pubmed[~prev_pubmed.isin(prev_pubmed['accession'])]
    pubmed_ids = prev_pubmed.append(new_pubmed).drop_duplicates()
    pubmed_ids.to_csv(args.pubmed_info_out, index=False)

    '''
    append new references. Use latest pubmed_ids and
    only if present in pubmed_ids dataframe
    '''
    new_refs = pandas.read_csv(args.new_references, dtype=str)
    prev_refs = pandas.read_csv(args.previous_references, dtype=str)
    prev_refs = prev_refs[~prev_refs['pubmed_id'].isin(new_refs['pubmed_id'])]
    refs = prev_refs.append(new_refs)
    refs = refs[refs['pubmed_id'].isin(pubmed_ids['pubmed_id'])]
    refs = refs.drop_duplicates()
    refs.to_csv(args.references_out, index=False)

    '''
    Filter the previous_versions with only existing records in ncbi. That way
    if any previous records appear again in a future esearch
    we can re-download.
    '''
    previous_versions = set(v.strip() for v in args.previous_versions)
    previous_versions = set(v for v in previous_versions if v in ncbi)

    '''
    output downloaded versions we do not want to download again. This
    includes records with no features and records that passed filtering
    up to this point
    '''
    versions = set(info['version'].tolist())
    versions |= set(n.strip() for n in args.no_features)
    versions |= previous_versions

    for v in sorted(versions):
        args.versions_out.write(v + '\n')


if __name__ == '__main__':
    main()
