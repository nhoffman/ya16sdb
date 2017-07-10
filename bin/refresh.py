#!/usr/bin/env python3
"""
Append new records with old records

Sequences that have no_features or passed all filtering into the info_out
are appended to the versions_out file.

Sequences that did not pass the vsearch or tax_id update_taxids steps will be
included in the versions_out file and consequently will be re-downloaded

Sequence id's in the fasta file are compared with the seq_info. If any seq.id
or seqname not present in either source all records by accession associated
will be dropped and re-downloaded at a later time.
"""

import argparse
import pandas

from Bio import SeqIO, Alphabet


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument(
        'new_records',
        help='records file in any Bio.SeqIO format')
    p.add_argument(
        'previous_records',
        help='records file in any Bio.SeqIO format')

    p.add_argument(
        'new_info',
        help='records file in csv format')
    p.add_argument(
        'previous_info',
        help='records file in csv format')

    p.add_argument(
        'new_pubmed_ids',
        help='pubmed_ids file')
    p.add_argument(
        'previous_pubmed_ids',
        help='pubmed_ids file')

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
        'ncbi',
        help='list of all the latest ncbi versions')
    p.add_argument(
        'vsearch',
        help='alignments for orientation information')

    p.add_argument(
        'records_out',
        type=argparse.FileType('w'),
        help='records file in same Bio.SeqIO format as input')
    p.add_argument(
        'info_out',
        help='seq_info file to output')
    p.add_argument(
        'pubmed_ids_out',
        help='pubmed_ids file to output')
    p.add_argument(
        'references_out',
        help='references file to output')
    p.add_argument(
        'versions_out',
        type=argparse.FileType('w'),
        help='list of all version downloaded')

    args = p.parse_args()

    '''
    remove old seqnames that have been replaced
    by the newly downloaded records (by accession)
    '''
    new_info = pandas.read_csv(args.new_info, dtype=str)

    vsearch = pandas.read_table(
        args.vsearch,
        header=None,
        names=['query', 'target', 'qstrand', 'id', 'tilo', 'tihi'],
        dtype=str)
    vsearch = vsearch[vsearch['target'] != '*']
    vsearch = vsearch.set_index('query')['qstrand'].to_dict()

    # ignore sequences not in 16s region
    new_info = new_info[new_info['seqname'].isin(vsearch)]

    # sort records and take the latest in case of multiple versions
    new_info['sort_date'] = pandas.to_datetime(new_info['date'])
    new_info = new_info.sort_values(by='sort_date')
    new_info = new_info.drop_duplicates(subset='seqname', keep='last')
    new_info = new_info.drop('sort_date', axis=1)

    # read in prev_info ignoring any accessions in the new data set
    prev_info = pandas.read_csv(args.previous_info, dtype=str)
    prev_info = prev_info[~prev_info['accession'].isin(new_info['accession'])]

    # combine info files and retain column order
    columns = prev_info.columns.tolist()
    for c in new_info.columns:
        if c not in columns:
            columns.append(c)
    info = prev_info.append(new_info)[columns]

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
    valid = set()
    invalid = set()

    new_recs = SeqIO.parse(
        args.new_records, 'fasta', Alphabet.IUPAC.ambiguous_dna)
    for r in new_recs:
        if r.id in info_names and Alphabet._verify_alphabet(r.seq):
            valid.add(r.id)
        else:
            invalid.add(r.id)
    for r in SeqIO.parse(args.previous_records, 'fasta'):
        if r.id in info_names:
            valid.add(r.id)
        else:
            invalid.add(r.id)

    # drop all invalid seq accessions, keep all valid seqnames
    invalid_acc = info[info['seqname'].isin(invalid)]['accession']
    info = info[~info['accession'].isin(invalid_acc)]
    info = info[info['seqname'].isin(valid)]

    '''
    * Write sequences in info starting with new sequences
    * Reverse complement new sequences in '-' orientation
    * Drop records that can exist between versions favoring newer recs
    '''
    seqnames = set(info['seqname'].tolist())
    new_recs = SeqIO.parse(
        args.new_records, 'fasta', Alphabet.IUPAC.ambiguous_dna)
    for r in new_recs:
        if r.id in seqnames:
            if vsearch.get(r.id, None) == '-':
                r.seq = r.seq.reverse_complement()
            SeqIO.write(r, args.records_out, 'fasta')
            seqnames.remove(r.id)
    for r in SeqIO.parse(args.previous_records, 'fasta'):
        if r.id in seqnames:
            SeqIO.write(r, args.records_out, 'fasta')
            seqnames.remove(r.id)

    info.to_csv(args.info_out, index=False)

    '''
    just like with the seq_info, replace old pubmed_ids with newly
    downloaded pubmed_ids, by accession
    '''
    new_pubmed = pandas.read_csv(args.new_pubmed_ids, dtype=str)
    prev_pubmed = pandas.read_csv(args.previous_pubmed_ids, dtype=str)
    prev_pubmed = prev_pubmed[~prev_pubmed.isin(prev_pubmed['accession'])]
    pubmed_ids = prev_pubmed.append(new_pubmed).drop_duplicates()
    pubmed_ids.to_csv(args.pubmed_ids_out, index=False)

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
