#!/usr/bin/env python3
"""
Parses Genbank records into csv format

Returns all sequences in the Plus, Forward, 5' - 3' orientation
"""
import argparse
import csv
import logging
import re
import sys

from Bio import SeqIO

log = logging.getLogger(__name__)

seq_info_columns = ['seqname', 'version', 'accession', 'name',
                    'description', 'gi', 'tax_id', 'date', 'source',
                    'keywords', 'organism', 'length', 'ambig_count',
                    'is_type', 'seq_start', 'seq_stop']

reference_columns = ['pubmed_id', 'medline_id', 'title',
                     'authors', 'journal', 'consrtm', 'comment']

pubmed_columns = ['pubmed_id', 'version', 'accession']

ACGT = frozenset('ACGT')
COORDINATES = re.compile('(?P<seq_start>\d+)..(?P<seq_stop>\d+)')


def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # ins
    parser.add_argument(
        'genbank',
        type=argparse.FileType('r'),
        default=sys.stdin,
        nargs='?',
        help='genbank record(s)')

    # outs
    parser.add_argument(
        'fasta_out',
        metavar='fasta',
        type=argparse.FileType('w'),
        help='sequence file')
    parser.add_argument(
        'info_out',
        type=argparse.FileType('w'),
        metavar='csv',
        help='write seq_info file')
    parser.add_argument(
        '--pubmed_ids',
        type=argparse.FileType('w'),
        metavar='csv',
        help=('csv with columns '
              '[version, accession, pubmed_id]'))
    parser.add_argument(
        '--references',
        type=argparse.FileType('w'),
        metavar='csv',
        help=('reference details'))

    return parser


def tax_of_genbank(gb):
    """
    Get the tax id from a genbank record,
    returning None if no taxonomy is present
    """
    # Check for bad name
    try:
        source = next(i for i in gb.features if i.type == 'source')
        taxon = next(i[6:] for i in source.qualifiers.get('db_xref', [])
                     if i.startswith('taxon:'))
        return taxon
    except StopIteration:
        return


def accession_version_of_genbank(record):
    """
    Return the accession and version of a Bio.SeqRecord.SeqRecord
    """
    annotations = record.annotations
    accession = annotations.get('accessions', [''])[0]
    if accession:
        if 'sequence_version' in annotations:
            version = int(annotations.get('sequence_version'))
        elif record.id.startswith(accession + '.'):
            version = int(record.id.split('.', 1)[1])
        else:
            version = 1
        version = '{}.{}'.format(accession, version)
    else:
        version = ''
    return accession, version


def count_ambiguous(seq):
    return sum(i not in ACGT for i in seq)


def is_type(record):
    """
    Returns a boolean indicating whether a sequence is a member of a type
    strain, as indicated by the presence of the string '(T)' within the
    record description.

    Note: NOT YET USED
    """
    type_keywords = ['(T)', 'ATCC', 'NCTC', 'NBRC', 'CCUG',
                     'DSM', 'JCM', 'NCDO', 'NCIB', 'CIP']

    for t in type_keywords:
        if t in record.description:
            return True

    return False


def parse_record(record):
    accession, version = accession_version_of_genbank(record)
    organism = record.annotations['organism']

    tax_id = tax_of_genbank(record)
    seq_start, seq_stop = parse_coordinates(record)

    if None not in (seq_start, seq_stop):
        seq_id = '{}_{}_{}'.format(accession, seq_start, seq_stop)
    else:
        seq_id = record.id

    info = dict(accession=accession,
                ambig_count=count_ambiguous(record.seq),
                length=len(record),
                seqname=seq_id,
                date=record.annotations['date'],
                description=record.description,
                gi=record.annotations.get('gi', ''),
                keywords=';'.join(record.annotations.get('keywords', [])),
                name=record.name,
                organism=organism,
                source=record.annotations['source'],
                tax_id=tax_id,
                seq_start=seq_start,
                seq_stop=seq_stop,
                version=version)
    return info


def parse_coordinates(record):
    accessions = record.annotations['accessions']
    if 'REGION:' in accessions:
        coordinates = accessions[accessions.index('REGION:') + 1]
        match = re.search(COORDINATES, coordinates)
        seq_start, seq_stop = match.group('seq_start'), match.group('seq_stop')
    else:
        seq_start, seq_stop = 1, len(record.seq)
    return seq_start, seq_stop


def parse_references(record):
    """
    Parse reference annotations that have a pubmed_id
    """
    references = []
    if 'references' in record.annotations:
        refs = [r for r in record.annotations['references'] if r.pubmed_id]
        for r in refs:
            references.append(
                dict(title=r.title,
                     authors=r.authors,
                     comment=r.comment,
                     consrtm=r.consrtm,
                     journal=r.journal,
                     medline_id=r.medline_id,
                     pubmed_id=r.pubmed_id))
    return references


def main():
    args = build_parser().parse_args()

    # setup csv output files
    info_out = csv.DictWriter(
        args.info_out,
        fieldnames=seq_info_columns,
        extrasaction='ignore')
    info_out.writeheader()

    if args.pubmed_ids:
        pubmed_ids_out = csv.DictWriter(
            args.pubmed_ids,
            fieldnames=pubmed_columns,
            extrasaction='ignore')
        pubmed_ids_out.writeheader()
    else:
        pubmed_ids_out = None

    if args.references:
        references_out = csv.DictWriter(
            args.references,
            fieldnames=reference_columns,
            extrasaction='ignore')
        references_out.writeheader()
    else:
        references_out = None

    for g in SeqIO.parse(args.genbank, 'genbank'):
        record = parse_record(g)
        args.fasta_out.write('>{}\n{}\n'.format(record['seqname'], g.seq))
        info_out.writerow(record)

        if pubmed_ids_out:
            references = parse_references(g)
            for r in references:
                pubmed_ids_out.writerow(
                    {'version': record['version'],
                     'accession': record['accession'],
                     'pubmed_id': r['pubmed_id']})
            if references_out:
                references_out.writerows(references)


if __name__ == '__main__':
    main()
