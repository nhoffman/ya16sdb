#!/usr/bin/env python3
"""
Removes up old and modified Genbank records
"""
import argparse
import datetime
import extract_genbank as eg
import hashlib
import io
from Bio import SeqIO
import sys
import ya16sdb


def parse_genbank(genbank):
    modifiedtr = b''
    for line in genbank:
        modifiedtr += line
        if line.strip() == b'//':
            yield modifiedtr
            modifiedtr = b''


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('genbank')
    p.add_argument(
        '--records', help='master list of accessions.version to keep')
    p.add_argument('--out', default=sys.stdout, type=argparse.FileType('wb'))
    args = p.parse_args()
    if args.records:
        versions = set(v for v in ya16sdb.open_clean(args.records))
    else:
        versions = set()
    modified = {}
    keep = {}
    with open(args.genbank, 'rb') as genbank:
        for i, gb in enumerate(parse_genbank(genbank)):
            sys.stderr.write('\rPass 1: ' + str(i))
            rec = SeqIO.read(io.StringIO(gb.decode()), format='genbank')
            if not versions or eg.parse_version(rec)[0] in versions:
                md5 = hashlib.md5(gb).hexdigest()
                acc = rec.annotations['accessions'][0]
                seq_start, seq_stop = eg.parse_coordinates(rec)
                if all([seq_start, seq_stop]):
                    seq_id = '{}_{}_{}'.format(acc, seq_start, seq_stop)
                else:
                    seq_id = rec.id
                date = rec.annotations['date']
                date = datetime.datetime.strptime(date, '%d-%b-%Y')
                if seq_id not in modified or modified[seq_id] <= date:
                    modified[seq_id] = date
                    keep[seq_id] = md5
    keep = set(keep.values())
    with open(args.genbank, 'rb') as genbank:
        for i, gs in enumerate(parse_genbank(genbank)):
            sys.stderr.write('\rPass 2: ' + str(i))
            md5 = hashlib.md5(gs).hexdigest()
            if md5 in keep:
                args.out.write(gs)
                keep.remove(md5)


if __name__ == '__main__':
    main()
