#!/usr/bin/env python3
# output a list of records accession.versions that have already been
# downloaded ignoring records not in a2t, have new tax_ids
# or have been modified
import argparse
import csv
import ya16sdb
import sys


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('a2t')
    p.add_argument('seq_info_cache')
    p.add_argument('modified')
    p.add_argument('cache', nargs=argparse.REMAINDER)
    p.add_argument('--out', type=argparse.FileType('w'), default=sys.stdout)
    args = p.parse_args()
    a2t = csv.DictReader(open(args.a2t))
    a2t = set((r['version'], r['tax_id']) for r in a2t)
    cseqinfo = csv.DictReader(open(args.seq_info_cache))
    cseqinfo = ((r['version'], r['tax_id']) for r in cseqinfo)
    unknown = set(r[0] for r in cseqinfo if r not in a2t)
    modified = set(ya16sdb.open_clean(args.modified))
    cache = set(ya16sdb.open_clean(args.cache))
    cache = cache - unknown - modified
    args.out.write('\n'.join(cache))


if __name__ == '__main__':
    main()
