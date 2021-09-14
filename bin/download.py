#!/usr/bin/env python3
import argparse
import csv
import sys
import ya16sdb


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('a2t')
    p.add_argument('cache')
    p.add_argument('do_not_download')
    p.add_argument('--out', type=argparse.FileType('w'), default=sys.stdout)
    args = p.parse_args()
    download = set(row['version'] for row in csv.DictReader(open(args.a2t)))
    cache = set(ya16sdb.open_clean(args.cache))
    do_not_download = set(ya16sdb.open_clean(args.do_not_download))
    download -= cache - do_not_download
    args.out.write('\n'.join(download))


if __name__ == '__main__':
    main()
