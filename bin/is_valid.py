#!/usr/bin/env python3
"""
Filter csv by tax_id validity
"""

import argparse
import csv
import sqlalchemy
import sys


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument(
        'csv',
        type=argparse.FileType('r'),
        help='file csv with column tax_id')

    p.add_argument(
        'url',
        help='Database string URI or filename.')

    db_parser = p.add_argument_group(title='database options')
    db_parser.add_argument(
        '--schema',
        help=('Name of SQL schema in database to query '
              '(if database flavor supports this).'))

    p.add_argument(
        '--out',
        default=sys.stdout,
        type=argparse.FileType('w'),
        help='list of all version downloaded')

    args = p.parse_args()

    cur = sqlalchemy.create_engine(args.url).raw_connection().cursor()

    cur.execute('select tax_id, is_valid from ' + args.schema + '.nodes')
    is_valid = dict(cur.fetchall())

    reader = csv.reader(args.csv, delimiter=',', quotechar='"')
    header = next(reader)
    tax_id = header.index('tax_id')

    rows = (r for r in reader if is_valid[r[tax_id]])

    writer = csv.writer(args.out)
    writer.writerow(header)
    writer.writerows(rows)


if __name__ == '__main__':
    main()
