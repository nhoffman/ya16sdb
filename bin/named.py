#!/usr/bin/env python3
"""
output all valid tax_ids in database
"""

import argparse
import sqlalchemy
import sys


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        'url',
        help='Database string URI or filename.')
    db_parser = p.add_argument_group(title='database options')
    db_parser.add_argument(
        '--schema',
        type=lambda x: x + '.',
        default='',
        help=('Name of SQL schema in database to query '
              '(if database flavor supports this).'))
    p.add_argument(
        '--out',
        default=sys.stdout,
        type=argparse.FileType('w'),
        help='list of all version downloaded')

    args = p.parse_args()
    cur = sqlalchemy.create_engine(args.url).raw_connection().cursor()
    cur.execute('select tax_id from ' + args.schema + 'nodes where is_valid')
    for i in cur.fetchall():
        args.out.write(i[0] + '\n')


if __name__ == '__main__':
    main()
