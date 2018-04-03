#!/usr/bin/env python3
"""
output all valid tax_ids in database
"""

import argparse
import configparser
import sqlalchemy
import sys


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        'url',
        help='Database string URI or filename.')
    p.add_argument(
        '--out',
        default=sys.stdout,
        type=argparse.FileType('w'),
        help='list of all version downloaded')
    args = p.parse_args()
    conf = configparser.SafeConfigParser(allow_no_value=True)
    conf.optionxform = str  # options are case-sensitive
    conf.read(args.url)
    url = conf.get('sqlalchemy', 'url')
    con = sqlalchemy.create_engine(url).raw_connection()
    result = con.execute('select tax_id from nodes where is_valid=1')
    for i in result.fetchall():
        args.out.write(i[0] + '\n')


if __name__ == '__main__':
    main()
