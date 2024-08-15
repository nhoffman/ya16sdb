#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--out', type=argparse.FileType('w'), default=sys.stdout)
    return parser.parse_args()


def main():
    args = get_args()
    if 'YA16SDB_VERSION' in os.environ:
        ver = os.environ['YA16SDB_VERSION']
    else:
        try:
            cmd = ['git', 'describe', '--tags', '--dirty']
            ver = subprocess.check_output(cmd, text=True).strip()
        except subprocess.CalledProcessError:
            ver = ''
    args.out.write(ver)


if __name__ == '__main__':
    sys.exit(main())
