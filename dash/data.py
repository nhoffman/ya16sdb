"""Reads a feather file containing input data for the Dash app.

When executed as a script, downsamples the full data set to a few
select species.

"""

import boto3
import datetime
import gzip
import io
import os
import pandas as pd
import sys
import argparse
import logging


log = logging.getLogger(__name__)


def read_feather(pathspec, aws_access_key_id=None, aws_secret_access_key=None,
                 get_data=True):
    """Reads gzip-compressed feather file into a pandas data_frame from s3
    bucket (if pathspec starts with 's3://') or file path. If get_data
    is False, return None for the data_frame and provide the
    modification time only.

    """

    if pathspec.startswith('s3://'):
        s3_bucket, s3_key = pathspec.replace('s3://', '').split('/', 1)
        if aws_access_key_id:
            log.info('using provided access keys')
            s3client = boto3.client(
                's3',
                aws_access_key_id=aws_access_key_id,
                aws_secret_access_key=aws_secret_access_key)
        else:
            # assume credentials are available in default credential chain
            log.info('assuming existing aws credentials')
            s3client = boto3.Session(region_name='us-west-2').client('s3')
        s3obj = s3client.head_object(Bucket=s3_bucket, Key=s3_key)
        last_modified = s3obj['LastModified']
        if get_data:
            log.warning('reading data from {}'.format(pathspec))
            s3obj = s3client.get_object(Bucket=s3_bucket, Key=s3_key)
            compressed = io.BytesIO(s3obj['Body'].read())
            with gzip.GzipFile(mode='rb', fileobj=compressed) as f:
                df = pd.read_feather(f)
        else:
            df = None
    else:
        last_modified = datetime.datetime.fromtimestamp(
            os.path.getmtime(pathspec)).isoformat()
        if get_data:
            log.warning('reading data from {}'.format(pathspec))
            if pathspec.endswith('.gz'):
                with gzip.open(pathspec) as f:
                    df = pd.read_feather(f)
            else:
                df = pd.read_feather(pathspec)
        else:
            df = None

    return df, last_modified


def main(arguments):

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--infile', default='filter_details.feather.gz')
    parser.add_argument('-o', '--outfile', help="filtered feather file")

    args = parser.parse_args(arguments)

    data, last_modified = read_feather(args.infile)
    print(f'last modified: {last_modified}')
    print(data.describe())

    # 1351: Enterococcus faecalis
    # 1352: Enterococcus faecium
    # 851: Fusobacterium nucleatum

    filtered = (data
                .loc[data['species'].isin(['1351', '851'])]
                .loc[data['dist_pct'] < 200]
                .reset_index())
    print(filtered.describe())

    if outfile := args.outfile:
        log.info(f'writing {outfile}')
        if outfile.endswith('.gz'):
            with gzip.open(outfile, 'wb') as f:
                filtered.to_feather(f)
        else:
            filtered.to_feather(outfile)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
