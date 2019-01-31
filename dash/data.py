import os
import gzip
from datetime import datetime
from io import BytesIO

import boto3

import numpy as np
import pandas as pd
import feather


def read_feather(pathspec, aws_access_key_id=None, aws_secret_access_key=None):
    """Reads gzip-compressed feather file into a pandas data_frame from s3
    bucket (if pathspec starts with 's3://') or file path.

    """

    print('reading data from {}'.format(pathspec))
    if pathspec.startswith('s3://'):
        s3_bucket, s3_key = pathspec.replace('s3://', '').split('/', 1)
        s3client = boto3.client(
            's3',
            aws_access_key_id=aws_access_key_id,
            aws_secret_access_key=aws_secret_access_key)
        s3obj = s3client.get_object(Bucket=s3_bucket, Key=s3_key)
        last_modified = s3obj['LastModified']
        compressed = BytesIO(s3obj['Body'].read())
        with gzip.GzipFile(mode='rb', fileobj=compressed) as f:
            df = feather.read_dataframe(f)
    else:
        last_modified = datetime.fromtimestamp(
            os.path.getmtime(pathspec)).isoformat()
        with gzip.open(pathspec) as f:
            df = feather.read_dataframe(f)

    return df, last_modified


if __name__ == '__main__':
    data, __ = read_feather('filter_details.feather.gz')
    print(data.head())
