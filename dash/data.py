import datetime
import gzip
import io
import os
import pandas

try:
    import boto3
except ImportError:
    pass


def read_feather(pathspec, aws_access_key_id=None, aws_secret_access_key=None,
                 get_data=True):
    """Reads gzip-compressed feather file into a pandas data_frame from s3
    bucket (if pathspec starts with 's3://') or file path. If get_data
    is False, return None for the data_frame and provide the
    modification time only.

    """

    print('reading data from {}'.format(pathspec))
    if pathspec.startswith('s3://'):
        s3_bucket, s3_key = pathspec.replace('s3://', '').split('/', 1)
        if aws_access_key_id:
            print('using provided access keys')
            s3client = boto3.client(
                's3',
                aws_access_key_id=aws_access_key_id,
                aws_secret_access_key=aws_secret_access_key)
        else:
            # assume credentials are available in default credential chain
            print('assuming existing aws credentials')
            s3client = boto3.Session(region_name='us-west-2').client('s3')
        s3obj = s3client.get_object(Bucket=s3_bucket, Key=s3_key)
        last_modified = s3obj['LastModified']
        if get_data:
            compressed = io.BytesIO(s3obj['Body'].read())
            with gzip.GzipFile(mode='rb', fileobj=compressed) as f:
                df = pandas.read_feather(f)
        else:
            df = None
    else:
        last_modified = datetime.datetime.fromtimestamp(
            os.path.getmtime(pathspec)).isoformat()
        if get_data:
            if pathspec.endswith('.gz'):
                with gzip.open(pathspec) as f:
                    df = pandas.read_feather(f)
            else:
                df = pandas.read_feather(pathspec)
        else:
            df = None

    return df, last_modified


if __name__ == '__main__':
    data, __ = read_feather('filter_details.feather.gz')
    print(data.head())
