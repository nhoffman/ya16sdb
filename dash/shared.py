'''
Shared global data for the Dash app
'''
import data
import logging
import os

log = logging.getLogger(__name__)

COLORS = ['blue', 'red', 'black', 'yellow',
          'gray', 'green', 'violet', 'silver']
DEFAULT_COLOR = 'is_out'
DEFAULT_SHAPE = 'confidence'
DEFAULT_GENUS = '1350'  # Enterococcus
DEFAULT_GENUS_NAME = 'Enterococcus'
DEFAULT_SPECIES = '1351'  # Enterococcus faecalis
DEFAULT_Y = 'y'
DEFAULT_X = 'x'
LEGEND_OTHER = 'other'
MAX_TABLE_RECORDS = 500
NO_DATA = 'No Data'
SEARCH_OPTS = ['seqname', 'accession', 'version',
               'species_name', 'species', 'genus']
SHAPES = ['circle', 'triangle-up', 'square', 'diamond',
          'pentagon', 'cross', 'star', 'hourglass']

FEATHER_FILE = os.environ.get('DATA_FILE', 'filter_details.feather.gz')
AWS_ACCESS_KEY_ID = os.environ.get('BUCKET_ACCESS_KEY')
AWS_SECRET_KEY_ID = os.environ.get('BUCKET_SECRET_KEY')

# GLOBAL VARS
df = None
genera = None
last_modified = None
seq_info = None
tax = None
species_to_genus = None
species_to_id = None


def set_global_data():
    global df, genera, last_modified, seq_info, species_to_genus, species_to_id, tax
    _, modified = data.read_feather(
        FEATHER_FILE,
        aws_access_key_id=AWS_ACCESS_KEY_ID,
        aws_secret_access_key=AWS_SECRET_KEY_ID,
        get_data=False)
    if last_modified != modified:
        df, last_modified = data.read_feather(
            FEATHER_FILE,
            aws_access_key_id=AWS_ACCESS_KEY_ID,
            aws_secret_access_key=AWS_SECRET_KEY_ID,
            get_data=True)
        df = df.astype({
            'modified_date': 'datetime64[ns]',
            'download_date': 'datetime64[ns]'
             })
        seq_info = df.copy()
        df = df[~df['x'].isna() & ~df['y'].isna()]
        df['genus_name'] = df['genus_name'].fillna('Unclassified')
        df['genus'] = df['genus'].fillna('')
        tcs = 'taxonomy-check-status'
        df[tcs] = df[tcs].fillna(NO_DATA)
        bmsn = 'best-match-species-name'
        df[bmsn] = df[bmsn].fillna(NO_DATA)
        tax = df[['genus', 'genus_name', 'species', 'species_name']]
        tax = tax.drop_duplicates().sort_values(
            by=['genus_name', 'species_name'])
        species_to_genus = dict(
            tax[['species', 'genus']].drop_duplicates().values)
        species_to_id = dict(
            tax[['species_name', 'species']].drop_duplicates().values)
        genera = tax.groupby(by='genus')


set_global_data()
