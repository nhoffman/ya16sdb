#!/usr/bin/env python

"""
"""

import argparse
import csv
import itertools
import gviz_api
import jinja2
import json
import logging
import os
import pandas as pd
import sys

import bokeh
from bokeh.plotting import figure, save, output_file, ColumnDataSource
from bokeh.models import CustomJS
from bokeh.models.widgets import DataTable, TableColumn
from bokeh.layouts import gridplot

# from bokeh import colors

logging.basicConfig(
    file=sys.stdout,
    format='%(levelname)s %(module)s %(lineno)s %(message)s',
    level=logging.WARNING)

log = logging


def get_color(name, alpha=None):

    col = getattr(bokeh.colors, name).copy()

    if alpha is not None:
        assert 0.0 <= alpha <= 1.0
        col.a = alpha

    return col


def paired_plots(data, title=None, text_cols=['seqname']):
    """Create an interactive scatterplot. ``data`` is a pandas dataframe
    with (at least) column 'dist' in addition to columns
    containing other features. ``text_cols`` is a list of columns to
    include in the DataTable.

    Returns a pair of plot objects (plt, tab).

    """

    red = get_color('red', alpha=0.5)
    black = get_color('black', alpha=0.5)

    data['i'] = range(len(data))
    data['color'] = data['is_out'].apply(lambda x: red if x else black)
    data = data.fillna('')
    data.loc[data['is_out'], 'size'] = 15
    data.loc[~data['is_out'], 'size'] = 10

    source = ColumnDataSource(data)

    # start with is_outs in table
    text_source = ColumnDataSource(data=data[data['is_out']])

    callback = CustomJS(
        args=dict(source=source, text_source=text_source),
        code="""
        var inds = cb_obj.get('selected')['1d'].indices;
        var data = source.get('data');
        var text_data = text_source.get('data');

        for (var column in text_data){
            text_data[column] = [];
            for (var i = 0; i < inds.length; i++) {
                var ind = inds[i];
                text_data[column].push(data[column][ind])
            }
        }

        source.trigger('change');
        text_source.trigger('change');
        """)

    tools = ['box_select', 'lasso_select', 'resize',
             'box_zoom', 'pan', 'reset', 'tap']

    step_plt = figure(
        title=title,
        plot_width=600,
        plot_height=600,
        tools=tools,
        x_axis_label='rank order',
        y_axis_label='% distance from centroid',
    )

    # underlying markers are visible only when selected
    step_plt.scatter(
        source=text_source,
        size='size',
        x='i',
        y='dist_from_cent',
        marker='circle',
        color='color',
    )

    step_plt.scatter(
        source=source,
        size='size',
        x='i',
        y='dist_from_cent',
        marker='circle',
        color='color',
    )

    pca_plt = figure(
        title=None,
        plot_width=600,
        plot_height=600,
        tools=tools,
        x_axis_label='distance',
        y_axis_label='distance'
    )

    # underlying markers are visible only when selected
    pca_plt.circle(
        source=text_source,
        size='size',
        x='x',
        y='y',
        color='color',
    )

    pca_plt.circle(
        source=source,
        size='size',
        x='x',
        y='y',
        color='color',
    )

    formatters = {
        'seqname': bokeh.models.widgets.tables.HTMLTemplateFormatter(
            template='<a href="https://www.ncbi.nlm.nih.gov/nuccore/'
                     '<%= version %>"><%= value %></a>'),
        'type strain hit': bokeh.models.widgets.tables.HTMLTemplateFormatter(),
        'description': bokeh.models.widgets.tables.HTMLTemplateFormatter(
            template='<div title="<%= value %>"><%= value %></div>'),
        'identity': bokeh.models.widgets.tables.NumberFormatter(format='0.0%'),
        'dist': bokeh.models.widgets.tables.NumberFormatter(format='0.00%')
    }

    titles = {
        'dist': 'distance from centroid'
    }

    tab = DataTable(
        source=text_source,
        columns=[TableColumn(
            field=col,
            title=titles.get(col, col),
            formatter=formatters.get(col, None))
            for col in text_cols],
        fit_columns=True,
        row_headers=False,
        width=1200,
        height=1200,
        sortable=True)

    # note that the callback is added to the source for the scatter plot only
    source.callback = callback

    return step_plt, pca_plt, tab


def main(arguments):

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        'details',
        help="details of filtering results")
    parser.add_argument(
        'taxonomy',
        help="taxonomy in csv format")
    parser.add_argument(
        'hits',
        help="output of annotate_hits.py")
    parser.add_argument(
        'hits_seq_info',
        help='seq_info for reference seqs')
    parser.add_argument(
        '--param',
        action='append',
        default=[],
        help='colon separated index.html information')
    parser.add_argument(
        '--index-template',
        default='data/index.html.jinja2',
        help='location for index.html template [%(default)s]')
    parser.add_argument(
        '--log',
        help="output of 'deenurp filter_is_outs --log'",
        type=argparse.FileType('r'))

    parser.add_argument(
        '-i',
        '--plot-index',
        default='index.html',
        help=('name of plot index file, assumed to '
              'have the same parent dir as plot_dir'))
    parser.add_argument(
        '-d',
        '--plot-dir',
        default='plot_details',
        help=('output directory for plots; index file '
              'will be named $plot_dir.html'))
    parser.add_argument(
        '--plot-map',
        metavar='CSV',
        type=argparse.FileType('w'),
        help='tax_id,outfile')
    parser.add_argument(
        '-N',
        type=int,
        help='max number of taxa to read')

    args = parser.parse_args(arguments)

    plot_dir = args.plot_dir
    plot_dir_name = os.path.basename(plot_dir)
    plot_index = args.plot_index
    json_data_file = 'index.json'

    try:
        os.makedirs(plot_dir)
    except OSError:
        pass

    if args.log:
        deenurp_args = {}
        for line in args.log:
            spl = line.split('args: ')
            if spl[-1].startswith('{'):
                deenurp_args = json.loads(spl[-1])
                break

    taxonomy = pd.read_csv(
        args.taxonomy,
        usecols=['tax_id', 'species', 'tax_name'],
        dtype={'species': str, 'tax_id': str})
    taxonomy = taxonomy.set_index('tax_id')

    species = taxonomy[['species']].join(
        taxonomy['tax_name'],
        on='species',
        how='inner')
    species = species.drop_duplicates().set_index('species')

    hits = pd.read_table(
        args.hits,
        header=None,
        dtype={'seqname': str, 'hit': str, 'pct_id': float},
        na_filter=False,
        names=['seqname', 'hit', 'pct_id'],
        usecols=['seqname', 'hit', 'pct_id'])
    hits = hits.set_index('seqname')

    hit_info = pd.read_csv(
        args.hits_seq_info,
        dtype=str,
        usecols=['seqname', 'tax_id', 'version'])
    hit_info = hit_info.set_index('seqname')

    cols = ['seqname', 'tax_id', 'ambig_count', 'centroid',
            'version', 'cluster', 'dist', 'is_out', 'species',
            'x', 'y', 'is_type', 'description']
    details = pd.read_csv(
        args.details,
        usecols=cols,
        dtype={'species': str, 'tax_id': str, 'version': str, 'dist': float})
    details = details.set_index('seqname')

    hits = hits.join(hit_info, on='hit')  # get tax_id
    hits = hits.join(taxonomy['species'], on='tax_id')  # get species_id
    hits = hits.join(species, on='species')  # get species tax name
    hits = hits.rename(
        columns={'species': 'hit_id',
                 'tax_name': 'type strain hit',
                 'version': 'hit_version'})
    hits = hits.drop('tax_id', axis=1)  # we only care about the species_id

    details = details.join(species, on='species')  # get species tax name
    details = details.join(hits)
    details = details.reset_index()  # move seqname back as a column

    # output file purposes
    details.loc[:, 'species_name'] = details['tax_name']

    # make the hit_name link back to an ncbi url to the actual record
    url = '<a href="https://www.ncbi.nlm.nih.gov/nuccore/{}">{}</a>'

    def hit_ncbi_link(row):
        if pd.isnull(row['type strain hit']):
            tag = ''
        else:
            tag = url.format(row['hit_version'], row['type strain hit'])
        return tag
    details.loc[:, 'type strain hit'] = details.apply(hit_ncbi_link, axis=1)

    details['outlier'] = details['is_out'].apply(
        lambda x: 'yes' if x else '')
    details['type strain'] = details['is_type'].apply(
        lambda x: 'yes' if x else '')
    # to retain sorting and use with NumberFormatter
    details['identity'] = details['pct_id'] / 100
    details['dist_from_cent'] = details['dist'] * 100

    # define table layout for index page
    table = gviz_api.DataTable(
        [('tax_id', 'string'),
         ('species', 'string'),
         ('records', 'number'),
         ('clustered', 'boolean'),
         ('is_outs', 'number'),
         ('pct out', 'number')]
    )

    if args.plot_map:
        map_out = csv.writer(args.plot_map)
        map_out.writerow(['tax_id', 'html'])

    table_data = []
    by_species = details.groupby(by=['species', 'species_name'], sort=False)
    by_species = itertools.islice(by_species, args.N)
    for (species_id, species_name), species in by_species:
        if species.x.isnull().all():
            continue

        species = species.copy().sort_values(by='dist')
        label = '_'.join(species_name.lower().split())

        step_plt, pca_plt, tab = paired_plots(
            species,
            title=species_name,
            text_cols=[
                'seqname',
                'description',
                'type strain',
                'dist',
                'outlier',
                'type strain hit',
                'identity'
            ])

        filename = '{}.html'.format(label)
        output_file(
            filename=os.path.join(args.plot_dir, filename),
            title=species_name)
        save(gridplot([[step_plt, pca_plt], [tab]]))

        if args.plot_map:
            map_out.writerow([species_id, filename])

        n_out = sum(species['is_out'])
        pct_out = (100.0 * n_out) / len(species)
        clustered = not (species.cluster == -1).all()

        table_data.append([
            species_id,
            '<a href="{}/{}">{}</a>'.format(
                plot_dir_name, filename, species_name),
            len(species),
            clustered,
            n_out if clustered else None,
            (pct_out, '%.01f%%' % pct_out) if clustered else None,
        ])

    with open(os.path.join(plot_dir, json_data_file), 'w') as fout:
        table.AppendData(table_data)
        fout.write(table.ToJSon())

    # plot index should refer to json data according to its relative
    # path, assuming that this file is in same directory as plot_dir
    with open(plot_index, 'w') as fout:
        index = jinja2.Template(open(args.index_template).read())
        fout.write(index.render(
            json_data=os.path.join(plot_dir_name, json_data_file),
            deenurp_args=deenurp_args,
            params=dict(p.split(':') for p in args.param)
        ))


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
