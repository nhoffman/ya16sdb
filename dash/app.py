#! /usr/bin/env python3
'''
Plotly Dash app exploring NCBI 16s records grouped by species taxonomy id
'''

import dash
import dash_core_components as dcc
import dash_html_components as html
import feather
import gzip
import os
import pandas
import urllib

from dash.dependencies import Input, State, Output

from data import read_feather

COLORS = ['blue', 'red', 'black', 'yellow',
          'gray', 'green', 'violet', 'silver']
DEFAULT_COLOR = 'is_out'
DEFAULT_SHAPE = 'confidence'
DEFAULT_GENUS = '1350'  # Enterococcus
DEFAULT_Y = 'y'
DEFAULT_X = 'x'
LEGEND_OTHER = 'other'
MAX_TABLE_RECORDS = 500
SEARCH_OPTS = ['seqname', 'accession', 'version',
               'species_name', 'species', 'genus']
SHAPES = ['circle', 'triangle-up', 'square', 'diamond',
          'pentagon', 'cross', 'star', 'hourglass']

app = dash.Dash()
app.title = 'Species Outlier Plots'

# access all environment variables here
FEATHER_FILE = os.environ.get('DATA_FILE', 'filter_details.feather.gz')
AWS_ACCESS_KEY_ID = os.environ.get('BUCKET_ACCESS_KEY')
AWS_SECRET_KEY_ID = os.environ.get('BUCKET_SECRET_KEY')

df, df_last_modified = read_feather(
    FEATHER_FILE,
    aws_access_key_id=AWS_ACCESS_KEY_ID,
    aws_secret_access_key=AWS_SECRET_KEY_ID)


df = df[~df['x'].isna() & ~df['y'].isna()]
df['genus_name'] = df['genus_name'].fillna('Unclassified')
df['genus'] = df['genus'].fillna('')
axes = ['confidence', 'dist_pct', 'x', 'y', 'match_pct',
        'match_species', 'match_version', 'rank_order']
tax = df[['genus', 'genus_name', 'species', 'species_name']]
tax = tax.drop_duplicates().sort_values(by=['genus_name', 'species_name'])
species_genus = dict(tax[['species', 'genus']].drop_duplicates().values)
species_id = dict(tax[['species_name', 'species']].drop_duplicates().values)
genera = tax.groupby(by='genus')
genus_opts = tax[['genus', 'genus_name']].drop_duplicates().values
genus_opts = [{'label': gn, 'value': gi} for gi, gn in genus_opts]

app.layout = html.Div(
    style={'width': 1175},
    children=[
        html.Div(id='state'),
        dcc.Location(id='url', refresh=False),
        dcc.Markdown(
            children=[
                str(df['download_date'].max().strftime('%A, %B %d, %Y'))
            ],
            containerProps={
                'style': {
                    'text-align': 'right',
                    'font-style': 'italic',
                    'height': 0,
                    'width': '100%'}}),
        html.Div(
            children=[
                dcc.Input(type='text', id='text-input'),
                html.Button(
                    id='submit-button',
                    n_clicks=0,
                    children='Search')]),
        dcc.Markdown(
            children=['**Genus**'],
            containerProps={
                'style': {
                    'display': 'inline-block',
                    'text-align': 'center',
                    'vertical-align': 'middle',
                    'width': '5%'}}),
        html.Div(
            children=[
                dcc.Dropdown(
                    id='genus-column',
                    options=genus_opts,
                    clearable=False)],
            style={
                'display': 'inline-block',
                'vertical-align': 'middle',
                'width': '40%'}),
        dcc.Markdown(
            children=['**Species**'],
            containerProps={
                'style': {
                    'display': 'inline-block',
                    'text-align': 'center',
                    'vertical-align': 'middle',
                    'width': '6%'}}),
        html.Div(
            children=[dcc.Dropdown(id='species-column', clearable=False)],
            style={
                'display': 'inline-block',
                'vertical-align': 'middle',
                'width': '48%'}),
        html.Div(
            children=[
                html.Div(style={'width': '11%', 'display': 'inline-block'}),
                dcc.Markdown(
                    children=['**Color**'],
                    containerProps={
                        'style': {
                            'width': '5%',
                            'display': 'inline-block'}}),
                dcc.Markdown(
                    children=['**Shape**'],
                    containerProps={
                        'style': {
                            'width': '10%',
                            'display': 'inline-block'}}),
                dcc.Markdown(
                     children=['**Selection**'],
                     containerProps={
                         'style': {
                             'width': '39%',
                             'display': 'inline-block'}}),
                dcc.Markdown(
                    children=['**Visibility**'],
                    containerProps={
                        'style': {
                            'width': '35%',
                            'display': 'inline-block'}}),
                html.Div(
                    children=[
                        dcc.Markdown(children=['**Outliers**']),
                        dcc.Markdown(children=['**Confidence**']),
                        dcc.Markdown(children=['**Match Species**']),
                        dcc.Markdown(children=['**Isolation Source**']),
                        ],
                    style={
                        'vertical-align': 'middle',
                        'width': '11%',
                        'display': 'inline-block',
                        }),
                dcc.RadioItems(
                    id='color-items',
                    options=[
                        {'value': 'is_out'},
                        {'value': 'confidence'},
                        {'value': 'match_species'},
                        {'value': 'isolation_source'}],
                    inputStyle={'height': 15, 'width': 15, 'margin': 11},
                    style={
                        'vertical-align': 'middle',
                        'width': '5%',
                        'display': 'inline-block'}),
                dcc.RadioItems(
                    id='shape-items',
                    options=[
                        {'value': 'is_out'},
                        {'value': 'confidence'},
                        {'value': 'match_species'},
                        {'value': 'isolation_source'}],
                    inputStyle={'height': 15, 'width': 15, 'margin': 11},
                    style={
                        'vertical-align': 'middle',
                        'width': '5%',
                        'display': 'inline-block'}),
                html.Div(
                    children=[
                        dcc.Dropdown(
                            id='outliers-selection',
                            multi=True),
                        dcc.Dropdown(
                            id='confidence-selection',
                            multi=True),
                        dcc.Dropdown(
                            id='match-species-selection',
                            multi=True),
                        dcc.Dropdown(
                            id='isolation-source-selection',
                            multi=True)],
                    style={
                        'vertical-align': 'middle',
                        'width': '39%',
                        'display': 'inline-block'}),
                html.Div(
                    children=[
                        dcc.Dropdown(
                            id='outliers-visibility',
                            multi=True),
                        dcc.Dropdown(
                            id='confidence-visibility',
                            multi=True),
                        dcc.Dropdown(
                            id='match-species-visibility',
                            multi=True),
                        dcc.Dropdown(
                            id='isolation-source-visibility',
                            multi=True)],
                    style={
                        'vertical-align': 'middle',
                        'width': '39%',
                        'display': 'inline-block'})],
            style={
                'border': 'thin lightgrey solid',
                'borderRadius': 5,
                'margin': 5,
                'padding': 10,
                'width': '97%'}),
        html.Div(
            children=[dcc.Slider(id='year--slider')],
            style={'margin': 15, 'width': '95%'}),
        dcc.Markdown(
            children=['**Axes**'],
            containerProps={
                'style': {
                    'display': 'inline-block',
                    'text-align': 'center',
                    'vertical-align': 'middle',
                    'width': '4%'}}),
        html.Div(
            children=[
                dcc.Dropdown(
                    clearable=False,
                    id='yaxis-column',
                    options=[{'label': i, 'value': i} for i in axes])],
            style={
                'width': '14%',
                'display': 'inline-block',
                'vertical-align': 'middle'}),
        html.Div(
            children=[
                dcc.Dropdown(
                    clearable=False,
                    id='xaxis-column',
                    options=[{'label': i, 'value': i} for i in axes])],
            style={
                'width': '14%',
                'display': 'inline-block',
                'vertical-align': 'middle'}),
        dcc.Graph(id='plot'),
        html.Table(
            id='table-div',
            style={
                'border': 'thin lightgrey solid',
                'borderRadius': 5,
                'display': 'inline-block',
                'height': 413,
                'margin': 5,
                'overflow-y': 'scroll',
                'padding': 10,
                'width': '97%'})])


def assign_hover_text(s):
    '''
    assign hover text to outliers for now
    '''
    text = None
    if s['is_out']:
        text = ('seqname: {seqname}<br>'
                'accession: {version}<br>'
                'modified_date: {modified_date}<br>'
                'match_species: {match_species}<br>'
                'isolation source: {isolation_source}'.format(**s))
    return text


def parse_search_input(dff, state, search, n_clicks, text):
    '''
    Determine where request came from url or text input
    '''
    request = None
    data = None
    if state is None:  # new instance of the web page
        args = urllib.parse.parse_qs(urllib.parse.urlparse(search).query)
        for o in SEARCH_OPTS:
            if o in args:
                request = o
                data = args[o][0]
    elif text is not None and state['n_clicks'] < n_clicks:  # button clicked
        data = text.strip()
        for o in SEARCH_OPTS:
            if text in dff[o].values:
                request = o
    return request, data


@app.callback(
    Output('yaxis-column', 'value'),
    [Input('url', 'search')])
def update_yaxis_value(search):
    args = urllib.parse.parse_qs(urllib.parse.urlparse(search).query)
    return args.get('y', [DEFAULT_Y])[0]


@app.callback(
    Output('xaxis-column', 'value'),
    [Input('url', 'search')])
def update_xaxis_value(search):
    args = urllib.parse.parse_qs(urllib.parse.urlparse(search).query)
    return args.get('x', [DEFAULT_X])[0]


# Setting the url search does not seem to work at the moment in addition
# to the circular logic this introduces:
# https://github.com/plotly/dash-core-components/issues/44
# @app.callback(
#     Output('url', 'search'),
#     [Input('species-column', 'value')])
# def update_url_search(value):
#     return '?' + urllib.parse.urlencode({'species_id': value})


@app.callback(
    Output('genus-column', 'value'),
    [Input('url', 'search'),
     Input('submit-button', 'n_clicks')],
    [State('text-input', 'value'),
     State('state', 'hidden')])
def update_genus_value(search, n_clicks, text, state):
    request, data = parse_search_input(df, state, search, n_clicks, text)
    if request is None:
        value = DEFAULT_GENUS
    elif request == 'species_name':
        value = species_genus.get(data, DEFAULT_GENUS)
    else:
        value = df[df[request] == data].iloc[0]['genus']
    return value


@app.callback(
    Output('species-column', 'options'),
    [Input('genus-column', 'value')])
def update_species_options(tax_id):
    group = genera.get_group(tax_id)
    group = group[['species', 'species_name']]
    return [{'label': sn, 'value': si} for si, sn in group.values]


@app.callback(
    Output('species-column', 'value'),
    [Input('species-column', 'options'),
     Input('url', 'search'),
     Input('submit-button', 'n_clicks')],
    [State('text-input', 'value'),
     State('state', 'hidden'),
     State('genus-column', 'value')])
def update_species_value(options, search, n_clicks, text, state, tax_id):
    dff = df[df['genus'] == tax_id]
    request, data = parse_search_input(dff, state, search, n_clicks, text)
    if request is None:  # no request
        value = options[0]['value']  # return first item in dropdown
    elif request == 'species_name':
        value = species_id.get(data, options[0]['value'])
    else:
        value = dff[dff[request] == data].iloc[0]['species']
    return value


@app.callback(
    Output('year--slider', 'value'),
    [Input('species-column', 'value'),
     Input('isolation-source-visibility', 'value')])
def update_slider_value(tax_id, iso_values):
    '''
    Reset the slider to generate the axis ranges with all data points.
    This will also help avoid confusion when users select a new species
    to view.
    '''
    dff = df[df['species'] == tax_id]
    if iso_values:
        dff = dff[dff['isolation_source'].isin(iso_values)]
    return dff['modified_date'].max().year


@app.callback(
    Output('year--slider', 'min'),
    [Input('species-column', 'value'),
     Input('isolation-source-visibility', 'value')])
def update_slider_min(tax_id, iso_values):
    '''
    reset min year to avoid None errors when drawing figure
    '''
    dff = df[df['species'] == tax_id]
    if iso_values:
        dff = dff[dff['isolation_source'].isin(iso_values)]
    return dff['modified_date'].min().year


@app.callback(
    Output('year--slider', 'max'),
    [Input('species-column', 'value'),
     Input('isolation-source-visibility', 'value')])
def update_slider_max(tax_id, iso_values):
    '''
    reset max year
    '''
    dff = df[df['species'] == tax_id]
    if iso_values:
        dff = dff[dff['isolation_source'].isin(iso_values)]
    return dff['modified_date'].max().year


@app.callback(
    Output('year--slider', 'marks'),
    [Input('species-column', 'value'),
     Input('isolation-source-visibility', 'value')])
def update_slider_marks(tax_id, iso_values):
    '''
    reset marks
    '''
    dff = df[df['species'] == tax_id]
    if iso_values:
        dff = dff[dff['isolation_source'].isin(iso_values)]
    modified_dates = dff['modified_date'].apply(lambda x: x.year).unique()
    return {str(year): str(year) for year in modified_dates}


@app.callback(
    Output('isolation-source-selection', 'options'),
    [Input('species-column', 'value')])
def update_isolation_selection(tax_id):
    '''
    '''
    dff = df[df['species'] == tax_id]
    iso = dff['isolation_source']
    options = []
    for k, v in iso.value_counts().items():
        options.append({'label': '{} ({})'.format(k, v), 'value': k})
    return options


@app.callback(
    Output('isolation-source-visibility', 'options'),
    [Input('species-column', 'value')])
def update_isolation_visiblity(tax_id):
    '''
    '''
    dff = df[df['species'] == tax_id]
    iso = dff['isolation_source']
    options = []
    for k, v in iso.value_counts().items():
        options.append({'label': '{} ({})'.format(k, v), 'value': k})
    return options


@app.callback(
    Output('match-species-selection', 'options'),
    [Input('species-column', 'value')])
def update_match_species_selection(tax_id):
    '''
    '''
    dff = df[df['species'] == tax_id]
    match_species = dff['match_species']
    options = []
    for k, v in match_species.value_counts().items():
        options.append({'label': '{} ({})'.format(k, v), 'value': k})
    return options


@app.callback(
    Output('match-species-visibility', 'options'),
    [Input('species-column', 'value')])
def update_match_species_visibility(tax_id):
    '''
    '''
    dff = df[df['species'] == tax_id]
    match_species = dff['match_species']
    options = []
    for k, v in match_species.value_counts().items():
        options.append({'label': '{} ({})'.format(k, v), 'value': k})
    return options


@app.callback(
    Output('outliers-selection', 'options'),
    [Input('species-column', 'value')])
def update_outliers_selection(tax_id):
    '''
    '''
    dff = df[df['species'] == tax_id]
    out = dff['is_out'].apply(lambda x: 'Yes' if x else 'No')
    options = []
    for k, v in out.value_counts().items():
        options.append({'label': '{} ({})'.format(k, v), 'value': k})
    return options


@app.callback(
    Output('outliers-visibility', 'options'),
    [Input('species-column', 'value')])
def update_outliers_visibility(tax_id):
    '''
    '''
    dff = df[df['species'] == tax_id]
    out = dff['is_out'].apply(lambda x: 'Yes' if x else 'No')
    options = []
    for k, v in out.value_counts().items():
        options.append({'label': '{} ({})'.format(k, v), 'value': k})
    return options


@app.callback(
    Output('confidence-selection', 'options'),
    [Input('species-column', 'value')])
def update_confidence_selection(tax_id):
    '''
    '''
    dff = df[df['species'] == tax_id]
    options = []
    for k, v in dff['confidence'].value_counts().items():
        options.append({'label': '{} ({})'.format(k, v), 'value': k})
    return options


@app.callback(
    Output('confidence-visibility', 'options'),
    [Input('species-column', 'value')])
def update_confidence_visibility(tax_id):
    '''
    '''
    dff = df[df['species'] == tax_id]
    options = []
    for k, v in dff['confidence'].value_counts().items():
        options.append({'label': '{} ({})'.format(k, v), 'value': k})
    return options


def parse_multi(state, search, option):
    value = None
    if state is None:
        args = urllib.parse.parse_qs(urllib.parse.urlparse(search).query)
        if option in args:
            value = args[option][0].split(',')
    return value


@app.callback(
    Output('isolation-source-selection', 'value'),
    [Input('species-column', 'value')],
    [State('state', 'hidden'),
     State('url', 'search')])
def update_isolation_source_selection_value(_, state, search):
    return parse_multi(state, search, 'selection_isolation_source')


@app.callback(
    Output('isolation-source-visibility', 'value'),
    [Input('species-column', 'value')],
    [State('state', 'hidden'),
     State('url', 'search')])
def update_isolation_source_visiblity_value(_, state, search):
    return parse_multi(state, search, 'visibility_isolation_source')


@app.callback(
    Output('match-species-selection', 'value'),
    [Input('species-column', 'value')],
    [State('state', 'hidden'),
     State('url', 'search')])
def update_match_species_selection_value(_, state, search):
    return parse_multi(state, search, 'selection_match_species')


@app.callback(
    Output('match-species-visibility', 'value'),
    [Input('species-column', 'value')],
    [State('state', 'hidden'),
     State('url', 'search')])
def update_match_species_visibility_value(_, state, search):
    return parse_multi(state, search, 'visibility_match_species')


@app.callback(
    Output('outliers-selection', 'value'),
    [Input('species-column', 'value')],
    [State('state', 'hidden'),
     State('url', 'search')])
def update_outliers_selection_value(_, state, search):
    return parse_multi(state, search, 'selection_is_out')


@app.callback(
    Output('outliers-visibility', 'value'),
    [Input('species-column', 'value')],
    [State('state', 'hidden'),
     State('url', 'search')])
def update_outliers_visibility_value(_, state, search):
    return parse_multi(state, search, 'visibility_is_out')


@app.callback(
    Output('confidence-selection', 'value'),
    [Input('species-column', 'value')],
    [State('state', 'hidden'),
     State('url', 'search')])
def update_confidence_selection_value(_, state, search):
    return parse_multi(state, search, 'selection_confidence')


@app.callback(
    Output('confidence-visibility', 'value'),
    [Input('species-column', 'value')],
    [State('state', 'hidden'),
     State('url', 'search')])
def update_confidence_visibility_value(_, state, search):
    return parse_multi(state, search, 'visibility_confidence')


@app.callback(
    Output('color-items', 'value'),
    [Input('url', 'search')])
def update_color_items(search):
    args = urllib.parse.parse_qs(urllib.parse.urlparse(search).query)
    return args.get('color', [DEFAULT_COLOR])[0]


@app.callback(
    Output('shape-items', 'value'),
    [Input('url', 'search')])
def update_symbol_items(search):
    args = urllib.parse.parse_qs(urllib.parse.urlparse(search).query)
    return args.get('shape', [DEFAULT_SHAPE])[0]


@app.callback(
    Output('plot', 'figure'),
    [Input('species-column', 'value'),
     Input('xaxis-column', 'value'),
     Input('yaxis-column', 'value'),
     Input('year--slider', 'value'),
     Input('isolation-source-visibility', 'value'),
     Input('isolation-source-selection', 'value'),
     Input('match-species-visibility', 'value'),
     Input('match-species-selection', 'value'),
     Input('outliers-visibility', 'value'),
     Input('outliers-selection', 'value'),
     Input('confidence-visibility', 'value'),
     Input('confidence-selection', 'value'),
     Input('color-items', 'value'),
     Input('shape-items', 'value'),
     Input('submit-button', 'n_clicks')],
    [State('state', 'hidden'),
     State('text-input', 'value'),
     State('url', 'search')])
def update_graph(tax_id, xaxis, yaxis, year_value,
                 viso_source, siso_source,
                 vmatch, smatch,
                 vout, sout,
                 vconf, sconf,
                 color, symbol,
                 n_clicks, state, text, search):
    dff = df[df['species'] == tax_id]
    dff = dff[dff['modified_date'] <= str(year_value+1)]

    # decide if we should allow plot to calculate axes ranges
    if state is None:
        x_range = None
        y_range = None
    elif tax_id != state['tax_id']:
        x_range = None
        y_range = None
    elif state['xaxis'] != xaxis:
        x_range = None
        y_range = None
    elif state['yaxis'] != yaxis:
        x_range = None
        y_range = None
    else:
        x_range = state['xrange']
        y_range = state['yrange']

    # decide visibile points
    if viso_source:
        dff = dff[dff['isolation_source'].isin(viso_source)]
    if vmatch:
        dff = dff[dff['match_species'].isin(vmatch)]
    if vout:
        dff = dff[dff['is_out'].isin(set(i == 'Yes' for i in vout))]
    if vconf:
        dff = dff[dff['confidence'].isin(vconf)]

    # decide selected points
    dff['selected'] = False
    request, data = parse_search_input(dff, state, search, n_clicks, text)
    if request in ['seqname', 'accession', 'version']:
        dff.loc[dff[request] == data, 'selected'] = True
    if siso_source:
        dff.loc[dff['isolation_source'].isin(siso_source), 'selected'] = True
    if smatch:
        dff.loc[dff['match_species'].isin(smatch), 'selected'] = True
    if sout:
        is_out = dff['is_out'].isin(set(i == 'Yes' for i in sout))
        dff.loc[is_out, 'selected'] = True
    if sconf:
        dff.loc[dff['confidence'].isin(sconf), 'selected'] = True

    dff['text'] = dff.apply(assign_hover_text, axis='columns')

    # assign symbols and colors
    for col, label, styles in [[symbol, 'symbol', SHAPES],
                               [color, 'color', COLORS]]:
        name = label + '_name'
        dff[name] = LEGEND_OTHER
        dff[label] = styles[0]
        if col == 'is_out':
            dff.loc[dff[col], name] = 'outlier'
            dff.loc[dff[col], label] = styles[1]
        else:
            style_count = len(styles) - 1  # minus styles[0]
            top = dff[col].value_counts().iloc[:style_count].keys()
            for i, c in enumerate(top, start=1):
                dff.loc[dff[col] == c, name] = c
                dff.loc[dff[col] == c, label] = styles[i]

    # create scatter plots for every symbol and color combination
    data = []
    by = ['symbol_name', 'color_name', 'symbol', 'color']
    for (sym_name, clr_name, _, _), d in dff.groupby(by=by):
        d = d.copy()
        d['iselected'] = range(len(d))
        if sym_name == clr_name:
            name = sym_name
        elif sym_name == LEGEND_OTHER:
            name = clr_name
        elif clr_name == LEGEND_OTHER:
            name = sym_name
        else:
            name = '{} and {}'.format(clr_name, sym_name)
        data.append({
            'customdata': d.index,
            'hoverinfo': 'text',
            'marker': {'symbol': d['symbol'], 'color': d['color'], 'size': 12},
            'mode': 'markers',
            'name': name,
            'legendgroup': clr_name,
            'selected': {'marker': {'size': 15, 'opacity': 0.7}},
            'selectedpoints': d[d['selected']]['iselected'],
            'type': 'scatter',
            'unselected': {'marker': {'size': 10, 'opacity': 0.4}},
            'hovertext': d['text'],
            'x': d[xaxis],
            'y': d[yaxis],
            })

    outliers = dff[dff['is_out']]  # for title denominator
    figure = {
        'data': data,
        'layout': {
            'dragmode': 'select',  # default tool from modebar
            'legend': {
                'orientation': 'v',
                },
            'height': 900,
            'hovermode': 'closest',
            'margin': {
                'b': 200 if xaxis == 'match_species' else None,
                'l': 200 if yaxis == 'match_species' else None
                },
            'showlegend': True,
            'title': '{}, tax_id {} outliers: {} ({:.0%}), total: {}'.format(
                dff.iloc[0]['species_name'],
                tax_id,
                len(outliers),
                (len(outliers) / len(dff)),
                len(dff)),
            'xaxis': {
                'scaleanchor': 'y' if xaxis == 'x' else None,
                'title': xaxis,
                'showgrid': True if xaxis == 'match_species' else False,
                'range': x_range
            },
            'yaxis': {
                'scaleanchor': 'x' if yaxis == 'y' else None,
                'title': yaxis,
                'showgrid': True if yaxis == 'match_species' else False,
                'range': y_range
            }
        }}
    return figure


@app.callback(
    Output('state', 'hidden'),
    [Input('submit-button', 'n_clicks'),
     Input('species-column', 'value'),
     Input('plot', 'figure'),
     Input('xaxis-column', 'value'),
     Input('yaxis-column', 'value')])
def update_state(n_clicks, tax_id, figure, xaxis, yaxis):
    '''
    preserve state on the client to help determine how actions are processed
    here on the server
    '''
    return {
        'n_clicks': n_clicks,  # determine if button was clicked
        'tax_id': tax_id,  # to check if tax_id has changed
        'xrange': figure['layout']['xaxis']['range'],  # preserve axes ranges
        'yrange': figure['layout']['yaxis']['range'],
        'xaxis': xaxis,
        'yaxis': yaxis}


@app.callback(
    Output('table-div', 'children'),
    [Input('plot', 'selectedData'),
     Input('isolation-source-selection', 'value'),
     Input('match-species-selection', 'value'),
     Input('outliers-selection', 'value'),
     Input('confidence-selection', 'value')],
    [State('submit-button', 'n_clicks'),
     State('species-column', 'value'),
     State('text-input', 'value'),
     State('url', 'search'),
     State('state', 'hidden')])
def update_table(selected, iso, match, outliers, confidence,
                 n_clicks, tax_id, text, search, state):
    dff = df[df['species'] == tax_id]
    dff = dff.sort_values(by='dist_pct', ascending=False)

    # parse selected points
    request, data = parse_search_input(dff, state, search, n_clicks, text)
    irows = pandas.Series(index=dff.index, data=False)
    if request is not None:
        irows |= dff[request] == data
    else:  # look for selected points
        if iso:
            irows |= dff['isolation_source'].isin(iso)
        if match:
            irows |= dff['match_species'].isin(match)
        if outliers:
            irows |= dff['is_out'].isin(i == 'Yes' for i in outliers)
        if confidence:
            irows |= dff['confidence'].isin(confidence)
        if selected is not None:
            idx = [i['customdata'] for i in selected['points']]
            if all(i in dff.index for i in idx):  # equivalent to new tax_id
                irows |= dff.index.isin(idx)
    if not irows.any():
        irows[:] = True

    rows = dff[irows].iloc[:MAX_TABLE_RECORDS].copy()  # pull selected rows

    # clean up boolean text and sort by dist
    rows['is_type'] = rows['is_type'].apply(lambda x: 'Yes' if x else '')
    rows['is_out'] = rows['is_out'].apply(lambda x: 'Yes' if x else '')
    rows['dist_pct'] = rows['dist_pct'].apply(lambda x: '{:.2f}'.format(x))

    TABLE_STYLE = {
        'width': 500,
        'padding': 10,
        'border-bottom': '1px solid #ddd',
        'text-align': 'left'}
    cols = ['seqname', 'description', 'is_type', 'dist_pct',
            'is_out', 'match_species', 'match_pct']
    trs = [html.Tr(children=[
        html.Th(children=[c], style=TABLE_STYLE) for c in cols])]
    for r in rows.to_dict('records'):
        tds = []
        for c in cols:
            if c == 'seqname':
                cell = html.A(
                    href=('https://www.ncbi.nlm.nih.gov/nuccore/' +
                          r['version']),
                    children=r[c],
                    target='_blank')
            elif c == 'match_species' and r[c] is not None:
                cell = html.A(
                    href=('https://www.ncbi.nlm.nih.gov/nuccore/' +
                          r['match_version']),
                    children=r[c],
                    target='_blank')
            else:
                cell = r[c]
            tds.append(html.Td(children=[cell], style=TABLE_STYLE))
        trs.append(html.Tr(children=tds))
    return trs


if __name__ == '__main__':
    app.run_server(debug=True)
