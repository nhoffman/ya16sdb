#! /usr/bin/env python3
import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas
import plotly.graph_objs as go
import urllib

from dash.dependencies import Input, State, Output

DEFAULT_GENUS = '547'  # Enterobacter
SEARCH_OPTS = ['seqname', 'accession', 'version',
               'species_name', 'species', 'genus']

app = dash.Dash()
app.title = 'Species Outlier Plots'

df = pandas.read_feather('filter_details.feather')
df = df[~df['x'].isna() & ~df['y'].isna()]

info = df[['x', 'y', 'match_species', 'dist', 'match_pct']].columns

tax = df[['genus', 'genus_name', 'species', 'species_name']]
tax = tax.drop_duplicates().sort_values(by='species_name')
species_genus = dict(tax[['species', 'genus']].drop_duplicates().values)
species_id = dict(tax[['species_name', 'species']].drop_duplicates().values)
genera = tax.groupby(by='genus')
genus_opts = tax[['genus', 'genus_name']].drop_duplicates().values
genus_opts = [{'label': gn, 'value': gi} for gi, gn in genus_opts]


app.layout = html.Div(
    style={'width': 1200},
    children=[
        html.Div(id='state'),
        dcc.Location(id='url', refresh=False),
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
                    'align': 'middle',
                    'display': 'inline-block',
                    'text-align': 'center',
                    'vertical-align': 'middle',
                    'width': '5%'}}),
        html.Div(
            children=[dcc.Dropdown(id='genus-column', options=genus_opts)],
            style={
                'display': 'inline-block',
                'vertical-align': 'middle',
                'width': '45%'}),
        dcc.Markdown(
            children=['**Species**'],
            containerProps={
                'style': {
                    'display': 'inline-block',
                    'text-align': 'center',
                    'vertical-align': 'middle',
                    'width': '6%'}}),
        html.Div(
            children=[dcc.Dropdown(id='species-column')],
            style={
                'vertical-align': 'middle',
                'display': 'inline-block',
                'width': '41%'}),
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
                        dcc.Markdown(children=['**Isolation Source**']),
                        dcc.Markdown(children=['**Match Species**']),
                        dcc.Markdown(children=['**Outliers**']),
                        dcc.Markdown(children=['**Type Strains**']),
                        dcc.Markdown(children=['**Published**']),
                        ],
                    style={
                        'vertical-align': 'middle',
                        'width': '11%',
                        'display': 'inline-block',
                        }),
                dcc.RadioItems(
                    options=[
                        {'value': 'colors-isolation'},
                        {'value': 'colors-match-species'},
                        {'value': 'colors-outlier'},
                        {'value': 'colors-type-strains'},
                        {'value': 'colors-published'}],
                    inputStyle={'height': 15, 'width': 15, 'margin': 11},
                    style={
                        'vertical-align': 'middle',
                        'width': '5%',
                        'display': 'inline-block'}),
                dcc.RadioItems(
                    options=[
                        {'value': 'symbols-isolation-source'},
                        {'value': 'symbols-match-species'},
                        {'value': 'symbols-outlier'},
                        {'value': 'symbols-type-strains'},
                        {'value': 'symbols-published'}],
                    inputStyle={'height': 15, 'width': 15, 'margin': 11},
                    style={
                        'vertical-align': 'middle',
                        'width': '5%',
                        'display': 'inline-block'}),
                html.Div(
                    children=[
                        dcc.Dropdown(
                            id='isolation-source-selection',
                            multi=True),
                        dcc.Dropdown(
                            id='match-species-selection',
                            multi=True),
                        dcc.Dropdown(
                            id='outliers-selection',
                            multi=True),
                        dcc.Dropdown(
                            id='type-strains-selection',
                            multi=True),
                        dcc.Dropdown(
                            id='published-selection',
                            multi=True),
                        ],
                    style={
                        'vertical-align': 'middle',
                        'width': '39%',
                        'display': 'inline-block',
                        }),
                html.Div(
                    children=[
                        dcc.Dropdown(
                            id='isolation-source-visibility',
                            multi=True),
                        dcc.Dropdown(
                            id='match-species-visibility',
                            multi=True),
                        dcc.Dropdown(
                            id='outliers-visibility',
                            multi=True),
                        dcc.Dropdown(
                            id='type-strains-visibility',
                            multi=True),
                        dcc.Dropdown(
                            id='published-visibility',
                            multi=True),
                        ],
                    style={
                        'vertical-align': 'middle',
                        'width': '39%',
                        'display': 'inline-block',
                        }),
                    ],
            style={
                'border': 'thin lightgrey solid',
                'borderRadius': 5,
                'margin': 5,
                'padding': 10,
                'width': '95%'}),
        html.Div(
            children=[
                dcc.Dropdown(
                    id='yaxis-column',
                    options=[{'label': i, 'value': i} for i in info],
                    value='y')],
            style={
                'width': '13%',
                'vertical-align': 'middle',
                'display': 'inline-block'}),
        html.Div(
            children=[
                dcc.Graph(id='plot'),
                html.Div(
                    children=[
                        dcc.Dropdown(
                            id='xaxis-column',
                            options=[{'label': i, 'value': i} for i in info],
                            value='x')],
                    style={
                       'text-align': 'left',
                       'display': 'inline-block',
                       'width': '17%'})],
            style={
                'display': 'inline-block',
                'text-align': 'center',
                'vertical-align': 'middle',
                'width': '77%'}),
        html.Div(
            children=[dcc.Slider(id='year--slider')],
            style={'padding': 15, 'display': 'inline-block', 'width': '95%'}),
        html.Table(
            id='table-div',
            style={
                'border': 'thin lightgrey solid',
                'borderRadius': 5,
                'display': 'inline-block',
                'height': 500,
                'margin': 5,
                'overflow-y': 'scroll',
                'padding': 10,
                'width': '97%'})])


@app.callback(
    Output('genus-column', 'value'),
    [Input('url', 'search'),
     Input('submit-button', 'n_clicks')],
    [State('text-input', 'value'),
     State('state', 'hidden')])
def update_genus_value(search, n_clicks, text, state):
    request, data = parse_search_input(df, state, search, n_clicks, text)
    if request is None:  # no request
        value = DEFAULT_GENUS  # return the default starting genus
    elif request == 'species_name':
        value = species_genus.get(data, DEFAULT_GENUS)
    else:
        value = value = df[df[request] == data].iloc[0]['genus']
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


def parse_search_input(dff, state, search, n_clicks, text):
    '''
    determine where the search input came from
    '''
    request = None
    data = None
    if state is None:  # new instance of the web page
        args = urllib.parse.parse_qs(urllib.parse.urlparse(search).query)
        for o in SEARCH_OPTS:
            if o in args:
                request = o
                data = args[o]
    elif text is not None and state['n_clicks'] < n_clicks:  # button clicked
        data = text.strip()
        for o in SEARCH_OPTS:
            if text in dff[o].values:
                request = o
    return request, data


# Setting the url search does not seem to work at the moment in addition
# to the circular logic this introduces:
# https://github.com/plotly/dash-core-components/issues/44
# @app.callback(
#     Output('url', 'search'),
#     [Input('species-column', 'value')])
# def update_url_search(value):
#     print(value)
#     return '?' + urllib.parse.urlencode({'species_id': value})


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
    Output('type-strains-selection', 'options'),
    [Input('species-column', 'value')])
def update_type_strains_selection(tax_id):
    '''
    '''
    dff = df[df['species'] == tax_id]
    types = dff['is_type'].apply(lambda x: 'Yes' if x else 'No')
    options = []
    for k, v in types.value_counts().items():
        options.append({'label': '{} ({})'.format(k, v), 'value': k})
    return options


@app.callback(
    Output('type-strains-visibility', 'options'),
    [Input('species-column', 'value')])
def update_type_strains_visiblity(tax_id):
    '''
    '''
    dff = df[df['species'] == tax_id]
    types = dff['is_type'].apply(lambda x: 'Yes' if x else 'No')
    options = []
    for k, v in types.value_counts().items():
        options.append({'label': '{} ({})'.format(k, v), 'value': k})
    return options


@app.callback(
    Output('published-selection', 'options'),
    [Input('species-column', 'value')])
def update_published_selection(tax_id):
    '''
    '''
    dff = df[df['species'] == tax_id]
    published = dff['is_type'].apply(lambda x: 'Yes' if x else 'No')
    options = []
    for k, v in published.value_counts().items():
        options.append({'label': '{} ({})'.format(k, v), 'value': k})
    return options


@app.callback(
    Output('published-visibility', 'options'),
    [Input('species-column', 'value')])
def update_published_visibility(tax_id):
    '''
    '''
    dff = df[df['species'] == tax_id]
    published = dff['is_type'].apply(lambda x: 'Yes' if x else 'No')
    options = []
    for k, v in published.value_counts().items():
        options.append({'label': '{} ({})'.format(k, v), 'value': k})
    return options


@app.callback(
    Output('isolation-source-selection', 'value'),
    [Input('species-column', 'value')])
def update_isolation_source_selection_value(tax_id):
    return None  # clear values


@app.callback(
    Output('isolation-source-visibility', 'value'),
    [Input('species-column', 'value')])
def update_isolation_source_visiblity_value(tax_id):
    return None  # clear values


@app.callback(
    Output('match-species-selection', 'value'),
    [Input('species-column', 'value')])
def update_match_species_selection_value(tax_id):
    return None  # clear values


@app.callback(
    Output('match-species-visibility', 'value'),
    [Input('species-column', 'value')])
def update_match_species_visibility_value(tax_id):
    return None  # clear values


@app.callback(
    Output('outliers-selection', 'value'),
    [Input('species-column', 'value')])
def update_outliers_selection_value(tax_id):
    return None  # clear values


@app.callback(
    Output('outliers-visibility', 'value'),
    [Input('species-column', 'value')])
def update_outliers_visibility_value(tax_id):
    return None  # clear values


@app.callback(
    Output('type-strains-selection', 'value'),
    [Input('species-column', 'value')])
def update_type_strain_selection_value(tax_id):
    return None  # clear values


@app.callback(
    Output('type-strains-visibility', 'value'),
    [Input('species-column', 'value')])
def update_type_strain_visibility_value(tax_id):
    return None  # clear values


@app.callback(
    Output('published-selection', 'value'),
    [Input('species-column', 'value')])
def update_published_selection_value(tax_id):
    return None  # clear values


@app.callback(
    Output('published-visibility', 'value'),
    [Input('species-column', 'value')])
def update_published_visibility_value(tax_id):
    return None  # clear values


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
     Input('type-strains-visibility', 'value'),
     Input('type-strains-selection', 'value'),
     Input('published-visibility', 'value'),
     Input('published-selection', 'value'),
     Input('submit-button', 'n_clicks')],
    [State('state', 'hidden'),
     State('text-input', 'value'),
     State('url', 'search')])
def update_graph(tax_id, xaxis_column_name, yaxis_column_name, year_value,
                 viso_source, siso_source,
                 vmatch, smatch,
                 vout, sout,
                 vtypes, stypes,
                 vpubs, spubs,
                 n_clicks, state, text, search):
    # pprint.pprint(locals())
    dff = df[df['species'] == tax_id]
    dff = dff[dff['modified_date'] <= str(year_value+1)]
    if viso_source:
        dff = dff[dff['isolation_source'].isin(viso_source)]
    if vmatch:
        dff = dff[dff['match_species'].isin(vmatch)]
    if vout:
        dff = dff[dff['is_out'].isin(set(i == 'Yes' for i in vout))]
    if vtypes:
        dff = dff[dff['is_type'].isin(set(i == 'Yes' for i in vtypes))]
    if vpubs:
        # pubmed_id
        dff = dff[dff['is_type'].isin(set(i == 'Yes' for i in vpubs))]

    # decide if we should allow plot to calculate axes ranges
    if state is None:
        x_range = None
        y_range = None
    elif tax_id != state['tax_id']:
        x_range = None
        y_range = None
    elif state['xaxis'] != xaxis_column_name:
        x_range = None
        y_range = None
    elif state['yaxis'] != yaxis_column_name:
        x_range = None
        y_range = None
    else:
        x_range = state['xrange']
        y_range = state['yrange']

    dff['text'] = dff.apply(
        lambda x: 'seqname: {seqname}<br>'
                  'accession: {version}<br>'
                  'modified_date: {modified_date}<br>'
                  'isolation source: {isolation_source}'.format(**x),
        axis='columns')

    # decide selected points
    dff['selected'] = False
    request, data = parse_search_input(dff, state, search, n_clicks, text)
    if request is not None:
        dff.loc[dff[request] == data, 'selected'] = True
    if siso_source:
        dff.loc[dff['isolation_source'].isin(siso_source), 'selected'] = True
    if smatch:
        dff.loc[dff['match_species'].isin(smatch), 'selected'] = True
    if sout:
        is_out = dff['is_out'].isin(set(i == 'Yes' for i in sout))
        dff.loc[is_out, 'selected'] = True
    if stypes:
        is_type = dff['is_type'].isin(set(i == 'Yes' for i in stypes))
        dff.loc[is_type, 'selected'] = True
    if spubs:
        is_pub = dff['isolation_source'].isin(set(i == 'Yes' for i in spubs))
        dff.loc[is_pub, 'selected'] = True

    inliers = ~dff['is_out'] & ~dff['is_type']
    outliers = dff['is_out'] & ~dff['is_type']
    types = dff['is_type']

    dff.loc[inliers | outliers, 'symbol'] = 'circle'
    dff.loc[types, 'symbol'] = 'triangle-up'
    dff.loc[inliers, 'color'] = 'lightblue'
    dff.loc[outliers, 'color'] = 'lightgreen'
    dff.loc[types, 'color'] = 'darkred'

    # sort markers to move non-circle markers in front
    cat_type = pandas.api.types.CategoricalDtype(
        categories=['circle', 'triangle-up'])
    dff['symbol'] = dff['symbol'].astype(cat_type)
    dff = dff.sort_values(by='symbol', ascending=True)

    dff['iselected'] = range(len(dff))

    figure = {
        'data': [
            go.Scatter(
                customdata=dff.index,
                marker={
                    'symbol': dff['symbol'],
                    'color': dff['color']},
                mode='markers',
                selected={'marker': {'size': 15, 'opacity': 0.7}},
                selectedpoints=dff[dff['selected']]['iselected'],
                unselected={'marker': {'size': 5, 'opacity': 0.4}},
                text=dff['text'],
                x=dff[xaxis_column_name],
                y=dff[yaxis_column_name]),
            ],
        'layout': go.Layout(
            xaxis={
                'scaleanchor': 'y',  # square grid ratio
                'title': xaxis_column_name,
                'type': 'linear',
                'showgrid': False,
                'range': x_range
            },
            yaxis={
                'scaleanchor': 'x',  # square grid ratio
                'title': yaxis_column_name,
                'type': 'linear',
                'showgrid': False,
                'range': y_range
            },
            dragmode='select',  # default tool from modebar
            showlegend=False,
            legend={'orientation': 'v'},
            hovermode='closest',
            height=900,
            width=900,
            title='{}, tax_id {} outliers: {} ({:.0%}), total: {}'.format(
                dff.iloc[0]['species_name'],
                tax_id,
                len(outliers),
                (len(outliers) / len(dff)),
                len(dff)),
        )
    }
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
        'tax_id': tax_id,  # if tax_id changes we reset everything
        'xrange': figure['layout']['xaxis']['range'],  # preserve axes ranges
        'yrange': figure['layout']['yaxis']['range'],
        'xaxis': xaxis,
        'yaxis': yaxis}


@app.callback(
    Output('table-div', 'children'),
    [Input('plot', 'figure'),
     Input('plot', 'selectedData')],
    [State('submit-button', 'n_clicks'),
     State('species-column', 'value'),
     State('text-input', 'value'),
     State('url', 'search'),
     State('state', 'hidden')])
def update_table(_, selected, n_clicks, tax_id, text, search, state):
    # pprint.pprint(locals())
    dff = df[df['species'] == tax_id]
    request, data = parse_search_input(dff, state, search, n_clicks, text)
    if request is not None:
        rows = dff[dff[request] == data]
    elif selected is not None:
        idx = [i['customdata'] for i in selected['points']]
        rows = dff.loc[idx]
    else:
        # raise exception?
        rows = dff[dff['is_out']]  # outliers

    # TODO: rebuild feather file to get match_version from seq_info file..
    rows = rows.copy()  # avoid SettingWithCopyWarning

    # clean up boolean text and sort by dist
    rows['is_type'] = rows['is_type'].apply(lambda x: 'Yes' if x else '')
    rows['is_out'] = rows['is_out'].apply(lambda x: 'Yes' if x else '')
    rows = rows.sort_values(by='dist', ascending=False)

    TABLE_STYLE = {
        'width': 500,
        'padding': 10,
        'border-bottom': '1px solid #ddd',
        'text-align': 'left'}

    cols = ['seqname', 'description', 'is_type', 'dist',
            'is_out', 'match_species', 'match_pct']

    trs = [html.Tr(children=[
        html.Th(children=[c], style=TABLE_STYLE) for c in cols])]
    for _, r in rows.iterrows():
        r = r.to_dict()
        tds = []
        for c in cols:
            if c == 'seqname':
                cell = html.A(
                    href=('https://www.ncbi.nlm.nih.gov/nuccore/' +
                          r['version']),
                    children=r[c],
                    target='_blank')
            elif c == 'match_species':
                cell = html.A(
                    href='https://www.ncbi.nlm.nih.gov/nuccore/' + r[c],
                    children=r[c],
                    target='_blank')
            else:
                cell = r[c]
            tds.append(html.Td(children=[cell], style=TABLE_STYLE))
        trs.append(html.Tr(children=tds))
    return trs


if __name__ == '__main__':
    app.run_server(port=8005)
