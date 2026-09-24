import dash
from dash.dependencies import Input, State, Output
import os
import pandas
import urllib

import shared

dash.register_page(__name__, path='/', name='Plot')

COLORS = shared.COLORS
DEFAULT_COLOR = shared.DEFAULT_COLOR
DEFAULT_SHAPE = shared.DEFAULT_SHAPE
DEFAULT_GENUS = shared.DEFAULT_GENUS
DEFAULT_SPECIES = shared.DEFAULT_SPECIES
DEFAULT_Y = shared.DEFAULT_Y
DEFAULT_X = shared.DEFAULT_X
LEGEND_OTHER = shared.LEGEND_OTHER
MAX_TABLE_RECORDS = shared.MAX_TABLE_RECORDS
NO_DATA = shared.NO_DATA
SEARCH_OPTS = shared.SEARCH_OPTS
SHAPES = shared.SHAPES


def layout():
    axes = ['confidence', 'dist_pct', 'x', 'y',
            'type_classification', 'rank_order']
    return dash.html.Div(
        style={'width': 1175},
        children=[
            dash.dcc.Store(id='state'),
            dash.dcc.Location(id='url', refresh=False),
            dash.dcc.Markdown(
                children='Loading...',
                id='download-date',
                style={
                    'textAlign': 'right',
                    'fontStyle': 'italic',
                    'height': 0,
                    'width': '100%'}),
            dash.html.Div(
                children=[
                    dash.dcc.Input(type='text', id='text-input'),
                    dash.html.Button(
                        id='submit-button',
                        n_clicks=0,
                        children='Search')]),
            dash.dcc.Markdown(
                children=['**Genus**'],
                style={
                    'display': 'inline-block',
                    'textAlign': 'center',
                    'verticalAlign': 'middle',
                    'width': '5%'}),
            dash.html.Div(
                children=[
                    dash.dcc.Dropdown(
                        id='genus-column',
                        clearable=False)],
                style={
                    'display': 'inline-block',
                    'verticalAlign': 'middle',
                    'width': '40%'}),
            dash.dcc.Markdown(
                children=['**Species**'],
                style={
                    'display': 'inline-block',
                    'textAlign': 'center',
                    'verticalAlign': 'middle',
                    'width': '6%'}),
            dash.html.Div(
                children=[
                    dash.dcc.Dropdown(id='species-column', clearable=False)],
                style={
                    'display': 'inline-block',
                    'verticalAlign': 'middle',
                    'width': '48%'}),
            dash.html.Div(
                children=[
                    dash.html.Div(
                        style={'width': '11%', 'display': 'inline-block'}),
                    dash.dcc.Markdown(
                        children=['**Color**'],
                        style={
                            'width': '5%',
                            'display': 'inline-block'}),
                    dash.dcc.Markdown(
                        children=['**Shape**'],
                        style={
                            'width': '10%',
                            'display': 'inline-block'}),
                    dash.dcc.Markdown(
                        children=['**Selection**'],
                        style={
                             'width': '39%',
                             'display': 'inline-block'}),
                    dash.dcc.Markdown(
                        children=['**Visibility**'],
                        style={
                                'width': '35%',
                                'display': 'inline-block'}),
                    dash.html.Div(
                        children=[
                            dash.dcc.Markdown(children=['**Outliers**']),
                            dash.dcc.Markdown(children=['**Confidence**']),
                            dash.dcc.Markdown(
                                children=['**Type Classification**']),
                            dash.dcc.Markdown(children=['**ANI Tax Check**']),
                            dash.dcc.Markdown(children=['**ANI Species**']),
                            dash.dcc.Markdown(
                                children=['**Isolation Source**']),
                            ],
                        style={
                            'verticalAlign': 'middle',
                            'width': '11%',
                            'display': 'inline-block',
                            }),
                    dash.dcc.RadioItems(
                        id='color-items',
                        options=[
                            {'label': '', 'value': 'is_out'},
                            {'label': '', 'value': 'confidence'},
                            {'label': '', 'value': 'type_classification'},
                            {'label': '', 'value': 'taxonomy-check-status'},
                            {'label': '', 'value': 'best-match-species-name'},
                            {'label': '', 'value': 'isolation_source'}],
                        inputStyle={'height': 15, 'width': 15, 'margin': 11},
                        style={
                            'verticalAlign': 'middle',
                            'width': '5%',
                            'display': 'inline-block'}),
                    dash.dcc.RadioItems(
                        id='shape-items',
                        options=[
                            {'label': '', 'value': 'is_out'},
                            {'label': '', 'value': 'confidence'},
                            {'label': '', 'value': 'type_classification'},
                            {'label': '', 'value': 'taxonomy-check-status'},
                            {'label': '', 'value': 'best-match-species-name'},
                            {'label': '', 'value': 'isolation_source'}],
                        inputStyle={'height': 15, 'width': 15, 'margin': 11},
                        style={
                            'verticalAlign': 'middle',
                            'width': '5%',
                            'display': 'inline-block'}),
                    dash.html.Div(
                        children=[
                            dash.dcc.Dropdown(
                                id='outliers-selection',
                                multi=True),
                            dash.dcc.Dropdown(
                                id='confidence-selection',
                                multi=True),
                            dash.dcc.Dropdown(
                                id='type-classification-selection',
                                multi=True),
                            dash.dcc.Dropdown(
                                id='taxcheck-selection',
                                multi=True),
                            dash.dcc.Dropdown(
                                id='ani-species-selection',
                                multi=True),
                            dash.dcc.Dropdown(
                                id='isolation-source-selection',
                                multi=True)],
                        style={
                            'verticalAlign': 'middle',
                            'width': '39%',
                            'display': 'inline-block'}),
                    dash.html.Div(
                        children=[
                            dash.dcc.Dropdown(
                                id='outliers-visibility',
                                multi=True),
                            dash.dcc.Dropdown(
                                id='confidence-visibility',
                                multi=True),
                            dash.dcc.Dropdown(
                                id='type-classification-visibility',
                                multi=True),
                            dash.dcc.Dropdown(
                                id='taxcheck-visibility',
                                multi=True),
                            dash.dcc.Dropdown(
                                id='ani-species-visibility',
                                multi=True),
                            dash.dcc.Dropdown(
                                id='isolation-source-visibility',
                                multi=True)],
                        style={
                            'verticalAlign': 'middle',
                            'width': '39%',
                            'display': 'inline-block'})],
                style={
                    'border': 'thin lightgrey solid',
                    'borderRadius': 5,
                    'margin': 5,
                    'padding': 10,
                    'width': '97%'}),
            dash.html.Div(
                children=[dash.dcc.Slider(id='year--slider', marks=None)],
                style={'margin': 15, 'width': '95%'}),
            dash.dcc.Markdown(
                children=['**Axes**'],
                style={
                    'display': 'inline-block',
                    'textAlign': 'center',
                    'verticalAlign': 'middle',
                    'width': '4%'}),
            dash.html.Div(
                children=[
                    dash.dcc.Dropdown(
                        clearable=False,
                        id='yaxis-column',
                        options=[{'label': i, 'value': i} for i in axes])],
                style={
                    'width': '14%',
                    'display': 'inline-block',
                    'verticalAlign': 'middle'}),
            dash.html.Div(
                children=[
                    dash.dcc.Dropdown(
                        clearable=False,
                        id='xaxis-column',
                        options=[{'label': i, 'value': i} for i in axes])],
                style={
                    'width': '14%',
                    'display': 'inline-block',
                    'verticalAlign': 'middle'}),
            dash.dcc.Graph(id='plot'),
            dash.html.Table(
                id='table-div',
                style={
                    'border': 'thin lightgrey solid',
                    'borderRadius': 5,
                    'display': 'inline-block',
                    'height': 413,
                    'margin': 5,
                    'overflowY': 'scroll',
                    'padding': 10,
                    'width': '97%'}),
            dash.html.Div(
                children=[
                    dash.html.Div(style={'flex': '1'}),
                    dash.html.A(
                        os.environ.get(
                            'VERSION',
                            'github.com/nhoffman/ya16sdb/dash'),
                        href='https://github.com/nhoffman/'
                             'ya16sdb/pkgs/container/ya16sdb-app',
                        style={
                            'textDecoration': 'none',
                            'textAlign': 'center',
                            'flex': '1'},
                        target='_blank'),
                    dash.html.Div(
                        dash.dcc.Link(
                            'Data Table',
                            href='/datatable',
                            style={'fontSize': '16px'}),
                        style={
                            'flex': '1',
                            'textAlign': 'right'}),
                ],
                style={
                    'display': 'flex',
                    'alignItems': 'center',
                    'padding': '10px 0',
                    'width': '97%'}),
                ])


def assign_hover_text(s):
    '''
    assign hover text to outliers for now
    '''
    text = None
    if s['is_out']:
        text = ('seqname: {seqname}<br>'
                'accession: {version}<br>'
                'modified_date: {modified_date}<br>'
                'type_classification: {type_classification}<br>'
                'ani_tax_check: {taxonomy-check-status}<br>'
                'ani_species: {best-match-species-name}<br>'
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


@dash.callback(
    Output('yaxis-column', 'value'),
    [Input('url', 'search')])
def update_yaxis_value(search):
    args = urllib.parse.parse_qs(urllib.parse.urlparse(search).query)
    return args.get('y', [DEFAULT_Y])[0]


@dash.callback(
    Output('xaxis-column', 'value'),
    [Input('url', 'search')])
def update_xaxis_value(search):
    args = urllib.parse.parse_qs(urllib.parse.urlparse(search).query)
    return args.get('x', [DEFAULT_X])[0]


@dash.callback(Output('download-date', 'children'), Input('state', 'data'))
def update_download_date(state):
    if state is None:
        shared.set_global_data()
    return str(shared.df['download_date'].max().strftime('%A, %B %d, %Y'))


@dash.callback(Output('genus-column', 'options'), Input('state', 'data'))
def update_genus_options(state):
    opts = shared.tax[['genus', 'genus_name']].drop_duplicates().values
    opts = [{'label': gn, 'value': gi} for gi, gn in opts]
    return opts


@dash.callback(
    Output('genus-column', 'value'),
    [Input('url', 'search'),
     Input('submit-button', 'n_clicks')],
    [State('text-input', 'value'),
     State('state', 'data')])
def update_genus_value(search, n_clicks, text, state):
    if state is None:
        shared.set_global_data()
    request, data = parse_search_input(
        shared.df, state, search, n_clicks, text)
    if request is None:
        value = DEFAULT_GENUS
    elif request == 'species_name':
        value = shared.species_to_genus.get(data, DEFAULT_GENUS)
    else:
        value = shared.df[shared.df[request] == data].iloc[0]['genus']
    return value


@dash.callback(
    Output('species-column', 'options'),
    [Input('genus-column', 'value')])
def update_species_options(tax_id):
    group = shared.genera.get_group(tax_id)
    group = group[['species', 'species_name']]
    return [{'label': sn, 'value': si} for si, sn in group.values]


@dash.callback(
    Output('species-column', 'value'),
    [Input('species-column', 'options'),
     Input('url', 'search'),
     Input('submit-button', 'n_clicks')],
    [State('text-input', 'value'),
     State('state', 'data'),
     State('genus-column', 'value')])
def update_species_value(options, search, n_clicks, text, state, tax_id):
    dff = shared.df[shared.df['genus'] == tax_id]
    request, data = parse_search_input(dff, state, search, n_clicks, text)
    if request is None:  # no request
        if state:
            value = options[0]['value']  # return first item in dropdown
        else:
            value = DEFAULT_SPECIES
    elif request == 'species_name':
        value = shared.species_to_id.get(data, options[0]['value'])
    else:
        value = dff[dff[request] == data].iloc[0]['species']
    return value


@dash.callback(
    Output('year--slider', 'value'),
    Output('year--slider', 'min'),
    Output('year--slider', 'max'),
    Output('year--slider', 'marks'),
    [Input('species-column', 'value'),
     Input('isolation-source-visibility', 'value')])
def update_slider(tax_id, iso_values):
    dff = shared.get_species(tax_id)
    if iso_values:
        dff = dff[dff['isolation_source'].isin(iso_values)]
    dates = dff['modified_date']
    min_year = dates.min().year
    max_year = dates.max().year
    years = dates.apply(lambda x: x.year).unique()
    marks = {str(y): str(y) for y in years}
    return max_year, min_year, max_year, marks


def _make_options(series):
    return [{'label': '{} ({})'.format(k, v), 'value': k}
            for k, v in series.value_counts().items()]


@dash.callback(
    Output('isolation-source-selection', 'options'),
    Output('isolation-source-visibility', 'options'),
    Output('type-classification-selection', 'options'),
    Output('type-classification-visibility', 'options'),
    Output('ani-species-selection', 'options'),
    Output('ani-species-visibility', 'options'),
    Output('taxcheck-selection', 'options'),
    Output('taxcheck-visibility', 'options'),
    Output('outliers-selection', 'options'),
    Output('outliers-visibility', 'options'),
    Output('confidence-selection', 'options'),
    Output('confidence-visibility', 'options'),
    [Input('species-column', 'value')])
def update_dropdown_options(tax_id):
    dff = shared.get_species(tax_id)
    iso = _make_options(dff['isolation_source'])
    tc = _make_options(dff['type_classification'])
    ani = _make_options(dff['best-match-species-name'])
    txc = _make_options(dff['taxonomy-check-status'])
    out = _make_options(dff['is_out'].apply(lambda x: 'Yes' if x else 'No'))
    conf = _make_options(dff['confidence'])
    return iso, iso, tc, tc, ani, ani, txc, txc, out, out, conf, conf


def parse_multi(state, search, option):
    value = None
    if state is None:
        args = urllib.parse.parse_qs(urllib.parse.urlparse(search).query)
        if option in args:
            value = args[option][0].split(',')
    return value


@dash.callback(
    Output('isolation-source-selection', 'value'),
    Output('isolation-source-visibility', 'value'),
    Output('type-classification-selection', 'value'),
    Output('type-classification-visibility', 'value'),
    Output('ani-species-selection', 'value'),
    Output('ani-species-visibility', 'value'),
    Output('taxcheck-selection', 'value'),
    Output('taxcheck-visibility', 'value'),
    Output('outliers-selection', 'value'),
    Output('outliers-visibility', 'value'),
    Output('confidence-selection', 'value'),
    Output('confidence-visibility', 'value'),
    [Input('species-column', 'value')],
    [State('state', 'data'),
     State('url', 'search')])
def update_dropdown_values(_, state, search):
    return (
        parse_multi(state, search, 'selection_isolation_source'),
        parse_multi(state, search, 'visibility_isolation_source'),
        parse_multi(state, search, 'selection_type_classification'),
        parse_multi(state, search, 'visibility_type_classification'),
        parse_multi(state, search, 'selection_ani_species'),
        parse_multi(state, search, 'visibility_ani_species'),
        parse_multi(state, search, 'selection_taxcheck'),
        parse_multi(state, search, 'visibility_taxcheck'),
        parse_multi(state, search, 'selection_is_out'),
        parse_multi(state, search, 'visibility_is_out'),
        parse_multi(state, search, 'selection_confidence'),
        parse_multi(state, search, 'visibility_confidence'),
    )


@dash.callback(
    Output('color-items', 'value'),
    [Input('url', 'search')])
def update_color_items(search):
    args = urllib.parse.parse_qs(urllib.parse.urlparse(search).query)
    return args.get('color', [DEFAULT_COLOR])[0]


@dash.callback(
    Output('shape-items', 'value'),
    [Input('url', 'search')])
def update_symbol_items(search):
    args = urllib.parse.parse_qs(urllib.parse.urlparse(search).query)
    return args.get('shape', [DEFAULT_SHAPE])[0]


@dash.callback(
    Output('plot', 'figure'),
    [Input('species-column', 'value'),
     Input('xaxis-column', 'value'),
     Input('yaxis-column', 'value'),
     Input('year--slider', 'value'),
     Input('isolation-source-visibility', 'value'),
     Input('isolation-source-selection', 'value'),
     Input('type-classification-visibility', 'value'),
     Input('type-classification-selection', 'value'),
     Input('taxcheck-visibility', 'value'),
     Input('taxcheck-selection', 'value'),
     Input('ani-species-visibility', 'value'),
     Input('ani-species-selection', 'value'),
     Input('outliers-visibility', 'value'),
     Input('outliers-selection', 'value'),
     Input('confidence-visibility', 'value'),
     Input('confidence-selection', 'value'),
     Input('color-items', 'value'),
     Input('shape-items', 'value'),
     Input('submit-button', 'n_clicks')],
    [State('state', 'data'),
     State('text-input', 'value'),
     State('url', 'search')])
def update_graph(tax_id, xaxis, yaxis, year_value,
                 viso_source, siso_source,
                 vmatch, smatch,
                 vtaxcheck, staxcheck,
                 vani, sani,
                 vout, sout,
                 vconf, sconf,
                 color, symbol,
                 n_clicks, state, text, search):
    dff = shared.get_species(tax_id)
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
        dff = dff[dff['type_classification'].isin(vmatch)]
    if vani:
        dff = dff[dff['best-match-species-name'].isin(vani)]
    if vout:
        dff = dff[dff['is_out'].isin(set(i == 'Yes' for i in vout))]
    if vconf:
        dff = dff[dff['confidence'].isin(vconf)]
    if vtaxcheck:
        dff = dff[dff['taxonomy-check-status'].isin(vtaxcheck)]

    # decide selected points
    dff['selected'] = False
    request, data = parse_search_input(dff, state, search, n_clicks, text)
    if request in ['seqname', 'accession', 'version']:
        dff.loc[dff[request] == data, 'selected'] = True
    if siso_source:
        dff.loc[dff['isolation_source'].isin(siso_source), 'selected'] = True
    if smatch:
        dff.loc[dff['type_classification'].isin(smatch), 'selected'] = True
    if sani:
        dff.loc[dff['best-match-species-name'].isin(sani), 'selected'] = True
    if sout:
        is_out = dff['is_out'].isin(set(i == 'Yes' for i in sout))
        dff.loc[is_out, 'selected'] = True
    if sconf:
        dff.loc[dff['confidence'].isin(sconf), 'selected'] = True
    if staxcheck:
        taxcheck = dff['taxonomy-check-status'].isin(staxcheck)
        dff.loc[taxcheck, 'selected'] = True

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
                'b': 200 if xaxis == 'type_classification' else None,
                'l': 200 if yaxis == 'type_classification' else None
                },
            'showlegend': True,
            'title': '{}, tax_id {} outliers: {} ({:.0%}), total: {}'.format(
                dff.iloc[0]['species_name'],
                tax_id,
                len(outliers),
                (len(outliers) / len(dff)),
                len(dff)),
            'xaxis': {
                'title': xaxis,
                'showgrid': True if xaxis == 'type_classification' else False,
                'range': x_range
            },
            'yaxis': {
                'title': yaxis,
                'showgrid': True if yaxis == 'type_classification' else False,
                'range': y_range
            }
        }}
    return figure


@dash.callback(
    Output('state', 'data'),
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


@dash.callback(
    Output('table-div', 'children'),
    [Input('plot', 'selectedData'),
     Input('isolation-source-selection', 'value'),
     Input('type-classification-selection', 'value'),
     Input('ani-species-selection', 'value'),
     Input('outliers-selection', 'value'),
     Input('confidence-selection', 'value')],
    [State('submit-button', 'n_clicks'),
     State('species-column', 'value'),
     State('text-input', 'value'),
     State('url', 'search'),
     State('state', 'data')])
def update_table(selected, iso, match, ani, outliers, confidence,
                 n_clicks, tax_id, text, search, state):
    dff = shared.get_species(tax_id)
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
            irows |= dff['type_classification'].isin(match)
        if ani:
            irows |= dff['best-match-species-name'].isin(ani)
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
        'borderBottom': '1px solid #ddd',
        'textAlign': 'left'}
    cols = ['seqname', 'is_type', 'dist_pct', 'outlier', 'type_classification',
            'ani_tax_check', 'ani_species', 'ani_pct']
    trs = [dash.html.Tr(children=[
        dash.html.Th(children=[c], style=TABLE_STYLE) for c in cols])]
    for r in rows.to_dict('records'):
        tds = []
        for c in cols:
            if c == 'seqname':
                cell = dash.html.A(
                    href=('https://www.ncbi.nlm.nih.gov/nuccore/' +
                          r['version']),
                    children=r[c],
                    target='_blank',
                    title=r['description'])
            elif c == 'ani_species':
                if pandas.notna(r['best-match-type-assembly']):
                    cell = dash.html.A(
                        href=('https://www.ncbi.nlm.nih.gov/assembly/' +
                              r['best-match-type-assembly']),
                        children=r['best-match-species-name'],
                        target='_blank')
                else:
                    cell = None
            elif c == 'ani_pct':
                cell = r['best-match-type-ANI']
            elif c == 'ani_tax_check':
                if r['taxonomy-check-status'] == NO_DATA:
                    cell = None
                else:
                    cell = r['taxonomy-check-status']
            elif c == 'outlier':
                cell = r['is_out']
            else:
                cell = r[c]
            tds.append(dash.html.Td(children=[cell], style=TABLE_STYLE))
        trs.append(dash.html.Tr(children=tds))
    return trs
