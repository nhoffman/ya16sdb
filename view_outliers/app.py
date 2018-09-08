#! /usr/bin/env python3
import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas
import plotly.graph_objs as go
import urllib

from dash.dependencies import Input, State, Output

DEFAULT_GENUS = '547'  # Enterobacter

app = dash.Dash()
app.title = 'Species Outlier Plots'

df = pandas.read_feather('filter_details.feather')
df = df[~df['x'].isna() & ~df['y'].isna()]

info = df[['x', 'y', 'dist']].columns

tax = df[['genus', 'genus_name', 'species', 'species_name']]
tax = tax.drop_duplicates().sort_values(by='species_name')
species_genus = dict(tax[['species', 'genus']].drop_duplicates().values)
species_id = dict(tax[['species_name', 'species']].drop_duplicates().values)
genera = tax.groupby(by='genus')
genus_opts = tax[['genus', 'genus_name']].drop_duplicates().values
genus_opts = [{'label': gn, 'value': gi} for gi, gn in genus_opts]


app.layout = html.Div(
    # style={'width': '1150'},
    children=[
        html.Div(id='state'),
        html.Div(
            children=[
                html.Div(
                    children=[
                        dcc.Location(id='url', refresh=False),
                        html.Div(
                            children=[
                                dcc.Input(
                                    type='text',
                                    id='text-input'),
                                html.Button(
                                    id='submit-button',
                                    n_clicks=0,
                                    children='Search')]),
                        dcc.Markdown('**Genus**'),
                        dcc.Dropdown(id='genus-column', options=genus_opts),
                        dcc.Markdown('**Species**'),
                        dcc.Dropdown(id='species-column'),
                        dcc.Markdown('**Isolation Source**'),
                        dcc.Dropdown(id='isolation-source-column', multi=True),
                        dcc.Markdown('**Axes**'),
                        dcc.Dropdown(
                            id='xaxis-column',
                            options=[{'label': i, 'value': i} for i in info],
                            value='x'),
                        dcc.Dropdown(
                            id='yaxis-column',
                            options=[{'label': i, 'value': i} for i in info],
                            value='y')],
                    style={'width': '30%', 'display': 'inline-block'}),
                html.Div(
                    children=[
                        dcc.Graph(id='plot'),
                        dcc.RadioItems(
                            id='radio-items',
                            options=[
                                {'label': 'All', 'value': 'all'},
                                {'label': 'Type Strain', 'value': 'is_type'},
                                {'label': 'Outliers', 'value': 'is_out'}],
                            value='all')],
                    style={
                        'width': '70%',
                        'display': 'inline-block',
                        'float': 'right'})]),
        html.Div(
            children=[dcc.Slider(id='year--slider')],
            style={
                'display': 'inline-block',
                'padding': '10',
                'width': '97%'}),
        html.Div(
            id='table-div',
            style={'height': '500', 'overflow-y': 'scroll'},
            children=[
                html.Table(
                    children=[
                        html.Tr(
                            children=[
                                html.Td()])])])
        ])


@app.callback(
    Output('genus-column', 'value'),
    [Input('url', 'search'),
     Input('submit-button', 'n_clicks')],
    [State('text-input', 'value'),
     State('state', 'hidden')])
def update_genus_value(search, n_clicks, text, state):
    if state and state['n_clicks'] < n_clicks:  # button clicked?
        text = text.strip()
        if text in df['seqname'].values:
            value = df[df['seqname'] == text].iloc[0]['genus']
        elif text in df['accession'].values:
            value = df[df['accession'] == text].iloc[0]['genus']
        elif text in df['version'].values:
            value = df[df['version'] == text].iloc[0]['genus']
        elif text in df['genus'].values:
            value = text
        elif text in df['species'].values:
            value = species_genus[text]
        elif text in df['species_name'].values:
            value = species_genus.get(species_id[text], DEFAULT_GENUS)
        else:
            value = DEFAULT_GENUS
    else:  # url
        args = urllib.parse.parse_qs(urllib.parse.urlparse(search).query)
        if 'seqname' in args:
            seqname = args['seqname'][0]
            value = df[df['seqname'] == seqname].iloc[0]['genus']
        elif 'accession' in args:
            acc = args['accession'][0]
            value = df[df['accession'] == acc].iloc[0]['genus']
        elif 'version' in args:
            acc = args['version'][0]
            value = df[df['version'] == acc].iloc[0]['genus']
        elif 'species_name' in args:
            species = species_id.get(args['species_name'][0], None)
            value = species_genus.get(species, DEFAULT_GENUS)
        elif 'species_id' in args:
            value = species_genus.get(args['species_id'][0], DEFAULT_GENUS)
        else:
            value = DEFAULT_GENUS
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
    if state is None:
        args = urllib.parse.parse_qs(urllib.parse.urlparse(search).query)
        if 'seqname' in args:
            seqname = args['seqname'][0]
            value = dff[dff['seqname'] == seqname].iloc[0]['species']
        elif 'accession' in args:
            acc = args['accession'][0]
            value = dff[dff['accession'] == acc].iloc[0]['species']
        elif 'species_name' in args:
            name = args['species_name'][0]
            value = species_id.get(name, options[0]['value'])
        elif 'species_id' in args:
            value = args['species_id'][0]
        else:
            value = options[0]['value']
    elif state['n_clicks'] < n_clicks:  # button clicked
        text = text.strip()
        if text in dff['seqname'].values:
            value = dff[dff['seqname'] == text].iloc[0]['species']
        elif text in dff['accession'].values:
            value = dff[dff['accession'] == text].iloc[0]['species']
        elif text in dff['version'].values:
            value = dff[dff['version'] == text].iloc[0]['species']
        elif text in dff['species'].values:
            value = text
        elif text in dff['species_name'].values:
            value = species_id.get(text, options[0]['value'])
        else:
            value = options[0]['value']
    else:
        value = options[0]['value']
    return value


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
     Input('isolation-source-column', 'value'),
     Input('radio-items', 'value')])
def update_slider_value(tax_id, iso_values, radio):
    '''
    Reset the slider to generate the axis ranges with all data points.
    This will also help avoid confusion when users select a new species
    to view.
    '''
    dff = df[df['species'] == tax_id]
    if radio == 'is_type':
        dff = dff[dff['is_type']]
    elif radio == 'is_out':
        dff = dff[dff['is_out']]
    if iso_values:
        dff = dff[dff['isolation_source'].isin(iso_values)]
    return dff['modified_date'].max().year


@app.callback(
    Output('year--slider', 'min'),
    [Input('species-column', 'value'),
     Input('isolation-source-column', 'value'),
     Input('radio-items', 'value')])
def update_slider_min(tax_id, iso_values, radio):
    '''
    reset min year to avoid None errors when drawing figure
    '''
    dff = df[df['species'] == tax_id]
    if radio == 'is_type':
        dff = dff[dff['is_type']]
    elif radio == 'is_out':
        dff = dff[dff['is_out']]
    if iso_values:
        dff = dff[dff['isolation_source'].isin(iso_values)]
    return dff['modified_date'].min().year


@app.callback(
    Output('year--slider', 'max'),
    [Input('species-column', 'value'),
     Input('isolation-source-column', 'value'),
     Input('radio-items', 'value')])
def update_slider_max(tax_id, iso_values, radio):
    '''
    reset max year
    '''
    dff = df[df['species'] == tax_id]
    if radio == 'is_type':
        dff = dff[dff['is_type']]
    elif radio == 'is_out':
        dff = dff[dff['is_out']]
    if iso_values:
        dff = dff[dff['isolation_source'].isin(iso_values)]
    return dff['modified_date'].max().year


@app.callback(
    Output('year--slider', 'marks'),
    [Input('species-column', 'value'),
     Input('isolation-source-column', 'value'),
     Input('radio-items', 'value')])
def update_slider_marks(tax_id, iso_values, radio):
    '''
    reset marks
    '''
    dff = df[df['species'] == tax_id]
    if radio == 'is_type':
        dff = dff[dff['is_type']]
    elif radio == 'is_out':
        dff = dff[dff['is_out']]
    if iso_values:
        dff = dff[dff['isolation_source'].isin(iso_values)]
    modified_dates = dff['modified_date'].apply(lambda x: x.year).unique()
    return {str(year): str(year) for year in modified_dates}


@app.callback(
    Output('isolation-source-column', 'options'),
    [Input('species-column', 'value'),
     Input('radio-items', 'value')])
def update_isolation_source_options(tax_id, radio):
    '''
    '''
    dff = df[df['species'] == tax_id]
    if radio == 'is_type':
        dff = dff[dff['is_type']]
    elif radio == 'is_out':
        dff = dff[dff['is_out']]
    iso = dff['isolation_source']
    options = []
    for k, v in iso.value_counts().items():
        options.append({'label': '{} ({})'.format(k, v), 'value': k})
    return options


@app.callback(
    Output('isolation-source-column', 'value'),
    [Input('species-column', 'value'),
     Input('radio-items', 'value')])
def update_isolation_source_value(tax_id, radio):
    return None  # clear values


@app.callback(
    Output('plot', 'figure'),
    [Input('species-column', 'value'),
     Input('xaxis-column', 'value'),
     Input('yaxis-column', 'value'),
     Input('year--slider', 'value'),
     Input('isolation-source-column', 'value'),
     Input('radio-items', 'value'),
     Input('submit-button', 'n_clicks')],
    [State('state', 'hidden'),
     State('text-input', 'value'),
     State('url', 'search')])
def update_graph(tax_id, xaxis_column_name, yaxis_column_name,
                 year_value, iso_values, radio, n_clicks,
                 state, text, search):
    # pprint.pprint(locals())
    dff = df[df['species'] == tax_id]
    if radio == 'is_type':
        dff = dff[dff['is_type']]
    elif radio == 'is_out':
        dff = dff[dff['is_out']]

    dff = dff[dff['modified_date'] <= str(year_value+1)]
    if iso_values:
        dff = dff[dff['isolation_source'].isin(iso_values)]

    # decide if we should allow plot to calculate axes ranges
    if state is not None and tax_id == state['tax_id']:
        x_range = state['xrange']
        y_range = state['yrange']
    else:
        x_range = None
        y_range = None

    dff['text'] = dff.apply(
        lambda x: 'seqname: {seqname}<br>'
                  'accession: {version}<br>'
                  'modified_date: {modified_date}<br>'
                  'isolation source: {isolation_source}'.format(**x),
        axis='columns')

    # decide selected points
    dff['selected'] = False
    if state is None:
        args = urllib.parse.parse_qs(urllib.parse.urlparse(search).query)
        if 'seqname' in args:
            seqname = args['seqname'][0]
            dff.loc[dff['seqname'] == seqname, 'selected'] = True
        elif 'version' in args:
            ver = args['version'][0]
            dff.loc[dff['version'] == ver, 'selected'] = True
        elif 'accession' in args:
            acc = args['accession'][0]
            dff.loc[dff['accession'] == acc, 'selected'] = True
        else:
            dff.loc[dff['is_out'] & ~dff['is_type'], 'selected'] = True
    elif state['n_clicks'] < n_clicks:  # button clicked
        if text in df['accession'].values:
            dff.loc[dff['accession'] == text, 'selected'] = True
        elif text in df['version'].values:
            dff.loc[dff['version'] == text, 'selected'] = True
        elif text in df['seqname'].values:
            dff.loc[dff['seqname'] == text, 'selected'] = True
        else:
            dff.loc[dff['is_out'] & ~dff['is_type'], 'selected'] = True
    else:
        dff.loc[dff['is_out'] & ~dff['is_type'], 'selected'] = True

    inliers = dff[~dff['is_out'] & ~dff['is_type']].copy()
    inliers['iselected'] = range(len(inliers))
    outliers = dff[dff['is_out'] & ~dff['is_type']].copy()
    outliers['iselected'] = range(len(outliers))
    types = dff[dff['is_type']].copy()
    types['iselected'] = range(len(types))

    figure = {
        'data': [
            go.Scatter(
                customdata=inliers.index,
                marker={
                    'size': 15,
                    'color': 'mediumblue',
                    'opacity': 0.5,
                    'line': {'width': 0.5, 'color': 'white'}},
                mode='markers',
                name='inliers',
                selected=go.scatter.Selected(marker={'size': 15}),
                selectedpoints=inliers[inliers['selected']]['iselected'],
                unselected=go.scatter.Unselected(marker={'size': 5}),
                text=inliers['text'],
                x=inliers[xaxis_column_name],
                y=inliers[yaxis_column_name]),
            go.Scatter(
                customdata=outliers.index,
                marker={
                    'size': 15,
                    'color': 'lightgreen',
                    'opacity': 0.5,
                    'line': {'width': 0.5, 'color': 'white'}},
                mode='markers',
                name='outliers',
                text=outliers['text'],
                selected=go.scatter.Selected(marker={'size': 15}),
                selectedpoints=outliers[outliers['selected']]['iselected'],
                unselected=go.scatter.Unselected(marker={'size': 5}),
                x=outliers[xaxis_column_name],
                y=outliers[yaxis_column_name]),
            go.Scatter(
                customdata=types.index,
                marker={
                    'size': 15,
                    'color': 'red',
                    'opacity': 0.5,
                    'symbol': 'triangle-up',
                    'line': {'width': 0.5, 'color': 'white'}},
                mode='markers',
                name='types',
                selected=go.scatter.Selected(marker={'size': 15}),
                selectedpoints=types[types['selected']]['iselected'],
                unselected=go.scatter.Unselected(marker={'size': 5}),
                text=types['text'],
                x=types[xaxis_column_name],
                y=types[yaxis_column_name]),
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
            showlegend=False,
            legend={'orientation': 'v'},
            hovermode='closest',
            height=750,
            width=750,
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
     Input('plot', 'figure')])
def update_state(n_clicks, tax_id, figure):
    '''
    preserve state on the client to help determine how actions are processed
    here on the server
    '''
    return {
        'n_clicks': n_clicks,  # determine if button was clicked
        'tax_id': tax_id,  # if tax_id changes we reset everything
        'xrange': figure['layout']['xaxis']['range'],  # preserve axes ranges
        'yrange': figure['layout']['yaxis']['range']}


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
    if state is None:
        args = urllib.parse.parse_qs(urllib.parse.urlparse(search).query)
        if 'seqname' in args:
            rows = dff[dff['seqname'] == args['seqname'][0]]
        elif 'version' in args:
            rows = dff[dff['version'] == args['version'][0]]
        elif 'accession' in args:
            rows = dff[dff['accession'] == args['accession'][0]]
        else:
            rows = dff[dff['is_out']]
    elif state['n_clicks'] < n_clicks:  # button was pressed
        rows = dff[(dff['seqname'] == text) |
                   (dff['accession'] == text) |
                   (dff['version'] == text)]
    elif tax_id != state['tax_id']:
        rows = dff[dff['is_out']]
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
        'padding': '10',
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
                    children=r[c])
            elif c == 'match_species':
                cell = html.A(
                    href='https://www.ncbi.nlm.nih.gov/nuccore/' + r[c],
                    children=r[c])
            else:
                cell = r[c]
            tds.append(html.Td(children=[cell], style=TABLE_STYLE))
        trs.append(html.Tr(children=tds))
    return [html.Table(children=trs)]


if __name__ == '__main__':
    app.run_server(port=8005)
