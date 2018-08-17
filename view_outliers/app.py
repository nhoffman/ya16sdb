#! /usr/bin/env python3
import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_table_experiments as dt
import pandas
import plotly.graph_objs as go
import pprint
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
        html.Div(
            children=[
                html.Div(
                    children=[
                        dcc.Location(id='url', refresh=False),
                        dcc.Markdown('**Search**'),
                        dcc.Input(value='', type='text', id='text-input'),
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
                    children=[dcc.Graph(id='indicator-graphic')],
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
            children=[
                dt.DataTable(
                    # Initialise the rows
                    rows=[{}],
                    editable=False,
                    row_selectable=True,
                    filterable=True,
                    resizable=True,
                    sortable=True,
                    id='datatable')],
            style={
                'display': 'inline-block',
                'padding': '10',
                'width': '97%'})])


@app.callback(
    Output('datatable', 'rows'),
    [Input('species-column', 'value'),
     Input('indicator-graphic', 'selectedData'),
     Input('indicator-graphic', 'clickData')])
def update_datatable(value, selected, clicked):
    if selected is not None:
        idx = [i['customdata'] for i in selected['points']]
        rows = df.loc[idx]
    elif clicked is not None:
        idx = [i['customdata'] for i in clicked['points']]
        rows = df.loc[idx]
    else:
        # outliers
        rows = df[(df['species'] == value) & (df['is_out'])]
    rows = rows.copy()  # avoid SettingWithCopyWarning
    rows['is_type'] = rows['is_type'].apply(lambda x: 'Yes' if x else '')
    rows['is_out'] = rows['is_out'].apply(lambda x: 'Yes' if x else '')
    rows = rows.sort_values(by='dist', ascending=False)
    rows = rows[['seqname', 'description', 'is_type', 'dist',
                 'is_out', 'match_species', 'match_pct']]
    return rows.to_dict('records')


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
    Output('genus-column', 'value'),
    [Input('url', 'search')])
def update_genus_value(search):
    args = urllib.parse.parse_qs(urllib.parse.urlparse(search).query)
    if 'species_name' in args:
        species = species_id.get(args['species_name'][0], None)
        value = species_genus.get(species, DEFAULT_GENUS)
    elif 'species_id' in args:
        value = species_genus.get(args['species_id'][0], DEFAULT_GENUS)
    else:
        value = DEFAULT_GENUS
    return value


@app.callback(
    Output('species-column', 'value'),
    [Input('species-column', 'options'),
     Input('url', 'search')])
def update_species_value(options, search):
    args = urllib.parse.parse_qs(urllib.parse.urlparse(search).query)
    if 'species_name' in args:
        value = species_id.get(args['species_name'][0], None)
    elif 'species_id' in args:
        value = args['species_id'][0]
    else:
        value = None
    if value is None or value not in set(o['value'] for o in options):
        value = options[0]['value']
    return value


@app.callback(
    Output('year--slider', 'value'),
    [Input('species-column', 'value'),
     Input('isolation-source-column', 'value')])
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
     Input('isolation-source-column', 'value')])
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
     Input('isolation-source-column', 'value')])
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
     Input('isolation-source-column', 'value')])
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
    Output('isolation-source-column', 'options'),
    [Input('species-column', 'value')])
def update_isolation_source_options(tax_id):
    iso = df[df['species'] == tax_id]['isolation_source']
    options = []
    for k, v in iso.value_counts().items():
        options.append({'label': '{} ({})'.format(k, v), 'value': k})
    return options


@app.callback(
    Output('isolation-source-column', 'value'),
    [Input('species-column', 'value')])
def update_isolation_source_value(tax_id):
    '''
    clear values
    '''
    return None


@app.callback(
    Output('species-column', 'options'),
    [Input('genus-column', 'value'),
     Input('text-input', 'value')])
def update_species_options(genus_id, text):
    group = genera.get_group(genus_id)
    group = group[['species', 'species_name']]
    return [{'label': sn, 'value': si} for si, sn in group.values]


@app.callback(
    Output('indicator-graphic', 'figure'),
    [Input('species-column', 'value'),
     Input('xaxis-column', 'value'),
     Input('yaxis-column', 'value'),
     Input('year--slider', 'value'),
     Input('text-input', 'value'),
     Input('isolation-source-column', 'value')],
    [State('indicator-graphic', 'figure')])
def update_graph(tax_id, xaxis_column_name, yaxis_column_name,
                 year_value, text_value, iso_values, figure):
    pprint.pprint(locals())
    dff = df[df['species'] == tax_id]
    dff = dff[dff['modified_date'] <= str(year_value+1)]
    if iso_values:
        dff = dff[dff['isolation_source'].isin(iso_values)]

    # use figure State to preserve max range
    if figure is not None and tax_id == figure['tax_id']:
        x_range = figure['layout']['xaxis']['range']
        y_range = figure['layout']['yaxis']['range']
    else:
        x_range = None
        y_range = None

    dff['text'] = dff.apply(
        lambda x: 'seqname: {seqname}<br>'
                  'accession: {version}<br>'
                  'modified_date: {modified_date}<br>'
                  'isolation source: {isolation_source}'.format(**x),
        axis='columns')

    inliers = dff[~dff['is_out'] & ~dff['is_type']]
    outliers = dff[dff['is_out'] & ~dff['is_type']]
    types = dff[dff['is_type']]

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
                # selected=go.scatter.Selected(marker={'size': 20}),
                unselected=go.scatter.Unselected(marker={'opacity': 0.1}),
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
                # selected=go.scatter.Selected(marker={'size': 20}),
                unselected=go.scatter.Unselected(marker={'opacity': 0.1}),
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
                name='outlier-types',
                # selected=go.scatter.Selected(marker={'size': 20}),
                unselected=go.scatter.Unselected(marker={'opacity': 0.1}),
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
        ),
        # great place to store State data
        'tax_id': tax_id,
    }
    return figure


if __name__ == '__main__':
    app.run_server(port=8005)
