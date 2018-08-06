#! /usr/bin/env python3
import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
import pprint
import urllib

DEFAULT_GENUS = '1350'
DEFAULT_SPECIES = '1351'

app = dash.Dash()
app.title = 'Species Outlier Plots'

df = pd.read_feather('filter_details.feather')

info = df.columns

modified_dates = df['modified_date'].apply(lambda x: x.year).unique()
tax = df[['genus', 'genus_name', 'species', 'species_name']]
tax = tax.drop_duplicates()
species_genus = dict(tax[['species', 'genus']].drop_duplicates().values)
species_id = dict(tax[['species_name', 'species']].drop_duplicates().values)
genera = tax.groupby(by='genus')

default_species_group = genera.get_group(DEFAULT_GENUS)
default_species_group = default_species_group[['species', 'species']]
species_opts = default_species_group.values
species_opts = [{'label': sn, 'value': si} for si, sn in species_opts]
genus_opts = tax[['genus', 'genus_name']].drop_duplicates().values
genus_opts = [{'label': gn, 'value': gi} for gi, gn in genus_opts]

app.layout = html.Div(children=[
    dcc.Location(id='url', refresh=False),
    dcc.Input(value='', type='text', id='text-input'),
    html.Div(children=[
        html.Div(children=[
            dcc.Dropdown(
                id='genus-column',
                options=genus_opts,
                value=DEFAULT_GENUS),
        ], style={'width': '48%', 'display': 'inline-block'}),

        html.Div(children=[
            dcc.Dropdown(
                id='species-column',
                options=species_opts,
                value=DEFAULT_SPECIES),
        ], style={'width': '48%', 'float': 'right', 'display': 'inline-block'})
    ]),
    html.Div(children=[
        html.Div(children=[
            dcc.Dropdown(
                id='xaxis-column',
                options=[{'label': i, 'value': i} for i in info],
                value='x'
            ),
        ], style={'width': '48%', 'display': 'inline-block'}),

        html.Div(children=[
            dcc.Dropdown(
                id='yaxis-column',
                options=[{'label': i, 'value': i} for i in info],
                value='y'
            ),
        ], style={'width': '48%', 'float': 'right', 'display': 'inline-block'})
    ]),
    dcc.Graph(id='indicator-graphic'),
    dcc.Slider(
        id='year--slider',
        min=modified_dates.min(),
        max=modified_dates.max(),
        value=modified_dates.max(),
        step=None,
        marks={str(year): str(year) for year in modified_dates}
    )],
    style={'width': '1200'})


@app.callback(
    dash.dependencies.Output('genus-column', 'value'),
    [dash.dependencies.Input('url', 'search')])
def update_genus_value(search):
    args = urllib.parse.parse_qs(urllib.parse.urlparse(search).query)
    if 'species_name' in args:
        species = species_id.get(args['species_name'][0], DEFAULT_SPECIES)
    elif 'species_id' in args:
        species = args['species_id']
    else:
        species = DEFAULT_SPECIES
    return species_genus[species]


@app.callback(
    dash.dependencies.Output('species-column', 'value'),
    [dash.dependencies.Input('species-column', 'options'),
     dash.dependencies.Input('url', 'search')])
def update_species_value(options, search):
    args = urllib.parse.parse_qs(urllib.parse.urlparse(search).query)
    if 'species_name' in args:
        value = species_id.get(args['species_name'][0], None)
    elif 'species_id' in args:
        value = args['species_id']
    else:
        value = None
    if value is None or value not in set(o['value'] for o in options):
        value = options[0]['value']
    return value


@app.callback(
    dash.dependencies.Output('species-column', 'options'),
    [dash.dependencies.Input('genus-column', 'value'),
     dash.dependencies.Input('text-input', 'value')])
def update_species_options(genus_id, text):
    group = genera.get_group(genus_id)
    group = group[['species', 'species_name']]
    return [{'label': sn, 'value': si} for si, sn in group.values]


@app.callback(
    dash.dependencies.Output('indicator-graphic', 'figure'),
    [dash.dependencies.Input('species-column', 'value'),
     dash.dependencies.Input('xaxis-column', 'value'),
     dash.dependencies.Input('yaxis-column', 'value'),
     dash.dependencies.Input('year--slider', 'value'),
     dash.dependencies.Input('text-input', 'value')],
    [dash.dependencies.State('indicator-graphic', 'relayoutData')])
def update_graph(species_id, xaxis_column_name, yaxis_column_name,
                 year_value, text_value, relayout_data):
    pprint.pprint(locals())
    dff = df[df['species'] == species_id]
    dff = dff[dff['modified_date'] <= str(year_value)]

    inliers = dff[~dff['is_out']]
    outliers = dff[dff['is_out']]

    dff.iloc[0]['species_name']
    figure = {
        'data': [
            go.Scatter(
                x=inliers[xaxis_column_name],
                y=inliers[yaxis_column_name],
                text='accession: ' + inliers['version'],
                mode='markers',
                marker={
                    'size': 15,
                    'color': None,
                    'opacity': 0.5,
                    'line': {'width': 0.5, 'color': 'white'}},
                name='inliers'),
            go.Scatter(
                x=outliers[xaxis_column_name],
                y=outliers[yaxis_column_name],
                text='accession: ' + outliers['version'],
                mode='markers',
                marker={
                    'size': 15,
                    'color': 'lightgreen',
                    'opacity': 0.5,
                    'line': {'width': 0.5, 'color': 'white'}},
                name='outliers')
            ],
        'layout': go.Layout(
            xaxis={
                'title': xaxis_column_name,
                'type': 'linear',
                'showgrid': False
            },
            yaxis={
                'title': yaxis_column_name,
                'type': 'linear',
                'showgrid': False
            },
            hovermode='closest',
            height=700,
            width=700,
            autosize=False,
            title='{}, tax_id {} outliers: {} ({:.0%}), total: {}'.format(
                dff.iloc[0]['species_name'],
                species_id,
                len(outliers),
                (len(outliers) / len(dff)),
                len(dff))
        )
    }
    return figure


if __name__ == '__main__':
    app.run_server(port=8005)
