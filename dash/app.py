#! /usr/bin/env python3
'''
Plotly Dash app exploring NCBI 16s records grouped by species taxonomy id
'''
import dash
import shared  # noqa: F401 - ensures global data is loaded before pages

app = dash.Dash(__name__, use_pages=True)
app.title = 'Species Outlier Plots'
app.config.suppress_callback_exceptions = True
server = app.server

app.layout = dash.html.Div([
    dash.page_container,
])


if __name__ == '__main__':
    app.run_server(debug=True)
