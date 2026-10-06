import dash
import pandas
from dash import dash_table
from dash.dependencies import Input, Output
from functools import reduce
from operator import or_

import shared

dash.register_page(__name__, path='/datatable', name='Data Table')

TABLE_PAGE_SIZE = 50
_STR_COLS = [c for c in shared.seq_info.columns
             if pandas.api.types.is_string_dtype(shared.seq_info[c])]


def layout():
    columns = [{'name': c, 'id': c} for c in shared.seq_info.columns]
    return dash.html.Div(
        style={'width': '98%', 'margin': '10px auto'},
        children=[
            dash.html.Div(
                dash.dcc.Link('Plot', href='/',
                              style={'fontSize': '16px'}),
                style={'margin': '10px 0'}),
            dash.dcc.Input(
                id='table-search',
                type='text',
                placeholder='Search all columns...',
                debounce=True,
                style={'width': '300px', 'margin': '10px 0'}),
            dash.dcc.Loading(
                dash_table.DataTable(
                    id='seq-info-table',
                    columns=columns,
                    page_current=0,
                    page_size=TABLE_PAGE_SIZE,
                    page_action='custom',
                    sort_action='custom',
                    sort_mode='multi',
                    filter_action='custom',
                    filter_query='',
                    style_table={'overflowX': 'auto'},
                    style_cell={'textAlign': 'left', 'padding': '5px'},
                    style_header={'fontWeight': 'bold'},
                ),
                type='circle'),
        ])


def _parse_filter_query(filter_query):
    """Parse Dash DataTable filter_query into (column, value) pairs."""
    filters = []
    if not filter_query:
        return filters
    import re
    for part in filter_query.split(' && '):
        match = re.match(
            r'\{(\S+)\}\s+(?:contains|scontains)\s+"?(.*?)"?$', part.strip())
        if match:
            filters.append((match.group(1), match.group(2)))
    return filters


@dash.callback(
    Output('seq-info-table', 'data'),
    Output('seq-info-table', 'page_count'),
    Input('seq-info-table', 'page_current'),
    Input('seq-info-table', 'page_size'),
    Input('seq-info-table', 'sort_by'),
    Input('seq-info-table', 'filter_query'),
    Input('table-search', 'value'))
def update_seq_info_table(page, page_size, sort_by, filter_query, search):
    dff = shared.seq_info
    if search:
        masks = (dff[c].str.contains(search, case=False, na=False)
                 for c in _STR_COLS)
        dff = dff[reduce(or_, masks)]
    for col, val in _parse_filter_query(filter_query):
        if col in dff.columns:
            dff = dff[dff[col].astype(str).str.contains(
                val, case=False, na=False)]
    if sort_by:
        dff = dff.sort_values(
            [s['column_id'] for s in sort_by],
            ascending=[s['direction'] == 'asc' for s in sort_by])
    page_count = max(1, -(-len(dff) // page_size))
    start = page * page_size
    return (
        dff.iloc[start:start + page_size].to_dict('records'),
        page_count,
    )
