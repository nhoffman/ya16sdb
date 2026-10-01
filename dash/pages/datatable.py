import dash
from dash import dash_table
from dash.dependencies import Input, Output

import shared

dash.register_page(__name__, path='/datatable', name='Data Table')

TABLE_PAGE_SIZE = 50


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
            dash_table.DataTable(
                id='seq-info-table',
                columns=columns,
                page_current=0,
                page_size=TABLE_PAGE_SIZE,
                page_action='custom',
                sort_action='custom',
                sort_mode='multi',
                style_table={'overflowX': 'auto'},
                style_cell={'textAlign': 'left', 'padding': '5px'},
                style_header={'fontWeight': 'bold'},
            ),
        ])


@dash.callback(
    Output('seq-info-table', 'data'),
    Output('seq-info-table', 'page_count'),
    Input('seq-info-table', 'page_current'),
    Input('seq-info-table', 'page_size'),
    Input('seq-info-table', 'sort_by'),
    Input('table-search', 'value'))
def update_seq_info_table(page, page_size, sort_by, search):
    dff = shared.seq_info
    if search:
        mask = dff.astype(str).apply(
            lambda col: col.str.contains(search, case=False, na=False)
        ).any(axis=1)
        dff = dff[mask]
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
