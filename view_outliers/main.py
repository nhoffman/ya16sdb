from itertools import zip_longest
from collections import Counter
from functools import partial
from itertools import groupby
import pprint

import pandas as pd

from bokeh.layouts import row, widgetbox, column
from bokeh.models import Select
from bokeh.palettes import Spectral5
from bokeh.plotting import curdoc, figure
from bokeh.models.widgets import Div, TextInput, Button, CheckboxButtonGroup

COLORS = Spectral5
discrete = ['cluster', 'is_out', 'is_type', 'match_species']
continuous = ['x', 'y', 'dist', 'ord', 'match_pct']
NO_POINTS_SELECTED = '(no points selected)'

# TODO: provide input file as an environment variable or command line
# argument (see 'bokeh server --args')
fname = 'output/filtered_details.feather'
DEFAULT_SPECIES = 'Enterococcus faecalis'
DEFAULT_GENUS = DEFAULT_SPECIES.split()[0]

global data
data = pd.read_feather(fname)

# TODO: move to prepare.py?
for colname in discrete + ['tax_name']:
    data[colname] = data[colname].fillna('<no value>')
    if colname.startswith('is_'):
        data[colname] = data[colname].astype(str)


def ifelse(items, if_t, if_f):
    return [if_t if x else if_f for x in items]


def get_species_data(tax_name=None):
    # globals: data (data_frame), organism (a selector widget)

    tax_name = tax_name or select_species.value
    tab = data.loc[data['tax_name'] == tax_name]
    tab = tab.sort_values(by=['dist']).reindex()
    tab['ord'] = pd.Series(range(len(tab.dist)), index=tab.index)
    return tab


def get_colors(values):
    counts = Counter(values).most_common()
    if len(counts) < len(COLORS):
        zfun = zip
    else:
        zfun = partial(zip_longest, fillvalue=COLORS[-1])

    d = {
        val: color
        for color, (val, count) in zfun(
            COLORS, Counter(values).most_common())
    }
    return [d[val] for val in values]


def create_figure(df=None):

    df = df if df is not None else get_species_data()

    xs = df[x.value].values
    ys = df[y.value].values
    x_title = x.value.title()
    y_title = y.value.title()

    kw = dict()
    if x.value in discrete:
        kw['x_range'] = sorted(set(xs))
    if y.value in discrete:
        kw['y_range'] = sorted(set(ys))

    n_out = sum(df['is_out'] == 'True')
    n_type = sum(df['is_type'] == 'True')

    plot_title = ('{tax_name}, tax_id {tax_id} '
                  'outliers: {n_out} ({pct_out:.01f}%) '
                  'type: {n_type} '
                  'total: {total}').format(
                      tax_name=df['tax_name'].iloc[0],
                      tax_id=df['species'].iloc[0],
                      n_out=n_out,
                      pct_out=100 * n_out / len(df.index),
                      n_type=n_type,
                      total=len(df.index),
                  )

    kw['title'] = plot_title

    p = figure(plot_height=600, plot_width=800,
               tools='pan,box_zoom,tap,reset,save',
               **kw)

    p.xaxis.axis_label = x_title
    p.yaxis.axis_label = y_title

    if x.value in discrete:
        p.xaxis.major_label_orientation = pd.np.pi / 4

    alpha_on, alpha_off = 0.6, 0.1
    sz_on, sz_off = 9, 1

    try:
        # global variable grid may not be defined yet
        selectors = grid.children[0].children[1].children
        selections = [(s.title, s.value) for s in selectors if s.value != 'All']
    except NameError:
        selections = None

    if selections:
        selected = list(pd.DataFrame(
            [df[col] == val for col, val in selections]).apply(all))
        print(selections)
        print('{} points selected'.format(sum(selected)))

        alpha = ifelse(selected, alpha_on, alpha_off)
        sz = ifelse(selected, sz_on, sz_off)
    else:
        alpha = alpha_on
        sz = sz_on

    c = "#31AADE"
    if color.value != 'None':
        c = get_colors(df[color.value].values)

    points = p.scatter(x=xs, y=ys,
                       marker='circle',
                       color=c,
                       size=sz,
                       alpha=alpha,
                       hover_color='white',
                       hover_alpha=0.5,
                       nonselection_alpha=0.3)

    points.data_source.on_change('selected', on_selection_change)

    return p


# updater callbacks refer to global variables!

def on_selection_change(attr, old, new):
    # 'new' contains indices into data source, eg
    # {'0d': {'get_view': {}, 'glyph': None, 'indices': []},
    #  '1d': {'indices': [6946]},
    #  '2d': {'indices': {}}}

    try:
        idx = new['1d']['indices'][0]
    except IndexError:
        select_msg.text = NO_POINTS_SELECTED
    else:
        df = get_species_data()
        d = df.iloc[idx].to_dict()
        pprint.pprint(d)

        url = 'https://www.ncbi.nlm.nih.gov/nuccore'

        msg = ('<a href="{url}/{accession}" style="color: white;" target=_blank>'
               '{accession}</a> {description} '
               '[{match_pct}% match to {match_species} {match_seqname}]').format(
                   url=url, **d)

        select_msg.text = msg


def update(attr, old, new):
    # see definition of grid object below to explain indexing
    grid.children[1].children[0] = create_figure()


def update_genus(attr, old, new):
    # reset list of species for the newly selected genus (second
    # element of the first widgetbox)
    grid.children[0].children[0].children[1].options = genera[new]

    # update species to first option
    select_species.value = genera[new][0]
    update_species(attr, '', genera[new][0])


def update_species(attr, old, new):
    """Get a df of species data based on current value of
    `select_species.value`

    """

    df = get_species_data()
    select_msg.text = NO_POINTS_SELECTED

    # clear selectors but not colors and axes
    grid.children[0].children[1] = widgetbox(get_selectors(df), width=200)
    grid.children[1].children[0] = create_figure(df)


def reset_species():

    # clear all selectors
    x.value = 'x'
    y.value = 'y'
    color.value = 'is_out'
    select_accession.value = None
    update_species(None, None, None)


def search_accession(attr, old, new):
    accession = new
    result = data.loc[data['accession'] == accession]
    if result.shape[0]:
        tax_name = list(result['tax_name'])[0]
        df = get_species_data(tax_name)

        select_genus.value = tax_name.split()[0]
        select_species.value = tax_name
        grid.children[1].children[0] = create_figure(
            df, select_col='accession', select_val=accession)

        # clear selectors
        grid.children[0].children[1] = widgetbox(get_selectors(df), width=200)
    else:
        select_msg.text = 'no match for accession "{}"'.format(accession)


def get_selectors(df, colnames=discrete):
    selectors = []
    for col in colnames:
        control = Select(
            title=col, value="All",
            options=['All'] + [val[0] for val in Counter(df[col]).most_common()])
        control.on_change('value', update)
        selectors.append(control)

    return selectors


x = Select(title='X-Axis', value='x', options=continuous + discrete)
x.on_change('value', update)

y = Select(title='Y-Axis', value='y', options=continuous + discrete)
y.on_change('value', update)

color = Select(title='Color', value='is_out', options=['None'] + discrete)
color.on_change('value', update)

select_msg = Div(text=NO_POINTS_SELECTED, width=900)

df = get_species_data(DEFAULT_SPECIES)

# get dict of genus and species names where there is plot data
select_accession = TextInput(title='Accession search')
select_accession.on_change('value', search_accession)

species_names = sorted(set(data['tax_name'].loc[data['x'].notna()]))
genera = {genus: list(species)
          for genus, species in groupby(species_names, lambda x: x.split()[0])}

select_genus = Select(
    title='Genus', value=DEFAULT_GENUS, options=sorted(genera.keys()))
select_genus.on_change('value', update_genus)

select_species = Select(
    title='Species', value=DEFAULT_SPECIES, options=genera[DEFAULT_GENUS])
select_species.on_change('value', update_species)

reset = Button(label='reset species', button_type='primary')
reset.on_click(reset_species)

# selectors
leftcol = column(
    widgetbox([select_genus, select_species, reset,
               select_accession, x, y, color], width=200),
    widgetbox(get_selectors(df), width=200),
)

rightcol = column(create_figure(df), select_msg)

# 'grid' is a global variable with references in update()
grid = row(leftcol, rightcol)

curdoc().add_root(grid)
curdoc().title = "Outliers"
