import collections
import functools
import itertools
import os
import pandas as pd

from pprint import pprint

from bokeh.layouts import row, widgetbox, column
from bokeh.models import Select
from bokeh.palettes import Spectral5
from bokeh.plotting import curdoc, figure
from bokeh.models.widgets import Div, TextInput, Button

COLORS = Spectral5
discrete = ['cluster', 'is_out', 'is_type', 'match_species']
continuous = ['x', 'y', 'dist', 'ord', 'match_pct']
NO_POINTS_SELECTED = '(no points selected)'
DEFAULT_SPECIES = 'Enterococcus faecalis'
DEFAULT_GENUS = DEFAULT_SPECIES.split()[0]

TABLE = ('<a href="https://www.ncbi.nlm.nih.gov/nuccore/{accession}" '
         'style="color: white;" target=_blank>'
         '{accession}</a> {description} '
         '[{match_pct}% match to '
         '{match_species} {match_seqname}]')

# TODO: provide input file as an environment variable or command line
# argument (see 'bokeh serve --args')
fname = 'filtered_details.feather'
if 'datadir' in os.environ:
    fname = os.path.join(os.environ['datadir'], 'dedup/1200bp/named', fname)

global data
data = pd.read_feather(fname)

# TODO: move to prepare.py?
for colname in discrete + ['tax_name']:
    data[colname] = data[colname].fillna('<no value>')
    if colname.startswith('is_'):
        data[colname] = data[colname].astype(str)


def ifelse(items, if_t, if_f):
    return [if_t if x else if_f for x in items]


def get_species_data(tax_name=None, tax_id=None):
    # globals: data (data_frame), organism (a selector widget)
    if tax_id is not None:
        tab = data.loc[data['tax_id'] == tax_id]
    elif tax_name is not None:
        tab = data.loc[data['tax_name'] == tax_name]
    else:
        tab = data.loc[data['tax_name'] == select_species.value]
    tab = tab.sort_values(by=['dist']).reindex()
    tab['ord'] = pd.Series(range(len(tab.dist)), index=tab.index)
    return tab


def get_colors(values):
    counts = collections.Counter(values).most_common()
    if len(counts) < len(COLORS):
        zfun = zip
    else:
        zfun = functools.partial(itertools.zip_longest, fillvalue=COLORS[-1])

    d = {
        val: color
        for color, (val, count) in zfun(
            COLORS, collections.Counter(values).most_common())
    }
    return [d[val] for val in values]


def create_figure(df=None, select_col=None, select_val=None):
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

    p = figure(plot_height=600,
               plot_width=800,
               toolbar_sticky=False,
               tools='pan,box_zoom,box_select,tap,reset,save,hover',
               **kw)

    p.xaxis.axis_label = x_title
    p.yaxis.axis_label = y_title

    if x.value in discrete:
        p.xaxis.major_label_orientation = pd.np.pi / 4

    alpha_on, alpha_off = 0.6, 0.1
    sz_on, sz_off = 9, 1

    selections = []
    if select_col is not None and select_val is not None:
        selections.append([select_col, select_val])
    else:
        try:
            # global variable layout may not be defined yet
            selectors = layout.children[0].children[1].children
            for s in selectors:
                if s.value != 'All':
                    selections.append([s.title, s.value])
        except NameError:
            pass

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
        pprint(d)
        select_msg.text = TABLE.format(**d)


def update(attr, old, new):
    # see definition of layout object below to explain indexing
    layout.children[1].children[0] = create_figure()


def update_genus(attr, old, new):
    # reset list of species for the newly selected genus (second
    # element of the first widgetbox)
    layout.children[0].children[0].children[1].options = genera[new]

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
    layout.children[0].children[1] = widgetbox(get_selectors(df), width=200)
    layout.children[1].children[0] = create_figure(df)


def reset_species():
    # clear all selectors
    x.value = 'x'
    y.value = 'y'
    color.value = 'is_out'
    select_accession.value = None
    update_species(None, None, None)


def search_accession(attr, old, new):
    df = None
    if new is not None:
        accession = new.strip()
        results = data.loc[data['accession'] == accession]
        if results.shape[0]:
            tax_name = list(results['tax_name'])[0]
            df = get_species_data(tax_name)

            select_genus.value = tax_name.split()[0]
            select_species.value = tax_name
            layout.children[1].children[0] = create_figure(
                df, select_col='accession', select_val=accession)
            msg = []
            for _, r in results.iterrows():
                pprint(r.to_dict())
                msg.append(TABLE.format(**r.to_dict()))
            select_msg.text = '<br>'.join(msg)

            # clear selectors
            layout.children[0].children[1] = widgetbox(
                get_selectors(df), width=200)
        else:
            select_msg.text = 'no match for accession "{}"'.format(accession)
    return df


def get_selectors(df, colnames=discrete):
    selectors = []
    for col in colnames:
        control = Select(
            title=col, value="All",
            options=['All'] + [val[0] for val in
                               collections.Counter(df[col]).most_common()])
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

get = curdoc().session_context.request.arguments

if 'tax_name' in get:
    df = get_species_data(tax_name=get['tax_id'][0].decode('ascii'))
elif 'tax_id' in get:
    df = get_species_data(tax_id=get['tax_id'][0].decode('ascii'))
else:
    df = get_species_data(tax_name=DEFAULT_SPECIES)

# get dict of genus and species names where there is plot data
select_accession = TextInput(title='Accession search')
select_accession.on_change('value', search_accession)

species_names = sorted(set(data['tax_name'].loc[data['x'].notna()]))
genera = {}
for genus, species in itertools.groupby(species_names, lambda x: x.split()[0]):
    genera[genus] = list(species)

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
               select_accession, x, y, color], width=250),
    widgetbox(get_selectors(df), width=250),
)

rightcol = column(create_figure(df), select_msg)

# 'layout' is a global variable with references in update()
layout = row(leftcol, rightcol)

curdoc().add_root(layout)
curdoc().title = "Outliers"

if 'accession' in get:
    search_accession(None, None, get['accession'][0].decode('ascii'))
