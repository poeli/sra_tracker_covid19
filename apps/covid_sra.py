import plotly.express as px
import plotly.graph_objs as go
import dash_bio as dashbio
import dash_bootstrap_components as dbc
import dash
from dash import Dash, dcc, html, State, Input, Output, dash_table
import pandas as pd
from datetime import datetime as dt
import numpy as np
import logging

from app import app

logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M',
)

## default file path

display_fields_dict = {
    'COUNTRY': 'Country',
    'SEQ_TYPE': 'Type', 
    'PRIMER': 'Protocol', 
    'SEQ_TECH': 'Platform', 
    'GENDER': 'Gender', 
    'VARIANT': 'VOC/VOI', 
    'ORG': 'Center', 
}

sra_agg_files = [
    'data/sars-cov-2.srr_seqmeta.tsv',
    'data/sars-cov-2.srr_age.tsv',
    'data/sars-cov-2.srr_date.tsv',
    'data/sars-cov-2.srr_gender.tsv',
    'data/sars-cov-2.srr_gisaid.tsv',
    'data/sars-cov-2.srr_loc.tsv',
    'data/sars-cov-2.srr_name.tsv',
    'data/sars-cov-2.srr_org.tsv',
    'data/sars-cov-2.srr_lineage.tsv',
    'data/sars-cov-2.srr_lat_lon.tsv',
    'data/sars-cov-2.srr_records.tsv',
]

file_alnstats = 'data/sars-cov-2.srr_ncbi_alnstats.tsv'
file_snp_pkl = 'data/SRR_EC19_project_df_snps.pkl'

## init data

df_sra = pd.DataFrame()

for filename in sra_agg_files:
    df = pd.read_csv(filename, sep='\t')

    if len(df_sra)>0:
        df_sra = df_sra.merge(df, how='left', on='ACC')
    else:
        df_sra = df
        
df_sra.loc[~df_sra.BIOPROJECT.str.startswith("PRJ"), 'BIOPROJECT'] = 'N/A'

for col in display_fields_dict:
    df_sra[col] = df_sra[col].fillna('N/A')

df_sra['COL_DATE'] = pd.to_datetime(df_sra['COL_DATE'], errors = 'coerce')    

df_sra_filtered = df_sra

# reading alnstats
df_ec19_alnstats = pd.read_csv(file_alnstats, sep="\t")
df_ec19_alnstats = df_ec19_alnstats.rename(columns={"Ref":"Sample", "Ref_GC%":"Ref_GC_pct",  "Ref_recovery%": "Ref_recovery_pct",  "Avg_fold(x)": "Avg_fold_x"})
df_ec19_alnstats = df_ec19_alnstats.drop(df_ec19_alnstats[df_ec19_alnstats.Sample=="Ref"].index) # removing additional header lines

df_ec19_alnstats = df_ec19_alnstats.astype({'Ref_GC_pct': 'float', 
                                            'Mapped_reads': 'int64', 
                                            'Ref_recovery_pct': 'float', 
                                            'Avg_fold_x': 'float',
                                            'Fold_std': 'float', 
                                            'Num_of_Gap': 'int64',
                                            'Total_Gap_bases': 'int64'})

# loading snps
df_snps = pd.read_pickle(file_snp_pkl)
df_snps.describe(include='all')

# prepare SRA table
sra_table_columns=[{'name': 'Run',        'id': 'ACC'},
                   {'name': 'Experiment', 'id': 'EXP'},
                   {'name': 'BioProject', 'id': 'BIOPROJECT'},
                   {'name': 'Type',       'id': 'SEQ_TYPE'},
                   {'name': 'Seq Tech',   'id': 'SEQ_TECH'},
                   {'name': 'Name',       'id': 'NAME'},
                   {'name': 'Collection', 'id': 'COL_DATE'},
                   {'name': 'Country',    'id': 'COUNTRY'},
                   {'name': 'Lineage',    'id': 'LINEAGE'},
                   {'name': 'Variant',    'id': 'VARIANT'},
                   {'name': 'Sequencing center', 'id': 'ORG'}]

sra_table = dash_table.DataTable(
                id='sra-table',
                data=df_sra_filtered.head(200).to_dict('records'),
                sort_action="native",
                columns=sra_table_columns,
                page_size=20,
                row_selectable='single',
                filter_action='native',
                style_as_list_view=True,
                style_cell={'fontSize':12, 'font-family':'Helvetica, Arial, sans-serif'},
                style_cell_conditional=[
                    {
                        'if': {'column_id': 'ORG'},
                        'textAlign': 'left'
                    }
                ],
                export_format='xlsx',
                export_headers='display',
)

# functions

def filter_data(df_sra, data):
    """
    Filtering SRA dataframe with input options

    :param df_sra: dataframe
    :param data: {}
    :return df: dataframe
    """
    global df_snps
    global df_ec19_alnstats
    
    idx = ~df_snps.SAMPLE.isna()
        
    if len(data['acc_selected']):
        search_idx = df_sra.ACC.isin(data['acc_selected']) | \
                        df_sra.STUDY.isin(data['acc_selected']) | \
                        df_sra.EXP.isin(data['acc_selected']) | \
                        df_sra.BIOPROJECT.isin(data['acc_selected']) | \
                        df_sra.GISAID.isin(data['acc_selected'])
        idx = idx & search_idx
        
    if len(data['seq_types']):
        idx = idx & df_sra.SEQ_TYPE.isin(data['seq_types'])

    if len(data['primer_types']):
        idx = idx & df_sra.PRIMER.isin(data['primer_types'])

    if len(data['tech_types']):
        idx = idx & df_sra.SEQ_TECH.isin(data['tech_types'])        
        
    if len(data['lineages']):
        idx = idx & df_sra.LINEAGE.isin(data['lineages'])

    if len(data['variants']):
        idx = idx & df_sra.VARIANT.isin(data['variants'])

    if len(data['mutations']):
        mut_idx = df_snps.Mutation.isin(data['mutations'])
        acc = df_snps[mut_idx].SAMPLE.unique()
        if len(acc):
            idx = idx & df_sra.ACC.isin(acc)

    if data['min_ref_cov'] > 0:
        idx = idx & (df_ec19_alnstats.Ref_recovery_pct >= data['min_ref_cov'])

    if data['min_depth'] > 0:
        idx = idx & (df_ec19_alnstats.Avg_fold_x >= data['min_depth'])

    return df_sra[idx]


def field_color_mapping(df, coloring_field):
    """
    Generating color mapping dict

    :param df: dataframe
    :param coloring_field: string
    :return dict:
    """
    if not coloring_field:
        coloring_field='COUNTRY'
    
    # add colors
    c = list(df[coloring_field].value_counts().sort_values(ascending=False).keys())

    colors = ['#AAAAAA'] + px.colors.qualitative.Dark24 + px.colors.qualitative.Light24
    
    if 'N/A' in c:
        c.remove('N/A')
        c = ['N/A']+c
    else:
        colors.remove('#AAAAAA')
    
    c_dict = dict(zip(c, colors[:len(c)]))
    
    return {'field': coloring_field, 'color_map': c_dict}


def generate_fig_metadata(df, color_data):
    """
    Generate figure for SRA metadata

    :param df: pd.Dataframe
    :param color_data: dict
    :return figure:
    """
    global display_fields_dict
    
    coloring_field = color_data['field']
    color_map = color_data['color_map']
    
    # get displaying fields and move the coloring field to the front
    display_fields = list(display_fields_dict.keys())
    display_fields.insert(0, display_fields.pop(display_fields.index(coloring_field)))

    # init df
    df = df[display_fields]

    df['COLORS'] = df[coloring_field].map(color_map)
    
    fig_metadata = px.parallel_categories(df, 
                                          dimensions=display_fields, 
                                          labels=display_fields_dict,
                                          color='COLORS')

    fig_metadata.update_traces(
        hoverinfo='count+probability',
    )
    
    fig_metadata.update_layout(
        margin=dict(l=20, r=60, t=20, b=20),
        height=800,
    )
    
    return fig_metadata


def generate_fig_sample_time(df, color_data):
    """
    Generate barplot over time scale for SRA records
    """

    coloring_field = color_data['field']
    color_map = color_data['color_map']
    
    df = df.groupby(['COL_DATE', coloring_field]).count()['ACC'].reset_index()
    
    fig = px.bar(
        df, 
        x='COL_DATE', 
        y='ACC',
        color=coloring_field,
        color_discrete_map=color_map,
        template='simple_white',
        log_y=True,
        custom_data=[coloring_field],
        labels={
             "COL_DATE": "Collection date",
             "ACC": "Num of runs"
         },
    )
    
    fig.update_traces(
         hovertemplate="Collection date: %{x}<br>" +
                       "Number of runs: %{y}<br><br>" +
                       "<i><b>[Click to display SRA list]</b></i>"
    )
    
    # Add range slider
    fig.update_layout(
        xaxis=dict(
            rangeselector=dict(
                buttons=list([
                    dict(count=1,
                         label="1m",
                         step="month",
                         stepmode="backward"),
                    dict(count=6,
                         label="6m",
                         step="month",
                         stepmode="backward"),
                    dict(count=1,
                         label="YTD",
                         step="year",
                         stepmode="todate"),
                    dict(count=1,
                         label="1y",
                         step="year",
                         stepmode="backward"),
                    dict(step="all")
                ])
            ),
            rangeslider=dict(
                visible=True
            ),
            type="date"
        )
    )

    fig.update_layout(
        font_size=11,
        title_font_family="Helvetica, Arial, sans-serif",
        height=500,
        legend=dict(
          orientation="h",
          xanchor="left",
          yanchor="bottom",
          x=0,
          y=-0.7,
        ),
        # plot_bgcolor='rgb(233,233,233, 0.1)',
    )

    return fig


def generate_snp_plot(df, color_data):
    """
    Generate SNP states
    """
    coloring_field = color_data['field']
    color_map = color_data['color_map']

    fig = px.bar(
        df, 
        x='SNP_position', 
        y='SAMPLE',
        color=coloring_field,
        color_discrete_map=color_map,
        custom_data=[coloring_field],
        template='simple_white',
        log_y=True,
        labels={
             "SNP_position": "Genome position",
             "SAMPLE": "Num of runs"
         },
    )

    fig.update_traces(
        width=5,
        marker=dict(line_width=0),
        opacity=1,
        hovertemplate="Position: %{x}<br>" +
                      "Number of runs: %{y}<br><br>" +
                      "<i><b>[Click to display SRA list]</b></i>"

    )

    fig.update_layout(
        barmode='stack',
        font_size=11,
        title_font_family="Helvetica, Arial, sans-serif",
        height=400,
        legend=dict(
          orientation="h",
          xanchor="left",
          yanchor="bottom",
          x=0,
          y=-0.4,
        ),
    )

    return fig

# Layouts

# info cards

card_1 = dbc.Card([
    dbc.CardHeader("Run(s)", class_name="small"),
    dbc.CardBody(
        [
            dcc.Loading(html.H5("NA", id='info-run', className="card-title text-center"))
        ]
    ),
], color="light"
)

card_2 = dbc.Card([
    dbc.CardHeader("Study(ies)", class_name="small"),
    dbc.CardBody(
        [
            dcc.Loading(html.H5("NA", id='info-exp', className="card-title text-center"))
        ]
    ),
], color="light"
)

card_3 = dbc.Card([
    dbc.CardHeader("BioProject(s)", class_name="small"),
    dbc.CardBody(
        [
            dcc.Loading(html.H5("NA", id='info-bioproj', className="card-title text-center"))
        ]
    ),
], color="light"
)

card_4 = dbc.Card([
    dbc.CardHeader("Countries", class_name="small"),
    dbc.CardBody(
        [
            dcc.Loading(html.H5("NA", id='info-country', className="card-title text-center"))
        ]
    ),
], color="light"
)

card_5 = dbc.Card([
    dbc.CardHeader("Seq Center(s)", class_name="small"),
    dbc.CardBody(
        [
            dcc.Loading(html.H5("NA", id='info-org', className="card-title text-center"))
        ]
    ),
], color="light"
)

card_6 = dbc.Card([
    dbc.CardHeader("Avg Depth", class_name="small"),
    dbc.CardBody(
        [
            dcc.Loading(html.H5("NA", id='info-avg-depth', className="card-title text-center"))
        ]
    ),
], color="light"
)

info_cards = dbc.Row(
    [
        dbc.Col(card_1, width=6, md=4, lg=2),
        dbc.Col(card_2, width=6, md=4, lg=2),
        dbc.Col(card_3, width=6, md=4, lg=2),
        dbc.Col(card_4, width=6, md=4, lg=2),
        dbc.Col(card_5, width=6, md=4, lg=2),
        dbc.Col(card_6, width=6, md=4, lg=2),
    ],
    className="mt-3"
)

info_cards = dbc.Row(
    [
        dbc.Col(card_1, width=6, lg=2),
        dbc.Col(card_2, width=6, lg=2),
        dbc.Col(card_3, width=6, lg=2),
        dbc.Col(card_4, width=6, lg=2),
        dbc.Col(card_5, width=6, lg=2),
        dbc.Col(card_6, width=6, lg=2),
    ],
    className="mt-3"
)

header_div = html.Div(
    [
        # html.Img(
        #     src="assets/favicon.ico", 
        #     className="mr-5",
        #     style={'height':'25px', 'width':'25px', 'display':'inline-block'}
        # ),    
        html.H4(
            'SARS-CoV-2 Sequence Read Archive',
            className="mr-3",
            style={'display':'inline-block'}
        ),
        html.P(
            '',
        )
    ]
)

# filters and controllers
acc_list = df_sra.ACC.dropna().unique().tolist() + \
            df_sra.EXP.dropna().unique().tolist() + \
            df_sra.STUDY.dropna().unique().tolist() + \
            df_sra.BIOPROJECT.dropna().unique().tolist() + \
            df_sra.GISAID.dropna().unique().tolist()

acc_list.remove('N/A')

controller = html.Div(
        id="control-card",
        children=[
            html.Span("Search ACC#"),
            dcc.Dropdown(
                id="acc-select",
                options=[{"label": i, "value": i} for i in acc_list],
                value=[],
                multi=True,
            ),
            html.Br(),
            html.Span("Select sequencing type"),
            dcc.Dropdown(
                id="seq-type-select",
                options=[{"label": i, "value": i} for i in df_sra.SEQ_TYPE.dropna().unique().tolist()],
                value=df_sra.SEQ_TYPE.unique().tolist(),
                multi=True,
            ),
            html.Br(),
            html.Span("Select primer/protocol"),
            dcc.Dropdown(
                id="primer-type-select",
                options=[{"label": i, "value": i} for i in df_sra.PRIMER.dropna().unique().tolist()],
                value=df_sra.PRIMER.unique().tolist(),
                multi=True,
            ),
            html.Br(),
            html.Span("Select sequencing tech"),
            dcc.Dropdown(
                id="tech-type-select",
                options=[{"label": i, "value": i} for i in df_sra.SEQ_TECH.dropna().unique().tolist()],
                value=df_sra.SEQ_TECH.unique().tolist(),
                multi=True,
            ),
            html.Br(),
            html.Span("Select lineage", className="mt-3"),
            dcc.Dropdown(
                id="lineage-select",
                options=[{"label": i, "value": i} for i in df_sra.LINEAGE.dropna().unique().tolist()],
                value=[],
                multi=True,
            ),
            html.Br(),
            html.Span("Select variant"),
            dcc.Dropdown(
                id="variant-select",
                options=[{"label": i, "value": i} for i in df_sra.VARIANT.dropna().unique().tolist()],
                value=[],
                multi=True,
            ),
            html.Br(),
            html.Span("Select protein mutation"),
            dcc.Dropdown(
                id="mutation-select",
                options=[{"label": i, "value": i} for i in df_snps.Mutation.dropna().unique()],
                value=[],
                multi=True,
            ),
        ],
    )

# the style arguments for the sidebar.
# SIDEBAR_STYLE = {
#     'position': 'fixed',
#     'top': 0,
#     'left': 0,
#     'bottom': 0,
#     "width": "17rem",
#     "padding": "2rem 1rem",
#     "background-color": "#f8f9fa",
#     "overflow-y": "scroll",
# }
SIDEBAR_STYLE = {
    'position': 'fixed',
    'top': '4em',
    'left': 0,
    'bottom': 0,
    "width": "18rem",
    "padding": "2rem 1rem",
    "background-color": "#f8f9fa",
    "overflow-y": "scroll",
}


# the style arguments for the main content page.
CONTENT_STYLE = {
    'margin-top': '4rem',
    "margin-left": "20rem",
    "margin-right": "2rem",
    "padding": "2rem 1rem",
}

TEXT_STYLE = {
    'textAlign': 'center',
    'color': '#191970'
}

CARD_TEXT_STYLE = {
    'textAlign': 'center',
    'color': '#0074D9'
}

sidebar = html.Div(
    [
        controller
    ],
    className="small",
    style=SIDEBAR_STYLE,
)

content_row_info = info_cards

content_row_controller = dbc.Row(
    [
        dbc.Col([
                     html.Span("Coloring the figure by "),
                     dcc.Dropdown(
                         id="color-field-select",
                         options=[{"label": label, "value": col} for col, label in display_fields_dict.items()],
                         searchable=False
                     ),
                ],width=12, lg=2
        ),
        dbc.Col([
                    html.Span(
                        'Minimal genome coverage (%)'
                    ),
                    dcc.Slider(
                        id='min-ref-cov',
                        min=0,
                        max=100,
                        value=0,
                        marks={
                            0: {'label': 'No limit'},
                            50: {'label': '50%'},
                            100: {'label': '100%'}
                        },
                        tooltip={"placement": "bottom", "always_visible": True},
                    )
                ],width=12, lg=4
        ),
        dbc.Col([
                    html.Span(
                        'Minimal depth (x)',
                    ),
                    dcc.Slider(
                        id='min-depth',
                        min=0,
                        max=500,
                        value=0,
                        marks={
                            0: {'label': 'No limit'},
                            100: {'label': '100x'},
                            200: {'label': '200x'},
                            300: {'label': '300x'},
                            400: {'label': '400x'},
                            500: {'label': '500x'}
                        },
                        tooltip={"placement": "bottom", "always_visible": True},
                    )
                ],width=12, lg=4
        ),
        dbc.Col([
                    dbc.Button(
                          "Display SRA list",
                          id="display-sra-btn",
                          color="primary",
                          class_name="btn-sm",
                    )
                ],width=12, lg=2, style={"text-align": "end"}
        ),
    ],
    class_name="mt-5 small",
    align="center",
)

content_row_metadata_fig = dbc.Row(
    [
        dbc.Col(
            dcc.Graph(id="fig-metadata"), md=12
        ),
    ]
)

content_row_time_fig_title = dbc.Row(
    [
        dbc.Col([
                  html.H5("Sample collection date"), 
                  html.P(id='num-sra-time', className='small'),
                ], 
                width=6
        ),
        dbc.Col([
                  dbc.Button(
                      "Display SRA in range",
                      id="display-sra-range-btn",
                      color="primary",
                      class_name="btn-sm",
                  )
                ], width=6, style={"text-align": "right"}
        ),
    ], class_name='mt-3'
)


content_row_time_fig = dbc.Row(
    [
        dbc.Col(
            dcc.Graph(id='fig-sra-time'), md=12
        ),
    ]
)


content_snp_fig_title = dbc.Row(
    [
        dbc.Col([
                  html.H5("SNP/InDels"), 
                  html.P(id='snp-note', className='small'),
                ], 
                width=6
        ),
    ], class_name='mt-3'
)


content_snp_fig = dbc.Row(
    [
        dbc.Col(
            dcc.Graph(id='fig-sra-snp'), md=12
        ),
    ]
)


content_row_modal = dbc.Row([
    html.Div([
        dbc.Modal(
            [
                dbc.ModalHeader(dbc.ModalTitle("SRA/ERA/DRA")),
                dbc.ModalBody([
                    dbc.Alert(
                        "EDGE COVID-19 data for selected SRA is unavailable.",
                        id="selected-sra-alert",
                        dismissable=True,
                        fade=True,
                        is_open=False,
                        duration=4000, 
                        color="danger",
                        className="w-100 small",
                    ),
                    dcc.Loading(sra_table),
                ]),
            ],
            id="modal-body-scroll",
            centered=True,
            scrollable=True,
            is_open=False,
            className='modal-xl',
        ),
    ])
])


content = html.Div(
    [
        # header_div,
        content_row_info,
        content_row_controller,
        content_row_metadata_fig,
        content_row_time_fig_title,
        content_row_time_fig,
        content_snp_fig_title,
        content_snp_fig,
        content_row_modal
    ],
    style=CONTENT_STYLE
)


# DASH app

# app = JupyterDash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
# app.title = "COVID19 SRA"
# app._favicon = "favicon.ico"

#app.layout = html.Div(
layout = html.Div(
    [
        dcc.Store(id='aggregate_data'),
        dcc.Store(id='color_mapping'),
        dcc.Store(id='clientside_callback_temp'),
        sidebar,
        content
    ]
)

# controllers -> aggregate_data.data
@app.callback(Output('aggregate_data', 'data'),
              [
                Input('acc-select', 'value'),
                Input('seq-type-select', 'value'),
                Input('primer-type-select', 'value'),
                Input('tech-type-select', 'value'),
                Input('mutation-select', 'value'),
                Input('lineage-select', 'value'),
                Input('variant-select', 'value'),
                Input('min-ref-cov', 'value'),
                Input('min-depth', 'value'),
                State('color-field-select', 'value'),
              ]
             )
def update_df_sra(acc_selected, seq_types, primer_types, tech_types, 
                  mutations, lineages, variants, min_ref_cov, min_depth, coloring_field):

    global df_sra_filtered
    global default_color_data
    
    data = dict(
        acc_selected = acc_selected,
        seq_types = seq_types, 
        primer_types = primer_types, 
        tech_types = tech_types,
        mutations = mutations, 
        lineages = lineages, 
        variants = variants, 
        min_ref_cov = min_ref_cov, 
        min_depth = min_depth, 
    )

    df_sra_filtered = filter_data(df_sra, data)
    default_color_data = field_color_mapping(df_sra_filtered, coloring_field)
    
    return data

    
# coloring field controllers -> color mapping
@app.callback(Output('color_mapping', 'data'),
              [
                Input('color-field-select', 'value'),
              ],
             )
def update_color_map(coloring_field):
    
    return field_color_mapping(df_sra_filtered, coloring_field)

    
# aggregate_data -> info cards
@app.callback(Output('info-run', 'children'),
              Output('info-exp', 'children'),
              Output('info-bioproj', 'children'),
              Output('info-country', 'children'),
              Output('info-org', 'children'),
              Output('info-avg-depth', 'children'),
              [
                Input('aggregate_data', 'data'),
              ]
             )
def update_info_text(data):
    global df_sra_filtered
    idx = df_ec19_alnstats.Sample.isin(df_sra_filtered.ACC)
    df_alnstats = df_ec19_alnstats[idx]

    avg_depth = np.average(df_alnstats.Avg_fold_x)
    
    return [f"{len(df_sra_filtered.ACC):,}",
            f"{len(df_sra_filtered.STUDY.unique()):,}",
            f"{len(df_sra_filtered.BIOPROJECT.unique()):,}",
            f"{len(df_sra_filtered.COUNTRY.unique()):,}",
            f"{len(df_sra_filtered.ORG.unique()):,}",
            f"{avg_depth:,.2f}x"]


# coloring field controllers -> fig_metadata
@app.callback(Output('fig-metadata', 'figure'),
              [
                Input('aggregate_data', 'data'),
                Input('color_mapping', 'data'),
              ],
              prevent_initial_call=True
             )
def update_metadata_graph(data, color_data):
    global default_color_data
    
    if not color_data: color_data = default_color_data

    return generate_fig_metadata(df_sra_filtered, color_data)


# coloring field controllers -> sra_time_graph
@app.callback(Output('fig-sra-time', 'figure'),
              [
                Input('aggregate_data', 'data'),
                Input('color_mapping', 'data'),
              ],
              prevent_initial_call=True
             )
def update_sra_time_graph(data, color_data):
    
    if not color_data: color_data = default_color_data

    return generate_fig_sample_time(df_sra_filtered, color_data)


# change range -> avail runs
@app.callback(Output('num-sra-time', 'children'),
              [
                Input('aggregate_data', 'data'),
                Input('color_mapping', 'data'),
                Input('fig-sra-time', 'relayoutData'),
              ],
              prevent_initial_call=True
             )
def update_sra_time_num(data, color_data, relayout_data):
    
    if not color_data: color_data = default_color_data
    
    coloring_field=color_data['field']

    df = df_sra_filtered[~df_sra_filtered.COL_DATE.isna()]
    start = None
    end = None
    
    if relayout_data:
        if 'xaxis.range' in relayout_data:
            (start, end) = relayout_data['xaxis.range']

        if 'xaxis.range[0]' in relayout_data:
            start = relayout_data['xaxis.range[0]']
            end = relayout_data['xaxis.range[1]']
    
    if start and end:
        df = df_sra_filtered.query(f'COL_DATE>="{start}" and COL_DATE<="{end}"')
    
    return f"Total {len(df):,} runs available in selected time period"


# hover sra time graph -> sra snp figure
@app.callback(Output('fig-sra-snp', 'figure'),
              Output('snp-note', 'children'),
              [
                Input('aggregate_data', 'data'),
                Input('color_mapping', 'data'),
                Input('fig-sra-time', 'hoverData'),
                Input('fig-sra-time', 'relayoutData'),
                State('fig-sra-time', 'figure'),
              ],
              prevent_initial_call=True
             )
def update_sra_snps(data, color_data, hover_data, relayout_data, fig_data):
    ctx = dash.callback_context

    coloring_field = color_data['field']
    target_df_sra = df_sra_filtered[~df_sra_filtered.COL_DATE.isna()]
    start = None
    end = None
    
    note = ""
    
    if "relayoutData" in ctx.triggered[0]['prop_id']:
        if 'xaxis.range' in relayout_data:
            (start, end) = relayout_data['xaxis.range']

        if 'xaxis.range[0]' in relayout_data:
            start = relayout_data['xaxis.range[0]']
            end = relayout_data['xaxis.range[1]']
            
    if "hoverData" in ctx.triggered[0]['prop_id']:
        date = hover_data['points'][0]['x']
        value = hover_data['points'][0]['customdata'][0]
        target_df_sra = target_df_sra.query(f'COL_DATE=="{date}" and {coloring_field}=="{value}"')[['ACC', coloring_field]]
        
        note = f"Samples collected on {date} / {display_fields_dict[coloring_field]} is {value}"

    if start and end:        
        # remove time from datetime
        start = start.split(' ')[0]
        end = end.split(' ')[0]

        target_df_sra = target_df_sra.query(f'COL_DATE>="{start}" and COL_DATE<="{end}"')
        note = f"Samples collected between {start} and {end}"
        
    target_df_snps = df_snps[df_snps.SAMPLE.isin(target_df_sra.ACC)]
    target_df_snps = target_df_snps.merge(target_df_sra, left_on='SAMPLE', right_on='ACC', how='left')

    df = target_df_snps.groupby(['SNP_position', coloring_field]).count()['SAMPLE'].reset_index()
    
    fig = generate_snp_plot(df, color_data)
        
    return [fig, note]


# click figure/buttons -> sra list
@app.callback(Output("modal-body-scroll", "is_open"),
              Output("sra-table", "data"),
              Output("sra-table", "selected_rows"),
              Output("sra-table", "selected_cells"),
              Output("sra-table", "active_cell"),
              [
                Input('fig-sra-time', 'clickData'),
                State('fig-sra-time', 'figure'),
                Input('fig-sra-snp', 'clickData'),
                State('snp-note', 'children'),
                Input('display-sra-btn', 'n_clicks'),
                Input('display-sra-range-btn', 'n_clicks'),
                State('color_mapping', 'data'),
                State("modal-body-scroll", "is_open"),
              ],
              prevent_initial_call=True
             )
def modal_sra_list(sra_time_click_data, figure_data, sra_snp_click_data, snp_note, sra_btn, sra_range_btn, color_map, is_open):

    ctx = dash.callback_context
    max_record = 1000
    df = df_sra_filtered
    coloring_field = color_map['field']

    triggered_id = None
    if ctx.triggered:
        triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]

    # clicking "sra list" button
    if triggered_id == 'display-sra-btn':
        # the default df
        pass
        
    # clicking "sra range" button
    elif triggered_id == 'display-sra-range-btn':
        start, end = figure_data['layout']['xaxis']['range']
        df = df.query(f'COL_DATE>="{start}" and COL_DATE<="{end}"')
    
    # clicking "sra range"
    elif triggered_id == 'fig-sra-time':
        
        date = ctx.triggered[0]['value']['points'][0]['x']
        value = ctx.triggered[0]['value']['points'][0]['customdata'][0]
        
        df = df.query(f'COL_DATE=="{date}" and {coloring_field}=="{value}"')

    # clicking "sra snp"
    elif triggered_id == 'fig-sra-snp':
        import re
        
        pos = ctx.triggered[0]['value']['points'][0]['x']
        value = ctx.triggered[0]['value']['points'][0]['customdata'][0]
        start, end, date = None, None, None
        
        idx = ~df.COL_DATE.isna()
        
        acc_list = df_snps.query(f'SNP_position=={pos}').SAMPLE.unique().tolist()
        if len(acc_list):
            idx = (idx & df.ACC.isin(acc_list))
        
        df = df[idx].query(f'{coloring_field}=="{value}"')
        
        if snp_note:
            if ' and ' in snp_note:
                start, end = re.search('(\S+) and (\S+)', snp_note).group(1,2)
                df = df.query(f'COL_DATE>="{start}" and COL_DATE<="{end}"')
            elif ' on ' in snp_note:
                date = re.search('on (\S+)', snp_note).group(1)
                df = df.query(f'COL_DATE=="{date}"')
    
    # reset selected row and cell every time the modal is opened

    return [not is_open, df.to_dict('records'), [], [], None]

@app.callback(Output("selected-sra-alert", "is_open"),
              [
                Input('sra-table', "derived_virtual_selected_rows"),
                Input('sra-table', "selected_rows"),
                Input('sra-table', "derived_viewport_selected_rows"),
                State('sra-table', "data"),
                State('sra-table', "derived_virtual_data"),
                State('sra-table', "derived_viewport_data"),
              ],
              prevent_initial_call=True
             )
def modal_sra_alert(derived_virtual_selected_rows, selected_rows, derived_viewport_selected_rows, data, derived_virtual_data, derived_viewport_data):
    from os.path import exists

    if derived_virtual_selected_rows:
        row = derived_virtual_selected_rows[0]
        acc = derived_virtual_data[row]['ACC']
        if not exists(f'assets/ec19_output/{acc}/Output/config.txt'):
            return True

    return False


# select a row on snp-table --> link to /sample?sra=acc
app.clientside_callback(
    """
    function(derived_virtual_selected_rows, selected_rows, derived_viewport_selected_rows, data, derived_virtual_data, derived_viewport_data) {
        /*
        console.log("derived_virtual_data=", derived_virtual_data)
        console.log("derived_viewport_data=", derived_viewport_data)
        console.log("derived_virtual_selected_rows=", derived_virtual_selected_rows)
        console.log("selected_rows=", selected_rows)
        console.log("derived_viewport_selected_rows=", derived_viewport_selected_rows)
        */

        row = derived_virtual_selected_rows[0]
        acc = derived_virtual_data[row]['ACC']

        var request = new XMLHttpRequest();  
        request.open('GET', `/assets/ec19_output/${acc}/Output/config.txt`, true);
        request.onreadystatechange = function(){
            if (request.readyState === 4){
                if (request.status === 200) {  
                    window.location.href = "/sample?sra="+acc
                }  
            }
        };
        request.send();    
    }
    """,
    Output('clientside_callback_temp', 'data'),
    Input('sra-table', "derived_virtual_selected_rows"),
    Input('sra-table', "selected_rows"),
    Input('sra-table', "derived_viewport_selected_rows"),
    State('sra-table', "data"),
    State('sra-table', "derived_virtual_data"),
    State('sra-table', "derived_viewport_data"),
    prevent_initial_call=True
)

if __name__ == '__main__':
    app.run_server(debug=True, port=8054)  # Turn off reloader if inside Jupyter