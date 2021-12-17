from datetime import date
import plotly.express as px
import plotly.graph_objs as go
import dash_bio as dashbio
import dash
from dash import Dash, dcc, html, Input, Output, State, dash_table, no_update
import dash_bootstrap_components as dbc
import logging
import pandas as pd
from apps.utilities import EC19_project
from apps.dashbio_pdb import SpikePdbData
import re
from glob import glob
from os.path import exists, basename

from app import app

logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M',
)

# init variables

# genes coordinates
gene = {'ORF1a'  : '266-13483',
        'ORF1ab' : '13468-21555',
        'S'      : '21563-25384',
        'ORF3a'  : '25393-26220',
        'E'      : '26245-26472',
        'M'      : '26523-27191',
        'ORF6'   : '27202-27387',
        'ORF7a'  : '27394-27759',
        'ORF7b'  : '27756-27887',
        'ORF8'   : '27894-28259',
        'N'      : '28274-29533',
        'ORF10'  : '29558-29674'}
        
# acc_num = 'SRR14135907'

# assets data url
assets_data_url = '/assets/data'

# assets data path
assets_data_path = 'assets/data'

assets_proj_url = None
proj_path = None

# for mol3d viewer
pdb         = f'{assets_data_path}/6VXX.pdb'
blast_xml   = f'{assets_data_path}/6VXX_blastp.xml'
spike_pdb   = SpikePdbData(pdb, blast_xml)
mol3d_spike = dashbio.Molecule3dViewer(id='dashbio-mol3d', modelData=spike_pdb.pdb_data)

# variant-mutation info
lineage_info_tsv = f'{assets_data_path}/lineage_info.tsv'
s_var_mut_tsv    = f'{assets_data_path}/S_protein_variant_mutation.tsv'
s_mut_category   = f'{assets_data_path}/S_mutation_category.tsv'

# filenames for igv 
bam_fn     = 'NC_045512.2.sort.bam'
bam_bai_fn = 'NC_045512.2.sort.bam.bai'
vcf_fn     = 'NC_045512.2.vcf'

# ec19 object
ec19_obj = None

# ec19 alnstats
alnstats = None
Avg_fold_x = None

# ec19 snp/snv/indels
composition_file = None # EC19 compositionlog
df_depth = None # EC19 compositionlog
df_snps = None
df_indels = None
df_depth_pos = None

# variant-mutation info
df_lineage_info = None
df_lineage_mutation = None
df_mutation_category = None
display_fields = None
df_var_mutation_plot_data = None


# functions

# snp table
def generate_snp_table():
    snp_table_columns=[{'name': 'Position', 'id': 'SNP_position'},
                       {'name': 'Ref', 'id': 'Ref_codon'},
                       {'name': 'Sub', 'id': 'Sub_codon'},
                       {'name': 'AA position', 'id': 'AA_pos'},
                       {'name': 'Mutation', 'id': 'Mutation'},
                       {'name': 'Depth', 'id': 'Count'},
                       {'name': 'Type', 'id': 'Type'},
                       {'name': 'F:R', 'id': 'F:R'}]

    snp_table = dash_table.DataTable(
                    id='snp-table',
                    sort_action="native",
                    columns=snp_table_columns,
                    data=[],
                    page_size=12,
                    style_as_list_view=True,
                    style_cell={'fontSize':12, 'font-family':'Helvetica, Arial, sans-serif'},
                    tooltip_header={
                        'SNP_position': 'Click a value to change IGV locus. Scroll down to view.',
                        'Mutation': 'Click a value to show on Mol3D',
                    },
                    tooltip_delay=1,
                    tooltip_duration=None,
                    style_data_conditional=[
                        {
                            'if': {
                                'filter_query': '{Product} = "S" && {Synonymous} = "No"',
                            },
                            'backgroundColor': '#FFD8D9',
                        },
                        {
                            'if': {
                                'column_id': 'Mutation'
                            },
                            'cursor': 'pointer',
                        },
                        {
                            'if': {
                                'column_id': 'SNP_position'
                            },
                            'cursor': 'pointer',
                        },
                    ],
                    style_header_conditional=[{
                        'if': {'column_id': col},
                        'textDecoration': 'underline',
                        'textDecorationStyle': 'dotted',
                    } for col in ['SNP_position', 'Mutation']],
    )

    return snp_table

# snp-table
snp_table = generate_snp_table()

def init_proj_data(acc_num):
    
    global ec19_obj
    global assets_proj_url
    global proj_path
    global alnstats
    global Avg_fold_x
    global df_depth
    global df_snps
    global df_indels
    global df_depth_pos
    global df_lineage_info
    global df_lineage_mutation
    global display_fields
    global df_var_mutation_plot_data
    global df_mutation_category
    
    # proj dir
    proj = f'assets/ec19_output/{acc_num}/Output'

    # assets EC19 url
    assets_proj_url = f'/assets/ec19_output/{acc_num}/Output/ReadsBasedAnalysis/readsMappingToRef/'

    # read EC19 output files
    proj_path = f'{proj}/ReadsBasedAnalysis/readsMappingToRef/'
    
    # for EC19 compositionlog
    composition_file = f'{proj_path}/NC_045512.2_consensus.compositionlog'
    
    # loading and processing EC19 project
    ec19_obj = EC19_project(proj, error_on_duplication=False, store_in_class=False)
    ec19_obj.add_project()
    
    # alnstats
    alnstats = ec19_obj.df_alnstats.query("Ref=='NC_045512_2'").reset_index().to_dict('records')[0]
    Avg_fold_x = alnstats['Avg_fold_x']

    # read EC19 compositionlog
    df_depth = pd.read_csv(composition_file, sep='\t')
    df_depth = df_depth.drop(columns=['#ID'])

    cols = ['A', 'C', 'G', 'T', 'non-ATCG']

    # cleanup 
    for col in cols:
        df_depth[col] = df_depth[col].str.replace(r' \(.*', '', regex=True)
        df_depth[col] = df_depth[col].astype(int)

    # calculate total non-ref-nuc
    df_depth['nonREF_prop'] = 0.0
    df_depth['nonREF_total'] = df_depth.apply(lambda x: x['TOTAL']-x[x['REF_NUC']]-x['non-ATCG'], axis=1)
    df_depth['nonREF'] = False
    idx = df_depth['nonREF_total']>0
    df_depth.loc[idx, 'nonREF'] = True

    # calculate non-ref-nuc proportion
    idx = df_depth['nonREF']==True
    df_depth.loc[idx, 'nonREF_prop'] = df_depth[idx].apply(lambda x: 1-x[x['REF_NUC']]/x['TOTAL'], axis=1)

    # snp/indels table
    df_snps = ec19_obj.df_snps.query("Chromosome=='NC_045512_2'")
    df_indels = ec19_obj.df_indels.query("Chromosome=='NC_045512_2'")

    df_depth_pos = df_depth.query("nonREF==True & nonREF_total>0")
    df_depth_pos = df_depth_pos.merge(df_snps[['SNP_position', 'Type']], left_on='POS', right_on='SNP_position', how='left')
    df_depth_pos = df_depth_pos.drop(columns=['SNP_position'])
    df_depth_pos.loc[df_depth_pos.Type.isna(), 'Type'] = 'SNV'

    df_snps = ec19_obj.df_snps.query("Chromosome=='NC_045512_2'")
    df_snps['id'] = df_snps['SNP_position'].astype(int)

    # variant-mutation info
    df_lineage_info = pd.read_csv(lineage_info_tsv, sep='\t')
    df_lineage_mutation = pd.read_csv(s_var_mut_tsv, sep='\t')
    df_lineage_mutation = df_lineage_mutation.merge(df_lineage_info, on='lineage')
    display_fields = ['lineage', 'mutation', 'type']
    df_lineage_mutation = df_lineage_mutation[display_fields]

    # mutation category
    df_mutation_category = pd.read_csv(s_mut_category, sep='\t')

# default mol3d stylish
def generate_mol3d_spike_default():
    
    global spike_pdb
    global df_snps
    
    mol3d_spike_pos = df_snps.query("Product=='S' and Synonymous=='No'").AA_pos.unique().tolist()
    mol3d_spike_pos = [int(float(x)) for x in mol3d_spike_pos]
    mol3d_spike_text = df_snps.query("Product=='S' and Synonymous=='No'").Mutation.unique().tolist()

    residue_indexes = spike_pdb.spike2pdb_residue_indexes(mol3d_spike_pos)
    
    # hightlight RBD
    opts = dict(highlight_bg_indexes = spike_pdb.spike2pdb_residue_indexes(range(319,541)))
    
    settings = dict(
        styles = spike_pdb.pdb_style(residue_indexes, **opts),
        labels = spike_pdb.residue_labels(mol3d_spike_pos, mol3d_spike_text, font_size=10),
        zoomTo = {"sel": {}, "animationDuration": 0, "fixedPath": False},
        zoom   = {"factor": 0.8, "animationDuration": 0, "fixedPath": False},
    )

    return settings

# prepare IGV tracks
def generate_igv_graph():

    global proj_path
    
    # check if required files exist for igv tracks
    igv_tracks = []
    
    if exists(f'{proj_path}/{vcf_fn}'):
        track = {
            'name': 'Variants',
            'type': 'variant',
            'format': 'vcf',
            'url': f'{assets_proj_url}/{vcf_fn}',
            'displayMode': 'squished',
        }
        igv_tracks.append(track)
    
    # bed paths
    bed_paths = glob(f'{proj_path}/*.bed')
    
    if len(bed_paths):
        for bed_file in bed_paths:
            # get baser name of bed file
            bed_fn = basename(bed_file)
            # convert to url
            bed = f'{assets_proj_url}/{bed_fn}'
            amplicon_track = bed_fn.replace(".bed", "")
            track = {
                'name': amplicon_track,
                'type': 'annotation',
                'format': 'bed',
                'url': bed,
                'displayMode': 'EXPANDED',
                'height': 100,
                'color': 'rgb(255,150,150)',
            }
            igv_tracks.append(track)
    
    if exists(f'{proj_path}/{bam_fn}') and exists(f'{proj_path}/{bam_bai_fn}'):
        track = {
            'name': 'Alignment',
            'type': 'alignment',
            'format': 'bam',
            'displayMode': 'SQUISHED',
            'squishedRowHeight': 10,
            'url': f'{assets_proj_url}/{bam_fn}',
            'indexURL': f'{assets_proj_url}/{bam_bai_fn}',
            'height': 500,
            # height: 200,
            # autoHeight: false
        }
        igv_tracks.append(track)

    if exists(f'{assets_data_path}/ref.gff3'):
        track = {
            'name': 'Annotations',
            'type': 'annotation',
            'url': f'{assets_data_url}/ref.gff3',
            'displayMode': 'EXPANDED',
            'filterTypes': ['CDS','mature_protein_region_of_CDS'],
            'nameField': 'gene',
            'height': 100,
            'color': 'rgb(176,141,87, 0.8)',
        }
        igv_tracks.append(track)

    igv = dashbio.Igv(
            id='reference-igv',
            genome='NC_045512_2',
            locus='NC_045512_2:1-1000',
            reference={
                'id': 'SARS-CoV-2',
                'name': 'NC_045512_2',
                'tracks': igv_tracks,
                'fastaURL': f'{assets_data_url}/ref.fa',
                'indexURL': f'{assets_data_url}/ref.fa.fai',
            },
        )
        
    return igv

# genome overview plot
def generate_overview_graph(min_nonref_pct, min_nonref_depth, pos_start, pos_end):
    
    global df_depth_pos
    global df_indels
    global Avg_fold_x
    
    df = df_depth_pos.query(f"nonREF_prop>={min_nonref_pct} & nonREF_total>={min_nonref_depth} & POS>={pos_start} & POS<={pos_end}")
    
    fig = px.scatter(df,
                     x='POS',
                     y='nonREF_total',
                     size='nonREF_prop',
                     size_max=20,
                     template='ggplot2',
                    #  title='Single-nucleotide variant overview',
                     custom_data=['REF_NUC','TOTAL','A','C','G','T','non-ATCG','nonREF','nonREF_prop','Type','POS'],
                     color='Type',
                     color_discrete_map={
                         "SNV" : "#BFBF71",
                         "Intergenic region" : "#55A05C",
                         "Synonymous" : "#B3C6E5",
                         "Non-synonymous" : "#EFA59D"
                     },
                     labels={
                         "POS": "Genome position (bp)",
                         "nonREF_total": "Depth of coverage (fold)"
                     },
    )

    fig.update_traces(
        marker=dict(
            # size=size,
            # sizemode='area',
            # sizeref=2.*max(size)/(40.**2),
            line_width=1,
            sizemin=5
        ),
        hovertemplate="""
    <b>Non-reference nucleotide (not %{customdata[0]}) %{customdata[8]:.2%}</b><br>
    <br>
    Position: %{customdata[10]}<br>
    non-Ref depth/total depth: %{y}x/%{customdata[1]}x<br>
    Reference base: %{customdata[0]}<br>
    Composition:<br>
    - A: %{customdata[2]}<br>
    - C: %{customdata[3]}<br>
    - G: %{customdata[4]}<br>
    - T: %{customdata[5]}<br>
    - Others: %{customdata[6]}<br>
    Non-reference: %{customdata[7]}<br>
    <i><b>[Click to view in IGV]</b></i>
    <extra>%{customdata[9]}</extra>"""
    )

    fig.add_trace(go.Scatter(
                    x=df_depth.query(f"POS>={pos_start} & POS<={pos_end}").POS, 
                    y=df_depth.query(f"POS>={pos_start} & POS<={pos_end}").TOTAL, 
                    line={'width': 0},
                    fillcolor='#BBBBBB', 
                    hoverinfo='skip',
                    name='Total depth',
                    fill='tozeroy'))

    fig.add_hline(y=Avg_fold_x, 
                  line_width=2, 
                  line_dash="dot",
                  annotation_text=f"Average depth {Avg_fold_x}x", 
                  annotation_position="bottom right")

    # add indels
    df_i = df_indels.query(f"INDEL_position>={pos_start} & INDEL_position<={pos_end}")
    
    fig.add_trace(go.Scatter(
                    x=df_i.INDEL_position, 
                    y=[-10]*len(df_i),
                    marker_symbol=['x']*len(df_i),
                    marker_color="#33477B",
                    marker_line_width=1,
                    marker_size=15,
                    text=df_i['Sequence'] if 'Sequence' in df_i else [],
                    hovertemplate=
                        "<b>Deletion: %{x}</b><br><br>" +
                        "Sequences: %{text}<br>" +
                        "<extra></extra>",
                    mode='markers',
                    name='Deletion'))
    
    for name in gene:
        (start, stop) = gene[name].split('-')
        label = name
        if '3' in name or '6' in name or '7' in name or '8' in name:
            label = ""
        
        fig.add_vline(x=int(start),
                      line_width=1,
                      line_dash="dot",
                      line_color="#415072",
                      annotation=dict(font_size=11),
                      annotation_text=label,
                      annotation_font_color="grey",
                      annotation_position="top right")

    fig.update_layout(
        # margin=dict(l=140, r=40, b=50, t=80),
        height=570,
        # hovermode='x unified',
        # clickmode='event+select',
        # ticks='outside',
        # xaxis= {
        #     'linecolor': 'black',
        #     'linewidth': 1,
        #     'mirror': True
        # },
        # yaxis= {
        #     'linecolor': 'black',
        #     'linewidth': 1,
        #     'mirror': True
        # },
        paper_bgcolor='white',
        plot_bgcolor='white',
    )
    
    # make data display in snp -> svn -> depth order
    data_depth = None
    data = []
    for d in fig.data:
        if d.name=='Total depth':
            data_depth = d
        else:
            data.append(d)
    
    data = [data_depth]+data
        
    fig.data = data
    fig.update_xaxes(range=[pos_start, pos_end])
    
    return fig

# variant(VOC/VOI) - mutation plot
def generate_var_mutation_plot(df):
    
    fig_var = px.scatter(df,
                         x='mutation',
                         y='lineage',
                         color='type',
                         color_discrete_map={
                             "VOC" : "#EFA59D",
                             "VOI" : "#BFBF71",
                             "Sample" : "#2C70F4",
                         },
                         height=500, 
                         size=[10]*len(df),
                         size_max=8,
                         template='ggplot2',
                         marginal_x="histogram", 
                         marginal_y="histogram",
                         symbol='category',
                         symbol_sequence=['circle', 'square', 'diamond','circle-x','square-x', 'diamond-x'],
                         labels={
                             "lineage": "Variant",
                             "mutation": "Protein mutation",
                             "type": "Variant type"
                         },
    )
    
    fig_var.update_traces(hovertemplate="%{x}: %{y}<br><i><b>[Click to view in Mol3D]</b></i>")
    fig_var.update_layout(barmode='stack',
                          hovermode = 'closest',
                          xaxis={
                              'categoryorder': 'array', 
                              'categoryarray': df.sort_values('pos').mutation.unique(),
                              'domain': [0, 0.9]
                          },
                          xaxis2={'domain': [0.9, 1]},
                          xaxis3={'domain': [0, 0.9]},
                          font_size=11,
                          legend=dict(
                              orientation="h",
                              xanchor="left",
                              yanchor="bottom",
                              x=0,
                              y=1,
                          ),
                          plot_bgcolor='rgb(233,233,233, 0.1)'
    )

    return fig_var

# generate info cards
def generate_layout_info_cards():
    card_1 = dbc.Card([
        dbc.CardHeader("Mapped reads", class_name="small"),
        dbc.CardBody(
            [
                html.H5(f"{alnstats['Mapped_reads']:,}", id='info-mapped-reads', className="card-title text-center"),
            ]
        ),
    ], color="light"
    )

    card_2 = dbc.Card([
        dbc.CardHeader("% Mapped reads", class_name="small"),
        dbc.CardBody(
            [
                html.H5(f"{alnstats['Mapped_reads_pct']}%", id='info-mapped-reads-pct', className="card-title text-center"),
            ]
        ),
    ], color="light"
    )

    card_3 = dbc.Card([
        dbc.CardHeader("% coverage", class_name="small"),
        dbc.CardBody(
            [
                html.H5(f"{alnstats['Ref_recovery_pct']:.2f}%", id='info-genome-coverage', className="card-title text-center"),
            ]
        ),
    ], color="light"
    )

    card_4 = dbc.Card([
        dbc.CardHeader("Avg depth", class_name="small"),
        dbc.CardBody(
            [
                html.H5(f"{alnstats['Avg_fold_x']}x", id='info-avg-depth', className="card-title text-center"),
            ]
        ),
    ], color="light"
    )

    card_5 = dbc.Card([
        dbc.CardHeader("SNPs", class_name="small"),
        dbc.CardBody(
            [
                html.H5(alnstats['Num_of_SNPs'], id='info-num-snps', className="card-title text-center"),
            ]
        ),
    ], color="light"
    )

    card_6 = dbc.Card([
        dbc.CardHeader("Gap region(s)", class_name="small"),
        dbc.CardBody(
            [
                html.H5(alnstats['Num_of_Gap'], id='info-num-gaps', className="card-title text-center"),
            ]
        ),
    ], color="light"
    )

    info_cards = [
            dbc.Col(card_1, width=6, md=4, xl=2),
            dbc.Col(card_2, width=6, md=4, xl=2),
            dbc.Col(card_3, width=6, md=4, xl=2),
            dbc.Col(card_4, width=6, md=4, xl=2),
            dbc.Col(card_5, width=6, md=4, xl=2),
            dbc.Col(card_6, width=6, md=4, xl=2),
    ]

    return info_cards



# variant-mutation comparison
dropdown_menu_variant = html.Div([
        html.P(
            'Select variant(s)',
        ),
        dcc.Dropdown(
            id='variant-selection',
            options=[],
            # options=[{"label": i, "value": i} for i in df_lineage_info.lineage],
            multi=True,
            # value=df_lineage_info.lineage,
        ),
])

dropdown_menu_mutation = html.Div([
        html.P(
            'Select mutation(s)',
        ),
        dcc.Dropdown(
            id='mutation-selection',
            # options=[{"label": i, "value": i} for i in df_lineage_mutation.mutation.unique()],
            multi=True,
            # value=[],
        ),
])

slider_num_variant = html.Div([
        html.P(
            'Mutation(s) found in variants'
        ),
        dcc.Slider(
            id='num-of-shared-variants',
            min=1,
            max=1,
            # max=len(df_lineage_info.lineage),
            value=1,
            # marks={},
            tooltip={"placement": "bottom", "always_visible": True},
        )
])


# Dash app conponents

# header
header_div = html.Div(
    [
        html.H4(
            id='header-title'
        ),
        html.H6(
            ['EDGE Bioinformatics SARS-CoV-2 analysis'],
            id='header-info'
        )
    ],
    className="mb-3 mt-3"
)

# filters and controllers
dropdown_menu_gene = html.Div([
        html.P(
            'Select a gene',
        ),
        dcc.Dropdown(
            id='gene-selection',
            options=[{"label": "All", "value": "all"}]+[{"label": i, "value": gene[i]} for i in gene],
            value='all'
        ),
])

slider_genome = html.Div([
        html.P(
            'Genome range'
        ),
        dcc.RangeSlider(
            id='genome-range',
            min=1,
            max=29903,
            value=[1, 29903],
            tooltip={"placement": "bottom", "always_visible": True},
            marks={
                266:   {'label': 'ORF1a' },
                13468: {'label': 'ORF1ab'},
                21563: {'label': 'S'     },
                # 25393: {'label': 'ORF3a' },
                # 26245: {'label': 'E'     },
                26523: {'label': 'M'     },
                # 27202: {'label': 'ORF6'  },
                # 27394: {'label': 'ORF7a' },
                # 27756: {'label': 'ORF7b' },
                # 27894: {'label': 'ORF8'  },
                28274: {'label': 'N'     },
                #29558: {'label': 'ORF10' }
            }
        )
])

slider_nonref_pct = html.Div([
        html.P(
            'Minimal non-ref %'
        ),
        dcc.Slider(
            id='min-nonref-pct',
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
])

slider_nonref_depth = html.Div([
        html.P(
            'Minimal non-ref depth',
        ),
        dcc.Slider(
            id='min-nonref-depth',
            min=0,
            max=100,
            value=0,
            marks={
                0: {'label': 'No limit'},
                50: {'label': '50x'},
                100: {'label': '100x'}
            },
            tooltip={"placement": "bottom", "always_visible": True},
        )
])

# alert for unavailable SRA
alert = dbc.Alert(
        "This is an alert for structure viewer",
        id="mol3d-alert",
        dismissable=True,
        fade=True,
        is_open=False,
        duration=4000,
        color="danger",
        className="w-100 inblock-alert",
        # style={"position": "fixed", "top": 80, "right": 15, "width": 500},
)

# warning modal for unavailable SRA
warning_modal = html.Div(
    [
        dbc.Modal(
            [
                dbc.Alert(
                    [
                        html.H2("Oooops..."),
                        html.Hr(),
                        html.P([
                            "We can't find the project data for SRA ",
                            html.Span(id="modal-sra-acc"),
                            ". Click ",
                            html.A(
                                "here", href="/browser",
                                # className="text-decoration-none"
                            ),
                            " to navigate back to the main page.",
                        ])
                    ],
                    dismissable=False,
                    fade=False,
                    is_open=True, 
                    color="danger",
                )
            ],
            id="warning-modal",
            centered=True,
            fullscreen=True,
            is_open=False,
        )
    ]
)

# Dash layout
layout = html.Div([
    dbc.Container([
        #dcc.Location(id='url', refresh=False),
        dcc.Store(id='aggregate-data'),
        dcc.Store(id='sra-acc'),
        dcc.Store(id="callback-holder"),
        warning_modal,
        html.Div(
            [
                header_div
            ],
        ),
        dbc.Row(
            id="info-cards",
            class_name="mt-5",
            align="start",
        ),
        dbc.Row(
            [
                dbc.Col(dcc.Graph(id="snp-overview-graph"), md=12),
            ],
            class_name="mt-3",
            align="start",
        ),
        dbc.Row(
            [
                dbc.Col(slider_genome, width=6, lg=3),
                dbc.Col(slider_nonref_pct, width=6, lg=3),
                dbc.Col(slider_nonref_depth, width=6, lg=3),
                dbc.Col(dropdown_menu_gene, width=6, lg=3),
            ],
            class_name="mt-3",
            align="start",
        ),
        dbc.Row(
            [
                dbc.Col([
                          html.H5("SNP/InDel"), 
                          html.P("SNPs and InDels are analyzed on EDGE COVID-19 pipeline using majority rule. The criteria are at least 5X depth of coverage, >50% allele frequency for SNPs, >50% for Indels from illumina, and >60% for Indels from ONT.", className='small'), 
                          alert,
                          snp_table
                        ], 
                        width=12, lg=8
                ),
                dbc.Col([ 
                          dcc.Store(id="dashbio-mol3d-default"),
                          html.Div([
                              html.H5("S protein structure ", className='ml-3 d-inline'),
                              dbc.Button("reset", id="mol3d-reset-btn", color="secondary", class_name="ml-3 btn-sm d-inline"),
                          ]),
                          html.Div([mol3d_spike], id="mol3d-div"),
                        ],
                        width=12, lg=4
                ),
            ],
            class_name="mt-5",
            align="start",
        ),
        dbc.Row(
            [
                dbc.Col([
                          html.H5("Variant / S protein mutation"),
                          html.P(["This plot summerizes the reference mutations of VOC/VOI and this sample. ",
                                  "The reference mutations for SARS-CoV-2 VOC/VOI are compiled by ",
                                  html.A("cov-lineages.org", href = "https://cov-lineages.org/constellations.html"),
                                  ". By clicking the bar on the right marginal plot, the protein structure viewer will highlight all mutations of the variant. ",
                                  "The single mutation will be highlighted when you click a marker.",
                                 ], 
                                 className='small'
                          ),
                        ], width=12)
            ],
            class_name="mt-5",
            align="start",
        ),
        dbc.Row(
            [
                dbc.Col(dropdown_menu_variant, width=12, lg=4),
                dbc.Col(dropdown_menu_mutation, width=12, lg=4),
                dbc.Col(slider_num_variant, width=12, lg=4),
            ],
            class_name="mt-1 small",
            align="start",
        ),
        dbc.Row(
            [
                dbc.Col([dcc.Graph(id="variant-mutation-graph")], width=12),
            ],
            class_name="mt-1",
            align="start",
        ),
        dbc.Row(
            [
                dbc.Col([
                          html.H5("Genome viewer ", className='d-lg-inline'),
                          dbc.Button(
                              "Sync range",
                              id="igv-range-btn",
                              color="primary",
                              class_name="btn-sm",
                          )
                        ]),
            ],
            id="igv-session",
            class_name="mt-5",
            align="start",
        ),
        dbc.Row(
            [
                dbc.Col(id='igv-graph', md=12),
            ],
            class_name="mt-3",
            align="start",
        ),
    ],
    class_name="p-5 m-3",
    fluid=True)
])

# Create callbacks

# checking SRA# in the URL
@app.callback(Output('modal-sra-acc', 'children'),
              Output('warning-modal', 'is_open'),
              [
                  Input('url', 'pathname'),
                  Input('url', 'search'),
              ]
)
def check_url(pathname, search):
    if '/sample' in pathname and '?sra=' in search:
        # parse sra# from URL
        acc_num = search.replace('?sra=', '')

        if not exists(f'assets/ec19_output/{acc_num}/Output/config.txt'):
            return [acc_num, True]
    
    return [dash.no_update, dash.no_update]


# init display page
@app.callback(Output('sra-acc', 'data'),
              Output('variant-selection', 'options'),
              Output('variant-selection', 'value'),
              Output('mutation-selection', 'options'),
              Output('num-of-shared-variants', 'max'),
              Output('num-of-shared-variants', 'marks'),
              Output('info-cards', 'children'),
              Output('header-title', 'children'),
              Output('igv-graph', 'children'),
              Output('dashbio-mol3d-default', 'data'),
              [
                  Input('url', 'search'),
              ]
)
def init_display_page(search):
    
    global df_lineage_info
    
    if not search:
        return [dash.no_update]*10
        
    # parse sra# from URL
    acc_num = search.replace('?sra=', '')
    
    # init data
    init_proj_data(acc_num)

    # variant-selection
    vs_options=[{"label": i, "value": i} for i in df_lineage_info.lineage]
    vs_value=df_lineage_info.lineage.tolist()

    # mutation-selection
    ms_options=[{"label": i, "value": i} for i in df_lineage_mutation.mutation.unique()]

    # num-of-shared-variants
    sv_max=len(df_lineage_info.lineage)
    sv_marks={
        1: {'label': 'No limit'},
        len(df_lineage_info): {'label': str(len(df_lineage_info))}
    }
    
    # info cards
    info_cards = generate_layout_info_cards()
    
    # igv-graph
    igv = generate_igv_graph()
    
    # dashbio-mol3d-default
    mol3d_default = generate_mol3d_spike_default()
    
    return [acc_num, vs_options, vs_value, ms_options, sv_max, sv_marks, info_cards, alnstats["SAMPLE"], igv, mol3d_default]


# Create callbacks
@app.callback(Output('aggregate-data', 'data'),
              [
                  Input('genome-range', 'value'),
                  Input('min-nonref-pct', 'value'),
                  Input('min-nonref-depth', 'value'),
                  Input('snp-table', 'selected_rows'),
                  Input('snp-overview-graph', 'clickData'),
                  Input('variant-mutation-graph', 'clickData'),
                  Input('header-title', 'children'),
              ],
              prevent_initial_call=True
            )
def update_aggregate_data(genome_range, min_nonref_pct, min_nonref_depth, selected_rows, so_click_data, vm_click_data, header_title):
    # need to convert selected_rows -> spike pos

    return {'start': genome_range[0], 
            'end': genome_range[1], 
            'min_nonref_pct': min_nonref_pct/100, 
            'min_nonref_depth': min_nonref_depth,
            'snp_overview_click_data': so_click_data,
            'var_mut_click_data': vm_click_data,
           }


# Gene drop down -> range
@app.callback(Output('genome-range', 'value'),
              [
                Input('gene-selection', 'value'),
                Input('snp-overview-graph', 'relayoutData'),
              ],
              prevent_initial_call=True
             )
def update_genome_range(gene_sel, relayout_data):
    
    ctx = dash.callback_context

    triggered_id = None
    if ctx.triggered:
        triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]
        
    if triggered_id == 'gene-selection':
        if gene_sel and '-' in gene_sel:
            (start, end) = gene_sel.split('-')
            return [int(start), int(end)]


    if triggered_id == 'snp-overview-graph':    
        if relayout_data and 'xaxis.range[0]' in relayout_data:
            return [relayout_data['xaxis.range[0]'], relayout_data['xaxis.range[1]']]

    return [1, 29903]

    
# aggregate-data -> snp-overview-graph
@app.callback(Output('snp-overview-graph', 'figure'),
              Output('snp-table', 'data'),
              [
                  Input('aggregate-data', 'data')
              ],
              prevent_initial_call=True
             )
def update_figures(data):
    
    fig = generate_overview_graph(data['min_nonref_pct'], data['min_nonref_depth'], data['start'], data['end'])
    
    if len(df_snps['Count'].isna())>0:
        df = df_snps.query(f"SNP_position>={data['start']} & SNP_position<={data['end']}").copy()
    else:
        df = df_snps.query(f"Count >= {data['min_nonref_depth']} & SNP_position>={data['start']} & SNP_position<={data['end']}").copy()

    return [fig, df.to_dict('records')]


# snp-table -> mol3d
@app.callback(Output('dashbio-mol3d', 'zoom'),
              Output('dashbio-mol3d', 'zoomTo'),
              Output('dashbio-mol3d', 'styles'),
              Output('dashbio-mol3d', 'labels'),
              Output('mol3d-alert', 'children'),
              Output('mol3d-alert', 'is_open'),
              [
                  Input('dashbio-mol3d-default', 'data'),
                  Input('snp-table', 'active_cell'),
                  Input('variant-mutation-graph', 'clickData'),
                  Input('mol3d-reset-btn', 'n_clicks'),
                  State('dashbio-mol3d', 'zoom'),
                  State('dashbio-mol3d', 'zoomTo'),
                  State('dashbio-mol3d', 'styles'),
                  State('dashbio-mol3d', 'labels'),
              ],
              prevent_initial_call=True
             )
def update_mol3d(mol3d_default, snp_selected_cell, 
                    vm_click_data, mol3d_reset,
                    zoom, zoomTo, styles, labels):
    
    global spike_pdb
    global df_var_mutation_plot_data
    global df_snps
    import re
    
    alert_content = ""
    alert_is_open = False
        
    ctx = dash.callback_context

    triggered_id = None
    if ctx.triggered:
        triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]

    # init mol3d to default settings
    if triggered_id == 'dashbio-mol3d-default' or triggered_id == 'mol3d-reset-btn':
        styles = mol3d_default['styles']
        labels = mol3d_default['labels']
        zoomTo = mol3d_default['zoomTo']
        zoom   = mol3d_default['zoom']
    
    # clicking "Mutation" on snp table
    if triggered_id == 'snp-table':
        
        if snp_selected_cell and snp_selected_cell['row_id']:            
            snp_pos = snp_selected_cell['row_id']    
            product = df_snps.set_index('SNP_position').loc[snp_pos,'Product']
            
            if product=='S' and snp_selected_cell['column_id']=='Mutation':
                pos = int(df_snps.set_index('SNP_position').loc[snp_pos, 'AA_pos'])
                # convert spike position to pdb_residue_index_chain
                residue_index = spike_pdb.spike2pdb_residue_index_chain(pos)
                
                zoomTo = {
                    "sel": {"chain": "A", "resi": residue_index},
                    "animationDuration": 0,
                    "fixedPath": False,
                }
            elif snp_selected_cell['column_id']=='SNP_position':
                # clicking on 'Position' field will trigger another callback to
                # change the locus to the clicked position.
                pass
            elif snp_selected_cell['column_id']=='Mutation':
                alert_content = "The protein structure viewer only supports S protein currently."
                alert_is_open = True
            else:
                alert_content = "Please click 'Position' and 'Mutation' field."
                alert_is_open = True            
    
    # Clicking on variant-mutation table
    #     -- highlight selections and the sample on Mol3D:
    # 1. Identifying the variant or SNP is being selected
    # 2. Getting spike pos and converting to residue_indexes
    # 3. Generating mol3d styles
    
    if triggered_id == 'variant-mutation-graph' and vm_click_data: 
        
        x = vm_click_data['points'][0]['x']
        y = vm_click_data['points'][0]['y']
        var_spike_pos = None
        
        # y-barchart variant selected
        if isinstance(x, int):
            # getting spike pos and converting to residue_indexes
            var_spike_pos = df_var_mutation_plot_data.query(f"lineage=='{y}'").pos.astype('int').tolist()
            var_spike_mut = df_var_mutation_plot_data.query(f"lineage=='{y}'").mutation.tolist()
        elif x.startswith("S:"):
            import re
            try:
                pos = re.search('(\d+)', x).group(1)
                var_spike_pos = [int(pos)]
                var_spike_mut = [x]
            except AttributeError:
                pass

        # var_spike_pos = [int(float(x)) for x in var_spike_pos]
        var_residue_indexes = spike_pdb.spike2pdb_residue_indexes(var_spike_pos)
        
        # current sample
        sample_spike_pos = df_var_mutation_plot_data.query("type=='Sample'").pos.astype('int').tolist()
        sample_spike_mut = df_var_mutation_plot_data.query("type=='Sample'").mutation.tolist()
        sample_residue_indexes = spike_pdb.spike2pdb_residue_indexes(sample_spike_pos)
        
        # generate residue_indexes color mapping dict
        residue_indexes_color_map = {}

        for x in var_residue_indexes:
            residue_indexes_color_map[x] = 'orange'

        for x in sample_residue_indexes:
            residue_indexes_color_map[x] = '#459DF8'
            if x in var_residue_indexes:
                 residue_indexes_color_map[x] = '#FD4546'
                    
        # generate s positions color mapping dict
        poss_color_map = {}
                
        for x in var_spike_pos:
            poss_color_map[x] = 'orange'

        for x in sample_spike_pos:
            poss_color_map[x] = '#459DF8'
            if x in var_spike_pos:
                 poss_color_map[x] = '#FD4546'
        
        # generate labels pos and text
        label_text = list(set(var_spike_mut+sample_spike_mut))
        label_pos  = [int(re.search('(\d+)', x).group(1)) for x in label_text]
        
        # hightlight RBD
        opts = dict(highlight_bg_indexes = spike_pdb.spike2pdb_residue_indexes(range(319,541)))
        
        styles = spike_pdb.pdb_style(residue_indexes_color_map=residue_indexes_color_map, **opts)
        labels = spike_pdb.residue_labels(label_pos, label_text, poss_color_map=poss_color_map, font_size=10)
        zoomTo = mol3d_default['zoomTo']

    return [zoom, zoomTo, styles, labels, alert_content, alert_is_open]


# select variant -> variant-mutation plot
@app.callback(
    Output('variant-mutation-graph', 'figure'),
    [
        Input('variant-selection', 'value'),
        Input('mutation-selection', 'value'),
        Input('num-of-shared-variants', 'value'),
    ],
    prevent_initial_call=True
)
def update_variant_mutation(variants, mutations, num_shared):
    
    global df_var_mutation_plot_data
    global df_lineage_mutation
    
    df_lineage_mutation_filtered = df_lineage_mutation
 
    # filter in variants
    if variants and len(variants):
        idx = df_lineage_mutation_filtered.lineage.isin(variants)
        df_lineage_mutation_filtered = df_lineage_mutation_filtered[idx]
    
    # filter in mutations
    if mutations and len(mutations):
        idx = df_lineage_mutation_filtered.mutation.isin(mutations)
        df_lineage_mutation_filtered = df_lineage_mutation_filtered[idx]
    
    # filter in shared mutations
    df_mut_lineage_cnt = df_lineage_mutation.groupby('mutation').count()[['lineage']].reset_index()
    shared_mutations = df_mut_lineage_cnt.query(f'lineage >= {num_shared}').mutation.to_list()
    
    if num_shared and len(shared_mutations):
        idx = df_lineage_mutation_filtered.mutation.isin(shared_mutations)
        df_lineage_mutation_filtered = df_lineage_mutation_filtered[idx]
    
    mock_var = df_snps.query('Synonymous=="No" & Product=="S"')[['SAMPLE', 'Mutation']]
    mock_var_rows = [{'lineage': v[0], 'mutation': v[1], 'type': 'Sample'} for idx, v in mock_var.iterrows()]
    df = pd.concat([df_lineage_mutation_filtered, pd.DataFrame(mock_var_rows)])
    df['pos'] = df.mutation.str.extract(r'(\d+)').astype(int) # add position to df for ordering x-axis labels
    
    # merging mutation categories
    df = df.merge(df_mutation_category, on='mutation', how='left')
    df = df.fillna(value='')
    
    # save data
    df_var_mutation_plot_data = df
    
    return generate_var_mutation_plot(df)


# aggregate-data -> igv
@app.callback(
    Output('reference-igv', 'locus'),
    Output('reference-igv', 'reference'),
    [
        Input('snp-table', 'active_cell'),
        Input('igv-range-btn', 'n_clicks'),
        Input('snp-overview-graph', 'clickData'),
        State('aggregate-data', 'data'),
        State('reference-igv', 'locus'),
        State('reference-igv', 'reference'),
    ],
    prevent_initial_call=True
)
def igv_update(snp_selected_cell, btn_clicks, so_click_data, data, locus, igv_reference):
    
    ctx = dash.callback_context

    if ctx.triggered:
        triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]
        
        if triggered_id == 'snp-overview-graph':
            point_x = so_click_data['points'][0]['x']
            locus = f"NC_045512_2:{point_x}-{point_x}"
            
        elif triggered_id == 'igv-range-btn':
            locus = f"NC_045512_2:{data['start']}-{data['end']}"

        elif triggered_id == 'snp-table':
            if snp_selected_cell and snp_selected_cell['column_id']=='SNP_position':                
                pos = snp_selected_cell['row_id'] # 'row_id' is 'SNP_position'
                locus = f"NC_045512_2:{pos}-{pos}"

        # elif triggered_id == 'header-title':
        #         igv_reference['tracks'] = generate_igv_graph()
    
    return [locus, igv_reference]

# click snp --> scrollIntoView IGV
app.clientside_callback(
    """
    function(clicks, elemid) {
        document.getElementById(elemid).scrollIntoView({
            behavior: 'smooth'
        });
    }
    """,
    Output('callback-holder', 'data'),
    [
        Input('snp-overview-graph', 'clickData'),
        State('igv-session', 'id'),
    ],
    prevent_initial_call=True
)

if __name__ == '__main__':
    app.run_server(mode='external', debug=True, port=8053, 
                   #threaded=True
                  )