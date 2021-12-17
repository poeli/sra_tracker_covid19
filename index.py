from re import S
from typing import Container
from dash import Input, Output, State, dcc, html
import dash_bootstrap_components as dbc

from app import app
# import all pages in the app
from apps import covid_sra, ec_sra

# building the navigation bar
# https://github.com/facultyai/dash-bootstrap-components/blob/master/examples/advanced-component-usage/Navbars.py
dropdown = dbc.DropdownMenu(
    id = "example-menu",
    children=[
        dbc.DropdownMenuItem("SRA Browser", href="/browser"),
        dbc.DropdownMenuItem("SRR14135907", href="/sample?sra=SRR14135907"),
        dbc.DropdownMenuItem("SRR12529581", href="/sample?sra=SRR12529581"),
        dbc.DropdownMenuItem("ERR5971761 (AY.4 ONT)", href="/sample?sra=ERR5971761"),
        dbc.DropdownMenuItem("ERR5669404 (B.1.1.7 ONT)", href="/sample?sra=ERR5669404"),

    ],
    nav = True,
    in_navbar = True,
    label = "example",
)

navbar = dbc.Navbar(
    dbc.Container(
        [
            html.A(
                # Use row and col to control vertical alignment of logo / brand
                dbc.Row(
                    [
                        dbc.Col(html.Img(src="/assets/virus-solid-white.png", height="27px")),
                        dbc.Col(dbc.NavbarBrand("SARS-CoV-2 SRA", className="ml-1")),
                    ],
                    align="center",
                ),
                href="/browser",
                className="text-decoration-none"
            ),
            dbc.NavbarToggler(id="navbar-toggler2"),
            dbc.Collapse(
                dbc.Nav(
                    # right align dropdown menu with ml-auto className
                    [dropdown], className="ml-auto", navbar=True
                ),
                id="navbar-collapse2",
                navbar=True,
            ),
            html.Div([
                dcc.Dropdown(
                    id="acc-select-nav",
                    options=[{'label': 'test', 'value': 'test'}],
                    placeholder="Select SRA/ERA/DRA",
                    style={'width': '15rem'},
                    className='small',
                )],
                className='d-flex'
            ),
        ]
    ),
    color="dark",
    dark=True,
    className='w-100 fixed-top container-fluid',
)

def toggle_navbar_collapse(n, is_open):
    if n:
        return not is_open
    return is_open

for i in [1, 2]:
    app.callback(
        Output(f"navbar-collapse{i}", "is_open"),
        [Input(f"navbar-toggler{i}", "n_clicks")],
        [State(f"navbar-collapse{i}", "is_open")],
    )(toggle_navbar_collapse)

# embedding the navigation bar
app.layout = html.Div([
    dcc.Store(id='clientside-callback-temp'),
    dcc.Location(id='url', refresh=False),
    navbar,
    html.Div(id='page-content')
])


@app.callback(Output('page-content', 'children'),
              [
                  Input('url', 'pathname'),
              ]
)
def display_page(pathname):
    if pathname == '/browser':
        return covid_sra.layout
    elif pathname == '/sample':
        #SRR14135907
        return ec_sra.layout
    else:
        return covid_sra.layout



# select a sra from nav bar --> link to /sample?sra=acc
app.clientside_callback(
    """
    function(acc) {
        window.location.href = "/sample?sra="+acc
    }
    """,
    Output('clientside-callback-temp', 'data'),
    Input('acc-select-nav', "value"),
    prevent_initial_call=True
)



if __name__ == '__main__':
    app.run_server(#port=8877,
                   threaded=True,
                   #debug=True
                  )
