import dash
import dash_bootstrap_components as dbc

external_stylesheets = [dbc.themes.BOOTSTRAP]

app = dash.Dash(__name__, 
    external_stylesheets = external_stylesheets,
    suppress_callback_exceptions = True
)
app.title = "SARS-CoV-2 SRA"
app._favicon = "favicon.ico"

# acc_num = 'SRR14135907'

server = app.server