from biomoni import Experiment, Yeast
from visualize_dash import visualize
import os 
import pathlib
from datetime import datetime
import numpy as np
import pandas as pd

import dash
from dash import dcc
from dash import html
from dash import dash_table
from dash.dependencies import Input, Output
from settings_dash import kwargs_experiment, kwargs_estimate, measurement_vars, simulated_vars, fitparams

from dash.long_callback import DiskcacheLongCallbackManager
import diskcache


from plotly.subplots import make_subplots
import plotly.graph_objs as go
import json



#external_stylesheets = ["https://codepen.io/chriddyp/pen/bWLwgP.css"]   #necessary?

all_vars = set([*measurement_vars, *simulated_vars])    #All variables only once


cache = diskcache.Cache("./cache")
long_callback_manager = DiskcacheLongCallbackManager(cache)

app = dash.Dash(__name__, long_callback_manager = long_callback_manager)     #external_stylesheets = external_stylesheets

Result_path = r"P:\Code\biomoni\Messdaten\OPCUA"        #pfad kann in settings.py
#Result_path = "/home/paul/Desktop/pCloudDrive/Code/biomoni/Messdaten/OPCUA"
Result_path = "/home/paul/Desktop/lel/Code/biomoni/Messdaten/OPCUA"


sub_paths = next(os.walk(Result_path))[1]       #yields the subsirectory in the given path
newest_results_dir = max([os.path.join(Result_path,i) for i in sub_paths], key=os.path.getmtime) #gives newest subdirectory

colors = {
    "background": "oxy",
    "text": "green",
    "settings" : "lightyellow",
    "table_header" : "red",
    "table_background" : "grey"
}

path = Result_path  #überflüssig
exp_dir_manual = newest_results_dir

#read the csv file with settings
Exp = Experiment(path = Result_path, exp_dir_manual = exp_dir_manual, **kwargs_experiment["online_est"])

y = Yeast()
y.estimate(Exp, tau = 1, max_nfev = 100) 
p_dict = {p : val.value for p, val in y.p.items()}
all_params = [i for i in y.p]

#dataframe for displaying experimental data in Table
df = Exp.dataset["on_CO2"]


#Available columns in the Dropdown menu
#available_indicators = [i for i in df.columns if i in measurement_vars]

def generate_table(dataframe, no_cols = []):          #Table to show measurement data numerically
    "Function to generate HTML table"
    columns = [col for col in dataframe.columns if col not in no_cols]
    return html.Table([
        html.Thead(
            html.Tr([html.Th(col) for col in columns])
        ),
        html.Tbody([
            html.Tr([
                html.Td(dataframe.iloc[i][col]) for col in columns
            ]) for i in range(len(dataframe))       #, style = {"color" : colors["text"]}
        ])
    ])

def make_list(input):
    "Function to assure that the selected variables are given in form of a list for the visualize function"
    assert isinstance(input, (list, str, type(None))), "colum types must be a str a list of str or None"
    if input is not None:
        if type(input) is list:
            pass
        elif type(input) is str:
            input = [input]
    else:
        input = []
    return input

def re_estimate(params):
    "Function to re-estimate params"
    for p in all_params:
        if p in params:
            y.change_params(p, vary= True)
        else:
            y.change_params(p, vary = False)

    y.estimate(Exp, tau = 1, max_nfev = 100) 


#Layout
app.layout = html.Div(style={"backgroundColor": colors["background"]}, children=[

    dcc.Interval(
    id = "interval",
    interval = 120 * 1000,      #milliseconds 
    n_intervals = 0

    ),

    dcc.Store(id = "data_store"),

    html.H1(
        children= "Biomonitoring Dashboard",
        style={
            "textAlign": "center",
            "color": colors["text"]
        }
    ),

    html.Div(children= ["A web application framework for your data.", 
        dcc.Markdown("""For more information visit [biomoni](https://github.com/PSenck/biomoni)""")     #htlm.A geht auch für link
        ], style={
        "textAlign": "center",
        "color": colors["text"]
    }),

    html.Br(),

    html.Div([
        
        html.Div([
            html.Div("Displayed variables of measurement data", style =  {"color" : colors["text"]}),
            dcc.Dropdown(
                id='meas_vars',
                options=[{'label': i, 'value': i} for i in measurement_vars],
                value = "CO2",
                multi = True,
        
            ),
            html.Div("Displayed variables of simulated data", style =  {"color" : colors["text"]}),
            dcc.Dropdown(
                id='sim_vars',
                options=[{'label': i, 'value': i} for i in simulated_vars],
                value = "CO2",
                multi = True,
        
            ),
            
            dcc.RadioItems(
                id='yaxis_type',
                options=[{'label': i, 'value': i} for i in ['linear', 'log']],
                value='linear',
                labelStyle={'display': 'inline-block'},
                style = {"color" : colors["text"]},
                
                
            )
        ], style={'width': '48%', 'display': 'inline-block'}),
        
        html.Div([
            html.Div("Variables to be displayed on the secondary_yaxis", style =  {"color" : colors["text"]}),
            dcc.Dropdown(
                id='secondary_yaxis',
                options=[{'label': i, 'value': i} for i in all_vars],
                value =  "CO2",
                multi = True
            ),
        

            dcc.RadioItems(
                id='secondary_yaxis_type',
                options=[{'label': i, 'value': i} for i in ['linear', 'log']],
                value='linear',
                labelStyle={'display': 'inline-block'},
                style = {"color" : colors["text"]}
            )
        ], style={'width': '48%', 'float': 'right', 'display': 'inline-block'})
    ], style = {"backgroundColor" : colors["settings"]}),

    dcc.Graph(
        id = "graph1",
        figure={} 

    ),
    html.Div([
        "Simulation time in hours: ",
        dcc.Input(id= "sim_time", value= 10, type='number')
    ], style = {"color": colors["text"]}),


    html.Div([html.P(id="paragraph_id", children=["Button not clicked"])]),
    html.Button(id="button_id", children="Estimate!"),

    dash_table.DataTable(
        id="table_params",
        columns= [],              #[{"name": i, "id": i} for i in p_dict.keys()],
        data= [],     #pd.DataFrame(p_dict, index = [0])
        style_header={
        "color" : colors["text"],
        'backgroundColor': colors["table_background"],
        'fontWeight': 'bold'
    },
    ),

    html.Div("This is before any iterations of dcc.Interval", id = "iteration_identifier")



])

@app.callback(
    Output("graph1", "figure"),
    Input("data_store", "data"),
    Input("meas_vars", "value"),
    Input("sim_vars", "value"),
    Input("secondary_yaxis", "value"),
    Input("yaxis_type", "value"),
    Input("secondary_yaxis_type", "value"))
def update_graph_1(jsonified_data, meas_vars, sim_vars, secondary_yaxis, yaxis_type, secondary_yaxis_type):
    if jsonified_data is not None:
        
        dataset = json.loads(jsonified_data)
        simulated_data = pd.read_json(dataset["simulated_data"], orient = "split").rename_axis("t")
        simulated_data = simulated_data.filter(items = sim_vars)
        measured_data = {}
        for typ, dat in dataset["measured_data"].items():
            measured_data[typ] = pd.read_json(dat, orient = "split").rename_axis("t")
            filter_cols = [col for col in  measured_data[typ].columns if col in meas_vars]
            measured_data[typ] = measured_data[typ].filter(items = filter_cols)

        #fig = visualize(measured_data, simulated_data)
        fig = visualize(measured_data, simulated_data,  secondary_y_cols= secondary_yaxis, yaxis_type = yaxis_type, sec_yaxis_type = secondary_yaxis_type )

    if jsonified_data is None:
        fig = {}

    return fig

@app.callback(
    Output("table_params", "columns"),
    Output("table_params", "data"),
    Input("data_store", "data")
)
def update_table_params(jsonified_data):
    dataset = json.loads(jsonified_data)
    params = dataset["params"]
    params = pd.DataFrame(params).T
    col_names = [{"name": i, "id": i} for i in params.columns]

    return col_names, params.to_dict("records")


@app.long_callback(
Output("data_store", "data"),
Output("iteration_identifier", "children"),
Output("paragraph_id", "children"),
Input("interval", "n_intervals"),
Input("sim_time", "value"),
Input("button_id", "n_clicks"),
running=[
    (Output("button_id", "disabled"), True, False),
],
)
def create_data(n_intervals, hours, n_clicks):
    ctx = dash.callback_context
    #last_input = ctx.triggered[0]["prop_id"].split(".")[0]

    Exp = Experiment(path = Result_path, exp_dir_manual = exp_dir_manual, **kwargs_experiment["online_est"])
    y.estimate(Exp, **kwargs_estimate["online_est"]) 
    dataset = Exp.dataset
    #p_dict = {p : val.value for p, val in y.p.items()}
    p_dict = {}
    for p, val in y.p.items():
        p_dict[p] = {
            "name": val.name,
            "value": val.value,
            "vary": val.vary,
            "min" : val.min,
            "max" : val.max,
            "standard_error" : val.stderr,
            "relative_error[%]" : val.stderr / val.value *100
            }
    t_grid = np.linspace(0,hours, round(hours*60)) 
    sim = y.simulate(Exp, t_grid = t_grid)

    dataset = {}
    for typ, df in Exp.dataset.items():
        dataset[typ] = df.to_json(orient = "split", date_format = "iso")
    
    all_data = {
    "measured_data" : dataset,
    "simulated_data" : sim.to_json(orient = "split", date_format = "iso"),
    "params" : p_dict
    }
    iteration_nr = "This is iteration: " + str(n_intervals)

    return json.dumps(all_data), iteration_nr, [f"Clicked {n_clicks} times"]


if __name__ == "__main__":
    app.run_server( debug=True, port = 7971)     #debug=True, 

# visit http://127.0.0.1:7971/ in your web browser.