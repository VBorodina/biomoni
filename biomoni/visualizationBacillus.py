import sys
sys.path.append("V:/biomoni/")  

from biomoni import Experiment

import pandas as pd
from plotly.subplots import make_subplots
import plotly.graph_objs as go


def visualizeBacillusFermentation(sim_data, Fdata= None, title= None,  suffix_1 = "", suffix_2 = "_fitted", mode_1 = "markers", mode_2 = "lines", secondary_y_cols = ["rP",("STIRR_2","Value"),"mu","rO2","Glc [g/L]","V",("BASET_2","Value"),"Acetate [g/L]"]):
    """X"""
    
    suffix = {1: suffix_2, 2: suffix_1}
    mode = {1: mode_1, 2: mode_2}
    data_both = {1: Fdata, 2: sim_data}
    data_dict ={}
    #sec_cols = ["rP",("STIRR_2","Value"),"mu","rO2","Glc [g/L]","V",("BASET_2","Value"),"Acetate [g/L]"]      
     
    col_dict = {("RF") : "rgb(239,170,33)", ("rP") : "rgb(158,112,21)", ("CO2") : "rgb(0,153,76)", ("RQ"):"rgb(0,51,25)",("STIRR_2","Value"):"rgb(0,153,153)", 
                ("CDW"):"rgb(0,204,102)", ("mu"): "rgb(128,255,0)", ("pO2_2","Value"):"rgb(0,0,153)", ("rO2"):"rgb(0,128,255)",("Glc [g/L]"):"rgb(0,0,102)", 
                ("rGlc"):"rgb(0,102,102)",("Feed_rate"):"rgb(0,102,0)", ("FOAMT_2" , "Value"):"rgb(0,51,25)", ("V"):"rgb(153,0,76)", ("F_out"):"rgb(255,0,127)", 
                ("Baserate"):"rgb(32,32,32)", ("pH_2", "Value"):"rgb(102,0,0)", ("BASET_2","Value"):"rgb(0,0,0)", ("Acetate [g/L]"):"rgb(255,255,0)",
                ("Peak 22.84"):"rgb(153,153,0)",("Peak 24.21"):"rgb(102,102,0)"}
    
    
    for i, data in data_both.items():
        assert isinstance(data, (pd.DataFrame, dict, list, Experiment, type(None))), "Data should be given either as an Experiment, a pd.DataFrame or as a list or dict with pd.DataFrame objects"

        if data is None:
            data_dict[i] = data
        
        elif type(data) is Experiment:
            data_dict[i] = list(data.dataset.values())


        elif type(data) is dict:
            assert all(isinstance(i, pd.DataFrame ) for i in data.values()) , "Each value in data_{0} has to be of type pd.DataFrame".format(i)
            data_dict[i] = list(data.values())
                

        elif type(data) is list:
            assert all(isinstance(i, pd.DataFrame ) for i in data) , "Each element in data_{0} hasS to be of type pd.DataFrame".format(i)
            data_dict[i] = data
                
        elif type(data) is pd.DataFrame: 
            data_dict[i] = [data]
      
    #fig = make_subplots(rows= 4, cols=2, shared_xaxes ="all", start_cell = "top-left", x_title = "F-time [h]", specs =[[{"secondary_y" : True},{"secondary_y" : True}],[{"secondary_y" : True},{"secondary_y" : True}],[{"secondary_y" : True},{"secondary_y" : True}],[{"secondary_y" : True},{"secondary_y" : True}]] )
    #fig = make_subplots(rows= 4, cols=2, specs =[[{"secondary_y": True},{"secondary_y": True}],[{"secondary_y": True},{"secondary_y": True}],[{"secondary_y": True},{"secondary_y": True}],[{"secondary_y": True},{"secondary_y": True}]] )
  
    fig = make_subplots(rows= 4, cols=2, specs =[[{"secondary_y": True},{"secondary_y": True}],
                                                 [{"secondary_y": True},{"secondary_y": True}],
                                                 [{"secondary_y": True},{"secondary_y": True}],
                                                 [{"secondary_y": True},{"secondary_y": True}]],
                        x_title = "F-time [h]")
    
    
    fig.show()
    
    for i, data_list in data_dict.items():
        if data_list is None:
            pass
        
        else: 
            for df in data_list:
                
                for col in df:
                    
                    secondary_y_flag = col in secondary_y_cols
                    if col in col_dict.keys():
                        if type(col) is tuple:
                            fig.add_trace(go.Scatter(x= df.index, y= df[col], mode = mode[i], marker = dict(color = col_dict[col], size = 5, symbol = "x")), 
                                          secondary_y= secondary_y_flag)
                        if type(col) is str:
        
                           fig.add_trace(go.Scatter(x= df.index, y= df[col], name = col+ suffix[i], mode = mode[i], marker = dict(color = col_dict[col], size = 5, symbol = "x")), 
                                          secondary_y= secondary_y_flag)
                        
                        



    cols_y1 = [col for col in col_dict.keys() if col not in secondary_y_cols]            
    cols_y2 = [col for col in col_dict.keys() if col in secondary_y_cols]
    #fig.update_yaxes(title_text= str(cols_y1), secondary_y=False , title_standoff = 20 )
    #fig.update_yaxes(title_text= str(cols_y2), secondary_y=True, title_standoff = 20)

    fig.update_xaxes(title_text = "Time [h]")
    fig.update_layout(title_text = title , title_x=0.5)
    fig.update_layout(legend=dict(x = 1))
    fig.show()
