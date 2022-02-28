import sys
sys.path.append("V:/biomoni/")   #without this line windows could not acces the package 
from plotly.subplots import make_subplots
import plotly.graph_objs as go


def visualizeBacillusFermentation(dict_for_graphs):
 """visualizeBacillusFermentation takes one dict as argument and creates a special overview of chosen fermentation parameters for bacillus experimental and simulated data.
 
 type dict_for_graphs: 
 :return: Graph with plotted and named data
 """
 dfg = dict_for_graphs
 
 col_dict = {"RF": "rgb(239,170,33)", "rP" : "rgb(158,112,21)", "CO2" : "rgb(0,153,76)", "RQ":"rgb(0,51,25)",("STIRR_2","Value"):"rgb(0,153,153)", 
                "CDW":"rgb(0,204,102)", "mu": "rgb(51,102,0)", ("pO2_2","Value"):"rgb(0,0,153)", "rO2":"rgb(0,128,255)","Glc [g/L]":"rgb(0,0,102)", 
                "rGlc":"rgb(0,102,102)","GlcFeed_rate":"rgb(0,102,0)", ("FOAMT_2" , "Value"):"rgb(0,51,25)", "V":"rgb(153,0,76)", "F_out":"rgb(255,0,127)", 
                "Baserate":"rgb(32,32,32)", "pH":"rgb(102,0,0)", ("BASET_2","Value"):"rgb(200,2,10)", "Acetate":"rgb(209,169,0)",
                "Peak 22.84":"rgb(100,200,0)","Peak 24.21":"rgb(100,50,0)"}


 for exp in dfg:
     fig = make_subplots(rows= 4, cols=2, specs =[[{"secondary_y": True},{"secondary_y": True}],
                                                  [{"secondary_y": True},{"secondary_y": True}],
                                                  [{"secondary_y": True},{"secondary_y": True}],
                                                  [{"secondary_y": True},{"secondary_y": True}]],
                         )
     #figure top left: RF and RF production rat
     fig.add_trace(go.Scatter(x=dfg[exp]["off"].index, y= dfg[exp]["off"]["RF [mg/L]"], mode="markers", name = "Riboflavin", marker= dict(color= col_dict["RF"], size= 8, symbol="x"), legendgroup ="group1"), row=1, col=1)
     fig.add_trace(go.Scatter(x=[dfg[exp]["off"].index], y= dfg[exp]["off"]["RF [mg/L]_rate"], mode="lines", name = "Riboflavin production rate", line= dict(color= col_dict["rP"], width= 1, dash= "solid"),legendgroup ="group1"),row=1, col=1, secondary_y=True)
     fig.add_trace(go.Scatter(x=dfg[exp]["simulated"].index, y=dfg[exp]["simulated"]["cP"], mode="lines", name = "Riboflavin_simulated", line= dict(color= col_dict["RF"], width= 2, dash="dot"), legendgroup ="group3"), row=1, col=1)
     fig.update_yaxes(title_text="Riboflavin [mg/L]", row=1,col=1, secondary_y =False, titlefont=dict(size=15), color=col_dict["RF"],showline=True, linewidth=2, linecolor=col_dict["RF"],showgrid=True, gridwidth=1, gridcolor="rgb(240,240,240)",zeroline=False,range = [0,5500])
     fig.update_yaxes(title_text="RF production rate [mg/L/h]", row=1,col=1, secondary_y =True, titlefont=dict(size=15), color=col_dict["rP"],showline=True, linewidth=2, linecolor=col_dict["rP"],showgrid=False, gridwidth=1, gridcolor='LightPink',zeroline=False, range = [0,150])

        #fig top right: CO2, RQ and Stirrer 
        
     if "CO2" in dfg[exp]:
         fig.add_trace(go.Scatter(x=dfg[exp]["CO2"].index, y= dfg[exp]["CO2"]["CO2"], mode="markers", name = "CO2", marker= dict(color= col_dict["CO2"], size= 3, symbol="0"),legendgroup ="group2"),row=1, col=2)
     else:
         fig.add_trace(go.Scatter(x=[], y= [], mode="markers", name = "NO CO2 available", marker= dict(color= col_dict["CO2"], size= 3, symbol="0"),legendgroup ="group2"),row=1, col=2)
        
        #fig.add_trace(go.Scatter(x=[0], y= [0], mode="markers", name = "RQ", marker= dict(color= col_dict["RQ"], size= 3, symbol="0"),legendgroup ="group2"),row=1, col=2)
     fig.add_trace(go.Scatter(x=dfg[exp]["on"].index, y= dfg[exp]["on"][("STIRR_2","Value")], mode="markers", name = "Stirrer", marker= dict(color= col_dict[("STIRR_2","Value")], size= 3, symbol="0"),legendgroup ="group2"),row=1, col=2,secondary_y=True)
     fig.update_yaxes(title_text="CO2 [%]", row=1,col=2, secondary_y =False, titlefont=dict(size=15), color=col_dict["CO2"],showline=True, linewidth=2, linecolor=col_dict["CO2"],showgrid=True, gridwidth=1, gridcolor="rgb(240,240,240)",zeroline=False,range = [0,8])
     fig.update_yaxes(title_text="Stirrer [rpm]", row=1,col=2, secondary_y =True, titlefont=dict(size=15), color=col_dict[("STIRR_2","Value")],showline=True, linewidth=2, linecolor=col_dict[("STIRR_2","Value")],showgrid=False, gridwidth=1, gridcolor='LightPink',zeroline=False,range = [500,1400])
        
        #fig 2nd row, left: CDW & mu 
     fig.add_trace(go.Scatter(x=dfg[exp]["off"].index, y=dfg[exp]["off"]["CDW_calc"], mode="markers", name = "CDW", marker= dict(color= col_dict["CDW"], size= 8, symbol="x"), legendgroup ="group3"), row=2, col=1)
     fig.add_trace(go.Scatter(x=dfg[exp]["simulated"].index, y=dfg[exp]["simulated"]["cX"], mode="lines", name = "CDW_simulated", line= dict(color= col_dict["CDW"], width= 2, dash="dot"), legendgroup ="group3"), row=2, col=1)
     fig.add_trace(go.Scatter(x=[dfg[exp]["off"].index], y= dfg[exp]["off"]["CDW_calc_rate"], mode="lines", name = "Growthrate", line= dict(color= col_dict["mu"], width= 3, dash="solid"),legendgroup ="group3"),row=2, col=1, secondary_y=True)
     fig.update_yaxes(title_text="CDW [g/L]", row=2,col=1, secondary_y =False, titlefont=dict(size=15), color=col_dict["CDW"],showline=True, linewidth=2, linecolor=col_dict["CDW"],showgrid=True, gridwidth=1, gridcolor="rgb(240,240,240)",zeroline=False,range = [0,30])
     fig.update_yaxes(title_text="Growthrate [g/L/h]", row=2,col=1, secondary_y =True, titlefont=dict(size=15), color=col_dict["mu"],showline=True, linewidth=2, linecolor=col_dict["mu"],showgrid=False, gridwidth=1, gridcolor='LightPink',zeroline=False,range = [-2,2])
        
        #fig 2nd row, right: pO2 and pO2 consumption rate 
     fig.add_trace(go.Scatter(x=dfg[exp]["on"].index, y= dfg[exp]["on"][("pO2_2","Value")], mode="markers", name = "first", marker= dict(color= col_dict[("pO2_2","Value")], size= 3, symbol="0"), legendgroup ="group4"), row=2, col=2)
     fig.add_trace(go.Scatter(x=dfg[exp]["on"].index, y=dfg[exp]["on"][("FOAMT_2","Value")], mode="markers", name = "Antifoam", marker= dict(color= col_dict[("FOAMT_2" , "Value")], size= 2, symbol="0"),legendgroup ="group6"),row=2, col=2, secondary_y=True)
       #fig.add_trace(go.Scatter(x=[0], y= [0], mode="markers", name = "Oxygen consumption rate", marker= dict(color= col_dict["rO2"], size= 5, symbol="0"),legendgroup ="group4"),row=2, col=2, secondary_y=True)
     fig.update_yaxes(title_text="DO[%]", row=2,col=2, secondary_y =False, titlefont=dict(size=15), color=col_dict[("pO2_2","Value")],showline=True, linewidth=2, linecolor=col_dict[("pO2_2","Value")],showgrid=True, gridwidth=1, gridcolor="rgb(240,240,240)",zeroline=False,range = [0,110])
     fig.update_yaxes(title_text="Antifoam [mL]", row=2,col=2, secondary_y =True, titlefont=dict(size=15), color=col_dict[("FOAMT_2" , "Value")],showline=True, linewidth=2, linecolor=col_dict[("FOAMT_2" , "Value")],showgrid=False, gridwidth=1, gridcolor='LightPink',zeroline=False,range = [0,3])
        
        #fig 3d row, left: Glc, Glc consumption rate, Feedrate
        
     fig.add_trace(go.Scatter(x=dfg[exp]["off"].index, y= dfg[exp]["off"]["Glucose [g/L]"], mode="markers", name = "Glucose", marker= dict(color= col_dict["Glc [g/L]"], size= 8, symbol="x"), legendgroup ="group5"), row=3, col=1)
     fig.add_trace(go.Scatter(x=dfg[exp]["off"].index, y= dfg[exp]["off"]["Glucose [g/L]_rate"], mode="lines", name = "Glc consumption rate", line= dict(color= col_dict["rGlc"], width= 1, dash="solid"),legendgroup ="group5"),row=3, col=1, secondary_y=True)
        #fig.add_trace(go.Scatter(x=[0], y= [0], mode="markers", name = "Glc Feed rate", marker= dict(color= col_dict["GlcFeed_rate"], size= 5, symbol="0"),legendgroup ="group5"),row=3, col=1, secondary_y=True)
     fig.add_trace(go.Scatter(x=dfg[exp]["simulated"].index, y=dfg[exp]["simulated"]["cS"], mode="lines", name = "Glucose_simulated", line= dict(color= col_dict["Glc [g/L]"], width= 2, dash="dot"), legendgroup ="group5"), row=3, col=1)
     fig.update_yaxes(title_text="Glucose [g/L]", row=3,col=1, secondary_y =False, titlefont=dict(size=15), color=col_dict["Glc [g/L]"],showline=True, linewidth=2, linecolor=col_dict["Glc [g/L]"],showgrid=True, gridwidth=1, gridcolor="rgb(240,240,240)",zeroline=False,range = [0,70])
     fig.update_yaxes(title_text="Glc consumption rate [g/L/h]", row=3,col=1, secondary_y =True, titlefont=dict(size=15), color=col_dict["rGlc"],showline=True, linewidth=2, linecolor=col_dict["rGlc"],showgrid=False, gridwidth=1, gridcolor='LightPink',zeroline=False,range = [-25,25])
       
        #fig 3d row, right: Volume, Fout
        
     fig.add_trace(go.Scatter(x=dfg[exp]["simulated"].index, y=dfg[exp]["simulated"]["V"], mode="lines", name = "Volume_simulated", line= dict(color= col_dict["V"], width= 2, dash = "dot"), legendgroup ="group6"), row=3, col=2,secondary_y=True)
       #fig.add_trace(go.Scatter(x=dfg[exp]["on"].index, y=dfg[exp]["on"][("FOAMT_2","Value")], mode="markers", name = "Antifoam", marker= dict(color= col_dict[("FOAMT_2" , "Value")], size= 2, symbol="0"),legendgroup ="group6"),row=3, col=2, secondary_y=False)
     fig.add_trace(go.Scatter(x=dfg[exp]["on"].index, y=dfg[exp]["on"][("SUBS_A2","Value")], mode="markers", name = "FeedPumpPower", marker= dict(color=col_dict["V"], size= 3, symbol="0"), legendgroup ="group7"), row=3, col=2,secondary_y=False)
     fig.add_trace(go.Scatter(x=dfg[exp]["simulated"].index, y=dfg[exp]["simulated"]["Fout"],  mode="lines", name = "Outflow from Sampling", line= dict(color= col_dict["Baserate"], width= 1, dash="solid"),legendgroup ="group6"),row=3, col=2, secondary_y=True)
     fig.update_yaxes(title_text="Feed Pump Power[%]", row=3, col=2, secondary_y =False, titlefont=dict(size=15), color=col_dict["F_out"],showline=True, linewidth=2, linecolor=col_dict["F_out"],showgrid=True, gridwidth=1, gridcolor="rgb(240,240,240)",zeroline=False,range = [0,10])
     fig.update_yaxes(title_text="Volume [L], Fout [L/h]", row=3, col=2, secondary_y =True, titlefont=dict(size=15), color=col_dict["V"],showline=True, linewidth=2, linecolor=col_dict["V"],showgrid=False, gridwidth=1, gridcolor='LightPink',zeroline=False,range = [0,1.5])
        
        #fig 4th row,left: Base demand(cummulated), Base consumption rate, pH
     fig.add_trace(go.Scatter(x=dfg[exp]["on"].index, y=dfg[exp]["on"][("BASET_2","Value")], mode="markers", name = "cumulated NH3 demand", marker= dict(color= col_dict[("BASET_2","Value")], size= 3, symbol="0"), legendgroup ="group7"), row=4, col=1)
     fig.add_trace(go.Scatter(x=dfg[exp]["on"].index, y=dfg[exp]["on"]["('BASET_2', 'Value')_rate"], mode="lines", name = "Base consumption rate", line= dict(color= col_dict["Baserate"], width= 1, dash="solid"),legendgroup ="group7"),row=4, col=1, secondary_y=True)
     #fig.add_trace(go.Scatter(x=dfg[exp]["on"].index, y= dfg[exp]["on"]["pH_2"], mode="markers", name = "pH", marker= dict(color= col_dict["pH"], size= 5, symbol="0"),legendgroup ="group7"),row=4, col=1, secondary_y=False)
     fig.update_yaxes(title_text="Base demand [mL]", row=4,col=1, secondary_y =False, titlefont=dict(size=15), color=col_dict[("BASET_2","Value")],showline=True, linewidth=2, linecolor=col_dict[("BASET_2","Value")],showgrid=True, gridwidth=1, gridcolor="rgb(240,240,240)",zeroline=False, range =[0,100])
     fig.update_yaxes(title_text="Base consumption rate[mL/h]", row=4,col=1, secondary_y =True, titlefont=dict(size=15), color=col_dict["Baserate"],showline=True, linewidth=2, linecolor=col_dict["Baserate"],showgrid=False, gridwidth=1, gridcolor='LightPink',zeroline=False,range = [0,50])
        
        #fig 4th row, right: Acetate, unknown Peaks 
     fig.add_trace(go.Scatter(x=dfg[exp]["off"].index, y=dfg[exp]["off"]["Acetate [g/L]"], mode="markers", name = "Acetate", marker= dict(color= col_dict["Acetate"], size= 8, symbol="x"), legendgroup ="group8"), row=4, col=2)
        
     if "Peak 22.87 min" in dfg[exp]["off"]:
         fig.add_trace(go.Scatter(x=dfg[exp]["off"].index, y= dfg[exp]["off"]["Peak 22.87 min"], mode="markers", name = "Peak at 22.87 min", marker= dict(color= col_dict["Peak 22.84"], size= 5, symbol="x"),legendgroup ="group8"),row=4, col=2, secondary_y=True)
     else:
         fig.add_trace(go.Scatter(x=[], y= [], mode="markers", name = "-", marker= dict(color= col_dict["Peak 22.84"], size= 5, symbol="x"),legendgroup ="group8"),row=4, col=2, secondary_y=True)
         
    
         
     if "Peak 24.21 min" in dfg[exp]["off"]:
         fig.add_trace(go.Scatter(x=dfg[exp]["off"].index, y= dfg[exp]["off"]["Peak 24.21 min"], mode="markers", name = "Peak at 24.21 min", marker= dict(color= col_dict["Peak 24.21"], size= 5, symbol="x"),legendgroup ="group8"),row=4, col=2, secondary_y=True)
     else:
         fig.add_trace(go.Scatter(x=[], y= [], mode="markers", name = "-", marker= dict(color= col_dict["Peak 24.21"], size= 5, symbol="x"),legendgroup ="group8"),row=4, col=2, secondary_y=True) 
            
     fig.update_yaxes(title_text="Acetate [g/L]", row=4,col=2, secondary_y =False, titlefont=dict(size=15), color=col_dict["Acetate"],showline=True, linewidth=2, linecolor=col_dict["Acetate"],showgrid=True, gridwidth=1, gridcolor="rgb(240,240,240)",zeroline=True,range = [0,15])
     fig.update_yaxes(title_text="Area of unknown Peaks [-]", row=4,col=2, secondary_y =True, titlefont=dict(size=15), color=col_dict["Peak 24.21"],showline=True, linewidth=2, linecolor=col_dict["Peak 24.21"],showgrid=False, gridwidth=1, gridcolor='LightPink',zeroline=False,range = [0,200000])
        
        
     fig.update_xaxes(title_text = "F-time [h]",matches = "x",showline=True, linewidth=2, linecolor='black',showgrid=True, gridwidth=1, gridcolor="rgb(240,240,240)",zeroline=False, range=[0,70])
     fig.update_layout(title_text = exp , title_x=0.5, title_y= 0.99, titlefont=dict(size=30))

     fig.update_layout(legend=dict(orientation="h",
                                   yanchor="bottom",
                                   y=1.01,
                                   xanchor="center",
                                   x=0.5,
                                   title_font_family="Times New Roman"))
     fig.update_layout(autosize=True,
                       width=1000,
                       height=1500,
                       template='plotly_white',
                       margin=dict( l=50,
                                   r=50,
                                   b=10,
                                   t=200,
                                   pad=4),
                       paper_bgcolor="rgb(234,242,242)")
     
     fig.show()