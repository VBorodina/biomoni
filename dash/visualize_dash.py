import pandas as pd
from plotly.subplots import make_subplots
import plotly.graph_objs as go
from biomoni.Experiment import Experiment

def visualize(data_1, data_2 = None, title = None
, column_dict = {"BASET_rate" : "cyan", "cX" : "red", "cS" : "green", "cE" : "blue", "CO2" : "orange"}
, suffix_1 = "", suffix_2 = "_fitted", mode_1 = "markers", mode_2 = "lines"
, secondary_y_cols = ["CO2", "BASET_rate"], yaxis_type = "linear", sec_yaxis_type = "linear"
): 
    """ Function to visualize experimental data and simulated data of one experiment either alone or in conjunction. It is necessary that the experimental data and the simulated data are from the same experiment to yield meaningful results. 
    The given data must be either a single pd.DataFrame: df or a list or a dict with pd.DataFrame objects, for example a dict with "off" : df1, "on" : df2 and "CO2" : df3. 


    :param data_1: Data to plot
    :type data_1: dict, list, pd.DataFrame
    :param data_2: Data to plot in conjunction with data_1
    :type data_2: dict, list, pd.DataFrame, optional
    :param title: title of plot, can be str, int (for experiment number) or None
    :type title: int, str, None, optional
    :param suffix_1: name appendix after column name in the legend for data_1
    :type suffix_1: str, optional
    :param suffix_2: name appendix after column name in the legend for data_2
    :type suffix_2: str, optional
    :param mode_1: type of marker (e.g. “lines”, “markers” or “lines+markers” for data_1
    :type mode_1: str, optional
    :param mode_2: type of marker (e.g. “lines”, “markers” or “lines+markers” for data_2
    :type mode_2: str, optional
    :param secondary_y_cols: list of column names which should be displayed on the secondary y axis
    :type secondary_y_cols: list, optional

    :return: Graph with plotted and named data


"""

    if secondary_y_cols is None:
        secondary_y_cols = []


    suffix = {1: suffix_1, 2: suffix_2}
    mode = {1: mode_1, 2: mode_2}

    data_both = {1: data_1, 2: data_2}
    data_dict ={}

    for i, data in data_both.items():
        assert isinstance(data, (pd.DataFrame, dict, list, Experiment, type(None))), "Data should be given either as an Experiment, a pd.DataFrame or as a list or dict with pd.DataFrame objects"

        if data is None:
            data_dict[i] = data
        
        elif type(data) is Experiment:
            data_dict[i] = list(data.dataset.values())


        elif type(data) is dict:
            assert all(isinstance(i, pd.DataFrame ) for i in data.values()) , "Each value in data_{0} have to be of type pd.DataFrame".format(i)
            data_dict[i] = list(data.values())
                

        elif type(data) is list:
            assert all(isinstance(i, pd.DataFrame ) for i in data) , "Each element in data_{0} have to be of type pd.DataFrame".format(i)
            data_dict[i] = data
                
        elif type(data) is pd.DataFrame: 
            data_dict[i] = [data]




    fig = make_subplots(specs = [[{"secondary_y" : True}]])

    for i, data_list in data_dict.items():
        if data_list is None:
            pass
    
        else:
            for df in data_list:
                
                for col in df:

                    secondary_y_flag = col in secondary_y_cols
                    if col in column_dict.keys():

                        fig.add_trace(
                            go.Scatter(x= df.index, y= df[col], name = col+ suffix[i], mode = mode[i], marker = dict(color = column_dict[col], size = 5, symbol = "x")
                            )
                            , secondary_y= secondary_y_flag
                        )
    

    cols_y1 = [col for col in column_dict.keys() if col not in secondary_y_cols]            
    cols_y2 = [col for col in column_dict.keys() if col in secondary_y_cols]
    fig.update_yaxes(title_text= str(cols_y1), secondary_y=False , title_standoff = 20 , type = yaxis_type)
    fig.update_yaxes(title_text= str(cols_y2), secondary_y=True, title_standoff = 20, type = sec_yaxis_type)

    fig.update_xaxes(title_text = "Time [h]")
    fig.update_layout(title_text = title , title_x=0.5)
    fig.update_layout(legend=dict(x = 1))
    #fig.show()

    return fig

