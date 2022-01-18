from biomoni.Experiment import Experiment
from biomoni.Yeast_variable_feedrate import Yeast_vf
import pandas as pd
import numpy as np

#Linux path
path = "/home/paul/pCloudDrive/Code/Messdaten" 

# #windows path
#path = r"P:\Code\biomoni\Messdaten"

experiment_dict = {exp : Experiment(path, exp) for exp in ["F4", "F5", "F6", "F7", "F8" ]}  #all experiments in a dictionary
experiment_dict["F8"].time_filter(dskey= "on", start = pd.to_datetime("14.12.2020  12:20:16")) #special time filter for experiment 8
[experiment_dict[exp].pop_dataframe("on") for exp in ["F4", "F5", "F6"]]   #delete online data in experiment 4,5,6 because of bad BASET_rate measurements
experiment_list = list(experiment_dict.values())   

a = Experiment(path, "F7")

y = Yeast_vf()
y.simulate(a)



y.estimate(experiment_dict)
y.report()