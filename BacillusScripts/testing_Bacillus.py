import sys
sys.path.append("V:/biomoni/")                      #without this line windows could not acces the package 
from BacillusScripts.Experiment import Experiment
from BacillusScripts.BacillusVariableFeedrate import Bacillus_vf
from BacillusScripts.visualizationBacillus import visualizeBacillusFermentation
import pandas as pd
import numpy as np
from IPython.display import display

path = r"V:/biomoni/BacillusData/Stamm185"

experiment_dict_for_estimation = {exp : Experiment(path, exp,endpoint = "F_end") for exp in ["F1", "F2", "F3", "F4", "F5" ]}  #all experiments in a dictionary   
experiment_dict_for_graphs = {exp : Experiment(path, exp, endpoint = "F_end") for exp in ["F1", "F2", "F3", "F4", "F5" ]} 


#Exp = experiment_dict 
b= Bacillus_vf()
b.estimate(experiment_dict_for_estimation)
b.report()

t_grid = np.linspace(0,65,101) 
sim_dict_all = {experiment.exp_id: b.simulate(experiment = experiment, t_grid = t_grid) for experiment in experiment_dict_for_estimation.values()} 

sim_list=list(sim_dict_all.values())
experiment_list=list(experiment_dict_for_graphs.values())
# put experimental data in dictionary in form of {"F1":{"off": ---DataFrame---,"on": ---DataFrame---,"CO2": ---DataFrame--- |etc.}}
ex ={}

for obj in experiment_list:
    ex[obj.exp_id]= obj.dataset 
        
# put simulated data in dictionary in form of {"F1":{"simulated": ---DataFrame--- |etc}}

sim_list=list(sim_dict_all.values())
sim ={}
for key in experiment_dict_for_graphs.keys():
    sim[key] = {}
 
sim_complete = {list(sim)[i]: {"simulated":sim_list[i]}for i in range(len(sim_list))}
   
# bring both dictionaries together in form of {"F1":{"simulated" ---DataFrame---,"off": ---DataFrame---,"on": ---DataFrame---,"CO2": ---DataFrame--- |etc.}}

dict_for_graphs = {}
dict_for_graphs = sim_complete
dfg = dict_for_graphs  
for d in dict_for_graphs.keys():
    dict_for_graphs[d].update(ex[d])         

visualizeBacillusFermentation(dfg)