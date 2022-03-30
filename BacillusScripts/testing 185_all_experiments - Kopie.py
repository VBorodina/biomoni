import sys
sys.path.append("V:/biomoni/")                      #without this line windows could not acces the package 
from biomoni import Experiment
from BacillusScripts.BacillusVariableFeedrate_copy import Bacillus_vf
from BacillusScripts.visualizationBacillus import visualizeBacillusFermentation
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from IPython.display import display

path = r"V:/biomoni/BacillusData/Stamm185"

experiment_dict_for_estimation = {exp : Experiment(path, exp,endpoint = "Sim_end") for exp in ["F1","F2","F3","F4","F5",]}  #all experiments in a dictionary   
experiment_dict_for_graphs = {exp : Experiment(path, exp, endpoint = "F_end") for exp in ["F1","F2","F3","F4","F5"]} 



#Exp = experiment_dict 
b= Bacillus_vf()
b.estimate(experiment_dict_for_estimation)
b.report()
print(b.p)
print(b.stat_single,b.stat_all)
t_grid = np.linspace(0,40,1001) 
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
   
#calculate Carbon recovery 
CRR_all ={}

for key in experiment_dict_for_graphs.keys():
    CRR_all[key]= {}

for key, dat in sim_dict_all.items():
    V_func= interp1d(x=sim_dict_all[key]["V"].index, y=sim_dict_all[key]["V"], fill_value = (sim_dict_all[key]["V"].iloc[0], sim_dict_all[key]["V"].iloc[-1]) , bounds_error= False)
    exp =Experiment(path, key, endpoint="F_end")
    c = b.create_controls(exp)
    y0 = b.create_y0(exp)
    #datasaet = exp.dataset
    
    CRR1, CRR2 = b.calc_CRR(experiment=exp,c=c,V=V_func,y0=y0)
    
    CRR_all[key]= {"CR1": CRR1,"CR2": CRR2}


   
# bring both dictionaries together in form of {"F1":{"simulated" ---DataFrame---,"off": ---DataFrame---,"on": ---DataFrame---,"CO2": ---DataFrame--- |etc.}}

dict_for_graphs = {}
dict_for_graphs = sim_complete
dfg = dict_for_graphs  
for id in dict_for_graphs.keys():
    dict_for_graphs[id].update(ex[id])
    dict_for_graphs[id].update(CRR_all[id])         
  

visualizeBacillusFermentation(dfg)