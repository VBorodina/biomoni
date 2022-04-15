import sys
sys.path.append("V:/biomoni/")                      #without this line windows could not acces the package 
from biomoni import Experiment
from BacillusScripts.BacillusVariableFeedrate import Bacillus_vf
from BacillusScripts.visualizationBacillus import visualizeBacillusFermentation
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from IPython.display import display

path = r"V:/biomoni/BacillusData/Stamm186"

experiment_dict_for_estimation = {exp : Experiment(path, exp,endpoint = "Sim_end") for exp in ["F6","F7"]}  #all experiments in a dictionary   
experiment_dict_for_graphs = {exp : Experiment(path, exp, endpoint = "F_end") for exp in ["F6","F7"]} 



#Exp = experiment_dict 
b= Bacillus_vf()
b.estimate(experiment_dict_for_estimation)
b.report()
print(b.p)

for key in b.stat_single.keys():
    df = pd.DataFrame(b.stat_single[key].values(), index = b.stat_single[key].keys())
    df_t=df.transpose()
    df_t["ID"]= key
    df_t["type"] = "single"
    print(df_t)
    


df_stat_all = pd.DataFrame(b.stat_all.values(), index = b.stat_all.keys())
df_stat_all_t= df_stat_all.transpose()
df_stat_all_t["ID"]= path
print(df_stat_all_t)



sim_dict_all= {}
for exp in experiment_dict_for_estimation.values():
    sim_end_h=exp.dataset["on"].index[-1]
    t_grid = np.linspace(0,sim_end_h,1001)
    sim_dict=  b.simulate(experiment = exp, t_grid = t_grid)
    sim_dict_all[exp.exp_id] = sim_dict



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