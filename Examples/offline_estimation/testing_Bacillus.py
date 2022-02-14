import sys
sys.path.append("V:/biomoni/")                      #without this line windows could not acces the package 
from biomoni.Experiment import Experiment
from biomoni.BacillusVariableFeedrate import Bacillus_vf
from biomoni.visualize import visualize
import pandas as pd
import numpy as np


path = r"V:/biomoni/BacillusDaten/Stamm185"

experiment_dict = {exp : Experiment(path, exp) for exp in ["F1", "F2", "F3", "F5" ]}  #all experiments in a dictionary
experiment_list = list(experiment_dict.values())
b= Bacillus_vf()
a = Experiment(path, "F4")

#y.simulate(a)
b.estimate(experiment_dict)
b.report()
t_grid = np.linspace(0,64,101) 
sim_dict_all = {experiment.exp_id: Bacillus_vf.simulate(experiment = experiment, t_grid = t_grid) for experiment in experiment_dict.values()}
for exp_id  in experiment_dict.keys():
    title = "Experiment {0}".format(exp_id) 
    visualize(experiment_dict[exp_id] , sim_dict_all[exp_id], title = title, suffix_1= "_experimental", suffix_2 = "_fitted")