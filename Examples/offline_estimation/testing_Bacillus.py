import sys
sys.path.append("V:/biomoni/")                      #without this line windows could not acces the package 
from biomoni.Experiment import Experiment
from biomoni.BacillusVariableFeedrate import Bacillus_vf
import pandas as pd
import numpy as np


path = r"V:/biomoni/BacillusDaten/Stamm185"

experiment_dict = {exp : Experiment(path, exp) for exp in ["F1", "F2", "F3", "F4", "F5" ]}  #all experiments in a dictionary
experiment_list = list(experiment_dict.values())
y = Bacillus_vf()
a = Experiment(path, "F4")

y.simulate(a)
y.estimate(experiment_dict)
y.report()