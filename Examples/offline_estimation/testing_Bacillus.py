import sys
sys.path.append("V:/biomoni/")
from biomoni.Experiment import Experiment
from biomoni.Yeast_variable_feedrate import Yeast_vf
import pandas as pd
import numpy as np


path = r"V:/biomoni/BacillusDaten/Stamm185"

experiment_dict = {exp : Experiment(path, exp) for exp in ["F1", "F2", "F3", "F4", "F5" ]}  #all experiments in a dictionary

y.estimate(experiment_dict)
y.report()