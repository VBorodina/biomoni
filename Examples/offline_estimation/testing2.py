from biomoni.Experiment import Experiment
from biomoni.Yeast_variable_feedrate import Yeast_vf
import pandas as pd
import numpy as np

#Linux path
path = "/home/paul/pCloudDrive/Code/Messdaten" 

Exp = Experiment(path, "F7")

print(Exp.path)