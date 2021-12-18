import path
import os
import pathlib
import time
import sys

import_path = "/home/paul/Desktop/pCloudDrive/param/BioMoni"
sys.path.append(import_path)        #to get modules from another directory

from Experiment import Experiment
from Yeast import Yeast
from settings_opcua import RESULT_PATH, data_name
from Yeast import Yeast
from visualize import visualize
from param_collection import p0, p1

#RESULT_PATH = r"P:\Code\Model_Martin_Version_ 08.09.2021\Testing_opcua\Results"    #MANUAL, DONT FORGET WEG

sub_paths = next(os.walk(RESULT_PATH))[1]       #yields the subsirectory in the given path
newest_results_folder = max([os.path.join(RESULT_PATH,i) for i in sub_paths], key=os.path.getmtime) #gives newest subdirectory
newest_results_folder

path = RESULT_PATH

typ1 = "on_CO2"
meta_path = "metadata_OPCUA.ods"



exp_id = "current_ferm"     #identifier to read correct data from the metadata



types = {typ1 : data_name}
index_ts = {typ1 : 0}
read_csv_settings = {typ1 : dict(sep=";",encoding= "unicode_escape",decimal=",", skiprows=[1,2] , skipfooter=1, usecols = None, engine="python")}
to_datetime_settings = {typ1 : dict(format = "%d.%m.%Y  %H:%M:%S", exact= False, errors = "coerce") }
calc_rate =(typ1, "BASET")
exp_dir_manual = pathlib.PurePath(newest_results_folder).name  #this is given because the subfoldername does not match the 
exp_dir_manual = newest_results_folder
read_excel_settings = dict(engine = "odf")

while True:
    #man könnte erst time sleep 10 min machen
    #2 verschiedene Fehler wenn man sofort nach beginn des datenerstellens sampled: einmal  File "p:\BioMoni\Experiment.py", line 225, in _calc_t df["t"] = (df["ts"] - df["ts"][0]) / pd.Timedelta(1,"h"), IndexError: index 0 is out of bounds for axis 0 with size 0
    #und einmal ein paar Messwerte später "p:\BioMoni\Model.py", line 103, in estimate,  result = minimize(self.residuals, p, args=(datasets_dict,  settings_dict, kwargs_solve_ivp), method = method, nan_policy = nan_policy, **fit_kws), \scipy\optimize\minpack.py", TypeError: Improper input: N=4 must not exceed M=3
    time.sleep(1)

    Exp = Experiment(path = RESULT_PATH, meta_path= meta_path, exp_id = exp_id, types = types, index_ts= index_ts
    , read_csv_settings= read_csv_settings, to_datetime_settings= to_datetime_settings, calc_rate= calc_rate, exp_dir_manual = exp_dir_manual, read_excel_settings = read_excel_settings)


    y = Yeast()
    y.set_params(p1)    #change initial parameters to test estimation
    y.estimate(Exp, tau = 1, max_nfev = 100)     #max function evaluations #, max_nfev = 100  



    y.report()
    print(y.stat_all)
    print(y.stat_single)
    #y.set_params()      #wegen online schätzung -> bessere default values
    print(y.p_default)
    print(y.p)

    sim = y.simulate(Exp)
    print(Exp.dataset)
    print(Exp.metadata)

    print("qsmax initial : ", p1["qsmax"])
    print("base_coef initial : ", p1["base_coef"])

    print("qsmax : ", y.p["qsmax"])
    print("base_coef : ", y.p["base_coef"])   

    visualize(Exp, sim, column_dict = {"BASET_rate" : "cyan", "CO2" : "orange"}, secondary_y_cols = ["CO2"])

    



