import os
import pathlib
import time

from biomoni.Experiment import Experiment
from biomoni.Yeast import Yeast
from biomoni.visualize import visualize
from settings import Result_path, data_name
from param_collection import p1


sub_paths = next(os.walk(Result_path))[1]       #yields the subsirectory in the given path
newest_results_dir = max([os.path.join(Result_path,i) for i in sub_paths], key=os.path.getmtime) #gives newest subdirectory


path = Result_path

exp_id = "current_ferm"     #identifier to read correct data from the metadata
typ1 = "on_CO2"     #type of data
meta_path = "metadata_OPCUA.ods"    #metadata location within path
types = {typ1 : data_name}
index_ts = {typ1 : 0}
read_csv_settings = {typ1 : dict(sep=";",encoding= "unicode_escape",decimal=",", skiprows=[1,2] , skipfooter=1, usecols = None, engine="python")}
to_datetime_settings = {typ1 : dict(format = "%d.%m.%Y  %H:%M:%S", exact= False, errors = "coerce") }
calc_rate =(typ1, "BASET")
exp_dir_manual = pathlib.PurePath(newest_results_dir).name  #this is given because the subfoldername does not match the 
exp_dir_manual = newest_results_dir
read_excel_settings = dict(engine = "odf")

while True:
    time.sleep(5*60)    #5 min sampling before each estimation

    Exp = Experiment(path = Result_path, meta_path= meta_path, exp_id = exp_id, types = types, index_ts= index_ts
    , read_csv_settings= read_csv_settings, to_datetime_settings= to_datetime_settings, calc_rate= calc_rate, exp_dir_manual = exp_dir_manual, read_excel_settings = read_excel_settings)


    y = Yeast()
    y.set_params(p1)    #change initial parameters to test estimation
    y.estimate(Exp, tau = 1, max_nfev = 100)     #max function evaluations #, max_nfev = 100  



    y.report()
    print(y.stat_all)
    print(y.stat_single)
    
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

    



