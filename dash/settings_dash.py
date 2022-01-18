#Settings required to run the dash app


measurement_vars = ["CO2", "BASET_rate"]
simulated_vars = ["cX", "cS", "cE", "CO2", "BASET_rate"]

fitparams = ["qsmax", "base_coef", "qO2max", "Yxs_ox"]


kwargs_experiment = {       #setting to create Experiment object

"online_est" : dict (


exp_id = "current_ferm",     #identifier to read correct data from the metadata
meta_path = "metadata_OPCUA.ods",    #metadata location within path
types = {"on_CO2" : "data_1"},
index_ts = {"on_CO2" : 0},
read_csv_settings = {"on_CO2" : dict(sep=";",encoding= "unicode_escape",decimal=",", skiprows=[1,2] , skipfooter=1, usecols = None, engine="python")},
to_datetime_settings = {"on_CO2" : dict(format = "%d.%m.%Y  %H:%M:%S", exact= False, errors = "coerce") },
calc_rate =("on_CO2", "BASET"),
read_excel_settings = dict(engine = "odf")

)
}

kwargs_estimate = {         #settings to estimate

"online_est" : dict (
tau = 1, 
max_nfev = 100      #max function evaluations
)


}
