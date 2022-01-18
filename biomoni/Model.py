from biomoni import Experiment
import time
from lmfit import Parameters, report_fit, minimize
import pandas as pd
import numpy as np

class Model:
    """
    :class Yeast: A class to create a basic model, able to estimate.
    
    """



    def __repr__(self):
        """Representation of the Model object in the print() call"""
    
        return  "Model()" 


    def __init__(self):
        """Initialization of an object of the Model class"""
        
        pass

    
    def set_params(self, p = None):
        """Set the parameters of the object"""
        self.p = p

    def chemical_balancing(self,p):
        """Calculate yield coefficients from several metabolic pathways """
        pass
    
    def create_settings(self, experiment):
        """Creates settings required for estimation. """

        pass
    
    def change_params(self, name, value = None, vary = None, min = None, max = None):    
        """Function to change existing parameter values
        :param name: name of the paramater to change 
        :type name: str 
        """
        
        try:  
            self.p[name]        #checking if variable exists
        except KeyError:
            raise NameError("There is no Parameter called {0}.".format(name))

            
        if value is not None:
            assert isinstance(value, (int, float)), "value must be int or float."
            self.p[name].value = value
        if vary is not None:
            assert isinstance(vary, bool), "vary must be True or False."
            self.p[name].vary = vary
        if min is not None:
            assert isinstance(min, (int, float)), "min must be int or float."
            self.p[name].min = min
        if max is not None:
            assert isinstance(max, (int, float)), "max must be int or float."
            self.p[name].max = max
        
        
        
        


    def estimate(self, experiments = None, method = "leastsq", nan_policy = "omit"
    ,  p_given = None, datasets_dict_given = None, settings_dict_given = None, tau = None
    ,kwargs_solve_ivp = dict(method= "Radau", first_step = 0.0000001, max_step= 0.1 )
    ,**fit_kws ):           
        """ Function which allows to estimate either directly from the experiments or from manually given data.

        :param experiments: Experiments from which the data are obtained, either single BioMoni.Experiment.Experiment object list or dict with Experiment objects
        :type experiments: BioMoni.Experiment.Experiment object, dict or list with Experiment objects
        :param method: Fitting method for lmfit.minimize 
        :type method: str, optional
        :param nan_policy: Decision on how to deal with nan values in the estimation, method from lmfit.minimize 
        :type nan_policy: str, optional
        :param p_given: Manually given structure with parameter values to be estimated, cf. lmfit.parameter.Parameters
        :type p_given: lmfit.parameter.Parameters object
        :param datasets_dict_given: Manually given dictionary with key: Experiment number, value: dataset with data from different measurement methods (e.g. off, on CO2 data)
        :type datasets_dict_given: dict

        :param settings_dict_given: Manually given dictionary with settings for the simulation for each given Experiment. Key: Experiment number, value: settings containing:

            * y0: Initial state vector 
            * c: Control variables
            * wf: Weighting factors

        :type settings_dict_given: dict
        :param tau: decay factor for a weighting factor decrase over time, required in the residuals function.
        :type tau: float, int or None
        :param kwargs_solve_ivp: Arguments for scipy.integrate.solve_ivp, in form of a dictionary
        :type kwargs_solve_ivp: dict, optional
        :param **fit_kws: optional Arguments for lmfit.minimize based on scipy.optimize
        :type **fit_kws: kwargs

        :return: Estimation outcomes (p and yields), estimation duration, simulated and measured data and residuals at the measurement time points, statistics (RMSE, BIAS, STDDEV).
        Everything is saved as object attribute.
    

        """
        
        if experiments is not None:

            self.prepare(experiments)

            datasets_dict = self.datasets_dict
            settings_dict = self.settings_dict

        
        else:   
            for i in [datasets_dict_given, settings_dict_given]:
                if i is None:
                    raise ValueError("If no Experiments are given, it is necessary to give datasets and corresponding settings")
        
        

        if p_given is not None:
            p = p_given
        else:
            p = self.p


        if datasets_dict_given is not None:
            datasets_dict = datasets_dict_given
        if settings_dict_given is not None:
            settings_dict = settings_dict_given


        start = time.time()
        result = minimize(self.residuals, p, args=(datasets_dict,  settings_dict, tau, kwargs_solve_ivp), method = method, nan_policy = nan_policy, **fit_kws)
        end = time.time()

        self.duration = end - start         
        self.result = result
        self.p = result.params
        self.yields = self.chemical_balancing(self.p)   #this may be not generic since the function chemical balancing may not be given in every submodel
        self.data_and_residuals = self.get_residuals(self.p, self.yields, datasets_dict, settings_dict, kwargs_solve_ivp)   #same prob with self.yields. there may be models which dont even need extra calculated yields in order to simulate 
        self.stat_single, self.stat_all = self.statistics(self.data_and_residuals)


    
    def report(self):
        """ Function to report fit results
        """
        report_fit(self.result)      


    def prepare(self, experiments):
        """Function that prepares the data for the estimation, datasets and corresponding settings are extracted from the experiment objects.

        :param experiments: Experiments from which the data are obtained, either single BioMoni.Experiment.Experiment object list or dict with Experiment objects.
        :type experiments: BioMoni.Experiment.Experiment object, dict or list with Experiment objects.

        :return: prepared dataset and settings saved as object attribute.

        """
        assert isinstance(experiments, (Experiment, dict, list)), "Given Experiments should either be of type dict or list with Experiment objects, or in case of a single Experiment of type Experiment"

        if type(experiments) is dict:
            assert all(isinstance(i, Experiment) for i in experiments.values()) , "Each value in experiments have to be of type Experiment" 
            experiment_list = list(experiments.values())
            experiment_dict = {experiment.exp_id : experiment for experiment in experiment_list}   

        elif type(experiments) is list:   
            assert all(isinstance(i, Experiment) for i in experiments) , "Each element in experiments have to be of type Experiment"    
            experiment_list = experiments
            experiment_dict = {experiment.exp_id : experiment for experiment in experiment_list}     

        elif type(experiments) is not list: 
            experiment_list = [experiments]    
            experiment_dict = {experiment.exp_id : experiment for experiment in experiment_list}   
        
        self.experiment_dict = experiment_dict

        self.settings_dict = {}
        self.datasets_dict = {}


        for exp_id, experiment in experiment_dict.items():
            self.datasets_dict[exp_id] = experiment.dataset
            self.settings_dict[exp_id] = self.create_settings(experiment)
            





    def get_residuals(self, p, yields, datasets_dict, settings_dict, kwargs_solve_ivp):
        """ Function that simulates on the measured data time points with given p and yields and calculates the residuals, The experimental data in datasets_dcit and settings must match.

        :param p: Structure with parameter values to be estimated, cf. lmfit.parameter.Parameters
        :type p: lmfit.parameter.Parameters object
        :param yields: Calculated yields coefficients
        :type yields: dict
        :param settings_dict_given: Manually given dictionary with settings for the simulation for each given Experiment. Key: Experiment number, value: settings containing:

            * y0: Initial state vector 
            * c: Control variables
            * wf: Weighting factors

        :type settings_dict: dict
        :param kwargs_solve_ivp: Arguments for scipy.integrate.solve_ivp, in form of a dictionary, should be the same as used in the parameter estimation
        :type kwargs_solve_ivp: dict 
        """
        
        

       
        residuals_dict_all = {}       #dict wich will contain residuals for all experiments 
        for exp_id, dataset in datasets_dict.items():
            settings = settings_dict[exp_id]
            y0 = settings["y0"] #extract y0,c and weighting factors for each experiment
            c = settings["c"]
            
            residuals_dict_single = {}     #dict which will 

            for dskey, dat in dataset.items():   #extract data for ("on", "off", "CO2") for each experiment
                     
                t_grid = dat.index.values           #extract time from dataframe
                sim_exp = self.simulate(None, t_grid, y0, p, c, yields, kwargs_solve_ivp)
                residuals_dict = {}     #dict which will contain residuals for one type of measurement ("on", "off", "CO2")

                for var in dat:     #loop over measured variables (columns)
                    if var in sim_exp.columns:
                        
                        res_var = (sim_exp[var] - dat[var]).values #residuals
                        residuals_dict[var] = dat[var].values
                        residuals_dict["{0}_fitted".format(var)] = sim_exp[var].values
                        residuals_dict["{0}_residuals".format(var)] = res_var   
                    else:
                        pass

                residuals_dict_single[dskey] = pd.DataFrame(residuals_dict, index = t_grid)     #residuals in form of pandas dataframe

            residuals_dict_all[exp_id] =residuals_dict_single
        
        return residuals_dict_all



    def statistics(self, residuals_dict_all):
        """ Function which calculates statistics (RMSE, BIAS, STTDEV) for each experiment individually and for all experiments together.

        :param residuals_dict_all: Dictionary containing residuals from experiments. Key: Experiment number, value: df with residuals.
        :type residuals_dict_all: dict
        :param data_types: type of measurement data e.g. off, on CO2 measurements
        :type data_types: list of strings

        :return: statistics_single_exp: statistics  for single each experiment, statistics_all_exp: statistics for all experiments simultaneously
        
        """

        #calculating statistics for every experiment
        statistics_single_exp = {}
        for exp_id, dataset in residuals_dict_all.items():
            RMSE = {}
            BIAS = {}
            STDDEV = {}
            stat = {}
            for dskey, df in dataset.items():
                for column in df:
                    if "residuals" in column:
                        name = column.split("_")[0]
                        RMSE[name] = np.sqrt(np.mean(df[column]**2))
                        BIAS[name] = np.mean(df[column])
                        STDDEV[name] = np.std(df[column])
            stat["RMSE"] = RMSE
            stat["BIAS"] = BIAS
            stat["STDDEV"] = STDDEV
            statistics_single_exp[exp_id] = stat


        #calculating statistics for all experiments simultaneously


        # first of all of we need every possible data type = dskey
        data_types_all = [residuals_dict_all[i].keys() for i in residuals_dict_all.keys()] # data types of each experiment
        data_types = max(data_types_all, key = len)    #data types with the highest length to be sure to get all different data types
        
        data = {}
        for i in data_types:
            data[i] = []

        for dataset in residuals_dict_all.values():
            for dskey, df in dataset.items():
                data[dskey].append(df)

        data_concat = {}
        for dskey in data.keys():
            data_concat[dskey] = pd.concat(data[dskey])


        statistics_all_exp = {}
        RMSE = {}
        BIAS = {}
        STDDEV = {}

        for df in data_concat.values():
            for column in df:
                if "residuals" in column:
                    name = column.split("_")[0]
                    RMSE[name] = np.sqrt(np.mean(df[column]**2))
                    BIAS[name] = np.mean(df[column])
                    STDDEV[name] = np.std(df[column])


            statistics_all_exp["RMSE"] = RMSE
            statistics_all_exp["BIAS"] = BIAS
            statistics_all_exp["STDDEV"] = STDDEV

        return statistics_single_exp, statistics_all_exp
