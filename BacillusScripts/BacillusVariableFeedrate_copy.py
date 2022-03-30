import sys
import xxlimited
sys.path.append("V:/biomoni/")  
from biomoni import Experiment
from biomoni import Model
import pandas as pd
import numpy as np
import scipy


from scipy.integrate import solve_ivp

import sympy as sp
from sympy import Eq
from sympy import symbols

from copy import deepcopy

from scipy.interpolate import interp1d

import time #in Model
from lmfit import Parameters, report_fit, minimize  #in Model




from biomoni import Experiment
from biomoni import Model
import pandas as pd
import numpy as np


from scipy.integrate import solve_ivp

import sympy as sp
from sympy import Eq
from sympy import symbols

from copy import deepcopy

from scipy.interpolate import interp1d

import time #in Model
from lmfit import Parameters, report_fit, minimize  #in Model




class Bacillus_vf(Model):     #Dependent on base class
    """
    :class Bacillus: 
    
    """


    def __repr__(self):
        """Representation of the Bacillus object in the print() call"""
        
        return  "Bacillus object: Bacillus_vf()" 

    def __init__(self):
        """Initialization of an object of the Bacillus_vf(Model) class

        """
        
        Model.__init__(self)        #in case if Base class has something generic for all classes
        self.set_params()
        self.p_default = deepcopy(self.p)       #deepcopy maybe not necessary
        self.yields_default = deepcopy(self.yields)

        


    def set_params(self, p = None):
        """Set the parameters of the object, if no parameters are given, the default parameters are set.

        :param p: Fixed parameters and fit parameters
        :type p: lmfit.parameter.Parameters object

        :return: Set self.p to given p or default p and execute chemical_balancing to calculate yields from p
        
        """
        if p is None:          

            p = Parameters()
           # p.add("qsmax", value= 1.61290497, min=0.01, max=5.0 , vary = False)     #Maximum glucose uptake rate [g/(g*h)]   #previous start val: 0.1  #estimated value: 1.61290497
           # p.add("rP", value=0.1, min=0.01, max=10, vary=True)    #Production rate [g/(g*h)]     

           # p.add("base_coef", value= 0.00733798 , min=0.0001, max =3, vary = True)     #base coefficient [mol/g]     #previous start val: 1 #estimated value: 0.00733798      
           # p.add("qO2max", value= 0.16473635, min=0.1, max=0.4, vary = True)    #Maximum oxygen uptake rate [g/(g*h)]  , min=0.1  0.1, 0.4  #previous start val: 0.255984, 1 #estimated value: 0.17451356
           # p.add("qm_max", value=0.01, min=0.0075, max=0.0125, vary = False)  #glucose uptake rate required for maintenance [g/(g*h)]  
            
            p.add("Km", value=0.1, min=0.01, max=0.5,  vary=False)   #Saturation constant, concentration of glucose at µ = 0.5 µmax [g/L]
           # p.add("Ke", value=0.1, min=0.01, max=1, vary=False)    #Saturation constant, concentration of ethanol at µ = 0.5 µmax [g/L]
           # p.add("Ki", value=0.1, min=0.01, max=1, vary=False)    #Inhibition constant, glucose inhibits uptake of ethanol [g/L]

           
            p.add("Yxs", value=0.1, min=0.001, max=0.5, vary=True)    #Yield biomass per glucose [g/g]
            p.add("Ypx", value= 100, min= 1, max= 250, vary= True)   #Yield RF per biomass [mg/g]
            #p.add("Yxsmain",value = 0.2, min= 0.00001, max=1, vary=True)
            
            #For single experiment fitting:
            #p.add("viab_f", value= 0.0009, min = 0.00001, max=1, vary=True)
            

            #Biomass composition paramaters, content H,O,N from paper sonnleitner

            p.add("HX", value=1.594, min=1.77, max=2.1, vary=False) #Stoichiometric hydrogen content of biomass [mol/mol]
            p.add("OX", value=0.387, min=0.54, max=0.63, vary=False) #Stoichiometric oxygen content of biomass [mol/mol]     
            p.add("NX", value=0.239, min=0.14, max=0.16, vary=False) #Stoichiometric nitrogen content of biomass [mol/mol]
            
            #Biomassgrowth
            p.add("mu_max", value= 0.2, min= 0.1, max= 1, vary = True)               #EINHEIT?! g/h?

            self.p = p
            
        else:

            assert isinstance(p, type(Parameters())), "Given parameters must be of type lmfit.parameter.Parameters"
            self.p = p
            
        self.yields = self.chemical_balancing(self.p)
        
       
   


    def create_settings(self, experiment):
        """Creates settings required for estimation.

        :param experiment: Experiment object 
        :type experiment: BioMoni.Experiment.Experiment
        :return: settings: Settings for single experiment as dictionary.
        
        """

        #create weighting factors
        wf = self.calc_weighting_factors(experiment)

        #create control variables
        c = self.create_controls(experiment)

        #create y0: Initial state vector
        y0 = self.create_y0(experiment)

        settings = { "wf": wf , "c" : c , "y0": y0 }

        return settings
        
   

    def calc_mean_pressure(self, experiment):
        """Calculates the CO2 mean pressure in the BlueSens sensor for all time points.
        
        :param experiment: Experiment object 
        :type experiment: BioMoni.Experiment.Experiment
        :return: A single value: CO2 mean pressure

        
        """
        if "CO2" in experiment.dataset.keys():      #not clean coded, assumes that data has the key CO2 and the column p in order to calculate pressure
            pressure = experiment.dataset["CO2"]["p"].mean()
        else:
            pressure = 1

        return pressure

    def calc_weighting_factors(self, experiment):
        """Calculates weighting factors from an single experiment with respect to the magnitude of values
         e.g. wf(i) = len("off") dataset / len(i) dataset for i = "CO2" or "on".
        
        :param experiment: Experiment object
        :type experiment: BioMoni.Experiment.Experiment
        :return: weighting factors for single experiment as dictionary
        """

        df_smallest_len = min(experiment.dataset.values() , key = len)   #return dataframe with smallest length
        weighting_factors = {}
        for dskey, df in experiment.dataset.items():
            weighting_factors[dskey] = {}
            for colname in df:
               x = ["Glucose [g/L]","RF [mg/L]","Acetat [g/L]","CDW_calc",("BASET_2","Value"),('pO2_2', 'Value'),"CO2"]
               
               if colname in x:
                   weighting_factors[dskey][colname] = (len(df_smallest_len) / len(experiment.dataset[dskey]))/experiment.dataset[dskey][colname].median()
                   #if colname == "Glucose [g/L]":
                       #weighting_factors[dskey][colname] = (len(df_smallest_len) / len(experiment.dataset[dskey])) 
                   #if colname == "RF [mg/L]":
                      # weighting_factors[dskey][colname] = (len(df_smallest_len) / len(experiment.dataset[dskey]))/1000
                   #if colname == "CDW_calc":
                       #weighting_factors[dskey][colname] = (len(df_smallest_len) / len(experiment.dataset[dskey]))
                   
               
        print(weighting_factors)
        return weighting_factors



    def create_controls(self, experiment):
        """Creates control variables for a single experiment and save them in a dictionary.

        :param experiment: Experiment object
        :type experiment: BioMoni.Experiment.Experiment
        :return: c: Control variables as dictionary.
        
        """

        c = {"feed_%_to_mL" : experiment.metadata.loc["feed_%_to_mL"],"feed_on" : experiment.metadata.loc["tsfeedOn"], "add_IPTG": experiment.metadata.loc["add_IPTG"], "feed_rate_mL_to_g" : experiment.metadata.loc["feed_rate_mL_to_g"], "csf" : experiment.metadata.loc["csf"]
        , "cNH3" : experiment.metadata.loc["cNH3"], "gas_flow_lpm" : experiment.metadata.loc["gas_flow"], "T" : experiment.metadata.loc["T"]}
        
        c["gas_flow"] = c["gas_flow_lpm"] * 60 #turn  Lpm to L/h

        c["pressure"] = self.calc_mean_pressure(experiment)

        feedpump_power = pd.to_numeric(experiment.dataset_raw["on"]["SUBS_A2","Value"] , downcast="float" , errors="coerce") #in some experiments the feed_pump from MFCS was recognized as string instead of float...
        #feedpump_power.fillna(0)       #nan values to zero
        feedrate = (feedpump_power * c["feed_%_to_mL"] * 60) / 1000                             #Feedrate in [L/h]
        feedrate_glc = feedrate * c["csf"]                                                      #Feedrate of pure Glucose in [g/h]
        #create interpolation function which can be called in the kinetics function at specific time points. if t < first value y assumes first value, if t > last value y assumes last value
        c["feedrate_glc"] = interp1d(x = feedrate_glc.index, y = feedrate_glc, fill_value = (feedrate_glc.iloc[0], feedrate_glc.iloc[-1]) , bounds_error= False)
        c["feedrate_volume"] = interp1d(x = feedrate.index, y = feedrate, fill_value = (feedrate.iloc[0], feedrate.iloc[-1]) , bounds_error= False)
        
        if "SampleVolume [g]" in experiment.dataset["off"]:
            c["Fout"] = self.calc_Fout(experiment)


        return c  

    def calc_Fout(self, experiment):
        """Calculates outflow Fout due to sampling from provided SampleVolume in offline data.
        It turns SampleVolume into Flowrate that is necessary to flow for a period of 10 min and result in the given amount of sample volume
        :param experiment: Experiment object
        :type experiment: BioMoni.Experiment.Experiment
        :return: Flow: inter1pd object describing outflow
        
        """
            
        if "SampleVolume [g]" in experiment.dataset["off"]:
                
                b = experiment.dataset["off"]["SampleVolume [g]"]

                #turning Sample Volume into Flowrate L/h that is necessry to flow for a period of 10 min to result the amount of sample volume
                #assumption: 1g = 1mL sample 
                c= pd.DataFrame(b*0.012)
                c = c.rename(columns={"SampleVolume [g]":"VolumeFlow [mL/h]"})

                e_after= [x + (5/60) for x in c.index]
                e_before = [x - (5/60) for x in c.index]
            
                for i in e_after:
                    c.loc[i]= 0

                for i in e_before:
                    c.loc[i]=0

                Flow_raw= c.sort_index(axis ="index")
                

                Flow = interp1d(x = Flow_raw.index, y = Flow_raw["VolumeFlow [mL/h]"] ,fill_value = (Flow_raw.iloc[0], Flow_raw.iloc[-1]), bounds_error= False)
                return Flow


    def create_y0(self, experiment):
        """Creates initial state vector.
            * y0[0] : initial Bio dry mass (mX) [g]
            * y0[1] : initial Substrate mass (mS) [g]  
            * y0[2] : initial Product mass(mP), in this case Riboflavin [mg]
            * y0[3] : initial Volume of fermentation broth [L]
        :param experiment: Experiment object
        :type experiment: BioMoni.Experiment.Experiment
        :return: y0: Initial state vector as list.
        
        """
       
        y0 = list(experiment.metadata.loc[["mX0","mS0","mP0", "V0"]].values)        #mX0 [g], mS0 [g], mP0 [mg], V0 [L]
        
    
        return y0

    def kinetics(self, t, y, p, c):

        """ r.h.s. (right hand side) function of the ODE model.

        :param t: Current time [h]
        :type t: float
        :param y: State vector:
            * y[0] : Bio dry mass (mX) [g]
            * y[1] : Substrate mass (mS) [g]  
            * y[2] : Product mass(mP), in this case Riboflavin [mg]
            * y[3] : Volume of fermentation broth [L]

        :type y:     array like
        :param p:    Structure with parameter values to be estimated, cf. lmfit.parameter.Parameters
        :type p:     lmfit.parameter.Parameters object
        :param c:    Dictionary with control values

           
            * c["feed_rate"]:           [L/h]
            * c["feedrate_function"]:   interpolation function of feed reate [L/h]
            * c["csf"]:                 Substrate concentration in feed [g/L]
            * c[cNH3"]:                 Base concentration, NH3 [g/L]
            * c["gas_flow"]:            Aeration rate [L/h]
            * c["T"]:                   Temperature in reactor [°C]
            * c["wf_on"]:               Weighting factor for online data [dimensionless]
            * c["wf_CO2"]:              Weighting factor for CO2 data [dimensionless]
            * c["pressure"]:            Mean pressure in CO2 sensor [bar]
            

        :type c: dict
        :param yields: Calculated yields coefficients
        :type yields: dict 
                
        :return: Fin, Fout, mu, rS, rP, rX, rP, V  kinetic rates at time t. 
        """

        # fit parameters
        mu_max = p["mu_max"].value
        Yxs = p["Yxs"].value
        Km = p["Km"].value
        Ypx = p["Ypx"].value
        #Yxsmain = p["Yxsmain"].value
        
        
        
        # controls
        #feed_on = c["feed_on"] # time point when feed was switched on [h]
        #feed_rate = c["feed_rate"] # feed rate [L/h]
        #Fin = feed_rate * (t > feed_on) # becomes 0 if t < feed_on    #this was used in Yeast.py when using a constant feedrate
        
        Fin = c["feedrate_volume"](t) #this is used when the feed rate changes over time: for Valeria 
        Fin = np.nan_to_num(Fin)        #if Nan oder infer values are generated, it cant be colved by solve_ivp anymore
        
        if "Fout" in c:
            Fout = c["Fout"](t)
            Fout = np.nan_to_num(Fout)
        else:
            Fout = 0  
        
        
        # masses and concentrations 
        mX, mS, mP, V = y
                   
        cX,cS,cP = [mX,mS,mP] / V 

        #kinetics [-]
        mu_Mo = cS / (Km + cS)
        

        #growth rates [1/h]
        mu = mu_max * mu_Mo 
        
        #Biomass growth [g/h]
        rX = mu * cX * V 
        
        #Substrate consummption for maintanance [g/h]
        #rmainS = rX *1/Yxsmain
        
        #SUbstrate consumption [g/h]
        rS = (rX * 1/Yxs) #+ rmainS
        
        #Product formation [mg/h]
        rP = rX * Ypx
        
        
        return Fin, Fout, mu, rS, rP, rX, rP,V 


    def model_rhs(self, t, y, p, c):
        """ This function calculates the time derivatives to be integrated by solve_ivp by processing the outcome of the kinetics function.

        :param t: Current time [h]
        :type t: float
        :param y: State vector:

            * y[1] : Substrate mass (mS) [g]  
            * y[0] : Bio dry mass (mX) [g]
            * y[3] : Product mass(mP), in this case Riboflavin [mg]
            * y[2] : Volume of fermentation broth [L]

        :type y: array like
        :param p: Structure with parameter values to be estimated, cf. lmfit.parameter.Parameters
        :type p: lmfit.parameter.Parameters object
        :param c:  Dictionary with control values

            * c["feed_on"]:             Time point when feed was switched on [h]
            * c["feed_rate"]:           Feed rate [L/h]
            * c["csf"]:                 Substrate concentration in feed [g/L]
            * c["feed_%_to_mL]: 
            * c["feed_rate_mL_to_g"]:
            * c[cNH3"]:                 Base concentration, NH3 [g/L]
            * c["gas_flow"]:            Aeration rate [L/h]
            * c["T"]:                   Temperature in reactor [°C]
            * c["wf_on"]:               Weighting factor for online data [dimensionless]
            * c["wf_CO2"]:              Weighting factor for CO2 data [dimensionless]
            * c["pressure"]:            Mean pressure in CO2 sensor [bar]

        :type c: dict
        :param yields: Calculated yields coefficients
        :type yields: dict 
                
        :return: dy_dt: Time derivative of state vector. State vector describes biomass, susbtrate, product and volume
        """
        
        csf = c["csf"] # substrate concentration in feed [g/L]  
        

        Fin, Fout, mu, rS, rP, rX, rP, V = self.kinetics(t, y, p, c)   
    
          
        dmX_dt = rX 
        dV_dt  = Fin - Fout       #later: insert here Fout = V Sample at distinct timepoints
        dmS_dt = -rS + (Fin * csf) 
        dmP_dt = rP 
        
        return dmX_dt, dmS_dt, dmP_dt, dV_dt

    
    
    def observation(self, t_grid, y0, p, c, **kwargs_solve_ivp):
        """ The observation function calculates simulated values by processing the outcome of the functions model_rhs and calc_CO2.

        :param t_grid: time grid with measurement time points [h]
        :type t_grid: array like
        :param y0: Initial state vector:

            
            * y[1] : Substrate mass (mS) [g]  
            * y[0] : Bio dry mass (mX) [g]
            * y[2] : Product mass(mP), in this case Riboflavin [mg]
            * y[3] : Volume of fermentation broth [L

        :type y: array like
        :param p: Structure with parameter values to be estimated, cf. lmfit.parameter.Parameters
        :type p: lmfit.parameter.Parameters object
        :param c:  Dictionary with control values

            * c["feed_on"]:             Time point when feed was switched on [h]
            * c["feed_rate"]:           Feed rate [L/h]
            * c["csf"]:                 Substrate concentration in feed [g/L]
            * c["feed_%_to_mL]: 
            * c["feed_rate_mL_to_g"]:
            * c[cNH3"]:                 Base concentration, NH3 [g/L]
            * c["gas_flow"]:            Aeration rate [L/h]
            * c["T"]:                   Temperature in reactor [°C]
            * c["wf_on"]:               Weighting factor for online data [dimensionless]
            * c["wf_CO2"]:              Weighting factor for CO2 data [dimensionless]
            * c["pressure"]:            Mean pressure in CO2 sensor [bar]

        :type c: dict
        :param yields: Calculated yields coefficients
        :type yields: dict

        :param **kwargs_solve_ivp: optional arguments for scipy.integrate.solve_ivp

        :return: sim_exp: Simulated variables for experiment.
        
        
        """
        y0_mod = deepcopy(y0) 
        
        # for viability factor 
        
        #y0_mod[0] = y0_mod[0]* p["viab_f"].value  
        #print("viab_f is "+ str(p["viab_f"].value))
        #print("mX0 changed to "+ str(y0_mod[0]))
        
        #result of IVP solver
        y_t = solve_ivp(self.model_rhs, [0, np.max(t_grid)], y0_mod, t_eval=t_grid, args = (p, c), **kwargs_solve_ivp).y.T
    

        # unpack solution into vectors
        mX, mS, mP, V= [y_t[:,i] for i in range(4)]
        cX, cS, cP  = [mX, mS, mP] / V
        

        #dmX_dt = np.array([self.model_rhs(t = t_grid[i], y = y_t[i,:], p = p, c = c, yields = yields) for i in range(len(t_grid))])[:,0]  #calc dmX_dt from model funcs
        #BASET_rate_wrong = dmX_dt * 1000  #base_rate in ml/h
        
        #for later to visualize 
        if "Fout" in c:
            Fout = c["Fout"](t_grid)
            Fout = np.nan_to_num(Fout)
        else:
            Fout = 0
        #CO2_percent = np.array([self.calc_CO2(t = t_grid[i], y = y_t[i,:], p = p, c = c, yields = yields) for i in range(len(t_grid))])     #CO2 in vol.% from calc_CO2
        
        #calculating carbon recovery 
        
        #CRR = self.calc_CRR()
        

                
        dict_sim_exp = {("t"): t_grid, ("V"): V, ("cX", "CDW_calc"):cX, 
                        ("BASET_rate","BASE"): 1, ("cS", "Glucose [g/L]"): cS,
                        ("cP","RF [mg/L]"): cP, ("Fout","F_out"): Fout}                 #necessary to extend the dict if col names are different in offline/online datasets for different Fermentation runs, warning: give BASET_rate column also another alternative name even if BASET_rate has no other names. If you dont: you will get each letter as a seperate column and t will be overwritten!  
        
        dict_sim_exp = {key: value for keys, value in dict_sim_exp.items() for key in keys}                             
        #simulated values
        sim_exp = pd.DataFrame(dict_sim_exp).set_index("t")           
        
        

        return sim_exp
        

    
    def simulate(self, experiment = None, t_grid = None, y0 = None, p = None,  c = None, kwargs_solve_ivp = dict(method= "Radau", first_step = 0.0000001, max_step= 0.1)       # kwargs_solve_ivp may have to be adapted
    , t_step = 0.001, endpoint = None, **kwargs_solve_ivp_given):  

        """
        Simulates values by calling the observation function, it is possible to simulate manually by using values for t_grid, p, y0 c and yields or via an Experiment object. 
        If using an Experiment object, it is still possible to give manual values for  t_grid, p, y0 c. 
        Values not given are created from the Experiment object."""

        if kwargs_solve_ivp_given:       #check if kwargs_direct are given to simulate, if yes use given kwargs if not use kwargs_solve_ivp handed over from estimate or default kwargs_solve_ivp
            kwargs_solve_ivp = kwargs_solve_ivp_given #it helps if you want to simulate manually with your own kwargs

        
          
        if experiment is not None:      #This if statement is to directly simulate from experiment, it is possible to give an experiment with modified t,y,p,c
            assert isinstance(experiment, Experiment), "Given experiment must be of type Experiment"

            if t_grid is None:

                if experiment.metadata["start"] is None:

                    ts_values_list = [df["ts"].values for df in experiment.dataset.values()]    #list with ts values for each measurement type
                    ts_values_lowest_ts = min(ts_values_list, key = np.min)     #ts values with lowest ts
                    start_ts = np.min(ts_values_lowest_ts)     #lowest ts

                else:

                    start_ts = experiment.metadata["start"]


                if endpoint is None:        #if the endpoint is not given, the same endpoint from Experiment creation is used
                    endpoint = experiment.endpoint
                else:                    
                    assert endpoint in experiment.metadata.index, "Given endpoint (\"end1\" or similar) must be in metadata columns"   #.index because of pd.series
                if experiment.metadata[endpoint] is None: 

                    t_grid_list = [df.index.values for df in experiment.dataset.values()]   #get t_grid for every typ in experiment as list
                    t_grid_largest_endpoint = max(t_grid_list, key = np.max)       #get t_grid with largest endpoint
                    end = np.max(t_grid_largest_endpoint) #largest endpoint               
                    
                else:

                    end = (experiment.metadata[endpoint] - start_ts) / pd.Timedelta(1,"h") 

                t_grid = np.linspace(0, end, t_step) #timeline for simulation, going from 0 to end. Simulated values start always at zero in contrast to real measurement values
            
            if y0 is None:
                y0 = deepcopy(self.create_y0(experiment))
               
                
            
            if p is None:
                p = self.p

            if c is None:
                c = self.create_controls(experiment)
            
            # if yields is None:
            #     yields = self.chemical_balancing(p)
            
              
            
        elif experiment is None:    #This if statement is to simulate with given t,y,p,c (used in the estimate process because it would require more computational power to calculate t, y and c within each estimate iteration)
            for i in [t_grid, p , y0, c]:
                if i is None:
                    raise ValueError("If no Experiment is given, it is necessary to give t_grid, p , y0, c and yields in order to simulate")

        
        assert isinstance(t_grid, (np.ndarray, list)), "t_grid must be either of type np.ndarray or list"
        assert isinstance(y0, (np.ndarray, list)), "y0 must be either of type np.ndarray or list"
        assert isinstance(p, type(Parameters())), "Given parameters must be of type lmfit.parameter.Parameters"
        assert isinstance(c, dict), "c must be a dictionary"
        #assert isinstance(yields, dict), "yields must be a dictionary"


        sim_exp = self.observation(t_grid, y0, p, c, **kwargs_solve_ivp)

        return sim_exp
   

    def residuals(self, p, datasets_dict, settings_dict, tau, kwargs_solve_ivp):
        """Function to calculate the residuals between simulated data and measurement data 

        :param p: Structure with parameter values to be estimated, cf. lmfit.parameter.Parameters
        :type p: lmfit.parameter.Parameters object
        :param datasets_dict: Dictionary with key: Experiment number, value: dataset with data from different measurement methods (e.g. off, on CO2 data)
        :type datasets_dict: dict
        :param settings_dict: Dictionary with settings for the simulation for each given Experiment. Key: Experiment number, value: settings containing:

            * y0: Initial state vector 
            * c: Control variables
            * wf: Weighting factors

        :type settings_dict: dict
        :param kwargs_solve_ivp: arguments for scipy.integrate.solve_ivp, coming from the estimate function
        :type kwargs_solve_ivp: dict

        _return: res_all: Super long residual vector containing the residuals of all experiments.
        """

        # yields = self.chemical_balancing(p)  #yields are calculated within the parameterestimation, thus influencing the parameter choice    

        res_all = []
        
        
        for exp_id, dataset in datasets_dict.items():
            
            res_single = np.array([]) #empty array which will contain residuals
            settings = settings_dict[exp_id]
            y0 = settings["y0"] #extract y0,c and weighting factors for each experiment
            
            
            
            c = settings["c"]
            weighting_factors = settings["wf"]
            

            for dskey, dat in dataset.items():   #extract data for ("on", "off", "CO2") for each experiment
                wf = weighting_factors[dskey]       #single weighting factor       
                t_grid = dat.index.values           #extract time from dataframe
                sim_exp = self.simulate(None, t_grid, y0, p, c, kwargs_solve_ivp = kwargs_solve_ivp)
                
                
                
                
                for var in dat:     #loop over measured variables (columns)
                    if var in sim_exp.columns:
                        
                        x = ["Glucose [g/L]","RF [mg/L]","Acetat [g/L]","CDW_calc",("BASET_2","Value"),('pO2_2', 'Value'),"CO2"]
                       
                        if var in x:
                    
                           #print("difference in " + str(var) + " is processed.")
                           res_var1 = (weighting_factors[dskey][var] * ((sim_exp[var] - dat[var]).values))
                           res_single = np.append(res_single, res_var1)
                            
                        
                            
                        #res_single = np.append(res_single, res_var)
            else:

                if tau is not None:     
                    weighting_decay = np.exp( (dat.index-np.max(dat.index) ) / tau)  #calculating exponential weighting decay 
                else:
                    weighting_decay = 1    #if no tau given: no additional weighting decay over time

                for var in dat:     #loop over measured variables (columns)
                    if var in sim_exp.columns:
                        print("unexpected entering of cycle")                       
                        res_var = (wf * ((sim_exp[var] - dat[var]).values)) * weighting_decay        # weighted residuals for this measured variable
                        res_single = np.append(res_single, res_var) # append to long residual vector
                else:
                    pass
            res_all = np.append(res_all, res_single)        #super long vector with all residuals
        
        
        return res_all
    
    def calc_CRR(self,experiment, V, c , y0):
        
        Mw_bm = 24.445          #g/mol molecular weight of biomass assuming composition of C=1,H=1.594,N=0.293,O=0.387,P=0.01,S=0005
        Mw_gluc = 180.156       #g/mol molecular weight of glucose C=6,H=12,O=6
        Mw_RF = 376.36          #g/mol molecular weight of Riboflavin C=17,H=20,N=4,O=6
        Mw_CO2 = 44.01
        Mw_Ac = 59.04
        R = 0.08314             #bar*l/mol*K
        T = c["T"] + 273.15     #Kelvin
        dV_gas_dt = c["gas_flow"]
        pressure = c["pressure"] #bar
        
        t_grid = list(experiment.dataset["off"].index.values)
        
        F_C_in = c["feedrate_glc"](t_grid)
        
        #initiate df where all relevant values for calculation will be stored
        df = pd.DataFrame(data=F_C_in,index=t_grid, columns=["feedrate pure glucose g/h"])
        
        # calculate cumulated amount of Glucose from Feed
        df["Glc g_cum"] = scipy.integrate.cumtrapz(df["feedrate pure glucose g/h"],df["feedrate pure glucose g/h"].index, initial=0)
        
        
        
        if "CO2" in experiment.dataset.keys():
            #creating column with calculcated Volume values for CO2 [L] from exhaust gas measurment
            experiment.dataset["CO2"]["V_CO2 L"] = (experiment.dataset["CO2"]["CO2"]- 0.04)/100 * dV_gas_dt   #30 L/h is F_exhaust or Outlet flow (assuming dV_gas_dt = inlet flow)

            ##creating column with cumulated volume values for CO2 [L] 
            experiment.dataset["CO2"]["cum_CO2 L"] = scipy.integrate.cumtrapz(experiment.dataset["CO2"]["V_CO2 L"],experiment.dataset["CO2"]["V_CO2 L"].index, initial=0)

            CO2_func = interp1d(x = experiment.dataset["CO2"].index, 
                                y = experiment.dataset["CO2"]["CO2"], fill_value = (experiment.dataset["CO2"]["CO2"].iloc[0], experiment.dataset["CO2"]["CO2"].iloc[-1]) , 
                                bounds_error= False)

            V_CO2_func = interp1d(x = experiment.dataset["CO2"].index, 
                                y = experiment.dataset["CO2"]["V_CO2 L"], fill_value = (experiment.dataset["CO2"]["V_CO2 L"].iloc[0], experiment.dataset["CO2"]["V_CO2 L"].iloc[-1]) ,
                                bounds_error= False)
            
            V_CO2_cum_func= interp1d(x = experiment.dataset["CO2"].index, 
                    y = experiment.dataset["CO2"]["cum_CO2 L"], fill_value = (experiment.dataset["CO2"]["cum_CO2 L"].iloc[0], experiment.dataset["CO2"]["cum_CO2 L"].iloc[-1]) ,
                    bounds_error= False)

            df["CO2 %"] = CO2_func(t_grid)
            df["V_CO2 L"] = V_CO2_func(t_grid)
            df["cum_CO2 L"] = V_CO2_cum_func(t_grid)
            
            #calculate cumulated amount of CO2 in exhaust gas 
            df["CO2 mol_cum"] = (df["cum_CO2 L"] * pressure) / (R * T)
            
            
            
            #calculation of recovery rate

            Carb_RR= {}
            Carb_RR["t"]= t_grid
            Carb_RR["given"]=[]
            Carb_RR["found"]=[]
            
            Carb_RR["mol glc"]=[]
            Carb_RR["mol bm"]=[]
            Carb_RR["mol RF"]=[]
            Carb_RR["mol Ac"]=[]
            Carb_RR["mol CO2"]=[]
            
            Carb_RR["RR"]=[]

            for t in t_grid:
                
                #mass balance of Carbon given (Substrate,Biomass,Product)
                
                nC_Sg = y0[1]/Mw_gluc *6 + df["Glc g_cum"].values[t_grid.index(t)]/Mw_gluc * 6   #initial Substrate + Substrate from Feed
                nC_Xg = float(experiment.metadata.loc[["mX0_ohne_viab_f"]].values)/Mw_bm
                nC_Pg = (y0[2]/1000)/Mw_RF *17
                
                amountC_given = nC_Sg + nC_Xg + nC_Pg 
                
                Volume = V(t)
                
                #mass balance of Carbon found (Substrate,Biomass,Riboflavin,Acetate,CO2)
                nC_Sf = (float(experiment.dataset["off"]["Glucose [g/L]"].loc[[t]])*Volume)/Mw_gluc * 6
                nC_Xf = (float(experiment.dataset["off"]["CDW_calc"].loc[[t]])*Volume)/ Mw_bm
                nC_RFf = (float(experiment.dataset["off"]["RF [mg/L]"].loc[[t]])*Volume /1000)/ Mw_RF * 17
                nC_Acf = (float(experiment.dataset["off"]["Acetate [g/L]"].loc[[t]])*Volume)/ Mw_Ac
                nC_CO2 = float(df["CO2 mol_cum"].values[t_grid.index(t)])

                amountC_found = nC_Sf + nC_Xf + nC_RFf + nC_Acf + nC_CO2
                
                RR = (amountC_found / amountC_given) * 100
                
                Carb_RR["given"].append(amountC_given)
                Carb_RR["found"].append(amountC_found)
                
                Carb_RR["mol glc"].append(nC_Sf)
                Carb_RR["mol bm"].append(nC_Xf)
                Carb_RR["mol RF"].append(nC_RFf)
                Carb_RR["mol Ac"].append(nC_Acf)
                Carb_RR["mol CO2"].append(nC_CO2)
            
                Carb_RR["RR"].append(RR)
            
            df2 = pd.DataFrame.from_dict(Carb_RR)
            

            df["RR"] = df2["RR"].values 
        else:
          df["RR"] = 0
          Carb_RR= {}
          df2 = pd.DataFrame.from_dict(Carb_RR)      

        return df, df2