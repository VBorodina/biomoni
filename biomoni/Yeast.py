from biomoni import Experiment
from biomoni import Model
import pandas as pd
import numpy as np

from scipy.integrate import solve_ivp

import sympy as sp
from sympy import Eq
from sympy import symbols

from copy import deepcopy


import time #in Model
from lmfit import Parameters, report_fit, minimize  #in Model




class Yeast(Model):     #Dependent on base class
    """
    :class Yeast: A class to create a yeast model, the class can simulate a yeast fermentation and estimate parameters based on lab data. This class is inherited from the base class: Model.

    """


    def __repr__(self):
        """Representation of the Yeast object in the print() call"""
        
        return  "Yeast()" 

    def __init__(self):
        """Initialization of an object of the Yeast(Model) class

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
            p.add("qsmax", value= 1.6, min=0.01, max=5.0 , vary = True)     #Maximum glucose uptake rate [g/(g*h)]   #previous start val: 0.1  #estimated value: 1.61290497
            p.add("qemax", value=0.2361, min=0.15, max=0.35, vary=False)    #Maximum ethanol uptake rate [g/(g*h)]       #value=0.2361, min=0.15, max=0.35, vary=True # value=0.05, min=0.001, max=0.5,  vary=True geht

            p.add("base_coef", value= 0.007395 , min=0.0001, max =3, vary = True)     #base coefficient [mol/g]     #previous start val: 1 #estimated value: 0.00733798      
            p.add("qO2max", value= 0.164745, min=0.1, max=0.4, vary = True)    #Maximum oxygen uptake rate [g/(g*h)]  , min=0.1  0.1, 0.4  #previous start val: 0.255984, 1 #estimated value: 0.17451356
            p.add("qm_max", value=0.01, min=0.0075, max=0.0125, vary = False)  #glucose uptake rate required for maintenance [g/(g*h)]  
            
            p.add("Ks", value=0.1, min=0.01, max=1,  vary=False)   #Saturation constant, concentration of glucose at µ = 0.5 µmax [g/L]
            p.add("Ke", value=0.1, min=0.01, max=1, vary=False)    #Saturation constant, concentration of ethanol at µ = 0.5 µmax [g/L]
            p.add("Ki", value=0.1, min=0.01, max=1, vary=False)    #Inhibition constant, glucose inhibits uptake of ethanol [g/L]

            p.add("Yxs_ox", value= 0.5389, min=0.3, max=0.6, vary=True) #Yield biomass per glucose (oxidative growth) [g/g] #previous start val: 49    #estimated value 0.52792125
            p.add("Yxs_red", value=0.05, min=0.01, max=0.8, vary=False)    #Yield biomass per glucose (reductive growth) [g/g] 
            p.add("Yxe_et", value=0.72, min=0.5, max=0.8, vary=False)      #Yield biomass per ethanol [g/g]
            p.add("Yxg_glyc", value=0.2, min=0.1, max=0.35, vary=False)    #Yield biomass per glycerol [g/g]              #value=0.2, min=0.1, max=0.35, vary=False

            #Biomass composition paramaters, content H,O,N from paper sonnleitner

            p.add("HX", value=1.79, min=1.77, max=2.1, vary=False) #Stoichiometric hydrogen content of biomass [mol/mol]
            p.add("OX", value=0.57, min=0.54, max=0.63, vary=False) #Stoichiometric oxygen content of biomass [mol/mol]     
            p.add("NX", value=0.15, min=0.14, max=0.16, vary=False) #Stoichiometric nitrogen content of biomass [mol/mol]

            p.add("g_e", value=0.21636282532087284, min=0.17624196916440354, max=0.2883951412907231, vary=False) #determined experimentally: formation glycerol per ethanol [g/g]

            self.p = p
            
        else:

            assert isinstance(p, type(Parameters())), "Given parameters must be of type lmfit.parameter.Parameters"
            self.p = p
        
        self.yields = self.chemical_balancing(self.p)

    
        
    def chemical_balancing(self,p):
        """Calculate yield coefficients from p for the following metabolic pathways:

        1. (ox) : Oxidative glucose consumption
        2. (red) : Reductive glucose consumption
        3. (et) : Ethanol consumption
        4. (glyc) : Glycerol consumption
        5. (m) : Maintenence metabolism

        :param p: Parameters 
        :type p: lmfit.parameter.Parameters object

        :return: yields: Yields coefficients of metabolic pathways
        """
        #Order: C,H,O,N numbers of atoms in the molecule e.g. glucose : C6H12O6N0
        gluc = np.array([6.0,12.0,6.0,0.0])
        O2 = np.array([0.0, 0.0, 2.0, 0.0])
        NH3 = np.array([0.0,3.0,0.0,1.0])
        biomass = np.array([1.0,p["HX"].value, p["OX"].value, p["NX"].value])
        CO2 = np.array([1.0,0.0,2.0,0.0])
        H2O = np.array([0.0,2.0,1.0,0.0])
        etoh = np.array([2.0,6.0,1.0,0.0])
        glyc = np.array([3.0,8.0,3.0,0.0])


        #Required yields to solve linear equation system taken from p
        yields = {"Yxs_ox": p["Yxs_ox"].value , "Yxs_red" : p["Yxs_red"].value      
        , "Yxe_et" : p["Yxe_et"].value, "Yxg_glyc" : p["Yxg_glyc"].value}

        


        MW_element_dict = {"C": 12.011, "H": 1.0079, "O": 15.999, "N": 14.007}        #molar masses of elements
        molecule = {"gluc": gluc, "O2": O2, "NH3" : NH3, "biomass": biomass, "CO2" : CO2, "H2O":  H2O, "etoh": etoh, "glyc": glyc} #molecules dict

        MW = {}           #MW = molar weights: dict with masses of molecules
        for key, mol in molecule.items():
            molecule_MW_array = ([])
            for vectorvalue, weight in zip (mol, MW_element_dict.values()):
                    vw = vectorvalue*weight
                    molecule_MW_array= np.append(molecule_MW_array, vw)
            MW[key] = sum(molecule_MW_array)

        NX1 = p["NX"].value #NX1 because NX is reserved for the symbol letter NX from sympy to solve the equation
        GE = (p["g_e"].value/MW["glyc"]) * MW["etoh"]   #from mass ratio(p["g_e"]) to molar ratio. Glycerol per ethanol

        #1. oxidative glucose consumption: gluc+ a*O2 + b*NX*NH3 = b*biomass + c*CO2 + d*H2O 
        a,b,c,d, NX = symbols("a b c d NX")         #set symbols to solve equation
        Yxs_ox = p["Yxs_ox"].value
        b1 = Yxs_ox* MW["gluc"]/MW["biomass"]     #calculate stoichiometric coefficient from mass related yield coefficient

        eqOx_list = []
        for num in range(3):        #Set up reaction equations three times because three unknowns are present
            eqOx = Eq(gluc[num]+ a*O2[num]+ b*NX*NH3[num], b*biomass[num]+ c*CO2[num]+ d*H2O[num])
            eqOx = eqOx.subs({b: b1, NX: NX1})
            eqOx_list.append(eqOx)
        
        solution_Ox = sp.solve(eqOx_list, (a, c, d), dict= True)        #solve equation system
        a1, c1, d1 = np.float(solution_Ox[0][a]), np.float(solution_Ox[0][c]), np.float(solution_Ox[0][d])  #assign results to variables

        YCO2x_ox = c1/b1 * MW["CO2"]/MW["biomass"]
        YCO2s_ox = c1/1 * MW["CO2"]/MW["gluc"]
        YO2s_ox = a1/1 * MW["O2"]/MW["gluc"]
        
        yields["YCO2x_ox"], yields["YCO2s_ox"], yields["YO2s_ox"] = YCO2x_ox, YCO2s_ox, YO2s_ox

        #2. reductive glucose consumption:  gluc+ g*NX*NH3 = g*biomass + h*CO2 + i*H2O + j*etOH + GE*j*glyc
        g,h,i,j, NX = symbols("g h i j NX")
        Yxs_red = p["Yxs_red"].value
        g1 = Yxs_red* MW["gluc"]/MW["biomass"]

        eqRed_list = []
        for num in range(3):
            eqRed = Eq(gluc[num]+ g*NX*NH3[num] , g*biomass[num]+ h*CO2[num]+ i*H2O[num]+ j*etoh[num]+ GE*j*glyc[num] )
            eqRed = eqRed.subs({g: g1, NX: NX1})
            eqRed_list.append(eqRed)
        
        solution_Red = sp.solve(eqRed_list, (h, i, j), dict= True)
        h1,i1,j1 = np.float(solution_Red[0][h]), np.float(solution_Red[0][i]), np.float(solution_Red[0][j])

        Yes_red = j1/1 * MW["etoh"]/MW["gluc"]
        Ygs_red = GE*j1/1 * MW["glyc"]/MW["gluc"]   
        YCO2x_red = h1/g1 * MW["CO2"]/MW["biomass"]
        YCO2s_red = h1/1 * MW["CO2"]/MW["gluc"]


        yields["Yes_red"], yields["Ygs_red"], yields["YCO2x_red"], yields["YCO2s_red"] = Yes_red, Ygs_red, YCO2x_red, YCO2s_red


        #3. oxidative ethanol consumption: etoh + k*O2 + l*NX*NH3 = l*biomass + m*CO2 + n*H2O
        k,l,m,n, NX = symbols("k l m n NX")
        Yxe_et = p["Yxe_et"].value      
        l1 = Yxe_et* MW["etoh"]/MW["biomass"]

        eqEt_list = []
        for num in range(3):
            eqEt = Eq(etoh[num]+ k*O2[num]+ l*NX*NH3[num], + l*biomass[num]+ m*CO2[num]+ n*H2O[num])
            eqEt = eqEt.subs({l: l1, NX: NX1})
            eqEt_list.append(eqEt)
        
        solution_Et = sp.solve(eqEt_list, (k, m, n), dict= True)
        k1, m1, n1 = np.float(solution_Et[0][k]), np.float(solution_Et[0][m]), np.float(solution_Et[0][n])

        YCO2x_et = m1/l1 * MW["CO2"]/MW["biomass"]
        YCO2e_et = m1/1 * MW["CO2"]/MW["etoh"]
        YO2e_et = k1 * MW["O2"]/MW["etoh"]


        yields["YCO2x_et"], yields["YCO2e_et"], yields["YO2e_et"] = YCO2x_et, YCO2e_et, YO2e_et

        #4. oxidative glycerol consumption: glyc+ u*O2 + v*NX*NH3 = v*biomass + w*CO2 + x*H2O
        u,v,w,x,NX = symbols("u v w x NX")
        Yxg_glyc = p["Yxg_glyc"].value
        v1 = Yxg_glyc * MW["glyc"]/MW["biomass"]

        eqGlyc_list = []
        for num in range(3):
            eqGlyc = Eq(glyc[num]+ u*O2[num]+ v*NX*NH3[num] , v*biomass[num]+ w*CO2[num]+ x*H2O[num])
            eqGlyc = eqGlyc.subs({v: v1, NX: NX1})
            eqGlyc_list.append(eqGlyc)
        
        solution_glyc = sp.solve(eqGlyc_list, (u, w, x), dict= True)
        u1, w1, x1 = np.float(solution_glyc[0][u]), np.float(solution_glyc[0][w]), np.float(solution_glyc[0][x])

        YCO2x_glyc = w1/v1 * MW["CO2"]/MW["biomass"]
        YCO2g_glyc = w1/1 * MW["CO2"]/MW["glyc"]
        YO2g_glyc = u1/1 * MW["O2"]/MW["etoh"]
        

        yields["YCO2x_glyc"], yields["YCO2g_glyc"], yields["YO2g_glyc"] = YCO2x_glyc, YCO2g_glyc, YO2g_glyc

        #5. maintenance metabolism : gluc + 6*O2 = 6*CO2 + 6*H2O

        YCO2s_m = 6 * MW["CO2"]/MW["gluc"]
  

        yields["YCO2s_m"] = YCO2s_m
        
        return yields


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
        """Calculates weighting factors from an single experiment: wf(i) = len("off") dataset / len(i) dataset for i = "CO2" or "on".
        
        :param experiment: Experiment object
        :type experiment: BioMoni.Experiment.Experiment
        :return: weighting factors for single experiment as dictionary
        """

        df_smallest_len = min(experiment.dataset.values() , key = len)   #return dataframe with smallest length
        weighting_factors = {}
        for dskey in experiment.dataset.keys():
    
            weighting_factors[dskey] = len(df_smallest_len) / len(experiment.dataset[dskey]) #weihting factor calculated by len smallest dataframe / len actual dataframe. This results in data frames with different lengths being weighted equally in the estimation.
      
        return weighting_factors



    def create_controls(self, experiment):
        """Creates controls variables for a single experiment and save them in a dictionary.

        :param experiment: Experiment object
        :type experiment: BioMoni.Experiment.Experiment
        :return: c: Control variables as dictionary.
        
        """

        c = {"feed_on" : experiment.metadata.loc["feed_on"], "feed_rate" : experiment.metadata.loc["feed_rate"], "csf" : experiment.metadata.loc["csf"]
        , "M_base" : experiment.metadata.loc["M_base"], "gas_flow" : experiment.metadata.loc["gas_flow"], "T" : experiment.metadata.loc["T"]}

        c["pressure"] = self.calc_mean_pressure(experiment)

        return c  

    def create_y0(self, experiment):
        """Creates initial state vector.

        :param experiment: Experiment object
        :type experiment: BioMoni.Experiment.Experiment
        :return: y0: Initial state vector as list.
        
        """

        y0 = list(experiment.metadata.loc[["mS0", "mX0", "mE0","V0"]].values)
        
        return y0

    def kinetics(self, t, y, p, c, yields):

        """ r.h.s. (right hand side) function of the ODE model.

        :param t: Current time [h]
        :type t: float
        :param y: State vector:

            * y[0] : Substrate mass (mS) [g]  
            * y[1] : Bio dry mass (mX) [g]
            * y[2] : Ethanol mass (mE) [g]
            * y[3] : Volume of fermentation broth [L]

        :type y: array like
        :param p: Structure with parameter values to be estimated, cf. lmfit.parameter.Parameters
        :type p: lmfit.parameter.Parameters object
        :param c:  Dictionary with control values

            * c["feed_on"]: Time point when feed was switched on [h]
            * c["feed_rate"]: Feed rate [L/h]
            * c["csf"]: Substrate concentration in feed [g/L]
            * c["M_base"]: Base concentration [mol/L]
            * c["gas_flow"]: Aeration rate [L/h]
            * c["T"]: Temperature in reactor [°C]
            * c["wf_on"]: Weighting factor for online data [dimensionless]
            * c["wf_CO2"]: Weighting factor for CO2 data [dimensionless]
            * c["pressure"]: Mean pressure in CO2 sensor [bar]

        :type c: dict
        :param yields: Calculated yields coefficients
        :type yields: dict 
                
        :return: cS, cX, cE, V, Fin, qsOx, qsRed, muTotal, qe, qCO2, qm: Substance concentrations and kinetic rates at time t. 
        """

        # (potential) fit parameters
        qsmax = p["qsmax"].value
        qemax = p["qemax"].value
        Ks    = p["Ks"].value
        Ke = p["Ke"].value
        Ki = p["Ki"].value
        
        qO2max = p["qO2max"].value
        qm_max = p["qm_max"].value
        g_e = p["g_e"].value
        YCO2s_m = yields["YCO2s_m"]
        
        Yxs_ox = yields["Yxs_ox"]
        Yxs_red = yields["Yxs_red"]
        Yxe_et = yields["Yxe_et"]
        Yxg_glyc = yields["Yxg_glyc"]  

        YCO2s_ox = yields["YCO2s_ox"]
        YCO2s_red = yields["YCO2s_red"]
        YCO2e_et = yields["YCO2e_et"]
        YCO2g_glyc = yields["YCO2g_glyc"]

        YO2s_ox = yields["YO2s_ox"]
        YO2e_et = yields["YO2e_et"]
        YO2g_glyc = yields["YO2g_glyc"]
        
        # controls
        feed_on = c["feed_on"] # time point when feed was switched on [h]
        feed_rate = c["feed_rate"] # feed rate [L/h]
        Fin = feed_rate * (t > feed_on) # becomes 0 if t < feed_on          ###muss vllt in rhs/derivatives function ####
       

        # masses and concentrations
        mS, mX, mE, V = y           ####vllt hier eher dict nehmen####
        cS, cX, cE = [mS, mX, mE] / V 

        #kinetics

        qs = qsmax * cS / (cS + Ks)

        if qs > qm_max:
            qm = qm_max
        else:
            qm = qs

        qs = qsmax * cS / (cS + Ks) - qm

        qO2_s_max = qs*YO2s_ox

        Monod_O2 = 1   #cO/(cO+Ko)

        qO2 = qO2max * Monod_O2


        if qO2_s_max <= qO2:
            
            qO2_s = qO2_s_max
            qsOx = qO2_s / YO2s_ox                  
            qsRed = 0
            
            qO2res = qO2 - qO2_s_max                       
            qe_Emax = qemax * cE/(cE + Ke) * Ki / (cS+ Ki)
            qO2_e = qe_Emax * YO2e_et
            qO2_g = g_e * qe_Emax * YO2g_glyc  

            if qO2_e + qO2_g <= qO2res:

                    qe = qe_Emax
            else:
                    qe = qO2res / (YO2e_et + g_e * YO2g_glyc)
                    
            qg = g_e*qe

        else:
            qO2_s = qO2 #before, there was qO2max
            qsOx = qO2_s / YO2s_ox         
            qsRed = qs - qsOx 
            qe = 0.0
            qg = 0.0

        #growth rates
        muOx = qsOx * Yxs_ox
        muRed = qsRed * Yxs_red
        muE = qe * Yxe_et
        muG = qg * Yxg_glyc

        muTotal = muOx+muRed+muE+muG
        qCO2 = qm*YCO2s_m + qsOx* YCO2s_ox + qsRed* YCO2s_red + qe* YCO2e_et + qg*YCO2g_glyc   #CO2 production rate

        return cS, cX, cE, V, Fin, qsOx, qsRed, muTotal, qe, qCO2, qm


    def model_rhs(self, t, y, p, c, yields):
        """ This function calculates the time derivatives to be integrated by solve_ivp by processing the outcome of the kinetics function.

        :param t: Current time [h]
        :type t: float
        :param y: State vector:

            * y[0] : Substrate mass (mS) [g]  
            * y[1] : Bio dry mass (mX) [g]
            * y[2] : Ethanol mass (mE) [g]
            * y[3] : Volume of fermentation broth [L]

        :type y: array like
        :param p: Structure with parameter values to be estimated, cf. lmfit.parameter.Parameters
        :type p: lmfit.parameter.Parameters object
        :param c:  Dictionary with control values

            * c["feed_on"]: Time point when feed was switched on [h]
            * c["feed_rate"]: Feed rate [L/h]
            * c["csf"]: Substrate concentration in feed [g/L]
            * c["M_base"]: Base concentration [mol/L]
            * c["gas_flow"]: Aeration rate [L/h]
            * c["T"]: Temperature in reactor [°C]
            * c["wf_on"]: Weighting factor for online data [dimensionless]
            * c["wf_CO2"]: Weighting factor for CO2 data [dimensionless]
            * c["pressure"]: Mean pressure in CO2 sensor [bar]

        :type c: dict
        :param yields: Calculated yields coefficients
        :type yields: dict 
                
        :return: dy_dt: Time derivative of state vector.
        """
        
        csf = c["csf"] # substrate concentration in feed [g/L]  
        Yes_red = yields["Yes_red"]          

        cS, cX, cE, V, Fin, qsOx, qsRed, muTotal, qe, qCO2, qm = self.kinetics(t, y, p, c, yields)   
    
        dmS_dt = cX *V*(- qsOx - qsRed -qm)+ csf *Fin    #-qm substraction ?
        dmX_dt = muTotal * cX * V
        dmE_dt = (qsRed * Yes_red - qe)*cX*V
        dV_dt  = + Fin

        return dmS_dt, dmX_dt, dmE_dt, dV_dt

    def calc_CO2(self, t, y, p, c, yields):
        """ This function calculates the CO2 pressure in the reactor in vol. % by processing the outcome of the kinetics function.
        
        :param t: Current time [h]
        :type t: float
        :param y: State vector:

            * y[0] : Substrate mass (mS) [g]  
            * y[1] : Bio dry mass (mX) [g]
            * y[2] : Ethanol mass (mE) [g]
            * y[3] : Volume of fermentation broth [L]

        :type y: array like
        :param p: Structure with parameter values to be estimated, cf. lmfit.parameter.Parameters
        :type p: lmfit.parameter.Parameters object
        :param c:  Dictionary with control values

            * c["feed_on"]: Time point when feed was switched on [h]
            * c["feed_rate"]: Feed rate [L/h]
            * c["csf"]: Substrate concentration in feed [g/L]
            * c["M_base"]: Base concentration [mol/L]
            * c["gas_flow"]: Aeration rate [L/h]
            * c["T"]: Temperature in reactor [°C]
            * c["wf_on"]: Weighting factor for online data [dimensionless]
            * c["wf_CO2"]: Weighting factor for CO2 data [dimensionless]
            * c["pressure"]: Mean pressure in CO2 sensor [bar]

        :type c: dict
        :param yields: Calculated yields coefficients
        :type yields: dict 
                
        :return: CO2_percent: CO2 pressure in vol. % 
    
        """

        #fix parameters 
        dV_gas_dt = c["gas_flow"]     # air flow in each experiment L/h
        R = 0.08314 #bar*l/mol*K
        T = c["T"] + 273.15     #Kelvin  
        pressure = c["pressure"] #bar
        M_CO2 = 44.01     #g/mol

        cS, cX, cE, V, Fin, qsOx, qsRed, muTotal, qe, qCO2, qm = self.kinetics(t, y, p, c, yields)

        #ideal gas quation to calculate CO2
        dmCO2_dt = qCO2 * cX * V
        dnCO2_dt = dmCO2_dt/M_CO2
        dvCO2_dt = (dnCO2_dt*R*T)/pressure
        CO2_percent = 100 * dvCO2_dt / dV_gas_dt

        return CO2_percent

    
    def observation(self, t_grid, y0, p, c, yields, **kwargs_solve_ivp):
        """ The observation function calculates simulated values by processing the outcome of the functions model_rhs and calc_CO2.

        :param t_grid: time grid with measurement time points [h]
        :type t_grid: array like
        :param y0: Initial state vector:

            * y[0] : Substrate mass (mS) [g]  
            * y[1] : Bio dry mass (mX) [g]
            * y[2] : Ethanol mass (mE) [g]
            * y[3] : Volume of fermentation broth [L]

        :type y0: array like
        :param p: Structure with parameter values to be estimated, cf. lmfit.parameter.Parameters
        :type p: lmfit.parameter.Parameters object
        :param c:  Dictionary with control values

            * c["feed_on"]: Time point when feed was switched on [h]
            * c["feed_rate"]: Feed rate [L/h]
            * c["csf"]: Substrate concentration in feed [g/L]
            * c["M_base"]: Base concentration [mol/L]
            * c["gas_flow"]: Aeration rate [L/h]
            * c["T"]: Temperature in reactor [°C]
            * c["wf_on"]: Weighting factor for online data [dimensionless]
            * c["wf_CO2"]: Weighting factor for CO2 data [dimensionless]
            * c["pressure"]: Mean pressure in CO2 sensor [bar]

        :type c: dict
        :param yields: Calculated yields coefficients
        :type yields: dict

        :param **kwargs_solve_ivp: optional arguments for scipy.integrate.solve_ivp

        :return: sim_exp: Simulated variables for experiment.
        
        
        """
        #result of IVP solver
        y_t = solve_ivp(self.model_rhs, [np.min(t_grid), np.max(t_grid)], y0, t_eval=t_grid, args = (p, c, yields), **kwargs_solve_ivp).y.T
    

        # unpack solution into vectors
        mS, mX, mE, V = [y_t[:,i] for i in range(4)]
        cS, cX, cE = [mS, mX, mE] / V

        dmX_dt = np.array([self.model_rhs(t = t_grid[i], y = y_t[i,:], p = p, c = c, yields = yields) for i in range(len(t_grid))])[:,1]  #calc dmX_dt from model funcs
        BASET_rate = p['base_coef'].value /c["M_base"] * dmX_dt * 1000  #base_rate in ml/h

        CO2_percent = np.array([self.calc_CO2(t = t_grid[i], y = y_t[i,:], p = p, c = c, yields = yields) for i in range(len(t_grid))])     #CO2 in vol.% from calc_CO2

        #simulated values
        sim_exp = pd.DataFrame(
            {"t": t_grid, "V": V,
            "cX": cX, "cS": cS, "cE" : cE, 
             "BASET_rate": BASET_rate,
            "CO2" : CO2_percent}
            ).set_index("t")

        return sim_exp
        

    
    def simulate(self, experiment = None, t_grid = None, y0 = None, p = None,  c = None, yields = None
    , kwargs_solve_ivp = dict(method= "Radau", first_step = 0.0000001, max_step= 0.1)       # kwargs_solve_ivp may have to be adapted
    , t_step = 1000, endpoint = None, **kwargs_solve_ivp_given):  

        """
        Simulates values by calling the observation function, it is possible to simulate manually by using values for t_grid, p, y0 c and yields or via an Experiment object. 
        If using an Experiment object, it is still possible to give manual values for  t_grid, p, y0 c. 
        Values not given are created from the Experiment object.

        :param experiment: Experiment object
        :type experiment: BioMoni.Experiment.Experiment

        :param t_grid: time grid with measurement time points [h]
        :type t_grid: array like
        :param y0: Initial state vector:

            * y[0] : Substrate mass (mS) [g] 
            * y[1] : Bio dry mass (mX) [g]
            * y[2] : Ethanol mass (mE) [g]
            * y[3] : Volume of fermentation broth [L]

        :type y0: array like
        :param p: Structure with parameter values to be estimated, cf. lmfit.parameter.Parameters
        :type p: lmfit.parameter.Parameters object
        :param c:  Dictionary with control values

            * c["feed_on"]: Time point when feed was switched on [h]
            * c["feed_rate"]: Feed rate [L/h]
            * c["csf"]: Substrate concentration in feed [g/L]
            * c["M_base"]: Base concentration [mol/L]
            * c["gas_flow"]: Aeration rate [L/h]
            * c["T"]: Temperature in reactor [°C]
            * c["wf_on"]: Weighting factor for online data [dimensionless]
            * c["wf_CO2"]: Weighting factor for CO2 data [dimensionless]
            * c["pressure"]: Mean pressure in CO2 sensor [bar]

        :type c: dict
        :param yields: Calculated yields coefficients
        :type yields: dict
        :param kwargs_solve_ivp: arguments for scipy.integrate.solve_ivp, either standard or coming from the estimate function
        :type kwargs_solve_ivp: dict
        
        :param t_step: steps for given time interval if t_grid is created from an Experiment object
        :type t_step: int, optional
        :param endpoint: Name of endpoint column in metadata, you may choose between different endpoints e.g. end1 and end2. Helpful if you have several endpoints in your metadata
        :type endpoint: str
        :param **kwargs_solve_ivp: optional arguments for scipy.integrate.solve_ivp
        
        :return: sim_exp: Simulated data.
        """
    
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
                    assert endpoint in experiment.metadata.index, "Given endpoint must be in metadata columns"   #.index because of pd.series

                if experiment.metadata[endpoint] is None: 

                    t_grid_list = [df.index.values for df in experiment.dataset.values()]   #get t_grid for every typ in experiment as list
                    t_grid_largest_endpoint = max(t_grid_list, key = np.max)       #get t_grid with largest endpoint
                    end = np.max(t_grid_largest_endpoint) #largest endpoint               
                    
                else:

                    end = (experiment.metadata[endpoint] - start_ts) / pd.Timedelta(1,"h") 

                t_grid = np.linspace(0, end, t_step) #timeline for simulation, going from 0 to end. Simulated values start always at zero in contrast to real measurement values
            
            if y0 is None:
                y0 = self.create_y0(experiment)
            
            if p is None:
                p = self.p

            if c is None:
                c = self.create_controls(experiment)
            
            if yields is None:
                yields = self.chemical_balancing(p)
                
            
        elif experiment is None:    #This if statement is to simulate with given t,y,p,c (used in the estimate process because it would require more computational power to calculate t, y and c within each estimate iteration)
            for i in [t_grid, p , y0, c, yields]:
                if i is None:
                    raise ValueError("If no Experiment is given, it is necessary to give t_grid, p , y0, c and yields in order to simulate")

        
        assert isinstance(t_grid, (np.ndarray, list)), "t_grid must be either of type np.ndarray or list"
        assert isinstance(y0, (np.ndarray, list)), "y0 must be either of type np.ndarray or list"
        assert isinstance(p, type(Parameters())), "Given parameters must be of type lmfit.parameter.Parameters"
        assert isinstance(c, dict), "c must be a dictionary"
        assert isinstance(yields, dict), "yields must be a dictionary"


        sim_exp = self.observation(t_grid, y0, p, c, yields, **kwargs_solve_ivp)

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

        yields = self.chemical_balancing(p)  #yields are calculated within the parameterestimation, thus influencing the parameter choice    

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
                sim_exp = self.simulate(None, t_grid, y0, p, c, yields, kwargs_solve_ivp = kwargs_solve_ivp)

                if tau is not None:     
                    weighting_decay = np.exp( (dat.index-np.max(dat.index) ) / tau)  #calculating exponential weighting decay 
                else:
                    weighting_decay = 1    #if no tau given: no additional weighting decay over time

                for var in dat:     #loop over measured variables (columns)
                    if var in sim_exp.columns:
                        res_var = wf * (sim_exp[var] - dat[var]).values * weighting_decay   # weighted residuals for this measured variable
                        res_single = np.append(res_single, res_var) # append to long residual vector
                    else:
                        pass
            res_all = np.append(res_all, res_single)        #super long vector with all residuals
        
        return res_all








        
            

