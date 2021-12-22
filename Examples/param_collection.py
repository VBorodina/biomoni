"""
Just a collection of different parameter values that can be loaded into the respective modules by the import function and passed to the yeast model object by the set_params function.

"""

from lmfit import Parameters


#Parameter for MFCS mimick simulation
p0 = Parameters()
p0.add("qsmax", value= 1.61290497, min=0.01, max=5.0 , vary = True)     #Maximum glucose uptake rate [g/(g*h)]   #previous start val: 0.1  #estimated value: 1.61290497
p0.add("qemax", value=0.2361, min=0.15, max=0.35, vary=False)    #Maximum ethanol uptake rate [g/(g*h)]       #value=0.2361, min=0.15, max=0.35, vary=True # value=0.05, min=0.001, max=0.5,  vary=True geht

p0.add("base_coef", value= 0.00733798 , min=0.0001, max =3, vary = True)     #base coefficient [mol/g]     #previous start val: 1 #estimated value: 0.00733798      
p0.add("qO2max", value= 0.16473635, min=0.1, max=0.4, vary = False)    #Maximum oxygen uptake rate [g/(g*h)]  , min=0.1  0.1, 0.4  #previous start val: 0.255984, 1 #estimated value: 0.17451356
p0.add("qm_max", value=0.01, min=0.0075, max=0.0125, vary = False)  #glucose uptake rate required for maintenance [g/(g*h)]  

p0.add("Ks", value=0.1, min=0.01, max=1,  vary=False)   #Saturation constant, concentration of glucose at µ = 0.5 µmax [g/L]
p0.add("Ke", value=0.1, min=0.01, max=1, vary=False)    #Saturation constant, concentration of ethanol at µ = 0.5 µmax [g/L]
p0.add("Ki", value=0.1, min=0.01, max=1, vary=False)    #Inhibition constant, glucose inhibits uptake of ethanol [g/L]

p0.add("Yxs_ox", value= 0.52792125, min=0.3, max=0.6, vary= False) #Yield biomass per glucose (oxidative growth) [g/g] #previous start val: 49    #estimated value 0.52792125
p0.add("Yxs_red", value=0.05, min=0.01, max=0.8, vary=False)    #Yield biomass per glucose (reductive growth) [g/g] 
p0.add("Yxe_et", value=0.72, min=0.5, max=0.8, vary=False)      #Yield biomass per ethanol [g/g]
p0.add("Yxg_glyc", value=0.2, min=0.1, max=0.35, vary=False)    #Yield biomass per glycerol [g/g]              #value=0.2, min=0.1, max=0.35, vary=False

#Biomass composition paramaters, content H,O,N from paper sonnleitner

p0.add("HX", value=1.79, min=1.77, max=2.1, vary=False) #Stoichiometric hydrogen content of biomass [mol/mol]
p0.add("OX", value=0.57, min=0.54, max=0.63, vary=False) #Stoichiometric oxygen content of biomass [mol/mol]     
p0.add("NX", value=0.15, min=0.14, max=0.16, vary=False) #Stoichiometric nitrogen content of biomass [mol/mol]

p0.add("g_e", value=0.21636282532087284, min=0.17624196916440354, max=0.2883951412907231, vary=False) #determined experimentally: formation glycerol per ethanol [g/g]


#Parameters to start with in laborloop.py
p1 = Parameters()
p1.add("qsmax", value= 1.63, min= 1.5, max=1.7 , vary = True)     #Maximum glucose uptake rate [g/(g*h)]   #previous start val: 0.1  #estimated value: 1.61290497
p1.add("qemax", value=0.2361, min=0.15, max=0.35, vary=False)    #Maximum ethanol uptake rate [g/(g*h)]       #value=0.2361, min=0.15, max=0.35, vary=True # value=0.05, min=0.001, max=0.5,  vary=True geht

p1.add("base_coef", value= 0.01 , min=0.0001, max =3, vary = True)     #base coefficient [mol/g]     #previous start val: 1 #estimated value: 0.00733798      
p1.add("qO2max", value= 0.16473635, min=0.1, max=0.4, vary = True)    #Maximum oxygen uptake rate [g/(g*h)]  , min=0.1  0.1, 0.4  #previous start val: 0.255984, 1 #estimated value: 0.17451356
p1.add("qm_max", value=0.01, min=0.0075, max=0.0125, vary = False)  #glucose uptake rate required for maintenance [g/(g*h)]  

p1.add("Ks", value=0.1, min=0.01, max=1,  vary=False)   #Saturation constant, concentration of glucose at µ = 0.5 µmax [g/L]
p1.add("Ke", value=0.1, min=0.01, max=1, vary=False)    #Saturation constant, concentration of ethanol at µ = 0.5 µmax [g/L]
p1.add("Ki", value=0.1, min=0.01, max=1, vary=False)    #Inhibition constant, glucose inhibits uptake of ethanol [g/L]

p1.add("Yxs_ox", value= 0.52792125, min=0.3, max=0.6, vary= True) #Yield biomass per glucose (oxidative growth) [g/g] #previous start val: 49    #estimated value 0.52792125
p1.add("Yxs_red", value=0.05, min=0.01, max=0.8, vary=False)    #Yield biomass per glucose (reductive growth) [g/g] 
p1.add("Yxe_et", value=0.72, min=0.5, max=0.8, vary=False)      #Yield biomass per ethanol [g/g]
p1.add("Yxg_glyc", value=0.2, min=0.1, max=0.35, vary=False)    #Yield biomass per glycerol [g/g]              #value=0.2, min=0.1, max=0.35, vary=False

#Biomass composition paramaters, content H,O,N from paper sonnleitner

p1.add("HX", value=1.79, min=1.77, max=2.1, vary=False) #Stoichiometric hydrogen content of biomass [mol/mol]
p1.add("OX", value=0.57, min=0.54, max=0.63, vary=False) #Stoichiometric oxygen content of biomass [mol/mol]     
p1.add("NX", value=0.15, min=0.14, max=0.16, vary=False) #Stoichiometric nitrogen content of biomass [mol/mol]

p1.add("g_e", value=0.21636282532087284, min=0.17624196916440354, max=0.2883951412907231, vary=False) #determined experimentally: formation glycerol per ethanol [g/g]


#Old param configurations from PAul Senck's Ma-Thesis configuration, not necessarily needed for anything, only saved in case something is no longer comprehensible.
p_old = Parameters()
p_old.add("qsmax", value=0.1, min=0.01, max=5.0 , vary = True)     #Maximum glucose uptake rate [g/(g*h)]
p_old.add("qemax", value=0.2361, min=0.15, max=0.35, vary=False)    #Maximum ethanol uptake rate [g/(g*h)]       

p_old.add("base_coef", value= 1, min=0.0001, max =3, vary = True)     #base coefficient [mol/g]      
p_old.add("qO2max", value= 0.255984, min=0.1, max=0.4, vary = True)    #Maximum oxygen uptake rate [g/(g*h)]  
p_old.add("qm_max", value=0.01, min=0.0075, max=0.0125, vary = False)  #glucose uptake rate required for maintenance [g/(g*h)]

p_old.add("Ks", value=0.1, min=0.01, max=1,  vary=False)   #Saturation constant, concentration of glucose at µ = 0.5 µmax [g/L]
p_old.add("Ke", value=0.1, min=0.01, max=1, vary=False)    #Saturation constant, concentration of ethanol at µ = 0.5 µmax [g/L]
p_old.add("Ki", value=0.1, min=0.01, max=1, vary=False)    #Inhibition constant, glucose inhibits uptake of ethanol [g/L]

p_old.add("Yxs_ox", value=0.49, min=0.3, max=0.6, vary=True) #Yield biomass per glucose (oxidative growth) [g/g]  
p_old.add("Yxs_red", value=0.05, min=0.01, max=0.8, vary=False)    #Yield biomass per glucose (reductive growth) [g/g] 
p_old.add("Yxe_et", value=0.72, min=0.5, max=0.8, vary=False)      #Yield biomass per ethanol [g/g]
p_old.add("Yxg_glyc", value=0.2, min=0.1, max=0.35, vary=False)    #Yield biomass per glycerol [g/g]             

#Biomass composition paramaters, content H,O,N from paper sonnleitner

p_old.add("HX", value=1.79, min=1.77, max=2.1, vary=False) #Stoichiometric hydrogen content of biomass [mol/mol]
p_old.add("OX", value=0.57, min=0.54, max=0.63, vary=False) #Stoichiometric oxygen content of biomass [mol/mol]     
p_old.add("NX", value=0.15, min=0.14, max=0.16, vary=False) #Stoichiometric nitrogen content of biomass [mol/mol]

p_old.add("g_e", value=0.21636282532087284, min=0.17624196916440354, max=0.2883951412907231, vary=False) #determined experimentally: formation glycerol per ethanol [g/g]



