#Settings required to build an OPC server, and an OPC client

url = "opc.tcp://127.0.0.1:4848"    #OPC server url

##windows
RESULT_PATH = r"P:\Code\Messdaten\OPCUA"    #path where Results (data from MFCS) are stored

#Linux
RESULT_PATH = "/home/paul/pCloudDrive/Code/Messdaten/OPCUA" #path where Results (data from MFCS) are stored



data_name = "data_1"        #name of the csv file in which data will be saved

sample_interval = 60    #Time in sec in which a measured value is generated (server) or a measured value is read (client). Measurement value generation and reading is equal is this case because otherwise timestamps would not fit the measurement values in a logical manner and would thus yield wrong derivatives when calculating the base_rate for example.

#Server settings
#Units identity
units_id = dict(
FM_1 = ('ns=2;s="FM_1"', "FM_1"),
FM_2 = ('ns=2;s="FM_2"', "FM_2"),
HM_3 = ('ns=2;s="HM_1"', "HM_1"),
HM_4 = ('ns=2;s="HM_4"', "HM_4")
)

#Variables_id
variables_id = dict(
BASET = ('ns=2;s="BASET_2"', "BASET_2"),
CO2 = ('ns=2;s="CO2"', "CO2"),
CO2_pressure = ('ns=2;s="CO2_p"', "CO2_pressure"),
ts = ('ns=2;s="PDatTime"', "PDatTime")
)

#Value_id
values_id = dict(
BASET = ('ns=2;s="Value_BASET_2"',"Value",0),        #default value is zero: gets changed anyway
CO2 = ('ns=2;s="Value_CO2"',"Value",0),
CO2_pressure = ('ns=2;s="CO2_pressure"',"Value",1.058),
ts = ('ns=2;s="Value_PDatTime"',"Value",0)
)



#Root identity: this is useful to get the values directly from the root node through root.get_child. This is used in the data collector.
#Martins Anmerkung: überall statt 0: 2: um die Anpassung von python UA Server auf OPCUA-server von MFCS zu machen
#Meine Anmerkung: überall steht Value+Variable nur bei BASET_2, CO2, CO2_pressure und PDatTime steht nur Value....

nID_HM_4 = {
    'ACIDT': ["0:Objects","0:Units","0:HM_4","0:Variables","0:BALANCE_B1","0:Value_BALANCE_B1"],
    'BALANCE_A1': ["0:Objects","0:Units","0:HM_4","0:Variables","0:BALANCE_A1","0:Value_BALANCE_A1"],
    'BALANCE_B1':["0:Objects","0:Units","0:HM_4","0:Variables","0:BALANCE_B1","0:Value_BALANCE_B1"],
    'BASET_2':["0:Objects","0:Units","0:HM_4","0:Variables","0:BASET_2","0:Value"],
    'CO2':["0:Objects","0:Units","0:HM_4","0:Variables","0:CO2","0:Value"],     #No idea why Martin named it only Value here and not Value_Variable
    'CO2_pressure':["0:Objects","0:Units","0:HM_4","0:Variables","0:CO2_pressure","0:Value"],
    'ECO2_1':["0:Objects","0:Units","0:HM_4","0:Variables","0:ECO2_1","0:Value_ECO2_1"],
    'EO2_1':["0:Objects","0:Units","0:HM_4","0:Variables","0:EO2_1","0:Value_EO2_1"],
    'EXT_A1':["0:Objects","0:Units","0:HM_4","0:Variables","0:EXT_A1","0:Value_EXT_A1"],
    'EXT_B1':["0:Objects","0:Units","0:HM_4","0:Variables","0:EXT_B1","0:Value_EXT_B1"],
    'FOAMT_1':["0:Objects","0:Units","0:HM_4","0:Variables","0:FOAMT_1","0:Value_FOAMT_1"],
    'MFC_B1':["0:Objects","0:Units","0:HM_4","0:Variables","0:MFC_B1","0:Value_MFC_B1"],
    'MFC_C1':["0:Objects","0:Units","0:HM_4","0:Variables","0:MFC_C1","0:Value_MFC_C1"],
    'MFC_D1':["0:Objects","0:Units","0:HM_4","0:Variables","0:MFC_D1","0:Value_MFC_D1"],
    'O2SPT_1':["0:Objects","0:Units","0:HM_4","0:Variables","0:O2SPT_1","0:Value_O2SPT_1"],
    'PUMPT_A1':["0:Objects","0:Units","0:HM_4","0:Variables","0:PUMPT_A1","0:Value_PUMPT_A1"],
    'REDOX_1':["0:Objects","0:Units","0:HM_4","0:Variables","0:REDOX_1","0:Value_REDOX_1"],
    'SUBST_A1':["0:Objects","0:Units","0:HM_4","0:Variables","0:SUBST_A1","0:Value_SUBST_A1"],
    'SUBST_B1':["0:Objects","0:Units","0:HM_4","0:Variables","0:SUBST_B1","0:Value_SUBST_B1"],
    'SUBST_C1':["0:Objects","0:Units","0:HM_4","0:Variables","0:SUBST_C1","0:Value_SUBST_C1"],
    'TURB_1':["0:Objects","0:Units","0:HM_4","0:Variables","0:TURB_1","0:Value_TURB_1"],
    'VALVET_A1':["0:Objects","0:Units","0:HM_4","0:Variables","0:VALVET_A1","0:Value_VALVET_A1"],
    'VALVET_C1':["0:Objects","0:Units","0:HM_4","0:Variables","0:VALVET_C1","0:Value_VALVET_C1"],
    'VALVET_D1':["0:Objects","0:Units","0:HM_4","0:Variables","0:VALVET_D1","0:Value_VALVET_D1"],
    'VALVE_A1':["0:Objects","0:Units","0:HM_4","0:Variables","0:VALVE_A1","0:Value_VALVE_A1"],
    'VALVE_C1':["0:Objects","0:Units","0:HM_4","0:Variables","0:VALVE_C1","0:Value_VALVE_C1"],
    'VALVE_D1':["0:Objects","0:Units","0:HM_4","0:Variables","0:VALVE_D1","0:Value_VALVE_D1"],
    'PDatTime':["0:Objects","0:Units","0:HM_4","0:Variables","0:PDatTime","0:Value"],           
    'pHo_1':["0:Objects","0:Units","0:HM_4","0:Variables","0:pHo_1","0:Value_pHo_1"],
    'pO2o_1':["0:Objects","0:Units","0:HM_4","0:Variables","0:pO2o_1","0:Value_pO2o_1"],
}
