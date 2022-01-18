from numpy import result_type
from opcua import Client
from settings  import Result_path, url,sample_interval, data_name, nID_HM_4 as nID
import time
import pandas as pd
import os
import csv
import errno
from datetime import datetime


Result_path = Result_path + "/Results" + datetime.now().strftime('_%Y-%m-%d_%H-%M-%S')      
try:            #Path + subfolder will be created, if path already exists, only new subfolder within path will be created
    os.makedirs(Result_path)
except OSError as exc:
    if exc.errno != errno.EEXIST:
        raise
    pass


client = Client(url)
client.connect()
print("client connected")

#The Server may yield BASET_2 from Herr Menta 4, but i named the original column BASET here to match with the BASET from the model in the YEast class
data = {"PDatTime":[],"BASET":[],"CO2":[],"CO2_pressure":[]}      #create empty df in which the values will be appended
cols = data.keys()
data = pd.DataFrame(data)

#create csv file in which the data will be stored similar to the MFCS outcome, if you start the script again the file will be overwritten
cols= data.keys()    #header row : column names
first_row = ["Value"] * len(cols); second_row = ["Unit"] * len(cols)        #first row after column names  was Value (or Setpoint or Mode) in MFCS and the second row contained the Units, then there was a third empty row [] aswell
csv_name = os.path.join(Result_path, data_name) 
with open(csv_name, 'w', newline = "") as f:       #csv file will be created #newline because open makes some extra lines, avoid by using newline = ""
    writer = csv.writer(f, delimiter = ";")
    [writer.writerow(i) for i in (cols, first_row, second_row, [])]       #[] for 1 empty row

    
while True:
    
     #root
    root = client.get_root_node()   #root_node

    #get the value from the root_node through all the child nodes by usind the complete node id path stored in settings_mimic
    Value_PDatTime = root.get_child(nID["PDatTime"]).get_value()
    Value_BASET= root.get_child(nID["BASET_2"]).get_value()   #This would be the same as getting the value directly from the node id of the value node: client.get_node(values_id["BASET"][0]).get_value())
    Value_CO2= root.get_child(nID["CO2"]).get_value()
    CO2_pressure=root.get_child(nID["CO2_pressure"]).get_value()

    appendix = pd.Series(            
        {
        "PDatTime" : Value_PDatTime,
        "BASET" : Value_BASET,
        "CO2" : Value_CO2,
        "CO2_pressure" : CO2_pressure
        }
    )   ##this will be appended to the data at every while loop iteration
    
    appendix = pd.Series(appendix)

    #append csv file row by row
    with open(csv_name, 'a', newline = "") as f:
        writer = csv.writer(f, delimiter = ";")
        writer.writerow(appendix)

    #just to show the data as pd.DataFrame
    data = data.append(appendix, ignore_index=True) #append data row by row
    print(data)

    time.sleep(sample_interval)







    