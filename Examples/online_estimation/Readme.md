# How to use the biomoni package in terms of online an estimation?
In this example, data from a Sartorius MFCS unit is automatically simulated and released via an OPC-UA server.
To do this, the following steps must be taken:
1. `open param_collection` and enter for `p0` the parameter vector for which the simulated data should be generated. `p1` is the parameter vector that will be used later just before the estimation in `laborloop`.
2. `open settings` and adjust the following:
    • `url` = OPC-UA Server url
    • `Simulation_path` = path in which the experiment to be simulated is located
    • `exp_id` = identifier of the Experiment to be simulated
    • `Result_path` = path where the results will be stored.
    • `data_name` = name of the csv file into which the data will be written.
    • `sample_interval` = in seconds
    • All the Node ids
3. Execute the module `MFCS_mmic` to simulate MFCS data as OPC-UA Server
4. Execute the module `Data_collector` which is an OPC-UA client and pulls data from the OPC-UA server and stores it as a csv file in the Results_path.
5. Start the module `labloop`, be sure that the metadata for the simulation contains the same values as the metadata for the estimation.
