#!/usr/bin/env python
# encoding: utf-8

import os
import csv
import numpy as np
from chemkin import ChemkinJob, getIgnitionDelay, getIgnitionDelayOH
import pandas as pd
from scipy.interpolate import CubicSpline


def find_nearest(array, value):
    array = np.asarray(np.abs(np.log(array)))
    value = np.abs(np.log(value))
    idx = (np.abs(array - value)).argmin()
    return idx,array[idx]

def fast_nearest_interp(xi, x, y):
    # Shift x points to centers
    spacing = np.diff(x) / 2
    x = x + np.hstack([spacing, spacing[-1]])
    # Append the last point in y twice for ease of use
    y = np.hstack([y, y[-1]])
    return y[np.searchsorted(x, xi)]

def ignitionDelay(df,pList, species,cond="max",specified_value ="None;None",exp_conc = "None"):
	if cond == "max":
		tau = df[species].idxmax()
	elif cond == "onset" and species != "p":
		time = df.index.to_numpy()
		conc = (df[species]).to_numpy()
		dtime = np.diff(df.index.to_numpy())
		dconc = np.diff((df[species]).to_numpy())
		slope = dconc/dtime
		intercept = int(np.diff(slope).argmax())
		tau = df.index.to_numpy()[intercept]
	elif cond == "onset" and species == "p":
		time = df.index.to_numpy()
		conc = np.asarray(pList)
		dtime = np.diff(df.index.to_numpy())
		dconc = np.diff(conc)
		slope = dconc/dtime
		intercept = int(np.diff(slope).argmax())
		tau = df.index.to_numpy()[intercept]
	elif cond == "dt-max" and species == "p":
		time = np.diff(df.index.to_numpy())
		conc = np.diff(np.asarray(pList))
		slope = conc/time
		index = int(slope.argmax())
		tau = df.index.to_numpy()[index]
	elif cond == "dt-max" and species != "p":
		time = np.diff(df.index.to_numpy())
		conc = np.diff(df[species].to_numpy())
		slope = conc/time
		index = int(slope.argmax())
		tau = df.index.to_numpy()[index]
	elif cond == "specific":
		if specified_value.split(";")[0] == None:
			raise Assertionerror("Input required for specified_value in ignition delay")
		else:
			target = float(specified_value.split(";")[0])
			unit = specified_value.split(";")[1]
			if unit == "molecule":
				avogadro = 6.02214E+23
				target = target/avogadro
				molecular_wt = gas.atomic_weight(species)
				target = molecular_wt*target	
				index,value = find_nearest(df[species],target)
				tau = df.index.to_numpy()[index[0]]	
			
			elif unit == "molecule/cm3":
				unit_conversion = 10**(6)
				time = df.index.to_numpy()
				conc = df[species].to_numpy()
				avogadro = 6.02214E+23
				target = (target/avogadro)*unit_conversion
				if exp_conc != "":
					exp_conc = (exp_conc/avogadro)*unit_conversion
					conc = conc*exp_conc
			
				f = CubicSpline(time,conc, bc_type='natural')
				time_new = np.arange(min(time),max(time),10**(-8))			
				conc_new = f(time_new)
				tau = fast_nearest_interp(target,conc_new,time_new)
				
			elif unit == "mole/cm3":
				time = df.index.to_numpy()
				conc = df[species].to_numpy()
				if exp_conc != "":
					conc = conc*exp_conc
				
				f = CubicSpline(time,conc, bc_type='natural')
				time_new = np.arange(min(time),max(time),10**(-8))
				conc_new = f(time_new)
				tau = fast_nearest_interp(target,conc_new,time_new)
	return tau


def convert_and_transpose_ckcsv(input_filepath, output_filepath):
    with open(input_filepath, 'r') as ckcsv_file:
        reader = csv.reader(ckcsv_file)

        # Skip the first row (metadata) and read the actual data
        rows = list(reader)
        data_rows = rows[1:-2]  # Skip metadata row
        
        # Extract variables, units, and time series data
        variables = [row[0] for row in data_rows]     # First column as variable names
        units = [row[1] for row in data_rows]         # Second column as units
        data_values = [row[2:] for row in data_rows]  # Time series data from the third column onward
        headers = [f"{row[0]} {row[1]}" for row in data_rows]
        #print(data_values[-1])    
        # Transpose the data values so each variable has its data in columns
        transposed_data = list(zip(*data_values))
        with open(output_filepath, 'w', newline='') as csv_file:
            writer = csv.writer(csv_file)

            # Write headers (variable names)
            writer.writerow(headers)

            # Write units row
            #writer.writerow(units)

            # Write transposed data rows
            writer.writerows(transposed_data)

    print(f"Conversion and transposition complete. Data saved to {output_filepath}")

#mechanism in yaml format: /data/TEST-THERMO-sens/FFCM_THERMO_UQ_OPT/PRS_INVESTIGATION_ID_1029_TIG_CHEMKIN/Perturbed_Mech/mechanism_0.yaml

# Parameters
reactorTemperature = 1944.0  # Replace with your desired temperature in Kelvin
reactorPressure = 112000.00000000001  # Replace with your desired pressure in atm
spec = {'H2': 0.0017, 'CO': 0.033, 'O2': 0.0347, 'AR': 0.9306}
conc_X = list(spec.items()) # Example species concentrations, adjust as needed

# Set up the Chemkin simulation
currentDir = os.path.dirname(__file__).strip("\n")
chemFile = os.path.join(currentDir, 'mechanism.inp')  # Update with the mechanism file path
tempDir = os.path.join(currentDir, 'output')


# Initialize Chemkin job
job = ChemkinJob(
    name=conc_X[0][0],
    chemFile=chemFile,
    tempDir=tempDir,
)

job.preprocess()

# Write Chemkin input for a single simulation
input_file = job.writeInputHomogeneousBatch(
    problemType='constrainVandSolveE',  # Problem type in Chemkin
    reactants=conc_X,
    temperature=reactorTemperature,
    pressure=reactorPressure,
    endTime= 0.04  # Simulation end time in seconds
)

# Run the Chemkin simulation
job.run(input_file, model='CKReactorGenericClosed', pro=True)
job.postprocess(sens=False, rop=False, all=True, transpose=False)

#tau = ignitionDelay(timeHistory,pressureList,"OH*","max","None",)

# Calculate ignition delay
try:
    tau = getIgnitionDelayOH(job.ckcsvFile)
except ValueError:
    print("Value Error")
    tau = None


#if saveAll:
#    convert_and_transpose_ckcsv(job.ckcsvFile, "time_history.csv")

#timeHistory = pd.read_csv("time_history.csv")
#pressureList = []
#tau = ignitionDelay(timeHistory,pressureList,"OH*","max","None",)

# Output the results
if tau:
    print(f"Ignition delay at {reactorTemperature} K and {reactorPressure} atm: {tau*1e3} ms")
    with open("output/tau.out", 'w') as tau_file:
        tau_file.write(f"#T(K)    tau(us)\n{reactorTemperature}    {tau * 1e6}\n")

# Save all output data if needed
saveAll = True  # Set to True if you want to save all time-history data
if saveAll:
    convert_and_transpose_ckcsv(job.ckcsvFile, "time_history.csv")
    #output_data.to_csv("time_history.csv", index=False)

import glob

def delete_files(directory):
    # Patterns for files to delete
    patterns = ["*.out", "*.asc", "*.zip"]
    exception_file = "tau.out"
    
    for pattern in patterns:
        # Get all files matching the pattern in the directory
        files = glob.glob(os.path.join(directory, pattern))
        
        for file in files:
            # Skip the exception file
            if os.path.basename(file).lower() == exception_file.lower():
                continue
            
            # Delete the file
            try:
                os.remove(file)
                print(f"Deleted: {file}")
            except Exception as e:
                print(f"Error deleting {file}: {e}")

# Specify the directory containing the files
directory_path = "output/"
delete_files(directory_path)
