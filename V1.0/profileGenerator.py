import re, os, sys
import subprocess
import numpy as np
import shutil
import yaml
from itertools import chain
def unitConverter(var,quant,unit):
	if var == "T":
		if unit == "K":
			output = quant
		elif unit == "mK":
			output = quant*0.001
		elif unit == "C":
			output = quant+273.15#0 C = 273.15 K
		else:
			raise AssertionError("Unit of temperature is {} instead of K".format(unit))
			exit()
	elif var == "P":
		if unit == "bar":
			output = quant*101350
		elif unit == "atm":
			output = quant*101350
		elif unit == "Pa":
			output = quant
		elif unit == "MPa":
			output = quant*1000000
		else:
			raise AssertionError("Unit of pressure is {} instead of Pa".format(unit))
			exit()
	return output

def dirichlet_str(species):
	str_syntax = "x->"
	dirichlet = ""
	for i in species:
		dirichlet += str_syntax+str(i)+"="+str(species[i])+"\n\t\t"
	return dirichlet
	
def concBounds_str(bound_type,bounds,exp_type="prem",inc=0.0005,inc_t = 25,To = None,From = None):
	if bound_type == "species":
		concBounds ="\nContBound = {}\nContInc = {}\nToSpecies is {}\nFromSpecies is {}".format(bounds,inc,To,From)
	elif bound_type == "temperature":
		concBounds = "\nContInc = {}\nContType is Temperature\nContBound = {}".format(inc_t,bounds)
	else:
		raise AssertionError("Invalid type used for concBounds")
	return concBounds
	
def dict_screening(imp,collection):
	dict1 = imp
	dict2 = collection
	keys = []
	for i in imp:
		keys.append(i)	
	for j in collection:
		keys.append(j)
	keys = set(keys)
	
	for i in keys:
		if i in imp:
			continue
		else:
			dict1[str(i)] = 0.0
	for i in keys:	
		if i in collection:
			continue
		else:
			dict2[str(i)] = 0.0
	return dict1,dict2		
def compare(fileA,fileB):
	g = open(fileA,"r").readlines()
	f = open(fileB,"r").readlines()
	for i in g:
		if i not in f:
			close = False
		else:
			close = True
	return close
	
def getNewProfile(inputs,phi,fuel,oxidizer,bathgas,globeRxn,P,T,controls):
	main.create_inputs_to_start(inputs,phi,fuel,oxidizer,bathgas,globeRxn,P,T,controls)
	subprocess.call("./run_generate")
	start_profile = None
	
	home_dir = os.getcwd()
	os.chdir("output")
	
	#if "previousProfile" in inputs:	
		#if inputs["previousProfile"] != inputs["initialProfile"]:
		#	os.remove(inputs["previousProfile"])
	l = list(os.listdir())
	if fuel["key"] == "Mono-Multi":
		available_fuel = fuel["From"]
		target_fuel = fuel["To"]["a"]
	elif fuel["key"] == "Multi-Multi":
		available_fuel = fuel["From"]["a"]
		target_fuel = fuel["To"]["a"]

	else:
		available_fuel = fuel["From"]
		target_fuel = fuel["To"]
	
	for files in l:
		if available_fuel in files or target_fuel in files:
			if files not in inputs["previousProfile"]:
				if files.endswith("noC"):
					#for i in os.listdir():
					#	os.remove(i)
					inputs["TimeFlag"] = "TRUE"
					os.chdir("..")
					main.create_inputs_to_start(inputs,phi,fuel,oxidizer,bathgas,globeRxn,P,T,controls)
					subprocess.call("./run_generate")
					os.chdir("output")
					l = list(os.listdir())
					for files in l:
						if available_fuel in files or target_fuel in files:
							if files.endswith("noC"):
								fileA = str(files).strip("noC")+"_t"
								if fileA in os.listdir():
									start_profile = os.path.abspath(fileA)
									os.remove(files)
									break
								else:
									raise AssertionError("startprofilegenerator failed!!!")
							#elif files in inputs["StartProfilesFile"]:
							#	bool_ = compare(inputs["StartProfilesFile"],os.path.abspath(files))
							#	if bool_ == True:
							#		continue
							#	else:
							#		start_profile = os.path.abspath(files)
							#		break
							#else:
							#	start_profile = os.path.abspath(files)
							#	break
				#elif files in inputs["StartProfilesFile"]:
				#	bool_ = compare(inputs["StartProfilesFile"],os.path.abspath(files))
				#	start_profile = os.path.abspath(files)
				#	break
				else:
					try:
						shutil.copy(files,home_dir)
					except:
						raise AssertionError("File already exist")
					start_profile = str(home_dir)+"/"+str(files)
					for i in os.listdir():
						os.remove(i)
					break
			else:
				continue
	#shutil.move(start_profile,home_dir)
	os.chdir(home_dir)
	return start_profile

def readjustment(spec,FromSpecies,ToSpecies):
	newSpecies = {}
	keys = list(FromSpecies.keys())
	values = []
	for i in FromSpecies:
		values.append(ToSpecies[i] - FromSpecies[i])
	max_inc = values.index(max(values))
	max_dec = values.index(min(values))
	 
	diff = 0.0
	for ind,key in enumerate(keys):
		if spec != key and ToSpecies[key] == FromSpecies[key]:
			newSpecies[key] = float(ToSpecies[key])
		elif spec == key:
			newSpecies[key] = float(ToSpecies[key])
			diff = newSpecies[key] - float(FromSpecies[key])
		elif spec !=key and diff<0 and ind == max_inc:
			newSpecies[keys[max_inc]] = float(FromSpecies[keys[max_inc]])+abs(diff)
		elif spec !=key and diff>0 and ind == max_dec:
			newSpecies[keys[max_dec]] = float(FromSpecies[keys[max_dec]])-abs(diff)
		
		else:
			newSpecies[key] = float(FromSpecies[key])
							
	#print(newSpecies)
	return newSpecies
					
if len(sys.argv) > 1:
	input_file = open(sys.argv[1],'r')
	text = input_file.read()
	print("Input file found\n")
else:
	print("Please enter a valid input file name as arguement. \n For details of preparing the input file, please see the UserManual\n\nProgram exiting")
	exit()


inputs = yaml.safe_load(text)["Inputs"]
profile = yaml.safe_load(text)["Boundary"]

binaryLoc = inputs["bin"]
flameis = inputs["Flame"]
#print(flameis)
mechLoc = inputs["MechanismFile"]
startprofileLoc = inputs["StartProfilesFile"]
inputs["initialProfile"] = startprofileLoc
if "ExpTempFile" in inputs :
	ExpTempFile = inputs["ExpTempFile"]
	if ExpTempFile == None:
		ExpTempFile = ""
		inputs["ExpTempFile"] = ExpTempFile
		
else:
	ExpTempFile = ""
	inputs["ExpTempFile"] = ExpTempFile

if "ComputeWithRadiation" in inputs:
	CoWin = inputs["ComputeWithRadiation"]
else:
	CoWin = "True"
	inputs["ComputeWithRadiation"] = CoWin

if "Thermodiffusion" in inputs:
	TheFusion = inputs["Thermodiffusion"]
else:
	TheFusion = "True"
	inputs["Thermodiffusion"] = TheFusion

if "phi" in profile:
	phi = profile["phi"]
else:
	phi = "#phi = "
	

fuel = profile["fuel"]
fuelType = fuel["key"]
oxidizer = profile["oxidizer"]
bathgas = profile["bathGas"]
globRxn = profile["globalReaction"]
pressureTo = profile["pressure"]
temperatureTo = profile["temperature"]
units = profile["units"]

#step: 1
#start the process of making flamemaster input files
sys.path.append(binaryLoc)
import make_input_file as main

if "UnstretchedPremixed" in flameis:
	if ExpTempFile == "":
		exp = "Fls"
		expType = "Laminar premixed"
	else:
		exp = "Flf"
		expType = "HeatFlux"
	inputs["exp_type"] = {"exp":exp,"type":expType,"sudo":"Fls"}

#else part will be developed later
#print(exp)
initialProfile = open(startprofileLoc).readlines()

for line in initialProfile:
	if "title" in line:
		start_profile_type = line.split("=")[1]
		continue
	if "mechanism" in line:
		start_profile_mechanism = line.split("=")[1]
		start_profile_mechanism = start_profile_mechanism.replace('"', '')
		continue
	if "author" in line:
		start_profile_author = line.split("=")[1]
		continue
	if "date" in line:
		start_profile_date = line.split("=")[1]
		continue
	if "ConstantLewisNumbers" in line:
		if "UnstretchedPremixed" in flameis:
			exp = "Flf"
			expType = "HeatFlux"
		else:
			exp = "Fls"
			expType = "Laminar premixed"
		inputs["exp_type"] = {"exp":exp,"type":expType,"sudo":"Fls"}
	if fuelType == "Mono-Multi":
		if "fuel" in line and "air" not in line:
			start_profile_fuel_name = line.split("=")[1]
			start_profile_fuel_name = start_profile_fuel_name.replace('"', '')
			continue
	
	elif fuelType == "Multi-Multi":
		start_profile_fuel_name = {}
		for i in fuel["From"]:
			if "fuel" in line and "air" not in line and fuel["From"][str(i)] in line:
				name = line.split("=")[1]
				start_profile_fuel_name[str(i)] = name.replace('"', '')
	else:
		if "fuel" in line and "air" not in line:
			start_profile_fuel_name = line.split("=")[1]
			start_profile_fuel_name = start_profile_fuel_name.replace('"', '')
			continue
							
	if "pressure" in line:
		start_profile_pressure = line.split("=")[1]
		if "bar" in start_profile_pressure.split("[")[1]:
			start_profile_pressure = float(start_profile_pressure.split("[")[0])*10**5
			continue
		elif "atm" in start_profile_pressure.split("[")[1]:
			start_profile_pressure = float(start_profile_pressure.split("[")[0])*10**5
			continue
		elif "Pa" in start_profile_pressure.split("[")[1]:
			start_profile_pressure = float(start_profile_pressure.split("[")[0])
			continue
		else:
			start_profile_pressure = float(start_profile_pressure.split("[")[0])*10**5
		continue
	if "fuel-air-equivalence-ratio" in line:
		start_profile_phi = line.split("=")[1]
		continue
	if "Tmax" in line:
		start_profile_Tmax = line.split("=")[1]
		continue
	if "T10mm" in line:
		start_profile_T10mm = line.split("=")[1]
		continue
	if "unburnt" in line:
		start = initialProfile.index(line)+1
		for j in initialProfile[initialProfile.index(line)+1:]:        
			if "end" in j:
				stop = initialProfile.index(j)+1
		boundary_string = " ".join(list(initialProfile[start:stop]))
		continue
	if "burningVelocity" in line:
		start_profile_fls = line.split("=")[1]
		continue
	if "FlameThickness" in line:
		start_profile_flt = line.split("=")[1]
		continue
	if "numOfSpecies" in line:
		start_profile_species_count = line.split("=")[1]
		continue
	if "gridPoints" in line:
		start_profile_gridPoints = line.split("=")[1]
	
		continue
if str(start_profile_fuel_name.strip()) != str(fuel["From"]):
	print("Error in start profile fuel name")
	raise AssertionError(start_profile_fuel_name)


for line in boundary_string.split("\n"):
	if bathgas["key"] == "Multi-Multi":
		bathgas["FromConc"] = {}
		for gas in bathgas["From"]:
			if "Molefraction-"+str(bathgas["From"][gas]).strip() in line:
				bathgas["FromConc"][gas] = float(line.split("=")[1])
				continue
	else:
		if "Molefraction-"+str(bathgas["From"]).strip() in line:
				bathgas["FromConc"] = float(line.split("=")[1])
				continue
	if fuel["key"] == "Multi-Multi":
		for gas in fuel["From"]:
			if "Molefraction-"+str(fuel["From"][gas]).strip() in line:
				fuel["FromConc"][gas] = float(line.split("=")[1])
				continue
	else:
		if "Molefraction-"+str(fuel["From"]).strip() in line:
				fuel["FromConc"] = float(line.split("=")[1])
				continue
			
	if "Molefraction-"+str(oxidizer["From"]) in line:
		oxidizer["FromConc"] = float(line.split("=")[1])
		continue
		
	elif "Temperature" in line:
		init_Temperature = float(line.split("=")[1].split("[")[0])
	elif "begin" in line:
		continue
	elif "end" in line:
		break
	else:
		continue

inputs["ContInc"] = 0.0005
ComputeWithRadiation = CoWin
Thermodiffusion = TheFusion
inputs["maxgridpoints"] = 300	
inputs["TimeFlag"] = "FALSE" #initially
inputs["initialgridpoints"] = start_profile_gridPoints
inputs["flowRate"] = profile["flowRate"]
species = {}
speciesTo = {}

if fuel["key"] == "Multi-Multi":
	for i in fuel["From"]:
		species[fuel["From"][i]] = fuel["FromConc"][i] 
	for i in fuel["To"]:
		speciesTo[fuel["To"][i]] = fuel["ToConc"][i] 
elif fuel["key"] == "Mono-Multi":
	species[fuel["From"]] = fuel["FromConc"] 
	for i in fuel["To"]:
		speciesTo[fuel["To"][i]] = fuel["ToConc"][i] 
elif fuel["key"] == "Multi-Mono":
	for i in fuel["From"]:
		species[fuel["From"][i]] = fuel["FromConc"][i] 
	speciesTo[fuel["To"]] = fuel["ToConc"]  
else:
	species[fuel["From"]] = fuel["FromConc"] 
	speciesTo[fuel["To"]] = fuel["ToConc"] 
	
species[oxidizer["From"]] = oxidizer["FromConc"] 
speciesTo[oxidizer["To"]] = oxidizer["ToConc"] 

if bathgas["key"] == "Multi-Multi":
	for i in bathgas["From"]:	
		species[bathgas["From"][i]] = bathgas["FromConc"][i] 
	for i in bathgas["To"]:	
		speciesTo[bathgas["To"][i]] = bathgas["ToConc"][i] 
elif bathgas["key"] == "Mono-Multi":
	species[bathgas["From"]] = bathgas["FromConc"] 
	for i in bathgas["To"]:	
		speciesTo[bathgas["To"][i]] = bathgas["ToConc"][i] 
elif bathgas["key"] == "Multi-Mono":
	for i in bathgas["From"]:	
		species[bathgas["From"][i]] = bathgas["FromConc"][i] 
	speciesTo[bathgas["To"]] = bathgas["ToConc"]
else:
	species[bathgas["From"]] = bathgas["FromConc"] 
	speciesTo[bathgas["To"]] = bathgas["ToConc"]

FromSpecies,ToSpecies = dict_screening(species,speciesTo)
initial_str = dirichlet_str(FromSpecies)
final_str = dirichlet_str(ToSpecies)
final_conc_bound = concBounds_str("temperature",temperatureTo)

#step: 2
#start the simulations
#inputs
#fuel
#oxidizer
#bathgas
#globalrxn
#P
#T
#phi

fuel["fuelIs"] = ""
if fuelType == "Multi-Multi":
	for i in fuel["From"]:
		if fuel["From"][i].strip() != "":
			fuel["fuelIs"] +="fuel is {}\n".format(str(fuel["From"][i]))
		else:
			continue
else:
	fuel["fuelIs"] +="fuel is {}\n".format(str(fuel["From"]))
oxidizer["oxidizerIs"] = oxidizer["From"]
controls = {}
controls["dirichlet"] = initial_str
controls["conc_bounds"] = concBounds_str("species",oxidizer["FromConc"],To = oxidizer["To"],From = oxidizer["From"])
print("Init start_profile is {}".format(inputs["StartProfilesFile"]))
inputs["previousProfile"] = ["hello"]
startProfileLoc = getNewProfile(inputs,phi,fuel,oxidizer,bathgas,globRxn["FromRxn"],start_profile_pressure,init_Temperature,controls)
#update the values after first iteration
inputs["previousProfile"].append(inputs["StartProfilesFile"])
inputs["StartProfilesFile"] = startProfileLoc

print("-1, start_profile is {}".format(inputs["StartProfilesFile"]))

for ind,spec in enumerate(speciesTo):
	inter_species = readjustment(spec,FromSpecies,ToSpecies)
	intermediate_str = dirichlet_str(inter_species)
	bound_str =	concBounds_str("species",speciesTo[spec],To = spec,From = spec)
	controls["dirichlet"] = intermediate_str
	controls["conc_bounds"] = bound_str
	startProfileLoc = getNewProfile(inputs,phi,fuel,oxidizer,bathgas,globRxn["FromRxn"],start_profile_pressure,init_Temperature,controls)
	inputs["previousProfile"].append(inputs["StartProfilesFile"])
	inputs["StartProfilesFile"] = startProfileLoc
	print("{}, start_profile is {}".format(ind,inputs["StartProfilesFile"]))
	FromSpecies = inter_species

#step:3
#All species are converted to their respective targets
#Convert the global rxn as well as the fuel
#switchOff the bounds
controls["dirichlet"] = final_str
controls["conc_bounds"] = concBounds_str("temperature",temperatureTo)
startProfileLoc = getNewProfile(inputs,phi,fuel,oxidizer,bathgas,globRxn["FromRxn"],start_profile_pressure,temperatureTo,controls)
inputs["previousProfile"].append(inputs["StartProfilesFile"])
inputs["StartProfilesFile"] = startProfileLoc

#step: 4
#Convert the globalRxn
if fuelType == "Mono-Multi" or fuelType == "Multi-Multi":
	for i in fuel["To"]:
		if fuel["To"][i].strip() != "":
			fuel["fuelIs"] +="fuel is {}\n".format(str(fuel["To"][i]))
		else:
			continue
else:
	fuel["fuelIs"] +="fuel is {}\n".format(str(fuel["To"]))
controls["conc_bounds"] = ""
startProfileLoc = getNewProfile(inputs,phi,fuel,oxidizer,bathgas,globRxn["ToRxn"],start_profile_pressure,temperatureTo,controls)
inputs["previousProfile"].append(inputs["StartProfilesFile"])
inputs["StartProfilesFile"] = startProfileLoc
#step: 5
#change pressure and achieve the final start profile
startProfileLoc = getNewProfile(inputs,phi,fuel,oxidizer,bathgas,globRxn["ToRxn"],pressureTo,temperatureTo,controls)
inputs["previousProfile"].append(inputs["StartProfilesFile"])
inputs["StartProfilesFile"] = startProfileLoc
#if inputs["exp_type"]["exp"] == "Flf":
#	inputs["exp_type"]["sudo"] = "Flf"
#	startProfileLoc = getNewProfile(inputs,phi,fuel,oxidizer,bathgas,globRxn["ToRxn"],start_profile_pressure,temperatureTo,controls)
#	inputs["startProfileLoc"] = startProfileLoc
print(str(inputs["CopyTo"])+"/"+str(inputs["UniqueID"]))
shutil.copy(startProfileLoc,str(inputs["CopyTo"])+"/"+str(inputs["UniqueID"]))
print("Task succesfully completed")
