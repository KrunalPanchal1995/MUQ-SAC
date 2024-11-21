#python default modules
import numpy as np
import scipy as sp
import scipy.stats as stats
from scipy.optimize import minimize
import os, sys, re, threading, subprocess, time
from sklearn.model_selection import train_test_split
#from sklearn.preprocessing import PolynomialFeatures
from collections import OrderedDict
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('figure', max_open_warning = 0)
from scipy.linalg import block_diag
from scipy import optimize as spopt
import json
import multiprocessing
import concurrent.futures
import asyncio
import pickle
import combustion_target_class
#import partial_PRS_system as P_PRS
#import yaml
#import ruamel.yaml as yaml
try:
    import ruamel_yaml as yaml
except ImportError:
    from ruamel import yaml
import pandas as pd
import yaml
sys.path.append('/parallel_yaml_writer.so')
#import parallel_yaml_writer
#print(dir(parallel_yaml_writer))
#raise AssertionError("Stop!")
####################################
##  Importing the sampling file   ##
##                                ##
####################################
#import sensitivity
#import test,Nominal,Optimal
#import generate
#import solution
from MechManipulator import Manipulator
#program specific modules
from copy import deepcopy
from MechanismParser import Parser
#import make_input_file
#import FlameMaster_in_parallel
#import combustion_dataset_class
#import combustion_variable_class
#from combustion_optimization_class import OptimizationTool
#from OptimizationTool import OptimizationTool as Optimizer
#import combustion_target_class
import data_management
#import data_management as dm
import simulation_manager2_0 as simulator
import Uncertainty as uncertainty
#import MechanismManipulator
#import plotter
#import Input_file_reader
#import statistics
#import ParallelWorkers as pk
from mpire import WorkerPool
import DesignMatrix as DM
#import ResponseSurface as PRS


from append_list import get_rxn_list


### KEY WORDS #######
optType = "optimization_type"
targets = "targets"
mech = "mechanism"
pre_file = "Initial_pre_file"
count = "Counts"
countTar = "targets_count"
home_dir = os.getcwd()
fuel = "fuel"
fuelClass = "fuelClass"
bin_solve = "solver_bin"
bin_opt = "bin"
globRxn = "global_reaction"
countThreads = "parallel_threads"
unsrt = "uncertainty_data"
thermoF = "thermo_file"
transF = "trans_file"
order = "Order_of_PRS"
startProfile = "StartProfilesData"
design = "Design_of_PRS"
countRxn = "total_reactions"
fT = "fileType"
add = "addendum"
###########
#open the input file and check for arguements
###########
global optInputs
if len(sys.argv) > 1:
	input_file = open(sys.argv[1],'r')
	optInputs = yaml.safe_load(input_file)
	import yaml
	print("Input file found\n")
else:
	print("Please enter a valid input file name as arguement. \n For details of preparing the input file, please see the UserManual\n\nProgram exiting")
	exit()

#!!!!!!! GET MECHANISM FILE , number of targets  from the input file !!!!!!!!!
iFile = str(os.getcwd())+"/"+str(sys.argv[1])
dataCounts = optInputs[count]
binLoc = optInputs["Bin"]
inputs = optInputs["Inputs"]
locations = optInputs["Locations"]
startProfile_location = optInputs[startProfile]
stats_ = optInputs["Stats"]
IMPORTANT_RXN_LIST = optInputs["Locations"]["SA_rxn_list"]
unsrt_location = locations[unsrt]
mech_file_location = locations[mech]
thermo_file_location = locations[thermoF]
trans_file_location = locations[transF]
fileType = inputs[fT]
samap_executable = optInputs["Bin"]["samap_executable"]
jpdap_executable = optInputs["Bin"]["jpdap_executable"]

if fileType == "chemkin":
	file_specific_input = "-f chemkin"
else:
	file_specific_input = ""
fuel = inputs[fuel]
global_reaction = inputs[globRxn]
design_type = stats_[design]
parallel_threads = dataCounts[countThreads]
targets_count = int(dataCounts["targets_count"])
rps_order = stats_[order]
PRS_type = stats_["PRS_type"]
#######################READ TARGET FILE ###################
print("Parallel threads are {}".format(parallel_threads))
targetLines = open(locations[targets],'r').readlines()
addendum = yaml.safe_load(open(locations[add],'r').read())
print(design_type)


####################################################
##  Unloading the target data	  		          ##
## TARGET CLASS CONTAINING EACH TARGET AS A	CASE  ##
####################################################


target_list = []
c_index = 0
string_target = ""

for target in targetLines[:targets_count]:
	if "#" in target:
		target = target[:target.index('#')]	
	add = deepcopy(addendum)
	t = combustion_target_class.combustion_target(target,add,c_index)
	string_target+=f"{t.dataSet_id}|{t.target}|{t.species_dict}|{t.temperature}|{t.pressure}|{t.phi}|{t.observed}|{t.std_dvtn}\n"
	c_index +=1
	target_list.append(t)
case_dir = range(0,len(target_list))

print("\n\n##########################################")
print("##\tNominal targets identified\t##")
print("##########################################\n")
target_file = open("target_data.txt","w")
target_file.write(string_target)
target_file.close()

############################################
##  Uncertainty Quantification            ##
##  					   ##
############################################

if "unsrt.pkl" not in os.listdir():
	UncertDataSet = uncertainty.uncertaintyData(locations,binLoc);
	############################################
	##   Get unsrt data from UncertDataSet    ##
	############################################

	unsrt_data = UncertDataSet.extract_uncertainty();
	# Save the object to a file
	with open('unsrt.pkl', 'wb') as file_:
		pickle.dump(unsrt_data, file_)
	print("Uncertainty analysis finished")

else:
	# Load the object from the file
	with open('unsrt.pkl', 'rb') as file_:
		unsrt_data = pickle.load(file_)
	print("Uncertainty analysis already finished")


######################################################################
##  CREATING A DICTIONARY CONTAINING ALL THE DATA FROM UNSRT CLASS  ##
##																    ##
######################################################################
manipulationDict = {}
selection = []
Cholesky_list = []
zeta_list = []
activeParameters = []
nominal_list = []
rxn_list = []
rIndex = []
rxn_list_a_fact = []

for rxn in unsrt_data:
	#print(i)
	data = {}
	data["temperatures"] = unsrt_data[rxn].temperatures
	data["uncertainties"] = unsrt_data[rxn].uncertainties
	data["Arrhenius"] = unsrt_data[rxn].nominal
	#print(rxn,data)
	rxn_list.append(rxn)
	rxn_list_a_fact.append(unsrt_data[rxn].rxn)
	selection.extend(unsrt_data[rxn].selection)
	Cholesky_list.append(unsrt_data[rxn].cholskyDeCorrelateMat)
	zeta_list.append(unsrt_data[rxn].perturb_factor)
	activeParameters.extend(unsrt_data[rxn].activeParameters)
	nominal_list.append(unsrt_data[rxn].nominal)
	rIndex.append(unsrt_data[rxn].rIndex)

manipulationDict['activeParameters'] = activeParameters
rxn_list_a_fact = set(rxn_list_a_fact)
rxn_a_fact = ""
for rxn in rxn_list_a_fact:
	rxn_a_fact+=f"{rxn}\n"
file_rxn_a_fact = open("rxn_list_a_fact.csv","w").write(rxn_a_fact)	

IMPORTANT_RXN_LIST = "rxn_list_a_fact.csv"

######################################
##  SENSITIVITY ANALYSIS            ##
##  				    ##
######################################

zeta_list = []
for rxn in unsrt_data:
	zeta_list.append(unsrt_data[rxn].perturb_factor)

"""
# Load the object from the file sensitivity_rx.pkl
if "sens_parameters.pkl" not in os.listdir():
		# Arguments to pass to script2.py
		args = [sys.argv[1], IMPORTANT_RXN_LIST,'&>SA.out']

		# Doing A factor sensitivity Analysis
		result = subprocess.run(['python3.9', binLoc["SA_tool"]] + args, capture_output=True, text=True)
		f = open("SA_A_FACTOR_tool.out","+a")
		f.write(result.stdout+"\n")
		print("\nSensitivity Analysis for A-factor is Done!!")
		# Printing the errors of script2.py, if any
		if result.stderr:
			f.write("Errors:\n"+result.stderr)	
		#	raise AssertionError("Sensitivity Analysis Done!!")
		with open('sens_parameters.pkl', 'rb') as file_:
			sensitivity_analysis = pickle.load(file_)
else:
		print("\n\tBrute-force Sensitivity analysis is alsready over!!\n")
		with open('sens_parameters.pkl', 'rb') as file_:
			sensitivity_analysis = pickle.load(file_)

with open('RXN_DICT.pkl', 'rb') as file_:
	RXN_DICT = pickle.load(file_)
#print(RXN_DICT)

"""
"""
Function Inputs:
	Unsrt Data:
Output:
	Desired Number of Samples - n_s for
		- Class A curves
		- Class B curves
		- Class C curves

"""
	
def getTotalUnknowns(N):
	n_ = 1 + 2*N + (N*(N-1))/2
	return int(n_)
	
def getSim(n,design):
	if "sampling_size" in optInputs["Stats"]:
		return int(optInputs["Stats"]["sampling_size"])
	
	else:
		n_ = getTotalUnknowns(n)
		if design == "A-facto":
			sim = 4*n_
		else:
			sim = 7*n_	
		return sim
no_of_sim = {}
if "DesignMatrix.csv" not in os.listdir():
	print(f"\nNo. of Simulations required: {getSim(len(manipulationDict['activeParameters']),design)}\n")
	no_of_sim_ = getSim(len(manipulationDict["activeParameters"]),design_type)
	design_matrix = DM.DesignMatrix(unsrt_data,design_type,getSim(len(manipulationDict["activeParameters"]),design_type),len(manipulationDict["activeParameters"])).getSamples()
	no_of_sim_ = len(design_matrix)
else:
	no_of_sim_ = getSim(len(manipulationDict["activeParameters"]),design_type)
	design_matrix_file = open("DesignMatrix.csv").readlines()
	design_matrix = []
	no_of_sim_ = len(design_matrix_file)
	for row in design_matrix_file:
		design_matrix.append([float(ele) for ele in row.strip("\n").strip(",").split(",")])
#raise AssertionError("Design Matrix created!!")
design_matrix_dict = {}
for case in case_dir:
	design_matrix_dict[case] = design_matrix
	no_of_sim[case] = no_of_sim_

SSM = simulator.SM(target_list,optInputs,unsrt_data,design_matrix)

if "Perturbed_Mech" not in os.listdir():
	os.mkdir("Perturbed_Mech")
	print("\nPerturbing the Mechanism files\n")
	#if "YAML_PERTUBED_FILES.pkl" not in os.listdir():
		
	#	yaml_list = SSM.getYAML_List()
	
		#with open('YAML_PERTUBED_FILES.pkl', 'wb') as file_:
		#	pickle.dump(yaml_list, file_)
	#else:
		#yaml_list = SSM.getYAML_List()
		#with open('YAML_PERTUBED_FILES.pkl', 'rb') as file_:
		
		#	yaml_list = pickle.load(file_)

	chunk_size = 500
	params_yaml = [design_matrix[i:i+chunk_size] for i in range(0, len(design_matrix), chunk_size)]
	count = 0
	yaml_loc = []
	for params in params_yaml:
		
		yaml_list = SSM.getYAML_List(params)
		#yaml_loc = []
		location_mech = []
		index_list = []
		for i,dict_ in enumerate(yaml_list):
			index_list.append(str(count+i))
			location_mech.append(os.getcwd()+"/Perturbed_Mech")
			yaml_loc.append(os.getcwd()+"/Perturbed_Mech/mechanism_"+str(count+i)+".yaml")
		count+=len(yaml_list)
		#gen_flag = False
		#SSM.getPerturbedMechLocation(yaml_list,location_mech,index_list)
		SSM.getPerturbedMechLocation(yaml_list,location_mech,index_list)
		print(f"\nGenerated {count} files!!\n")
	print("\nGenerated the YAML files required for simulations!!\n")
else:
	print("\nYAML files already generated!!")
	yaml_loc = []
	location_mech = []
	index_list = []
	for i,sample in enumerate(design_matrix):
		index_list.append(i)
		location_mech.append(os.getcwd()+"/Perturbed_Mech")
		yaml_loc.append(os.getcwd()+"/Perturbed_Mech/mechanism_"+str(i)+".yaml")

selected_params = []
for params in activeParameters:
	selected_params.append(1)

selected_params_dict = {}
design_matrix_dict = {}
yaml_loc_dict = {}
for case in case_dir:
	yaml_loc_dict[case] = yaml_loc
	design_matrix_dict[case] = design_matrix
	selected_params_dict[case] = selected_params
##############################################################
##     Doing simulations using the design matrix            ## 
##     The goal is to test the design matrix                ##
##############################################################

if "SIMULATION_FOLDER" not in os.listdir():
	os.mkdir("SIMULATION_FOLDER")
	os.chdir("SIMULATION_FOLDER")
	optDir = os.getcwd()
	os.mkdir("Data")
	os.chdir("Data")
	os.mkdir("Simulations")
	os.mkdir("ResponseSurface")
	os.chdir("..")
else:
	os.chdir("SIMULATION_FOLDER")
	optDir = os.getcwd()

if os.path.isfile("progress") == False:
	FlameMaster_Execution_location = SSM.make_dir_in_parallel(yaml_loc_dict)
	#raise AssertionError("The Target class, Uncertainty class, Design Matrix and Simulations")
else:
	print("Progress file detected")
	progress = open(optDir+"/progress",'r').readlines()
	FlameMaster_Execution_location = []
	with open(optDir+"/locations") as infile:
		for line in infile:
			FlameMaster_Execution_location.append(line)
	#FlameMaster_Execution_location = open(optDir+"/locations",'r').readlines()
	#missing_location = data_management.find_missing_location(FlameMaster_Execution_location, progress)

#############################################
##     Extracting simulation results       ##
##       (The module if validated          ##
#############################################
temp_sim_opt = {}
for case in case_dir:	
	os.chdir("Data/Simulations")
	if "sim_data_case-"+str(case)+".lst" in os.listdir():
		ETA = [float(i.split("\t")[1]) for i in open("sim_data_case-"+str(case)+".lst").readlines()]
		temp_sim_opt[str(case)] = ETA
		os.chdir(optDir)
		#print(ETA)
		#raise AssertionError("Generating ETA list for all cases")	
	else:
		os.chdir(optDir)
		#print(os.getcwd())
		os.chdir("case-"+str(case))	
		data_sheet,failed_sim, ETA = data_management.generate_target_value_tables(FlameMaster_Execution_location, target_list, case, fuel)
		#print(data_sheet)
		#raise AssertionError("!STOP")
		temp_sim_opt[str(case)] = ETA
		f = open('../Data/Simulations/sim_data_case-'+str(case)+'.lst','w').write(data_sheet)
		g = open('../Data/Simulations/failed_sim_data_case-'+str(case)+'.lst','w').write(failed_sim)
		#f.write(data_sheet)
		#f.close()
		os.chdir(optDir)


raise AssertionError(f"Simulations Finished: {int(optInputs['Stats']['sampling_size'])} simulations for {len(target_list)} targets")

###############################################
##     Finding the      ##
##                                           ##
###############################################

ResponseSurfaces = {}
selected_PRS = {}
for case_index,case in enumerate(temp_sim_opt):
	yData = np.asarray(temp_sim_opt[case]).flatten()[0:no_of_sim[case_index]]
	xData = np.asarray(design_matrix_dict[case_index])[0:no_of_sim[case_index]]
	#print(np.shape(xData))
	#print(np.shape(yData))
	#raise AssertionError("Stop!!")



#print(sens_data)
new_rxn_dict, new_sc_dict = get_rxn_list(sensitivity_analysis , unsrt_data) 

if "UPPER_BOUND" not in os.listdir() and "LOWER_BOUND" not in os.listdir():
	os.mkdir("UPPER_BOUND")
	os.mkdir("LOWER_BOUND")

originalMech = Parser(mech_file_location).mech
copy_of_mech = deepcopy(originalMech)
zeta_list = np.asarray(zeta_list)

YAML_UPPER_DICT = {}
YAML_LOWER_DICT = {}
for case in new_sc_dict:
	new_sc_list = np.asarray(new_sc_dict[case])
	#print(new_sc_list,zeta_list)
	positive_zeta_list = zeta_list *  np.sign(new_sc_list)
	negative_zeta_list = -1*zeta_list* np.sign(new_sc_list)
	####### Positive sign zeta list
	new_mechanism_upper,a = Manipulator(copy_of_mech,unsrt_data,positive_zeta_list.flatten(),rxn_dict = RXN_DICT).doPerturbation()
	string = yaml.dump(new_mechanism_upper,default_flow_style=False)
	f = open(f"UPPER_BOUND/new_mech_{case}.yaml","w").write(string)
	YAML_UPPER_DICT[case] = Parser(f"UPPER_BOUND/new_mech_{case}.yaml").mech
	####### Negative sign zeta list
	new_mechanism_lower,a = Manipulator(copy_of_mech,unsrt_data,negative_zeta_list.flatten(),rxn_dict = RXN_DICT).doPerturbation()
	string = yaml.dump(new_mechanism_lower,default_flow_style=False)
	f = open(f"LOWER_BOUND/new_mech_{case}.yaml","w").write(string)
	YAML_LOWER_DICT[case] = Parser(f"LOWER_BOUND/new_mech_{case}.yaml").mech
###################################
#### Nominal Simulations ##########
###################################

##############################################################
##     Doing simulations using the design matrix            ## 
##     The goal is to test the design matrix                ##
##############################################################

if "NOMINAL" not in os.listdir():
	os.mkdir("NOMINAL")
	os.mkdir("Plot_NOMINAL")
	os.mkdir("Plot_NOMINAL/Dataset")
	os.mkdir("Plot_NOMINAL/Dataset/Tig")
	os.mkdir("Plot_NOMINAL/Dataset/Fls")
	os.mkdir("Plot_NOMINAL/Dataset/JSR")
	os.mkdir("Plot_NOMINAL/Dataset/RCM")
	os.chdir("NOMINAL")
	nominalDir = os.getcwd()
	os.mkdir("Data")
	os.chdir("Data")
	os.mkdir("Simulations")
	os.chdir("..")
else:
	os.chdir("NOMINAL")
	nominalDir = os.getcwd()

unsrt_data = {}
design_matrix = []
if os.path.isfile("progress") == False:
	FlameMaster_Execution_location = simulator.SM(target_list,optInputs,unsrt_data,design_matrix).make_nominal_dir_in_parallel()
	#print(FlameMaster_Execution_location)
	#raise AssertionError("The Target class, Uncertainty class, Design Matrix and Simulations")
else:
	print("Progress file detected")
	progress = open(nominalDir+"/progress",'r').readlines()
	FlameMaster_Execution_location = []
	with open(nominalDir+"/locations") as infile:
		for line in infile:
			FlameMaster_Execution_location.append(line)

#############################################
##     Extracting simulation results       ##
##       (The module if validated          ##
#############################################
temp_sim_opt = {}
#dataset = "#DS_ID,T,Obs(us),Nominal"
for case in case_dir:	
	os.chdir("Data/Simulations")
	if "sim_data_case-"+str(case)+".lst" in os.listdir():
		ETA = [np.exp(float(i.split("\t")[1]))/10 for i in open("sim_data_case-"+str(case)+".lst").readlines()]
		temp_sim_opt[str(case)] = ETA
		os.chdir(nominalDir)
	
	else:
		os.chdir(nominalDir)
		#print(os.getcwd())
		os.chdir("case-"+str(case))	
		data_sheet,failed_sim, ETA = data_management.generate_target_value_tables(FlameMaster_Execution_location, target_list, case, fuel)
		temp_sim_opt[str(case)] = np.exp(ETA)/10
		f = open('../Data/Simulations/sim_data_case-'+str(case)+'.lst','w').write(data_sheet)
		g = open('../Data/Simulations/failed_sim_data_case-'+str(case)+'.lst','w').write(failed_sim)
		os.chdir(nominalDir)

#for d_set in dataset:
def listToString(s):
 
    # initialize an empty string
    str1 = ""
 
    # traverse in the string
    for ele in s:
        str1 += f"{ele},"
 
    # return string
    return str1

dataset = []
for t in target_list:
	dataset.append(t.dataSet_id)
dataset = set(dataset)
for d_set in dataset:
	print(d_set)
	string_1 = "DS_ID,T,Obs(us),Nominal\n"
	string_2 = "DS_ID,T,P,Phi,Fuel,Ox,BathGas,Obs(us),Nominal\n"
	flag = None
	for case,target in enumerate(target_list):
		if target.target == "Tig" or target.target == "JSR":
			#print(d_set)
			#print(target.target)
			if target.dataSet_id == d_set:
				folder = target.target
				#print(d_set)
				#print(temp_sim_opt[str(target.index)])
				flag =True
				string_1 += f"{target.uniqueID},{target.temperature},{target.observed},{temp_sim_opt[str(target.index)][0]}\n"
		elif target.target == "Fls":
			if target.dataSet_id == d_set:
				#print(target.fuel_dict,target.BG_dict)
				#raise AssertionError("Stop!")
				flag =False	
				string_2 += f"{target.uniqueID},{target.temperature},{target.pressure},{target.phi},{target.fuel_dict},{target.oxidizer_x},{target.BG_dict},{target.observed},{temp_sim_opt[str(target.index)][0]}\n"
		#else
	if flag == True:
		folder_path = f"../Plot/Dataset/{folder}/"
		os.makedirs(folder_path, exist_ok=True) 
		g = open(f"../Plot_NOMINAL/Dataset/{folder}/"+d_set+".csv","+w").write(string_1)
	else:
		g = open("../Plot_NOMINAL/Dataset/Fls/"+d_set+".csv","+w").write(string_2)
#raise AssertionError("Stop!!")

os.chdir("..")






##############################################################
##     Doing simulations using the design matrix            ## 
##     The goal is to test the design matrix                ##
##############################################################

if "NOMINAL_UPPER" not in os.listdir():
	os.mkdir("NOMINAL_UPPER")
	os.mkdir("Plot_UPPER")
	os.mkdir("Plot_UPPER/Dataset")
	os.mkdir("Plot_UPPER/Dataset/Tig")
	os.mkdir("Plot_UPPER/Dataset/Fls")
	os.mkdir("Plot_UPPER/Dataset/JSR")
	os.mkdir("Plot_UPPER/Dataset/RCM")
	os.chdir("NOMINAL_UPPER")
	nominalDir = os.getcwd()
	os.mkdir("Data")
	os.chdir("Data")
	os.mkdir("Simulations")
	os.chdir("..")
else:
	os.chdir("NOMINAL_UPPER")
	nominalDir = os.getcwd()

unsrt_data = {}
design_matrix = []
if os.path.isfile("progress") == False:
	FlameMaster_Execution_location = simulator.SM(target_list,optInputs,unsrt_data,design_matrix).make_nominal_dir_in_parallel(yaml_list = YAML_UPPER_DICT)
	#print(FlameMaster_Execution_location)
	#raise AssertionError("The Target class, Uncertainty class, Design Matrix and Simulations")
else:
	print("Progress file detected")
	progress = open(nominalDir+"/progress",'r').readlines()
	FlameMaster_Execution_location = []
	with open(nominalDir+"/locations") as infile:
		for line in infile:
			FlameMaster_Execution_location.append(line)

#############################################
##     Extracting simulation results       ##
##       (The module if validated          ##
#############################################
temp_sim_opt = {}
#dataset = "#DS_ID,T,Obs(us),Nominal"
for case in case_dir:	
	os.chdir("Data/Simulations")
	if "sim_data_case-"+str(case)+".lst" in os.listdir():
		ETA = [np.exp(float(i.split("\t")[1]))/10 for i in open("sim_data_case-"+str(case)+".lst").readlines()]
		temp_sim_opt[str(case)] = ETA
		os.chdir(nominalDir)
	
	else:
		os.chdir(nominalDir)
		#print(os.getcwd())
		os.chdir("case-"+str(case))	
		data_sheet,failed_sim, ETA = data_management.generate_target_value_tables(FlameMaster_Execution_location, target_list, case, fuel)
		temp_sim_opt[str(case)] =  np.exp(ETA)/10
		f = open('../Data/Simulations/sim_data_case-'+str(case)+'.lst','w').write(data_sheet)
		g = open('../Data/Simulations/failed_sim_data_case-'+str(case)+'.lst','w').write(failed_sim)
		os.chdir(nominalDir)
#for d_set in dataset:
def listToString(s):
 
    # initialize an empty string
    str1 = ""
 
    # traverse in the string
    for ele in s:
        str1 += f"{ele},"
 
    # return string
    return str1

dataset = []
for t in target_list:
	dataset.append(t.dataSet_id)
dataset = set(dataset)

for d_set in dataset:
	print(d_set)
	string_1 = "DS_ID,T,Obs(us),Nominal\n"
	string_2 = "DS_ID,T,P,Phi,Fuel,Ox,BathGas,Obs(us),Nominal\n"
	flag = None
	for case,target in enumerate(target_list):
		if target.target == "Tig" or target.target == "JSR":
			#print(d_set)
			#print(target.target)
			if target.dataSet_id == d_set:
				folder = target.target
				#print(d_set)
				#print(temp_sim_opt[str(target.index)])
				flag =True
				string_1 += f"{target.uniqueID},{target.temperature},{target.observed},{temp_sim_opt[str(target.index)][0]}\n"
		elif target.target == "Fls":
			if target.dataSet_id == d_set:
				#print(target.fuel_dict,target.BG_dict)
				#raise AssertionError("Stop!")
				flag =False	
				string_2 += f"{target.uniqueID},{target.temperature},{target.pressure},{target.phi},{target.fuel_dict},{target.oxidizer_x},{target.BG_dict},{target.observed},{temp_sim_opt[str(target.index)][0]}\n"
		#else
	if flag == True:
		g = open(f"../Plot_UPPER/Dataset/{folder}/"+d_set+".csv","+w").write(string_1)
	else:
		g = open("../Plot_UPPER/Dataset/Fls/"+d_set+".csv","+w").write(string_2)
#raise AssertionError("Stop!!")


os.chdir("..")
##############################################################
##     Doing simulations using the design matrix            ## 
##     The goal is to test the design matrix                ##
##############################################################

if "NOMINAL_LOWER" not in os.listdir():
	os.mkdir("NOMINAL_LOWER")
	os.mkdir("Plot_LOWER")
	os.mkdir("Plot_LOWER/Dataset")
	os.mkdir("Plot_LOWER/Dataset/Tig")
	os.mkdir("Plot_LOWER/Dataset/Fls")
	os.mkdir("Plot_LOWER/Dataset/JSR")
	os.mkdir("Plot_LOWER/Dataset/RCM")
	os.chdir("NOMINAL_LOWER")
	nominalDir = os.getcwd()
	os.mkdir("Data")
	os.chdir("Data")
	os.mkdir("Simulations")
	os.chdir("..")
else:
	os.chdir("NOMINAL_LOWER")
	nominalDir = os.getcwd()

unsrt_data = {}
design_matrix = []
if os.path.isfile("progress") == False:
	FlameMaster_Execution_location = simulator.SM(target_list,optInputs,unsrt_data,design_matrix).make_nominal_dir_in_parallel(yaml_list = YAML_LOWER_DICT)
	#print(FlameMaster_Execution_location)
	#raise AssertionError("The Target class, Uncertainty class, Design Matrix and Simulations")
else:
	print("Progress file detected")
	progress = open(nominalDir+"/progress",'r').readlines()
	FlameMaster_Execution_location = []
	with open(nominalDir+"/locations") as infile:
		for line in infile:
			FlameMaster_Execution_location.append(line)

#############################################
##     Extracting simulation results       ##
##       (The module if validated          ##
#############################################
temp_sim_opt = {}
#dataset = "#DS_ID,T,Obs(us),Nominal"
for case in case_dir:	
	os.chdir("Data/Simulations")
	if "sim_data_case-"+str(case)+".lst" in os.listdir():
		ETA = [np.exp(float(i.split("\t")[1]))/10 for i in open("sim_data_case-"+str(case)+".lst").readlines()]
		temp_sim_opt[str(case)] = ETA
		os.chdir(nominalDir)
	
	else:
		os.chdir(nominalDir)
		#print(os.getcwd())
		os.chdir("case-"+str(case))	
		data_sheet,failed_sim, ETA = data_management.generate_target_value_tables(FlameMaster_Execution_location, target_list, case, fuel)
		temp_sim_opt[str(case)] =  np.exp(ETA)/10
		f = open('../Data/Simulations/sim_data_case-'+str(case)+'.lst','w').write(data_sheet)
		g = open('../Data/Simulations/failed_sim_data_case-'+str(case)+'.lst','w').write(failed_sim)
		os.chdir(nominalDir)
#for d_set in dataset:
def listToString(s):
 
    # initialize an empty string
    str1 = ""
 
    # traverse in the string
    for ele in s:
        str1 += f"{ele},"
 
    # return string
    return str1

dataset = []
for t in target_list:
	dataset.append(t.dataSet_id)
dataset = set(dataset)

for d_set in dataset:
	print(d_set)
	string_1 = "DS_ID,T,Obs(us),Nominal\n"
	string_2 = "DS_ID,T,P,Phi,Fuel,Ox,BathGas,Obs(us),Nominal\n"
	flag = None
	for case,target in enumerate(target_list):
		if target.target == "Tig" or target.target == "JSR":
			#print(d_set)
			#print(target.target)
			if target.dataSet_id == d_set:
				folder = target.target
				#print(d_set)
				#print(temp_sim_opt[str(target.index)])
				flag =True
				string_1 += f"{target.uniqueID},{target.temperature},{target.observed},{temp_sim_opt[str(target.index)][0]}\n"
		elif target.target == "Fls":
			if target.dataSet_id == d_set:
				#print(target.fuel_dict,target.BG_dict)
				#raise AssertionError("Stop!")
				flag =False	
				string_2 += f"{target.uniqueID},{target.temperature},{target.pressure},{target.phi},{target.fuel_dict},{target.oxidizer_x},{target.BG_dict},{target.observed},{temp_sim_opt[str(target.index)][0]}\n"
		#else
	if flag == True:
		folder_path = f"../Plot_LOWER/Dataset/{folder}/"
		os.makedirs(folder_path, exist_ok=True) 
		g = open(f"../Plot_LOWER/Dataset/{folder}/"+d_set+".csv","+w").write(string_1)
	else:
		g = open("../Plot_LOWER/Dataset/Fls/"+d_set+".csv","+w").write(string_2)
#raise AssertionError("Stop!!")


os.chdir("..")
print("code run succesfully..return 0")





















raise AssertionError("stop")	

