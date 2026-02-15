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
from copy import deepcopy
#import yaml
#import ruamel.yaml as yaml
try:
    import ruamel_yaml as yaml
except ImportError:
    from ruamel import yaml
import pandas as pd
import yaml
####################################
##  Importing the sampling file   ##
##                                ##
####################################

#sys.path.append('/home/krithika/Desktop/MUQ-SAC/MUQSAC')
#print(sys.path)
#import sensitivity
#import test,Nominal,Optimal
#import generate
#import solution
from MechManipulator2_0 import Manipulator
#program specific modules
from MechanismParser import Parser
#import make_input_file
#import FlameMaster_in_parallel
#import combustion_dataset_class
#import combustion_variable_class
#from combustion_optimization_class import OptimizationTool
from OptimizationTool import OptimizationTool as Optimizer
import combustion_target_class
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
import ResponseSurface as PRS
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
if len(sys.argv) > 1:
	input_file = open(sys.argv[1],'r')
	optInputs = yaml.safe_load(input_file)
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
	print(target)
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

#UncertDataSet = uncertainty.uncertaintyData(locations,binLoc);

############################################
##   Get unsrt data from UncertDataSet    ##
############################################

#unsrt_data = UncertDataSet.extract_uncertainty();
#unsrt_data,rxnUnsrt_data,plogUnsrt_data,plog_interpolated_data,focUnsrt_data,tbdUnsrt_data,thermoUnsrt_data,transportUnsrt_data, reaction_index,plog_boundary_index,plog_index,fallOffCurve_index, thirdBody_index, thermo_index, transport_index = UncertDataSet.extract_uncertainty();
#print("Uncertainty analysis finished")

###########################################################
#### Printing the reactions and their index in the file ###
###########################################################

#string_Rxn = ""
#for i in unsrt_data:
#	string_Rxn += f"{i}\t{unsrt_data[i].index}\n"
#file_rxn = open("RXN.csv","w").write(string_Rxn)


######################################################################
##  CREATING A DICTIONARY CONTAINING ALL THE DATA FROM UNSRT CLASS  ##
## 																    ##
######################################################################
#manipulationDict = {}
#selection = []
#Cholesky_list = []
#zeta_list = []
#activeParameters = []
#nominal_list = []
#rxn_list = []
#rIndex = []
#for rxn in unsrt_data:
	#print(i)
#	rxn_list.append(rxn)
#	selection.extend(unsrt_data[rxn].selection)
#	Cholesky_list.append(unsrt_data[rxn].cholskyDeCorrelateMat)
#	zeta_list.append(unsrt_data[rxn].perturb_factor)
#	activeParameters.extend(unsrt_data[rxn].activeParameters)
#	nominal_list.append(unsrt_data[rxn].nominal)
#	rIndex.append(unsrt_data[rxn].rIndex)
	
#manipulationDict["selection"] = deepcopy(selection)#.deepcopy()
#manipulationDict["Cholesky"] = deepcopy(Cholesky_list)#.deepcopy()
#manipulationDict["zeta"] = deepcopy(zeta_list)#.deepcopy()
#manipulationDict["activeParameters"] = deepcopy(activeParameters)#.deepcopy()
#manipulationDict["nominal"] = deepcopy(nominal_list)#.deepcopy()

#print(manipulationDict["activeParameters"])

##################################################################
##  Use the unsrt data to sample the Arrhenius curves           ##
##  MUQ-SAC: Method of Uncertainty Quantification and           ##
##	Sampling of Arrhenius Curves                                ##
##################################################################

"""
Function Inputs:
	Unsrt Data:
Output:
	Desired Number of Samples - n_s for
		- Class A curves
		- Class B curves
		- Class C curves

"""


#def getTotalUnknowns(N):
#	n_ = 1 + 2*N + (N*(N-1))/2
#	return int(n_)
	

#if "DesignMatrix.csv" not in os.listdir():
#	design_matrix = DM.DesignMatrix(unsrt_data,design_type,7*getTotalUnknowns(len(manipulationDict["activeParameters"]))).getSamples()
#else:
#	design_matrix_file = open("DesignMatrix.csv").readlines()
#	design_matrix = []
#	for row in design_matrix_file:
#		design_matrix.append([float(ele) for ele in row.strip("\n").strip(",").split(",")])
#raise AssertionError("Design Matrix created!!")

##############################################################
##     Doing simulations using the design matrix            ## 
##     The goal is to test the design matrix                ##
##############################################################

if "NOMINAL" not in os.listdir():
	os.mkdir("NOMINAL")
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
		if target.target == "Tig" or target.target == 'RCM' or target.target == "JSR":
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
		g = open(f"../Plot/Dataset/{folder}/"+d_set+".csv","+w").write(string_1)
	else:
		g = open("../Plot/Dataset/Fls/"+d_set+".csv","+w").write(string_2)
print("code run succesfully..return 0")








