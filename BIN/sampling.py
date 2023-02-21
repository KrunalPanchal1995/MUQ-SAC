#python default modules
import numpy as np
import scipy as sp
import scipy.stats as stats
from scipy.optimize import minimize
import os, sys, re, threading, subprocess, time
from sklearn.preprocessing import PolynomialFeatures
from collections import OrderedDict
import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('figure', max_open_warning = 0)
from scipy.linalg import block_diag
from scipy import optimize as spopt
import json
import yaml
import pandas as pd
import csv
#program specific modules
import ParallelWorkers as pk
import Uncertainty 
import combustion_target_class
import simulation_manager
import ResponseSurface as PRS
import statistics
### KEY WORDS #######
optType = "optimization_type"
targets = "targets"
mech = "mechanism"
pre_file = "Initial_pre_file"
count = "Counts"
countTar = "targets_count"
home_dir = os.getcwd()
fuel = "fuel/s"
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
targetLines = open(locations[targets],'r').readlines()
addendum = yaml.safe_load(open(locations[add],'r').read())
target_list = []
c_index = 0
string_target = ""
for target in targetLines[:targets_count]:
	if "#" in target:
		target = target[:target.index('#')]	
	t = combustion_target_class.combustion_target(target,addendum,c_index)
	string_target+=f"{t.dataSet_id}|{t.target}|{t.species_dict}|{t.temperature}|{t.pressure}|{t.phi}|{t.observed}|{t.std_dvtn}\n"
	c_index +=1
	target_list.append(t)

case_dir = range(0,len(target_list))

#print(targetLines)
##!!!  MAKE A LIST OF TARGET CLASS CONTAINING EACH TARGET AS A CA
##########Obtain uncertainty data from user input file!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

UncertDataSet = Uncertainty.uncertaintyData(locations);
unsrt_data,rxnUnsrt_data,focUnsrt_data,tbdUnsrt_data,thermoUnsrt_data,transportUnsrt_data, reaction_index,fallOffCurve_index, thirdBody_index, thermo_index, transport_index = UncertDataSet.extract_uncertainty();
print("\n\nUncertainty data collected \n\n")

BranchingRxns = []
activeParameters = []
unsrtDatabase = {}
thirdBody_dict = {}
total_m_params = []                                                   
rIndex = []
rIndex_dict = {}
for index,rxn in enumerate(reaction_index):
	branch_bool = rxnUnsrt_data[rxn].branching
	if branch_bool.strip() == "True":
		BranchingRxns.append(rxnUnsrt_data[rxn].rxn) 
	activeParameters.append(rxn)
	unsrtDatabase[rxn] = unsrt_data[rxn].getDtList()
	rIndex.append(unsrt_data[rxn].rIndex)
	rIndex_dict[unsrt_data[rxn].rIndex] = rxn
for i in fallOffCurve_index:
	activeParameters.append(i)
	unsrtDatabase[i] = unsrt_data[i].getDtList()

for i in thirdBody_index:
	activeParameters.append(i)
	unsrtDatabase[i] = unsrt_data[i].getDtList()
	temp = []
	for j in unsrt_data[i].branches.split(","):
		temp.append(j)
	thirdBody_dict[i] = temp
	total_m_params.append(temp)
	
for i in thermo_index:
	activeParameters.append(i)
	unsrtDatabase[i] = unsrt_data[i].getDtList()

for i in transport_index:
	activeParameters.append(i)
	unsrtDatabase[i] = unsrt_data[i].getDtList()

print("Uncertainty data acquired \n")
#print(rIndex_dict)
selectedParams = []
for rxn in rxnUnsrt_data:
	#print(rxn)
	selectedParams.extend(rxnUnsrt_data[rxn].activeParameters)

zeta_file = optInputs["Locations"]["zeta_file"]
"""
Get the samples from the simulation manager
"""
sample_length = 50
#pk.Sample(rxnUnsrt_data,reaction_index,sample_length).get_SAMAP(samap_executable,jpdap_executable)
pk.Sample(rxnUnsrt_data,reaction_index,sample_length).plotSamples(zeta_file)

