#python default modules
import json
import pickle
import asyncio
import numpy as np
import scipy as sp
import pandas as pd
import multiprocessing
import matplotlib as mpl
import concurrent.futures
import parallel_yaml_writer
import scipy.stats as stats
import matplotlib.pyplot as plt
try:
    import ruamel_yaml as yaml
except ImportError:
    from ruamel import yaml
import yaml 

import partial_PRS_system as P_PRS
from scipy.optimize import minimize
from collections import OrderedDict
from scipy.linalg import block_diag
from scipy import optimize as spopt
mpl.rc('figure', max_open_warning = 0)
import os, sys, re, threading, subprocess, time
from sklearn.model_selection import train_test_split
sys.path.append('/parallel_yaml_writer.so')
from scipy.interpolate import Akima1DInterpolator
####################################
##  Importing the sampling file   ##
##                                ##
####################################
#program specific modules
import data_management
from copy import deepcopy
import DesignMatrix as DM
from mpire import WorkerPool
import ResponseSurface as PRS
import combustion_target_class
import Uncertainty as uncertainty
from MechanismParser import Parser
import simulation_manager2_0 as simulator
from MechManipulator3_0_A_factor import Manipulator
from OptimizationTool import OptimizationTool as Optimizer
sys.path.append('/yamlwriter.so')  # Adjust this path to the correct build directory
#print(sys.path)
import yamlwriter

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

################################################
# Open the input file and check for arguements #
################################################

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
global A_fact_samples

A_fact_samples = stats_["Sampling_of_PRS"]
if "sensitive_parameters" not in stats_:
	stats_["sensitive_parameters"] = "Principle_SubMatrix"
	optInputs["Stats"]["sensitive_parameters"] = "Principle_SubMatrix"
if "Arrhenius_Selection_Type" not in stats_:
	stats_["Arrhenius_Selection_Type"] = "some"
	optInputs["Stats"]["Arrhenius_Selection_Type"] = "some"
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
gr = inputs[globRxn]
global_reaction = gr

design_type = stats_[design]
parallel_threads = dataCounts[countThreads]
targets_count = int(dataCounts["targets_count"])
rps_order = stats_[order]
PRS_type = stats_["PRS_type"]
#######################READ TARGET FILE ###################

print("\nParallel threads are {}".format(parallel_threads))
targetLines = open(locations[targets],'r').readlines()
addendum = yaml.safe_load(open(locations[add],'r').read())

####################################################
##  Unloading the target data	  	           ##
## TARGET CLASS CONTAINING EACH TARGET AS A CASE  ##
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
print(case_dir)
print("\n\nOptimization targets identified.\nStarted the MUQ process.......\n")
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
#for rxn in unsrt_data:
#	print(unsrt_data[rxn].linked_rIndex)
#raise AssertionError("Pull over")
###########################################################
#### Printing the reactions and their index in the file ###
###########################################################

string_Rxn = ""
str_rxn = ""
for i in unsrt_data:
	string_Rxn += f"{i}\t{unsrt_data[i].index}\n"
	str_rxn += f"{i}\n"
file_rxn = open("RXN.csv","w").write(string_Rxn)
file_rxn_list = open("rxn_list.csv","w").write(str_rxn)


######################################################################
##  CREATING A DICTIONARY CONTAINING ALL THE DATA FROM UNSRT CLASS  ##
##		                                                    ##
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

rxn_list_a_fact = set(rxn_list_a_fact)
rxn_a_fact = ""
for rxn in rxn_list_a_fact:
	rxn_a_fact+=f"{rxn}\n"
file_rxn_a_fact = open("rxn_list_a_fact.csv","w").write(rxn_a_fact)	
manipulationDict["selection"] = deepcopy(selection)#.deepcopy()
manipulationDict["Cholesky"] = deepcopy(Cholesky_list)#.deepcopy()
manipulationDict["zeta"] = deepcopy(zeta_list)#.deepcopy()
manipulationDict["activeParameters"] = deepcopy(activeParameters)#.deepcopy()
manipulationDict["nominal"] = deepcopy(nominal_list)#.deepcopy()
print(f"\n\n################################################\nThe total reactions choosen for this study: {len(manipulationDict['activeParameters'])}\n\tThe list is as follows:\n\t")
print("\t"+f"{manipulationDict['activeParameters']}")
print("\n################################################\n\n")


"""
For Partial PRS (p-PRS), we need to do a second level of sensitivity analysis:	
	
	sensitivity_analysis: a dictionary that stores the sensitivity coefficients of all the active
			      parameters (selected in the first set of sensitivity analysis) reaction wise
			      as in unsrt_data dictionary
"""


##################################################
##        Optimization Procedure                ##
##   Inputs: Traget list and Response Surfaces  ## 
##################################################
if "solution_zeta.save" not in os.listdir():

	#opt, opt_zeta,posterior_cov = Optimizer(target_list).run_optimization_with_selected_PRS(unsrt_data,ResponseSurfaces,optInputs)
	opt, opt_zeta,posterior_cov = Optimizer(target_list).run_optimization_with_MLE_no_PRS(unsrt_data,optInputs)

	string_save = ""
	for i, j in enumerate(activeParameters):
		print("\n{}=\t{}".format(activeParameters[i], opt_zeta[i]))
		string_save+="{}=\t{}\n".format(activeParameters[i], opt_zeta[i])
	save = open("solution_zeta.save","w").write(string_save)


	#string_save = ""
	#for i, j in enumerate(activeParameters):
	#	print("\n{}=\t{}".format(activeParameters[i], opt[i]))
	#	string_save+="{}=\t{}\n".format(activeParameters[i], opt[i])
	#save = open("solution.save","w").write(string_save)


	originalMech = Parser(mech_file_location).mech
	copy_of_mech = deepcopy(originalMech)#.deepcopy()
	new_mechanism,a = Manipulator(copy_of_mech,unsrt_data,opt_zeta).doPerturbation()

	#new_mechanism,a,b,c = self.MechManipulator.GeneratePerturbedMechanism(self.target_list[case],self.beta_[i],np.ones(len(selectedParams)),reactionList,self.simulation,"False",extra_arg = self.activeIndexDict)
	yamlwriter.dump_to_yaml(os.getcwd(),f"Opt_mechanism.yaml",new_mechanism)
	#string = yaml.dump(new_mechanism,default_flow_style=False)
	#f = open("new_mech.yaml","w").write(string)
else:
	
	save = open("solution_zeta.save","r").readlines()
	
	opt_zeta = []
	for i in save:
		opt_zeta.append(float(i.split("=")[1].strip()))
	originalMech = Parser(mech_file_location).mech
	copy_of_mech = deepcopy(originalMech)#.deepcopy()
	new_mechanism,a = Manipulator(copy_of_mech,unsrt_data,opt_zeta).doPerturbation()

	#new_mechanism,a,b,c = self.MechManipulator.GeneratePerturbedMechanism(self.target_list[case],self.beta_[i],np.ones(len(selectedParams)),reactionList,self.simulation,"False",extra_arg = self.activeIndexDict)
	
	yamlwriter.dump_to_yaml(os.getcwd(),f"Opt_mechanism.yaml",new_mechanism)
	#string = yaml.dump(new_mechanism,default_flow_style=False)
	#f = open("new_mech.yaml","w").write(string)



raise AssertionError("The Target class, Uncertainty class, Design Matrix and Simulations and Response surface")
