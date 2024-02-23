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
import subprocess
import time
import sys
import concurrent.futures
import asyncio

#import yaml
#import ruamel.yaml as yaml
try:
    import ruamel_yaml as yaml
except ImportError:
    from ruamel import yaml
import pandas as pd

####################################
##  Importing the sampling file   ##
##                                ##
####################################
#import sensitivity
#import test,Nominal,Optimal
#import generate
#import solution
from MechManipulator2_0 import Manipulator
#program specific modules
from copy import deepcopy
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
print("optimization targets identified\n")
target_file = open("target_data.txt","w")
target_file.write(string_target)
target_file.close()


############################################
##  Uncertainty Quantification  		  ##
##  									  ##
############################################

UncertDataSet = uncertainty.uncertaintyData(locations,binLoc);

############################################
##   Get unsrt data from UncertDataSet    ##
############################################

unsrt_data = UncertDataSet.extract_uncertainty();
#unsrt_data,rxnUnsrt_data,plogUnsrt_data,plog_interpolated_data,focUnsrt_data,tbdUnsrt_data,thermoUnsrt_data,transportUnsrt_data, reaction_index,plog_boundary_index,plog_index,fallOffCurve_index, thirdBody_index, thermo_index, transport_index = UncertDataSet.extract_uncertainty();
print("Uncertainty analysis finished")

###########################################################
#### Printing the reactions and their index in the file ###
###########################################################

string_Rxn = ""
for i in unsrt_data:
	string_Rxn += f"{i}\t{unsrt_data[i].index}\n"
file_rxn = open("RXN.csv","w").write(string_Rxn)


######################################################################
##  CREATING A DICTIONARY CONTAINING ALL THE DATA FROM UNSRT CLASS  ##
## 																    ##
######################################################################
manipulationDict = {}
selection = []
Cholesky_list = []
zeta_list = []
activeParameters = []
nominal_list = []
rxn_list = []
rIndex = []
for rxn in unsrt_data:
	#print(i)
	rxn_list.append(rxn)
	selection.extend(unsrt_data[rxn].selection)
	Cholesky_list.append(unsrt_data[rxn].cholskyDeCorrelateMat)
	zeta_list.append(unsrt_data[rxn].perturb_factor)
	activeParameters.extend(unsrt_data[rxn].activeParameters)
	nominal_list.append(unsrt_data[rxn].nominal)
	rIndex.append(unsrt_data[rxn].rIndex)
	
manipulationDict["selection"] = deepcopy(selection)#.deepcopy()
manipulationDict["Cholesky"] = deepcopy(Cholesky_list)#.deepcopy()
manipulationDict["zeta"] = deepcopy(zeta_list)#.deepcopy()
manipulationDict["activeParameters"] = deepcopy(activeParameters)#.deepcopy()
manipulationDict["nominal"] = deepcopy(nominal_list)#.deepcopy()

print(manipulationDict["activeParameters"])

import os
import numpy as np
import itertools as it
import more_itertools as mit
dict_ = {}
for rxn in unsrt_data:
	dict_[rxn] = []
file_open = open("zeta_guess_values.txt","r").readlines()

for line in file_open:
	list_1 = line.strip("\n").split(",")
	list_1 = [float(i) for i in list_1]
	n = 3
	iterable = list_1
	list_2 = list(mit.chunked(iterable, n))
	list_2 = [np.asarray(i) for i in list_2]
	for ind,rxn in enumerate(unsrt_data):
		dict_[rxn].append(list_2[ind])
	#print(dict_)
	#print(list_2)
	#raise AssertionError("Stop!!")
delta_n = {}
p = {}
V_opt = {}
V = {}#Populating V for unshuffled portion of the design matrix
ch = {}
nominal = {}
p_max = {}
p_min = {}
theta = {}
Temp = {}
d_n = {}
for rxn in dict_:
	ch[rxn] = unsrt_data[rxn].cholskyDeCorrelateMat
	nominal[rxn] = unsrt_data[rxn].nominal
	p = unsrt_data[rxn].nominal
	p_max[rxn] = unsrt_data[rxn].P_max
	p_min[rxn] = unsrt_data[rxn].P_min
	Temp[rxn] = unsrt_data[rxn].temperatures
	Tp = unsrt_data[rxn].temperatures
	Theta_p = np.array([Tp/Tp,np.log(Tp),-1/(Tp)])
	kmax = Theta_p.T.dot(p_max[rxn])
	kmin = Theta_p.T.dot(p_min[rxn])
	ka_o = Theta_p.T.dot(nominal[rxn])
	zeta_list = dict_[rxn]
	fig = plt.figure()
	plt.title(str(rxn))
	plt.xlabel(r"1000/T\K$^{-1}$")
	plt.ylabel(r"$log_{10}(k)$ / $s^{-1}$ or $log_{10}$(k) / $cm^{3}\,molecule^{-1}\,s^{-1}$")
	plt.plot(1/Tp,kmax,'k--',label="Uncertainty limits")
	plt.plot(1/Tp,kmin,'k--')
	plt.plot(1/Tp,ka_o,'b-',label='Prior rate constant')
	for zeta in zeta_list:
		p_zet = p + np.asarray(np.dot(ch[rxn],zeta)).flatten();
		k =  Theta_p.T.dot(p_zet)
		plt.plot(1/Tp,k,'r-',linewidth=0.55)
	plt.legend()
	plt.savefig("Plots/reaction_"+str(rxn)+".png",bbox_inches='tight')
