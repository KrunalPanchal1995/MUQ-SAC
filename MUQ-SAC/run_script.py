#python default modules
import numpy as np
import scipy as sp
import scipy.stats as stats
from scipy.optimize import minimize
import os, sys, re, threading, subprocess, time
from sklearn.preprocessing import PolynomialFeatures
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
import combustion_target_class
#import data_management
#import data_management as dm
#import simulation_manager
import Uncertainty as uncertainty
#import MechanismManipulator
#import plotter
#import Input_file_reader
#import statistics
#import ParallelWorkers as pk
from mpire import WorkerPool
#import ResponseSurface as PRS
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
	t = combustion_target_class.combustion_target(target,addendum,c_index)
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
















raise AssertionError("The Target class and Uncertainty class")
"""
Selection of the parameters
"""
print(activeParameters)
selectedParams = deepcopy(activeParameters)#.deepcopy()

unsrtDatabase = {}
for index,rxn in enumerate(reaction_index):
	print(rxn)
	unsrtDatabase[rxn] = unsrt_data[rxn].getDtList()

for index,rxn in enumerate(plog_boundary_index):
	unsrtDatabase[rxn] = unsrt_data[rxn].getDtList()

for index,rxn in enumerate(plog_index):
	unsrtDatabase[rxn] = unsrt_data[rxn].getDtList()
		
for i in fallOffCurve_index:
	unsrtDatabase[i] = unsrt_data[i].getDtList()

for i in thirdBody_index:
	unsrtDatabase[i] = unsrt_data[i].getDtList()
	temp = []
	for j in unsrt_data[i].branches.split(","):
		temp.append(j)
	thirdBody_dict[i] = temp
	total_m_params.append(temp)
	
for i in thermo_index:
	unsrtDatabase[i] = unsrt_data[i].getDtList()

for i in transport_index:
	unsrtDatabase[i] = unsrt_data[i].getDtList()	
	

raise AssertionError("The Target class and Uncertainty class")
#
"""
Get the samples from the simulation manager
"""
#sample_length = 100
#pk.Sample(rxnUnsrt_data,reaction_index,sample_length).get_SAMAP(samap_executable,jpdap_executable)
#pk.Sample(rxnUnsrt_data,reaction_index,sample_length).getSamples()


"""
Get the samples from the csv file
"""

cm = 1/2.54 
ZetaFile = open(optInputs["Locations"]["zeta_file"],"r").readlines()
x_data = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in ZetaFile])

y_data = []
"""
Send the samples to simulator
"""
#print(activeParameters)
activeIndexDict = {'H+O2<=>O+OH': 1, 'O+H2<=>H+OH:A': 0, 'O+H2<=>H+OH:B': 0, 'OH+H2<=>H+H2O': 1, '2OH<=>O+H2O': 1, 'O+H+M<=>OH+M:High': 1, 'H+O2(+M)<=>HO2(+M):Low': 1, 'H+O2(+M)<=>HO2(+M):High': 1, 'HO2+H<=>H2+O2': 1, 'HO2+H<=>2OH': 1, 'HO2+OH<=>H2O+O2:A': 1, 'HO2+OH<=>H2O+O2:B': 1, '2HO2<=>H2O2+O2:A': 0, '2HO2<=>H2O2+O2:B': 0, 'H2O2+H<=>OH+H2O': 0, 'H2O2+H<=>HO2+H2': 1, 'CO+O(+M)<=>CO2(+M):Low': 0, 'CO+O(+M)<=>CO2(+M):High': 0, 'CO+OH<=>H+CO2:A': 1, 'CO+OH<=>H+CO2:B': 1, 'CO+HO2<=>OH+CO2': 0, 'HCO+M<=>H+CO+M:High': 1, 'HCO+H<=>H2+CO': 0, 'HCO+O2<=>HO2+CO': 0, 'H2O2+OH<=>H2O+HO2:A': 0, 'H2O2+OH<=>H2O+HO2:B': 0}
{0: {'R1_A': 1, 'R1_n': 1, 'R1_Ea': 1, 'R3a_A': 0, 'R3a_n': 0, 'R3a_Ea': 0, 'R3b_A': 0, 'R3b_n': 0, 'R3b_Ea': 0, 'R4_A': 1, 'R4_n': 1, 'R4_Ea': 1, 'R5_A': 1, 'R5_n': 1, 'R5_Ea': 1, 'R6:High_A': 1, 'R6:High_n': 1, 'R6:High_Ea': 1, 'R2:Low_A': 1, 'R2:Low_n': 1, 'R2:Low_Ea': 1, 'R2:High_A': 1, 'R2:High_n': 1, 'R2:High_Ea': 1, 'R7:B1_A': 1, 'R7:B1_n': 1, 'R7:B1_Ea': 1, 'R7:B2_A': 1, 'R7:B2_n': 1, 'R7:B2_Ea': 1, 'R8a_A': 1, 'R8a_n': 1, 'R8a_Ea': 1, 'R8b_A': 1, 'R8b_n': 1, 'R8b_Ea': 1, 'R9a_A': 0, 'R9a_n': 0, 'R9a_Ea': 0, 'R9b_A': 0, 'R9b_n': 0, 'R9b_Ea': 0, 'R10:B1_A': 0, 'R10:B1_n': 0, 'R10:B1_Ea': 0, 'R10:B2_A': 1, 'R10:B2_n': 1, 'R10:B2_Ea': 1, 'R12:High_A': 0, 'R12:High_n': 0, 'R12:High_Ea': 0, 'R12:Low_A': 0, 'R12:Low_n': 0, 'R12:Low_Ea': 0, 'R13a_A': 1, 'R13a_n': 1, 'R13a_Ea': 1, 'R13b_A': 1, 'R13b_n': 1, 'R13b_Ea': 1, 'R14_A': 0, 'R14_n': 0, 'R14_Ea': 0, 'R15:High_A': 1, 'R15:High_n': 1, 'R15:High_Ea': 1, 'R16_A': 0, 'R16_n': 0, 'R16_Ea': 0, 'R17_A': 0, 'R17_n': 0, 'R17_Ea': 0, 'R11a_A': 0, 'R11a_n': 0, 'R11a_Ea': 0, 'R11b_A': 0, 'R11b_n': 0, 'R11b_Ea': 0}}

activeReactions = {0: {'H+O2<=>O+OH': 1, 'O+H2<=>H+OH:A': 1, 'O+H2<=>H+OH:B': 1, 'OH+H2<=>H+H2O': 1, '2OH<=>O+H2O': 1, 'O+H+M<=>OH+M:High': 1, 'H+O2(+M)<=>HO2(+M):Low': 1, 'H+O2(+M)<=>HO2(+M):High': 1, 'HO2+H<=>H2+O2': 1, 'HO2+H<=>2OH': 1, 'HO2+OH<=>H2O+O2:A': 1, 'HO2+OH<=>H2O+O2:B': 1, '2HO2<=>H2O2+O2:A': 1, '2HO2<=>H2O2+O2:B': 1, 'H2O2+H<=>OH+H2O': 1, 'H2O2+H<=>HO2+H2': 1, 'CO+O(+M)<=>CO2(+M):Low': 1, 'CO+O(+M)<=>CO2(+M):High': 1, 'CO+OH<=>H+CO2:A': 1, 'CO+OH<=>H+CO2:B': 1, 'CO+HO2<=>OH+CO2': 1, 'HCO+M<=>H+CO+M:High': 1, 'HCO+H<=>H2+CO': 1, 'HCO+O2<=>HO2+CO': 1, 'H2O2+OH<=>H2O+HO2:A': 1, 'H2O2+OH<=>H2O+HO2:B': 1}, 1: {'H+O2<=>O+OH': 1, 'O+H2<=>H+OH:A': 1, 'O+H2<=>H+OH:B': 1, 'OH+H2<=>H+H2O': 1, '2OH<=>O+H2O': 1, 'O+H+M<=>OH+M:High': 1, 'H+O2(+M)<=>HO2(+M):Low': 1, 'H+O2(+M)<=>HO2(+M):High': 1, 'HO2+H<=>H2+O2': 1, 'HO2+H<=>2OH': 1, 'HO2+OH<=>H2O+O2:A': 1, 'HO2+OH<=>H2O+O2:B': 1, '2HO2<=>H2O2+O2:A': 1, '2HO2<=>H2O2+O2:B': 1, 'H2O2+H<=>OH+H2O': 1, 'H2O2+H<=>HO2+H2': 1, 'CO+O(+M)<=>CO2(+M):Low': 1, 'CO+O(+M)<=>CO2(+M):High': 1, 'CO+OH<=>H+CO2:A': 1, 'CO+OH<=>H+CO2:B': 1, 'CO+HO2<=>OH+CO2': 1, 'HCO+M<=>H+CO+M:High': 1, 'HCO+H<=>H2+CO': 1, 'HCO+O2<=>HO2+CO': 1, 'H2O2+OH<=>H2O+HO2:A': 1, 'H2O2+OH<=>H2O+HO2:B': 1}}

"""
This routine is to test the performance of response surface in optimization code
"""

def run_sampling(sample,data,generator,length):
	A = Uncertainty.UncertaintyExtractor(data)
	#print(generator)
	a1 = generator[0]
	a2 = generator[1]
	
	#a1 = np.random.uniform(-1,1,1)
	#a2 = np.random.uniform(-1,1,1)
	#a3 = np.random.uniform(-1,1,1)
	#A.generator.append([a1[0],a2[0],a3[0]])
	#if a1 == float(0) and a2 == float(0) and a3 == float(0):
	#	zeta = [0,0,0]
	#else:
	A.populateValues(a1,a2)
	A.getCovariance(flag=False)
	A.getUnCorrelated(flag = False)
	#zeta = A.getLinearCombinationZeta()
	zeta = A.getB2Zeta(flag=True)
	#print(data)
	#zeta =UncertaintyExtractor(data).getExtremeCurves(sample)	
	del A
	#print(zeta)
	
	return (sample,generator,zeta,length)
	
def run_sampling_c(sample,data,generator,length):
	A = Uncertainty.UncertaintyExtractor(data)
	#print(generator)
	a1 = generator[0]
	a2 = generator[1]
	#print(a1,a2,a3)
	#a1 = np.random.uniform(-1,1,1)
	#a2 = np.random.uniform(-1,1,1)
	#a3 = np.random.uniform(-1,1,1)
	#A.generator.append([a1[0],a2[0],a3[0]])
	#if a1 == float(0) and a2 == float(0) and a3 == float(0):
	#	zeta = [0,0,0]
	#else:
	A.populateValues(a1,a2)
	A.getCovariance(flag=False)
	A.getUnCorrelated(flag = False)
	#zeta = A.getLinearCombinationZeta()
	zeta = A.getC2Zeta(flag=True)
	#print(data)
	#zeta =UncertaintyExtractor(data).getExtremeCurves(sample)	
	del A
	#print(zeta)
	
	return (sample,generator,zeta,length)

def run_sampling_b(sample,data,generator,length):
	A = Uncertainty.UncertaintyExtractor(data)
	#print(generator)
	a1 = generator[0]
	a2 = generator[1]
	#print(a1,a2)
	#a1 = np.random.uniform(-1,1,1)
	#a2 = np.random.uniform(-1,1,1)
	#a3 = np.random.uniform(-1,1,1)
	#A.generator.append([a1[0],a2[0],a3[0]])
	#if a1 == float(0) and a2 == float(0) and a3 == float(0):
	#	zeta = [0,0,0]
	#else:
	A.populateValues(a1,a2)
	A.getCovariance(flag=False)
	A.getUnCorrelated(flag = False)
	#zeta = A.getLinearCombinationZeta()
	zeta = A.getB2Zeta(flag=True)
	#print(data)
	#zeta =UncertaintyExtractor(data).getExtremeCurves(sample)	
	del A
	#print(zeta)
	
	return (sample,generator,zeta,length)

def mpire_run_sampling_b(sample,data,generator,length):
	#print(d)
	#sample = d["sample"]
	#data = d['data']
	#generator = d['generator']
	#length = d['length']
	A = Uncertainty.UncertaintyExtractor(data)
	#print(generator)
	a1 = generator[0]
	a2 = generator[1]
	#print(a1,a2)
	#a1 = np.random.uniform(-1,1,1)
	#a2 = np.random.uniform(-1,1,1)
	#a3 = np.random.uniform(-1,1,1)
	#A.generator.append([a1[0],a2[0],a3[0]])
	#if a1 == float(0) and a2 == float(0) and a3 == float(0):
	#	zeta = [0,0,0]
	#else:
	A.populateValues(a1,a2)
	A.getCovariance(flag=False)
	A.getUnCorrelated(flag = False)
	#zeta = A.getLinearCombinationZeta()
	zeta = A.getB2Zeta(flag=True)
	#print(data)
	#zeta =UncertaintyExtractor(data).getExtremeCurves(sample)	
	del A
	#print(zeta)
	
	return (sample,generator,zeta,length)

def run_sampling_direct(sample,rxn,data,generator,length):
	A = Uncertainty.UncertaintyExtractor(data)
	#print(generator)
	a1 = generator[0]
	a2 = generator[1]
	#print(a1,a2,a3)
	#a1 = np.random.uniform(-1,1,1)
	#a2 = np.random.uniform(-1,1,1)
	#a3 = np.random.uniform(-1,1,1)
	#A.generator.append([a1[0],a2[0],a3[0]])
	#if a1 == float(0) and a2 == float(0) and a3 == float(0):
	#	zeta = [0,0,0]
	#else:
	A.populateValues(a1,a2)
	A.getCovariance(flag=False)
	A.getUnCorrelated(flag = False)
	zeta = A.getB2Zeta(flag=True)
	#print(data)
	#zeta =UncertaintyExtractor(data).getExtremeCurves(sample)	
	del A
	#print(zeta)
	
	return (sample,rxn,generator,zeta,length)

def run_sampling_kappa_direct(sample,rxn,data,generator,length):
	A = Uncertainty.UncertaintyExtractor(data)
	A.getCovariance(flag=False)
	A.getUnCorrelated(flag = False)
	zeta = A.getZeta_typeA(generator)
	del A	
	return (sample,rxn,generator,zeta,length)

def run_generate_dir(location):
	os.mkdir(location)
	os.mkdir(location+"/output")
	
def run_executable_files(args):
	os.chdir(args[0])
	#file_paste = open("run","w")
	#file_paste.write(string)
	#file_paste.close()
	#print(location[i])	
	#subprocess.check_call([str(location[i])+"/"+file_name])
	subprocess.call(["./"+args[1]])
	return (args[0],args[2])

def run_map(params):
	location = str(params[4])
	yaml_string = yaml.safe_dump(params[0],default_flow_style=False)
	yamlfile = open(location+"/mechanism.yaml","w").write(yaml_string)
	sim = open(location+"/cantera_.py",'w').write(params[1])
	sim = open(location+"/FlameMaster.input",'w').write(params[1])
	extract = open(location+"/extract.py",'w').write(params[5])
	runConvertorScript = open(location+"/run_convertor",'w').write(params[2])
	runScript = open(location+"/run","w").write(params[3])
	subprocess.call(["chmod","+x",location+"/run_convertor"])
	subprocess.call(["chmod","+x",location+"/run"])
	del yaml_string
	del location
	#return (params[4])

def run_direct_map(file_dict):
	location = str(file_dict["mkdir"])
	#print(location)
	betaFile = open(location+"/beta.txt",'w').write(str(file_dict["beta"]))
	#mechFile = open(location+"/mechanism.mech",'w').write(file_dict["mechanism"])
	yaml_string = yaml.safe_dump(file_dict["mechanism"],default_flow_style=False)
	yamlfile = open(location+"/mechanism.yaml","w").write(yaml_string)
	del yaml_string
	#thermoFile = open(location+"/thermo.therm",'w').write(file_dict["thermo"])
	#transportFile = open(location+"/transport.trans",'w').write(file_dict["transport"])
	sim = open(location+"/cantera_.py",'w').write(file_dict["simulationInputString"])
	sim = open(location+"/FlameMaster.input",'w').write(file_dict["simulationInputString"])
	extract = open(location+"/extract.py",'w').write(file_dict["extractString"])
	runConvertorScript = open(location+"/run_convertor",'w').write(file_dict["file_convertor_script"])
	runScript = open(location+"/run","w").write(file_dict["run_script"])
	#print(location+"/run_convertor")
	subprocess.call(["chmod","+x",location+"/run_convertor"])
	subprocess.call(["chmod","+x",location+"/run"])
	del location
	
def run_copy_files(i,sim_dict,mech_dict,thermo_dict,trans_dict,run_convert_dict,run_dict,locations,n):
	#os.mkdir(location)
	location = locations[i]
	#print(location)
	os.chdir(location)
	NewMechanismString = mech_dict[i]
	NewTransportString = trans_dict[i]
	NewThermoString = thermo_dict[i]
	NewRunConvertorScript = run_convert_dict[i]
	NewRunScript = run_dict[i]
	NewSimScript = sim_dict[i]
	sim = open("cantera_.py",'w')
	sim.write(NewSimScript)
	sim.close()
	mechFile = open("mechanism.mech",'w')
	mechFile.write(NewMechanismString)
	mechFile.close()
	thermoFile = open("thermo.therm",'w')
	thermoFile.write(NewThermoString)
	thermoFile.close()
	transportFile = open("transport.trans",'w')
	transportFile.write(NewTransportString)
	transportFile.close()
	runConvertorScript = open("run_convertor",'w')
	runConvertorScript.write(NewRunConvertorScript)
	runConvertorScript.close()
	runScript = open("run","w")
	runScript.write(NewRunScript)
	runScript.close()
	subprocess.call(["chmod","+x","run_convertor"])
	subprocess.call(["chmod","+x","run"])
	os.mkdir("output")
	return(i,n)
	
def create_files_for_directSimulation(file_dict,n):
	#os.mkdir(location)
	location = file_dict["mkdir"]
	#print(location)
	os.chdir(location)
	NewMechanismString = file_dict["mechanism"]
	NewTransportString = file_dict["transport"]
	NewThermoString = file_dict["thermo"]
	NewRunConvertorScript = file_dict["file_convertor_script"]
	NewRunScript = file_dict["run_script"]
	NewSimScript = file_dict["simulationInputString"]
	sim = open("cantera_.py",'w')
	sim.write(NewSimScript)
	sim.close()
	mechFile = open("mechanism.mech",'w')
	mechFile.write(NewMechanismString)
	mechFile.close()
	thermoFile = open("thermo.therm",'w')
	thermoFile.write(NewThermoString)
	thermoFile.close()
	transportFile = open("transport.trans",'w')
	transportFile.write(NewTransportString)
	transportFile.close()
	runConvertorScript = open("run_convertor",'w')
	runConvertorScript.write(NewRunConvertorScript)
	runConvertorScript.close()
	runScript = open("run","w")
	runScript.write(NewRunScript)
	runScript.close()
	subprocess.call(["chmod","+x","run_convertor"])
	subprocess.call(["chmod","+x","run"])
	os.mkdir("output")
	return(dictonary,n)
	
#class Ray():
#	def __init__(self,workers):
		

class Worker():
	def __init__(self, workers):
		self.pool_mpire = WorkerPool(n_jobs=workers)
		self.pool = multiprocessing.Pool(processes=workers)
		self.progress = []
		self.parallized_zeta = []
		self.generator = []
		self.parallel_zeta_dict = {}
		
	def callback_run(self, result):
		#print(result)
		string=""
		for i in result:
			self.progress.append(i[0])
			string+=f"{i[0]}/run\n"
			total = i[1]
		#string+="\n"
		sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(self.progress)/float(total)*100))
		sys.stdout.flush()
		f = open('../progress','+a')
		f.write(string)
		f.close()
		
	def callback_runConvertor(self, result):
		self.progress.append(result[0])
		sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(self.progress)/float(result[-1])*100))
		sys.stdout.flush()
		#self.pool.terminate()
	def callback_create(self, result):
		self.progress.append(result[0])
		sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(self.progress)/float(result[-1])*100))
		sys.stdout.flush()
		
	def callback_create_1(self, result):
		#string=""
		#for i in result:
		#	self.progress.append(i[0])
			#string+=f"{i[0]}\n"
			#total = i[1]
		#self.progress.append(result[0])
		sys.stdout.write("\t\t\r{} is complete".format(len(result)))
		sys.stdout.flush()
		#self.pool.terminate()

	def callback(self,result):
		#print("Entered callback\n")
		self.progress.append(result[0])
		self.generator.append(result[1])
		self.parallized_zeta.append(result[2])
		sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(self.progress)/float(result[-1])*100))
		sys.stdout.flush()
	
	def callback_direct(self,result):
		#print("Entered callback\n")
		self.progress.append(result[0])
		self.generator.append(result[2])
		self.parallized_zeta.append(result[3])
		self.parallel_zeta_dict[result[1]] = result[3]
		sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(self.progress)/float(result[-1])*100))
		sys.stdout.flush()
	
	def callback_kappa_direct(self,result):
		#print("Entered callback\n")
		self.progress.append(result[0])
		self.generator.append(result[2])
		self.parallized_zeta.append(result[3])
		self.parallel_zeta_dict[result[1]] = result[3]
		sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(self.progress)/float(result[-1])*100))
		sys.stdout.flush()
		
	def callback_error(self,result):
    		print('error', result)
	def do_job_async_unsrt_direct(self,data,generator,beta):
		#print("Entered Async\n")
		#print("Entered async\n")
		for args in data:
			self.pool.apply_async(run_sampling_direct, 
				  args=(1,args,data[args],generator[args],beta,), 
				  callback=self.callback_direct)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
		return self.parallel_zeta_dict
	
	def do_job_async_unsrt_kappa_direct(self,data,generator,beta):
		#print("Entered Async\n")
		#print("Entered async\n")
		for args in data:
			self.pool.apply_async(run_sampling_kappa_direct, 
				  args=(1,args,data[args],generator[args],beta,), 
				  callback=self.callback_kappa_direct,error_callback=self.callback_error)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
		return self.parallel_zeta_dict
	
	def do_unsrt(self,data,sampling_points):
		#generators = data["generators"]
		for args in range(sampling_points):
			self.pool.apply_async(run_sampling, 
				  args=(1,data,data["generators"][args],sampling_points,), 
				  callback=self.callback)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
		return self.generator,self.parallized_zeta
	
	def do_unsrt_c(self,data,sampling_points):
		#generators = data["generators"]
		for args in range(sampling_points):
			self.pool.apply_async(run_sampling_c, 
				  args=(1,data,data["generators_c"][args],sampling_points,), 
				  callback=self.callback,error_callback=self.callback_error)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
		return self.generator,self.parallized_zeta
	
	
	def MPIRE_do_unsrt_b(self,data,sampling_points):
		list_compre = [({"sample": 1,"data":data,"generator":data["generators_b"][args],"length":sampling_points}) for args in range(sampling_points)]
		for result in self.pool_mpire.imap(mpire_run_sampling_b,list_compre):
			#print(result)
			self.progress.append(result[0])
			self.generator.append(result[1])
			self.parallized_zeta.append(result[2])
			sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(self.progress)/float(result[-1])*100))
			sys.stdout.flush()
		self.pool_mpire.terminate()
		return self.generator,self.parallized_zeta
				
	
	def do_unsrt_b(self,data,sampling_points):
		for args in range(sampling_points):
			self.pool.apply_async(run_sampling_b, 
				  args=(1,data,data["generators_b"][args],sampling_points,), 
				  callback=self.callback)#,error_callback=self.callback_error)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
		return self.generator,self.parallized_zeta
	
	def do_job_async_convertor(self,args):
			#self.pool.apply_async(run_executable_files, 
			#	  args=(args,file_name,len(location)), 
			#	  callback=self.callback_runConvertor,error_callback=self.callback_error)
		self.pool.map_async(run_executable_files,args,callback=self.callback_run)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
	
	def do_job_async(self,args):

			#self.pool.apply_async(run_executable_files, 
			#	  args=(arg,file_name,len(location)), 
			#	  callback=self.callback_run,error_callback=self.callback_error)
		self.pool.map_async(run_executable_files,args,callback=self.callback_run,error_callback=self.callback_error)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
		
	def do_job_create_async(self,sim_list,sim_dict,mech_dict,thermo_dict,trans_dict,run_convert_dict,run_dict,locations):
		for args in sim_list:
			self.pool.apply_async(run_copy_files, 
				  args=(args,sim_dict,mech_dict,thermo_dict,trans_dict,run_convert_dict,run_dict,locations,len(sim_list)), 
				  callback=self.callback_create,error_callback=self.callback_error)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
	
	def do_job_direct_async(dictionary_list):
		for args in dictionary_list:
			self.pool.apply_async(create_files_for_directSimulation, 
				  args=(args,len(dictionary_list)), 
				  callback=self.callback_create,error_callback=self.callback_error)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
	
	def do_job_direct_map(self,params):
		self.pool.map_async(run_direct_map,params)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
		
	def do_job_map(self,locations):
		self.pool.map_async(run_generate_dir,locations)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
	
	def do_job_map_create(self,params):
		self.pool.map_async(run_map,params,callback=self.callback_create_1)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
	
	def do_job_executor(self,locations):
		with concurrent.futures.ProcessPoolExecutor() as executor:
			executor.map(run_generate_dir,locations)
		
	






sim_type = "opt"
simulator = simulation_manager.SM(optInputs,iFile,case_dir,rps_order,activeParameters, reaction_index,fallOffCurve_index,thirdBody_index,thermo_index,transport_index, mech_file_location,fileType,unsrt_data, rxnUnsrt_data, focUnsrt_data, tbdUnsrt_data, thermoUnsrt_data, transportUnsrt_data, target_list, fuel, global_reaction, thermo_file_location, trans_file_location,startProfile_location,design_type,sim_type,parallel_threads,selectedParams,manipulationDict,activeReactions,extra_arg = activeIndexDict)
#sample_length = 100
#pk.Sample(rxnUnsrt_data,reaction_index,sample_length).get_SAMAP(samap_executable,jpdap_executable)
case1 = target_list[0]
case2 = target_list[1]
iter_number = 0
objective = 1
cases = {}
l = []
for ind,i in enumerate(target_list):
	l.append(ind)
	cases[ind] = target_list[ind]
if "yData.txt" not in os.listdir():
	#y_data = simulator.getSimulatedValues(x_data)
	if "DirectSimulation" in os.listdir():
		os.chdir("DirectSimulation")
		for i,ele in enumerate(x_data):
			#print(ele)
			#print(len(x_data))
			#raise AssertionError("!00")
			os.chdir(f"{i}")
			print(i)
			temp_output = []
			for j,target in enumerate(target_list):
				eta, out_path = dm.extract_output(target,"N/A",f"case-{j}/output/",i)#index is dummy, fuel is also not useful in this script, path and case are important
				temp_output.append(eta)
			os.chdir("..")
			y_data.append(temp_output)
		os.chdir("..")
			
		iString = ""
		for i in y_data:
			for j in i:
				iString+=f"{float(j)},"
			iString+="\n"
		write_yData = open("yData.txt","w").write(iString)
	
	else:
		mkdir_list = []
		dictionary_list = []
		args_convert = []
		args = []
		for x in x_data[0:5]:
			#print(x)
			y_data.append(simulator.getSimulatedValues(x,l,cases,iter_number,objective))
			#args_c,args_r = simulator.getSimulatedValues(x,l,cases,iter_number,objective)
			#mkdir_list.extend(mkdir)
			#dictionary_list.extend(dict_lst)
			#args_convert.extend(args_c)
			#args.extend(args_r)
			#args.extend(simulator.getSimulatedValues(x,l,cases,iter_number,objective))
			iter_number+=1
		
		#W = Worker(int(parallel_threads))
		#W.do_job_map(mkdir_list)
		#print("\tDirectories for direct simulations are generated\n")
		
		#V = Worker(int(parallel_threads))
		#print(params[0][6])
		#V.do_job_direct_map(dictionary_list)
		#print("\tRequired files for {} iteration is generated\n".format(iter_number))
		
		#U = Worker(int(parallel_threads))
		#U.do_job_async_convertor(args_convert)
		
		#print("\tFiles converted to standard input files for iteration \n".format(iter_number))
		#del U
		
		#start_time = time.time()
		#X = Worker(int(parallel_threads))
		#X.do_job_async(args)
		
		#print("\tSimulations for {} iteration is Done!!".format(iter_number))
			
		#dt = int(time.time() - start_time)
		#hours = dt/3600
		#minutes = (dt%3600)/60
		#seconds = dt%60
		#os.system("clear")
		#print("Performed {} Simulations....".format(len(args)))
		#print("Time for performing simulations : {h} hours,  {m} minutes, {s} seconds\n................................................ \n".format(h = hours, m = minutes, s =seconds))
		#del X
		
		#os.chdir("DirectSimulation")
		#for i,ele in enumerate(x_data):
		#	os.chdir(f"{i}")
		#	temp_output = []
		#	for j,target in enumerate(target_list):
		#		eta, out_path = dm.extract_output(target,"N/A",f"case-{j}/output/",i)#index is dummy, fuel is also not useful in this script, path and case are important
		#		temp_output.append(eta)
		#	os.chdir("..")
		#	y_data.append(temp_output)
		#os.chdir("..")
		iString = ""
		for i in y_data:
			iString+=f"{i}\n"
		
		write_yData = open("yData.txt","w").write(iString)
else:
	yDATA = open("yData.txt","r").readlines()
	y_data = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in yDATA])
	#for i in yDATA:
	#	y_data.append(float(i))
#print(len(y_data))
y_dict = {}

def statistics_(error):
	h = list(error)
	h.sort()
	hmean = np.mean(h)
	hstd = np.std(h)
	return h,hmean,hstd


for i in range(targets_count):
	temp = []
	for j in range(len(y_data)):
		temp.append(y_data[j][i])
	y_dict[i] = temp
#print(y_dict)
"""
Now test the respones surface 
"""
string = "Target,training_max,training_mean,testing_max,testing_mean,search_max,search_mean\n"
for z in range(targets_count):
	cases = target_list[z].target
	X_train = optInputs["Locations"]["xData"]+"/Beta_list_case-"+str(z)+".csv"

	Y_train = optInputs["Locations"]["yData"]+"/sim_data_case-"+str(z)+".lst"
	
	X_test = optInputs["Locations"]["xTest"]+"/Beta_list_case-"+str(z)+".csv"

	Y_test = optInputs["Locations"]["yTest"]+"/sim_data_case-"+str(z)+".lst"
	
	file_x = open(X_train,"r").readlines()
	file_y = open(Y_train,"r").readlines()
	#print(len(file_y))
	x_train_data = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in file_x])
	y_train_data = np.exp(np.asarray([float(i.strip("\n").split()[1]) for i in file_y]))
	
	#y_train_data = np.log(np.asarray([float(i.strip("\n").split()[1]) for i in file_y]))
	
	
	file_x_test = open(X_test,"r").readlines()
	file_y_test = open(Y_test,"r").readlines()

	x_test_data = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in file_x_test])
	y_test_data = np.exp((np.asarray([float(i.strip("\n").split()[1]) for i in file_y_test])))
	
	y_data_test = np.exp((y_dict[z]))
	
	#y_test_data = np.log(np.asarray([float(i.strip("\n").split()[1]) for i in file_y_test]))
	
	#y_data_test = np.log(y_dict[z])

	x_data_test = x_data[0:5]
	#for row in x_data: 
	#	count = 0
	#	temp = []
	#	for i in activeParameters:
	#		if activeReactions[0][i] == 1:
	#			temp.append(row[count])
	#			temp.append(row[count+1])
	#			temp.append(row[count+2])
	#			count+=3
	#		else:
	#			count+=3
	#	x_data_test.append(np.asarray(temp))


	#print(np.shape(x_train_data))
	#print(np.shape(x_train_data))
	#print(np.shape(x_data_test))
	#print(np.shape(y_data_test))

	#print(len(x_train_data[0]))
	#print(len(x_data_test[0]))
	model = PRS.ResponseSurface(x_train_data,y_train_data,2)
	
	"""
	SVR_response = model.create_SVR_response_surface()
	### testing of the response surface #####
	svr_yfit_test = model.svr.predict(x_data_test)
	error_svr = [abs(model.svr_yfit[i_i]-y_data[i_i])*100/(y_data[i_i]) for i_i,i in enumerate(y_data)]
	max_error_svr = max(error_svr)
	mean_error_svr = statistics.mean(error_svr)

	error_test_svr = [abs(svr_yfit_test[i_i]-i)*100/(i) for i_i,i in enumerate(y_data_test)]
	max_error_test_svr = max(error_test_svr)
	mean_error_test_svr = statistics.mean(error_test_svr)
	"""
	######
	#####   2] General PRS
	######
	######
	######

	general_response = model.create_response_surface()
	y_prediction = model.resFramWrk
	error_general = [abs(model.resFramWrk[i_i]-i)*100/(i) for i_i,i in enumerate(y_train_data)]
	max_error_general = max(error_general)
	#print(error_general)
	mean_error_general = statistics.mean(error_general)
	
	"""
	Testing the response surface using the zetas generated from optimization algorithm
	"""
	general_yfit_test = []
	for i in x_data_test:
		general_yfit_test.append(model.evaluate(i))
	#print(len(general_yfit_test))
	#print(len(y_data_test))
	error_test_general = [abs(general_yfit_test[i_i]-i)*100/(i) for i_i,i in enumerate(y_data_test)]
	residual_error_general = [abs(general_yfit_test[i_i]-i) for i_i,i in enumerate(y_data_test)]
	max_error_test_general = max(error_test_general)
	mean_error_test_general = statistics.mean(error_test_general)
	
	max_error_residual_test_general = max(residual_error_general)
	mean_error_residual_test_general = statistics.mean(residual_error_general)

	
	"""
	Testing the response surface using the zetas generated from function
	"""
	general_yfit = []
	for i in x_test_data:
		general_yfit.append(model.evaluate(i))
	#print(general_yfit_test)
	error_test = [abs(general_yfit[i_i]-i)*100/(i) for i_i,i in enumerate(y_test_data)]
	residual_error_test = [abs(general_yfit[i_i]-i) for i_i,i in enumerate(y_test_data)]
	max_error_test = max(error_test)
	mean_error_test = statistics.mean(error_test)
	
	max_error_residual_test = max(residual_error_test)
	mean_error_residual_test = statistics.mean(residual_error_test)
	
	
	"""
	HuberRegression
	"""
	"""
	
	Huber_response = model.create_Isotonic_response_surface()
	### testing of the response surface #####
	huber_yfit_test = model.huber.predict(x_data_test)
	#print(model.huber_yfit)
	error_huber = [abs(model.huber_yfit[i_i]-y_train_data[i_i])*100/(y_train_data[i_i]) for i_i,i in enumerate(y_train_data)]
	#print(error_huber)
	max_error_huber = max(error_huber)
	mean_error_huber = statistics.mean(error_huber)

	error_test_huber = [abs(huber_yfit_test[i_i]-i)*100/(i) for i_i,i in enumerate(y_data_test)]
	max_error_test_huber = max(error_test_huber)
	mean_error_test_huber = statistics.mean(error_test_huber)
		
	
	general_yfit_huber_test = []
	general_yfit_huber_test = model.huber.predict(x_data_test)
	#for i in x_data_test:
	#	general_yfit_huber_test.append(model.huber.predict(i))
	#print(general_yfit_test)
	error_huber_test_general = [abs(general_yfit_huber_test[i_i]-i)*100/(i) for i_i,i in enumerate(y_data_test)]
	max_error_huber_test_general = max(error_huber_test_general)
	mean_error_huber_test_general = statistics.mean(error_huber_test_general)
	
	
	general_huber_yfit = []
	general_huber_yfit = model.huber.predict(x_test_data)
	#for i in x_test_data:
	#	general_huber_yfit.append(model.huber.predict(i))
	#print(general_yfit_test)
	error_huber_test = [abs(general_huber_yfit[i_i]-i)*100/(i) for i_i,i in enumerate(y_test_data)]
	max_error_huber_test = max(error_huber_test)
	mean_error_huber_test = statistics.mean(error_huber_test)
	
	"""
	#####
	### testing of the response surface #####
	#####
	x = np.linspace(-500,500,100)
	fig = plt.figure()
	#plt.title("Training")
	#ax = fig.add_subplot(111)
	plt.plot(x,x,"-",linewidth=0.32)
	plt.ylim(abs(min(general_yfit))*0.95,abs(max(general_yfit))*1.05)
	plt.xlim(abs(min(general_yfit))*0.95,abs(max(general_yfit))*1.05)
	plt.plot(np.asarray(general_yfit_test),np.asarray(y_data_test),"r^",ms=4,label=f"Search Iterations (max error = {max_error_test_general:.3f}%), mean error = {mean_error_test_general:.3f}%)")
	#plt.scatter(np.asarray(model.resFramWrk), np.asarray(y_train_data), color="none", edgecolor="green")
	plt.scatter(np.asarray(general_yfit), np.asarray(y_test_data),s=10,marker="o",color="black", edgecolor="black",label=f"Testing dataset (max error = {max_error_test:.3f}%, mean error = {mean_error_test:.3f}%)")
	
	#plt.scatter(np.asarray(model.resFramWrk), np.asarray(y_train_data),s=8,marker="s",color="none", edgecolor="green",label=f"Training dataset (max error = {max_error_general:.3f}%, mean error = {mean_error_general:.3f}%)")

	plt.xlabel("Direct Simulations",fontsize=10)
	plt.ylabel("Response surface predictions",fontsize=10)
	#plt.plot(y_data,model.svr_yfit,"g.",label="SVR")
	#plt.plot(y_train_data,model.resFramWrk,"r.",label="Traning")
	#plt.plot(y_data_test,svr_yfit_test,"g.",label ="SVR")
	#plt.plot(y_data_test,general_yfit_test,"g.",label="Testing")
	plt.legend()
	#plt.show()
	#pk.Sample(rxnUnsrt_data,reaction_index,sample_length).getSamples()
	plt.savefig("Testing"+str(z)+".png",bbox_inches="tight")
	
	string+=f"{cases},{max_error_general},{mean_error_general},{max_error_test},{mean_error_test},{max_error_test_general},{mean_error_test_general}\n"
	"""
	x = np.linspace(-500,500,100)
	fig = plt.figure()
	#plt.title("Training")
	#ax = fig.add_subplot(111)
	plt.plot(x,x,"-",linewidth=0.32)
	plt.ylim(abs(min(general_yfit))*0.95,abs(max(general_yfit))*1.05)
	plt.xlim(abs(min(general_yfit))*0.95,abs(max(general_yfit))*1.05)
	plt.plot(np.asarray(general_yfit_huber_test),np.asarray(y_data_test),"r^",ms=4,label=f"Search Iterations  (max error = {max_error_huber_test_general:.3f}%, mean error = {mean_error_huber_test_general:.3f}%)")
	#plt.scatter(np.asarray(model.resFramWrk), np.asarray(y_train_data), color="none", edgecolor="green")
	plt.scatter(np.asarray(general_huber_yfit), np.asarray(y_test_data),s=9,marker="o",color="black", edgecolor="black",label=f"Testing dataset (max error = {max_error_huber_test:.3f}%, mean error = {mean_error_huber_test:.3f}%)")
	plt.xlabel("Black box simulations")
	plt.ylabel("Response surface predictions (Huber)")
	#plt.plot(y_data,model.svr_yfit,"g.",label="SVR")
	#plt.plot(y_train_data,model.resFramWrk,"r.",label="Traning")
	#plt.plot(y_data_test,svr_yfit_test,"g.",label ="SVR")
	#plt.plot(y_data_test,general_yfit_test,"g.",label="Testing")
	plt.legend()
	#plt.show()
	#pk.Sample(rxnUnsrt_data,reaction_index,sample_length).getSamples()
	plt.savefig("Testing_huber"+str(z)+".pdf",bbox_inches="tight")
	"""
	#string+="\n"
	fig = plt.figure()
	h,hmean,hstd = statistics_(list(y_train_data))
	h_,hmean_,hstd_ = statistics_(list(y_prediction))
	pdf = stats.norm.pdf(h, hmean, hstd)
	pdf_ = stats.norm.pdf(h_, hmean_, hstd_)
	y = np.arange(0,max(pdf),100)
	x = hmean*np.ones(len(y))
	
	fig = plt.figure()
	plt.plot(h, pdf,label="Simulations")
	plt.axvline(x = hmean, color = 'b', label = f'simulation mean ({hmean}) ')
	plt.plot(h_, pdf_,label="Response surface")
	plt.axvline(x = hmean_, color = 'b', label = f'Prediction mean ({hmean_})')
	
	#plt.plot(x,y,"--",label= "Mean")
	plt.legend()
	plt.savefig("Testing_pdf"+str(z)+".pdf",bbox_inches="tight")
file_print =open("New_FLS_cm_percentage.csv","w").write(string)
