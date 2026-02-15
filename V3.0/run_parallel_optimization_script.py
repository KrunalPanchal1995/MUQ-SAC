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
from MechManipulator2_0 import Manipulator
from ParallelOptimizationTool import OptimizationTool as Optimizer


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
MSO = "Multi_Stage_Optimization"
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
MSO_data = optInputs[MSO]
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
#targetLines = open(locations[targets],'r').readlines()
addendum = yaml.safe_load(open(locations[add],'r').read())

############################################################


####################################################
##  FINDING THE LOCATIONS FOR PO OPTIMIZATION	  ##
##                                                ##
####################################################
HTT_OPT_FOLDER = MSO_data["Opt_HTT"]
LTT_OPT_FOLDER = MSO_data["Opt_LTT"]
HTC_file = MSO_data["HTC"]
LTC_file = MSO_data["LTC"]
HTT_file = MSO_data["HTT"]
LTT_file = MSO_data["LTT"] 
targetLinesLTT = open(LTT_file,'r').readlines()
targetLinesHTT = open(HTT_file,'r').readlines()
####################################################
##  Unloading the target data	  	           ##
## TARGET CLASS CONTAINING EACH TARGET AS A CASE  ##
####################################################

"""
Targets for LTT regime
"""

targets_LTT = []
c_index_LTT = 0
string_target_LTT = ""

for target in targetLinesLTT:
	if "#" in target:
		target = target[:target.index('#')]	
	add = deepcopy(addendum)
	t = combustion_target_class.combustion_target(target,add,c_index_LTT)
	string_target_LTT+=f"{t.dataSet_id}|{t.target}|{t.species_dict}|{t.temperature}|{t.pressure}|{t.phi}|{t.observed}|{t.std_dvtn}\n"
	c_index_LTT +=1
	targets_LTT.append(t)
case_dir_LTT = range(0,len(targets_LTT))
print(f"No. of LTT targets: {case_dir_LTT}\n")
target_file = open("target_data_LTT.txt","w")
target_file.write(string_target_LTT)
target_file.close()

"""
Targets for HTT regime
"""


targets_HTT = []
c_index_HTT = 0
string_target_HTT = ""

for target in targetLinesHTT:
	if "#" in target:
		target = target[:target.index('#')]	
	add = deepcopy(addendum)
	t = combustion_target_class.combustion_target(target,add,c_index_HTT)
	string_target_HTT+=f"{t.dataSet_id}|{t.target}|{t.species_dict}|{t.temperature}|{t.pressure}|{t.phi}|{t.observed}|{t.std_dvtn}\n"
	c_index_HTT +=1
	targets_HTT.append(t)
case_dir_HTT = range(0,len(targets_HTT))
print(f"No. of HTT targets: {case_dir_HTT}\n")
print("\n\nOptimization targets identified.\nStarted the MUQ process.......\n")
target_file = open("target_data_HTT.txt","w")
target_file.write(string_target_HTT)
target_file.close()

#################################################################################
opt_dir = os.getcwd()
#################################################################################

############################################
##  Uncertainty Quantification for LTC    ##
##  					                  ##
############################################
os.chdir(LTT_OPT_FOLDER)
# Load the object from the file
with open('unsrt.pkl', 'rb') as file_:
	unsrt_data_LTC = pickle.load(file_)
print("Uncertainty analysis for LTC is already finished\n")


######################################################################
##  CREATING A DICTIONARY CONTAINING ALL THE DATA FROM UNSRT CLASS  ##
##		                                                    ##
######################################################################

manipulationDictLTC = {}
selection_LTC = []
Cholesky_list_LTC = []
zeta_list_LTC = []
activeParameters_LTC = []
nominal_list_LTC = []
rxn_list_LTC = []
rIndex_LTC = []
rxn_list_a_fact_LTC = []

for rxn in unsrt_data_LTC:
	#print(i)
	data_LTC = {}
	data_LTC["temperatures"] = unsrt_data_LTC[rxn].temperatures
	data_LTC["uncertainties"] = unsrt_data_LTC[rxn].uncertainties
	data_LTC["Arrhenius"] = unsrt_data_LTC[rxn].nominal
	#print(rxn,data)
	rxn_list_LTC.append(rxn)
	rxn_list_a_fact_LTC.append(unsrt_data_LTC[rxn].rxn)
	selection_LTC.extend(unsrt_data_LTC[rxn].selection)
	Cholesky_list_LTC.append(unsrt_data_LTC[rxn].cholskyDeCorrelateMat)
	zeta_list_LTC.append(unsrt_data_LTC[rxn].perturb_factor)
	activeParameters_LTC.extend(unsrt_data_LTC[rxn].activeParameters)
	nominal_list_LTC.append(unsrt_data_LTC[rxn].nominal)
	rIndex_LTC.append(unsrt_data_LTC[rxn].rIndex)

rxn_list_a_fact_LTC = set(rxn_list_a_fact_LTC)
rxn_a_fact_LTC = ""
for rxn in rxn_list_a_fact_LTC:
	rxn_a_fact_LTC+=f"{rxn}\n"
file_rxn_a_fact_LTC = open("rxn_list_a_fact.csv","w").write(rxn_a_fact_LTC)	
manipulationDictLTC["selection"] = deepcopy(selection_LTC)#.deepcopy()
manipulationDictLTC["Cholesky"] = deepcopy(Cholesky_list_LTC)#.deepcopy()
manipulationDictLTC["zeta"] = deepcopy(zeta_list_LTC)#.deepcopy()
manipulationDictLTC["activeParameters"] = deepcopy(activeParameters_LTC)#.deepcopy()
manipulationDictLTC["nominal"] = deepcopy(nominal_list_LTC)#.deepcopy()
print(f"\n\n################################################\nThe total reactions choosen for parallel optimization (LTC): {len(manipulationDictLTC['activeParameters'])}\n\tThe list is as follows:\n\t")
print("\t"+f"{manipulationDictLTC['activeParameters']}")
print("\n################################################\n\n")

############################################
##  Uncertainty Quantification for HTC    ##
##  					                  ##
############################################


os.chdir(HTT_OPT_FOLDER)
# Load the object from the file
with open('unsrt.pkl', 'rb') as file_:
	unsrt_data_HTC = pickle.load(file_)
print("Uncertainty analysis for LTC is already finished\n")


######################################################################
##  CREATING A DICTIONARY CONTAINING ALL THE DATA FROM UNSRT CLASS  ##
##		                                                    ##
######################################################################

manipulationDictHTC = {}
selection_HTC = []
Cholesky_list_HTC = []
zeta_list_HTC = []
activeParameters_HTC = []
nominal_list_HTC = []
rxn_list_HTC = []
rIndex_HTC = []
rxn_list_a_fact_HTC = []

for rxn in unsrt_data_HTC:
	#print(i)
	data_HTC = {}
	data_HTC["temperatures"] = unsrt_data_HTC[rxn].temperatures
	data_HTC["uncertainties"] = unsrt_data_HTC[rxn].uncertainties
	data_HTC["Arrhenius"] = unsrt_data_HTC[rxn].nominal
	#print(rxn,data)
	rxn_list_HTC.append(rxn)
	rxn_list_a_fact_HTC.append(unsrt_data_HTC[rxn].rxn)
	selection_HTC.extend(unsrt_data_HTC[rxn].selection)
	Cholesky_list_HTC.append(unsrt_data_HTC[rxn].cholskyDeCorrelateMat)
	zeta_list_HTC.append(unsrt_data_HTC[rxn].perturb_factor)
	activeParameters_HTC.extend(unsrt_data_HTC[rxn].activeParameters)
	nominal_list_HTC.append(unsrt_data_HTC[rxn].nominal)
	rIndex_HTC.append(unsrt_data_HTC[rxn].rIndex)

rxn_list_a_fact_HTC = set(rxn_list_a_fact_HTC)
rxn_a_fact_HTC = ""
for rxn in rxn_list_a_fact_HTC:
	rxn_a_fact_HTC+=f"{rxn}\n"
file_rxn_a_fact_HTC = open("rxn_list_a_fact.csv","w").write(rxn_a_fact_HTC)	
manipulationDictHTC["selection"] = deepcopy(selection_HTC)#.deepcopy()
manipulationDictHTC["Cholesky"] = deepcopy(Cholesky_list_HTC)#.deepcopy()
manipulationDictHTC["zeta"] = deepcopy(zeta_list_HTC)#.deepcopy()
manipulationDictHTC["activeParameters"] = deepcopy(activeParameters_HTC)#.deepcopy()
manipulationDictHTC["nominal"] = deepcopy(nominal_list_HTC)#.deepcopy()
print(f"\n\n################################################\nThe total reactions choosen for parallel optimization (HTC): {len(manipulationDictHTC['activeParameters'])}\n\tThe list is as follows:\n\t")
print("\t"+f"{manipulationDictHTC['activeParameters']}")
print("\n################################################\n\n")

###########################################################
###########################################################
### ResponseSurfaceLTT with reaction classification (LTC)
### 
###########################################################
###########################################################
os.chdir(LTT_OPT_FOLDER)
"""
For Partial PRS (p-PRS), we need to do a second level of sensitivity analysis:	
	
	sensitivity_analysis: a dictionary that stores the sensitivity coefficients of all the active
			      parameters (selected in the first set of sensitivity analysis) reaction wise
			      as in unsrt_data dictionary
"""

if PRS_type == "Partial":
	####################################################
	## Sensitivity Analysis of the selected reactions ##
	####################################################
	
	print("\nPartial Polynomial Response Surface is choosen as a Solution Mapping Technique\n\n\tRecollecting the Sensitivity Analysis: \n\t\tKindly be Patient\n")
	if "A-facto" in design_type:
		status_ = "Pending"
		print("\n\tBrute-force Sensitivity analysis is alsready over!!\n")
		with open('sens_parameters.pkl', 'rb') as file_:
			sensitivity_analysis_LTT = pickle.load(file_)
	#print(sensitivity_analysis)
	#raise AssertionError("Stop")
	partialPRS_Object = []
	selected_params_dict_LTT = {}
	design_matrix_dict = {}
	yaml_loc_dict = {}
	no_of_sim = {}
	print("\n################################################\n###  Starting to generate Design Matrix      ###\n###  for all targets                        ###\n################################################\n\n")
	for case in case_dir_LTT:
		#print(case)
		#PPRS_system = P_PRS.PartialPRS(sensitivity_analysis[str(case)],unsrt_data,optInputs,target_list,case,activeParameters,design_type)
		PPRS_system = P_PRS.PartialPRS(sensitivity_analysis_LTT[str(case)],unsrt_data_LTC,optInputs,targets_LTT,str(case),activeParameters_LTC,design_type,status=status_)
		yaml_loc,design_matrix,selected_params = PPRS_system.partial_DesignMatrix()
		#print(len(design_matrix))
		partialPRS_Object.append(PPRS_system)
		no_of_sim[case] = int(PPRS_system.no_of_sim)
		yaml_loc_dict[case] = yaml_loc
		design_matrix_dict[case] = design_matrix
		selected_params_dict_LTT[case] = selected_params
##################################################################
##  Use the unsrt data to sample the Arrhenius curves           ##
##  MUQ-SAC: Method of Uncertainty Quantification and           ##
##	Sampling of Arrhenius Curves                                ##
##################################################################
	print("\n\nSensitivity Analysis is already finished!!\n")
	SSM = simulator.SM(targets_LTT,optInputs,unsrt_data_LTC,design_matrix_dict,tag="Partial")
	
else:
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
		n_ = getTotalUnknowns(n)
		if design == "A-facto":
			sim = int(A_fact_samples)*n_
		else:
			sim = 7*n_	
		return sim
	no_of_sim = {}
	
	no_of_sim_ = getSim(len(manipulationDictLTC["activeParameters"]),design_type)
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
	
	SSM = simulator.SM(targets_LTT,optInputs,unsrt_data_LTC,design_matrix)

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
	for params in activeParameters_LTC:
		selected_params.append(1)
	
	selected_params_dict_LTT = {}
	design_matrix_dict = {}
	yaml_loc_dict = {}
	for case in case_dir:
		yaml_loc_dict[case] = yaml_loc
		design_matrix_dict[case] = design_matrix
		selected_params_dict_LTT[case] = selected_params


##############################################################
##     Doing simulations using the design matrix            ## 
##     The goal is to test the design matrix                ##
##############################################################

if "Opt" not in os.listdir():
	os.mkdir("Opt")
	os.chdir("Opt")
	optDir = os.getcwd()
	os.mkdir("Data")
	os.chdir("Data")
	os.mkdir("Simulations")
	os.mkdir("ResponseSurface")
	os.chdir("..")
else:
	os.chdir("Opt")
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
for case in case_dir_LTT:	
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
		data_sheet,failed_sim, ETA = data_management.generate_target_value_tables(FlameMaster_Execution_location, targets_LTT, case, fuel,input_=optInputs)
		#print(data_sheet)
		#raise AssertionError("!STOP")
		temp_sim_opt[str(case)] = ETA
		f = open('../Data/Simulations/sim_data_case-'+str(case)+'.lst','w').write(data_sheet)
		g = open('../Data/Simulations/failed_sim_data_case-'+str(case)+'.lst','w').write(failed_sim)
		#f.write(data_sheet)
		#f.close()
		os.chdir(optDir)


###############################################
##      Generating the response surface      ##
##                                           ##
###############################################

ResponseSurfacesLTT = {}
selected_PRS_LTT = {}
for case_index,case in enumerate(temp_sim_opt):
	yData = np.asarray(temp_sim_opt[case]).flatten()#[0:no_of_sim[case_index]]
	xData = np.asarray(design_matrix_dict[case_index])#[0:no_of_sim[case_index]]
	#print(np.shape(xData))
	#print(np.shape(yData))
	#raise AssertionError("Stop!!")
	
	xTrain,xTest,yTrain,yTest = train_test_split(xData,yData,
									random_state=104, 
                                	test_size=0.1, 
                                   	shuffle=True)
	#print(np.shape(xTest))
	#print(np.shape(yTrain))
	Response = PRS.ResponseSurface(xTrain,yTrain,case,case_index,prs_type=PRS_type,selected_params=selected_params_dict_LTT[case_index])
	Response.create_response_surface()
	Response.test(xTest,yTest)
	Response.plot(case_index)
	#Response.DoStats_Analysis() #Generates stastical analysis report
	#print(Response.case)
	ResponseSurfacesLTT[case_index] = Response
	#print(Response.selection)
	del xTrain,xTest,yTrain,yTest
#raise AssertionError("The Target class, Uncertainty class, Design Matrix and Simulations and Response surface")
os.chdir("..")

###########################################################
###########################################################
### ResponseSurfaceHTT with reaction classification (HTC)
### 
###########################################################
###########################################################
os.chdir(HTT_OPT_FOLDER)
"""
For Partial PRS (p-PRS), we need to do a second level of sensitivity analysis:	
	
	sensitivity_analysis: a dictionary that stores the sensitivity coefficients of all the active
			      parameters (selected in the first set of sensitivity analysis) reaction wise
			      as in unsrt_data dictionary
"""

if PRS_type == "Partial":
	####################################################
	## Sensitivity Analysis of the selected reactions ##
	####################################################
	
	print("\nPartial Polynomial Response Surface is choosen as a Solution Mapping Technique\n\n\tRecollecting the Sensitivity Analysis: \n\t\tKindly be Patient\n")
	if "A-facto" in design_type:
		status_ = "Pending"
		print("\n\tBrute-force Sensitivity analysis is alsready over!!\n")
		with open('sens_parameters.pkl', 'rb') as file_:
			sensitivity_analysis_HTT = pickle.load(file_)
	#print(sensitivity_analysis)
	#raise AssertionError("Stop")
	partialPRS_Object = []
	selected_params_dict_HTT = {}
	design_matrix_dict = {}
	yaml_loc_dict = {}
	no_of_sim = {}
	print("\n################################################\n###  Starting to generate Design Matrix      ###\n###  for all targets                        ###\n################################################\n\n")
	for case in case_dir_HTT:
		#PPRS_system = P_PRS.PartialPRS(sensitivity_analysis[str(case)],unsrt_data,optInputs,target_list,case,activeParameters,design_type)
		PPRS_system = P_PRS.PartialPRS(sensitivity_analysis_HTT[str(case)],unsrt_data_HTC,optInputs,targets_HTT,str(case),activeParameters_HTC,design_type,status=status_)
		yaml_loc,design_matrix,selected_params = PPRS_system.partial_DesignMatrix()
		#print(len(design_matrix))
		partialPRS_Object.append(PPRS_system)
		no_of_sim[case] = int(PPRS_system.no_of_sim)
		yaml_loc_dict[case] = yaml_loc
		design_matrix_dict[case] = design_matrix
		selected_params_dict_HTT[case] = selected_params
##################################################################
##  Use the unsrt data to sample the Arrhenius curves           ##
##  MUQ-SAC: Method of Uncertainty Quantification and           ##
##	Sampling of Arrhenius Curves                                ##
##################################################################
	print("\n\nSensitivity Analysis is already finished!!\n")
	SSM = simulator.SM(targets_HTT,optInputs,unsrt_data_HTC,design_matrix_dict,tag="Partial")
	
else:
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
		n_ = getTotalUnknowns(n)
		if design == "A-facto":
			sim = int(A_fact_samples)*n_
		else:
			sim = 7*n_	
		return sim
	no_of_sim = {}
	
	no_of_sim_ = getSim(len(manipulationDictHTC["activeParameters"]),design_type)
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
	
	SSM = simulator.SM(targets_HTT,optInputs,unsrt_data_HTC,design_matrix)

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
	for params in activeParameters_HTC:
		selected_params.append(1)
	
	selected_params_dict_HTT = {}
	design_matrix_dict = {}
	yaml_loc_dict = {}
	for case in case_dir:
		yaml_loc_dict[case] = yaml_loc
		design_matrix_dict[case] = design_matrix
		selected_params_dict_HTT[case] = selected_params


##############################################################
##     Doing simulations using the design matrix            ## 
##     The goal is to test the design matrix                ##
##############################################################

if "Opt" not in os.listdir():
	os.mkdir("Opt")
	os.chdir("Opt")
	optDir = os.getcwd()
	os.mkdir("Data")
	os.chdir("Data")
	os.mkdir("Simulations")
	os.mkdir("ResponseSurface")
	os.chdir("..")
else:
	os.chdir("Opt")
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
for case in case_dir_HTT:	
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
		data_sheet,failed_sim, ETA = data_management.generate_target_value_tables(FlameMaster_Execution_location, targets_HTT, case, fuel,input_=optInputs)
		#print(data_sheet)
		#raise AssertionError("!STOP")
		temp_sim_opt[str(case)] = ETA
		f = open('../Data/Simulations/sim_data_case-'+str(case)+'.lst','w').write(data_sheet)
		g = open('../Data/Simulations/failed_sim_data_case-'+str(case)+'.lst','w').write(failed_sim)
		#f.write(data_sheet)
		#f.close()
		os.chdir(optDir)


###############################################
##      Generating the response surface      ##
##                                           ##
###############################################

ResponseSurfacesHTT = {}
selected_PRS_HTT = {}
for case_index,case in enumerate(temp_sim_opt):
	yData = np.asarray(temp_sim_opt[case]).flatten()#[0:no_of_sim[case_index]]
	xData = np.asarray(design_matrix_dict[case_index])#[0:no_of_sim[case_index]]
	#print(np.shape(xData))
	#print(np.shape(yData))
	#raise AssertionError("Stop!!")
	
	xTrain,xTest,yTrain,yTest = train_test_split(xData,yData,
									random_state=104, 
                                	test_size=0.1, 
                                   	shuffle=True)
	#print(np.shape(xTest))
	#print(np.shape(yTrain))
	Response = PRS.ResponseSurface(xTrain,yTrain,case,case_index,prs_type=PRS_type,selected_params=selected_params_dict_HTT[case_index])
	Response.create_response_surface()
	Response.test(xTest,yTest)
	Response.plot(case_index)
	#Response.DoStats_Analysis() #Generates stastical analysis report
	#print(Response.case)
	ResponseSurfacesHTT[case_index] = Response
	#print(Response.selection)
	del xTrain,xTest,yTrain,yTest
#raise AssertionError("The Target class, Uncertainty class, Design Matrix and Simulations and Response surface")
os.chdir(opt_dir)


##################################################
##        Optimization Procedure                ##
##   Inputs: Traget list and Response Surfaces  ## 
##################################################
if "solution_zeta.save" not in os.listdir():

	opt, opt_zeta,posterior_cov = Optimizer(targets_LTT,targets_HTT).run_multi_stage_parallel_optimization_with_selected_PRS(unsrt_data_LTC,unsrt_data_HTC,ResponseSurfacesLTT,ResponseSurfacesHTT,optInputs)
	
	opt_zeta_LTC = opt_zeta[:len(activeParameters_LTC)]
	opt_zeta_HTC = opt_zeta[len(activeParameters_LTC):]
	string_save_LTC = ""
	for i, j in enumerate(activeParameters_LTC):
		print("\n{}=\t{}".format(activeParameters_LTC[i], opt_zeta_LTC[i]))
		string_save_LTC+="{}=\t{}\n".format(activeParameters_LTC[i], opt_zeta_LTC[i])
	save = open("solution_zeta_LTC.save","w").write(string_save_LTC)
	
	
	string_save_HTC = ""
	for i, j in enumerate(activeParameters_HTC):
		print("\n{}=\t{}".format(activeParameters_HTC[i], opt_zeta_HTC[i]))
		string_save_HTC+="{}=\t{}\n".format(activeParameters_HTC[i], opt_zeta_HTC[i])
	save = open("solution_zeta_HTC.save","w").write(string_save_HTC)
	
	originalMech = Parser(mech_file_location).mech
	copy_of_mech = deepcopy(originalMech)#.deepcopy()
	new_mechanism_LTC,a = Manipulator(copy_of_mech,unsrt_data_LTC,opt_zeta_LTC).doPerturbation()
	new_mechanism_HTC,a = Manipulator(new_mechanism_LTC,unsrt_data_HTC,opt_zeta_HTC).doPerturbation()
	
	string = yaml.safe_dump(new_mechanism_HTC,default_flow_style=False)
	f = open("new_mech.yaml","w").write(string)
else:
	
	save_LTC = open("solution_zeta_LTC.save","r").readlines()
	save_HTC = open("solution_zeta_HTC.save","r").readlines()
	opt_zeta_LTC = []
	for i in save_LTC:
		opt_zeta_LTC.append(float(i.split("=")[1].strip()))
	opt_zeta_HTC = []
	for i in save_HTC:
		opt_zeta_HTC.append(float(i.split("=")[1].strip()))
	
	originalMech = Parser(mech_file_location).mech
	copy_of_mech = deepcopy(originalMech)#.deepcopy()
	new_mechanism_LTC,a = Manipulator(copy_of_mech,unsrt_data_LTC,opt_zeta_LTC).doPerturbation()
	new_mechanism_HTC,a = Manipulator(new_mechanism_LTC,unsrt_data_HTC,opt_zeta_HTC).doPerturbation()
	
	string = yaml.safe_dump(new_mechanism_HTC,default_flow_style=False)
	f = open("new_mech.yaml","w").write(string)



raise AssertionError("The Target class, Uncertainty class, Design Matrix and Simulations and Response surface")





