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
import partial_PRS_system as P_PRS
#import yaml
#import ruamel.yaml as yaml
try:
    import ruamel_yaml as yaml
except ImportError:
    from ruamel import yaml
import pandas as pd

sys.path.append('/parallel_yaml_writer.so')
import parallel_yaml_writer
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
PRS_type = stats_["PRS_type"]
#######################READ TARGET FILE ###################
print("Parallel threads are {}".format(parallel_threads))
targetLines = open(locations[targets],'r').readlines()
addendum = yaml.safe_load(open(locations[add],'r').read())

print(design_type)

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
print("optimization targets identified\n")
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
	#unsrt_data,rxnUnsrt_data,plogUnsrt_data,plog_interpolated_data,focUnsrt_data,tbdUnsrt_data,thermoUnsrt_data,transportUnsrt_data, reaction_index,plog_boundary_index,plog_index,fallOffCurve_index, thirdBody_index, thermo_index, transport_index = UncertDataSet.extract_uncertainty();
	print("Uncertainty analysis finished")

else:
	# Load the object from the file
	with open('unsrt.pkl', 'rb') as file_:
		unsrt_data = pickle.load(file_)
	print("Uncertainty analysis already finished")

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
print("\nFollowing list is the choosen reactions\n")
print(manipulationDict["activeParameters"])

if PRS_type == "Partial":
	####################################################
	## Sensitivity Analysis of the selected reactions ##
	####################################################
	print("\nPartial Polynomial Response Surface is choosen as a Solution Mapping Technique\n\n\tStarting the Sensitivity Analysis: \n\t\tKindly be Patient\n")
	#if "A-facto" in design_type:
	#if "SA" not in os.listdir():
	# Arguments to pass to script2.py
	args = [sys.argv[1], 'rxn_list_a_fact.csv','&>SA.out']

	# Doing A factor sensitivity Analysis
	result = subprocess.run(['python3.9', binLoc["SA_tool"]] + args, capture_output=True, text=True)

	# Printing the output of script2.py
	#print("\nOutput of sens.py:\n")
	#print(result.stdout)
	print("\nSensitivity Analysis for A-factor is Done!!")
	# Printing the errors of script2.py, if any
	#if result.stderr:
	#	print("Errors:")
	#	print(result.stderr)
			
	#	raise AssertionError("Sensitivity Analysis Done!!")
	with open('sens_parameters.pkl', 'rb') as file_:
		sensitivity_analysis = pickle.load(file_)
	
	#else:
	#if "SA_3P" not in os.listdir():

	args = [sys.argv[1], 'rxn_list.csv','&>SA_3p.out']

	result = subprocess.run(['python3.9', binLoc["SA_tool_3p"]] + args, capture_output=True, text=True)

	# Printing the output of script2.py
	#print("\nOutput of sens.py:\n")
	#print(result.stdout)
	print("\nSensitivity Analysis for 3-param is Done!!")
	# Printing the errors of script2.py, if any
	#if result.stderr:
	#	print("Errors:")
	#	print(result.stderr)
	
	raise AssertionError("Sensitivity Analysis Done!!")
	with open('sens_3p_parameters.pkl', 'rb') as file_:
		sensitivity_analysis_3p = pickle.load(file_)
	
	partialPRS_Object = []
	selected_params_dict = {}
	design_matrix_dict = {}
	yaml_loc_dict = {}
	no_of_sim = {}
	for case in case_dir:
		PPRS_system = P_PRS.PartialPRS(sensitivity_analysis_3p[str(case)],unsrt_data,optInputs,target_list,case,activeParameters,design_type)
		yaml_loc,design_matrix,selected_params = PPRS_system.partial_DesignMatrix()
		partialPRS_Object.append(PPRS_system)
		no_of_sim[case] = int(PPRS_system.no_of_sim)
		yaml_loc_dict[case] = yaml_loc
		design_matrix_dict[case] = design_matrix
		selected_params_dict[case] = selected_params
##################################################################
##  Use the unsrt data to sample the Arrhenius curves           ##
##  MUQ-SAC: Method of Uncertainty Quantification and           ##
##	Sampling of Arrhenius Curves                                ##
##################################################################
	print("\n\nSensitivity Analysis is already finished!!\n")
	SSM = simulator.SM(target_list,optInputs,unsrt_data,design_matrix_dict)
	
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
	for case in target_list:
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
	for case in target_list:
		yaml_loc_dict[case] = yaml_loc
		design_matrix_dict[case] = design_matrix
		selected_params_dict[case] = selected_params
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

###############################################
##      Generating the response surface      ##
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
	
	xTrain,xTest,yTrain,yTest = train_test_split(xData,yData,
									random_state=104, 
                                	test_size=0.1, 
                                   	shuffle=True)
	#print(np.shape(xTest))
	#print(np.shape(yTrain))
	Response = PRS.ResponseSurface(xTrain,yTrain,case,case_index,prs_type=PRS_type,selected_params=selected_params_dict[case_index])
	Response.create_response_surface()
	Response.test(xTest,yTest)
	Response.plot(case_index)
	#print(Response.case)
	ResponseSurfaces[case_index] = Response
	#print(Response.selection)
	del xTrain,xTest,yTrain,yTest
#raise AssertionError("The Target class, Uncertainty class, Design Matrix and Simulations and Response surface")
os.chdir("..")
##################################################
##        Optimization Procedure                ##
##   Inputs: Traget list and Response Surfaces  ## 
##################################################

opt, opt_zeta,posterior_cov = Optimizer(target_list).run_optimization_with_selected_PRS(unsrt_data,ResponseSurfaces,optInputs)

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
string = yaml.dump(new_mechanism,default_flow_style=False)
f = open("new_mech.yaml","w").write(string)



raise AssertionError("The Target class, Uncertainty class, Design Matrix and Simulations and Response surface")









