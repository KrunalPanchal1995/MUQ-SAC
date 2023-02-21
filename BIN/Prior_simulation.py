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
import yaml
import pandas as pd
import sensitivity
import test,Nominal,Optimal
import generate
import solution
#program specific modules

import make_input_file
import FlameMaster_in_parallel
import combustion_dataset_class
import combustion_variable_class
from combustion_optimization_class import OptimizationTool
import combustion_target_class
import data_management
import simulation_manager
import Uncertainty as uncertainty
import MechanismManipulator
import plotter
import Input_file_reader
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
#print(targetLines)
addendum = yaml.safe_load(open(locations[add],'r').read())
##!!!  MAKE A LIST OF TARGET CLASS CONTAINING EACH TARGET AS A CASE
def filter_list(List):
	temp = []
	for i in List:
		if i != "":
			temp.append(i)
	return temp
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
##########Obtain uncertainty data from user input file!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

UncertDataSet = uncertainty.uncertaintyData(locations,binLoc);
unsrt_data,rxnUnsrt_data,focUnsrt_data,tbdUnsrt_data,thermoUnsrt_data,transportUnsrt_data, reaction_index,fallOffCurve_index, thirdBody_index, thermo_index, transport_index = UncertDataSet.extract_uncertainty();
print("Uncertainty analysis finished")
BranchingRxns = []
#pDeptRxns=[]
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

import multiprocessing

print("Number of cpu : ", multiprocessing.cpu_count())
### DO prior simulations #####

#testDir,manipulation_dict,sim_dict,response_surface_dict = Nominal.simulations(optInputs,iFile,case_dir,rps_order,activeParameters, reaction_index,fallOffCurve_index, thirdBody_index, thermo_index, transport_index, mech_file_location, fileType, rxnUnsrt_data, focUnsrt_data,tbdUnsrt_data, thermoUnsrt_data, transportUnsrt_data, target_list, fuel, global_reaction, thermo_file_location, trans_file_location,startProfile_location,design_type,parallel_threads,file_specific_input,rIndex,unsrt_data,manipulation_dict,sim_dict,response_surface_dict,selectedParams,activeIndexDict)

manipulation_dict = {}
sim_dict = {}
response_surface_dict = {}
activeIndexDict = {}
activeReactions = {}

originalDir,manipulation_dict,sim_dict,response_surface_dict = Nominal.simulations(optInputs,iFile,case_dir,rps_order,activeParameters, reaction_index,fallOffCurve_index, thirdBody_index, thermo_index, transport_index, mech_file_location, fileType, rxnUnsrt_data, focUnsrt_data,tbdUnsrt_data, thermoUnsrt_data, transportUnsrt_data, target_list, fuel, global_reaction, thermo_file_location, trans_file_location,startProfile_location,design_type,parallel_threads,file_specific_input,rIndex,unsrt_data,manipulation_dict,sim_dict,response_surface_dict,selectedParams,activeIndexDict,activeReactions)

### Plot prior results ######

#############################################################
### 	      Building maga database 	   					#
###			for parameters and targets						#
#############################################################
#print(manipulation_dict)

parametric_space = {}
for i in activeParameters:
	temp_database = {}
	temp_database["unsrtDatabase"] = unsrtDatabase[str(i)]
	temp_database["sensitivityManipulations"] = None
	temp_database["optManipulations"] = None
	#temp_database["originalManipulations"] = manipulation_dict["Original"][str(i)]
	temp_database["optimizedParameter"] = None
	temp_database["posteriorCovariance"] = None
	parametric_space[str(i)] = temp_database


DataSet_Tig = []
DataSet_Fls = []
DataSet_Flw = []
DataSet_Flf = []
target_space = {}
for case in case_dir:
	if target_list[case].target == "Tig":
		DataSet_Tig.append(target_list[case].d_set)
	if target_list[case].target == "Fls":
		DataSet_Fls.append(target_list[case].d_set)
	if target_list[case].target == "Flw":
		DataSet_Flw.append(target_list[case].d_set)
	if target_list[case].target == "Flf":
		DataSet_Flf.append(target_list[case].d_set)
target_space["Tig"] = list(OrderedDict.fromkeys(DataSet_Tig))
target_space["Fls"] = list(OrderedDict.fromkeys(DataSet_Fls))
target_space["Flw"] = list(OrderedDict.fromkeys(DataSet_Flw))
target_space["Flf"] = list(OrderedDict.fromkeys(DataSet_Flf))

dataSet = {}
#populating dataSet

for d_types in target_space:
	points = {}
	for set_index in target_space[d_types]:
		set_list = []
		for case in case_dir:
			set_data = {}
			if target_list[case].d_set == set_index:
				set_data["Target_database"] = target_list[case]
				set_data["SA_simulations"] = None
				set_data["Original_simulations"] = sim_dict["Original"][str(case)]
				if "PRS" in optInputs["Type"]["optimization_type"]:
					#set_data["Opt_simulations"] = sim_dict["Opt"][str(case)]
					#set_data["SA_PRS"] = response_surface_dict["sa"][str(case)]
					#set_data["Opt_PRS"] = response_surface_dict["Opt"][str(case)]
					set_data["optimization_type"]= "PRS"
				else:
					set_data["optimization_type"] = "Direct"
				set_data["Opt_simulations"] = None
				set_data["SA_PRS"] = None
				set_data["Opt_PRS"] = None
				set_data["Prior_model_unsrt"] = None
				set_data["Prior_model_response"] = None
				set_data[str(case)] = None
				set_list.append(set_data)
			else:
				continue
		points[str(set_index)] = set_list
	dataSet[str(d_types)]= points
#print(dataSet)
#parameter_list = []
#for i in parametric_space:
#	p = combustion_variable_class.combustion_variable(i,parametric_space[i])
#	parameter_list.append(p)

target_sets = []
for targetClass in dataSet:
	for index,sets in enumerate(dataSet[targetClass]):
		tar = combustion_dataset_class.combustion_dataset(sets,dataSet[targetClass][sets],[])
		target_sets.append(tar)

#plotter.plot(activeParameters,parameter_list,case_dir,target_sets)
#plotter.plot_errors(target_list,opt.x)
#Theta = []
#for i in reaction_index:
#	T1,T2 = unsrt_data[i].getTempRange()
#	T = np.linspace(T1,T2,1000)
#	theta = []
#	for i in T:
#		Th1 = 1
#		Th2 = np.log(i)
#		Th3 = -1/i
#		theta.append([Th1,Th2,Th3])
#	theta = np.asarray(theta)
#	Theta.append(theta)
#plotter.plot_errors(target_list, opt.x)
#plotter.plot_vector(reaction_index, opt.x)
print("generating prior files..... \n\n")
#plotter.plot_graphs(target_list,init_guess,posterior_cholesky,opt.x)
#plotter.plot_rxnUnsrt(mech_file_location,unsrt_location,thermo_file_location,trans_file_location,target_list,sim_type,Theta,Posterior_cov,opt.x,"True")
#IFR_unopt = Input_file_reader.MechParsing(mech_file_location)
#IFR_opt = Input_file_reader.MechParsing("./mechanism_opt.mech")

#for index,rxn_list in enumerate(BranchingRxns):
#	unOptRatios = IFR_unopt.getRatioData(rxn_list,"unOpt")
	
#	OptRatios = IFR_opt.getRatioData(rxn_list,"Opt")
	
#	fig = plt.figure()

#	T = np.linspace(500,2500,len(np.asarray(unOptRatios)[:,0]))
#	for i in range(len(unOptRatios[0])):
#		plt.plot(T,np.asarray(unOptRatios)[:,i],"--")
	
	#	plt.plot(T,np.asarray(OptRatios)[:,i],"-")
	#	plt.xlabel("Temperature (K)")
	#	plt.ylabel("Branching ratios")
	#	plt.title("Branching ratios for {} \nbefore and after optimization".format(rxn_list))
	#	plt.savefig("./BranchRatio/"+str(rxn_list)+".png")
#plotter.print_cov(Posterior_cov)
#plotter.print_cov(pcov_rxn,"reaction")
#plotter.print_cov(pcov_foc,"FoC")
#plotter.print_cov(pcov_hcp,"Hcp")
#plotter.print_cov(pcov_m,"collision_efficiency")
#plotter.print_cov(pcov_enthalpy,"enthalpy")
#plotter.print_cov(pcov_enthalpy,"entropy")

for ind,sets in enumerate(target_sets):
	index = sets.set_id[0]
	#print(index)
	add = sets.addendum[0]
	plotter.plot_optResults(index,add,sets,opt_tag="off")

print(">>>>>> Program Exiting <<<<<<")

