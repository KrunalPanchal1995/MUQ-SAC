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
#import yaml
#import ruamel.yaml as yaml
try:
    import ruamel_yaml as yaml
except ImportError:
    from ruamel import yaml
import pandas as pd
import sensitivity
import test,Nominal,Optimal
import generate
import solution
from MechManipulator2_0 import Manipulator
#program specific modules
from copy import deepcopy
from MechanismParser import Parser
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
unsrt_data,rxnUnsrt_data,plogUnsrt_data,plog_interpolated_data,focUnsrt_data,tbdUnsrt_data,thermoUnsrt_data,transportUnsrt_data, reaction_index,plog_boundary_index,plog_index,fallOffCurve_index, thirdBody_index, thermo_index, transport_index = UncertDataSet.extract_uncertainty();
print("Uncertainty analysis finished")
string_Rxn = ""
for i in unsrt_data:
	string_Rxn += f"{i}\t{unsrt_data[i].index}\n"
file_rxn = open("RXN.csv","w").write(string_Rxn)

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

raise AssertionError("First stop!!")

for i in rxnUnsrt_data:
	print(i)

"""
Selection of the parameters
"""
print(activeParameters)
selectedParams = deepcopy(activeParameters)#.deepcopy()


#beta = np.ones(3*len(selectedParams))
#print(selectedParams)
#print(beta)
#selection = np.ones(len(selectedParams))
###
# Manipulation test ###
###
#mech = open(mech_file_location,"r")
#mechanism_test = yaml.safe_load(mech)
#print(mechanism_test["reactions"][0])

#test1 = Manipulator(mechanism_test,unsrt_data,beta,rxn_list)
#new_mechanism = test1.doPerturbation()


#string = yaml.dump(new_mechanism,default_flow_style=False)#,explicit_start=True,width=70)
#f = open("new_mech.yaml","w").write(string)
#raise AssertionError("Plz wait the code is incomplete yet!!")
unsrtDatabase = {}
reac = []
for index,rxn in enumerate(unsrt_data):
	reac.append(rxn)
	unsrtDatabase[rxn] = unsrt_data[rxn].getDtList()
	
"""
BranchingRxns = []
#pDeptRxns=[]

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
"""
"""
selectedParams = []
for rxn in rxnUnsrt_data:
	#print(rxn)
	selectedParams.extend(rxnUnsrt_data[rxn].activeParameters)
"""
import multiprocessing

print("Number of cpu : ", multiprocessing.cpu_count())


manipulation_dict = {}
sim_dict = {}
response_surface_dict = {}
########### DO SENSITIVITY ANALYSIS ######################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#print(activeParameters)
"""
sensDir,manipulation_dict,sim_dict,response_surface_dict,sensitivityDict,activeIndexDict,optimization_bounds,optimization_init_guess,activeReactions = sensitivity.analysis(optInputs,
										iFile,case_dir,
										rps_order,
										activeParameters,
										reaction_index,
										fallOffCurve_index,
										thirdBody_index, 
										thermo_index, 
										transport_index, 
										mech_file_location, 
										fileType, 
										rxnUnsrt_data, 
										focUnsrt_data,
										tbdUnsrt_data, 
										thermoUnsrt_data, 
										transportUnsrt_data, 
										target_list, 
										fuel, 
										global_reaction, 
										thermo_file_location,
										trans_file_location,
										startProfile_location,
										design_type,
										parallel_threads,
										file_specific_input,
										rIndex,
										unsrt_data,
										manipulation_dict,
										sim_dict,
										response_surface_dict,selectedParams,manipulationDict)



#print(activeIndexDict)
#print("\n===activeReactions===\n")
#print(activeReactions)
#raise AssertionError("Stop!!")

"""

sensDir = {}
activeReactions = {}
activeIndexDict = {}

if "PRS" in optInputs["Type"]["optimization_type"]: 
	testDir,manipulation_dict,sim_dict,response_surface_dict = test.simulations(optInputs,iFile,case_dir,rps_order,activeParameters, reaction_index,fallOffCurve_index, thirdBody_index, thermo_index, transport_index, mech_file_location, fileType, rxnUnsrt_data, focUnsrt_data,tbdUnsrt_data, thermoUnsrt_data, transportUnsrt_data, target_list, fuel, global_reaction, thermo_file_location,trans_file_location,startProfile_location,design_type,parallel_threads,file_specific_input,rIndex,unsrt_data,manipulation_dict,sim_dict,response_surface_dict,selectedParams,manipulationDict,activeIndexDict,activeReactions)
	########### DO Original ANALYSIS ######################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	originalDir,manipulation_dict,sim_dict,response_surface_dict = Nominal.simulations(optInputs,iFile,case_dir,rps_order,activeParameters, reaction_index,fallOffCurve_index, thirdBody_index, thermo_index, transport_index, mech_file_location, fileType, rxnUnsrt_data, focUnsrt_data,tbdUnsrt_data, thermoUnsrt_data, transportUnsrt_data, target_list, fuel, global_reaction, thermo_file_location, trans_file_location,startProfile_location,design_type,parallel_threads,file_specific_input,rIndex,unsrt_data,manipulation_dict,sim_dict,response_surface_dict,selectedParams,manipulationDict,activeIndexDict,activeReactions)
	
	
	########### DO Original ANALYSIS ######################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	optDir,manipulation_dict,sim_dict,response_surface_dict = Optimal.simulations(optInputs,iFile,case_dir,rps_order,activeParameters, reaction_index,fallOffCurve_index, thirdBody_index, thermo_index, transport_index, mech_file_location, fileType, rxnUnsrt_data, focUnsrt_data,tbdUnsrt_data, thermoUnsrt_data, transportUnsrt_data, target_list, fuel, global_reaction, thermo_file_location,trans_file_location,startProfile_location,design_type,parallel_threads,file_specific_input,rIndex,unsrt_data,manipulation_dict,sim_dict,response_surface_dict,selectedParams,manipulationDict,activeIndexDict,activeReactions)

	############################################       
	#	RESPONSE SURFACE DEVELOPMENT       #
	#	and Testing	 		    #
	############################################
	selectedPRS,bound,init_guess = generate.response_surface(sensDir,testDir,optDir,optInputs,case_dir,response_surface_dict,target_list,activeParameters,unsrt_data,stats_,selectedParams,activeIndexDict)
	####################################################################################
	####									############
	####  Build the objective function and perform scipy.optimize routine.!!!!!!!!!!!!!!
	####									############
	####################################################################################
	sim_type = ""
	simulator = simulation_manager.SM(optInputs,iFile,case_dir,rps_order,activeParameters, reaction_index,fallOffCurve_index,thirdBody_index,thermo_index,transport_index, mech_file_location,fileType,unsrt_data, rxnUnsrt_data, focUnsrt_data, tbdUnsrt_data, thermoUnsrt_data, transportUnsrt_data, target_list, fuel, global_reaction, thermo_file_location, trans_file_location,startProfile_location,design_type,sim_type,parallel_threads,selectedParams,manipulationDict,activeReactions,extra_arg = activeIndexDict)
	#selectedPRS = {}
	#for i in case_dir:
	#	selectedPRS[i] = 1		
	#print(selectedPRS)
	var = len(selectedParams)
	#print(var)
	mean = np.zeros(var)
	init_variance = np.ones(var)
	#init_guess = np.random.normal(mean,init_variance,var)
	init_guess = np.zeros(var)
	init_covariance = 0.25*np.eye(var)
	inv_cov_var = np.linalg.inv(init_covariance)
	#print("\n\nAlgorithm used for optimization : {}\n\n".format(opt.get_algorithm_name()))
	t = time.time()
	#opt = OptimizationTool(target_list).run_optimization(simulator,selectedPRS,init_guess,initial_covariance=inv_cov_var)
	#bnds = tuple((-i,i) for i in bound)
	#bound = 1000*np.ones(var)
	#bound = optimization_bounds
	#init_guess = optimization_init_guess
	
	bnds = tuple((-1,1) for _ in bound)
	
	
	#opt,opt_zeta = OptimizationTool(target_list).run_optimization_with_selected_PRS(simulator,
	#									selectedPRS,
	#									method=optInputs["Type"]["Algorithm_type"],
	#									algorithm=optInputs["Type"]["Algorithm"],
	#									initial_guess=init_guess,
	#									bounds=bnds,
	#									bounds_array=bound,
	#									initial_covariance=inv_cov_var)
	
	opt,opt_zeta = OptimizationTool(target_list).run_optimization_with_selected_PRS(simulator,
										selectedPRS,
										method=optInputs["Type"]["Algorithm_type"],
										algorithm=optInputs["Type"]["Algorithm"],
										initial_guess=init_guess,
										bounds=bnds,
										bounds_array=bound,
										initial_covariance=inv_cov_var)
	
	
	#print(len(opt_zeta))
	string_save = ""
	for i, j in enumerate(selectedParams):
		print("\n{}=\t{}".format(selectedParams[i], opt_zeta[i]))
		string_save+="{}=\t{}\n".format(selectedParams[i], opt_zeta[i])
	save = open("solution_zeta.save","w").write(string_save)
	
	J = []
	W = []
	E_y = []
	E_d = []
	for i,case in enumerate(target_list):
		if selectedPRS[i] == 1:
			#J.append(case.getJacobian(opt))
			W.append(case.d_weight)
			E_y.append(case.std_dvtn**2)
			#E_d.append(float(case.calculated_target_value(opt)-case.observed)**2)
			
	#W = np.diag(W)
	#E_y = np.diag(E_y)
	#E_d = np.diag(E_d)
	#J = np.asarray(J)
	#J = np.matrix(J)
	#iE_y = np.diag(1/E_y)
	#S = np.add(E_y,E_d)
	#Matrix_B = np.linalg.inv(J.T.dot(W.dot(iE_y.dot(J))))
	#Matrix_A = J.T.dot(W.dot(iE_y))
	#Matrix_C = Matrix_B.dot(Matrix_A.T)
	#posterior_cov = Matrix_C.dot(S.dot(Matrix_C.T))
	#print(np.shape(posterior_cov))
	#print(np.shape(init_covariance))
	dt = int(time.time() - t)
	hours = dt/3600
	minutes = (dt%3600)/60
	seconds = dt%60
	cov_opt = init_covariance
	print("Time for performing Optimization: {h} hours,  {m} minutes, {s} seconds\n................................................ ".format(h = hours, m = minutes, s =seconds))
	print(">>>>>>>>>>>>>>>\n\nOptimized Vectors\n\n>>>>>>>>>>>>>>>>>")


elif "Direct" in optInputs["Type"]["optimization_type"]:
	manipulation_dict["Opt"] = {}
	#manipulation_dict["Original"] = {}
	originalDir,manipulation_dict,sim_dict,response_surface_dict = Nominal.simulations(optInputs,iFile,case_dir,rps_order,activeParameters, reaction_index,fallOffCurve_index, thirdBody_index, thermo_index, transport_index, mech_file_location, fileType, rxnUnsrt_data, focUnsrt_data,tbdUnsrt_data, thermoUnsrt_data, transportUnsrt_data, target_list, fuel, global_reaction, thermo_file_location, trans_file_location,startProfile_location,design_type,parallel_threads,file_specific_input,rIndex,unsrt_data,manipulation_dict,sim_dict,response_surface_dict,selectedParams,manipulationDict,activeIndexDict)
	"""
	This script is for direct optimization
	without the use of Polynomial Response Surface
	
	"""
	sim_type = ""
	simulator = simulation_manager.SM(optInputs,iFile,case_dir,rps_order,activeParameters, reaction_index,fallOffCurve_index,thirdBody_index,thermo_index,transport_index, mech_file_location,fileType, rxnUnsrt_data, focUnsrt_data, tbdUnsrt_data, thermoUnsrt_data, transportUnsrt_data, target_list, fuel, global_reaction, thermo_file_location, trans_file_location,startProfile_location,design_type,sim_type,parallel_threads,selectedParams,manipulationDict,extra_arg = activeIndexDict)

	"""
	This is a test script:
		input value is random
		need to check if all the files are generated as desired
	"""
	var = len(selectedParams)
	mean = np.zeros(var)
	init_variance = np.ones(var)
	init_guess = np.random.normal(mean,init_variance,var)
	
	init_covariance = 0.25*np.eye(var)
	inv_cov_var = np.linalg.inv(init_covariance)
	#print("\n\nAlgorithm used for optimization : {}\n\n".format(opt.get_algorithm_name()))
	#init_guess = np.zeros(var)
	t = time.time()
	selectedPRS = {}
	for i in case_dir:
		selectedPRS[i] = 0 # as no PRS are selected
	bnds = tuple((-1,1) for i in range(var))
	opt,opt_zeta = OptimizationTool(target_list).run_optimization_with_selected_PRS(simulator,
										selectedPRS,
										method=optInputs["Type"]["Algorithm_type"],
										initial_guess=init_guess,
										bounds=bnds,
										initial_covariance=inv_cov_var)
	string_save = ""
	for i, j in enumerate(selectedParams):
		print("\n{}=\t{}".format(selectedParams[i], opt[i]))
		string_save+="{}=\t{}\n".format(selectedParams[i], opt[i])
	save = open("solution.save","w").write(string_save)
	cov_opt = init_covariance
	dt = int(time.time() - t)
	hours = dt/3600
	minutes = (dt%3600)/60
	seconds = dt%60
	print("Time for performing Optimization: {h} hours,  {m} minutes, {s} seconds\n................................................ ".format(h = hours, m = minutes, s =seconds))
	print(">>>>>>>>>>>>>>>\n\nOptimized Vectors\n\n>>>>>>>>>>>>>>>>>")
	

#opt,cov_opt = OptimizationTool(target_list).run_optimization(init_guess,initial_covariance=inv_cov_var)

string_save = ""
for i, j in enumerate(selectedParams):
	print("\n{}=\t{}".format(selectedParams[i], opt[i]))
	string_save+="{}=\t{}\n".format(selectedParams[i], opt[i])
save = open("solution.save","w").write(string_save)
#for i, j in enumerate(reaction_index):
#	print("\n{}\t{}".format(j, opt[i]))
#print("\n\n\nInitial objective function value : {}\n".format(obj_function(init_guess)))
#print("Optimized objective function value : {}\n".format(obj_function(opt)))
###########################################
#    Test response surface again	   #
#					   #
###########################################
"""
for k,case in enumerate(target_list):
	xData_testing, yData_testing = data_management.getTestingData(sensDir,case_dir[target_list.index(case)])
	
	for index,data in enumerate(yData_testing):
		c.append(np.log(case.std_dvtn/1000))
		actual_error_testing.append(((case.calculated_target_value(np.asarray(xData_testing[index]))-data)/data)*100)
		temp_relative = abs((case.calculated_target_value(np.asarray(xData_testing[index]))-data)/data)*100
		actual_value = case.calculated_target_value(np.asarray(xData_testing[index]))
		Sim_value_testing.append(actual_value)
		error_testing.append(abs((case.calculated_target_value(np.asarray(xData_testing[index]))-data)))
		string +="{},{},{}\n".format(data,actual_value,temp_relative)

"""
####################################
###		Posterior covariance	####
####################################

model_priorResponse = {}
model_priorUnsrt = {}


model_posteriorResponse = {}
model_posteriorUnsrt = {}

for case in target_list:
	model_priorResponse[case] = case.evaluate(np.zeros(var),stats_[order])#case.calculated_target_value(np.zeros(var))
	
	#model_priorResponse[case], model_priorUnsrt[case] = case.evaluateResponse(np.zeros(var),cov_x=np.eye(var))
	#model_posteriorResponse[case],model_posteriorUnsrt[case] = case.evaluateResponse(opt,cov_x=cov_opt)
	model_posteriorResponse[case] = case.evaluate(opt,stats_[order])
###################################################      
########                        	 ###########              
########   	 Post Data Analysis      ###########           
########                        	 ###########              
####################################################       

	
####Build the optimized mechanism!!!!!!!!!!!!!!!!!!!         

originalMech = Parser(mech_file_location).mech
copy_of_mech = deepcopy(originalMech)#.deepcopy()
new_mechanism,a = Manipulator(copy_of_mech,unsrt_data,opt,reaction_index,np.ones(len(opt))).doPerturbation()

#new_mechanism,a,b,c = self.MechManipulator.GeneratePerturbedMechanism(self.target_list[case],self.beta_[i],np.ones(len(selectedParams)),reactionList,self.simulation,"False",extra_arg = self.activeIndexDict)
string = yaml.dump(new_mechanism,default_flow_style=False)
f = open("new_mech.yaml","w").write(string)
#raise AssertionError("Stop!!")
#MechManipulator = MechanismManipulator.MechanismManipulator(mech_file_location,fileType,thermo_file_location,trans_file_location,reaction_index,fallOffCurve_index,thirdBody_index,thermo_index,transport_index,rxnUnsrt_data, focUnsrt_data,tbdUnsrt_data, thermoUnsrt_data,transportUnsrt_data,selectedParams,activeIndexDict,design_type=design_type)
#a,b,c,param_opt_dict = MechManipulator.GeneratePerturbedMechanism(target_list,opt_zeta,np.ones(len(selectedParams)),reaction_index,"Opt","True")

#originalDir,manipulation_dict,sim_dict,response_surface_dict = Nominal.simulations(optInputs,iFile,case_dir,rps_order,activeParameters, reaction_index,fallOffCurve_index, thirdBody_index, thermo_index, transport_index, locations["opt_mech"], fileType, rxnUnsrt_data, focUnsrt_data,tbdUnsrt_data, thermoUnsrt_data, transportUnsrt_data, target_list, fuel, global_reaction, locations["opt_thermo"], locations["opt_trans"],startProfile_location,design_type,parallel_threads,file_specific_input,rIndex,unsrt_data,manipulation_dict,sim_dict,response_surface_dict,selectedParams,activeIndexDict)



#### Posterior covariance matrix

J = []
W = []
Ey = []
dY = []
'''
for case in target_list:
	#print("The index for case-{}\n".format(target_list.index(case)))
	if case.target.strip() == "Tig":
		dY.append(case.evaluate(opt,2) - (case.observed+case.std_dvtn))
		J.append(case.Jacobian(opt,2))
		W.append(case.d_weight)
		Ey.append((np.log(case.std_dvtn))**2)
	elif case.target.strip() == "Fls":
		dY.append(case.evaluate(opt,2) - case.observed)
		J.append(case.Jacobian(opt,2))
		W.append(case.d_weight)
		Ey.append((np.log(case.std_dvtn))**2)
	elif case.target.strip() == "Flw":
		dY.append(case.evaluate(opt,2) - case.observed)
		J.append(case.Jacobian(opt,2))
		W.append(case.d_weight)
		Ey.append((np.log(case.std_dvtn))**2)

#######################################################################
#####								#######
#####Posterior Cov Method of David Sheen and Turanyi, Nagy et at#######
#####								#######
#######################################################################
Ey = np.matrix(np.diag(np.asarray(Ey)))
#3Ed = np.matrix(np.outer(np.asarray(dY),np.asarray(dY)))
#S = Ed+np.matrix(Ey)
J = np.matrix(J)
#print(np.shape(J))
#print(np.shape(J.T))
#W = np.diag(np.asarray(W))
Ey_ = np.linalg.inv(np.matrix(Ey))
#print("Ed is {}\n".format(Ed))
#print("Ey_ is {}\n".format(Ey_))
#print("J is {}\n".format(J))
#print("W is {}\n".format(W))
#print("S is {}\n".format(S))
#print("Ey is {}\n".format(Ey))
#A_o = np.dot(J.T,np.dot(W,Ey_))
#print("A_o shape is {}\n".format(np.shape(A_o)))
#print("J shape is {}\n".format(np.shape(J)))
#A_j = np.dot(A_o,J)
#print("A_j shape is {}\n".format(np.shape(A_j)))
#print("A_o is {}\n".format(A_o))
#A_o_ = np.linalg.inv(A_j)
#B_o = np.dot(A_o_,A_o)
#cov = np.linalg.inv(np.dot(B_o,np.dot(S,B_o.T)))
#posterior_cholesky = np.linalg.cholesky(cov)
#cov_inverse = np.dot(J.T,np.dot(Ey_,J)) + np.identity(len(opt),dtype=float)
#cov = np.linalg.inv(cov_inverse)
#posterior_cholesky = np.linalg.cholesky(cov)
#print("shape of cov is  {}\n".format(np.shape(posterior_cov)))
#posterior_cholesky = np.linalg.cholesky(posterior_cov)
#print("Cov has shape {}\n".format(shape(cov)))
#print("B_o is {}\n".format(B_o))x

#print("cov is {}\n".format(cov))
#print("cholesky is {}\n".format(posterior_cholesky.T))
'''
unsrtPosterior = {}
count = 0
#dict_active_params

#for i,rxn in enumerate(reaction_index):
#	m = dict_active_params[rxn]
#	pcov_rxn = plotter.extract_block_diag(cov_opt[count:count+m,count:count+m],M=m)
#	unsrtPosterior[str(rxn)] = pcov_rxn
#	count+=m

#pcov_foc = plotter.extract_block_diag(posterior_cholesky[count:count+4*len(fallOffCurve_index),count:count+4*len(fallOffCurve_index)],M=1)
#for i,foc in enumerate(fallOffCurve_index):
#	unsrtPosterior[str(foc)] = list(pcov_foc)[i]

#count=count+len(fallOffCurve_index)
#pcov_m = plotter.extract_block_diag(cov_opt[count_Arrhenius_params:,count_Arrhenius_params:],M=1)
#for i,m in enumerate(thirdBody_index):
#	unsrtPosterior[str(m)] = pcov_m[i]
	#start = end
#print(unsrtPosterior)
#count = count+len(np.asarray(total_m_params).flatten())
#if len(thermo_index) != 0:
#	pcov_thermo = plotter.extract_block_diag(posterior_cholesky[count:count+7*len(thermo_index),count:count+7*len(thermo_index)],M=7)
#	for i,th in enumerate(thermo_index):
#		unsrtPosterior[str(th)] = list(pcov_thermo)[i]
#	count = count+7*len(thermo_index)
#if len(transport_index) != 0:
#	pcov_transport = plotter.extract_block_diag(posterior_cholesky[count:count+2*len(transport_index),count:count+2*len(transport_index)],M=2)
#	for i,trans in enumerate(transport_index):
#		unsrtPosterior[str(trans)] = list(pcov_transport)[i]
#Finding the posterior uncertainty for reactions

#############################################################
### 	      Building maga database 	   					#
###			for parameters and targets						#
#############################################################
#print(manipulation_dict)

parametric_space = {}
for i in reac:
	temp_database = {}
	temp_database["unsrtDatabase"] = unsrtDatabase[str(i)]
	#temp_database["sensitivityManipulations"] = manipulation_dict["sa"][str(i)]
	#temp_database["optManipulations"] = manipulation_dict["Opt"][str(i)]
	#temp_database["originalManipulations"] = manipulation_dict["Original"][str(i)]
	temp_database["sensitivityManipulations"] = manipulation_dict["sa"]
	temp_database["optManipulations"] = manipulation_dict["Opt"]
	temp_database["originalManipulations"] = manipulation_dict["Original"]
	
	temp_database["optimizedParameter"] = param_opt_dict[str(i)]
	#temp_database["posteriorCovariance"] = unsrtPosterior[str(i)]
	temp_database["posteriorCovariance"] = unsrtPosterior
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
				set_data["SA_simulations"] = sim_dict["sa"][str(case)]
				set_data["Original_simulations"] = sim_dict["Original"][str(case)]
				if "PRS" in optInputs["Type"]["optimization_type"]:
					set_data["Opt_simulations"] = sim_dict["Opt"][str(case)]
					set_data["SA_PRS"] = response_surface_dict["sa"][str(case)]
					set_data["Opt_PRS"] = response_surface_dict["Opt"][str(case)]
					set_data["optimization_type"]= "PRS"
					#set_data["Optimized_simulations"] = sim_dict["Optimized"][str(case)]
					
				else:
					set_data["optimization_type"] = "Direct"
					#set_data["Optimized_simulations"] = sim_dict["Optimized"][str(case)]
				
				set_data["Prior_model_unsrt"] = 0				
				#set_data["Prior_model_unsrt"] = model_priorUnsrt[target_list[case]]
				#set_data["Posterior_model_unsrt"] = model_posteriorUnsrt[target_list[case]]
				set_data["Prior_model_response"] = model_priorResponse[target_list[case]]
				set_data["Posterior_model_response"] = model_posteriorResponse[target_list[case]]
				#print(model_unsrt[target_list[case]])
				#set_data[str(case)] = temp_targets
				set_list.append(set_data)
			else:
				continue
		points[str(set_index)] = set_list
	dataSet[str(d_types)]= points
#print(dataSet)
parameter_list = []
for i in parametric_space:
	p = combustion_variable_class.combustion_variable(i,parametric_space[i])
	parameter_list.append(p)

target_sets = []
for targetClass in dataSet:
	for index,sets in enumerate(dataSet[targetClass]):
		tar = combustion_dataset_class.combustion_dataset(sets,dataSet[targetClass][sets],opt)
		target_sets.append(tar)

plotter.plot(activeParameters,parameter_list,case_dir,target_sets)
#plotter.plot_errors(target_list,opt)
Theta = []
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
#plotter.plot_errors(target_list, opt)
#plotter.plot_vector(reaction_index, opt)
print("generating log files and optimized mechanism file..... \n\n")
#plotter.plot_graphs(target_list,init_guess,posterior_cholesky,opt)
#plotter.plot_rxnUnsrt(mech_file_location,unsrt_location,thermo_file_location,trans_file_location,target_list,sim_type,Theta,Posterior_cov,opt,"True")
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
for sets in target_sets:
	index = sets.set_id[0]
	add = sets.addendum[0]
	plotter.plot_optResults(index,add,sets,opt_tag="on")
print(">>>>>> Program Exiting <<<<<<")
