#python default modules
import numpy as np
import scipy as sp
import scipy.stats as stats
from scipy.optimize import minimize
import os, sys, re, threading, subprocess, time
from sklearn.preprocessing import PolynomialFeatures
from collections import OrderedDict
import matplotlib.pyplot as plt

#program specific modules
import make_input_file
import FlameMaster_in_parallel
import combustion_target_class
import data_management
import simulation_manager
import uncertainty
import MechanismManipulator
import plotter
import Input_file_reader
### KEY WORDS #######
mechanism = "mechanism"
targets_count = "targets_count"
targets = "targets:"
home_dir = os.getcwd()
fuel_name = "fuel"
global_reaction_eqn = "global_reaction"
parallel_thread_count = "parallel_threads"
plus = "plus"
minus = "minus"
unsrt = "uncertainty_data"
thermo_file = "thermo_file"
start_profile = "start_profile"
trans_file = "trans_file"
Sim_type = "sim_count"
order = "Order_of_PRS"
design = "Design_of_PRS"
rxn_total = "total_reactions"
###########
#open the input file and check for arguements
###########

if len(sys.argv) > 1:
	input_file = open(sys.argv[1],'r')
	lines = input_file.readlines()
	print("Input file found\n")
else:
	print("Please enter a valid input file name as arguement. \n For details of preparing the input file, please see the UserManual\n\nProgram exiting")
	exit()


#!!!!!!! GET MECHANISM FILE , number of targets  from the input file !!!!!!!!!

for line in lines:
	######check for comments!!!!!!!!!!
	if "#" in line:
		line = line[:line.index('#')]
		
	word = line.split()
	if mechanism in word:
		mech_file_location = os.path.abspath(word[word.index(mechanism) + 2])
			
	if thermo_file in word:
		thermo_file_location = os.path.abspath(word[word.index(thermo_file) + 2])
	
	if trans_file in word:
		trans_file_location = os.path.abspath(word[word.index(trans_file) + 2])
		
	if start_profile in word:
		startProfile_location = os.path.abspath(word[word.index(start_profile) + 2])
	
	if targets_count in word:
		no_of_targets = int(word[word.index(targets_count) + 2])
	
	if rxn_total in word:
		no_of_reactions = int(word[word.index(rxn_total) + 2])
	
	if fuel_name in word:
		fuel = word[word.index(fuel_name) + 2]
		
	if global_reaction_eqn in word:
		global_reaction = word[word.index(global_reaction_eqn) + 2:]
		
	if parallel_thread_count in word:
		parallel_threads = int(word[word.index(parallel_thread_count) + 2])+1
		
	if unsrt in word:
		unsrt_location = os.path.abspath(word[word.index(unsrt) + 2])
	
	
	if Sim_type in word:
		sim_count = word[word.index(Sim_type) + 2:]

	if design in word:
		design_type = word[word.index(design) + 2]

	if order in word:
		rps_order = int(word[word.index(order) + 2])

	if targets in word:
		start_line = lines.index(line) + 1
		stop_line = start_line + no_of_targets



##!!!  MAKE A LIST OF TARGET CLASS CONTAINING EACH TARGET AS A CASE


target_list = []
c_index = 0
for target in lines[start_line:stop_line]:
	if "#" in target:
		target = target[:target.index('#')]
		
	t = combustion_target_class.combustion_target(target,c_index)
	c_index +=1
	target_list.append(t)
case_dir = range(0,len(target_list))
print("optimization targets identified\n")

##########Obtain uncertainty data from user input file!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

UncertDataSet = uncertainty.uncertaintyData(unsrt_location,mech_file_location,thermo_file_location,trans_file_location);
unsrt_data, reaction_index,fallOffCurve_index, thirdBody_index, thermo_index, transport_index = UncertDataSet.extract_uncertainty();
BranchingRxns = []
pDeptRxns=[]
activeParameters = []
for index in reaction_index:
	rxn_type, IsBranching, branchRxn = unsrt_data[index].getRxnType()
	activeParameters.append(index)
	if IsBranching == "True":
		rxn_list = []
		rxn_list.append(index)
		for i in branchRxn:
			rxn_list.append(i)
		BranchingRxns.append(rxn_list)
	if rxn_type == "p-Dependent":
		pDeptRxns.append(index)

#print("Branching reactions found\n\t{}\n".format(BranchingRxns))
#print("Pressure dependent reactions found\n\t{}\n".format(pDeptRxns))

for i in fallOffCurve_index:
	activeParameters.append(i)
for i in thirdBody_index:
	activeParameters.append(i)
for i in thermo_index:
	activeParameters.append(i)
for i in transport_index:
	activeParameters.append(i)
print("Uncertainty data acquired \n")

sim_type = 'sa'
#print(reaction_index)
########### DO SENSITIVITY ANALYSIS ######################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# CREATE reaction list
#rn = []
#string = ""
#for i in range(1,no_of_reactions+1):
#	string = '{}f'.format(i)
#	rn.append(string)
#uncertainty_data = None


if "SA" not in os.listdir():
	os.mkdir("SA")
	os.chdir("SA")
	os.mkdir("Data")
	os.chdir("Data")
	os.mkdir("NumericalAnalysis")
	os.mkdir("ResponseSurface")
	os.mkdir("Simulations")
	os.chdir("..")
	os.mkdir("Plots")
	os.chdir("Plots")
	os.mkdir("SingularValues")
	os.chdir("..")
else:
	os.chdir("SA")
	if os.path.isfile("progress") == False:
		sim_type = 'opt'
		FlameMaster_Execution_location = simulation_manager.make_directories_for_simulation(case_dir,rps_order, reaction_index,fallOffCurve_index,thirdBody_index,thermo_index,transport_index, mech_file_location, unsrt_data, target_list, fuel, global_reaction, thermo_file_location, trans_file_location,startProfile_location,design_type,sim_type)
		locations = open("locations",'w')
		for i in FlameMaster_Execution_location:
			locations.write(i+"\n")
		locations.close()
		FlameMaster_in_parallel.run_FlameMaster_parallel(FlameMaster_Execution_location, parallel_threads, thermo_file_location, trans_file_location,startProfile_location)
	else:
		progress = open("progress",'r').readlines()
		FlameMaster_Execution_location = open("locations",'r').readlines()
		missing_location = data_management.find_missing_location(FlameMaster_Execution_location, progress)
		FlameMaster_in_parallel.run_FlameMaster_parallel(missing_location, parallel_threads, thermo_file_location, trans_file_location,startProfile_location)


	SA_dir = os.getcwd()
	for case in case_dir:	
		os.chdir("case-"+str(case))
		f = open('../Data/Simulations/sim_data_case-'+str(case)+'.lst','w')
		data_sheet = data_management.generate_target_value_tables(FlameMaster_Execution_location, target_list, case, fuel,sim_type)
		f.write(data_sheet)
		f.close()
		os.chdir(SA_dir)
	
	for case in target_list:
		case.create_response_surface(activeParameters,unsrt_data,1)
		s_a = []
		s_n = []
		s_ea = []
		count = 1
		print(len(case.resCoef[1:]))
		for i in range(len(case.resCoef[1:])//3):
			if count<=len(case.resCoef[1:]):
				#print(count)
				s_a.append(case.resCoef[count])
				s_n.append(case.resCoef[count+1])
				s_ea.append(case.resCoef[count+2])
				count+=3
		fig = plt.figure()
		y_pos = range(0,len(s_a))
		plt.barh(y_pos, s_a, align='center', alpha=0.5)
		plt.yticks(y_pos, activeParameters)
		plt.xlabel('Sensitivity of pre-exponential factor (A)')
		plt.title('Sensitivity Analysis using Response surface method')
		plt.savefig('./Plots/sensitivity_A_'+str(target_list.index(case))+'.png')
	
		fig = plt.figure()
		y_pos = range(0,len(s_n))
		plt.barh(y_pos, s_n, align='center', alpha=0.5)
		plt.yticks(y_pos, activeParameters)
		plt.xlabel('Sensitivity of temperature exponent (n)')
		plt.title('Sensitivity Analysis using Response surface method')
		plt.savefig('./Plots/sensitivity_n_'+str(target_list.index(case))+'.png')
	
		fig = plt.figure()
		y_pos = range(0,len(s_ea))
		plt.barh(y_pos, s_ea, align='center', alpha=0.5)
		plt.yticks(y_pos, activeParameters)
		plt.xlabel('Sensitivity of Activation Energy (Ea)')
		plt.title('Sensitivity Analysis using Response surface method')
		plt.savefig('./Plots/sensitivity_Ea_'+str(target_list.index(case))+'.png')
	
		sorted_Params = []
		Params = case.resCoef[1:]
		Max_param = max(abs(np.asarray(Params)))
		"""
		criteria to select the active reactions
		"""
		
		criteria = 0.001
		index = []
		for i in Params:
			if i>criteria:
				index.append(case.resCoef.index(i))
		print(index)		
	os.chdir("..")

# call simulation manager
# generate beta list and perturbe the mech file accordingly
# DO simulations and find the sensitivity using standard formula (both temp. independent as well as temp dependent)
#Do varification

###[BONUS TASK]create a list of top 20 reactions based on optimization potential
print("Sensitivity analysis completed for all the targets\n")

