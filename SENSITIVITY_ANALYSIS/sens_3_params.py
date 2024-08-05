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
#import yaml
#import ruamel.yaml as yaml
try:
    import ruamel_yaml as yaml
except ImportError:
    from ruamel import yaml
import pandas as pd

sys.path.append('/parallel_yaml_writer.so')
import parallel_yaml_writer

####################################
##  Importing the sampling file   ##
##                                ##
####################################
import reaction_selection as rs
from MechManipulator2_0 import Manipulator
#program specific modules
from copy import deepcopy
from MechanismParser import Parser

import combustion_target_class
import data_management
#import data_management as dm
import simulation_manager2_0 as simulator
import Uncertainty as uncertainty
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
#########################################
###    Reading the input file        ####
#########################################
if len(sys.argv) > 2:
	input_file = open(sys.argv[1],'r')
	optInputs = yaml.safe_load(input_file)
	rxn_list =[i.strip("\n") for i in open(sys.argv[2],"r").readlines()]
	#print(rxn_list)
	print("\n\t########################\n\tInput file and List of reactions are found\n\t########################\n")
	#raise AssertionError("Stop")
elif len(sys.argv) > 1:
	input_file = open(sys.argv[1],'r')
	optInputs = yaml.safe_load(input_file)
	rxn_list = []
	print("\n\t########################\n\tInput file found\n\t########################\n")
else:
	print("Please enter a valid input file name as arguement. \n Two arguments can be passed:\n\t1. Traget opt file\n\t2. List of reactions\nThe code will still work by passing only the first argument\n\nProgram exiting")
	exit()

count = "Counts"
targets = "targets"
countThreads = "parallel_threads"
add = "addendum"
design = "Design_of_PRS"
fuel = "fuel"
dataCounts = optInputs[count]
parallel_threads = dataCounts[countThreads]
print("\n\tParallel threads are {}".format(parallel_threads))
locations = optInputs["Locations"]
targets_count = int(dataCounts["targets_count"])
stats_ = optInputs["Stats"]
design_type = stats_[design]
fuel = optInputs["Inputs"][fuel]

#########################################
###    Reading the mechanism file    ####
#########################################

MECH = locations["mechanism"]
carbon_number = optInputs["SA"]["carbon_number"]

with open(MECH,'r') as file_:
	yaml_mech = file_.read()

mechanism = yaml.safe_load(yaml_mech)
species = mechanism['phases'][0]["species"]
species_data = mechanism["species"]
reactions = mechanism["reactions"]
if len(rxn_list) == 0:
	if carbon_number != 0:
		# get the species list greater than or equal to the predetermined carbon number
		selected_species = rs.species_selection(species,species_data,carbon_number)

		# get the reaction list containing the selected species
		selected_reactions = rs.reaction_selection(selected_species,reactions)

		#get the reaction index for the selected reactions
		reaction_dict = rs.reaction_index(selected_reactions,reactions)

	else:	
		selected_species = species
		selected_reactions = []
		reaction_dict = {}
		for index,rxn in enumerate(reactions):
			selected_reactions.append(rxn["equation"])
			reaction_dict[index+1] = rxn["equation"]
else:
	selected_reactions = rxn_list
	reaction_dict = rs.reaction_index(selected_reactions,reactions)

rxn_type = rs.getRxnType(mechanism,selected_reactions)
string_f = ""
string_g = ""
index_dict = {}
for index in reaction_dict:
	index_dict[reaction_dict[index]] = index
	string_f+=f"{index}\t{reaction_dict[index]}\n"

for rxn in rxn_type:
	string_g+=f"{rxn}\t{rxn_type[rxn]}\n"
f = open("Reaction_dict.txt","w").write(string_f)
g = open("Reaction_type.txt","w").write(string_g)
rxn_dict = {}
rxn_dict["reaction"] = reaction_dict
rxn_dict["type"] = rxn_type
rxn_dict["data"] = rs.getRxnDetails(mechanism,selected_reactions)
string_reaction = ""
for index in reaction_dict:
	string_reaction+=f"{index}\t{reaction_dict[index]}\n"
g = open("selected_rxn.txt","+w").write(string_reaction)
#raise AssertionError("Selected reaction")
####################################################
##  Unloading the target data	  		   ##
## TARGET CLASS CONTAINING EACH TARGET AS A CASE  ##
####################################################

targetLines = open(locations[targets],'r').readlines()
addendum = yaml.safe_load(open(locations[add],'r').read())

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

print("\n\toptimization targets identified\n")
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
print("\nFollowing list is the choosen reactions\n")
print(manipulationDict["activeParameters"])

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



#########################################
###    Creating Design Matrix for    ####
###    sensitivity analysis          ####
#########################################

"""
For sensitivity analysis we create two design matrix
	- for one, we multiply all reactions by a factor of 2
	- for second, we devide all reactions by a factor of 0.5
"""
ap = len(manipulationDict["activeParameters"])
if "DesignMatrix_x0.csv" not in os.listdir():
	design_matrix_x0 = DM.DesignMatrix(unsrt_data,design_type,ap,getSim(ap,design)).getNominal_samples()
	s =""
	for row in design_matrix_x0:
		for element in row:
			s+=f"{element},"
		s+="\n"
	ff = open('DesignMatrix_x0.csv','w').write(s)
	#design_matrix_x0 = DM.DesignMatrix(selected_reactions,design_type,len(reaction_dict)).getNominal_samples()
else:
	design_matrix_file = open("DesignMatrix_x0.csv").readlines()
	design_matrix_x0 = []
	for row in design_matrix_file:
		design_matrix_x0.append([float(ele) for ele in row.strip("\n").strip(",").split(",")])
		
if "DesignMatrix_x2.csv" not in os.listdir():
	design_matrix_x2 = DM.DesignMatrix(unsrt_data,design_type,ap,500).getSamples()
	s =""
	for row in design_matrix_x2:
		for element in row:
			s+=f"{element},"
		s+="\n"
	ff = open('DesignMatrix_x2.csv','w').write(s)
else:
	design_matrix_file = open("DesignMatrix_x2.csv").readlines()
	design_matrix_x2 = []
	for row in design_matrix_file:
		design_matrix_x2.append([float(ele) for ele in row.strip("\n").strip(",").split(",")])

yaml_loc_nominal = []
yaml_loc_nominal.append(MECH)
SSM = simulator.SM(target_list,optInputs,unsrt_data,design_matrix_x2)

if "Perturbed_Mech_SA_3P" not in os.listdir():
	os.mkdir("Perturbed_Mech_SA_3P")
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
	params_yaml = [design_matrix_x2[i:i+chunk_size] for i in range(0, len(design_matrix_x2), chunk_size)]
	count = 0
	yaml_loc = []
	for params in params_yaml:
		
		yaml_list = SSM.getYAML_List(params)
		#yaml_loc = []
		location_mech = []
		index_list = []
		for i,dict_ in enumerate(yaml_list):
			index_list.append(str(count+i))
			location_mech.append(os.getcwd()+"/Perturbed_Mech_SA_3P")
			yaml_loc.append(os.getcwd()+"/Perturbed_Mech_SA_3P/mechanism_"+str(count+i)+".yaml")
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
	for i,sample in enumerate(design_matrix_x2):
		index_list.append(i)
		location_mech.append(os.getcwd()+"/Perturbed_Mech_SA_3P")
		yaml_loc.append(os.getcwd()+"/Perturbed_Mech_SA_3P/mechanism_"+str(i)+".yaml")


#########################################
###    Creating Simulation Field     ####
#########################################
print("\t\t#########################################\n\t\t###    Creating Simulation Field     ####\n\t\t#########################################")

if "SA_3P" not in os.listdir():
	os.mkdir("SA_3P")
	os.chdir("SA_3P")
	os.mkdir("multiply")# we multiply reactions by 2 in this folder
	os.mkdir("divide")# we divide the reactions by 2 in this folder
	os.mkdir("nominal")
	os.mkdir("Data")
	os.chdir("Data")
	os.mkdir("Simulations")
	os.chdir("Simulations")
	os.mkdir("Multiply")
	os.mkdir("Divide")
	os.mkdir("Nominal")
	os.chdir("..")
	os.mkdir("ResponseSurface")
	os.chdir("..")
	os.chdir("multiply")
	SADir = os.getcwd()
else:
	os.chdir("SA_3P")
	os.chdir("multiply")
	SADir = os.getcwd()

#####################################################
#### Multiplying the reactions by a factor of 2   ###
#####################################################
print("\n\t\t#####################################################\n\t\t#### Multiplying the reactions for 3-Params UQ   ###\n\t\t#####################################################")


if os.path.isfile("progress") == False:
	FlameMaster_Execution_location_x2 = simulator.SM(target_list,optInputs,rxn_dict,design_matrix_x2).make_dir_in_parallel(yaml_loc)
	
else:
	print("\t\tProgress file detected")
	progress = open(SADir+"/progress",'r').readlines()
	FlameMaster_Execution_location_x2 = []
	with open(SADir+"/locations") as infile:
		for line in infile:
			FlameMaster_Execution_location_x2.append(line)

os.chdir("..")
os.chdir("nominal")
SADir = os.getcwd()

#############################
#### Nominal simulations  ###
#############################
print("\n\t\t#############################\n\t\t#### Nominal simulations ###\n\t\t#############################")

if os.path.isfile("progress") == False:
	FlameMaster_Execution_location_x0 = simulator.SM(target_list,optInputs,rxn_dict,design_matrix_x0).make_dir_in_parallel(yaml_loc_nominal)
	
else:
	print("\t\tProgress file detected")
	progress = open(SADir+"/progress",'r').readlines()
	FlameMaster_Execution_location_x0 = []
	with open(SADir+"/locations") as infile:
		for line in infile:
			FlameMaster_Execution_location_x0.append(line)
os.chdir("..")
SAdir = os.getcwd()
########################################################
#### collecting sensitivity data from the simulation ###
########################################################
##################################
##### From Nominal Folder #######
##################################
temp_sim_opt_x0 = {}
for case in case_dir:	
	os.chdir("Data/Simulations/Nominal")
	if "sim_data_case-"+str(case)+".lst" in os.listdir():
		ETA_x0 = [float(i.split("\t")[1]) for i in open("sim_data_case-"+str(case)+".lst").readlines()]
		
		folderName_x0 = [float(i.split("\t")[0]) for i in open("sim_data_case-"+str(case)+".lst").readlines()]
		temp_sim_opt_x0[str(case)] = {}
		temp_sim_opt_x0[str(case)]["ETA"] = ETA_x0
		temp_sim_opt_x0[str(case)]["index"] = folderName_x0
		os.chdir(SAdir)
		#print(ETA)
		#raise AssertionError("Generating ETA list for all cases")	
	else:
		os.chdir(SAdir)
		os.chdir("nominal/case-"+str(case))	
		data_sheet_x0,failed_sim_x0,index_x0, ETA_x0,eta_x0 = data_management.generate_SA_target_value_tables(FlameMaster_Execution_location_x0, target_list, case, fuel)
		#print(data_sheet)
		#raise AssertionError("!STOP")
		temp_sim_opt_x0[str(case)] = {}
		temp_sim_opt_x0[str(case)]["ETA"] = eta_x0
		temp_sim_opt_x0[str(case)]["index"] = index_x0
		f = open('../../Data/Simulations/Nominal/sim_data_case-'+str(case)+'.lst','w').write(data_sheet_x0)
		g = open('../../Data/Simulations/Nominal/failed_sim_data_case-'+str(case)+'.lst','w').write(failed_sim_x0)
		#f.write(data_sheet)
		#f.close()
		os.chdir(SAdir)


##################################
##### From Multiply Folder #######
##################################
temp_sim_opt_x2 = {}
for case in case_dir:	
	os.chdir("Data/Simulations/Multiply")
	if "sim_data_case-"+str(case)+".lst" in os.listdir():
		ETA_x2 = [float(i.split("\t")[1]) for i in open("sim_data_case-"+str(case)+".lst").readlines()]
		folderName_x2 = [float(i.split("\t")[0]) for i in open("sim_data_case-"+str(case)+".lst").readlines()]
		temp_sim_opt_x2[str(case)] = {}
		temp_sim_opt_x2[str(case)]["ETA"] = ETA_x2
		temp_sim_opt_x2[str(case)]["index"] = folderName_x2
		os.chdir(SAdir)
		#print(ETA)
		#raise AssertionError("Generating ETA list for all cases")	
	else:
		os.chdir(SAdir)
		os.chdir("multiply/case-"+str(case))	
		data_sheet_x2,failed_sim_x2,index_x2, ETA_x2,eta_x2 = data_management.generate_SA_target_value_tables(FlameMaster_Execution_location_x2, target_list, case, fuel)
		#print(data_sheet)
		#raise AssertionError("!STOP")
		temp_sim_opt_x2[str(case)] = {}
		temp_sim_opt_x2[str(case)]["ETA"] = eta_x2
		temp_sim_opt_x2[str(case)]["index"] = index_x2
		f = open('../../Data/Simulations/Multiply/sim_data_case-'+str(case)+'.lst','w').write(data_sheet_x2)
		g = open('../../Data/Simulations/Multiply/failed_sim_data_case-'+str(case)+'.lst','w').write(failed_sim_x2)
		#f.write(data_sheet)
		#f.close()
		os.chdir(SAdir)


###############################################
##      Generating the response surface      ##
##                                           ##
###############################################

ResponseSurfaces = {}
selected_PRS = {}

for case_index,case in enumerate(temp_sim_opt_x2):
	T = int(target_list[case_index].temperature)
	yData = np.asarray(temp_sim_opt_x2[case]["ETA"]).flatten()
	xData = np.asarray(design_matrix_x2)
	xTrain,xTest,yTrain,yTest = train_test_split(xData,yData,
									random_state=104, 
                                	test_size=0.2, 
                                   	shuffle=True)
	Response = PRS.ResponseSurface(xTrain,yTrain,case,case_index,1)
	Response.create_response_surface()
	Response.test(xTest,yTest)
	Response.plot(case_index)
	#print(Response.case)
	ResponseSurfaces[case_index] = Response
	#print(Response.selection)
	del xTrain,xTest,yTrain,yTest
	if design_type == "A-facto":
		rxn_Sa = {}
		count = 0
		for ind,rxn in enumerate(rxn_list):
			temp = Response.coeff[1:][count]
			count+=1
			rxn_Sa[rxn] = temp
			
		SA_dict = dict(sorted(rxn_Sa.items(), key=lambda item: abs(item[1]),reverse = True))
		selected_PRS[str(case_index)] = rxn_Sa
		sort_rlist = []
		sort_alist = []
		ticks = []
		for ind,rxn in enumerate(SA_dict):
			sort_rlist.append(rxn)
			sort_alist.append(SA_dict[rxn])
			ticks.append(ind)
			
		
		fig = plt.figure()
		y_pos = range(0,len(sort_alist))
		#print(y_pos)
		plt.barh(y_pos,sorted(sort_alist,key =abs), alpha=0.51)
		plt.yticks(y_pos, sort_rlist)
		plt.xlabel(r'normalized sensitivities $S_i = \frac{\partial ln(S_{u}^{o})}{\partial x_i}}$')
		#plt.title('Sensitivity Analysis using Response surface method')
		plt.savefig('../Plots_SA/sensitivity_A_'+str(case_index)+'.png',bbox_inches="tight")
		plt.close()
			
	else:
		rxn_Sa = {}
		count = 0
		for ind,rxn in enumerate(rxn_list):
			temp = []
			temp.append(Response.coeff[1:][count])
			temp.append(Response.coeff[1:][count+1])
			temp.append(Response.coeff[1:][count+2])
			count+=3
			rxn_Sa[rxn] = temp
			
		SA_dict = dict(sorted(rxn_Sa.items(), key=lambda item: abs(item[1][0]),reverse = True))
		selected_PRS[str(case_index)] = rxn_Sa
		sort_rlist = []
		sort_alist = []
		sort_nlist = []
		sort_ealist = []
		ticks = []
		for ind,rxn in enumerate(SA_dict):
			sort_rlist.append(rxn)
			sort_alist.append(SA_dict[rxn][0])
			sort_nlist.append(SA_dict[rxn][1])
			sort_ealist.append(SA_dict[rxn][2])
			ticks.append(ind)
			
		
		fig = plt.figure()
		y_pos = range(0,len(sort_alist))
		#print(y_pos)
		plt.barh(y_pos,sorted(sort_alist,key =abs), alpha=0.51)
		plt.yticks(y_pos, sort_rlist)
		plt.xlabel(r'normalized sensitivities $S_i = \frac{\partial ln(S_{u}^{o})}{\partial x_i}}$')
		#plt.title('Sensitivity Analysis using Response surface method')
		plt.savefig('../Plots_SA/sensitivity_A_'+str(case_index)+'.png',bbox_inches="tight")
		plt.close()
		"""
		Plotting the sensitivity in same fig
		"""
		
		#print(sort_rlist)
		fake_data = pd.DataFrame({"index": list(sort_rlist), 0: sort_alist , 1: sort_nlist, 2: np.asarray(sort_ealist)*10})
		fake_data.set_index("index",drop=False)
		fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(8, 8), frameon=False)
		
		fake_data[0].plot.barh(ax=ax1)
		fake_data[1].plot.barh(ax=ax2)
		fake_data[2].plot.barh(ax=ax3)
		
		ax1.set_yticks(ticks,sort_rlist)
		ax1.set_xlabel(r'$\partial ln(\eta) / \partial \zeta_{\alpha}$')
		ax2.set_xlabel(r'$\partial ln(\eta) / \partial \zeta_{n}$')
		ax3.set_xlabel(r'$\partial ln(\eta) / \partial \zeta_{\epsilon} (\times 10^{-1})$')
		fig.savefig('../Plots_SA/sensitivity_'+str(case_index)+'.png',bbox_inches="tight")
		plt.close()
	
	string_FM = "Sensitivity Analysis (using Cantera) Tig:\n"
	for ind,rxn in enumerate(sort_rlist):
		string_FM +=f"\t{sort_alist[ind]:.8f}\t{index_dict[rxn.split(':')[0]]}\t{rxn}\n"

	g = open(f"Data/FM_sensitivity_T_{T}_case_{case_index}.txt","w").write(string_FM)
	
os.chdir("..")
if "sens_3p_parameters.pkl" not in os.listdir():
	with open('sens_3p_parameters.pkl', 'wb') as file_:
			pickle.dump(selected_PRS, file_)

raise AssertionError("HOHOHO, 3-Param unsrt analysis done!!")
for case in case_dir:
	#print(case)
	T = target_list[case].temperature
	index = temp_sim_opt_x2[str(case)]["index"]
	#index_0 = temp_sim_opt_x0[str(case)]["index"]
	multiply = np.asarray(temp_sim_opt_x2[str(case)]["ETA"])
	#divide = np.asarray(temp_sim_opt_x0_5[str(case)]["ETA"])
	nominal = np.asarray(temp_sim_opt_x0[str(case)]["ETA"])
	#case_sens = np.log(multiply/divide)/np.log(4)
	#case_sens_0 = np.log(multiply/nominal)/np.log(2)
	case_FM = (multiply - nominal)/nominal
		
	sens_FM = {}
	for i,ind in enumerate(index):
		sens_FM[ind] = case_FM[i]

	sensitivity_FM = dict(sorted(sens_FM.items(), key=lambda item: abs(item[1]),reverse = True))

	string_FM = "Sensitivity Analysis (using Cantera) Tig:\n"
	for ind,sens in enumerate(sensitivity_FM):
		
		string_FM +=f"\t{sensitivity_FM[sens]:.8f}\t{int(int(sens)+1)}\t{reaction_dict[sens]}\n"

	g = open(f"Data/FM_sensitivity_T_{T}.txt","w").write(string_FM)
	
print("\n\t################################\n\tSENSITIVITY ANALYSIS DONE!!\n\t################################\n\t")	
