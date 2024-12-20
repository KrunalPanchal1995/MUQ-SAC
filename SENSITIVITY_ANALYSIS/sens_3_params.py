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
import yaml
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
if "sensitive_parameters" not in stats_:
	stats_["sensitive_parameters"] = "Arrhenius"
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

########################################################################
## Analysis for temperature-exponent (n) and activation energies (Ea) ##
##                                                                    ##
########################################################################
#Finding sigma_n and sigma_Ea from uncertainty analysis
#sigma_n = []
#sigma_Ea = []
#for rxn in unsrt_data:
	

######################################################################
##  CREATING A DICTIONARY CONTAINING ALL THE DATA FROM UNSRT CLASS  ##
##																    ##
######################################################################
manipulationDict = {}
selection = []
Cholesky_list = []
zeta_list = []
activeParameters = []
P_nominal_list = []
P_upper = []
P_lower = []
P_multiply_A = []
P_multiply_n = []
P_multiply_Ea = []
rxn_list = []
rIndex = []
sigma_n = []
sigma_Ea = []
zeta_list_A = []
zeta_list_B = []
zeta_list_C = []
for rxn in unsrt_data:
	activeParameters.extend(unsrt_data[rxn].activeParameters)
ap = len(activeParameters)
zeta_b_max = DM.DesignMatrix(unsrt_data,design_type,ap,1).getB_TYPE_fSAC(1000,tag="MAX_N")
zeta_c_max,gen = DM.DesignMatrix(unsrt_data,design_type,ap,1).getClassC_Curves(1,generator=np.array([[1,-1]]))


if "zeta" in stats_["sensitive_parameters"]:
	z_a = np.asarray([1,0,0])
	z_n = np.asarray([0,1,0])
	z_e = np.asarray([0,0,1])
	p_a = p_n = p_e = np.asarray([1,1,1])
else:
	z_a = z_n = z_e = np.asarray([1,1,1])
	p_a = np.asarray([1,0,0])
	p_n = np.asarray([0,1,0])
	p_e = np.asarray([0,0,1])
	
for rxn in unsrt_data:
	rxn_list.append(rxn)
	selection.extend(unsrt_data[rxn].selection)
	Cholesky_list.append(unsrt_data[rxn].cholskyDeCorrelateMat)
	covariance_mat = np.dot(unsrt_data[rxn].L,unsrt_data[rxn].L.T)
	sigma_vector = np.diag(covariance_mat)
	cov = unsrt_data[rxn].L
	zeta = np.asarray(unsrt_data[rxn].perturb_factor)
	Po = np.asarray(unsrt_data[rxn].nominal)
	zeta_A = z_a*zeta
	zeta_B = z_n*zeta_b_max[rxn]
	zeta_C = z_e*zeta_c_max[rxn][0]
	zeta_list_A.append(zeta_A)
	zeta_list_B.append(zeta_B)
	zeta_list_C.append(zeta_C)
	P_upper.append(np.asarray(Po + np.array([1,1,1])*np.asarray(cov.dot(zeta)).flatten()))
	P_lower.append(np.asarray(Po - np.array([1,1,1])*np.asarray(cov.dot(zeta)).flatten()))
	P_multiply_A.append(np.asarray(Po + p_a*np.asarray(cov.dot(zeta_A)).flatten()))
	P_multiply_n.append(np.asarray(Po + p_n*np.asarray(cov.dot(zeta_B)).flatten()))
	P_multiply_Ea.append(np.asarray(Po + p_e*np.asarray(cov.dot(zeta_C)).flatten()))

	sigma_n.append(sigma_vector[1])
	sigma_Ea.append(sigma_vector[2])
	zeta_list.append(unsrt_data[rxn].perturb_factor)
	
	#activeParameters.extend(unsrt_data[rxn].activeParameters)
	P_nominal_list.append(unsrt_data[rxn].nominal)
	rIndex.append(unsrt_data[rxn].rIndex)
#print(zeta_list)
#raise AssertionError("Stop")
manipulationDict["selection"] = deepcopy(selection)#.deepcopy()
manipulationDict["Cholesky"] = deepcopy(Cholesky_list)#.deepcopy()
manipulationDict["zeta"] = deepcopy(zeta_list)#.deepcopy()
manipulationDict["activeParameters"] = deepcopy(activeParameters)#.deepcopy()
manipulationDict["nominal"] = deepcopy(P_nominal_list)#.deepcopy()
print("\nFollowing list is the choosen reactions\n")
print(manipulationDict["activeParameters"])


############################################
### Plotting all the samples taken for    ##
### Sensitivity analysis                  ##
############################################
global T, theta
T = np.linspace(300,2500,100)
theta = np.array([T/T,np.log(T),-1/T])

def getUnsrtLimit(Po,P_u,P_l):
	K_o = np.asarray([ i.dot(Po) for i in theta.T ]).flatten()
	K_u = np.asarray([ i.dot(P_u) for i in theta.T ]).flatten()
	K_l = np.asarray([ i.dot(P_l) for i in theta.T ]).flatten()
	return K_o,K_u,K_l
	
def getKappa(P):
	K = np.asarray([ i.dot(P) for i in theta.T ]).flatten()
	return np.exp(K)


for index,rxn in enumerate(unsrt_data):
	Kappa_o, Kappa_U, Kappa_L = getUnsrtLimit(P_nominal_list[index],P_upper[index],P_lower[index])
	Kappa_A = getKappa(P_multiply_A[index])
	Kappa_n = getKappa(P_multiply_n[index])
	Kappa_Ea = getKappa(P_multiply_Ea[index])
	fig = plt.figure()
	plt.plot(1/T,Kappa_o,'b-',label="Nominal Curve")
	plt.plot(1/T,Kappa_U,'k--',label="Unsrt Limits")
	plt.plot(1/T,Kappa_L,'k--')
	plt.plot(1/T,Kappa_A,'r-',label="A-factor perturbation")
	plt.plot(1/T,Kappa_n,'g-',label="n perturbation")
	plt.plot(1/T,Kappa_Ea,'y-',label="Ea perturbation")
	plt.legend()
	plt.savefig(f"Plots/{rxn}.pdf",bbox_inches='tight')
	plt.close()
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

if "DesignMatrix_A.csv" not in os.listdir():
	flag = "A"
	design_matrix_A = DM.DesignMatrix(unsrt_data,design_type,ap).getSA3P_samples(zeta_list_A,flag)
	s =""
	for row in design_matrix_A:
		for element in row:
			s+=f"{element},"
		s+="\n"
	ff = open('DesignMatrix_A.csv','w').write(s)
else:
	design_matrix_file = open("DesignMatrix_A.csv").readlines()
	design_matrix_A = []
	for row in design_matrix_file:
		design_matrix_A.append([float(ele) for ele in row.strip("\n").strip(",").split(",")])
		
if "DesignMatrix_n.csv" not in os.listdir():
	flag = "n"
	design_matrix_n = DM.DesignMatrix(unsrt_data,design_type,ap).getSA3P_samples(zeta_list_B,flag)
	s =""
	for row in design_matrix_n:
		for element in row:
			s+=f"{element},"
		s+="\n"
	ff = open('DesignMatrix_n.csv','w').write(s)
else:
	design_matrix_file = open("DesignMatrix_n.csv").readlines()
	design_matrix_n = []
	for row in design_matrix_file:
		design_matrix_n.append([float(ele) for ele in row.strip("\n").strip(",").split(",")])

if "DesignMatrix_Ea.csv" not in os.listdir():
	flag = "Ea"
	design_matrix_Ea = DM.DesignMatrix(unsrt_data,design_type,ap).getSA3P_samples(zeta_list_C,flag)
	s =""
	for row in design_matrix_Ea:
		for element in row:
			s+=f"{element},"
		s+="\n"
	ff = open('DesignMatrix_Ea.csv','w').write(s)
else:
	design_matrix_file = open("DesignMatrix_Ea.csv").readlines()
	design_matrix_Ea = []
	for row in design_matrix_file:
		design_matrix_Ea.append([float(ele) for ele in row.strip("\n").strip(",").split(",")])
		
yaml_loc_nominal = []
yaml_loc_nominal.append(MECH)
SSM = simulator.SM(target_list,optInputs,unsrt_data,design_matrix_x2)
if "Perturbed_Mech_SA_3P_BruteForce" not in os.listdir():
	os.mkdir("Perturbed_Mech_SA_3P_BruteForce")
	os.mkdir("Perturbed_Mech_SA_3P_BruteForce/A_factor")
	os.mkdir("Perturbed_Mech_SA_3P_BruteForce/n")
	os.mkdir("Perturbed_Mech_SA_3P_BruteForce/Ea")
	print("\nPerturbing the Mechanism files for 3-Parameter brute force sensitivity analysis\n")
	chunk_size = 500
	params_yaml_A = [design_matrix_A[i:i+chunk_size] for i in range(0, len(design_matrix_A), chunk_size)]
	params_yaml_n = [design_matrix_n[i:i+chunk_size] for i in range(0, len(design_matrix_n), chunk_size)]
	params_yaml_Ea = [design_matrix_Ea[i:i+chunk_size] for i in range(0, len(design_matrix_Ea), chunk_size)]
	count = 0
	yaml_loc_A = []
	yaml_loc_n = []
	yaml_loc_Ea = []
	for params in params_yaml_A:
		yaml_list = SSM.getYAML_List(params,"A")
		#yaml_loc = []
		location_mech = []
		index_list = []
		for i,dict_ in enumerate(yaml_list):
			index_list.append(str(count+i))
			location_mech.append(os.getcwd()+"/Perturbed_Mech_SA_3P_BruteForce/A_factor")
			yaml_loc_A.append(os.getcwd()+"/Perturbed_Mech_SA_3P_BruteForce/A_factor/mechanism_"+str(count+i)+".yaml")
		count+=len(yaml_list)
		#gen_flag = False
		#SSM.getPerturbedMechLocation(yaml_list,location_mech,index_list)
		SSM.getPerturbedMechLocation(yaml_list,location_mech,index_list)
		print(f"\n\tGenerated {count} files!!\n")
	for params in params_yaml_n:
		yaml_list = SSM.getYAML_List(params,"n")
		#yaml_loc = []
		location_mech = []
		index_list = []
		for i,dict_ in enumerate(yaml_list):
			index_list.append(str(count+i))
			location_mech.append(os.getcwd()+"/Perturbed_Mech_SA_3P_BruteForce/n")
			yaml_loc_n.append(os.getcwd()+"/Perturbed_Mech_SA_3P_BruteForce/n/mechanism_"+str(count+i)+".yaml")
		count+=len(yaml_list)
		#gen_flag = False
		#SSM.getPerturbedMechLocation(yaml_list,location_mech,index_list)
		SSM.getPerturbedMechLocation(yaml_list,location_mech,index_list)
		print(f"\n\tGenerated {count} files!!\n")
	for params in params_yaml_Ea:
		yaml_list = SSM.getYAML_List(params,"Ea")
		#yaml_loc = []
		location_mech = []
		index_list = []
		for i,dict_ in enumerate(yaml_list):
			index_list.append(str(count+i))
			location_mech.append(os.getcwd()+"/Perturbed_Mech_SA_3P_BruteForce/Ea")
			yaml_loc_Ea.append(os.getcwd()+"/Perturbed_Mech_SA_3P_BruteForce/Ea/mechanism_"+str(count+i)+".yaml")
		count+=len(yaml_list)
		#gen_flag = False
		#SSM.getPerturbedMechLocation(yaml_list,location_mech,index_list)
		SSM.getPerturbedMechLocation(yaml_list,location_mech,index_list)
		print(f"\n\tGenerated {count} files!!\n")
	print("\n\tGenerated the YAML files required for simulations!!\n")
else:
	print("\nYAML files for 3-Parameter Bruteforce analysis is already generated!!")
	yaml_loc_A = []
	yaml_loc_n = []
	yaml_loc_Ea = []
	location_mech_A = []
	location_mech_n = []
	location_mech_Ea = []
	index_list_A = []
	index_list_n = []
	index_list_Ea = []
	for i,sample in enumerate(design_matrix_A):
		index_list_A.append(i)
		location_mech_A.append(os.getcwd()+"/Perturbed_Mech_SA_3P_BruteForce/A_factor")
		yaml_loc_A.append(os.getcwd()+"/Perturbed_Mech_SA_3P_BruteForce/A_factor/mechanism_"+str(i)+".yaml")
	for i,sample in enumerate(design_matrix_n):
		index_list_n.append(i)
		location_mech_n.append(os.getcwd()+"/Perturbed_Mech_SA_3P_BruteForce/n")
		yaml_loc_n.append(os.getcwd()+"/Perturbed_Mech_SA_3P_BruteForce/n/mechanism_"+str(i)+".yaml")
	for i,sample in enumerate(design_matrix_Ea):
		index_list_Ea.append(i)
		location_mech_Ea.append(os.getcwd()+"/Perturbed_Mech_SA_3P_BruteForce/Ea")
		yaml_loc_Ea.append(os.getcwd()+"/Perturbed_Mech_SA_3P_BruteForce/Ea/mechanism_"+str(i)+".yaml")

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
	os.mkdir("multiply_A")
	os.mkdir("multiply_n")
	os.mkdir("multiply_Ea")
	os.mkdir("multiply")# we multiply reactions by 2 in this folder
	os.mkdir("divide")# we divide the reactions by 2 in this folder
	os.mkdir("nominal")
	os.mkdir("Data")
	os.chdir("Data")
	os.mkdir("Simulations")
	os.chdir("Simulations")
	os.mkdir("Multiply")
	os.mkdir("Multiply_A")
	os.mkdir("Multiply_n")
	os.mkdir("Multiply_Ea")
	os.mkdir("Divide")
	os.mkdir("Nominal")
	os.chdir("..")
	os.mkdir("ResponseSurface")
	os.chdir("..")
	os.chdir("multiply_A")
	SADir = os.getcwd()
else:
	os.chdir("SA_3P")
	os.chdir("multiply_A")
	SADir = os.getcwd()

################################################################################
#### Multiplying the A-factor of reactions within the uncertainty limits     ###
################################################################################
print("\n\t\t#####################################################\n\t\t#### Multiplying the A-factor of all reactions for 3-Params UQ   ###\n\t\t#####################################################")


if os.path.isfile("progress") == False:
	FlameMaster_Execution_location_A = simulator.SM(target_list,optInputs,rxn_dict,design_matrix_A).make_dir_in_parallel(yaml_loc_A)
	
else:
	print("\t\tProgress file detected")
	progress = open(SADir+"/progress",'r').readlines()
	FlameMaster_Execution_location_A = []
	with open(SADir+"/locations") as infile:
		for line in infile:
			FlameMaster_Execution_location_A.append(line)

os.chdir("..")
os.chdir("multiply_n")
SADir = os.getcwd()

####################################################################################
#### Multiplying the "n" parameter of reactions within the uncertainty limits    ###
####################################################################################
print("\n\t\t#####################################################\n\t\t#### Multiplying the n of all reactions for 3-Params UQ   ###\n\t\t#####################################################")


if os.path.isfile("progress") == False:
	FlameMaster_Execution_location_n = simulator.SM(target_list,optInputs,rxn_dict,design_matrix_n).make_dir_in_parallel(yaml_loc_n)
	
else:
	print("\t\tProgress file detected")
	progress = open(SADir+"/progress",'r').readlines()
	FlameMaster_Execution_location_n = []
	with open(SADir+"/locations") as infile:
		for line in infile:
			FlameMaster_Execution_location_n.append(line)

os.chdir("..")
os.chdir("multiply_Ea")
SADir = os.getcwd()

#######################################################################################
#### Multiplying the "Ea" paramteter of reactions within the uncertainty limits     ###
#######################################################################################
print("\n\t\t#####################################################\n\t\t#### Multiplying the Ea of all reactions for 3-Params UQ   ###\n\t\t#####################################################")


if os.path.isfile("progress") == False:
	FlameMaster_Execution_location_Ea = simulator.SM(target_list,optInputs,rxn_dict,design_matrix_Ea).make_dir_in_parallel(yaml_loc_Ea)
	
else:
	print("\t\tProgress file detected")
	progress = open(SADir+"/progress",'r').readlines()
	FlameMaster_Execution_location_Ea = []
	with open(SADir+"/locations") as infile:
		for line in infile:
			FlameMaster_Execution_location_Ea.append(line)

os.chdir("..")
os.chdir("multiply")
SADir = os.getcwd()


#####################################################
#### Multiplying the reactions                    ###
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
		temp_sim_opt_x0[str(case)]["ETA"] = ETA_x0
		temp_sim_opt_x0[str(case)]["index"] = index_x0
		f = open('../../Data/Simulations/Nominal/sim_data_case-'+str(case)+'.lst','w').write(data_sheet_x0)
		g = open('../../Data/Simulations/Nominal/failed_sim_data_case-'+str(case)+'.lst','w').write(failed_sim_x0)
		#f.write(data_sheet)
		#f.close()
		os.chdir(SAdir)

####################################
##### From Multiply_A Folder #######
####################################
temp_sim_opt_A = {}
for case in case_dir:	
	os.chdir("Data/Simulations/Multiply_A")
	if "sim_data_case-"+str(case)+".lst" in os.listdir():
		ETA_A = [float(i.split("\t")[1]) for i in open("sim_data_case-"+str(case)+".lst").readlines()]
		folderName_A = [float(i.split("\t")[0]) for i in open("sim_data_case-"+str(case)+".lst").readlines()]
		temp_sim_opt_A[str(case)] = {}
		temp_sim_opt_A[str(case)]["ETA"] = ETA_A
		temp_sim_opt_A[str(case)]["index"] = folderName_A
		os.chdir(SAdir)
		#print(ETA)
		#raise AssertionError("Generating ETA list for all cases")	
	else:
		os.chdir(SAdir)
		os.chdir("multiply_A/case-"+str(case))	
		data_sheet_A,failed_sim_A,index_A, ETA_A,eta_A = data_management.generate_SA_target_value_tables(FlameMaster_Execution_location_A, target_list, case, fuel)
		#print(data_sheet)
		#raise AssertionError("!STOP")
		temp_sim_opt_A[str(case)] = {}
		temp_sim_opt_A[str(case)]["ETA"] = ETA_A
		temp_sim_opt_A[str(case)]["index"] = index_A
		f = open('../../Data/Simulations/Multiply_A/sim_data_case-'+str(case)+'.lst','w').write(data_sheet_A)
		g = open('../../Data/Simulations/Multiply_A/failed_sim_data_case-'+str(case)+'.lst','w').write(failed_sim_A)
		#f.write(data_sheet)
		#f.close()
		os.chdir(SAdir)

####################################
##### From Multiply_n Folder #######
####################################
temp_sim_opt_n = {}
for case in case_dir:	
	os.chdir("Data/Simulations/Multiply_n")
	if "sim_data_case-"+str(case)+".lst" in os.listdir():
		ETA_n = [float(i.split("\t")[1]) for i in open("sim_data_case-"+str(case)+".lst").readlines()]
		folderName_n = [float(i.split("\t")[0]) for i in open("sim_data_case-"+str(case)+".lst").readlines()]
		temp_sim_opt_n[str(case)] = {}
		temp_sim_opt_n[str(case)]["ETA"] = ETA_n
		temp_sim_opt_n[str(case)]["index"] = folderName_n
		os.chdir(SAdir)
		#print(ETA)
		#raise AssertionError("Generating ETA list for all cases")	
	else:
		os.chdir(SAdir)
		os.chdir("multiply_n/case-"+str(case))	
		data_sheet_n,failed_sim_n,index_n, ETA_n,eta_n = data_management.generate_SA_target_value_tables(FlameMaster_Execution_location_n, target_list, case, fuel)
		#print(data_sheet)
		#raise AssertionError("!STOP")
		temp_sim_opt_n[str(case)] = {}
		temp_sim_opt_n[str(case)]["ETA"] = ETA_n
		temp_sim_opt_n[str(case)]["index"] = index_n
		f = open('../../Data/Simulations/Multiply_n/sim_data_case-'+str(case)+'.lst','w').write(data_sheet_n)
		g = open('../../Data/Simulations/Multiply_n/failed_sim_data_case-'+str(case)+'.lst','w').write(failed_sim_n)
		#f.write(data_sheet)
		#f.close()
		os.chdir(SAdir)

#####################################
##### From Multiply_Ea Folder #######
#####################################
temp_sim_opt_Ea = {}
for case in case_dir:	
	os.chdir("Data/Simulations/Multiply_Ea")
	if "sim_data_case-"+str(case)+".lst" in os.listdir():
		ETA_Ea = [float(i.split("\t")[1]) for i in open("sim_data_case-"+str(case)+".lst").readlines()]
		folderName_Ea = [float(i.split("\t")[0]) for i in open("sim_data_case-"+str(case)+".lst").readlines()]
		temp_sim_opt_Ea[str(case)] = {}
		temp_sim_opt_Ea[str(case)]["ETA"] = ETA_Ea
		temp_sim_opt_Ea[str(case)]["index"] = folderName_Ea
		os.chdir(SAdir)
		#print(ETA)
		#raise AssertionError("Generating ETA list for all cases")	
	else:
		os.chdir(SAdir)
		os.chdir("multiply_Ea/case-"+str(case))	
		data_sheet_Ea,failed_sim_Ea,index_Ea, ETA_Ea,eta_Ea = data_management.generate_SA_target_value_tables(FlameMaster_Execution_location_Ea, target_list, case, fuel)
		#print(data_sheet)
		#raise AssertionError("!STOP")
		temp_sim_opt_Ea[str(case)] = {}
		temp_sim_opt_Ea[str(case)]["ETA"] = ETA_Ea
		temp_sim_opt_Ea[str(case)]["index"] = index_Ea
		f = open('../../Data/Simulations/Multiply_Ea/sim_data_case-'+str(case)+'.lst','w').write(data_sheet_Ea)
		g = open('../../Data/Simulations/Multiply_Ea/failed_sim_data_case-'+str(case)+'.lst','w').write(failed_sim_Ea)
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
		ETA_x2 = [ np.log(float(i.split("\t")[1])*10) for i in open("sim_data_case-"+str(case)+".lst").readlines()]
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
if "SensitivityCoeffs_PRS" not in os.listdir("Data/"):
	os.mkdir("Data/SensitivityCoeffs_PRS")
if "SensitivityCoeffs" not in os.listdir("Data/"):
	os.mkdir("Data/SensitivityCoeffs")
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

	g = open(f"Data/SensitivityCoeffs_PRS/FM_sensitivity_T_{T}_case_{case_index}.txt","w").write(string_FM)


selected_BRUTE_FORCE_PARAMETERS = {}
for case_index,case in enumerate(temp_sim_opt_A):
	if design_type == "A-facto":
		rxn_Sa = {}
		rxn_Sa_1 = {}
		count = 0
		T = target_list[case_index].temperature
		multiply_A = np.asarray(temp_sim_opt_A[str(case)]["ETA"])
		nominal = np.asarray(temp_sim_opt_x0[str(case)]["ETA"])
		index = temp_sim_opt_A[str(case)]["index"]
		SA_coeff = []
		for rxn_index,rxn in enumerate(unsrt_data):
			k_perturbed.append(float(np.exp(P_multiply_A[rxn_index][0])))
			k_o = float(np.exp(P_nominal_list[rxn_index][0]))
			SA_coeff.append((k_o/(k_perturbed-k_o))*(multiply_A[rxn_index] - nominal)/nominal)
			#print((k_o/(k_perturbed-k_o)),(multiply_A - nominal)/nominal)
			SA_coeff_without_k_perturbed = (multiply_A[rxn_index] - nominal)/nominal
		for ind,rxn in enumerate(rxn_list):
			temp = SA_coeff[count]
			count+=1
			rxn_Sa[rxn] = temp
		
		#for ind,rxn in enumerate(rxn_list):
		#	temp = SA_coeff_without_k_perturbed[count]
		#	count+=1
		#	rxn_Sa_1[rxn] = temp
		
		
		SA_dict = dict(sorted(rxn_Sa.items(), key=lambda item: abs(item[1]),reverse = True))
		#SA_dict_1 = dict(sorted(rxn_Sa_1.items(), key=lambda item: abs(item[1]),reverse = True))
		selected_BRUTE_FORCE_PARAMETERS[str(case_index)] = rxn_Sa
		#sort_rlist = []
		#sort_alist_1 = []
		#ticks = []
		#for ind,rxn in enumerate(SA_dict):
		#	sort_rlist.append(rxn)
		#	sort_alist_1.append(SA_dict_1[rxn])
		#	ticks.append(ind)
		
		
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
		rxn_Sa_1 = {}
		count = 0
		index = temp_sim_opt_A[str(case)]["index"]
		multiply_A = np.asarray(temp_sim_opt_A[str(case)]["ETA"])
		nominal = np.asarray(temp_sim_opt_x0[str(case)]["ETA"])
		
		index_n = temp_sim_opt_n[str(case)]["index"]
		multiply_n = np.asarray(temp_sim_opt_n[str(case)]["ETA"])
		
		index_ea = temp_sim_opt_Ea[str(case)]["index"]
		multiply_Ea = np.asarray(temp_sim_opt_Ea[str(case)]["ETA"])
		#Sensitivity analysis for A parameter
		T_ = float(target_list[case_index].temperature)
		SA_coeff_A = []
		SA_coeff_n = []
		SA_coeff_Ea = []
		SA_coeff_without_k_perturbed = []
		for rxn_index, rxn in enumerate(unsrt_data):
			k_perturbed = getKappa(P_multiply_A[rxn_index])
			k_o = getKappa(P_nominal_list[rxn_index])
			normalized_A = np.asarray((k_o/(k_perturbed-k_o))).flatten()
			fact_A = np.asarray((multiply_A[rxn_index] - nominal)/nominal).flatten()
			#SA_coeff_A.append((k_o/(k_perturbed-k_o))*((multiply_A[rxn_index] - nominal)/nominal))
			SA_coeff_A.append(np.asarray(max(list(normalized_A))*fact_A)[0])
			SA_coeff_without_k_perturbed.append(np.asarray((multiply_A[rxn_index] - nominal)/nominal).flatten()[0])
			
			#Sensitivity analysis for n parameter
			#T = target_list[case_index].temperature
			k_perturbed = getKappa(P_multiply_n[rxn_index])
			k_o = getKappa(P_nominal_list[rxn_index])
			#print(len(k_o),len(k_perturbed),len(multiply_n))
			#print(multiply_n)
			fact_n = np.asarray(((multiply_n[rxn_index] - nominal)/nominal)).flatten()
			normalized_n = np.asarray((k_o/(k_perturbed-k_o))).flatten()
			max_morm_n = max(list(normalized_n))
			#SA_coeff_n.append((k_o/(k_perturbed-k_o))*((multiply_n[rxn_index] - nominal)/nominal))
			SA_coeff_n.append(max_morm_n*fact_n[0])
			#Sensitivity analysis for Ea parameter
			#T = target_list[case_index].temperature
			k_perturbed = getKappa(P_multiply_Ea[rxn_index])
			k_o = getKappa(P_nominal_list[rxn_index])
			fact_ea = np.asarray(((multiply_Ea[rxn_index] - nominal)/nominal)).flatten()
			normalized_ea = np.asarray((k_o/(k_perturbed-k_o))).flatten()
			max_morm_ea = max(list(normalized_ea))
			#SA_coeff_Ea.append((k_o/(k_perturbed-k_o))*(multiply_Ea[rxn_index] - nominal)/nominal)
			SA_coeff_Ea.append(max_morm_ea*fact_ea[0])
		

		
		for ind,rxn in enumerate(rxn_list):
			temp = []
			temp_1 = []
			temp.append(SA_coeff_A[count])
			temp_1.append(SA_coeff_without_k_perturbed[count])
			temp.append(SA_coeff_n[count])
			temp.append(SA_coeff_Ea[count])
			count+=1
			rxn_Sa[rxn] = temp
			rxn_Sa_1[rxn] = temp_1
			
		SA_dict = dict(sorted(rxn_Sa.items(), key=lambda item: abs(item[1][0]),reverse = True))
		SA_dict_1 = dict(sorted(rxn_Sa_1.items(), key=lambda item: abs(item[1][0]),reverse = True))
		#print(SA_dict)
		selected_BRUTE_FORCE_PARAMETERS[str(case_index)] = rxn_Sa
		sort_rlist = []
		sort_alist = []
		sort_alist_1 = []
		sort_nlist = []
		sort_ealist = []
		ticks = []
		for ind,rxn in enumerate(SA_dict):
			sort_rlist.append(rxn)
			sort_alist_1.append(SA_dict_1[rxn][0])
			sort_alist.append(SA_dict[rxn][0])
			sort_nlist.append(SA_dict[rxn][1])
			sort_ealist.append(SA_dict[rxn][2])
			ticks.append(ind)
			
		
		#print(sort_alist)
		fig = plt.figure()
		#y_pos = [i for i in range(0,len(sort_alist))]
		y_pos = range(0,len(sort_alist))
		#print(y_pos)
		plt.barh(y_pos,sorted(sort_alist,key =abs), alpha=0.51)
		#plt.barh(y_pos,sort_alist, alpha=0.51)
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
	string_FM_1 = "Sensitivity Analysis (using Cantera) Tig:\n"
	for ind,rxn in enumerate(sort_rlist):
		string_FM +=f"\t{sort_alist[ind]:.8f}\t{index_dict[rxn.split(':')[0]]}\t{rxn}\n"
		string_FM_1 +=f"\t{sort_alist_1[ind]:.8f}\t{index_dict[rxn.split(':')[0]]}\t{rxn}\n"
		
	string_FM_n = "Sensitivity Analysis (using Cantera) Tig:\n"
	for ind,rxn in enumerate(sort_rlist):
		string_FM_n +=f"\t{sort_nlist[ind]:.8f}\t{index_dict[rxn.split(':')[0]]}\t{rxn}\n"
	
	string_FM_Ea = "Sensitivity Analysis (using Cantera) Tig:\n"
	for ind,rxn in enumerate(sort_rlist):
		string_FM_Ea +=f"\t{sort_ealist[ind]:.8f}\t{index_dict[rxn.split(':')[0]]}\t{rxn}\n"

	g = open(f"Data/SensitivityCoeffs/FM_sensitivity_T_{T_}_case_{case_index}.txt","w").write(string_FM)
	g = open(f"Data/SensitivityCoeffs/FM_sensitivity_T_{T_}_case_{case_index}_1.txt","w").write(string_FM_1)
	g_n = open(f"Data/SensitivityCoeffs/FM_sensitivity_T_{T_}_case_{case_index}_n.txt","w").write(string_FM_n)
	g_ea = open(f"Data/SensitivityCoeffs/FM_sensitivity_T_{T_}_case_{case_index}_ea.txt","w").write(string_FM_Ea)

	
os.chdir("..")
if "sens_3p_parameters_using_PRS.pkl" not in os.listdir():
	with open('sens_3p_parameters_using_PRS.pkl', 'wb') as file_:
			pickle.dump(selected_PRS, file_)
if "sens_3p_parameters.pkl" not in os.listdir():
	with open('sens_3p_parameters.pkl', 'wb') as file_:
			pickle.dump(selected_BRUTE_FORCE_PARAMETERS, file_)

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
