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
import statistics
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


ZetaFile = open(optInputs["Locations"]["zeta_file"],"r").readlines()
x_data = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in ZetaFile])
y_data = []
if "DesignMatrix.csv" not in os.listdir():
	design_matrix = DM.DesignMatrix(unsrt_data,design_type,2000,len(manipulationDict["activeParameters"])).getSamples()
else:
	design_matrix_file = open("DesignMatrix.csv").readlines()
	design_matrix = []
	for row in design_matrix_file:
		design_matrix.append([float(ele) for ele in row.strip("\n").strip(",").split(",")])
sim_type = "opt"
simulator_func = simulator.SM(target_list,optInputs,unsrt_data,design_matrix)

iter_number = 0
objective = 1
cases = {}
l = []
for ind,i in enumerate(target_list):
	l.append(ind)
	cases[ind] = target_list[ind]
if "yData.txt" not in os.listdir():
	#y_data = simulator.getSimulatedValues(x_data)
	for x in x_data:
		y_data.append(simulator_func.getSimulatedValues(x,l,cases,iter_number,objective))
		iter_number+=1
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
#raise AssertionError("Ydata is generated!!")
"""
Now test the respones surface 
"""
for z in range(targets_count):
	X_train = optInputs["Locations"]["xData"]+"/Beta_list_case-"+str(z)+".csv"

	Y_train = optInputs["Locations"]["yData"]+"/sim_data_case-"+str(z)+".lst"
	
	X_test = optInputs["Locations"]["xTest"]+"/Beta_list_case-"+str(z)+".csv"

	Y_test = optInputs["Locations"]["yTest"]+"/sim_data_case-"+str(z)+".lst"
	
	file_x = open(X_train,"r").readlines()
	file_y = open(Y_train,"r").readlines()

	x_train_data = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in file_x])
	y_train_data = (np.asarray([float(i.strip("\n").split()[1]) for i in file_y]))
	
	#y_train_data = np.log(np.asarray([float(i.strip("\n").split()[1]) for i in file_y]))
	
	
	file_x_test = open(X_test,"r").readlines()
	file_y_test = open(Y_test,"r").readlines()

	x_test_data = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in file_x_test])
	y_test_data = ((np.asarray([float(i.strip("\n").split()[1]) for i in file_y_test])))
	
	y_data_test = ((y_dict[z]))
	
	#y_test_data = np.log(np.asarray([float(i.strip("\n").split()[1]) for i in file_y_test]))
	
	#y_data_test = np.log(y_dict[z])

	x_data_test = []
	for row in x_data: 
		count = 0
		temp = []
		for i in unsrt_data:
			temp.append(row[count])
			temp.append(row[count+1])
			temp.append(row[count+2])
			count+=3
		x_data_test.append(np.asarray(temp))


	#print(np.shape(x_train_data))
	#print(np.shape(x_train_data))
	#print(np.shape(x_data_test))
	#print(np.shape(y_data_test))

	#print(len(x_train_data[0]))
	#print(len(x_data_test[0]))
	model = PRS.ResponseSurface(x_train_data,y_train_data,z,2)
	
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
	train_error_general = [abs(model.resFramWrk[i_i]-i) for i_i,i in enumerate(y_train_data)]
	max_error_general = max(error_general)
	train_max_error_general = max(train_error_general)
	
	#print(error_general)
	mean_error_general = statistics.mean(error_general)
	
	"""
	Testing the response surface using the zetas generated from optimization algorithm
	"""
	general_yfit_test = []
	for i in x_data_test:
		general_yfit_test.append((model.evaluate(i)))
	#print(general_yfit_test)
	#print(len(general_yfit_test))
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
		general_yfit.append((model.evaluate(i)))
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
	plt.plot(np.asarray(general_yfit_test),np.asarray(y_data_test),"r^",ms=4,label=f"Search Iterations (max error = {max_error_test_general:.3f}%)")#, \nmean error = {mean_error_test_general:.3f}%)")
	
	plt.scatter(np.asarray(general_yfit), np.asarray(y_test_data),s=10,marker="o",color="black", edgecolor="black",label=f"Testing dataset (max error = {max_error_test:.3f}%)")#, \nmean error = {mean_error_test:.3f}%)")
	
	#plt.scatter(np.asarray(model.resFramWrk), np.asarray(y_train_data),s=8,marker="s",color="none", edgecolor="green",label=f"Training dataset (max error = {max_error_general:.3f}%, mean error = {mean_error_general:.3f}%)")
	
	
	#plt.plot(np.asarray(general_yfit_test),np.asarray(y_data_test),"r^",ms=4,label=f"Search Iterations (max error = {max_error_residual_test_general:.3f}cm/s)")#, \nmean error = {mean_error_residual_test_general:.3f}cm/s)")
	
	#plt.scatter(np.asarray(general_yfit), np.asarray(y_test_data),s=10,marker="o",color="black", edgecolor="black",label=f"Testing dataset (max error = {max_error_residual_test:.3f}cm/s)")#, \nmean error = {mean_error_residual_test:.3f}cm/s)")
			
	#plt.scatter(np.asarray(model.resFramWrk), np.asarray(y_train_data),s=8,marker="s",color="none", edgecolor="green",label=f"Training dataset (max error = {train_max_error_general:.3f}cm/s")#, mean error = {mean_error_general:.3f}%)")
	

	plt.xlabel("Black box simulations")
	plt.ylabel("Response surface predictions")
	#plt.plot(y_data,model.svr_yfit,"g.",label="SVR")
	#plt.plot(y_train_data,model.resFramWrk,"r.",label="Traning")
	#plt.plot(y_data_test,svr_yfit_test,"g.",label ="SVR")
	#plt.plot(y_data_test,general_yfit_test,"g.",label="Testing")
	plt.legend()
	#plt.show()
	#pk.Sample(rxnUnsrt_data,reaction_index,sample_length).getSamples()
	plt.savefig("Testing"+str(z)+".pdf",bbox_inches="tight")
	
	
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
	"""
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
	"""



