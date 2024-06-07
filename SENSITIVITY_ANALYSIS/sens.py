try:
    import ruamel_yaml as yaml
except ImportError:
    from ruamel import yaml

import numpy as np
import sys,os
import reaction_selection as rs
import simulation_manager as simulator
import DesignMatrix as DM
from copy import deepcopy
import combustion_target_class
import data_management

#########################################
###    Reading the input file        ####
#########################################

if len(sys.argv) > 1:
	input_file = open(sys.argv[1],'r')
	optInputs = yaml.safe_load(input_file)
	print("\n\t########################\n\tInput file found\n\t########################\n")
else:
	print("Please enter a valid input file name as arguement. \n For details of preparing the input file, please see the UserManual\n\nProgram exiting")
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
		reaction_dict[index] = rxn["equation"]


rxn_type = rs.getRxnType(mechanism,selected_reactions)
string_f = ""
string_g = ""
for index in reaction_dict:
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
##  Unloading the target data	  		          ##
## TARGET CLASS CONTAINING EACH TARGET AS A	CASE  ##
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

#########################################
###    Creating Design Matrix for    ####
###    sensitivity analysis          ####
#########################################

"""
For sensitivity analysis we create two design matrix
	- for one, we multiply all reactions by a factor of 2
	- for second, we devide all reactions by a factor of 0.5
"""
if "DesignMatrix_x0.csv" not in os.listdir():
	design_matrix_x0 = DM.DesignMatrix(selected_reactions,design_type,len(reaction_dict)).getNominal_samples()
else:
	design_matrix_file = open("DesignMatrix_x0.csv").readlines()
	design_matrix_x0 = []
	for row in design_matrix_file:
		design_matrix_x0.append([float(ele) for ele in row.strip("\n").strip(",").split(",")])
		
if "DesignMatrix_x2.csv" not in os.listdir():
	design_matrix_x2 = DM.DesignMatrix(selected_reactions,design_type,len(reaction_dict)).getSA_samples(np.log(2.0))
else:
	design_matrix_file = open("DesignMatrix_x2.csv").readlines()
	design_matrix_x2 = []
	for row in design_matrix_file:
		design_matrix_x2.append([float(ele) for ele in row.strip("\n").strip(",").split(",")])

#if "DesignMatrix_x0_5.csv" not in os.listdir():#
#	design_matrix_x0_5 = DM.DesignMatrix(selected_reactions,design_type,len(selected_reactions)).getSA_samples(np.log(0.5))
#else:
#	design_matrix_file = open("DesignMatrix_x0_5.csv").readlines()
#	design_matrix_x0_5 = []
#	for row in design_matrix_file:
#		design_matrix_x0_5.append([float(ele) for ele in row.strip("\n").strip(",").split(",")])

#########################################
###    Creating Simulation Field     ####
#########################################
print("\t\t#########################################\n\t\t###    Creating Simulation Field     ####\n\t\t#########################################")

if "SA" not in os.listdir():
	os.mkdir("SA")
	os.chdir("SA")
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
	os.chdir("SA")
	os.chdir("multiply")
	SADir = os.getcwd()

#####################################################
#### Multiplying the reactions by a factor of 2   ###
#####################################################
print("\n\t\t#####################################################\n\t\t#### Multiplying the reactions by a factor of 2   ###\n\t\t#####################################################")


if os.path.isfile("progress") == False:
	FlameMaster_Execution_location_x2 = simulator.SM(target_list,optInputs,rxn_dict,design_matrix_x2).make_dir_in_parallel()
	
else:
	print("\t\tProgress file detected")
	progress = open(SADir+"/progress",'r').readlines()
	FlameMaster_Execution_location_x2 = []
	with open(SADir+"/locations") as infile:
		for line in infile:
			FlameMaster_Execution_location_x2.append(line)

#os.chdir("..")
#os.chdir("divide")
#SADir = os.getcwd()	
#####################################################
#### Multiplying the reactions by a factor of 0.5 ###
#####################################################
#print("\n\t\t#####################################################\n\t\t#### Multiplying the reactions by a factor of 0.5   ###\n\t\t#####################################################")
#raise AssertionError("Making directories completed !!!")
#if os.path.isfile("progress") == False:
#	FlameMaster_Execution_location_x0_5 = simulator.SM(target_list,optInputs,rxn_dict,design_matrix_x0_5).make_dir_in_parallel()
	
#else:
#	print("\t\tProgress file detected")
#	progress = open(SADir+"/progress",'r').readlines()
#	FlameMaster_Execution_location_x0_5 = []
#	with open(SADir+"/locations") as infile:
#		for line in infile:
#			FlameMaster_Execution_location_x0_5.append(line)

os.chdir("..")
os.chdir("nominal")
SADir = os.getcwd()

#############################
#### Nominal simulations  ###
#############################
print("\n\t\t#############################\n\t\t#### Nominal simulations ###\n\t\t#############################")

if os.path.isfile("progress") == False:
	FlameMaster_Execution_location_x0 = simulator.SM(target_list,optInputs,rxn_dict,design_matrix_x0).make_dir_in_parallel()
	
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
		data_sheet_x0,failed_sim_x0,index_x0, ETA_x0 = data_management.generate_SA_target_value_tables(FlameMaster_Execution_location_x0, target_list, case, fuel)
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
		data_sheet_x2,failed_sim_x2,index_x2, ETA_x2 = data_management.generate_SA_target_value_tables(FlameMaster_Execution_location_x2, target_list, case, fuel)
		#print(data_sheet)
		#raise AssertionError("!STOP")
		temp_sim_opt_x2[str(case)] = {}
		temp_sim_opt_x2[str(case)]["ETA"] = ETA_x2
		temp_sim_opt_x2[str(case)]["index"] = index_x2
		f = open('../../Data/Simulations/Multiply/sim_data_case-'+str(case)+'.lst','w').write(data_sheet_x2)
		g = open('../../Data/Simulations/Multiply/failed_sim_data_case-'+str(case)+'.lst','w').write(failed_sim_x2)
		#f.write(data_sheet)
		#f.close()
		os.chdir(SAdir)

##################################
##### From Divide Folder   #######
##################################
#temp_sim_opt_x0_5 = {}
#for case in case_dir:	
#	os.chdir("Data/Simulations/Divide")
#	if "sim_data_case-"+str(case)+".lst" in os.listdir():
#		ETA_x0_5 = [float(i.split("\t")[1]) for i in open("sim_data_case-"+str(case)+".lst").readlines()]
#		folderName_x0_5 = [float(i.split("\t")[0]) for i in open("sim_data_case-"+str(case)+".lst").readlines()]
#		temp_sim_opt_x0_5[str(case)] = {}
#		temp_sim_opt_x0_5[str(case)]["ETA"] = ETA_x0_5
#		temp_sim_opt_x0_5[str(case)]["index"] = folderName_x0_5
#		os.chdir(SAdir)
		#print(ETA)
		#raise AssertionError("Generating ETA list for all cases")	
#	else:
#		os.chdir(SAdir)
#		os.chdir("divide/case-"+str(case))	
#		data_sheet_x0_5,failed_sim_x0_5, index_x0_5, ETA_x0_5 = data_management.generate_SA_target_value_tables(FlameMaster_Execution_location_x0_5, target_list, case, fuel)
		#print(data_sheet)
		#raise AssertionError("!STOP")
#		temp_sim_opt_x0_5[str(case)] = {}
#		temp_sim_opt_x0_5[str(case)]["ETA"] = ETA_x0_5
#		temp_sim_opt_x0_5[str(case)]["index"] = index_x0_5
#		f = open('../../Data/Simulations/Divide/sim_data_case-'+str(case)+'.lst','w').write(data_sheet_x0_5)
#		g = open('../../Data/Simulations/Divide/failed_sim_data_case-'+str(case)+'.lst','w').write(failed_sim_x0_5)
#		#f.write(data_sheet)
		#f.close()
#		os.chdir(SAdir)

#print("\n\t##############################\n\tSensitivity Analysis Done\n\tCalculating sensitivity coefficients\n\t##############################\n")
for case in case_dir:
	T = target_list[case].temperature
	index = temp_sim_opt_x2[str(case)]["index"]
	#index_0 = temp_sim_opt_x0[str(case)]["index"]
	multiply = np.asarray(temp_sim_opt_x2[str(case)]["ETA"])
	#divide = np.asarray(temp_sim_opt_x0_5[str(case)]["ETA"])
	nominal = np.asarray(temp_sim_opt_x0[str(case)]["ETA"])
	#case_sens = np.log(multiply/divide)/np.log(4)
	#case_sens_0 = np.log(multiply/nominal)/np.log(2)
	case_FM = (multiply - nominal)/nominal
	##### Using new formula: S = ln(Tig_2/Tig_0.5)/ln(k_2/k_0.5)
	#sens = {}
	#for i,ind in enumerate(index):
	#	sens[ind] = case_sens[i]

	#sensitivity = dict(sorted(sens.items(), key=lambda item: abs(item[1]),reverse = True))

	#string = "Sensitivity Analysis (using Cantera) Tig:\n"
	#for ind,sens in enumerate(sensitivity):
	#	string +=f"\t{sensitivity[sens]:,.8f}\t{str(int(sens+1))}\t{reaction_dict[sens]}\n"

	#g = open(f"Data/sensitivity_case_{case}.txt","w").write(string)
	##### Using new formula: S = ln(Tig_2/Tig_0)/ln(k_2/k_0)
	#sens_0 = {}
	#for i,ind in enumerate(index):
	#	sens_0[ind] = case_sens_0[i]

	#sensitivity_0 = dict(sorted(sens_0.items(), key=lambda item: abs(item[1]),reverse = True))

	#string_0 = "Sensitivity Analysis (using Cantera) Tig:\n"
	#for ind,sens in enumerate(sensitivity_0):
	#	string_0 +=f"\t{sensitivity_0[sens]:,.8f}\t{str(int(sens+1))}\t{reaction_dict[sens]}\n"

	#g = open(f"Data/Nom_sensitivity_case_{case}.txt","w").write(string_0)
	##### Using new formula: S = (Tig_2 - Tig_0)/Tig_0
	
	sens_FM = {}
	for i,ind in enumerate(index):
		sens_FM[ind] = case_FM[i]

	sensitivity_FM = dict(sorted(sens_FM.items(), key=lambda item: abs(item[1]),reverse = True))

	string_FM = "Sensitivity Analysis (using Cantera) Tig:\n"
	for ind,sens in enumerate(sensitivity_FM):
		string_FM +=f"\t{sensitivity_FM[sens]:.8f}\t{str(int(sens+1))}\t{reaction_dict[sens]}\n"

	g = open(f"Data/FM_sensitivity_T_{T}.txt","w").write(string_FM)
	
print("\n\t################################\n\tSENSITIVITY ANALYSIS DONE!!\n\t################################\n\t")	
