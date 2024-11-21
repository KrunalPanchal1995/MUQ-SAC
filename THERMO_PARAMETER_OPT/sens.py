try:
    import ruamel_yaml as yaml
except ImportError:
    from ruamel import yaml
import yaml
import numpy as np
import sys,os
from copy import deepcopy
import pickle

# custom modlues .................................................
import reaction_selection as rs # called first in line 69
import combustion_target_class  # first called in 132
import DesignMatrix as DM           # line 154
import simulation_manager as simulator # line 225
import data_management
import create_parameter_dictionary as create_dict



#########################################
###    Reading the input file        ####
#########################################
print(len(sys.argv))
if len(sys.argv) > 2:
	input_file = open(sys.argv[1],'r')
	optInputs = yaml.safe_load(input_file)
	rxn_list =[i.strip("\n") for i in open(sys.argv[2],"r").readlines()]
	species_list = [i.strip("\n") for i in open(sys.argv[3],"r").readlines()]
	#print(rxn_list)
	print("\n\t########################\n\tInput file , List of reactions, List of Species are found\n\t########################\n")
	#raise AssertionError("Stop")
elif len(sys.argv) > 1:
	input_file = open(sys.argv[1],'r')
	optInputs = yaml.safe_load(input_file)
	rxn_list = []
	species_list = []
	print("\n\t########################\n\tInput file found\n\t########################\n")
else:
	print("Please enter a valid input file name as arguement. \n Two arguments can be passed:\n\t1. Traget opt file\n\t2. List of reactions\nThe code will still work by passing only the first argument\n\nProgram exiting")
	exit()

#print("list of reactions ::::",rxn_list)

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
species = mechanism['phases'][0]["species"] # .................takes all the species availabel in mechanism file....................
species_data = mechanism["species"] # ................ takes data corresponding to each species like..........composition, thermo,transport.................
#print(species_data)
#raise AssertionError
reactions = mechanism["reactions"] #....... takes all the reactions present in mechanism file ....#3984 in case of butonate- MB2D.yaml file.....................
######...........we can choose thermo data here instead of reactions for sensitivity analysis of thermo data .........................#########


#print(reaction_dict)
#print(selected_species)
#raise AssertionError("stop")


analysis_type = optInputs['Inputs'].get('AnalysisType', None)
print("Analysis Type is :" ,analysis_type)
parameter_dict,selected_species,selected_reactions = create_dict.dictionary_creator(analysis_type, mechanism,carbon_number,rxn_list,species_list)
print(f"Total species selected:  {len(selected_species)}\n")
#print(parameter_dict)
#raise AssertionError
#raise AssertionError("STOP")
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
	
#########################################
###    Creating Design Matrix for    ####
###    sensitivity analysis          ####
#########################################
if analysis_type == "reaction":
    """
    For sensitivity analysis of reactions we create two design matrices:
        - For one, we multiply all reactions by a factor of 2
        - For the second, we divide all reactions by a factor of 0.5
    """
    if "DesignMatrix_x0_a_fact.csv" not in os.listdir():
        design_matrix_x0 = DM.DesignMatrix(selected_reactions, design_type, len(parameter_dict)).getNominal_samples(flag=analysis_type)
        s = ""
        for row in design_matrix_x0:
            for element in row:
                s += f"{element},"
            s += "\n"
        with open('DesignMatrix_x0_a_fact.csv', 'w') as ff:
            ff.write(s)
    else:
        design_matrix_file = open("DesignMatrix_x0_a_fact.csv").readlines()
        design_matrix_x0 = []
        for row in design_matrix_file:
            design_matrix_x0.append([float(ele) for ele in row.strip("\n").strip(",").split(",")])

    if "DesignMatrix_x2_a_fact.csv" not in os.listdir():
        design_matrix_x2 = DM.DesignMatrix(selected_reactions, design_type, len(parameter_dict)).getSA_samples(flag=analysis_type)
        s = ""
        for row in design_matrix_x2:
            for element in row:
                s += f"{element},"
            s += "\n"
        with open('DesignMatrix_x2_a_fact.csv', 'w') as ff:
            ff.write(s)
    else:
        design_matrix_file = open("DesignMatrix_x2_a_fact.csv").readlines()
        design_matrix_x2 = []
        for row in design_matrix_file:
            design_matrix_x2.append([float(ele) for ele in row.strip("\n").strip(",").split(",")])

elif analysis_type == "thermo":
    
    if "DesignMatrix_x0_a_fact.csv" not in os.listdir():
        design_matrix_x0 = DM.DesignMatrix(selected_species, design_type, len(parameter_dict)).getNominal_samples(flag=analysis_type)
        s = ""
        for row in design_matrix_x0:
            for element in row:
                s += f"{element},"
            s += "\n"
        with open('DesignMatrix_x0_a_fact.csv', 'w') as ff:
            ff.write(s)
    else:
        design_matrix_file = open("DesignMatrix_x0_a_fact.csv").readlines()
        design_matrix_x0 = []
        for row in design_matrix_file:
            design_matrix_x0.append([float(ele) for ele in row.strip("\n").strip(",").split(",")])

    if "DesignMatrix_2x.csv" not in os.listdir():
        design_matrix_x2 = DM.DesignMatrix(selected_species, design_type, len(parameter_dict)).getSA_samples(flag=analysis_type)
        s = ""
        for row in design_matrix_x2:
            for element in row:
                s += f"{element},"
            s += "\n"
        with open('DesignMatrix_2x.csv', 'w') as ff:
            ff.write(s)
    else:
        design_matrix_file = open("DesignMatrix_2x.csv").readlines()
        design_matrix_x2 = []
        for row in design_matrix_file:
            design_matrix_x2.append([float(ele) for ele in row.strip("\n").strip(",").split(",")])

#########################################
###    Creating Simulation Field     ####
#########################################
print("\t\t#########################################\n\t\t###    Creating Simulation Field     ####\n\t\t#########################################")
if "SA" not in os.listdir():
	list_dir = ["multiply","divide","nominal","Data","Data/Simulations","Data/Simulations/Multiply","Data/Simulations/Divide","Data/Simulations/Nominal","Data/ResponseSurface"]
	os.mkdir("SA")
	os.chdir("SA")
	for i in list_dir:
		if i not in os.listdir():
			os.mkdir(i)
	os.chdir("multiply")
	SADir = os.getcwd()
else:
	os.chdir("SA")
	list_dir = ["multiply","divide","nominal","Data","Data/Simulations","Data/Simulations/Multiply","Data/Simulations/Divide","Data/Simulations/Nominal","Data/ResponseSurface"]
	for i in list_dir:
		if os.path.exists(i):
			continue
		else:
			os.mkdir(i)
	os.chdir("multiply")
	SADir = os.getcwd()

#####################################################
#### Multiplying the reactions by a factor of 2   ###
#####################################################
print("\n\t\t#####################################################\n\t\t#### Multiplying the reactions by a factor of 2   ###\n\t\t#####################################################")

#raise AssertionError
if os.path.isfile("progress") == False:
	FlameMaster_Execution_location_x2 = simulator.SM(target_list,optInputs,parameter_dict,design_matrix_x2,flag=analysis_type).make_dir_in_parallel()
	
else:
	print("\t\tProgress file detected")
	progress = open(SADir+"/progress",'r').readlines()
	FlameMaster_Execution_location_x2 = []
	with open(SADir+"/locations") as infile:
		for line in infile:
			FlameMaster_Execution_location_x2.append(line)
#raise AssertionError
#############################
#### Nominal simulations  ###
#############################

os.chdir("..")
os.chdir("nominal")
SADir = os.getcwd()


print("\n\t\t#############################\n\t\t#### Nominal simulations ###\n\t\t#############################")

if os.path.isfile("progress") == False:
	FlameMaster_Execution_location_x0 = simulator.SM(target_list,optInputs,parameter_dict,design_matrix_x0,flag=analysis_type).make_dir_in_parallel()
	
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
		folderName_x0 = [i.split("\t")[0] for i in open("sim_data_case-"+str(case)+".lst").readlines()]
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


##################################
##### From Multiply Folder #######
##################################
temp_sim_opt_x2 = {}
for case in case_dir:	
	os.chdir("Data/Simulations/Multiply")
	if "sim_data_case-"+str(case)+".lst" in os.listdir():
		ETA_x2 = [float(i.split("\t")[1]) for i in open("sim_data_case-"+str(case)+".lst").readlines()]
		folderName_x2 = [i.split("\t")[0] for i in open("sim_data_case-"+str(case)+".lst").readlines()]
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
		#print(ETA_x2)
		temp_sim_opt_x2[str(case)] = {}
		temp_sim_opt_x2[str(case)]["ETA"] = ETA_x2
		temp_sim_opt_x2[str(case)]["index"] = index_x2
		f = open('../../Data/Simulations/Multiply/sim_data_case-'+str(case)+'.lst','w').write(data_sheet_x2)
		g = open('../../Data/Simulations/Multiply/failed_sim_data_case-'+str(case)+'.lst','w').write(failed_sim_x2)
		#f.write(data_sheet)
		#f.close()
		os.chdir(SAdir)

selected_PRS = {}
for case in case_dir:
	T = int(target_list[case].temperature)
	index = temp_sim_opt_x2[str(case)]["index"]
	multiply = np.asarray(temp_sim_opt_x2[str(case)]["ETA"])
	nominal = np.asarray(temp_sim_opt_x0[str(case)]["ETA"])
	case_FM = (multiply - nominal)/nominal
	#print(reaction_dict,index)
	sens_FM = {}
	param_Sa = {}
if(analysis_type =="reaction"):
	for i,ind in enumerate(index):
		param_Sa[parameter_dict[int(ind)]] = case_FM[i]
		sens_FM[ind] = case_FM[i]
	selected_PRS[str(case)] = param_Sa
	sensitivity_FM = dict(sorted(sens_FM.items(), key=lambda item: abs(item[1]),reverse = True))
	
	string_FM = "Sensitivity Analysis (using Cantera) Tig:\n"
	for ind,sens in enumerate(sensitivity_FM):
		sens = str(sens)
		string_FM +=f"\t{sensitivity_FM[sens]:.8f}\t{int(int(sens))}\t{reaction_dict[int(sens)]}\n"

	g = open(f"Data/FM_sensitivity_T_{T}_case_{case}.txt","w").write(string_FM)

	if "sens_parameters.pkl" not in os.listdir(".."):
		with open('../sens_parameters.pkl', 'wb') as file_:
			pickle.dump(selected_PRS, file_)	
	print("\n\t################################\n\tSENSITIVITY ANALYSIS DONE!!\n\t################################\n\t")	

elif(analysis_type =="thermo"):
	for i,ind in enumerate(index):
		param_Sa[(ind)] = case_FM[i]
		sens_FM[ind] = case_FM[i]
	selected_PRS[str(case)] = param_Sa
	sensitivity_FM = dict(sorted(sens_FM.items(), key=lambda item: abs(item[1]),reverse = True))
	
	string_FM = "Sensitivity Analysis (using Cantera) Tig:\n"
	for ind,sens in enumerate(sensitivity_FM):
		sens = str(sens)
		string_FM +=f"\t{sensitivity_FM[sens]:.8f}\t{sens}\n"

	g = open(f"Data/FM_case_{case}.txt","w").write(string_FM)

	if "sens_parameters.pkl" not in os.listdir(".."):
		with open('../sens_parameters.pkl', 'wb') as file_:
			pickle.dump(selected_PRS, file_)	
	print("\n\t################################\n\tSENSITIVITY ANALYSIS DONE!!\n\t################################\n\t")	



	
