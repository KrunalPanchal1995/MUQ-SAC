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
import Uncertainty as uncertainty
from sample_curve_plotter  import sample_plot

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
iFile = str(os.getcwd())+"/"+str(sys.argv[1])
dataCounts = optInputs[count]
binLoc = optInputs["Bin"]
inputs = optInputs["Inputs"]
locations = optInputs["Locations"]
startProfile_location = optInputs[startProfile]
stats_ = optInputs["Stats"]
global A_fact_samples

A_fact_samples = stats_["Sampling_of_PRS"]
if "sensitive_parameters" not in stats_:
	stats_["sensitive_parameters"] = "zeta"
	optInputs["Stats"]["sensitive_parameters"] = "zeta"
if "Arrhenius_Selection_Type" not in stats_:
	stats_["Arrhenius_Selection_Type"] = "all"
	optInputs["Stats"]["Arrhenius_Selection_Type"] = "all"
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

print("\nParallel threads are {}".format(parallel_threads))
targetLines = open(locations[targets],'r').readlines()
addendum = yaml.safe_load(open(locations[add],'r').read())

#########################################
###    Reading the mechanism file    ####
#########################################

MECH = locations["mechanism"]
carbon_number = optInputs["SA"]["carbon_number"]

with open(MECH,'r') as file_:
	yaml_mech = file_.read()

"""
#mechanism = yaml.safe_load(yaml_mech)
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
"""
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
selected_species = []
#parameter_dict,selected_species,selected_reactions = create_dict.dictionary_creator(analysis_type, mechanism,carbon_number,rxn_list,species_list)
for species in unsrt_data:
	selected_species.append(species)
	#print(species ,dir(unsrt_data[species]) )

print(f"Total species selected:  {len(selected_species)}\n")

analysis_type = optInputs['Inputs'].get('AnalysisType', None)


design_matrix = DM.DesignMatrix(unsrt_data,"Tcube",10,10).get_thermo_samples()
#print(design_matrix)

sample_plot(unsrt_data , design_matrix) # generating sample plots in the folder "sample_plots"





raise AssertionError("UQ analysis is Done!!")
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
            


