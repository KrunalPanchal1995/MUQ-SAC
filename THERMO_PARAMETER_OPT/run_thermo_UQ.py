try:
    import ruamel_yaml as yaml
except ImportError:
    from ruamel import yaml
import yaml
import numpy as np
import sys,os
from copy import deepcopy
import pickle
from sklearn.model_selection import train_test_split
# custom modlues .................................................
import ResponseSurface as PRS
import reaction_selection as rs # called first in line 69
import combustion_target_class  # first called in 132
import DesignMatrix as DM           # line 154
import simulation_manager as simulator # line 225
import data_management
import create_parameter_dictionary as create_dict
import Uncertainty as uncertainty
from sample_curve_plotter  import sample_plot
from OptimizationTool import OptimizationTool as Optimizer
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
##  Unloading the target data                 ##
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
##                         ##
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
#raise AssertionError
analysis_type = optInputs['Inputs'].get('AnalysisType', None)

"""
For PLotting purpose

design_matrix = DM.DesignMatrix(unsrt_data,"Tcube",5,10).get_thermo_samples()
#print(design_matrix)

sample_plot(unsrt_data , design_matrix) # generating sample plots in the folder "sample_plots"

"""
"""
def getTotalUnknowns(N):
    n_ = 1 + 5*N + (5*N*(5*N-1))/2
    return int(n_)
    
def getSim(n,design):
    n_ = getTotalUnknowns(n)
    if design == "A-facto":
        sim = int(A_fact_samples)*n_
    elif design == "Tcube":
        sim = 7*n_
    else:
        sim = 7*n_    
    return sim
  

design_matrix = DM.DesignMatrix(unsrt_data,"Tcube",5,len_active_params).get_thermo_samples()
sample_plot(unsrt_data , design_matrix) # generating sample plots in the folder "sample_plots"
raise AssertionError("UQ analysis is Done!!")
"""
len_active_params = 5*len(selected_species)
thermo_design = "Tcube"
#no_sim = getSim(len_active_params,thermo_design)
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

    '''
    if "DesignMatrix_x0_a_fact.csv" not in os.listdir():
        design_matrix_x0 = design_matrix
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
    '''
    def getTotalUnknowns(N):
        n_ = 1 + 2*N + (N*(N-1))/2
        return int(n_)
        
    def getSim(n,design):
        n_ = getTotalUnknowns(n)
        if design == "A-facto":
            sim = int(A_fact_samples)*n_
        elif design == "Tcube":
            sim = 3*n_
        else:
            sim = 7*n_    
        return sim
    
    no_of_sim = {}
    if "DesignMatrix.csv" not in os.listdir():
        no_of_sim_ = getSim(len_active_params,thermo_design)
        print(no_of_sim_)
        design_matrix = DM.DesignMatrix(unsrt_data,thermo_design,no_of_sim_,len_active_params).get_thermo_samples()
        no_of_sim_ = len(design_matrix)
    else:
        no_of_sim_ = getSim(len_active_params,thermo_design)
        design_matrix_file = open("DesignMatrix.csv").readlines()
        design_matrix = []
        no_of_sim_ = len(design_matrix_file)
        for row in design_matrix_file:
            design_matrix.append([float(ele) for ele in row.strip("\n").strip(",").split(",")])
    #raise AssertionError("Design Matrix created!!")
        #sample_plot(unsrt_data , design_matrix)
        #raise AssertionError("Stop!")
        design_matrix_dict = {}
        for case in case_dir:
            design_matrix_dict[case] = design_matrix
            no_of_sim[case] = no_of_sim_
############ ok upto here ##################################
#print(design_matrix)
#print("1st sample is plotted")
#sample_plot(unsrt_data , design_matrix) # generating sample plots in the folder "sample_plots"

###################################################
SSM = simulator.SM(target_list,optInputs,unsrt_data,design_matrix,flag = analysis_type)

if "Perturbed_Mech" not in os.listdir():
    os.mkdir("Perturbed_Mech")
    print("\nPerturbing the Mechanism files\n")

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
activeParameters = []
activeParameters.extend(unsrt_data[species].activeParameters)
for params in activeParameters:
    selected_params.append(1)

selected_params_dict = {}
design_matrix_dict = {}
yaml_loc_dict = {}
for case in case_dir:
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
        file_name = "sim_data_case-"+str(case)+".lst"
        with open(file_name, "r") as file_:
            unique_lines = file_.readlines()        
        # Process the unique lines to extract ETA
        ETA = [float(line.split("\t")[1]) for line in unique_lines]
        temp_sim_opt[str(case)] = ETA
        os.chdir(optDir)
        #print(ETA)
        #raise AssertionError("Generating ETA list for all cases")    
    else:
        os.chdir(optDir)
        #print(os.getcwd())
        os.chdir("case-"+str(case))    
        data_sheet,failed_sim, ETA = data_management.generate_target_value_tables(FlameMaster_Execution_location, target_list, case, fuel,input_=optInputs)
        #print(data_sheet)
        #raise AssertionError("!STOP")
        temp_sim_opt[str(case)] = ETA
        f = open('../Data/Simulations/sim_data_case-'+str(case)+'.lst','w').write(data_sheet)
        g = open('../Data/Simulations/failed_sim_data_case-'+str(case)+'.lst','w').write(failed_sim)
        #f.write(data_sheet)
        #f.close()
        os.chdir(optDir)
#raise AssertionError("only for checking the search iteration of optimization proces.. take the repose_sueface file from -opt/data/responsesurface- to the current directory.")
###############################################
##      Generating the response surface      ##
##                                           ##
###############################################

ResponseSurfaces = {}
selected_PRS = {}
for case_index,case in enumerate(temp_sim_opt):
    yData = np.asarray(temp_sim_opt[case]).flatten()#[0:no_of_sim[case_index]]
    xData = np.asarray(design_matrix_dict[case_index])#[0:no_of_sim[case_index]]
    #print(np.shape(xData))
    #print(np.shape(yData))
    #raise AssertionError("Stop!!")
    
    xTrain,xTest,yTrain,yTest = train_test_split(xData,yData,
                                    random_state=104, 
                                    test_size=0.9, 
                                       shuffle=True)
    #print(np.shape(xTest))
    #print(np.shape(yTrain))
    Response = PRS.ResponseSurface(xTrain,yTrain,case,case_index,prs_type=PRS_type,selected_params=selected_params_dict[case_index])
    Response.create_response_surface()
    if Response.generated == False:
    	Response.test(xTest,yTest)
    	Response.plot(case_index)
    else:
    	Response.test(xTest,yTest)
    	Response.plot(case_index)
    #Response.DoStats_Analysis() #Generates stastical analysis report
    #print(Response.case)
    ResponseSurfaces[case_index] = Response
    #print(Response.selection)
    del xTrain,xTest,yTrain,yTest
raise AssertionError("The Target class, Uncertainty class, Design Matrix and Simulations and Response surface")
os.chdir("..")
##################################################
##        Optimization Procedure                ##
##   Inputs: Traget list and Response Surfaces  ## 
##################################################
if "solution_zeta.save" not in os.listdir():

    opt, opt_zeta,posterior_cov = Optimizer(target_list,opt_method="GD").run_optimization_with_selected_PRS_thermo(unsrt_data,ResponseSurfaces,optInputs)
    
    with open('Z_star.pkl', 'wb') as file_:
        pickle.dump(opt, file_)
    
    with open('V_star.pkl', 'wb') as file_:
        pickle.dump(opt_zeta, file_)
    
    with open('posterior_cov.pkl', 'wb') as file_:
        pickle.dump(posterior_cov, file_)
    
    string_save = ""
    for i, j in enumerate(activeParameters):
        print("\n{}=\t{}".format(activeParameters[i], opt_zeta[i]))
        string_save+="{}=\t{}\n".format(activeParameters[i], opt_zeta[i])
    save = open("solution_zeta.save","w").write(string_save)


    #string_save = ""
    #for i, j in enumerate(activeParameters):
    #    print("\n{}=\t{}".format(activeParameters[i], opt[i]))
    #    string_save+="{}=\t{}\n".format(activeParameters[i], opt[i])
    #save = open("solution.save","w").write(string_save)


    originalMech = Parser(mech_file_location).mech
    copy_of_mech = deepcopy(originalMech)#.deepcopy()
    new_mechanism,a = Manipulator(copy_of_mech,unsrt_data,opt_zeta).doPerturbation()

    #new_mechanism,a,b,c = self.MechManipulator.GeneratePerturbedMechanism(self.target_list[case],self.beta_[i],np.ones(len(selectedParams)),reactionList,self.simulation,"False",extra_arg = self.activeIndexDict)
    string = yaml.dump(new_mechanism,default_flow_style=False)
    f = open("new_mech.yaml","w").write(string)
else:
    
    save = open("solution_zeta.save","r").readlines()
    with open('Z_star.pkl', 'rb') as file_:
        opt = pickle.load(file_)
    
    with open('V_star.pkl', 'rb') as file_:
        opt_zeta = pickle.load(file_)
    
    with open('posterior_cov.pkl', 'rb') as file_:
        posterior_cov = pickle.load(file_)
    opt_zeta = []
    for i in save:
        opt_zeta.append(float(i.split("=")[1].strip()))
    originalMech = Parser(mech_file_location).mech
    copy_of_mech = deepcopy(originalMech)#.deepcopy()
    new_mechanism,a = Manipulator(copy_of_mech,unsrt_data,opt_zeta).doPerturbation()

    #new_mechanism,a,b,c = self.MechManipulator.GeneratePerturbedMechanism(self.target_list[case],self.beta_[i],np.ones(len(selectedParams)),reactionList,self.simulation,"False",extra_arg = self.activeIndexDict)
    string = yaml.dump(new_mechanism,default_flow_style=False)
    f = open("new_mech.yaml","w").write(string)



raise AssertionError("The Target class, Uncertainty class, Design Matrix and Simulations and Response surface")
raise AssertionError()

