import numpy as np
import os, sys
from copy import deepcopy
import yaml
import data_management
import combustion_target_class
import simulation_manager2_0 as simulator
from ruamel.yaml import YAML
# Load input file
if len(sys.argv) > 1:
    with open(sys.argv[1], 'r') as input_file:
        optInputs = yaml.safe_load(input_file)
        print("Input file found\n")
else:
    print("Please provide a valid input file as argument.")
    exit()

home_dir = os.getcwd()
dataCounts = optInputs["Counts"]
inputs = optInputs["Inputs"]
locations = optInputs["Locations"]
stats_ = optInputs["Stats"]

mech_file_location = locations["mechanism"]
fuel = inputs["fuel"]
parallel_threads = dataCounts["parallel_threads"]
targets_count = int(dataCounts["targets_count"])
targetLines = open(locations["targets"], 'r').readlines()
addendum = yaml.safe_load(open(locations["addendum"], 'r'))

# Parse targets
target_list, string_target = [], ""
for c_index, target in enumerate(targetLines[:targets_count]):
    if "#" in target:
        target = target[:target.index('#')]
    add = deepcopy(addendum)
    t = combustion_target_class.combustion_target(target, add, c_index)
    string_target += f"{t.dataSet_id}|{t.target}|{t.species_dict}|{t.temperature}|{t.pressure}|{t.phi}|{t.observed}|{t.std_dvtn}\n"
    target_list.append(t)

with open("target_data.txt", "w") as f:
    f.write(string_target)

# Prepare NOMINAL simulation directory
if "NOMINAL" not in os.listdir():
    os.mkdir("NOMINAL")
os.chdir("NOMINAL")
nominalDir = os.getcwd()
os.makedirs("Data/Simulations", exist_ok=True)

# Run nominal simulations
unsrt_data = {}
design_matrix = []
progress_file = "progress"
if not os.path.isfile(progress_file):
    FlameMaster_Execution_location = simulator.SM(target_list, optInputs, unsrt_data, design_matrix).make_nominal_dir_in_parallel()
else:
    with open("locations") as infile:
        FlameMaster_Execution_location = [line.strip() for line in infile]

# Extract simulation results
temp_sim_opt = {}
case_dir = range(len(target_list))

for case in case_dir:
    os.chdir("Data/Simulations")
    sim_file = f"sim_data_case-{case}.lst"
    if sim_file in os.listdir():
        ETA = [np.exp(float(i.split("\t")[1])) / 10 for i in open(sim_file).readlines()]
        temp_sim_opt[str(case)] = ETA
    else:
        os.chdir(nominalDir)
        os.chdir(f"case-{case}")
        data_sheet, failed_sim, ETA = data_management.generate_target_value_tables(FlameMaster_Execution_location, target_list, case, fuel)
        temp_sim_opt[str(case)] = np.exp(ETA) / 10
        with open(f'../Data/Simulations/{sim_file}', 'w') as f:
            f.write(data_sheet)
        with open(f'../Data/Simulations/failed_sim_data_case-{case}.lst', 'w') as f:
            f.write(failed_sim)
    os.chdir(nominalDir)

# Save target data CSVs
dataset = set([t.dataSet_id for t in target_list])
for d_set in dataset:
    string_1 = "DS_ID,T,Obs(us),Nominal\n"
    string_2 = "DS_ID,T,P,Phi,Fuel,Ox,BathGas,Obs(us),Nominal\n"
    folder = None
    for target in target_list:
        if target.dataSet_id != d_set:
            continue
        if target.target in ["Tig", "RCM", "JSR"]:
            folder = target.target
            string_1 += f"{target.uniqueID},{target.temperature},{target.observed},{temp_sim_opt[str(target.index)][0]}\n"
        elif target.target == "Fls":
            folder = "Fls"
            string_2 += f"{target.uniqueID},{target.temperature},{target.pressure},{target.phi},{target.fuel_dict},{target.oxidizer_x},{target.BG_dict},{target.observed},{temp_sim_opt[str(target.index)][0]}\n"

    out_dir = f"../Plot/Dataset/{folder}/"
    os.makedirs(out_dir, exist_ok=True)
    file_path = os.path.join(out_dir, f"{d_set}.csv")
    with open(file_path, "w") as f:
        f.write(string_1 if folder != "Fls" else string_2)
os.chdir(home_dir)
print("Nominal simulations done.. return 0")
