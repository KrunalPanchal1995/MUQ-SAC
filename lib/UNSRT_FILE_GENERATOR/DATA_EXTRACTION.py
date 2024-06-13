import numpy as np
import os,sys, pandas as pd
import matplotlib.pyplot as plt
import reaction_selection as rs

#Read all the files and extract the reaction index as well as unsrt data
if len(sys.argv) > 1:
	sens_files_location = sys.argv[1]
	rxn_data = pd.read_csv(sys.argv[2])
	df = pd.DataFrame(rxn_data)
	rxn_list = df["Rxn"].to_list()
	Unsrt = np.log(df["Unsrt"].to_numpy())
	#	optInputs = yaml.safe_load(input_file)
	print("Sensitivity files found\n")
else:
	print("Please enter a valid input folder dir as arguement. \n Program exiting")
	exit()

start = os.getcwd()

os.chdir(sens_files_location)
sens_dir = os.getcwd()
list_files = os.listdir()
file_name = []
Temperature = []
selected_reactions = []
sensitivity_dict = {}
for file_ in list_files:
	file_name.append(file_)
	Temperature.append(float(file_.strip(".txt").split("_")[3]))
	T = float(file_.strip(".txt").split("_")[3])
	file_data = open(file_,"r").readlines()
	for rxn_id in rxn_list:
		#print(type(rxn_id))
		for line in file_data[1:]:
			line_split = line.strip(" \t").split("\t")
			#print(line_split)
			if rxn_id == int(line_split[1]):
				#print(rxn_id)
				selected_reactions.append(line_split[2])
				sensitivity_dict[line_split[2]] = {}
				sensitivity_dict[line_split[2]][T] = float(line_split[1][0])
				break
	
#print(sort(file_name,Temperature))
#print(list(set(selected_reactions)))
raise AssertionError("Stop!")
file_name = ""
for i in os.listdir():
	if 'FinalSorted' in i:
		file_name = i

file_data = open(file_name,"r").readlines()

top_50 = []
for i in file_data:
	if len(top_50) < 50:
		if "=" in i:
			top_50.append(i)

print(top_50)

write_data = ""
for i in top_50:
	temp = i.split("\t")
	write_data+=f"{temp[3].strip()}\t{temp[2]}\n"

file_write = open("top_50.csv","w").write(write_data)
				
