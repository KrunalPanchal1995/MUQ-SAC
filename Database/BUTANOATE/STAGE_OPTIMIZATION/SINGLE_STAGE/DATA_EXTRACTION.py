import numpy as np
import os,sys, pandas as pd
import matplotlib.pyplot as plt
import reaction_selection as rs
import yaml
#yaml = yaml.YAML(typ="safe",pure=True)
#Read all the files and extract the reaction index as well as unsrt data
if len(sys.argv) > 1:
	sens_files_location = sys.argv[1]
	rxn_data = pd.read_csv(sys.argv[2])
	mechanism = yaml.safe_load(open(sys.argv[3],"r").read())
	df = pd.DataFrame(rxn_data)
	rxn_id_list = [i for i in np.sort(df["Rxn"].to_numpy())]
	Unsrt = np.log(df["Unsrt"].to_numpy())
	data_type = df["data_type"].to_list()
	#	optInputs = yaml.safe_load(input_file)
	print("Sensitivity files found\n")
else:
	print("Please enter a valid input folder dir as arguement. \n Program exiting")
	exit()

#print(rxn_id_list)
rxn_unsrt_data = {}
for i,rxn_id in enumerate(rxn_id_list):
	rxn_unsrt_data[str(rxn_id)] = f"{Unsrt[i]},{data_type[i]}"
start = os.getcwd()
#data_type = "constant;end_points"
os.chdir(sens_files_location)
sens_dir = os.getcwd()
list_files = os.listdir()
file_name = []
Temperature = []
selected_reactions = {}
selected_id = {}
sensitivity_dict = {}
for rxn_id in rxn_id_list:
	sensitivity_dict[str(rxn_id)] = {}
for file_ in list_files:
	file_name.append(file_)
	Temperature.append(float(file_.strip(".txt").split("_")[3]))
	T = float(file_.strip(".txt").split("_")[3])
	file_data = open(file_,"r").readlines()
	for rxn_id in rxn_id_list:
		#print(type(rxn_id))
		for line in file_data[1:]:
			line_split = line.strip(" \t").split("\t")
			if rxn_id == int(line_split[1]):
				#print(rxn_id)
				selected_reactions[line_split[1]] = line_split[2].strip("\n")
				selected_id[line_split[2].strip("\n")] = line_split[1]
				#sensitivity_dict[line_split[2]] = {}
				sensitivity_dict[line_split[1]][T] = float(line_split[0])
							
				break
os.chdir("..")
#print(sensitivity_dict['2266'])
#raise AssertionError("Stop!1")
###
Temperature = np.sort(np.asarray(Temperature))
rxn_list = [selected_reactions[key].strip("\n") for key in selected_reactions]
#print(len(rxn_list))
RXN_ID_LIST = [key.strip("\n") for key in selected_reactions]
rxn_type = rs.getRxnType_unsrt(mechanism,rxn_list)
rxn_count = {}
for rxn in rxn_type:
	if rxn_type[rxn] == "PLOG":
		rxn_count[rxn] = rs.getRxnCount(mechanism,rxn)
	elif rxn_type[rxn] == "PLOG-Duplicate":
		rxn_count[rxn] = rs.getRxnCount(mechanism,rxn)
#raise AssertionError("Stop!")
###Plot sensitivity data


fig = plt.figure()
for rxn_id in sensitivity_dict:
	#temperature_data = []
	sensitivity_data = []
	for T in Temperature:
		#temperature_data.append(float(data))
		#print(rxn_id,T)
		sensitivity_data.append(float(sensitivity_dict[rxn_id][T]))
	plt.plot(Temperature,sensitivity_data,"-",label=f"{rxn_id}")
plt.plot(Temperature,np.repeat(0.05,len(Temperature)),"--",label="Cut-off")
plt.plot(Temperature,np.repeat(-0.05,len(Temperature)),"--",label="Cut-off")
plt.legend()
plt.ylabel("Sensitivity Coeffifient")
plt.xlabel("Temperature (K)")
plt.savefig("SA.pdf",bbox_inches="tight")

string_g = ""
#print(len(rxn_type),len(RXN_ID_LIST))
for i,rxn in enumerate(rxn_type):
	string_g+=f"{rxn}\t{RXN_ID_LIST[i]}\t{rxn_type[rxn]}\n"#{RXN_ID_LIST[i]}
g = open("rxn.data","w").write(string_g)
#raise AssertionError("Stop!!")

def getElementaryRxn(rxn,rxn_id,data):
	#data_type = 
	data_type = data.strip("\n").split(",")[1]
	un = data.strip("\n").split(",")[0]
	if data_type == "constant;end_points":
		unsrt = f"{un},{un}"
		temp = "300,2500"
	rxn = rxn.replace("=","&#61;")
	rxn = rxn.replace("<","&#60;")
	rxn = rxn.replace(">","&#62;")
	rxn = rxn.replace(" ","&#032;")
	string = f"""\t<reaction rxn = '{rxn}' no = '{rxn_id}'>
	 <class>Elementary</class>
		<type>pressure_independent</type>
		<sub_type name = "forward">
			<multiple> False </multiple>
			<branches>N/A</branches>
			<pressure_limit>N/A</pressure_limit>
			<common_temp>N/A</common_temp>
			<temp_limit>N/A</temp_limit>
		</sub_type>
		<perturbation_type>all</perturbation_type>
		<data_type>{data_type}</data_type>
		<unsrt>{unsrt}</unsrt>
		<temp>{temp}</temp>
		<file></file>	
	</reaction>\n"""
	return string
def getDuplicateRxn(rxn,rxn_id,data):
	data_type = data.strip("\n").split(",")[1]
	un = data.strip("\n").split(",")[0]
	if data_type == "constant;end_points":
		unsrt = f"{un},{un}"
		temp = "300,2500"
	rxn = rxn.replace("=","&#61;")
	rxn = rxn.replace("<","&#60;")
	rxn = rxn.replace(">","&#62;")
	rxn = rxn.replace(" ","&#032;")
	string = f"""\t<reaction rxn = '{rxn}' no = '{rxn_id}a'>
	 <class>Duplicate</class>
		<type>pressure_independent</type>
		<sub_type name = "duplicate">
			<multiple> True </multiple>
			<branches>A</branches>
			<pressure_limit>N/A</pressure_limit>
			<common_temp>N/A</common_temp>
			<temp_limit>N/A</temp_limit>
		</sub_type>
		<perturbation_type>all</perturbation_type>
		<data_type>{data_type}</data_type>
		<unsrt>{unsrt}</unsrt>
		<temp>{temp}</temp>
		<file></file>	
	</reaction>\n\t<reaction rxn = '{rxn}' no = '{rxn_id}b'>
	 <class>Duplicate</class>
		<type>pressure_independent</type>
		<sub_type name = "duplicate">
			<multiple> True </multiple>
			<branches>B</branches>
			<pressure_limit>N/A</pressure_limit>
			<common_temp>N/A</common_temp>
			<temp_limit>N/A</temp_limit>
		</sub_type>
		<perturbation_type>all</perturbation_type>
		<data_type>{data_type}</data_type>
		<unsrt>{unsrt}</unsrt>
		<temp>{temp}</temp>
		<file></file>	
	</reaction>\n"""
	return string
def getTROERxn():
	string = f""""""
	return string
def getFallOffRxn():
	string = f""""""
	return string
def getPLOGRxn(rxn,rxn_id,data,rxn_count):
	rxn_count = int(rxn_count)-2
	data_type = data.strip("\n").split(",")[1]
	un = data.strip("\n").split(",")[0]
	if data_type == "constant;end_points":
		unsrt = f"{un},{un}"
		temp = "300,2500"
	rxn = rxn.replace("=","&#61;")
	rxn = rxn.replace("<","&#60;")
	rxn = rxn.replace(">","&#62;")
	rxn = rxn.replace(" ","&#032;")
	string = f"""	<PLOG rxn ='{rxn}' no = '{rxn_id}:High'>
	 <class>PLOG</class>
		<type>pressure_dependent</type>
		<sub_type name = "forward">
			<multiple> True </multiple>
			<branches></branches>
			<pressure_limit>High</pressure_limit>
			<common_temp>N/A</common_temp>
			<temp_limit>N/A</temp_limit>
		</sub_type>
		<perturbation_type>all</perturbation_type>
		<data_type>{data_type}</data_type>
		<unsrt>{unsrt}</unsrt>
		<temp>{temp}</temp>
		<file></file>	
	</PLOG>
	<PLOG rxn = '{rxn}' no = '{rxn_id}:Low'>
	 <class>PLOG</class>
		<type>pressure_dependent</type>
		<sub_type name = "forward">
			<multiple> True </multiple>
			<branches></branches>
			<pressure_limit>Low</pressure_limit>
			<common_temp>N/A</common_temp>
			<temp_limit>N/A</temp_limit>
		</sub_type>
		<perturbation_type>all</perturbation_type>
		<data_type>{data_type}</data_type>
		<unsrt>{unsrt}</unsrt>
		<temp>{temp}</temp>
		<file></file>	
	</PLOG>
	
	<PLOG-Interconnectedness>
	  <class>PLOG</class>
	  <InterConnectedRxns>{rxn_id}:High,{rxn_id}:Low</InterConnectedRxns>
	  <RxnCount>{rxn_count}</RxnCount>
	  <unsrtQuantType>Interpolation</unsrtQuantType>
	</PLOG-Interconnectedness>	\n"""
	return string

def getPLOG_DuplicateRxn(rxn,rxn_id,data,rxn_count):
	rxn_count = int(rxn_count)-2
	data_type = data.strip("\n").split(",")[1]
	un = data.strip("\n").split(",")[0]
	if data_type == "constant;end_points":
		unsrt = f"{un},{un}"
		temp = "300,2500"
	rxn = rxn.replace("=","&#61;")
	rxn = rxn.replace("<","&#60;")
	rxn = rxn.replace(">","&#62;")
	rxn = rxn.replace(" ","&#032;")
	string = f"""	<PLOG rxn = '{rxn}' no = '{rxn_id}a:High'>
	 <class>PLOG-Duplicate</class>
		<type>pressure_dependent</type>
		<sub_type name = "forward">
			<multiple>True</multiple>
			<branches>A</branches>
			<pressure_limit>High</pressure_limit>
			<common_temp>N/A</common_temp>
			<temp_limit>N/A</temp_limit>
		</sub_type>
		<perturbation_type>all</perturbation_type>
		<data_type>{data_type}</data_type>
		<unsrt>{unsrt}</unsrt>
		<temp>{temp}</temp>
		<file></file>	
	</PLOG>
	<PLOG rxn = '{rxn}' no = '{rxn_id}a:Low'>
	 <class>PLOG-Duplicate</class>
		<type>pressure_dependent</type>
		<sub_type name = "forward">
			<multiple>True</multiple>
			<branches>A</branches>
			<pressure_limit>Low</pressure_limit>
			<common_temp>N/A</common_temp>
			<temp_limit>N/A</temp_limit>
		</sub_type>
		<perturbation_type>all</perturbation_type>
		<data_type>{data_type}</data_type>
		<unsrt>{unsrt}</unsrt>
		<temp>{temp}</temp>
		<file></file>	
	</PLOG>
	<PLOG-Interconnectedness>
	  <class>PLOG-Duplicate</class>
	  <InterConnectedRxns>{rxn_id}a:High,{rxn_id}a:Low</InterConnectedRxns>
	  <RxnCount>{rxn_count}</RxnCount>
	  <unsrtQuantType>Interpolation</unsrtQuantType>
	</PLOG-Interconnectedness>
	<PLOG rxn = '{rxn}' no = '{rxn_id}b:High'>
	 <class>PLOG-Duplicate</class>
		<type>pressure_dependent</type>
		<sub_type name = "forward">
			<multiple>True</multiple>
			<branches>B</branches>
			<pressure_limit>High</pressure_limit>
			<common_temp>N/A</common_temp>
			<temp_limit>N/A</temp_limit>
		</sub_type>
		<perturbation_type>all</perturbation_type>
		<data_type>{data_type}</data_type>
		<unsrt>{unsrt}</unsrt>
		<temp>{temp}</temp>
		<file></file>	
	</PLOG>
	<PLOG rxn = 'rxn' no = '{rxn_id}b:Low'>
	 <class>PLOG-Duplicate</class>
		<type>pressure_dependent</type>
		<sub_type name = "forward">
			<multiple>True</multiple>
			<branches>B</branches>
			<pressure_limit>Low</pressure_limit>
			<common_temp>N/A</common_temp>
			<temp_limit>N/A</temp_limit>
		</sub_type>
		<perturbation_type>all</perturbation_type>
		<data_type>{data_type}</data_type>
		<unsrt>{unsrt}</unsrt>
		<temp>{temp}</temp>
		<file></file>	
	</PLOG>
	<PLOG-Interconnectedness>
	  <class>PLOG-Duplicate</class>
	  <InterConnectedRxns>{rxn_id}b:High,{rxn_id}b:Low</InterConnectedRxns>
	  <RxnCount>{rxn_count}</RxnCount>
	  <unsrtQuantType>Interpolation</unsrtQuantType>
	</PLOG-Interconnectedness>\n"""
	return string

string = '<?xml version="1.0" encoding="UTF-8"?>\n<uncertainty>\n'
for rxn in rxn_list:
	if rxn_type[rxn] == "Elementary":
		string+=getElementaryRxn(rxn,selected_id[rxn],rxn_unsrt_data[selected_id[rxn]])
	
	elif rxn_type[rxn] == "Duplicate":
		string+=getElementaryRxn(rxn,selected_id[rxn],rxn_unsrt_data[selected_id[rxn]])
	
	elif rxn_type[rxn] == "PLOG":
		string+=getPLOGRxn(rxn,selected_id[rxn],rxn_unsrt_data[selected_id[rxn]],rxn_count[rxn])
	
	elif rxn_type[rxn] == "PLOG-Duplicate":
		#print(selected_id[rxn])
		#print(rxn_unsrt_data[selected_id[rxn]])
		#print(rxn_count[rxn])
		string+=getPLOG_DuplicateRxn(rxn,selected_id[rxn],rxn_unsrt_data[selected_id[rxn]],rxn_count[rxn])
		
	elif rxn_type[rxn] == "ThirdBody":
		string+=getElementaryRxn(rxn,selected_id[rxn],rxn_unsrt_data[selected_id[rxn]])


string+='</uncertainty>'
g = open("new_unsrt.xml","w").write(string)
raise AssertionError("Stop!")	
#1] Tasks 
#	- Create plots for the sensitivity data
#	- Get reaction type
#	- Create xml file for uncertainty
				
