import os.path
import re, os, sys
import subprocess
import numpy as np
import shutil
import yaml
import pandas as pd
from itertools import chain

#keywords
if len(sys.argv) > 1:
	os.chdir(sys.argv[1])
	print("Input directory found\n")
else:
	print("Please enter a valid input directory as arguement. \n For details of preparing the input file, please see the UserManual\n\nProgram exiting")
	
def filter_list(List):
	temp = []
	for i in List:
		if i != "":
			temp.append(i)
	return temp

DATA = "DATA"
path = os.getcwd();
input_files = [f for f in os.listdir(path) if f.endswith('.csv')]
input_headers = [f for f in os.listdir(path) if f.endswith('.header')]
if "output" not in os.listdir():
	os.mkdir("output")

for i in input_files:
	csv_file = i
	if str(i.split(".")[0])+".header" in input_headers:
		header_file = input_headers[input_headers.index(str(i.split(".")[0])+".header")]
	else:
		continue
	file_head = open(header_file,"r").read()
	file_data = open(csv_file,"r").read()
	head = yaml.safe_load(file_head)
	data = pd.read_csv(i,sep=',', comment='#')
	fuel_type = len(data[head["Data_structure"][0]])*[head["Fuel_type"]]
	oxidizer = len(data[head["Data_structure"][0]])*[head["Oxidizer"]]
	bath_gas = head["Bath_gas"]
	target = head["Target"]
	#print(head)
	if "phi" in head["Data_structure"]:
		phi = data["phi"]
	elif "Phi" in head:
		phi = np.ones(len(data[head["Data_structure"][0]]))*float(head["Phi"])
	else:
		phi =  len(data[head["Data_structure"][0]])*["N/A"]
	if "res" in head["Data_structure"]:
		res = data["res"]
	
	if "Pi" in head["Data_structure"]:
		pressure_i = data["Pi"]
	
	if "P" in head["Data_structure"]:
		pressure = data["P"]
	else:
		pressure = np.ones(len(data[head["Data_structure"][0]]))*float(head["Pressure"])
		
	if "Ti" in head["Data_structure"]:
		temperature_i = data["Ti"]
		
	if "T" in head["Data_structure"]:
		temperature = data["T"]
	else:
		temperature = np.ones(len(data[head["Data_structure"][0]]))*float(head["Temperature"])
	
	if "Sr" in head["Data_structure"]:
		serial_id = data["Sr"]
		dataSetIndex = []
		for i in serial_id:
			string = head["DataSet"]+"_"+str(i)
			dataSetIndex.append(string)
	else:
		dataSetIndex = len(data[head["Data_structure"][0]])*head["DataSet"]
	
	if "Multi" in fuel_type[0]:
		fuel = len(data[head["Data_structure"][0]])*[head["Fuel"]]
		if "a" in head["Data_structure"]:
			fuel_x = []
			for i in range(len(data[head["Data_structure"][0]])):
				temp = {}
				for j in fuel[0]:
					temp[j] = data[j][i]
				fuel_x.append(temp)
		else:
			fuel_x = np.ones(len(data[head["Data_structure"][0]]))*float(head["Fuel_x"])
	else:
		fuel = len(data[head["Data_structure"][0]])*[head["Fuel"]]
		if "Fuel" in head["Data_structure"]:
			fuel_x = data["Fuel"]
		else:
			fuel_x = np.ones(len(data[head["Data_structure"][0]]))*float(head["Fuel_x"])
		
	if "Ox" in head["Data_structure"]:
		oxidizer_x = data["Ox"]
	else:
		oxidizer_x = np.ones(len(data[head["Data_structure"][0]]))*float(head["Oxidizer_x"])
		
	if "bg1" in head["Data_structure"]:
		bath_gas_x = []
		for i in range(len(data[head["Data_structure"][0]])):
			temp = {}
			for j in bath_gas:
				temp[j] = data[j][i]
			bath_gas_x.append(temp)	
	else:
		bath_gas_x = len(data[head["Data_structure"][0]])*[head["Bath_gas_x"]]
	
	if target in head["Data_structure"]:
		target_value = data[target]
	else:
		target_value = np.ones(len(data[head["Data_structure"][0]]))*float(head[target])
		
	if "Sim" in head["Data_structure"]:
		sim_type = data["Sim"]
	else:
		sim_type = len(data[head["Data_structure"][0]])*[head["Simulation_type"]]
	
	if "Measure_type" in head["Data_structure"]:
		Measure_type = data["Measure_type"]
	else:
		Measure_type = len(data[head["Data_structure"][0]])*[head["Measurnment_type"]]
		
	
	if "flame_type" in head["Data_structure"]:
		flame_type = data["flame_type"]
	elif "Flame_type" in head:
		flame_type = len(data[head["Data_structure"][0]])*[head["Flame_type"]]
	else:
		flame_type =  len(data[head["Data_structure"][0]])*["N/A"]
		
	if "reactor_type" in head["Data_structure"]:
		reactor_type = data["reactor_type"]
	elif "Reactor_type" in head:
		reactor_type = len(data[head["Data_structure"][0]])*[head["Reactor_type"]]
	else:
		reactor_type =  len(data[head["Data_structure"][0]])*["N/A"]
	
	if "ign_mode" in head["Data_structure"]:
		ign_type = data["ign_mode"]
	elif "Ignition_mode" in head:
		ign_type = len(data[head["Data_structure"][0]])*[head["Ignition_mode"]]
	else:
		ign_type =  len(data[head["Data_structure"][0]])*["N/A"]
	
	if "start_pro" in head["Data_structure"]:
		start_pro = data["start_pro"]
	elif "startprofile" in head:
		start_pro = len(data[head["Data_structure"][0]])*[head["startprofile"]]
	else:
		start_pro =  len(data[head["Data_structure"][0]])*["N/A"]
		
	if "flow_rate" in head["Data_structure"]:
		flow_rate = data["flow_rate"]
	elif "Flow_rate" in head:
		flow_rate = np.ones(len(data[head["Data_structure"][0]]))*float(head["Flow_rate"])
	else:
		flow_rate =  len(data[head["Data_structure"][0]])*["N/A"]
			
	if "sigma" in head["Data_structure"]:
		unsrt = []
		unsrt_type = head["Unsrt"]['type']
		sigma = data["sigma"]
		if unsrt_type == "absolute":
			for i in range(len(data[head["Data_structure"][0]])):
				temp = float(sigma[i])
				unsrt.append(temp)
		else:
			for i in range(len(data[head["Data_structure"][0]])):
				temp = float(target_value[i])*float(sigma[i])
				unsrt.append(temp)
	else:
		unsrt = len(data[head["Data_structure"][0]])*[head["Unsrt"]]
	
	if "obs_unit" in head["Data_structure"]:
		units = len(data[head["Data_structure"][0]])*[head["Unsrt"]]
		obs_unit = data["obs_unit"]
	else:
		units = len(data[head["Data_structure"][0]])*[head["Unit"]]
		obs_unit = len(data[head["Data_structure"][0]])*[""]
	string = ""
	weight = len(data[head["Data_structure"][0]])
	bath_gas = weight*[head["Bath_gas"]]
	target = weight*[target]
	for i in range(len(dataSetIndex)):
		if "Tig" in target:
			string+="{}\t|{}\t| target -- {}\t| simulation -- {}\t| measurnment_type -- {}\t|Ignition_mode -- {}\t|Flame_type -- {}\t|Reactor_type -- {}\t|Fuel_type -- {}\t| Fuel -- x->{} = {}\t| Oxidizer -- x->{} = {}\t| Bath_gas -- x->{} = {}\t|T -- {}\t| P -- {}\t| flow_rate -- {}\t|Phi -- {}\t| observed -- {}\t|obs_unit -- {}\t| deviation -- {}\t|data_weight -- {}\t|units -- {}\n".format(i+1,dataSetIndex[i],target[i],sim_type[i],Measure_type[i],ign_type[i],flame_type[i],reactor_type[i],fuel_type[i],fuel[i],fuel_x[i],oxidizer[i],oxidizer_x[i],bath_gas[i],bath_gas_x[i],temperature[i],pressure[i],flow_rate[i],phi[i],target_value[i],obs_unit[i],unsrt[i],weight,units[i])
		elif "RCM" in target:
			string+="{}\t|{}\t| target -- {}\t| simulation -- {}\t| measurnment_type -- {}\t|Ignition_mode -- {}\t|Flame_type -- {}\t|Reactor_type -- {}\t|Fuel_type -- {}\t| Fuel -- x->{} = {}\t| Oxidizer -- x->{} = {}\t| Bath_gas -- x->{} = {}\t| Ti -- {}\t| T -- {}\t| Pi -- {}\t| P -- {}\t| flow_rate -- {}\t|Phi -- {}\t| observed -- {}\t|obs_unit -- {}\t| deviation -- {}\t|data_weight -- {}\t|units -- {}\n".format(i+1,dataSetIndex[i],target[i],sim_type[i],Measure_type[i],ign_type[i],flame_type[i],reactor_type[i],fuel_type[i],fuel[i],fuel_x[i],oxidizer[i],oxidizer_x[i],bath_gas[i],bath_gas_x[i],temperature_i[i],temperature[i],pressure_i[i],pressure[i],flow_rate[i],phi[i],target_value[i],obs_unit[i],unsrt[i],weight,units[i])
		elif "Fls" in target:
			string+="{}\t|{}\t| target -- {}\t| simulation -- {}\t| measurnment_type -- {}\t|Ignition_mode -- {}\t|Flame_type -- {}\t|Reactor_type -- {}\t|Fuel_type -- {}\t| Fuel -- x->{} = {}\t| Oxidizer -- x->{} = {}\t| Bath_gas -- x->{} = {}\t|T -- {}\t| P -- {}\t| flow_rate -- {}\t|Phi -- {}\t| observed -- {}\t|obs_unit -- {}\t| deviation -- {}\t|data_weight -- {}\t|units -- {}\n".format(i+1,dataSetIndex[i],target[i],sim_type[i],Measure_type[i],ign_type[i],flame_type[i],reactor_type[i],fuel_type[i],fuel[i],fuel_x[i],oxidizer[i],oxidizer_x[i],bath_gas[i],bath_gas_x[i],temperature[i],pressure[i],flow_rate[i],phi[i],target_value[i],obs_unit[i],unsrt[i],weight,units[i])
		else:
			string+="{}\t|{}\t| target -- {}\t| simulation -- {}\t| measurnment_type -- {}\t|Ignition_mode -- {}\t|Flame_type -- {}\t|Reactor_type -- {}\t|Fuel_type -- {}\t| Fuel -- x->{} = {}\t| Oxidizer -- x->{} = {}\t| Bath_gas -- x->{} = {}\t|T -- {}\t| P -- {}\t| flow_rate -- {}\t|Phi -- {}\t| res -- {}\t| observed -- {}\t|obs_unit -- {}\t| deviation -- {}\t|data_weight -- {}\t|units -- {}\n".format(i+1,dataSetIndex[i],target[i],sim_type[i],Measure_type[i],ign_type[i],flame_type[i],reactor_type[i],fuel_type[i],fuel[i],fuel_x[i],oxidizer[i],oxidizer_x[i],bath_gas[i],bath_gas_x[i],temperature[i],pressure[i],flow_rate[i],phi[i],res[i],target_value[i],obs_unit[i],unsrt[i],weight,units[i])
	os.chdir("output")
	outFile = open(csv_file.split(".")[0]+'.out','+w')
	outFile.write(string)
	outFile.close()
	os.chdir("..")
	
	
	
	
	
