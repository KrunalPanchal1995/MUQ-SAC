import os
import os.path
import sys

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
for j in range (0,len(input_files)):	
	f = open(input_files[j],'r');
	lines = f.readlines()
	print("{} file found\n".format(input_files[j]));
	fuel = []
	oxidizer = []
	bath_gas1 = []
	bath_gas2 = []
	bath_gas3 = [] 
	temperature = []
	pressure = []
	target = []
	time = []
	simulation = []
	phi = []
	std_dvtn = []
	response_unsrt = []
	measurnment_type = []
	ig_mode = []
	unit = []
	data_structure = []
	tig = []
	fls = []
	scp = []
	time = []
	fuel_value = []
	oxidizer_value = []
	bath_gas1_value = []
	bath_gas2_value = []
	bath_gas3_value = [] 
	stop_line = 0
	start_profile = []
	k = len(lines)
	for line in lines:
		if "!" not in line:
			key_and_value = line.split(":")
			if len(key_and_value) == 1:
				continue	
				
			key = key_and_value[0]
			content = key_and_value[1].split("\n")[0]
			
			if "Target" in key:
				target.append(content)
			if "Fuel"in key:
				fuel.append(content)
			if "Oxidizer" in key:
				oxidizer.append(content)
			if  "Bath_gas1" in key:
				bath_gas1.append(content)
			if "Bath_gas2" in key:
				bath_gas2.append(content)	
			if "Bath_gas3" in key:
				bath_gas3.append(content)		
			if "Pressure" in key:
				pressure.append(content.strip())
				
			if  "Temperature" in key:
				temperature.append(content.strip())
			if "Phi" in key:
				phi.append(content)
			if "Unsrt_factor" in key:
				std_dvtn.append(float(content))
			if "Measurnment_type" in key:
				measurnment_type.append(content)
			if "Ignition_mode" in key:
				ig_mode.append(content)
			if "Data_structure" in key:
				data_structure.append(content)
			if "Unit" in key:
				unit.append(content)
			if "Simulation_type" in key:
				simulation.append(content);
			if "Start_profile" in key:
				start_profile.append(content)
			if "DATA" in key:
				start = (int(lines.index(line)+1));
	# standard units is [mol,Pa,K,ms]	
	
	fuel_species = fuel[0].split("=")[0]
	fuel_conc = fuel[0].split("=")[1]
	
	oxidizer_species = oxidizer[0].split("=")[0]
	oxidizer_conc = oxidizer[0].split("=")[1]
	
	if bath_gas1[0] != "":
		bath_gas1_species = bath_gas1[0].split("=")[0]
		bath_gas1_conc = bath_gas1[0].split("=")[1]
	else:
		bath_gas1_species = "N/N"
		bath_gas1_conc = "N/N"
	if bath_gas2[0] != "":
		bath_gas2_species = bath_gas2[0].split("=")[0]
		bath_gas2_conc = bath_gas2[0].split("=")[1]
	else:
		bath_gas2_species = "N/N"
		bath_gas2_conc = "N/N"
	
	if bath_gas3[0] != "":
		bath_gas3_species = bath_gas3[0].split("=")[0]
		bath_gas3_conc = bath_gas3[0].split("=")[1]
	else:
		bath_gas3_species = "N/N"
		bath_gas3_conc = "N/N"
	
	
	
	Target_units = unit[0].split("\t") 
	if Target_units[0] == "mol":
		cons_factor = 0
	if Target_units[1] == "Pa":
		pressure_factor = 1
	elif Target_units[1] == "kPa":
		 pressure_factor = 1000
	elif Target_units[1] == "MPa":
		 pressure_factor = 10**6
	elif Target_units[1] == "mPa":
		 pressure_factor = 10**(-3)
	elif Target_units[1] == "atm":
		 pressure_factor = 10**(5)
	elif Target_units[1] == "bar":
		 pressure_factor = 10**(5)
	if Target_units[2] == "K":
		temperature_factor = 0
	elif Target_units[2] == "1000/K":
		temperature_factor = 1000
	elif Target_unit[2] == "10000/K":
		temperature_factor = 10000
	elif Target_unit[2] == "1/K":
		temperature_factor = 1
	if Target_units[3] == "ms":
		tig_factor = 0
	elif Target_units[3] == "us":
		tig_factor = 1000
	elif Target_units[3] == "s":
		tig_factor = 0.001 
	if Target_units[4] == "cm/s":
		fls_factor = 0
	
	if pressure[0] !="":
		Target_pressure = float(pressure_factor)*float(pressure[0])
	
	#print(tig_factor)
	Target_data = data_structure[0].split("\t")
	#if Target_data[0].strip() == "T":
		#print("yes")
	Target_columns = len(Target_data)
	#print(Target_columns)

	for i in range(Target_columns):
		for k in lines[start:len(lines)]:
			string = k.split("\n")
			content = string[0].split("\t")	
			if len(content) !=1:
				if Target_data[i].strip() == "F":
					print("entered")
					fuel_value.append(content[i])
				if Target_data[i].strip() == "Ox":
					oxidizer_value.append(content[i])
				if Target_data[i].strip() == "BG1":
					bath_gas1_value.append(content[i])
				if Target_data[i].strip() == "BG2":
					bath_gas2_value.append(content[i])
				if Target_data[i].strip() == "BG3":
					bath_gas3_value.append(content[i])
				if Target_data[i].strip() == "P":
					pressure.append(pressure_factor*float(content[i]))
					
				if Target_data[i].strip() == "T":
					#print(float(content[i]))
					#print(temperature_factor)
					if temperature_factor == 0:
						t = float(content[i])
						temperature.append(float(t))
					elif temperature_factor == 1:
						t = float(content[i])
						temperature.append(float(1/t))
					elif temperature_factor == 1000:
						t = float(content[i])
						temperature.append(float(1000/t))
					elif temperature_factor == 10000:
						t = float(content[i])
						temperature.append(float(10000/t))	
					#print(temperature)
				if Target_data[i].strip() == "Time":
					time.append(content[i])
				if Target_data[i].strip() == "Phi":
					phi.append(content[i].strip())
				if Target_data[i].strip() == "Tig":
					#print(float(content[i]))
					if tig_factor == 0:
						tig.append(float(content[i]))
					else:
						tig.append(tig_factor*(float(content[i])))
					#print(tig)
				if Target_data[i].strip() == "Fls":
					fls.append(content[i])
				if Target_data[i].strip() == "Scp":
					scp.append(content[i])
		
		
	tig = filter_list(tig)
	scp = filter_list(scp)
	fls = filter_list(fls)
	temperature = filter_list(temperature)
	#print(temperature)
	pressure = filter_list(pressure)
	phi = filter_list(phi)
	time = filter_list(time)
	#print(Target_pressure)
	count = 1
	if len(phi) == 0:
		Target_phi = "N/N"
	else:
		Target_phi = float(phi[0])
	
	if target[0].strip() == "Tig":
		for i in range(len(tig)):
			response_unsrt.append(std_dvtn[0]*float(tig[i]))
		g = open("targets.opt","+a");
		for i in range (len(tig)):
			if len(fuel_value)==0 and len(pressure)==1:
				g.write("{},{}, target : {}, simulation : {}, measurnment_type : {},Ignition_mode : {}, Fuel : {}, Oxidizer : {},BG1 : {}, BG2 : {},BG3 : {}, T : {}, P : {}, Phi : {}, observed : {}, deviation : {},data_weight = {}\n".format(count,input_files[j],target[0],simulation[0],measurnment_type[0],ig_mode[0],fuel[0],oxidizer[0],bath_gas1[0],bath_gas2[0],bath_gas3[0],temperature[i],Target_pressure,Target_phi,tig[i],response_unsrt[i],len(tig)));
			
			elif len(fuel_value)==0 and len(pressure)>1:
				g.write("{},{}, target : {}, simulation : {}, measurnment_type : {},Ignition_mode : {}, Fuel : {}, Oxidizer : {},BG1 : {}, BG2 : {},BG3 : {}, T : {}, P : {}, Phi : {}, observed : {}, deviation : {},data_weight = {}\n".format(count,input_files[j],target[0],simulation[0],measurnment_type[0],ig_mode[0],fuel[0],oxidizer[0],bath_gas1[0],bath_gas2[0],bath_gas3[0],temperature[i],pressure[i],Target_phi,tig[i],response_unsrt[i],len(tig)));
			
			elif len(fuel_value)!=0 and len(pressure)==1:
				g.write("{},{}, target : {}, simulation : {}, measurnment_type : {},Ignition_mode : {}, Fuel : x->{}={}, Oxidizer : x->{}={},BG1 : x->{}={}, BG2 : x->{}={},BG3 : x->{}={}, T : {}, P : {}, Phi : {}, observed : {}, deviation : {},data_weight = {}\n".format(count,input_files[j],target[0],simulation[0],measurnment_type[0],ig_mode[0],fuel_species,fuel_conc,oxidizer_species,oxidizer_cons,bath_gas1_species,bath_gas1_cons,bath_gas2_species,bath_gas2_conc,bath_gas3_species,bath_gas3_conc,temperature[i],Target_pressure,Target_phi,tig[i],response_unsrt[i],len(tig)));
			
			elif len(fuel_value)!=0 and len(pressure)!=1:
				g.write("{},{}, target : {}, simulation : {}, measurnment_type : {},Ignition_mode : {}, Fuel : x->{}={}, Oxidizer : x->{}={},BG1 : x->{}={}, BG2 : x->{}={},BG3 : x->{}={}, T : {}, P : {}, Phi : {}, observed : {}, deviation : {},data_weight = {}\n".format(count,input_files[j],target[0],simulation[0],measurnment_type[0],ig_mode[0],fuel_species,fuel_conc,oxidizer_species,oxidizer_cons,bath_gas1_species,bath_gas1_cons,bath_gas2_species,bath_gas2_conc,bath_gas3_species,bath_gas3_conc,temperature[i],Target_pressure[i],Target_phi,tig[i],response_unsrt[i],len(tig)));
		
			#print(count)
			count+=1
			
	elif target[0].strip() == "Scp":
		for i in range(len(scp)):
			response_unsrt.append(std_dvtn[0]*float(i))
	elif target[0].strip() == "Fls":
		for i in range(len(fls)):
			response_unsrt.append(std_dvtn[0]*float(i))

	#units = 
	#print("units")
	#for data in lines[10:len(lines)]:
	#	word = data.split();
	#	temperature.append(word[0]);
	#	ignition_delay.append(word[1]);
	#	std_dvtn.append(word[2]);
	#print(Pressure,target,Phi,simulation,temperature,ignition_delay,std_dvtn)
	
	#save_path = '/home/krunal/Downloads/Project work/Objectives/Main Burner (Uncertainity Quantification)/Data sets/Olm_Methanol/Exp data points/Shock_tube_Ignition_delay/DATA_SET_TXT/Generated_target_dataset'
	#name_of_file = input("target_input")
	#completeName = os.path.join(save_path, name_of_file+text_files[j])	
	
	#g = open("target_input_"+text_files[j],"+w");
		
	# All data are converted into a standard format
	
	#for i in range (0,len(temperature)):
	#	g.write("{}, target = {}, simulation = {}, measurnment_type = {}, Fuel = {}, O2 = {},BG1 = {}, BG2 = {},BG3 = {}, T = {}, P = {}, phi = {}, observed = {}, deviation = {}\n".format(input_files[j],target[0],simulation[0],methanol[0],oxygen[0],nitrogen[0],argon[0],temperature[i],Pressure[0],Phi[0],ignition_delay[i],std_dvtn[i]));
		#print("i+1, target = "<target[i]<", simulation ="simulation[i]<", T = "<temperature[i]<", P = "<Pressure[i]<"%f, phi = "<Phi[i]<", observed = "<ignition_delay[i]<"#Burke16");

