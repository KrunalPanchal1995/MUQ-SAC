import os
import os.path
#keywords
DATA = "DATA"
path = os.getcwd();
text_files = [f for f in os.listdir(path) if f.endswith('.txt')]
count = 1
target_string =""
for j in range (0,len(text_files)):	
	f = open(text_files[j],'r');
	lines = f.readlines()
	#print(len(lines));
	print("{} file found\n".format(text_files[j]));
	methanol = [];
	argon = [];
	time = [];
	oxygen = [];
	nitrogen = [];
	phi = [];
	speciesCONC = []
	temperature = [];
	Pressure = [];
	target = [];
	simulation = [];
	Phi = [];
	ignition_delay = [];
	std_dvtn = [];
	target_list = [];
	start = [];
	tend = [];
	species = [];
	dataset = []
	stop_line = 0;
	k = len(lines)
	for line in lines[1:13]:
		data_set = line.split();		
		key_and_value = line.split()
		if len(key_and_value) == 1:
			continue			
		key = key_and_value[0]
		content = key_and_value[2]
		if key == "CH3OH":
			methanol.append(content);
		if key == "AR":
			argon.append(content)
		if key == "O2":
			oxygen.append(content)
		if key == "N2":
			nitrogen.append(content);		
		if key == "P":
			Pressure.append(content)
		if key == "species":
			species.append(content)
		if key == "Phi":
			phi.append(content)
		if key == "target":
			target.append(content);
		if key == "Simulation":
			simulation.append(content);
		if key =="T":
			temperature.append(content)
		if key == "Tend":
			tend.append(content)
		if key =="Dataset":
			dataset.append(content)
		if DATA in data_set:
			start.append(lines.index('Data')+1);
		
	for data in lines[14:len(lines)]:
		if len(data.split(",")) == 1:
			continue
		else:
			word = data.split(",");
			time.append(word[0]);
			speciesCONC.append(word[1].split()[0]);
			std_dvtn.append(0.13*float(word[1].split()[0]));
	weight = len(time)
	#print(Pressure,target,Phi,simulation,temperature,ignition_delay,std_dvtn)
	#save_path = '/home/krunal/Downloads/Project work/Objectives/Main Burner (Uncertainity Quantification)/Data sets/Olm_Methanol/Exp data points/Shock_tube_Ignition_delay/DATA_SET_TXT/Generated_target_dataset'
	#name_of_file = input("target_input")
	#completeName = os.path.join(save_path, name_of_file+text_files[j])	
	g = open("target_input_"+text_files[j],"+w");
	temp = ""
	for i in range (0,len(time)):
		if i<len(time):			
			temp += "scp{}, target = {}, simulation = {},species = {}, Fuel = {}, Tend = {}, AR = {}, O2 = {}, N2 = {}, Phi = {}, T = {}, P = {}, time = {}, observed = {}, deviation = {}, DataSet = {}, weight = {}  #2013_{}\n".format(count,target[0],simulation[0],species[0],methanol[0],tend[0],argon[0],oxygen[0],nitrogen[0],phi[0],temperature[0],Pressure[0],time[i],speciesCONC[i],std_dvtn[i],dataset[0],weight,text_files[j]);
			count+=1	
	g.write(temp)
	target_string+=temp
	g.close()

target = open("target_input.opt","w")
target.write(target_string)
target.close()
		#print("i+1, target = "<target[i]<", simulation ="simulation[i]<", T = "<temperature[i]<", P = "<Pressure[i]<"%f, phi = "<Phi[i]<", observed = "<ignition_delay[i]<"#Burke16");
