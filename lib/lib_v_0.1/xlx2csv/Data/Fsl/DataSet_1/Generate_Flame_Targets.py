import os
import os.path
#keywords
DATA = "DATA"
path = os.getcwd();
list_sp = os.listdir("startprofile")
print()
text_files = [f for f in os.listdir(path) if f.endswith('.txt')]
for j in range (0,len(text_files)):	
	f = open(text_files[j],'r');
	lines = f.readlines()
	print(len(lines));
	print("{} file found\n".format(text_files[j]));
	equivalence_ratio = [];
	Pressure = [];
	target = [];
	startprofile = [];
	simulation = [];
	temperature = [];
	flame_speed = [];
	std_dvtn = [];
	target_list = [];
	start = [];
	stop_line = 0
	Dataset = 1
	k = len(lines)
	for line in lines[1:len(lines)-1]:
		data_set = line.split();		
		key_and_value = line.split()
		if len(key_and_value) == 1:
			continue			
		key = key_and_value[0]
		content = key_and_value[2]
		if key == "P":
			Pressure.append(content);
		if key == "Temp":
			temperature.append(content);
		if key == "target":
			target.append(content);
		if key == "simulation":
			simulation.append(content);
		if DATA in data_set:
			start.append(lines.index('Data')+1);
	for data in lines[7:len(lines)]:
		word = data.split(',');
		equivalence_ratio.append(word[0]);
		flame_speed.append(word[1]);
		std_dvtn.append(word[2].split()[0]);
	weight = len(equivalence_ratio)
	
	#print(Pressure,target,Phi,simulation,temperature,ignition_delay,std_dvtn)
	#save_path = '/home/krunal/Downloads/Project work/Objectives/Main Burner (Uncertainity Quantification)/Data sets/Olm_Methanol/Exp data points/Shock_tube_Ignition_delay/DATA_SET_TXT/Generated_target_dataset'
	#name_of_file = input("target_input")
	#completeName = os.path.join(save_path, name_of_file+text_files[j])	
	g = open("target_input_"+text_files[j],"w");
	text = ""
	for i in range (0,len(equivalence_ratio)):
		text += "fsl{}, target = {}, simulation = {}, startProfile = {}, T = {}, P = {}, phi = {}, observed = {}, deviation = {}, DataSet = {}, weight = {} #_{} \n".format(i+1,target[0],simulation[0],list_sp[i],temperature[0],Pressure[0],equivalence_ratio[i],flame_speed[i],std_dvtn[i],Dataset,weight,text_files[j]);
	g.write(text)
	g.close
			#print("i+1, target = "<target[i]<", simulation ="simulation[i]<", T = "<temperature[i]<", P = "<Pressure[i]<"%f, phi = "<Phi[i]<", observed = "<ignition_delay[i]<"#Burke16");
