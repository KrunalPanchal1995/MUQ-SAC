import os
import os.path
#keywords
DATA = "DATA"
path = os.getcwd();
text_files = [f for f in os.listdir(path) if f.endswith('.csv')]
for j in range (0,len(text_files)):	
	f = open(text_files[j],'r');
	lines = f.readlines()
	print(len(lines));
	print("{} file found\n".format(text_files[j]));
	methanol = [];
	oxygen = [];
	nitrogen = [];
	argon = [];
	equivalence_ratio = [];
	Pressure = [];
	target = [];
	simulation = [];
	temperature = [];
	flame_speed = [];
	std_dvtn = [];
	target_list = [];
	start = [];
	stop_line = 0;
	k = len(lines)
	for line in lines[1:len(lines)-1]:
		data_set = line.split();		
		key_and_value = line.split()
		if len(key_and_value) == 1:
			continue			
		key = key_and_value[0]
		content = key_and_value[2]
		if key == "CH3OH":
			methanol.append(content);
		if key == "O2":
			oxygen.append(content);
		if key == "N2":
			nitrogen.append(content);
		if key == "AR":
			argon.append(content);		
		if key == "H2":
			hydrogen.append(content);
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
		std_dvtn.append(word[2]);
	#print(Pressure,target,Phi,simulation,temperature,ignition_delay,std_dvtn)
	#save_path = '/home/krunal/Downloads/Project work/Objectives/Main Burner (Uncertainity Quantification)/Data sets/Olm_Methanol/Exp data points/Shock_tube_Ignition_delay/DATA_SET_TXT/Generated_target_dataset'
	#name_of_file = input("target_input")
	#completeName = os.path.join(save_path, name_of_file+text_files[j])	
	g = open("target_input_"+text_files[j],"+w");
	for i in range (0,len(equivalence_ratio)):
		g.write("{}, target = {}, simulation = {}, Fuel = {}, O2 = {},N2 = {}, AR = {}, H2 = {}, T = {}, P = {}, phi = {}, observed = {}, deviation = {}#Burke16_{}\n".format(i+1,target[0],simulation[0],methanol[0],oxygen[0],nitrogen[0],argon[0],hydrogen[0],temperature[0],Pressure[0],equivalence_ratio[i],flame_speed[i],std_dvtn[i],text_files[j]));
		#print("i+1, target = "<target[i]<", simulation ="simulation[i]<", T = "<temperature[i]<", P = "<Pressure[i]<"%f, phi = "<Phi[i]<", observed = "<ignition_delay[i]<"#Burke16");
	g.close()
