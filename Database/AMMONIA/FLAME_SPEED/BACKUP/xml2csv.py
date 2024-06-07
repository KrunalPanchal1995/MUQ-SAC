#!/usr/bin/env python
# coding: utf-8

# In[15]:
#target_line
#1,Tig0010.csv, target :  Tig, simulation :  isochor, measurnment_type :  OHEX ;max,Ignition_mode :  reflected, Fuel :  x->H2 = 0.0587, Oxidizer :  x->O2 = 0.0295,BG1 :  x->Ar = 0.9118, BG2 : ,BG3 : , T : 1160.0, P : 412000.0, Phi : N/N, observed : 12.195121951219512, deviation : 1.2195121951219514,data_weight = 13
#progress
#




import xml.etree.ElementTree as ET
import numpy as np
import os
fuel_type = "Multi"
fuel_list = ["H2","CO"]
input_fuel = {'a':'H2','b':'CO'}
input_oxidizer = 'O2'
files = [i for i in os.listdir() if i.endswith("xml")]
count = 0
#print(len(files))
add = ""

#raise AssertionError("Stop")
for file_ in files:
	Target = None
	Target_kind = None
	Target_mode = None
	bibliography = None
	commonProp = None
	dataGroup = None
	add_mod = []
	print(file_)
	dataset = file_.strip(".xml")
	add+=f"{dataset}:\n"
	
	tree = ET.parse(file_)
	root = tree.getroot()
	for i in root.findall(".//*"):
		if "experimentType" in i.tag:
			Target = root.find(".//experimentType").text
	
		if "kind" in i.tag:
			Target_kind = root.find(".//kind").text
		
		if "mode" in i.tag:
			Target_mode = root.find(".//mode").text
			Target_add = root.findall(".//mode")
			for i in Target_add:		
				add_mod.append(i.text)
			
		if "bibliographyLink" in i.tag:
			bibliography = root.find(".//bibliographyLink")

		
		if "commonProperties" in i.tag:

			commonProp  = root.find(".//commonProperties")
			
		if "dataGroup" in i.tag:
			dataGroup = root.findall(".//dataGroup")
		
		
	
	count = 0
	for group in dataGroup:
		for child in group:
			if child.tag == "property" and child.attrib["name"] == "composition":
				count +=1
	if count == 0:
		add += " solver: cantera\n type: phi\n\n"
	else: 
		add += " solver: cantera\n\n"
	#raise AssertionError("Stop")
	#Target_kind = root.find(".//kind").text
	#Target_mode = root.find(".//mode").text
	#bibliography = root.find(".//bibliographyLink")
	#commonProp  = root.find(".//commonProperties")
	#dataGroup = root.findall(".//dataGroup")
	 
	
		
	for i in root.findall("*"):
		if "ignitionType" in i.tag:
		    MeasurnmentType = root.find(".//ignitionType").attrib["target"]+root.find(".//ignitionType").attrib["type"]
		else:
		    MeasurnmentType = ""
		

	# In[18]:


	properties = {}
	concentartion_Dict = {}
	unsrt = {}
	datagroupDict = {}


	# In[31]:
	
	for child in commonProp:
		
		#print(child.[@attribName]
		if child.tag == "property":
		    if child.attrib["name"] == "initial composition":
		        for zero in child:       
		            temp = {}
		            for minusOne in zero:
		                if minusOne.tag == "speciesLink":
		                    temp["name"] = minusOne.attrib["preferredKey"]
		                    temp["CAS"] = minusOne.attrib["CAS"]
		                    temp["InChI"] = minusOne.attrib["InChI"]
		                    temp["SMILES"] = minusOne.attrib["SMILES"]
		                    temp["chemName"] = minusOne.attrib["chemName"]
		                elif minusOne.tag == "amount":
		                    value = minusOne.text
		                    temp["value"] = value
		                    temp["unit"] = minusOne.attrib["units"]
		            concentartion_Dict[temp["name"]] = temp
		        properties["initial composition"] = concentartion_Dict
		    if child.attrib["name"] == "evaluated standard deviation" or child.attrib["name"] == "uncertainty":
		        
		        type_ = child.attrib["reference"]
		        unsrt_Dict = {}
		        
		        if type_ == "composition":
		        	for i in child:
		        		if i.tag == "speciesLink":
		        			unsrt_Dict["specification"] = i.attrib
		        			type_ = child.attrib["reference"]+"_"+i.attrib["preferredKey"]
		       
		        unsrt_Dict["type"] = child.attrib["kind"]
		        if "method" in child.attrib:
		            unsrt_Dict["method"] = child.attrib["method"]
		        else:
		            unsrt_Dict["method"] = "instrumental error"
		            unsrt_Dict["bound"] = child.attrib["bound"]
		            unsrt_Dict["label"] = child.attrib["label"]
		        unsrt_Dict["sourcetype"] = child.attrib["sourcetype"]
		        unsrt_Dict["units"] = child.attrib["units"]
		        for zero in child:
		            #print(zero.text)
		            unsrt_Dict["value"] = zero.text

		        unsrt[type_] = unsrt_Dict
		        print(unsrt_Dict)
		        properties["Unsrt"] = unsrt
		    if child.attrib["name"] == "pressure":
		        temp = {}
		        temp["label"] = child.attrib["label"]
		        temp["unit"] = child.attrib["units"]
		        temp["sourcetype"] = child.attrib["sourcetype"]
		        for zero in child:
		            if zero.tag == "value":
		                temp["value"] = zero.text
		        properties["pressure"] = temp
		    if child.attrib["name"] == "temperature":
		        temp = {}
		        temp["label"] = child.attrib["label"]
		        temp["unit"] = child.attrib["units"]
		        temp["sourcetype"] = child.attrib["sourcetype"]
		        for zero in child:
		            if zero.tag == "value":
		                temp["value"] = zero.text
		        properties["temperature"] = temp
	
		

	for group in dataGroup:
		Dict = {}
		points = {}
		tempID = {}
		temp = {}
		tempID["id"] = group.attrib["id"]
		#print(group.attrib["id"])
		for ind,child in enumerate(group):
		    if child.tag == "property":
		        if child.attrib["name"] == "composition":
		            tem = {}
		            for zero in child:
		                tem["specification"] = zero.attrib
		            tem["name"] = child.attrib["name"]
		            tem["id"] = child.attrib["id"]
		            tem["label"] = child.attrib["label"]
		            tem["sourcetype"] = child.attrib["sourcetype"]
		            tem["units"] = child.attrib["units"]
		            temp[child.attrib["id"]] = tem		        	
		        elif child.attrib["name"] == "uncertainty" and child.attrib["units"] == "cm/s" or child.attrib["name"] == "evaluated standard deviation" and child.attrib["units"] == "cm/s":
		        	if "laminar burning velocity" in unsrt:
		        		unsrt["laminar burning velocity"]["var"] = child.attrib["id"]
		        	else:
		        		unsrt["laminar burning velocity"] = {}
		        		unsrt["laminar burning velocity"]["var"] = child.attrib["id"]
		        		unsrt["laminar burning velocity"]["type"] = child.attrib["kind"]
		        		unsrt["laminar burning velocity"]["units"] = child.attrib["units"]
		        		
		        	
		        else:
		            temp[child.attrib["id"]] = child.attrib
		    if child.tag == "dataPoint":
		        data = {}
		        for zero in child:
		            #print(zero.tag)
		            data[zero.tag]=zero.text
		        points[str(ind-2)] = data
		Dict["Variables"] = temp
		Dict["DataPoints"] = points
		datagroupDict[tempID["id"]] = Dict
	
	#no of data points
	ndpt = len(datagroupDict["dg1"]["DataPoints"])
	
	if "laminar burning velocity" in unsrt:
		if "var" in unsrt["laminar burning velocity"]:
			unsrt["laminar burning velocity"]["value"] = []
			for points in datagroupDict["dg1"]["DataPoints"]:
				unsrt["laminar burning velocity"]["value"].append(datagroupDict["dg1"]["DataPoints"][points][unsrt["laminar burning velocity"]["var"]])
	#print(unsrt)
	#print(ndpt)
	#print(datagroupDict["dg1"])
	
	
	# In[30]:
	
	if "volume" in properties:
		volume = float(properties["colume"]["value"])
		volume_unit = properties["flow rate"]["unit"]
	else:
		volume = ""
		volume_unit = ""
	if "residence time" in properties:
		res_unit = properties["flow rate"]["unit"]
		if res_unit == "ms":
			res_time = 0.001*float(properties["colume"]["value"])
			res_unit = "s"
		elif res_unit == "us":
			res_time = 0.000001*float(properties["colume"]["value"])
			res_unit = "s" 
	else:	
		res_time = ""
		res_unit = ""
	
	# Extracting pressure
	
	if "pressure" in properties:
		p_unit = properties["pressure"]["unit"]
		if p_unit == "mbar":
			fact = 0.001
			p_unit = "bar"
		elif p_unit == "bar" or p_unit == "atm":
			fact = 1.0
			p_unit = "bar"
		elif p_unit == "Pa":
			fact = 1.0
			p_unit == "Pa"
		elif p_unit == "MPa":
			fact = 1e6
			p_unit == "Pa"
		else:
			raise AssertionError("New pressure unit found in {}".format(file_))
		pressure = fact*float(properties["pressure"]["value"])
		p_sigma = 0 #by default
		pressure = np.repeat(pressure,ndpt)
		if "pressure" in unsrt:
		    p_sigma = unsrt["pressure"]["value"]
	
	else:
	
		pressure = []
		for var in datagroupDict["dg1"]["Variables"]:
			#print(var)
			#print(datagroupDict["dg1"]["DataPoints"])
			if datagroupDict["dg1"]["Variables"][var]["name"] == "pressure":
				p_unit = datagroupDict["dg1"]["Variables"][var]["units"] 
				
				if p_unit == "mbar":
					fact = 0.001
					p_unit = "bar"
				elif p_unit == "bar" or p_unit == "atm":
					fact = 1.0
					
				elif p_unit == "Pa":
					fact = 1.0
					p_unit == "Pa"
				elif p_unit == "MPa":
					fact = 1e6
					p_unit = "Pa"
				else:
					raise AssertionError("New pressure unit found in {}".format(file_))
				for point in datagroupDict["dg1"]["DataPoints"]:
			    		pressure.append(fact*float(datagroupDict["dg1"]["DataPoints"][point][var]))

		if "pressure" in unsrt:
		    p_sigma = unsrt["pressure"]["value"]


		
	if "flow rate" in properties:
		
		flowRate = float(properties["flow rate"]["value"])
		flow_unit = properties["flow rate"]["unit"]
		flow_sigma = 0 #by default
		if "flow rate" in unsrt:
		    flow_sigma = unsrt["flow rate"]["value"]
	
	elif "flow rate" in datagroupDict:
	
		flowRate = []
		for var in datagroupDict["dg1"]["Variables"]:
		    if datagroupDict["dg1"]["Variables"][var]["name"] == "flow rate":
		        flow_unit = datagroupDict["dg1"]["Variables"][var]["units"] 
		        for point in datagroupDict["dg1"]["DataPoints"]:
		            flowRate.append(datagroupDict["dg1"]["DataPoints"][point][var])
		
		if "flow rate" in unsrt:
		    flow_sigma = unsrt["flow rate"]["value"]
	else:
		flow_unit = ""
		flowRate = np.repeat(flow_unit,ndpt)
	
	# In[22]:

	
#print(datagroupDict)
	
	#Extracting initial composition
	if "initial composition" in properties:
		#print(properties["initial composition"])
		List_fuel=[]
		List_fuel_x= []
		input_bath_gas = {}
		bathgas = []
		bathgas_conc = []
		ox = []
		ox_x = []
		bathgasDict = {}
		for i in properties["initial composition"]:
			if i in fuel_list:
				if fuel_type == "Multi":
					input_fuel_x = {}
					for fuel in input_fuel:
						temp = input_fuel[fuel]
						if temp == i: 
							input_fuel_x[str(fuel)] = float(properties["initial composition"][i]["value"])
				else:

					input_fuel_x = float(properties["initial composition"][i]["value"])
			elif i == "O2":
				input_oxidizer_x = float(properties["initial composition"][i]["value"])
			else:
				input_bath_gas[i] = float(properties["initial composition"][i]["value"])
		#print(properties)
		bath_gas = {}
		bath_gas_x = {}
		list_a = []
		for j in range(97,97+len(input_bath_gas)):
			list_a.append(chr(j))
		
		for ind,i in enumerate(input_bath_gas):
			bath_gas[list_a[ind]] = i
			bath_gas_x[list_a[ind]] = input_bath_gas[i]
		
		List_fuel = np.repeat(input_fuel,ndpt)
		List_fuel_x = np.repeat(input_fuel_x,ndpt)
		bathgas = np.repeat(bath_gas,ndpt)
		bathgas_conc = np.repeat(bath_gas_x,ndpt)
		ox = np.repeat(input_oxidizer,ndpt)
		ox_x = np.repeat(input_oxidizer_x,ndpt)
		
	else:
		input_fuel_x = {}
		List_fuel=[]
		List_fuel_x= []
		bathgas = []
		bathgas_conc = []
		ox = []
		ox_x = []
		bathgasDict = {}
		for var in datagroupDict["dg1"]["Variables"]:
			if datagroupDict["dg1"]["Variables"][var]["name"] == "composition":
				if datagroupDict["dg1"]["Variables"][var]["specification"]["preferredKey"] in fuel_list:
					for fuel in fuel_list:
						fuel_conc = []	
						if datagroupDict["dg1"]["Variables"][var]["specification"]["preferredKey"] == fuel:
							#fuel = np.repeat(datagroupDict["dg1"]["Variables"][var]["specification"]["preferredKey"],ndpt)
							#fuel_unit = datagroupDict["dg1"]["Variables"][var]["units"]
							for point in datagroupDict["dg1"]["DataPoints"]:
								fuel_conc.append(float(datagroupDict["dg1"]["DataPoints"][point][var]))
							input_fuel_x[fuel] = fuel_conc
							continue
				elif datagroupDict["dg1"]["Variables"][var]["specification"]["preferredKey"] == "O2":
					ox = np.repeat(datagroupDict["dg1"]["Variables"][var]["specification"]["preferredKey"],ndpt)
					#oxidizer_unit = datagroupDict["dg1"]["Variables"][var]["units"]

					for point in datagroupDict["dg1"]["DataPoints"]:
						ox_x.append(float(datagroupDict["dg1"]["DataPoints"][point][var]))
					continue
				else:
					component = datagroupDict["dg1"]["Variables"][var]["specification"]["preferredKey"]
					#bathgas_unit = datagroupDict["dg1"]["Variables"][var]["units"]
					temp_list = []
					#temp_list[component] = []
					for point in datagroupDict["dg1"]["DataPoints"]:
						temp_list.append(float(datagroupDict["dg1"]["DataPoints"][point][var]))
					bathgasDict[component] = temp_list
		
		List_fuel = None
		List_fuel_x = []
		bath_gas = {}
		bathgas = None
		bathgas_conc = []
		list_a = []
		for j in range(97,97+len(bathgasDict)):
			list_a.append(chr(j))
		list_b = []
		for j in range(97,97+len(input_fuel_x)):
			list_b.append(chr(j))
		for ind,i in enumerate(bathgasDict):
			bath_gas[list_a[ind]] = i 
		bathgas = np.repeat(bath_gas,ndpt)
		for i in range(ndpt):
			bath_gas_x = {}
			for ind,j in enumerate(bathgasDict):
				bath_gas_x[list_a[ind]] = bathgasDict[j][i]
			bathgas_conc.append(bath_gas_x)
		temp_fuel = {}
		for ind,i in enumerate(input_fuel_x):
			temp_fuel[list_b[ind]] = i
		List_fuel = np.repeat(temp_fuel,ndpt)
		for i in range(ndpt):
			temp_fuel_x = {}
			for ind,j in enumerate(input_fuel_x):
				temp_fuel_x[list_b[ind]] = input_fuel_x[j][i]
			List_fuel_x.append(temp_fuel_x)	
			
		#print(ox_x)
	#Extracting temperature

	for ind,i in enumerate(List_fuel):
		for j in i:
			if j not in List_fuel_x[ind]:
				List_fuel_x[ind][j] = 0.0
	if "temperature" in properties:
		t_unit = properties["temperature"]["unit"]
		temperature = np.repeat(float(properties["temperature"]["value"]),ndpt)
		t_sigma = 0 #by default
		if "temperature" in unsrt:
		    t_sigma = unsrt["temperature"]["value"]
	else:
		temperature = []
		#print(datagroupDict['dg1']["Variables"])
		for var in datagroupDict["dg1"]["Variables"]:
		    if datagroupDict["dg1"]["Variables"][var]["name"] == "temperature":
		        t_unit = datagroupDict["dg1"]["Variables"][var]["units"] 
		        for point in datagroupDict["dg1"]["DataPoints"]:
		            temperature.append(float(datagroupDict["dg1"]["DataPoints"][point][var]))
		if "temperature" in unsrt:
		    t_sigma = unsrt["temperature"]["value"]

	# In[34]:

	#Extracting experimental targets
	exp_target = []
	if Target_kind == "flame" and "HFM" not in add_mod and Target_mode == "premixed" or Target_mode == "spherical" or Target_mode == "modified Bunsen burner" or  Target_mode == "unstretched" or Target_mode == "slot burner" or Target_mode == "laminar" or Target_mode == "premixed" or Target_mode == "outwardly propagating spherical flame" or Target_mode == "constant volume combustion chamber" or Target_mode == "extrapolation method to zero stretch: LS" or Target_mode == "extrapolation method to zero stretch: NQ":
		exp = 'laminar burning velocity'
		exp_tag = "Fls"
		Simulation_type = "UnstretchedPremixed"
		reactorType = ""
		ignitionMode = ""
		flameType = "Unstretched"
	elif Target_kind == "shock tube":
		exp = "ignition delay"
		exp_tag = "Tig"
		Simulation_type = "Isochor Homo Reactor"
		reactorType = ""
		ignitionMode = "reflected"
		flameType = ""
	elif Target_kind == "stirred reactor" or Target_kind == "flow reactor":
		exp = "composition"
		exp_tag = "Flw"
		Simulation_type = "Isobar Homo Reactor"
		MeasurnmentType = "CONC"
		reactorType = "FlowReactor"
		ignitionMode = ""
		flameType = ""
		#if ""
	elif Target_kind == "flame" and Target_mode == "burner-stabilized":
		exp = "composition"
		exp_tag = "Flf"
		Simulation_type = "UnstretchedPremixed"
		MeasurnmentType = "CONC"
		reactorType = ""
		ignitionMode = ""
		flameType = "HeatFlux"
	elif Target_kind == "flame" and "HFM" in add_mod:
		exp = "laminar burning velocity"
		exp_tag = "Flf"
		Simulation_type = "UnstretchedPremixed"
		MeasurnmentType = "velocity"
		reactorType = ""
		ignitionMode = ""
		flameType = "HeatFlux"
	
	elif Target_kind == "rapid compression machine":
		exp = "ignition delay"
		exp_tag = "RCM"
		volume_profile = datagroupDict["dg2"]
		print(volume_profile)
	else:
		raise AssertionError("New exp found !!! {}".format(file_))
	#print(exp)
	specific_cond = {}
	exp_composition = {}
	for var in datagroupDict['dg1']['Variables']:
		if exp_tag == "Flw":
			if "residence time" in datagroupDict['dg1']['Variables'][var]['name']:
				specific_cond["residence time"] = []
				x_unit = datagroupDict['dg1']['Variables'][var]['units']
				x_id = datagroupDict['dg1']['Variables'][var]['id']
				if x_unit == "s":
					fact = 1.0
				elif x_unit == "ms":
					fact = 0.001
				elif x_unit == "us":
					fact = 0.000001
				else:
					raise AssertionError("New units found for time in FLW {}".format(file_))
					
				for point in datagroupDict["dg1"]["DataPoints"]:
					specific_cond["residence time"].append(fact*float(datagroupDict["dg1"]["DataPoints"][point][var]))
					
			if "composition" in datagroupDict['dg1']['Variables'][var]['name']:
				x_id = datagroupDict['dg1']['Variables'][var]['id']
				exp_unit = datagroupDict['dg1']['Variables'][var]['units']
				species = datagroupDict['dg1']['Variables'][var]["specification"]["preferredKey"]
				exp_composition[species] = {}
				exp_composition[species]["observation"] = []
				for point in datagroupDict["dg1"]["DataPoints"]:
					exp_composition[species]["observation"].append(float(datagroupDict["dg1"]["DataPoints"][point][var]))
				exp_ = exp+"_"+species
				if exp_ in unsrt:
					exp_sigma = float(unsrt[exp_]["value"])
					exp_sigma_type = unsrt[exp_]["type"]
					exp_std_dvnt_unit = unsrt[exp_]["units"]
				
				if unsrt[exp_]["units"] == "ppm":
					exp_sigma = float(unsrt[exp_]["value"])/10**6
				elif unsrt[exp_]["units"] == "mole fraction":
					exp_sigma = float(unsrt[exp_]["value"])
				else:
					raise AssertionError("New type of unsrt {}".format(file_))

				if exp_sigma_type == "relative":
					exp_std_dnvt = exp_sigma*np.asarray(exp_composition[species]["observation"])
				elif exp_sigma_type == "absolute":
					exp_std_dnvt = np.repeat(exp_sigma,ndpt)
				exp_composition[species]["std_unit"] = exp_std_dvnt_unit
				exp_composition[species]["std"] = exp_std_dnvt
			
			if "time" in datagroupDict['dg1']['Variables'][var]['name']:
				print("yes")
				specific_cond["time"] = []
				x_unit = datagroupDict['dg1']['Variables'][var]['units']
				x_id = datagroupDict['dg1']['Variables'][var]['id']
				
				if x_unit == "s":
					fact = 1.0
				elif x_unit == "ms":
					fact = 0.001
				elif x_unit == "us":
					fact = 0.000001
				else:
					raise AssertionError("New units found for time in FLW {}".format(file_))
				for point in datagroupDict["dg1"]["DataPoints"]:
					specific_cond["time"].append(fact*float(datagroupDict["dg1"]["DataPoints"][point][var]))
			
			if "temperature" in datagroupDict['dg1']['Variables'][var]['name']:
				x_id = datagroupDict['dg1']['Variables'][var]['id']
				specific_cond["temperature"] = []
				for point in datagroupDict["dg1"]["DataPoints"]:
					specific_cond["temperature"].append(float(datagroupDict["dg1"]["DataPoints"][point][var]))
		
		
		
			
		else:
			if exp in datagroupDict['dg1']['Variables'][var]['name']:
				print("yes")
				exp_unit = datagroupDict['dg1']['Variables'][var]['units']
				for point in datagroupDict["dg1"]["DataPoints"]:
					exp_target.append(float(datagroupDict["dg1"]["DataPoints"][point][var]))
				
			
			if exp in unsrt:
				if "var" in unsrt[exp]:
					exp_sigma = np.asarray(unsrt[exp]["value"])
					exp_sigma_type = unsrt[exp]["type"]
					exp_std_dvnt_unit = unsrt[exp]["units"]
				else:
				
					exp_sigma = unsrt[exp]["value"]
					exp_sigma_type = unsrt[exp]["type"]
					exp_std_dvnt_unit = unsrt[exp]["units"]
					if unsrt[exp]["units"] == "ppm":
						exp_sigma = unsrt[exp]["value"]/10**6
			else:
				print(unsrt)
				raise AssertionError("New type of unsrt {}".format(file_))

			if exp_sigma_type == "relative":
				#print(exp_sigma)
				#print(exp_target)
				exp_std_dnvt = exp_sigma
			elif exp_sigma_type == "absolute":
				if len(exp_sigma) == ndpt:
					exp_std_dnvt = exp_sigma
				else:
					exp_std_dnvt = np.repeat(exp_sigma,ndpt)

	if exp_tag == "Flw":
		if "time" not in specific_cond:
			specific_cond["time"] = np.repeat("",ndpt)
		if "temperature" not in specific_cond:
			specific_cond["temperature"] = np.repeat("",ndpt)
		if "residence time" not in specific_cond:
			specific_cond["residence time"] = np.repeat("",ndpt)
	
	# In[35]:
	#print(exp_std_dnvt)
	#print(exp_target)
	if exp_sigma_type == "relative":
		#print(exp_sigma)
		#print(exp_target)
		exp_std_dnvt = np.asarray(exp_target)*float(exp_sigma)
	
	
	Input_units = {}
	Input_units["conc"] = 'mol'
	Input_units["P"] = p_unit
	Input_units["T"] = t_unit
	Input_units["observed"] = exp_unit
	Input_units["flow"] = flow_unit
	#Input_data_structure = ""
	#for i in datagroupDict['dg1']["Variables"]:
	#	Input_data_structure +=datagroupDict['dg1']["Variables"][i]["label"].strip("[]")+"\t"


	opt_target_file = open("target.opt","a+")
	string = ""

	#dataUnits = [fuel_unit,t_unit,p_unit,exp_unit,flow_unit]
	if exp_tag == "Flw":
		for each in exp_composition:
			
			for i in range(ndpt):
				
				string+="{}\t|{}_{}\t| target -- {}\t| simulation -- {}\t| measurnment_type -- {}\t|Ignition_mode -- {}\t| Flame_type -- {}\t| Reactor_type -- {}\t| Fuel_type -- {}\t| Fuel -- x->{} = {}\t| Oxidizer -- x->{} = {}\t| Bath_gas -- x->{}={}\t|T -- {}\t| P -- {}\t| flow_rate -- {}\t|Phi -- {}\t| volume -- {}\t|volume_unit -- {}\t| residence_time -- {}\t|observed_quantity -- {}\t| specific_details -- {}\t| observation_type -- {}\t|observed -- {}\t| obs_unit --{} \t| deviation -- {}\t|data_weight -- {}\t|units -- {}\n".format(count+int(i+1),file_.split(".")[0],each,exp_tag,Simulation_type,MeasurnmentType,ignitionMode,flameType,reactorType,fuel_type,List_fuel[i],List_fuel_x[i],ox[i],ox_x[i],bathgas[i],bathgas_conc[i],temperature[i],pressure[i],flowRate[i],"N/A",volume,volume_unit,res_time,each,{"time":specific_cond["time"][i],"residence time":specific_cond["residence time"][i],"temperature":specific_cond["temperature"][i]},"profile",exp_composition[each]["observation"][i],exp_composition[each]["std_unit"],exp_composition[each]["std"][i],ndpt,Input_units)
	
	else:
		for i in range(ndpt):
			string+="{}\t|{}\t| target -- {}\t| simulation -- {}\t| measurnment_type -- {}\t|Ignition_mode -- {}\t| Flame_type -- {}\t| Reactor_type -- {}\t| Fuel_type -- {}\t| Fuel -- x->{} = {}\t| Oxidizer -- x->{} = {}\t| Bath_gas -- x->{}={}\t|T -- {}\t| P -- {}\t| flow_rate -- {}\t|Phi -- {}\t| observed -- {}\t| obs_unit --{} \t| deviation -- {}\t|data_weight -- {}\t|units -- {}\n".format(count+int(i+1),file_.split(".")[0],exp_tag,Simulation_type,MeasurnmentType,ignitionMode,flameType,reactorType,fuel_type,List_fuel[i],List_fuel_x[i],ox[i],ox_x[i],bathgas[i],bathgas_conc[i],temperature[i],pressure[i],flowRate[i],"N/A",exp_target[i],exp_std_dvnt_unit,exp_std_dnvt[i],ndpt,Input_units)

	count = count+ndpt
	print(string)
	opt_target_file.write(string)
	opt_target_file.close()
	# In[ ]:

add_file = open("addundem.add","+w").write(add)

