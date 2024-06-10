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
input_fuel = "H2"
input_oxidizer = "O2"
files = [i for i in os.listdir() if i.endswith("xml")]

count = 0
for file_ in files:
	print(file_)
	tree = ET.parse(file_)
	root = tree.getroot()
	Target = root.find(".//experimentType").text
	Target_kind = root.find(".//kind").text
	Target_mode = root.find(".//mode").text
	bibliography = root.find(".//bibliographyLink")
	commonProp  = root.find(".//commonProperties")
	dataGroup = root.findall(".//dataGroup")
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
	Dict = {}
	points = {}


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
	#print(unsrt)


	# In[30]:


	for group in dataGroup:
		tempID = {}
		temp = {}
		tempID["id"] = group.attrib["id"]
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
	print(ndpt)
	
	
	
	# Extracting pressure
	if "pressure" in properties:
		
		pressure = float(properties["pressure"]["value"])
		p_unit = properties["pressure"]["unit"]
		p_sigma = 0 #by default
		if "pressure" in unsrt:
		    p_sigma = unsrt["pressure"]["value"]
	
	else:
	
		pressure = []
		for var in datagroupDict["dg1"]["Variables"]:
		    if datagroupDict["dg1"]["Variables"][var]["name"] == "pressure":
		        p_unit = datagroupDict["dg1"]["Variables"][var]["units"] 
		        for point in datagroupDict["dg1"]["DataPoints"]:
		            pressure.append(datagroupDict["dg1"]["DataPoints"][point][var])
		
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


	#Extracting initial composition
	if "initial composition" in properties:
		
		bathgas = []
		bathgas_conc = []
		bathgas_unit = []
		#print(properties)
		for component in properties["initial composition"]:
		    if component == input_fuel:
		        fuel = np.repeat(properties["initial composition"][component]["name"],ndpt)
		        fuel_conc = np.repeat(float(properties["initial composition"][component]["value"]),ndpt)
		        fuel_unit = properties["initial composition"][component]["unit"]
		        continue
		    elif component == input_oxidizer:
		        oxidizer = np.repeat(properties["initial composition"][component]["name"],ndpt)
		        oxidizer_conc = np.repeat(float(properties["initial composition"][component]["value"]),ndpt)
		        oxidizer_unit = properties["initial composition"][component]["unit"]
		        continue
		    else:
		        bathgas.append(properties["initial composition"][component]["name"])
		        bathgas_conc.append(float(properties["initial composition"][component]["value"]))
		        bathgas_unit.append(properties["initial composition"][component]["unit"])
		
		if len(bathgas) == 1:
			empty = ""
			bg1_gas = np.repeat(bathgas[0],ndpt)
			bg1_gas_conc = np.repeat(bathgas_conc[0],ndpt)
			bg2_gas = np.repeat(empty,ndpt)
			bg2_gas_conc = np.repeat(empty,ndpt)
			bg3_gas = np.repeat(empty,ndpt)
			bg3_gas_conc = np.repeat(empty,ndpt)
		
		elif len(bathgas) == 2:
			empty = ""
			bg1_gas = np.repeat(bathgas[0],ndpt)
			bg1_gas_conc = np.repeat(bathgas_conc[0],ndpt)
			bg2_gas = np.repeat(bathgas[1],ndpt)
			bg2_gas_conc = np.repeat(bathgas_conc[1],ndpt)
			bg3_gas = np.repeat(empty,ndpt)
			bg3_gas_conc = np.repeat(empty,ndpt)
		
		elif len(bathgas) == 3:
			bg1_gas = np.repeat(bathgas[0],ndpt)
			bg1_gas_conc = np.repeat(bathgas_conc[0],ndpt)
			bg2_gas = np.repeat(bathgas[1],ndpt)
			bg2_gas_conc = np.repeat(bathgas_conc[1],ndpt)
			bg3_gas = np.repeat(bathgas[2],ndpt)
			bg3_gas_conc = np.repeat(bathgas_conc[2],ndpt)
			
		else:
			print("bathgas exceeded the max limit")
	else:
		
		fuel_conc = []
		oxidizer_conc = []
		bathgas = []
		bathgas_conc = []
		
		for var in datagroupDict["dg1"]["Variables"]:
		    
		    if datagroupDict["dg1"]["Variables"][var]["name"] == "composition":
		        if datagroupDict["dg1"]["Variables"][var]["specification"]["preferredKey"] == "H2":
		            fuel = np.repeat(datagroupDict["dg1"]["Variables"][var]["specification"]["preferredKey"],ndpt)
		            fuel_unit = datagroupDict["dg1"]["Variables"][var]["units"]
		            for point in datagroupDict["dg1"]["DataPoints"]:
		                fuel_conc.append(float(datagroupDict["dg1"]["DataPoints"][point][var]))
		        elif datagroupDict["dg1"]["Variables"][var]["specification"]["preferredKey"] == "O2":
		            oxidizer = np.repeat(datagroupDict["dg1"]["Variables"][var]["specification"]["preferredKey"],ndpt)
		            oxidizer_unit = datagroupDict["dg1"]["Variables"][var]["units"]

		            for point in datagroupDict["dg1"]["DataPoints"]:
		                oxidizer_conc.append(float(datagroupDict["dg1"]["DataPoints"][point][var]))
		        else:
		            bathgas.append(datagroupDict["dg1"]["Variables"][var]["specification"]["preferredKey"])
		            bathgas_unit = datagroupDict["dg1"]["Variables"][var]["units"]

		            temp_list = []
		            for point in datagroupDict["dg1"]["DataPoints"]:
		                temp_list.append(float(datagroupDict["dg1"]["DataPoints"][point][var]))
		            bathgas_conc.append(temp_list)

		
		if len(bathgas) == 1:
			bg1_gas = np.repeat(bathgas[0],ndpt)
			bg1_gas_conc = bathgas_conc[0]
			bg2_gas = np.repeat("",ndpt)
			bg2_gas_conc = np.repeat("",ndpt)
			bg3_gas = np.repeat("",ndpt)
			bg3_gas_conc = np.repeat("",ndpt)
		
		elif len(bathgas) == 2:
			bg1_gas = np.repeat(bathgas[0],ndpt)
			bg1_gas_conc = bathgas_conc[0]
			bg2_gas = np.repeat(bathgas[1],ndpt)
			bg2_gas_conc = bathgas_conc[1]
			bg3_gas = np.repeat("",ndpt)
			bg3_gas_conc = np.repeat("",ndpt)
		
		elif len(bathgas) == 3:
			bg1_gas = np.repeat(bathgas[0],ndpt)
			bg1_gas_conc = bathgas_conc[0]
			bg2_gas = np.repeat(bathgas[1],ndpt)
			bg2_gas_conc = bathgas_conc[1]
			bg3_gas = np.repeat(bathgas[2],ndpt)
			bg3_gas_conc = bathgas_conc[2]
			
		else:
			print("bathgas exceeded the max limit")
		print(bg1_gas_conc)
	#Extracting temperature
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
	if Target_kind == "flame" and Target_mode == "premixed":
		exp = 'laminar burning velocity'
		exp_tag = "fls"
		Simulation_type = "UnstretchedPremixed"
		
	elif Target_kind == "shock tube":
		exp = "ignition delay"
		exp_tag = "Tig"
		Simulation_type = "isochor"
		
	for var in datagroupDict['dg1']['Variables']:
		if exp in datagroupDict['dg1']['Variables'][var]['name']:
		    exp_unit = datagroupDict['dg1']['Variables'][var]['units']
		    for point in datagroupDict["dg1"]["DataPoints"]:
		        exp_target.append(float(datagroupDict["dg1"]["DataPoints"][point][var]))
	
	if exp in unsrt:
		exp_sigma = float(unsrt[exp]["value"])
		exp_sigma_type = unsrt[exp]["type"]
	
	if exp_sigma_type == "relative":
		exp_std_dnvt = exp_sigma*np.asarray(exp_target)
	elif exp_sigma_type == "absolute":
		exp_std_dnvt = np.repeat(exp_sigma,ndpt)


	
	# In[35]:


	Input_units = "{}\t{}\t{}\t{}".format(fuel_unit,p_unit,t_unit,exp_unit)
	Input_data_structure = ""
	for i in datagroupDict['dg1']["Variables"]:
		Input_data_structure +=datagroupDict['dg1']["Variables"][i]["label"].strip("[]")+"\t"

	pressure = np.repeat(pressure,ndpt)

	opt_target_file = open("target.opt","a+")
	string = ""
	
	
	dataUnits = [fuel_unit,t_unit,p_unit,exp_unit,flow_unit]
	
	for i in range(ndpt):
		print(pressure[i])
		string+="{}\t,{}\t, target : {}\t, simulation : {}\t, measurnment_type : {}\t,Ignition_mode : {}\t, Fuel : x->{} = {}\t, Oxidizer : x->{} = {}\t, BG1 : x->{} = {}\t, BG2 : x->{} = {}\t,BG3 : x->{} = {}\t,T : {}\t, P : {}\t, flow_rate : {}\t,Phi : {}\t, observed : {}\t, deviation : {}\t,data_weight : {}\t,units : {}\n".format(int(count+i),file_,Target,Simulation_type,MeasurnmentType,Target_mode,fuel[i],fuel_conc[i],oxidizer[i],oxidizer_conc[i],bg1_gas[i],bg1_gas_conc[i],bg2_gas[i],bg2_gas_conc[i],bg3_gas[i],bg3_gas_conc[i],temperature[i],pressure[i],flowRate[i],"N/A",exp_target[i],exp_std_dnvt[i],ndpt,dataUnits)

	count = count+ndpt
	print(string)
	opt_target_file.write(string)
	opt_target_file.close()
	# In[ ]:




