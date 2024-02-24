import os
import math
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import PolynomialFeatures
from sklearn import linear_model
from numpy import linalg as LA
import json
from scipy.interpolate import InterpolatedUnivariateSpline
import statistics
import pandas as pd
import time
from scipy.optimize import minimize
class combustion_target():
###Class definition and acquisition of input parameters.	
	def __init__(self,data,addendum,index):
		self.molecularWt = {}
		self.molecularWt["CO"] = 28
		self.molecularWt["H2"] = 2
		self.molecularWt["O2"] = 32
		self.molecularWt["co"] = 28
		self.molecularWt["h2"] = 2
		self.molecularWt["o2"] = 32
		self.molecularWt["CO2"] = 0
		self.molecularWt["AR"] = 0
		self.molecularWt["HE"] = 0
		self.molecularWt["CO2"] = 0
		self.molecularWt["Ar"] = 0
		self.molecularWt["He"] = 0
		self.molecularWt["N2"] = 0
		self.molecularWt["H2O"] = 0
		self.molecularWt["NC7H16"] = 100.21
		self.molecularWt["MB-C5H10O2"] = 86
		self.stoichiometry = {}
		self.stoichiometry["H2"] = 0.5
		self.stoichiometry["CO"] = 0.5
		self.stoichiometry["O2"] = 0.0
		self.stoichiometry["N2"] = 0.0
		self.stoichiometry["h2"] = 0.5
		self.stoichiometry["co"] = 0.5
		self.stoichiometry["o2"] = 0.0
		self.stoichiometry["n2"] = 0.0
		self.stoichiometry["HE"] = 0.0
		self.stoichiometry["AR"] = 0.0
		self.stoichiometry["He"] = 0.0
		self.stoichiometry["Ar"] = 0.0
		self.stoichiometry["H2O"] = 0.0
		self.stoichiometry["CO2"] = 0.0
		self.stoichiometry["NC7H16"] = 11
		self.stoichiometry["MB-C5H10O2"] = 6.5
		self.data = data
		parameters = self.data.split('|')
		self.dataSet_id = parameters[1].strip("\t")
		self.calculated = 0
		self.x = None
		self.index = index
		self.case_index = "case-"+str(index)
		self.units = self.f_unit = self.t_unit = self.temperature_factor = self.p_unit = self.pressure_factor = self.target_unit = self.target_factor = self.flow_unit  = self.flow_rate = self.target_key = self.input_file = self.target = self.simulation = self.temperature = self.fuel = self.oxygen = self.nitrogen = self.argon = self.pressure = self.phi = self.observed = self.d_weight = self.d_set = self.Tend = self.species = self.s_p_name = None
		self.add = {}	
		
		#raise ValueError("addendum")
		if self.dataSet_id in addendum:
			self.add = addendum[self.dataSet_id]
			if self.add == None:
				self.add = {}
				self.add["solver"] = "FlameMaster"
			#print(self.add)
		else:
			self.add = {}
			
			self.add["solver"] = "FlameMaster"
		
		self.uniqueID = str(self.dataSet_id)+"_"+str(parameters[0].strip("\t"))
		print(self.uniqueID)
		#By default  
#		default = {}
#		cantera_def = {}
#		cantera_def["saveAll"] = '#' 
#		cantera_def["ign_species"] = 'OH'
#		cantera_def["ign_cond"] = 'max'
#		cantera_def["ign_specific_cond"] = "None;None"
#		FlameMan_def = {}
#		default["cantera"] = cantera_def
#		default["FlameMaster"] = FlameMan_def
		self.species_dict = {}
		self.ignition_type = "reflected"
		self.flame = "laminar"
		self.reactor = "JSR" 
		self.std_dvtn = 1.0
		
		for param in parameters:
			#print(param)
			key_and_value = param.split("--")
			
			if len(key_and_value) == 1:
				continue
			
			key = key_and_value[0].strip()
			#print(key)
			content = key_and_value[1].strip()
			
			if "input_file" in key:
				self.input_file = os.path.abspath(content)
				continue
					
			if "target" in key:
				self.target = content
#				print(content)
#				if content.strip(" ") == "ignition delay measurement" or "Tig":
#					self.target = "Tig"
#				
#				elif content.strip(" ") == "laminar burning velocity measurement" or "Fsl" or "Sl":
#					self.target = "Fsl"
#					
#				elif content.strip(" ") == "Heat Flux Burner" or "Flf" or "HeatFlux":
#					self.target = "Flf"
#					
#				elif content.strip(" ") == "laminar burning velocity measurement" or "Fsl" or "Sl":
#					self.target = "Flw"
#					
#				else:
#					raise AssertionError("Invalid target type")
				#print(self.target)
			if "simulation" in key :
				self.simulation = content
			
			if "measurnment_type" in key:
				self.measure_type = content
			
			if "Flame_type" in key:
				self.flame_type = content
			
			if "Ignition_mode" in key:
				self.ig_mode = content.strip("[]")
			
			if "Reactor_type" in key:
				self.reactor = content
			
			if "Residence_time" in key:
				self.residenceTime = content
			
			if "startProfile" in key:
				self.s_p_name = content
				
			if "Fuel_type" in key:
				self.fuel_type = content
	
			if key.strip() == "Fuel":
				#print(content.split("->"))
				#print(content)
				if "Multi" in self.fuel_type:
					self.fuel_id = json.loads(content.split("=")[0].split("->")[1].replace("'",'"'))
					self.fuel_x = json.loads(content.split("=")[1].replace("'",'"'))	
					self.fuel_is = ""
					for i in self.fuel_id:
						self.species_dict[str(self.fuel_id[i])] = float(self.fuel_x[i])
						self.fuel_is+="fuel is "+str(self.fuel_id[i])+"\n"
				else:
					self.fuel_id = content.split("->")[1].split("=")[0]
					self.fuel_x = float(content.split("->")[1].split("=")[1])
					self.species_dict[str(self.fuel_id)] = float(self.fuel_x)
					
			
			if "Oxidizer" in key:
				self.oxidizer = content.split("->")[1].split("=")[0].strip(" ")
				self.oxidizer_x = float(content.split("->")[1].split("=")[1])# to convert to percentage
				self.species_dict[str(self.oxidizer)] = float(self.oxidizer_x)
				actual_oxidizer = self.species_dict[self.oxidizer]*self.molecularWt[self.oxidizer]
			
			
			if "Bath_gas" in key:
				self.bath_gas_id = "Multi"
				self.bath_gas = json.loads(content.split("=")[0].split("->")[1].replace("'",'"'))
				self.bath_gas_x = json.loads(content.split("=")[1].replace("'",'"'))
	
				for i in self.bath_gas:
					self.species_dict[str(self.bath_gas[i])] = float(self.bath_gas_x[i])
						
			if "BG1" in key:
				if content.split("=")[1].strip() != '':
					self.bath_gas1 = content.split("->")[1].split("=")[0]
					self.bath_gas1_x = float(content.split("->")[1].split("=")[1])# to convert to percentage
					self.bg1_tag = ""
				else:
					self.bath_gas1 = ""
					self.bath_gas1_x = ""# to convert to percentage
					self.bg1_tag = "#"
				
			if "BG2" in key:#mole %/100
				if content.split("=")[1].strip() != '':
					self.bath_gas2 = content.split("->")[1].split("=")[0]
					self.bath_gas2_x = float(content.split("->")[1].split("=")[1])# to convert to percentage
					self.bg2_tag = ""
				else:
					self.bath_gas2 = ""
					self.bath_gas2_x = ""# to convert to percentage
					self.bg2_tag = "#"
				
			if "BG3" in key:#mole %/100
				
				if content.split("=")[1].strip() != '':
					self.bath_gas3 = content.split("->")[1].split("=")[0]
					self.bath_gas3_x = float(content.split("->")[1].split("=")[1])# to convert to percentage
					self.bg3_tag = ""
				else:
					self.bath_gas3 = ""
					self.bath_gas3_x = ""# to convert to percentage
					self.bg3_tag = "#"
			
			if  key == "T":#K
				self.temperature = float(content)
			
			if key.strip() == "units":
				self.units = json.loads(content.replace("'",'"'))
			
			if key.strip() == "flow_rate":
				self.flow_rate = content
			
			if key.strip() == "Fuel_units":		
				self.f_unit = content
			
			if key.strip() == "T_units":
				self.t_unit = content
							
			if key.strip() == "P_units":
				self.p_unit = content
			
			if key.strip() == "obs_unit":
				self.target_unit = content
			
			if key.strip() == "flow_units":
				self.flow_unit = content
				
			if key.strip() == "P":
				self.pressure = float(content) #atm
				
			if key.strip() == "Phi":
				self.phi = str(content)
					
			if key.strip() == "observed":
				self.observed = float(content);
				
			if "deviation" in key:
				self.std_dvtn = float(content);
			
			if "dataSet" in key:
				self.d_set = content;#index of a data set
			else:
				self.d_set = self.dataSet_id

			if "data_weight" in key:
				self.d_weight = 1/float(content.strip());#weight of each data point in a perticular data set
			
			if "species" in key:
				self.species = content
				#print(content)
			
			if "time" in key:
				self.Tend = float(content)
				#print(content)	
			
			#print(self.target)		
		
		if self.units["P"].strip(" ") == "atm" or self.units["P"].strip(" ") == "bar":
			self.pressure = float(self.pressure)*101350
		if self.units["P"].strip("	") == "Pa":
			self.pressure = float(self.pressure)
		if "solver" not in self.add:
			self.add["solver"] = "FlameMaster"
		if "Thermodiffusion" not in self.add:
			self.add["Thermodiffusion"] = True
		if "ComputeWithRadiation" not in self.add:
			self.add["ComputeWithRadiation"] = True
		if "saveAll" not in self.add:
			self.add["saveAll"] = "#"
		if "BoundaryLayer" not in self.add:
			self.add["BoundaryLayer"] = False
		if "ign_delay_def" not in self.add:
			self.add["ign_delay_def"] = "OH"
		if "ign_cond" not in self.add:
			self.add["ign_cond"] = "max"
		if "specific_cond" not in self.add:
			self.add["specific_cond"] = "None;None"
		if "rcm_vt_profile" not in self.add:
			self.add["rcm_vt_profile"] = ""
		if "exp_conc" not in self.add:
			self.add["exp_conc"] = {}
			self.add["exp_conc"][float(self.temperature)] = ""
		else:
			if int(self.temperature) not in self.add["exp_conc"]:
				self.add["exp_conc"][float(self.temperature)] = "" 
		if "width" not in self.add:
			self.add["width"] = 0.014 
		if "ratio" not in self.add:
			self.add["ratio"] = 3
		if "slope" not in self.add:
			self.add["slope"] = 0.015
		if "curve" not in self.add:
			self.add["curve"] = 0.015
		if "loglevel" not in self.add:
			self.add["loglevel"] = 1		
		if "auto" not in self.add:
			self.add["auto"] = "True"
		if "width_bsf" not in self.add:
			self.add["width_bsf"] = 0.5 #m 
		if "slope_bsf" not in self.add:
			self.add["slope_bsf"] = 0.05
		if "phi" not in self.add:
			self.add["phi"] = "N/A"
		if "flw_length" not in self.add:
			self.add["flw_length"] = 1 #m 
		if "flow_velocity" not in self.add:
			self.add["flow_velocity"] = 0.06  #m/s 
		if "crossSectionalArea" not in self.add:
			self.add["crossSectionalArea"] = 0.000314 # m**2,pi,r=0.01 m
		if "reactorVolume" not in self.add:
			self.add["reactorVolume"] = 0.000314 #m**3,area*l	
		if "residenceTime" not in self.add:
			self.add["residenceTime"] = 1
		if "total_time" not in self.add:
			self.add["total_time"] = 10 #s
		if "time_step" not in self.add:
			self.add["time_step"] = 2000 #steps		
		if "flw_method" not in self.add:
			self.add["flw_method"] = "slope"
		if "flw_species" not in self.add:
			self.add["flw_species"] = "H2" #species 
		if "flw_limits" not in self.add:
			self.add["flw_limits"] = (40,60) #two points between which slope is to be found 
		if "anchor" not in self.add:
			self.add["anchor"] = 50 # anchor point 
		if "limit_units" not in self.add:
			self.add["limit_units"] = "percentage"
		if "heat" not in self.add:
			self.add["heat"] = ""
		else:
			self.add["heat"] = "HeatTransCoeff is {}".format(self.add["heat"])
		if "flow_rate" not in self.add:
			self.add["flow_rate"] = "0.06"
		if "flf_target" not in self.add:
			self.add["flf_target"] = "H" #m**3	
		if "flf_cond" not in self.add:
			self.add["flf_cond"] = "max" #m**3	
		if "T_amb" not in self.add:
			self.add["T_ambient"] = 298
		else:
			self.add["T_amb"] = "AmbientTemp is {}".format(self.add["T_ambient"])
		if "isIsotherm" not in self.add:
			self.add["isIsotherm"] = ""
		else:
			self.add["isIsotherm"] = "Isothermal is {}".format(self.add["isIsotherm"])
		if "transport_model" not in self.add:
			self.add["transport_model"] = "Multi"
			self.add["solve_bsf"] = "loglevel,refine_grid"
			self.add["group"] = "multi"
			self.add["description"] = "solution with multicomponent transport"
		else:
			if "Mix" in self.add["transport_model"]:
				self.add["solve_bsf"] = "loglevel,auto=True"
				self.add["group"] = "mix"
				self.add["description"] = "solution with mixture-averaged transport"
			
		self.initialCond = self.getInitialCond()
		if "Flf" in self.target:
		
			if "ExpTempFile" not in self.add:
				self.add["ExpTempFile"] = ""
				self.burner_temp = 0
			else:
				# Populating the expTempProfile
				# Using spline function in scipy
				read = open(str(self.add["ExpTempFile"]),'r').readlines()
				data = []
				for i in read:
					if "#" not in i:
						data.append(i)
				header = self.add["file_header"]
				x,y,xx,yy,self.burner_temp = self.getBurnerTemp(data,header)
				string = ""
				for ind,val in enumerate(xx):
					if self.add["solver"] == 'cantera':
						string+="{},{}".format(val,yy[ind])
					else:
						string+="{}\t{}".format(val,yy[ind])
					if ind<len(xx):
						string+="\n"
					else:
						continue
				
				for ind,i in enumerate(x):
					if float(i)>1.0:
						if self.add["solver"] == 'cantera':
							string+="{},{}".format(i,y[ind])
						else:
							string+="{}\t{}".format(i,y[ind])
						if ind<len(x):
							string+="\n"
						else:
							continue
					else:
						continue
				
				self.add["flf_grid"] = x[-1]
				file_dump = open(self.add["ExpTempFile"],"w+")
				file_dump.write(string)
				file_dump.close()
		    
		if "Fsl" in self.target:
			if key == "start_profile":
				self.startProfile_location = content	
		if "N/A" in self.phi:
			self.phi = self.add["phi"]
		
		
		
		actual_fuel = []
		st_fuel = []
		st_ox = []
		for i in self.species_dict:
			actual_fuel.append(self.molecularWt[i]*self.species_dict[i])
			if self.species_dict[i] != 0:
				st_fuel.append(self.molecularWt[i])
				st_ox.append(self.stoichiometry[i]*self.molecularWt['O2'])
		#print(actual_fuel)
		#print(actual_oxidizer)
		Fast = 0
		Fact = 0
		for i,ele in enumerate(st_fuel):
			if st_ox[i] !=0:
				Fast += st_fuel[i]/st_ox[i]	
				Fact += actual_fuel[i]/actual_oxidizer	
		#print(Fast)
		#print(Fact)
		
		self.phi = float(Fact/Fast)
		if "Fls" in self.target:
			print("Target {} has phi of {}\n".format(self.uniqueID,self.phi))
		#self.equivalence_ratio = (F/A)act/(F/A)st
#This function goes through a file previously created by the "generate_target_value_tables" function in data management module 
#and extract target value information for each case and sort them based on the extend of perturbation.
	def getBurnerTemp(self,data,structure):
		x = []
		y = []
		for i in data:
			if len(i)>1:
				for ind,s in enumerate(structure):
					if "distance" in s:
						if self.add["solver"] == 'cantera':
							x.append(float(i.split(",")[ind].strip("\n")))
						else:
							x.append(float(i.split("	")[ind].strip("\n")))
					if "T" in s:
						if self.add["solver"] == 'cantera':
							y.append(float(i.split(",")[ind].strip("\n")))
						else:
							y.append(float(i.split("	")[ind].strip("\n")))	
		x = np.asarray(x)
		y = np.asarray(y)
		xo = 0
		order = 1
		s = InterpolatedUnivariateSpline(x, y, k=order)
		inter = 4
		xx = np.arange(x[0],1,0.02)
		ss = InterpolatedUnivariateSpline(x, y, k=inter)
		yy = ss(xx)
		yo = s(xo)
		return x,y,xx,yy,float(yo)
		
	def getInitialCond(self):
		species = self.species_dict
		str_syntax = "x->"
		dirichlet = ""
		for i in species:
			dirichlet += str_syntax+str(i)+"="+str(species[i])+"\n\t\t"
		return dirichlet
	
	def make_eta_lists(self,reaction_index,unsrt):
		xdata = []
		self.xdata = []		
		self.ydata = []
		home_dir = os.getcwd()
		os.chdir(self.case_index)
		eta_values = open(home_dir+'/Data/Simulations/sim_data_'+self.case_index+'.lst','r').readlines()
		for eta in eta_values:
			params = eta.split("\t")
			self.ydata.append(float(params[1]))
		
		#beta = open(home_dir+'/Data/Simulations/Generator_list_'+self.case_index+'.csv','r').readlines()		
		beta = open(home_dir+'/Data/Simulations/Beta_list_'+self.case_index+'.csv','r').readlines()		
		#beta = open(home_dir+'/Data/Simulations/NormalizedBeta_list_'+self.case_index+'.csv','r').readlines()	
		for i in beta[0:len(self.ydata)]:
			X = []			
			line = i.split(',')
			for j in line:
				if j != '\n':				
					X.append(float(j))
			self.xdata.append(X)	
		#print(np.shape(self.xdata))
		os.chdir(home_dir)			

	def MatPolyFitTransform(self,order):
		BTrsMatrix = []
		for i in self.xdata:
			tow = order
			row = i
			row_ = []
			row_.append(1)
			if tow > 0:		
				for i in row:				
					row_.append(i)
				tow = tow - 1

			if tow > 0:
				for i,j in enumerate(row): 
					for k in row[i:]:					
						row_.append(j*k)
				tow = tow - 1

			if tow > 0:		
				for i in row:
					for j in row[row.index(i):]: 
						for k in row[row.index(j):]:					
							row_.append(i*j*k)
				tow = tow - 1					

			if tow > 0:
				for i,j in enumerate(row):
					for j in row[row.index(i):]: 
						for k in row[row.index(j):]:					
							for l in row[row.index(k):]:
								row_.append(i*j*k*l)
				tow = tow - 1

			if tow > 0:
				for i in row:
					for j in row[row.index(i):]: 
						for k in row[row.index(j):]:					
							for l in row[row.index(k):]:
								for m in row[row.index(l):]:
									row_.append(i*j*k*l*m)
				tow = tow - 1
			BTrsMatrix.append(row_)
		return BTrsMatrix

	def evaluate(self,BZeta,order):
		tow = order
		row = BZeta
		#print(len(row))
		coeff = self.resCoef
		#print(len(coeff))
		#print(len(row))
		#if self.PRS_type == "Partial":
		#	if self.activeIndexDict is not None:
		#		row = []
		#		for i,ele in enumerate(self.activeIndexDict[self.index]):
		#			if self.activeIndexDict[self.index][ele] == 1:
		#				row.append(BZeta[i])
				#coeff = self.resCoef_full
		val = coeff[0]
		count = 1
		if tow > 0:		
			for i in row:				
				val+=coeff[count]*i
				count +=1
			tow = tow - 1

		if tow > 0:
			for i,j in enumerate(row): 
				for k in row[i:]:
					if count < len(coeff):					
						val+=coeff[count]*j*k
						count +=1
			tow = tow - 1

		if tow > 0:		
			for i,j in enumerate(row):
				for k in row[row.index(j):]: 
					for m in row[row.index(k):]:
						if count < len(coeff):					
							val+=coeff[count]*j*k*m
							count +=1
			tow = tow - 1					

		if tow > 0:
			for i,j in enumerate(row):
				for k,l in enumerate(row[i:]): 
					for m,n in enumerate(row[k:]):					
						for o in row[m:]:
							if count < len(coeff):
								val+=coeff[count]*j*l*n*o
								count +=1
			tow = tow - 1

		if tow > 0:
			for i,j in enumerate(row):
				for k,l in enumerate(row[i:]): 
					for m,n in enumerate(row[k:]):					
						for o,p in enumerate(row[m:]):
							for q in row[o:]:
								if count<len(coeff):
									val+=coeff[count]*q*p*n*l*j
									count +=1
			tow = tow - 1
			
		return val
	
	def calculated_target_value(self,BZeta):	
		if self.PRS_type == "Partial":
			if self.activeIndexDict is not None:
				row = []
				for i,ele in enumerate(self.activeIndexDict[self.index]):
					if self.activeIndexDict[self.index][ele] == 1:
						row.append(BZeta[i])
				#coeff = self.resCoef_full
				
		else:
			row = BZeta
		tow = self.order
		coeff = self.resCoef
		val = coeff[0]
		count = 1
		if tow > 0:		
			for i in BZeta:				
				val+=coeff[count]*i
				count +=1
			tow = tow - 1

		if tow > 0:
			for i,j in enumerate(BZeta): 
				for k in BZeta[i:]:
					if count < len(coeff):					
						val+=coeff[count]*j*k
						count +=1
			tow = tow - 1

		if tow > 0:		
			for i,j in enumerate(BZeta):
				for k in BZeta[BZeta.index(j):]: 
					for m in row[row.index(k):]:
						if count < len(coeff):					
							val+=coeff[count]*j*k*m
							count +=1
			tow = tow - 1					

		if tow > 0:
			for i,j in enumerate(BZeta):
				for k,l in enumerate(BZeta[i:]): 
					for m,n in enumerate(BZeta[k:]):					
						for o in BZeta[m:]:
							if count < len(coeff):
								val+=coeff[count]*j*l*n*o
								count +=1
			tow = tow - 1

		if tow > 0:
			for i,j in enumerate(BZeta):
				for k,l in enumerate(BZeta[i:]): 
					for m,n in enumerate(BZeta[k:]):					
						for o,p in enumerate(BZeta[m:]):
							for q in BZeta[o:]:
								if count<len(coeff):
									val+=coeff[count]*q*p*n*l*j
									count +=1
			tow = tow - 1
			
		return val
		
	def evaluateResponse(self,x,cov_x=None):
		#print(len(x))
		if self.PRS_type == "Partial":
			if self.activeIndexDict is not None:
				count = 0
				x_ = []
				for i in self.activeIndexDict[self.index]:
					if self.activeIndexDict[self.index][i] == 1:
						x_.append(x[count])
						count+=1
					else:
						continue
						count+=1 
			x = np.asarray(x_)
		val = self.zero
		val += np.dot(self.a,x)
		#print(len(self.a))
		#print(len(x))
		if cov_x is not None:
			a_times_cov = np.dot(self.a,cov_x)
			variance = np.dot(self.a,a_times_cov.T)
        
		if self.b is not None:
			b_times_x = np.asarray(np.dot(self.b,x)).flatten()
			val += np.dot(b_times_x.T,x)
			if cov_x is not None:
				b_times_cov = np.dot(self.b,cov_x)
				variance += 2*np.trace(np.dot(b_times_cov,b_times_cov))
         
		if cov_x is not None:
			computed_unc = math.sqrt(variance)
			return val,computed_unc
		return val
		 
	def estimate(self,x):
		val = self.zero
		val += np.dot(self.a,x)
		
		response_grad = self.a
		
		if not(self.b is None):
			b_times_x = np.asarray(np.dot(self.b,x)).flatten()
			val += np.dot(b_times_x.T,x)
			response_grad += 2*b_times_x

		return val,response_grad
		
	def resCoeffTransform(self,order):
		x = self.x
		coeff = self.resCoef
		#if self.PRS_type == "Partial":
		#	if self.activeIndexDict is not None:
		#		coeff = self.resCoef_full
		tow = order
		row = self.xdata[0]
		row_transformed = row
		#if self.PRS_type == "Partial":
		#	if self.activeIndexDict is not None:
		#		count = 0
		#		for i in self.activeIndexDict[self.index]:
		#			if self.activeIndexDict[self.index][i] == 1:
		#				row_transformed.append(row[count])
		#				count+=1
		#			else:
		#				row_transformed.append(0)
		#else:
		#	row_transformed = row
		zero = []
		a = []
		b = []
		zero.append(coeff[0])
		count = 1
		if tow > 0:		
			for _ in row_transformed:				
				a.append(coeff[count])
				count +=1
			tow = tow - 1

		if tow > 0:
			for i,j in enumerate(row_transformed): 
				temp = []
				for k,l in enumerate(row_transformed[i:]):
					if count < len(coeff):				
						temp.append(coeff[count])
						count +=1
				b.append(list(np.zeros(len(row_transformed)-len(temp)))+temp)					
			tow = tow - 1
		return float(x[0]),np.asarray(a),np.matrix(b)
	
	def getFullResCoef(self,BZeta):
		coeff = self.resCoef
		tow = self.order
		row = list(BZeta)
		#print(len(row))
		bzeta_transformed = []
		count = 0
		for i in self.activeIndexDict[self.index]:
			if self.activeIndexDict[self.index][i] == 1:
				bzeta_transformed.append(row[count])
				count+=1
			else:
				bzeta_transformed.append(0)
		#print(count)
		new_coefficients = []
		val = coeff[0]
		count = 1
		new_coefficients.append(val)
		if tow > 0:		
			for i in bzeta_transformed:				
				if i == 0:
					new_coefficients.append(0)
				else:
					new_coefficients.append(coeff[count])
					count+=1
			tow = tow - 1

		if tow > 0:
			for i,j in enumerate(bzeta_transformed): 
				for k in bzeta_transformed[i:]:					
					temp = j*k
					if temp == 0:
						new_coefficients.append(0)
					else:
						new_coefficients.append(coeff[count])
						count+=1

			tow = tow - 1
		##After order 2, response surface is not developed
		if tow > 0:		
			for i in row:
				for j in row[row.index(i):]: 
					for k in row[row.index(j):]:					
						row_.append(i*j*k)
			tow = tow - 1					

		if tow > 0:
			for i,j in enumerate(row):
				for j in row[row.index(i):]: 
					for k in row[row.index(j):]:					
						for l in row[row.index(k):]:
							row_.append(i*j*k*l)
			tow = tow - 1

		if tow > 0:
			for i in row:
				for j in row[row.index(i):]: 
					for k in row[row.index(j):]:					
						for l in row[row.index(k):]:
							for m in row[row.index(l):]:
								row_.append(i*j*k*l*m)
			tow = tow - 1
		return new_coefficients
	
	
	
	def model_response_uncertainty(self,x,cov,order):
		coeff = self.resCoef
		tow = order
		val = 0
		a = 0
		b_ii = 0
		b_ij = 0
		row = x
		count = 1
		if tow > 0:		
			for i in BZeta:				
				a += coeff[count]**2
				count +=1
			tow = tow - 1

		if tow > 0:
			for i,j in enumerate(BZeta): 
				for k,l in enumerate(BZeta[i:]):
					if count < len(coeff):
						if i == k:
							b_ii+=coeff[count]**2
						else:
							b_ij+=coeff[count]**2					
						count +=1
			tow = tow - 1
		
		#print(a)
		#print(b_ii)
		#print(b_ij)
		#val = a + 2*b_ii+b_ij
		
		val = a + 2*b_ii+b_ij
		
		#val = np.linalg.norm(np.dot(np.asarray(a),cov))**2+ np.linalg.norm(np.dot(cov.T,np.dot(np.matrix(b),cov)),'fro')**2
		### Calculating the variance for each response surface
		#print("Len of x mdres is {}".format(len(self.resCoef)))
		#print("Count is {}\n".format(count))	
		return np.sqrt(val)
	
	def jacobian_element(self,coeff,BZeta,x,ind,resp_order):
		"""
		J = a + bx
		"""
		tow = resp_order
		row = BZeta   
		val = 0
		count = 1
		index=[]
		"""
		find a
		"""
		if tow > 0:
			for i,j in enumerate(BZeta):
				
				if i == ind:
					val+=coeff[count]
					index.append(count)
					count +=1
				
				else:
					count +=1
			tow = tow - 1

		"""
		find b*x
		"""
		if tow > 0:
			for i,j in enumerate(BZeta):
				l = i 
				for k in BZeta[i:]:
					#if count < len(coeff):
					if i == ind and l == ind:
						val += 2*coeff[count]*k
						index.append(count)
						#print("{}{}\t{}\n".format(i,l,count))
						count+=1
						l+=1
					elif i == ind and l!=ind:
						val += coeff[count]*k
						index.append(count)
						#print("{}{}\t{}\n".format(i,BZeta.index(k),count))
						count+=1
						l+=1
					elif i!= ind and l==ind:
						val += coeff[count]*j
		#				print("{}{}\t{}\n".format(i,BZeta.index(k),count))
						index.append(count)
						count+=1
						l+=1
					else:
						#print("uncounted {}{}\t{}\n".format(i,BZeta.index(k),count))
						count+=1
						l+=1
			tow = tow - 1
		#print(val)
		#print("\t\tThe index for {} in x[{}]\n{}\n".format(x,ind,index))
		return val
	
	def getJacobian(self,x):
		"""
		Jacobian = a + 2bx
		
		"""
		J = np.asarray(self.a)+2*np.asarray(self.b.dot(x))
		return J
		
	def Jacobian(self,x,resp_order):
		
		j = []
		x = list(x)
		#print("\tx{} (type-{}) in response surface is\n".format(x,type(x)))
		for i,opt in enumerate(x):
			j.append(self.jacobian_element(self.resCoef,x,opt,i,resp_order))
		return j 
	
#This function creates a set of coefficients from the calculated eta values collected by the function make_eta _lists of combustion_target_class
#These co-efficients define the response surface for each target hence they're connected to the respective instance of the class.	
	
	#def getSAB_coefficients(self,reaction_index,unsrt,order):
	def maxColumn(self,L):    
    		return list(map(max, zip(*L)))
	def objective(self,z):
		"""
		(Ax-y) + penalty = 0
		"""
		A = self.BTrsMatrix
		prediction = A.dot(z)
		simulation = self.actualValue
		residual = (prediction - simulation)
		obj = np.dot(residual,residual.T)
		for i in z:
			obj+=i**2
		return obj
		
	def create_response_surface(self,reaction_index,unsrt,order,selected=None,activeParams=None,PRS_type = "Full"):
		self.PRS_type = PRS_type
		self.selectedParams = selected
		self.activeIndexDict = activeParams
		self.make_eta_lists(reaction_index,unsrt)			
		self.order = order
		##Write code for scipy curve fit
		X = np.matrix(self.xdata)
		bnd = np.asarray(X)
		df = pd.DataFrame(bnd)
		bounds = np.asarray(df.abs().max(axis=0))
		init_guess = bnd[0]
		#print(np.shape(X))
		
		Y = np.array(self.ydata)
		#print(np.shape(Y))
		BTrsMatrix = self.MatPolyFitTransform(int(order))
		print(np.shape(BTrsMatrix))
		self.BTrsMatrix = np.matrix(BTrsMatrix)
		#print(np.matrix(BTrsMatrix))
		self.actualValue = Y
 		#print(len(selected))
		#print(len(X[0]))
		stringX = ""
		#Singular Value analysis		
		#u, s, vh = np.linalg.svd(BTrsMatrix, full_matrices=True)
		#sv = open(os.getcwd()+'/Data/NumericalAnalysis/Singular_values_'+self.case_index+'.csv','w')
		#string = ""
		#for i in s:
		#	string+="{}\n".format(float(i))
		#sv.write(string)
		#sv.close()
		#S = s/s[0]
		significant_singular_values = len(self.ydata)
		#for i in range(len(S)):
		#	if S[i] < 10**(-3):
		#		significant_singular_values = i
		#	break
		
		#PLOTTING SINGULAR VALUES
		#x = np.linspace(0,len(s),len(s))
		#fig = plt.figure()
		#plt.title('Singular values of the Transformed Matrix\n for Response Surface generation with 18 Rxn\n(Values from Largest to smallest) for case'+self.case_index)
		#plt.plot(x,np.log(s),'-',label = 'Singular Values')
		#plt.xlabel('Index')
		#plt.ylabel('logarithmic Singular Values')
		#plt.legend()
		#plt.savefig(os.getcwd()+'/Plots/SingularValues/Logarithmic SingularValue_18Rnx'+self.case_index+'.png')

			
		#SOLVING FOR THE CO-EFFICIENTS
		start = os.getcwd()
		#coef_dir = str(start)+'/Data/ResponseSurface/'
		#os.chdir(coef_dir)
		#coef_file = "responsecoef_"+str(self.case_index)+".csv"
		"""
		if coef_file in os.listdir():
			#print("PRS already generated")
			x_file = open(coef_file,"r").readlines()
			x = []
			for i in x_file[1:]:
				x.append(float(i.strip("\n").strip("")))
			self.x = np.asarray(x)
			os.chdir(start)
			norm_Ax_b_0 = 1
			norm_b_0 = 1
			norm_Ax_b_1 = 1
			norm_b_1 = 1
			norm_Ax_b_2 = 1
			norm_b_2 = 1
		else:
		"""
		"""
		Doing least square using QR decomposition
		"""
		
		Q, R = np.linalg.qr(BTrsMatrix)
	#print(self.ydata)
		y = np.dot(np.transpose(Q),self.ydata)
		self.x = np.linalg.solve(R,np.transpose(y))
		
		"""
		coeff_size = np.asarray(self.BTrsMatrix[0]).flatten()
		guess = np.ones(len(coeff_size))
		strat = time.time()
		x = minimize(self.objective,guess,method="SLSQP")
		stop = time.time()
		print(f"time taken is {stop-strat}")
		print(x)
		self.x = x.x
		"""
		"""
		Using the penalty function to fit the response surface
		"""		
	
		#norm_Ax_b_0 = np.linalg.norm(abs(np.dot(R,self.x)-y),ord = 0)
		#norm_b_0 = np.linalg.norm(abs(y),ord = 0)
		#norm_Ax_b_1 = np.linalg.norm(abs(np.dot(R,self.x)-y),ord = 1)
		#norm_b_1 = np.linalg.norm(abs(y),ord = 1)
		#norm_Ax_b_2 = np.linalg.norm(abs(np.dot(R,self.x)-y),ord = 2)
		#norm_b_2 = np.linalg.norm(abs(y),ord = 2)
		#norm_Ax_b_inf = np.linalg.norm(np.dot(R,x)-y,ord = 'inf',axis = 1)
		#norm_b_inf = np.linalg.norm(y,ord = 'inf',axis = 1)
		
		# Populating the a vector and b matrix:


		#Generating Response Co-efficients
		self.resCoef = []
		rr = open(os.getcwd()+'/Data/ResponseSurface/responsecoef_'+self.case_index+'.csv','w')		
		res = "Coefficients\n"
		for i in self.x:
			self.resCoef.append(float(i))
			res +='{}\n'.format(float(i))
		rr.write(res)
		rr.close()
		self.resFramWrk = []		
		#print(self.resCoef)
		#self.zero,self.a,self.b = self.resCoeffTransform(2)
		if self.PRS_type == "Partial":
			if activeParams is not None:
				self.resCoef_full = np.asarray(self.getFullResCoef(self.xdata[0]))
				#print(self.resCoef_full)
		
		
		self.zero,self.a,self.b = self.resCoeffTransform(2)
		

		#simVSresp  = "FlameMaster,Response Surface\n"
		for i in self.xdata:
		    self.resFramWrk.append(self.evaluate(i,int(order)))
		fileResponse = open(os.getcwd()+'/Data/ResponseSurface/FlaMan_Response_comparison_'+self.case_index+'.csv','w')
		simVSresp  = "FlameMaster,Response Surface,Error_(ln(tau)),Error(Tau)\n"
		self.RSM_error = []
		self.RSM_error_relative = []
		TraError = []
		for i in range(len(self.resFramWrk)):
			self.RSM_error.append(abs(self.ydata[i]-self.resFramWrk[i]))
			self.RSM_error_relative.append((abs(self.ydata[i]-self.resFramWrk[i])/self.ydata[i])*100)
			TraError.append(1-np.exp(-(Y[i]-self.resFramWrk[i])))
			simVSresp +='{},{},{},{}\n'.format(Y[i],self.resFramWrk[i],(self.ydata[i]-self.resFramWrk[i])/self.ydata[i],1-np.exp(-(Y[i]-self.resFramWrk[i])))
		    #print('{},{},{}\n'.format(Y[i],self.resFramWrk[i],100*(Y[i]-self.resFramWrk[i])/Y[i]))

		
		
		#Simulation report		
		fileResponse.write(simVSresp)		
		fileResponse.close() 
		ff = open(os.getcwd()+'/Data/ResponseSurface/report_'+self.case_index+'.csv','w')
		self.MaxError = max(self.RSM_error)
		self.MeanError = statistics.mean(self.RSM_error)
		
		self.MaxError_relative = max(self.RSM_error_relative)
		self.MeanError_relative = statistics.mean(self.RSM_error_relative)
		
		
		MaxTauError = max(np.abs(TraError))
		string = 'Significant singular values\n'
		string += '{}\n'.format(significant_singular_values)
		string = 'Max relative error in ln(Tau)\n'
		string += '{}\n'.format(self.MaxError)
		string = 'Max % relative error in ln(Tau)\n'
		string += '{}\n'.format(self.MaxError*100)
		string = 'Max relative error in Tau\n'
		string += '{}\n'.format(MaxTauError)
		string = 'Max % relative error in Tau\n'
		string += '{}\n'.format(MaxTauError*100)
		string += 'Conditioning Number of BZetaMatrix\n'
		string += '{}\n'.format(LA.cond(self.xdata))
		string += 'Conditioning Number of BTrsMatrix\n'
		string += '{}\n'.format(LA.cond(BTrsMatrix))
		string += 'Relative error based on 0-norm for ln(Tau)\n,'
		#string += '{}\n'.format(abs((norm_Ax_b_0-norm_b_0))/norm_b_0)
		#string += 'Relative error based on 1-norm ln(Tau)\n'
		#string += '{}\n'.format(abs(norm_Ax_b_1-norm_b_1)/norm_b_1)
		#string += 'Relative error based on 2-norm ln(Tau)\n'
		#string += '{}\n'.format(abs(norm_Ax_b_2-norm_b_2)/norm_b_2)
		#string += 'Relative error based on inf-norm ln(Tau)\n'
		#string += '{}\n'.format((norm_Ax_b_inf-norm_b_inf)/norm_b_inf)
		ff.write(string)
		ff.close()
		#print("---------------------------------------------")
		#print("Conditioning number of the response surface matrix X is {} \n".format(LA.cond(self.xdata)))
		#print("Conditioning number of the response surface matrix X is {} \n".format(LA.cond(BTrsMatrix)))
		#print("---------------------------------------------")		
		return bounds,init_guess
		
#This function calculates the target value for any set of normalised vectiors, values of all x must be less than 1 and greater than -1. Function returns a second order polynomial. 		
		
	def calculated_target_value(self,x):
		#print(x)
		x_ = []
		if self.PRS_type == "Partial":
			if self.activeIndexDict is not None:
				for i,ele in enumerate(self.activeIndexDict):
					if self.activeIndexDict[ele] == 1:
						x_.append(x[i])
					
		else:
			x_= x
		target_value = self.resCoef[0]
		count = 1
		#print(len(x))
		#print(len(self.resCoef))
		for i in x_:
			target_value += self.resCoef[count]*i
			count +=1
		for i in x_:
			for j in x[x.tolist().index(i):]:	
				if count<=len(self.resCoef)-1:		
					target_value += self.resCoef[count]*i*j
					count +=1
		#print("Len of x is {}".format(len(self.resCoef)))
		#print("Count is {}\n".format(count))		
		return target_value
								
			
		
