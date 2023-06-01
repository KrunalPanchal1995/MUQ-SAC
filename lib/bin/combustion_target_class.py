import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import PolynomialFeatures
from sklearn import linear_model
from numpy import linalg as LA
import json

class combustion_target():
###Class definition and acquisition of input parameters.	
	def __init__(self,data,addendum,index):
		self.data = data
		parameters = self.data.split('|')
		self.calculated = 0
		self.x = None
		self.case_index = "case-"+str(index)
		self.units = self.f_unit = self.t_unit = self.temperature_factor = self.p_unit = self.pressure_factor = self.target_unit = self.target_factor = self.flow_unit  = self.flow_rate = self.target_key = self.input_file = self.target = self.simulation = self.temperature = self.fuel = self.oxygen = self.nitrogen = self.argon = self.pressure = self.phi = self.observed = self.d_weight = self.d_set = self.Tend = self.species = self.s_p_name = None
		
		
		#By default  
		default = {}
		cantera_def = {}
		cantera_def["saveAll"] = '#' 
		cantera_def["ign_species"] = 'OH'
		cantera_def["ign_cond"] = 'max'
		cantera_def["ign_specific_cond"] = "None;None"
		FlameMan_def = {}
		
		
		self.ignition_type = "reflected"
		self.flame = "laminar"
		self.reactor = "JSR" 
		self.std_dvtn = 1.0
		
		for param in parameters:
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
				self.target_key = content
				if content == "ignition delay measurement" or "Tig":
					self.target = "Tig"
				
				elif content == "laminar burning velocity measurement" or "Fsl" or "Sl":
					self.target = "Fsl"
			
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
					
				else:
					self.fuel_id = content.split("->")[1].split("=")[0]
					self.fuel_x = float(content.split("->")[1].split("=")[1])

			if "Oxidizer" in key:
				self.oxidizer = content.split("->")[1].split("=")[0]
				self.oxidizer_x = float(content.split("->")[1].split("=")[1])# to convert to percentage
			if "Bath_gas" in key:
				self.bath_gas = json.loads(content.split("=")[0].split("->")[1].replace("'",'"'))
				self.bath_gas_x = json.loads(content.split("=")[1].replace("'",'"'))
			
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
				self.pressure = float(content) #bar
				
			if key.strip() == "Phi":
				self.phi = content
					
					
			if key.strip() == "observed":
				self.observed = float(content);
				
			if "deviation" in key:
				self.std_dvtn = float(content);
			
			if "dataSet" in key:
				self.d_set = content;#index of a data set

			if "data_weight" in key:
				self.d_weight = 1/float(content.strip());#weight of each data point in a perticular data set
			
			if "species" in key:
				self.species = content
				#print(content)
			
			if "time" in key:
				self.Tend = float(content)
				#print(content)	
			
			#print(self.target)
		if "Fsl" in self.target:
			if key == "start_profile":
				self.startProfile_location = content

		self.species_dict = {}
		if self.
		
#This function goes through a file previously created by the "generate_target_value_tables" function in data management module 
#and extract target value information for each case and sort them based on the extend of perturbation.
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
		
		beta = open(home_dir+'/Data/Simulations/Beta_list_'+self.case_index+'.csv','r').readlines()		
		for i in beta:
			X = []			
			line = i.split(',')
			for j in line:
				if j != '\n':				
					X.append(float(j))
			self.xdata.append(X)	
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

	
	def evaluate(self,coeff, BZeta,order):
		tow = order
		row = BZeta
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
				for k in BZeta[BZeta.index(i):]: 
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

	def model_response_uncertainty(self,BZeta,cov,order):
		coeff = self.resCoef
		tow = order
		row = BZeta
		val = 0
		a = []
		b = []
		count = 1
		if tow > 0:		
			for i in BZeta:				
				a.append(coeff[count])
				count +=1
			tow = tow - 1

		if tow > 0:
			for i,j in enumerate(BZeta): 
				row = []
				for k,l in enumerate(BZeta):
				
					#if count <= len(coeff)-1:
					if k<i:
						row.append(float(0))
					else:
						row.append(coeff[count])
						count+=1
				b.append(row)
					
						#if i == k:						
						#	val+=2*coeff[count]**2
						#	count +=1
							#print(count)
						#else:	
						#	val+=coeff[count]**2
						#	count +=1
							#print(count)
			tow = tow - 1
		val = np.linalg.norm(np.dot(np.asarray(a),cov))**2+2*np.linalg.norm(np.dot(cov.T,np.dot(np.matrix(b),cov)),'fro')**2
		#print("Len of x mdres is {}".format(len(self.resCoef)))
		#print("Count is {}\n".format(count))	
		return val
	def jacobian_element(self,coeff,BZeta,x,ind,resp_order):
		tow = resp_order
		row = BZeta   
		val = 0
		count = 1
		index=[]
		
		if tow > 0:
			for i,j in enumerate(BZeta):
				
				if i == ind:
					val+=coeff[count]
					index.append(count)
				#	print("{}\t{}\n".format(i,count))
					count +=1
				
				else:
					count +=1
			tow = tow - 1

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

	def Jacobian(self,x,resp_order):
		j = []
		x = list(x)
		#print("\tx{} (type-{}) in response surface is\n".format(x,type(x)))
		for i,opt in enumerate(x):
			j.append(self.jacobian_element(self.resCoef,x,opt,i,resp_order))
		return j 
		
#This function creates a set of coefficients from the calculated eta values collected by the function make_eta _lists of combustion_target_class
#These co-efficients define the response surface for each target hence they're connected to the respective instance of the class.	
	
	
	def create_response_surface(self,reaction_index,unsrt,order):
		self.make_eta_lists(reaction_index,unsrt)			
		##Write code for scipy curve fit
		X = np.matrix(self.xdata)
		Y = np.array(self.ydata)
		BTrsMatrix = self.MatPolyFitTransform(int(order))
		
		stringX = ""
		#Singular Value analysis		
		u, s, vh = np.linalg.svd(BTrsMatrix, full_matrices=True)
		sv = open(os.getcwd()+'/Data/NumericalAnalysis/Singular_values_'+self.case_index+'.csv','w')
		string = ""
		for i in s:
			string+="{}\n".format(float(i))
		sv.write(string)
		sv.close()
		S = s/s[0]
		significant_singular_values = len(self.ydata)
		for i in range(len(S)):
			if S[i] < 10**(-3):
				significant_singular_values = i
			break
		
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
		Q, R = np.linalg.qr(BTrsMatrix)
		#print(self.ydata)
		y = np.dot(np.transpose(Q),self.ydata)
		self.x = np.linalg.solve(R,np.transpose(y))
		norm_Ax_b_0 = np.linalg.norm(abs(np.dot(R,self.x)-y),ord = 0)
		norm_b_0 = np.linalg.norm(abs(y),ord = 0)
		norm_Ax_b_1 = np.linalg.norm(abs(np.dot(R,self.x)-y),ord = 1)
		norm_b_1 = np.linalg.norm(abs(y),ord = 1)
		norm_Ax_b_2 = np.linalg.norm(abs(np.dot(R,self.x)-y),ord = 2)
		norm_b_2 = np.linalg.norm(abs(y),ord = 2)
		#norm_Ax_b_inf = np.linalg.norm(np.dot(R,x)-y,ord = 'inf',axis = 1)
		#norm_b_inf = np.linalg.norm(y,ord = 'inf',axis = 1)
		



		#Generating Response Co-efficients
		self.resCoef = []
		rr = open(os.getcwd()+'/Data/ResponseSurface/responsecoef_'+self.case_index+'.csv','w')		
		res = "Coefficients\n"
		for i in self.x:
			self.resCoef.append(float(i))
			res +='{}\n'.format(float(i))
		rr.write(res)
		rr.close()
		resFramWrk = []		

		#simVSresp  = "FlameMaster,Response Surface\n"
		for i in self.xdata:
		    resFramWrk.append(self.evaluate(self.resCoef,i,int(order)))
		fileResponse = open(os.getcwd()+'/Data/ResponseSurface/FlaMan_Response_comparison_'+self.case_index+'.csv','w')
		simVSresp  = "FlameMaster,Response Surface,Error_(ln(tau)),Error(Tau)\n"
		error = []
		TraError = []
		for i in range(len(resFramWrk)):
			error.append((self.ydata[i]-resFramWrk[i])/Y[i])
			TraError.append(1-np.exp(-(Y[i]-resFramWrk[i])))
			simVSresp +='{},{},{},{}\n'.format(Y[i],resFramWrk[i],(self.ydata[i]-resFramWrk[i])/self.ydata[i],1-np.exp(-(Y[i]-resFramWrk[i])))
		    #print('{},{},{}\n'.format(Y[i],resFramWrk[i],100*(Y[i]-resFramWrk[i])/Y[i]))


		
		#Simulation report		
		fileResponse.write(simVSresp)		
		fileResponse.close() 
		ff = open(os.getcwd()+'/Data/ResponseSurface/report_'+self.case_index+'.csv','w')
		self.MaxError = max(np.abs(error))
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
		string += '{}\n'.format(abs((norm_Ax_b_0-norm_b_0))/norm_b_0)
		string += 'Relative error based on 1-norm ln(Tau)\n'
		string += '{}\n'.format(abs(norm_Ax_b_1-norm_b_1)/norm_b_1)
		string += 'Relative error based on 2-norm ln(Tau)\n'
		string += '{}\n'.format(abs(norm_Ax_b_2-norm_b_2)/norm_b_2)
		#string += 'Relative error based on inf-norm ln(Tau)\n'
		#string += '{}\n'.format((norm_Ax_b_inf-norm_b_inf)/norm_b_inf)
		ff.write(string)
		ff.close()
		#print("---------------------------------------------")
		#print("Conditioning number of the response surface matrix X is {} \n".format(LA.cond(self.xdata)))
		#print("Conditioning number of the response surface matrix X is {} \n".format(LA.cond(BTrsMatrix)))
		#print("---------------------------------------------")		
		
		
#This function calculates the target value for any set of normalised vectiors, values of all x must be less than 1 and greater than -1. Function returns a second order polynomial. 		
		
	def calculated_target_value(self,x):
		#print(x)
		target_value = self.resCoef[0]
		count = 1
		for i in x:
			target_value += self.resCoef[count]*i
			count +=1
		for i in x:
			for j in x[x.tolist().index(i):]:	
				if count<=len(self.resCoef)-1:		
					target_value += self.resCoef[count]*i*j
					count +=1
		#print("Len of x is {}".format(len(self.resCoef)))
		#print("Count is {}\n".format(count))		
		return target_value
								
			
		
