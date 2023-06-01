import xml.etree.ElementTree as ET
import scipy as sp
import numpy as np
from scipy.optimize import minimize
from scipy.optimize import curve_fit
import os,re
import Input_file_reader as IFR
import matplotlib.pyplot as plt
from scipy.linalg import block_diag

#############################################################
###       Uncertainty for arrhenius parameters         ######
###       of elementary reactions                      ######
#############################################################

class reaction:
	def __init__(self, Element,mechPath):
		self.rxn = self.classification = self.type = self.sub_type = self.exp_data_type = self.temperatures = self.uncertainties = self.branching  = self.branches = self.pressure_limit = self.common_temp = self.temp_limit = None
		self.tag = Element.tag
		self.rxn = str(Element.attrib["rxn"])
		self.rxn_dict = IFR.MechParsing(mechPath).getKappa(self.rxn)
		for item in Element:
			if item.tag == "class":
				self.classification = item.text
			if item.tag == "type":
				self.type = item.text
			if item.tag == "sub_type":
				self.sub_type = item.attrib["name"]
				for subitem in item:
					if subitem.tag == "multiple":
						self.branching = subitem.text
					if subitem.tag == "branches":
						self.branches = subitem.text
					if subitem.tag == "pressure_limit":
						self.pressure_limit = subitem.text
					if subitem.tag == "common_temp":
						self.common_temp = subitem.text
					if subitem.tag == "temp_limit":
						self.temp_limit = subitem.text
			
			if item.tag == "data_type":
				self.exp_data_type = item.text
			if item.tag == "temp":
				#print(item.text)
				self.temperatures = np.asarray([float(i) for i in item.text.split(",")])
			if item.tag == "unsrt":
				self.uncertainties = np.asarray([float(i) for i in item.text.split(",")])
			
		if self.exp_data_type.split(";")[0] == "constant":
			if self.exp_data_type.split(";")[1] == "array":	
				self.temperatures = self.temperatures
				self.uncertainties = self.uncertainties
				
			elif self.exp_data_type.split(";")[1] == "end_points":
				self.temperatures = np.linspace(self.temperatures[0],self.temperatures[1],20)
				self.uncertainties = np.linspace(self.uncertainties[0],self.uncertainties[1],20)
		
		
		if len(self.temperatures) != len(self.uncertainties):
			print("Error in unsrt data for {}".format(self.rxn))
	
		
		if self.type == "pressure_dependent" and self.sub_type.strip() != "":
			self.nametag = self.rxn+":"+self.sub_type
		else:
			self.nametag = self.rxn
		uncertPriorGuess = np.array([-100.0,-100.0,0.5,20.0,10.0,5.0]);
		
		#Solving using scipy.minimize
		self.solution = minimize(self.uncertPriorObjective,uncertPriorGuess);
		
		#this module can be added when u require it
		L1 = np.array([[self.solution.x[0],0,0],[self.solution.x[1],self.solution.x[2],0],[self.solution.x[3],self.solution.x[4],self.solution.x[5]]]);#cholesky lower triangular matric
		self.L11 = L1[0][0];
		self.L21 = L1[1][0];
		self.L22 = L1[1][1];
		self.L31 = L1[2][0];
		self.L32 = L1[2][1];
		self.L33 = L1[2][2];
		self.Sigma_p = np.dot(L1,np.transpose(L1));#L dot L_transpose
		self.sigma_alpha = np.sqrt(self.Sigma_p[0][0]);# std. deviation of alpha
		self.sigma_n = np.sqrt(self.Sigma_p[1][1]);#std. deviation of n
		self.sigma_epsilon = np.sqrt(self.Sigma_p[2][2]);#std. deviation of epsilon
		self.r_alpha_n = self.Sigma_p[0][1]/(self.sigma_alpha*self.sigma_n);
		self.r_alpha_epsilon = self.Sigma_p[0][1]/(self.sigma_alpha*self.sigma_epsilon);
		self.r_n_epsilon = self.Sigma_p[0][1]/(self.sigma_n*self.sigma_epsilon);
		self.cholskyDeCorrelateMat = L1
		
		
		#Solving for optimized zeta values
		z = np.array([1,1,1]);#guess value to get optimized zeta values
		con1 = {'type': 'ineq', 'fun': self.normRandZetaCons1}
		con2 = {'type': 'ineq', 'fun': self.normRandZetaCons2}
		self.cons = [con1, con2]
		self.zeta = minimize(self.normRandZetaObjective,z,constraints=self.cons)
				
		
		# Plots for Temperature dependent uncertainity versus their respective Temperatures 
		self.T = self.temperatures;
		#self.theta = np.array([self.T/self.T,np.log(self.T),-1/self.T])
		self.f = self.getUncertainty(self.T)		
		#self.zetaUnsrt = np.asarray(np.dot(self.theta.T,np.array(np.dot(self.cholskyDeCorrelateMat,self.zeta.x)).flatten())).flatten()
		#print(self.f,self.zetaUnsrt,self.uncertainties)
		'''
		fig = plt.figure()
		plt.xlabel('Temperature (K)')
		plt.ylabel('Uncertainity ($\sigma$)')
		plt.title('Reaction {}'.format(self.name), fontsize = 10)		
		plt.plot(self.temperatures,self.uncertainties,'o',label='exp uncertainties');
		plt.plot(self.T,self.f,'-',label='Cholesky decomposition')
		plt.ylim(0,2*max(self.uncertainties[0],self.uncertainties[-1]))	
		plt.legend()
		my_path = os.getcwd()
		plt.savefig(my_path+'/Plots/TDeptUnsrt/'+self.name+'.png')
		
		# Plots zeta values
		
		fig = plt.figure()
		plt.xlabel('Temperature (K)')
		plt.ylabel('zeta function ($\sigma$)')
		plt.title('zeta_values for {}\n: {}'.format(self.name,self.zeta.x), fontsize = 10)	
		plt.plot(self.temperatures,self.uncertainties,'o',label='exp uncertainties');
		plt.plot(self.T,self.fit_zeta(self.T,self.zeta.x),'-',label='zeta_function')
		plt.ylim(0,2*max(self.uncertainties[0],self.uncertainties[-1]))	
		plt.legend()
		my_path = os.getcwd()
		plt.savefig(my_path+'/Plots/Zeta/'+self.name+'.png')
		plt.close()
		
		'''
		
#public function to get the uncertainity values for discrete temperatures
	def getDtList(self):
		self.Unsrt_dict = {}
		key = ["Tag","Solution","Class","Type","Sub_type","Branch_boolen","Branches","Pressure_limit","Common_temp","temp_list","Nominal","Exp_input_data_type","priorCovariance","Basis_vector","Uncertainties","Temperatures","unsrt_func","Data_key"]
		values = [self.tag,self.solution,self.classification,self.type,self.sub_type,self.branching,self.branches,self.pressure_limit,self.common_temp,self.temp_limit,self.rxn_dict,self.exp_data_type,self.cholskyDeCorrelateMat,self.zeta.x,self.uncertainties,self.temperatures,self.f,""]
		for i,element in enumerate(key):
			self.Unsrt_dict[element] = values[i]
		return self.Unsrt_dict
	def getAllData(self):
		if len(self.branches.split(","))==1:
			b1 = self.branches.split(",")[0]
			b2 = ""
			b3 = ""
		elif len(self.branches.split(",")) >1:
			b1 = self.branches.split(",")[0]
			b2 = self.branches.split(",")[1]
			if len(self.branches.split(",")) >2:
				b3 = self.branches.split(",")[3]
			else:
				b3 = ""
		else:
			b1 = ""
			b2 = ""
			b3 = ""
		exp_data_type = self.exp_data_type.split(";")[0]
		exp_format = self.exp_data_type.split(";")[1]
		Log_string = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.classification,self.type,self.sub_type,self.branching,b1,b2,b3,self.pressure_limit,self.common_temp,self.temp_limit,exp_data_type,exp_format,self.nametag)
		exp_unsrt_string = ""
		solver_log = "######################\n{}######################\n\t\t{}\n\n".format(self.nametag,self.solution)
		calc_cholesky = self.cholskyDeCorrelateMat
		zeta_string = "{}".format(self.zeta)
		for i in range(len(self.temperatures)):
			exp_unsrt_string += "{}\t{}\n".format(self.temperatures[i],self.uncertainties[i])
		return Log_string,exp_unsrt_string,solver_log,calc_cholesky,zeta_string
	def getUncertainty(self,T): 
		M = 3/(np.log(10.0));# Truncation at 3 sigma
		
		return M*np.sqrt((self.solution.x[0]**2+2*self.solution.x[0]*self.solution.x[1]*np.log(T)+(self.solution.x[1]*np.log(T))**2-2*self.solution.x[0]*self.solution.x[3]/T-2*self.solution.x[1]*self.solution.x[3]*np.log(T)/T+(self.solution.x[3]/T)**2+(self.solution.x[2]*np.log(T))**2-2*self.solution.x[4]*self.solution.x[2]*np.log(T)/T+(self.solution.x[4]/T)**2+(self.solution.x[5]/T)**2))
 
#public function to get the temperature range for the uncertainity of a perticular reaction
	def getTempRange(self):
		return self.temperatures[0], self.temperatures[-1]
	def getTemperature(self):
		return self.temperatures
	def getRxnType(self):
		return self.type,self.branching,self.branches

#public function to get the zeta values for a perticulat reaction
	def getData(self):
		
		return self.zeta.x
	
	def zetaValues(self):
		return self.zeta.x
#public function to get the cholesky decomposed matrix for normalization of variables
	
	
	def getCholeskyMat(self):
		return self.cholskyDeCorrelateMat

	def fit_zeta(self,T,z):
		func = (self.L11*z[0]+(self.L21*z[0]+self.L22*z[1])*np.log(T)-(self.L31*z[0]+self.L32*z[1]+self.L33*z[2])*(1/T))
		return func
	def uncertPriorObjective (self,uncertPriorGuess):
		M = 3/(np.log(10.0));# Truncation at 3 sigma
		y = np.array(self.uncertainties)#experimental uncertainity
		T = np.array(self.temperatures)#Temperature domain
		f = np.zeros(len(T))
		f = ((y) - M*np.sqrt((uncertPriorGuess[0]**2+2*uncertPriorGuess[0]*uncertPriorGuess[1]*np.log(T)+(uncertPriorGuess[1]*np.log(T))**2-2*uncertPriorGuess[0]*uncertPriorGuess[3]/T-2*uncertPriorGuess[1]*uncertPriorGuess[3]*np.log(T)/T+(uncertPriorGuess[3]/T)**2+(uncertPriorGuess[2]*np.log(T))**2-2*uncertPriorGuess[4]*uncertPriorGuess[2]*np.log(T)/T+(uncertPriorGuess[4]/T)**2+(uncertPriorGuess[5]/T)**2)))/(y/M);
		obj = np.dot(f,f);
		return obj


	def uncertPriorCons(self,uncertPriorGuess):
		M = 3/(np.log(10.0));
		y = np.array(self.uncertainties)#experimental uncertainity
		T = np.array(self.temperatures)#Temperature domain
		f = np.zeros(len(T));
		f = (y) - M*np.sqrt((uncertPriorGuess[0]**2+2*uncertPriorGuess[0]*uncertPriorGuess[1]*np.log(T)+(uncertPriorGuess[1]*np.log(T))**2-2*uncertPriorGuess[0]*uncertPriorGuess[3]/T-2*uncertPriorGuess[1]*uncertPriorGuess[3]*np.log(T)/T+(uncertPriorGuess[3]/T)**2+(uncertPriorGuess[2]*np.log(T))**2-2*uncertPriorGuess[4]*uncertPriorGuess[2]*np.log(T)/T+(uncertPriorGuess[4]/T)**2+(uncertPriorGuess[5]/T)**2))
		return np.amin(f)


	def normRandZetaObjective(self,z):
		T =np.linspace(self.temperatures[0],self.temperatures[-1],1000)
		dd = self.getUncertainty(T)
		fdiff = (dd-(self.L11*z[0]+(self.L21*z[0]+self.L22*z[1])*np.log(T)-(self.L31*z[0]+self.L32*z[1]+self.L33*z[2])*(1/T)))
		obj = np.dot(fdiff,fdiff)
		return obj


	def normRandZetaCons1(self,z):
		T =np.linspace(self.temperatures[0],self.temperatures[-1],1000)
		dd = self.getUncertainty(T)
		f = (dd-(self.L11*z[0]+(self.L21*z[0]+self.L22*z[1])*np.log(T)-(self.L31*z[0]+self.L32*z[1]+self.L33*z[2])*(1/T)))
		return np.amin(f)


	def normRandZetaCons2(self,z):
		T =np.linspace(self.temperatures[0],self.temperatures[-1],1000)
		f = (self.L11*z[0]+(self.L21*z[0]+self.L22*z[1])*np.log(T)-(self.L31*z[0]+self.L32*z[1]+self.L33*z[2])*(1/T))
		return np.amin(f)

#############################################################
###       center broadening factors for showing        ######
###          decomposition and recombination of        ###### 
###          pressure dependent reactions              ######
#############################################################

class fallOffCurve:
	def __init__(self, Element,mechPath):
		self.rxn = self.classification = self.type = self.sub_type = self.exp_data_type = self.temperatures = self.uncertainties = None
		
		self.rxn = Element.attrib["rxn"]
		self.tag = Element.tag	
		self.foc_dict = IFR.MechParsing(mechPath).getFocData(self.rxn)[0]
		#print(self.rxn_dict)
		for item in Element:
			if item.tag == "class":
				self.classification = item.text
			if item.tag == "type":
				self.type = item.text
			if item.tag == "sub_type":
				self.sub_type = item.attrib["name"]
				for subitem in item:
					if subitem.tag == "multiple":
						self.branching = subitem.text
					if subitem.tag == "branches":
						self.branches = subitem.text
					if subitem.tag == "pressure_limit":
						self.pressure_limit = subitem.text
					if subitem.tag == "common_temp":
						self.common_temp = subitem.text
					if subitem.tag == "temp_limit":
						self.temp_limit = subitem.text
			
			if item.tag == "data_type":
				self.exp_data_type = item.text
			if item.tag == "temp":
				#print(item.text)
				self.temperatures = np.asarray([float(i) for i in item.text.split(",")])
			if item.tag == "unsrt":
				self.uncertainties = np.asarray([float(i) for i in item.text.split(",")])
			
		if self.exp_data_type.split(";")[0] == "constant":
			if self.exp_data_type.split(";")[1] == "array":	
				self.temperatures = self.temperatures
				self.uncertainties = self.uncertainties
				
			elif self.exp_data_type.split(";")[1] == "end_points":
				self.temperatures = np.linspace(self.temperatures[0],self.temperatures[1],20)
				self.uncertainties = np.linspace(self.uncertainties[0],self.uncertainties[1],20)
		
		
#		for item in Element:
#			if item.tag == "class":
#				self.classification = item.text
#			if item.tag == "type":
#				self.type = item.text
#			if item.tag == "sub_type":
#				self.sub_type = item.attrib["name"]
#				for subitem in item:
#					if subitem.tag == "branch":
#						self.branching = subitem.text
#					if subitem.tag == "branches":
#						self.branches = subitem.text
#					if subitem.tag == "pressure_limit":
#						self.pressure_limit = subitem.text
#					if subitem.tag == "common_temp":
#						self.common_temp = subitem.text
#					if subitem.tag == "temp_limit":
#						self.temp_limit = subitem.text
#			
#			if item.tag == "data_type":
#				self.exp_data_type = item.text
#			if item.tag == "temp":
#				self.temperatures = np.asarray([float(i) for i in item.text.split(",")])
#			if item.tag == "unsrt":
#				self.uncertainties = np.asarray([float(i) for i in item.text.split(",")])
#			
#		if self.exp_data_type.split(";")[0] == "constant":
#			if self.exp_data_type.split(";")[1] == "array":
#				self.temperatures = self.temperatures
#				self.uncertainties = self.uncertainties
#			elif self.exp_data_type.split(";")[1] == "end_points":
#				self.temperatures = np.linspace(self.temperatures[0],self.temperatures[1],20)
#				self.uncertainties = np.linspace(self.uncertainties[0],self.uncertainties[1],20)
		
		self.nametag = self.rxn+":"+self.sub_type

		L11 = self.uncertainties[0]/3
		L = np.array([L11]) #Cholesky_lower_triangular_matrix
		self.cholskyDeCorrelateMat = L
		#print(L)
		self.zeta = 1
		
		self.T = self.temperatures
		self.f = self.getUncertainty(self.T)		
		self.solution = self.cholskyDeCorrelateMat
#public function to get the uncertainity values for discrete temperatures
	def getAllData(self):
		Log_string = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.classification,self.type,self.sub_type,self.branching,self.branches,self.pressure_limit,self.common_temp,self.temp_limit,self.exp_data_type,self.nametag)
		exp_unsrt_string = ""
		solver_log = "{}\n{}\n".format(self.nametag,self.solution)
		calc_cholesky = self.cholskyDeCorrelateMat
		zeta_string = "{}".format(self.zeta)
		for i in range(len(self.temperatures)):
			exp_unsrt_string += "{}\t{}\n".format(self.temperatures[i],self.uncertainties[i])
		return Log_string,exp_unsrt_string,solver_log,calc_cholesky,zeta_string

	def getDtList(self):
		self.Unsrt_dict = {}
		key = ["Tag","Solution","Class","Type","Sub_type","Branch_boolen","Branches","Pressure_limit","Common_temp","temp_list","Nominal","Exp_input_data_type","priorCovariance","Basis_vector","Uncertainties","Temperatures","unsrt_func","Data_key"]
		values = [self.tag,self.solution,self.classification,self.type,self.sub_type,self.branching,self.branches,self.pressure_limit,self.common_temp,self.temp_limit,self.foc_dict,self.exp_data_type,self.cholskyDeCorrelateMat,self.zeta,self.uncertainties,self.temperatures,self.f,""]
		for i,element in enumerate(key):
			self.Unsrt_dict[element] = values[i]
		return self.Unsrt_dict
	
	def getUncertainty(self,T):
		L11 = self.uncertainties[0]/3
		
		Foc_unsrt = 3*np.sqrt((L11*(T/T))**2)
		return Foc_unsrt
 
#public function to get the temperature range for the uncertainity of a perticular reaction
	def getTempRange(self):
		return self.temperatures[0], self.temperatures[-1]
	def getTemperature(self):
		return self.temperatures
	def getRxnType(self):
		
		return self.type,self.IsBranching,self.branchRxn
	def fit_zeta(self,T,z):
		L11 = self.uncertainties[0]/3
		func =z*L11*(T/T)
		
		return func

	def getZeta(self):
		#printZeta = "{}\t{}\t{}\n".format(self.zeta.x[0],self.zeta.x[1],self.zeta.x[2])
		#printL ="{}".format(self.cholskyDeCorrelateMat)
		#printL+="\n"
		#fileZeta = open('../zetaValues.txt','a')
		#fileL = open('../Cholesky.txt','a')
		#fileZeta.write(printZeta)
		#fileL.write(printL)
		#fileZeta.close()
		#fileL.close()
		#print("\n{}\n".format(self.zeta.x));
		return self.zeta
	
	def zetaValues(self):
		return self.zeta
#public function to get the cholesky decomposed matrix for normalization of variables
	
	def getCholeskyMat(self):
		return self.cholskyDeCorrelateMat


#############################################################
###       Uncertainty for heat capacities              ######
###       of kinetic species                           ######
#############################################################


class thermodynamic:
	def __init__(self, Element,thermo_loc):
		self.species = self.classification = self.type = self.sub_type =  self.branching = self.branches  = self.pressure_limit  = self. common_temp = self.temp_limit= None
		self.tag = Element.tag
		IFRT = IFR.ThermoParsing(thermo_loc)
				
		self.exp_data_type = {}
		self.temperatures = {}
		self.uncertainties = {}
		self.cholskyDeCorrelateMat = {}
		self.zeta = {}
		self.species = Element.attrib["species"]
		self.nominal = {}
		
		
		for item in Element:
			#print(item.tag)
			if item.tag == "class":
				self.classification = item.text
				continue
			if item.tag == "type":
				self.type = item.text
				continue
			if item.tag == "sub_type":
				self.sub_type = item.attrib["name"]
				
				for subitem in item:
					if "multiple" in subitem.tag:
						self.branching = str(subitem.text)
						continue
					if "branches" in subitem.tag:
						self.branches = str(subitem.text)
						continue
						
					if "pressure_limit" in subitem.tag:
						self.pressure_limit = str(subitem.text)
						continue
					if "common_temp" in subitem.tag:
						self.common_temp = str(subitem.text)
						continue
					if  "temp_limit" in subitem.tag :
						self.temp_limit = str(subitem.text)
						#print(self.temp_limit)
						if self.temp_limit == "Low":
		
							self.thermo_dict = IFR.ThermoParsing(thermo_loc).getThermoLow(self.species)
						else:
							self.thermo_dict = IFR.ThermoParsing(thermo_loc).getThermoHigh(self.species)
						continue
				continue
				
		for i in self.branches.split(","):
			
			for item in Element:
				if item.tag == "data_type_"+str(i):
					self.exp_data_type[str(i)] = item.text
					continue
					
				if item.tag == "temp_"+str(i):
					#print(item.text)
					self.temperatures[str(i)] = np.asarray([float(j) for j in item.text.split(",")])
					continue
						
				if item.tag == "unsrt_"+str(i):
					self.uncertainties[str(i)] = np.asarray([float(i) for i in item.text.split(",")])
					continue
			
		for i in self.exp_data_type:
			if self.exp_data_type[str(i)].split(";")[0] ==  "constant":
				if self.exp_data_type[str(i)].split(";")[1] == "array":
					continue
				elif self.exp_data_type[str(i)].split(";")[1] == "end_points":
					self.temperatures[str(i)] = np.linspace(self.temperatures[str(i)][0],self.temperatures[str(i)][1],20)
					self.uncertainties[str(i)] = np.linspace(self.uncertainties[str(i)][0],self.uncertainties[str(i)][1],20)
					continue

			elif self.exp_data_type[str(i)].split(";")[0] == "percentage":
				if self.exp_data_type[str(i)].split(";")[1] == "array":
					a = self.thermo_dict[str(i)]
					func = IFRT.function(str(i),a,self.temperatures[str(i)])
					y = self.uncertainties[str(i)]
					#print(y)						
					self.uncertainties[str(i)] = np.asarray(np.dot(y,func)).flatten()
					continue

				elif self.exp_data_type[str(i)].split(";")[1] == "end_points":

					self.temperatures[str(i)] = np.linspace(self.temperatures[str(i)][0],self.temperatures[str(i)][1],20)
					self.uncertainties[str(i)] = np.linspace(self.uncertainties[str(i)][0],self.uncertainties[str(i)][1],20)

					a = self.thermo_dict[str(i)]
					func = IFRT.function(str(i),a,self.temperatures[str(i)])
					#print(str(i),a,self.temperatures[str(i)])
					#print(func)
					y = self.uncertainties[str(i)]						
					#print(y*func)
					self.uncertainties[str(i)] = np.asarray((y*func)).flatten()
					continue
			
		
			
		self.doUnsrtAnalysis()
		for i in self.branches.split(","):
			self.cholskyDeCorrelateMat[str(i)],self.zeta[str(i)] = self.unsrt(str(i))
		self.nametag = self.species+":"+self.temp_limit	
			
		self.corelate_block = block_diag(*(self.cholskyDeCorrelateMat["Hcp"],self.cholskyDeCorrelateMat["h"],self.cholskyDeCorrelateMat["e"]))
		
		self.f = self.getUncertainty(self.temperatures)
		
		
	def getAllData(self):
		b1 = self.branches.split(",")[0]
		b2 = self.branches.split(",")[1]
		if len(self.branches.split(",")) >2:
			b3 = self.branches.split(",")[3]
		else:
			b3 = ""
		exp_data_type = self.exp_data_type["Hcp"].split(";")[0]
		exp_format = self.exp_data_type["Hcp"].split(";")[1]
		Log_string = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.classification,self.type,self.sub_type,self.branching,b1,b2,b3,self.pressure_limit,self.common_temp,self.temp_limit,exp_data_type,exp_format,self.nametag)
		exp_unsrt_string = ""
		solver_log = "######################\n{}######################\n\t\t{}\n\n".format(self.nametag,self.solution)
		calc_cholesky = self.corelate_block
		zeta_string = "{}\t{}\t{}\n".format(self.zeta["Hcp"],self.zeta["h"],self.zeta["e"])
		#print(self.temperatures)
		#print(self.uncertainties)
		for i in range(len(self.temperatures["Hcp"])):
			exp_unsrt_string += "{}\t{}\t{}\t{}\t{}\t{}\n".format(self.temperatures["Hcp"][i],self.uncertainties["Hcp"][i],self.temperatures["h"][i],self.uncertainties["h"][i],self.temperatures["e"][i],self.uncertainties["e"][i])
		string_2 = "temp\tHcp\ttemp\th\ttemp\te\t\n"
		file_unsrt = open("./Data/"+self.nametag+"_usrtData.log","w")
		file_unsrt.write(string_2+"\n"+exp_unsrt_string)
		file_unsrt.close()
		return Log_string,exp_unsrt_string,solver_log,calc_cholesky,zeta_string
	
	
	def getDtList(self):
		self.Unsrt_dict = {}
		key = ["Tag","Solution","Class","Type","Sub_type","Branch_boolen","Branches","Pressure_limit","Common_temp","temp_list","Nominal","Exp_input_data_type","priorCovariance","Basis_vector","Uncertainties","Temperatures","unsrt_func","Data_key"]
		values = [self.tag,self.solution,self.classification,self.type,self.sub_type,self.branching,self.branches,self.pressure_limit,self.common_temp,self.temp_limit,self.thermo_dict,self.exp_data_type,self.cholskyDeCorrelateMat,self.zeta,self.uncertainties,self.temperatures,self.f,"Hcp,h,e"]		
		
		for i,element in enumerate(key):
			self.Unsrt_dict[element] = values[i]
		return self.Unsrt_dict
		
	def doUnsrtAnalysis(self):
		T = self.temperatures[self.branches.split(",")[0]]
		guess = 0.001*np.ones(17)
		self.solution = minimize(self.uncertPriorObjective,guess,method="Nelder-Mead",options={'maxiter': 100000, 'maxfev': 100000, 'disp': False, 'return_all': False, 'initial_simplex': None, 'xatol': 1E-05, 'fatol': 1E-05, 'adaptive': True})
		Tscale = 5000
		#print(self.solution)
		L11 = self.solution.x[0]
		L21 = self.solution.x[1]/Tscale
		L22 = self.solution.x[2]/Tscale
		L31 = self.solution.x[3]/Tscale**2
		L32 = self.solution.x[4]/Tscale**2
		L33 = self.solution.x[5]/Tscale**2
		L41 = self.solution.x[6]/Tscale**3
		L42 = self.solution.x[7]/Tscale**3
		L43 = self.solution.x[8]/Tscale**3
		L44 = self.solution.x[9]/Tscale**3
		L51 = self.solution.x[10]/Tscale**4
		L52 = self.solution.x[11]/Tscale**4
		L53 = self.solution.x[12]/Tscale**4
		L54 = self.solution.x[13]/Tscale**4
		L55 = self.solution.x[14]/Tscale**4
		L66 = self.solution.x[15]
		L77 = self.solution.x[16]
		self.Lcp = np.array([[L11,L21,L31,L41,L51],[0,L22,L32,L42,L52],[0,0,L33,L43,L53],[0,0,0,L44,L54],[0,0,0,0,L55]])
		self.LH = np.array([L66])
		self.LS = np.array([L77])
		
		self.cholskyDeCorrelateMat_cp = np.matrix(self.Lcp.T)
		self.cholskyDeCorrelateMat_H  = np.matrix(self.LH.T)
		self.cholskyDeCorrelateMat_S  = np.matrix(self.LS.T)
		
		
		theta_cp = np.array([T/T,T,T**2,T**3,T**4])
		theta_H = np.array([T/T])
		theta_S = np.array([T/T])
		
		#Find zeta values
		guess_zeta = 0.01*np.array([1,1,1,1,1,1,1])
		self.zeta = minimize(self.obj_zeta,guess_zeta)
		self.zeta_cp = self.zeta.x[0:5]
		self.zeta_h = self.zeta.x[5]
		self.zeta_s = self.zeta.x[6]
	
	def unsrt(self,index):
		if index == "Hcp":
			L = self.Lcp
			z = self.zeta_cp
		if index == "h":
			L = self.LH
			z = self.zeta_h
		if index == "e":
			L = self.LS
			z = self.zeta_s
		return L,z
		
		
	def func_4(self,T,L11,L12,L22,L13,L23,L33,L14,L24,L34,L44,L15,L25,L35,L45,L55):
		unsrt = 9*((L55*T**4)**2+(L44*T**3+L55*T**4)**2+(L33*T**2+L34*T**3+L35*T**4)**2+(L22*T+L23*T**2+L24*T**3+L25*T**4)**2+(L11+L12*T+L13*T**2+L14*T**3+L15*T**4)**2)
		return unsrt
	def func_zeta(self,T,L0,L1,L2,L3,L4):
		z = np.array([L0,L1,L2,L3,L4])
		L11 = self.solution[0]
		L12 = self.solution[1]
		L22 = self.solution[2]
		L13 = self.solution[3]
		L23 = self.solution[4]
		L33 = self.solution[5]
		L14 = self.solution[6]
		L24 = self.solution[7]
		L34 = self.solution[8]
		L44 = self.solution[9]
		L15 = self.solution[10]
		L25 = self.solution[11]
		L35 = self.solution[12]
		L45 = self.solution[13]
		L55 = self.solution[14]
		fdiff = ((z[0]*L11)+(z[0]*L12+z[1]*L22)*T+(z[0]*L13+z[1]*L23+z[2]*L33)*T**2+(z[0]*L14+z[1]*L24+z[2]*L34+z[3]*L44)*(T**3)+(L15*z[0]+L25*z[1]+L35*z[2]+L45*z[3]+L55*z[4])*T**4)    
		return fdiff

	def getUncertainty(self,T): 
		unsrt = {}
		Tscale = 5000
		T = T["Hcp"]
		L11 = self.solution.x[0]
		L21 = self.solution.x[1]/Tscale
		L22 = self.solution.x[2]/Tscale
		L31 = self.solution.x[3]/Tscale**2
		L32 = self.solution.x[4]/Tscale**2
		L33 = self.solution.x[5]/Tscale**2
		L41 = self.solution.x[6]/Tscale**3
		L42 = self.solution.x[7]/Tscale**3
		L43 = self.solution.x[8]/Tscale**3
		L44 = self.solution.x[9]/Tscale**3
		L51 = self.solution.x[10]/Tscale**4
		L52 = self.solution.x[11]/Tscale**4
		L53 = self.solution.x[12]/Tscale**4
		L54 = self.solution.x[13]/Tscale**4
		L55 = self.solution.x[14]/Tscale**4
		L66 = self.solution.x[15]
		L77 = self.solution.x[16]
		# Truncation at 3 sigma
		unsrt_cp = 3*np.sqrt((L55*T**4)**2+(L44*T**3+L55*T**4)**2+(L33*T**2+L43*T**3+L53*T**4)**2+(L22*T+L32*T**2+L42*T**3+L52*T**4)**2+(L11+L21*T+L31*T**2+L41*T**3+L51*T**4)**2)
		unsrt_H = 3*np.sqrt((L66*(T/T))**2)
		unsrt_S = 3*np.sqrt((L77*(T/T))**2)
		unsrt["Hcp"] = unsrt_cp
		unsrt["h"] = unsrt_H
		unsrt["e"] = unsrt_S
		return unsrt
	
	def uncertPriorObjective (self,guess):
		Tscale = 5000
		R = 8.314
		Z= guess
		#z = np.array([Z[15],Z[16],Z[17],Z[18],Z[19]])
		L11 =Z[0]
		L21 =Z[1]
		L22 =Z[2]
		L31 =Z[3]
		L32 =Z[4]
		L33 =Z[5]
		L41 =Z[6]
		L42 =Z[7]
		L43 =Z[8]
		L44 =Z[9]
		L51 =Z[10]
		L52 =Z[11]
		L53 =Z[12]
		L54 =Z[13]
		L55 =Z[14]
		L66=Z[15]
		L77=Z[16]
		Lcp = np.array([[L11,L21,L31,L41,L51],[0,L22,L32,L42,L52],[0,0,L33,L43,L53],[0,0,0,L44,L54],[0,0,0,0,L55]])
		Lh = np.array([L66])
		Ls = np.array([L77])

		if "h" in self.sub_type:
			Y_h = self.uncertainties["h"]
			thetaH = np.array([T/T])
			sigma_H = 9*(np.dot(Lh,thetaH))**2

		if "e" in self.sub_type:
			Y_s = self.uncertainties["e"]
			thetaS = np.array([T/T])
			sigma_S = 9*(np.dot(Ls,thetaS))**2

		if "Hcp" in self.sub_type:
			Y_cp = self.uncertainties["Hcp"]
			thetaCP = np.array([T/T,T,T**2,T**3,T**4])
			unsrt = 9*((L55*T**4)**2+(L44*T**3+L55*T**4)**2+(L33*T**2+L43*T**3+L53*T**4)**2+(L22*T+L32*T**2+L42*T**3+L52*T**4)**2+(L11+L21*T+L31*T**2+L41*T**3+L51*T**4)**2)

		T = self.temperatures[self.branches.split(",")[0]]/Tscale


		unsrt = 9*((L55*T**4)**2+(L44*T**3+L55*T**4)**2+(L33*T**2+L43*T**3+L53*T**4)**2+(L22*T+L32*T**2+L42*T**3+L52*T**4)**2+(L11+L21*T+L31*T**2+L41*T**3+L51*T**4)**2)


		if "Hcp" in self.sub_type:
			residual_cp = (Y_cp-np.sqrt(unsrt))/(Y_cp/3)
		else:
			residual_cp = np.array([0])
		if "h" in self.sub_type:
			residual_h =(Y_h-np.sqrt(sigma_H))/(Y_h/3)
		else:
			residual_h =  np.array([0])
		if "e" in self.sub_type:
			residual_s =(Y_s-np.sqrt(sigma_S))/(Y_s/3)
		else:
			residual_s = np.array([0])
	
		obj = np.dot(residual_cp.T,residual_cp)+np.dot(residual_h.T,residual_h)+np.dot(residual_s.T,residual_s)
		return obj
	
	def obj_zeta(self,guess):
		T = self.temperatures[self.branches.split(",")[0]]
		#T = np.linspace(1000,5000,100)
		Tscale =5000
		z = np.ones(7)
		L11 = self.solution.x[0]
		L21 = self.solution.x[1]/Tscale
		L22 = self.solution.x[2]/Tscale
		L31 = self.solution.x[3]/Tscale**2
		L32 = self.solution.x[4]/Tscale**2
		L33 = self.solution.x[5]/Tscale**2
		L41 = self.solution.x[6]/Tscale**3
		L42 = self.solution.x[7]/Tscale**3
		L43 = self.solution.x[8]/Tscale**3
		L44 = self.solution.x[9]/Tscale**3
		L51 = self.solution.x[10]/Tscale**4
		L52 = self.solution.x[11]/Tscale**4
		L53 = self.solution.x[12]/Tscale**4
		L54 = self.solution.x[13]/Tscale**4
		L55 = self.solution.x[14]/Tscale**4
		L66 = self.solution.x[15]
		L77 = self.solution.x[16]
		z[0] = guess[0]
		z[1] = guess[1]
		z[2] = guess[2]
		z[3] = guess[3]
		z[4] = guess[4]
		z[5] = guess[5]
		z[6] = guess[6]
		
		zetaFunc_cp = ((z[0]*L11)*(T/T)+(z[0]*L21+z[1]*L22)*T+(z[0]*L31+z[1]*L32+z[2]*L33)*T**2+(z[0]*L41+z[1]*L42+z[2]*L43+z[3]*L44)*(T**3)+(L51*z[0]+L52*z[1]+L53*z[2]+L54*z[3]+L55*z[4])*T**4)
		zetaFunc_H = (T/T)*z[5]*L66
		zetaFunc_S = (T/T)*z[6]*L77
		
		
		if "Hcp" in self.sub_type:
			residual_cp =self.uncertainties["Hcp"]-zetaFunc_cp
		else:
			residual_cp = 0
		if "h" in self.sub_type:
			residual_H = self.H_uncertainties["h"]-zetaFunc_H
		else:
			residual_H = 0
		if "e" in self.sub_type:
			residual_S = self.S_uncertainties["e"]-zetaFunc_S
		else:
			residual_S = 0	
		
		obj = np.dot(residual_cp,residual_cp)+np.dot(residual_H,residual_H)+np.dot(residual_S,residual_S)
		return obj
		
#############################################################
###       Uncertainty analysis for fallOffCurves,      ######
###       enthalpy, entropy and transport              ######
###       properties for kinetic species               ######
#############################################################		
#Temperature independent uncetainties
#zeta = 3
#sigma_p = L11
	
class transport:
	def __init__(self, Element,transport_loc):
		self.species = self.classification = self.type = self.sub_type  = self. branching = self.branches  = self.pressure_limit  = self. common_temp = self.temp_limit =  None
		self.tag = Element.tag
		self.exp_data_type = {}
		self.temperatures = {}
		self.uncertainties = {}
		self.cholskyDeCorrelateMat = {}
		self.zeta = {} 
		self.species = Element.attrib["species"]
		
		self.trans_dict = IFR.TransportParsing(transport_loc).getTransportData(self.species)
		
		for item in Element:
			#print(item.tag)
			if item.tag == "class":
				self.classification = item.text
				continue
			if item.tag == "type":
				self.type = item.text
				continue
			if item.tag == "sub_type":
				self.sub_type = item.attrib["name"]
				
				for subitem in item:
					if "multiple" in subitem.tag:
						self.branching = str(subitem.text)
						continue
					if "branches" in subitem.tag:
						self.branches = str(subitem.text)
						continue
						
					if "pressure_limit" in subitem.tag:
						self.pressure_limit = str(subitem.text)
						continue
					if "common_temp" in subitem.tag:
						self.common_temp = str(subitem.text)
						continue
					if  "temp_limit" in subitem.tag :
						self.temp_limit = str(subitem.text)
						
						continue
				continue
					
		for i in self.branches.split(","):
			
			for item in Element:
				if item.tag == "data_type_"+str(i):

					self.exp_data_type[str(i)] = item.text

					continue
					
				if item.tag == "temp_"+str(i):
					#print(item.text)
					self.temperatures[str(i)] = np.asarray([float(j) for j in item.text.split(",")])
					continue
						
				if item.tag == "unsrt_"+str(i):
					self.uncertainties[str(i)] = np.asarray([float(i) for i in item.text.split(",")])
					continue
		
		
		for i in self.exp_data_type:
			
			if self.exp_data_type[str(i)].split(";")[0] ==  "constant":
				
				if self.exp_data_type[str(i)].split(";")[1] == "array":
				
					continue
				
				elif self.exp_data_type[str(i)].split(";")[1] == "end_points":
				
					self.temperatures[str(i)] = np.linspace(self.temperatures[str(i)][0],self.temperatures[str(i)][1],20)
				
					self.uncertainties[str(i)] = np.linspace(self.uncertainties[str(i)][0],self.uncertainties[str(i)][1],20)
				
					continue

			elif self.exp_data_type[str(i)].split(";")[0] == "percentage":
				
				if self.exp_data_type[str(i)].split(";")[1] == "array":
				
					a = self.trans_dict[str(i)]
					y = self.uncertainties[str(i)]						
					self.uncertainties[str(i)] = y*a
					
					continue

				elif self.exp_data_type[str(i)].split(";")[1] == "end_points":

					self.temperatures[str(i)] = np.linspace(self.temperatures[str(i)][0],self.temperatures[str(i)][1],20)
					self.uncertainties[str(i)] = np.linspace(self.uncertainties[str(i)][0],self.uncertainties[str(i)][1],20)

					a = self.trans_dict[str(i)]
					y = self.uncertainties[str(i)]						
					self.uncertainties[str(i)] = y*a
					continue
				
			
			
		
		for i in self.branches.split(","):
			self.cholskyDeCorrelateMat[str(i)],self.zeta[str(i)] = self.unsrt(str(i))
			
			
		self.nametag = self.species	
		self.f = self.getUncertainty(self.temperatures)
		'''
		fig = plt.figure()
		plt.xlabel('Temperature (K)')
		plt.ylabel('Uncertainity ($\sigma$)')
		plt.title( string+'{}'.format(self.name), fontsize = 10)		
		plt.plot(self.temperatures,self.uncertainties,'o',label='exp uncertainties');
		plt.ylim(0,2*max(self.uncertainties[0],self.uncertainties[-1]))	
		plt.legend()
		my_path = os.getcwd()
		plt.savefig(my_path+'/Plots/TDeptUnsrt/'+self.name+'.png')
		
		'''	
		self.solution = self.cholskyDeCorrelateMat
		self.corelate_block = block_diag(*(self.cholskyDeCorrelateMat["LJe"],self.cholskyDeCorrelateMat["LJs"]))
	
	
	
	def getDtList(self):
		self.Unsrt_dict = {}
		key = ["Tag","Solution","Class","Type","Sub_type","Branch_boolen","Branches","Pressure_limit","Common_temp","temp_list","Nominal","Exp_input_data_type","priorCovariance","Basis_vector","Uncertainties","Temperatures","unsrt_func","Data_key"]
		values = [self.tag,self.solution,self.classification,self.type,self.sub_type,self.branching,self.branches,self.pressure_limit,self.common_temp,self.temp_limit,self.trans_dict,self.exp_data_type,self.cholskyDeCorrelateMat,self.zeta,self.uncertainties,self.temperatures,self.f,"LJe,LJs"]
		
		for i,element in enumerate(key):
			self.Unsrt_dict[element] = values[i]
		return self.Unsrt_dict
		
	def getAllData(self):
		b1 = self.branches.split(",")[0]
		b2 = self.branches.split(",")[1]
		if len(self.branches.split(",")) >2:
			b3 = self.branches.split(",")[3]
		else:
			b3 = ""
		exp_data_type = self.exp_data_type["LJe"].split(";")[0]
		exp_format = self.exp_data_type["LJe"].split(";")[1]
		Log_string = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.classification,self.type,self.sub_type,self.branching,b1,b2,b3,self.pressure_limit,self.common_temp,self.temp_limit,exp_data_type,exp_format,self.nametag)
		exp_unsrt_string = ""
		solver_log = "######################\n{}######################\n\t\t{}\n\n".format(self.nametag,self.solution)
		calc_cholesky = self.corelate_block
		zeta_string = "{}\t{}\n".format(self.zeta["LJe"],self.zeta["LJs"])
		#print(self.temperatures)
		#print(self.uncertainties)
		for i in range(len(self.temperatures["LJe"])):
			exp_unsrt_string += "{}\t{}\t{}\t{}\n".format(self.temperatures["LJe"][i],self.uncertainties["LJe"][i],self.temperatures["LJs"][i],self.uncertainties["LJs"][i])
		string_2 = "temp\tLJe\ttemp\tLJs\n"
		file_unsrt = open("./Data/"+self.nametag+"_usrtData.log","w")
		file_unsrt.write(string_2+"\n"+exp_unsrt_string)
		file_unsrt.close()
		return Log_string,exp_unsrt_string,solver_log,calc_cholesky,zeta_string
	
	
	def getUncertainty(self,T): 
		unsrt = {}
		T1 = T["LJe"]
		T2 = T["LJs"]
		unsrt_LJe = np.sqrt((float(self.uncertainties["LJe"][0])*(T1/T1))**2)
		unsrt_LJs = np.sqrt((float(self.uncertainties["LJs"][0])*(T2/T2))**2)
		unsrt["LJe"] = unsrt_LJe
		unsrt["LJs"] = unsrt_LJs
		return unsrt
	
	def unsrt(self,index):
		if index == "LJe":
			L = self.uncertainties[index][0]
			z = 1
		if index == "LJs":
			L = self.uncertainties[index][0]
			z = 1
		#print(L)
		return L,z		
	

class collision:
	def __init__(self, Element,string,mechanism_loc):
		self.rxn = self.classification = self.type = self.sub_type = self. branching = self.branches  = self.pressure_limit  = self. common_temp = self.temp_limit= None
		
		self.tag = Element.tag
		self.cholskyDeCorrelateMat = {}
		self.zeta = {}
		self.exp_data_type = {}
		self.temperatures = {} 
		self.uncertainties = {}
		self.nominal = {}
		self.rxn = Element.attrib["rxn"]
		
		
		for item in Element:
			#print(item.tag)
			if item.tag == "class":
				self.classification = item.text
				continue
			if item.tag == "type":
				self.type = item.text
				continue
			if item.tag == "sub_type":
				self.sub_type = item.attrib["name"]
				
				for subitem in item:
					if "multiple" in subitem.tag:
						self.branching = str(subitem.text)
						continue
					if "branches" in subitem.tag:
						self.branches = str(subitem.text)
						self.m_dict = IFR.MechParsing(mechanism_loc).getThirdBodyCollisionEff(self.rxn,self.branches)
						continue
						
					if "pressure_limit" in subitem.tag:
						self.pressure_limit = str(subitem.text)
						continue
					if "common_temp" in subitem.tag:
						self.common_temp = str(subitem.text)
						continue
					if  "temp_limit" in subitem.tag :
						self.temp_limit = str(subitem.text)
						
						continue
				continue
					
		for i in self.branches.split(","):
			
			for item in Element:
				if item.tag == "data_type_"+str(i):

					self.exp_data_type[str(i)] = item.text

					continue
					
				if item.tag == "temp_"+str(i):
					#print(item.text)
					self.temperatures[str(i)] = np.asarray([float(j) for j in item.text.split(",")])
					continue
						
				if item.tag == "unsrt_"+str(i):
					self.uncertainties[str(i)] = np.asarray([float(i) for i in item.text.split(",")])
					continue
		
		
		for i in self.exp_data_type:
			
			if self.exp_data_type[str(i)].split(";")[0] ==  "constant":
				
				if self.exp_data_type[str(i)].split(";")[1] == "array":
				
					continue
				
				elif self.exp_data_type[str(i)].split(";")[1] == "end_points":
				
					self.temperatures[str(i)] = np.linspace(self.temperatures[str(i)][0],self.temperatures[str(i)][1],20)
				
					self.uncertainties[str(i)] = np.linspace(self.uncertainties[str(i)][0],self.uncertainties[str(i)][1],20)
				
					continue

			elif self.exp_data_type[str(i)].split(";")[0] == "percentage":
				
				if self.exp_data_type[str(i)].split(";")[1] == "array":
				
					a = self.trans_dict[str(i)]
					y = self.uncertainties[str(i)]						
					self.uncertainties[str(i)] = y*a
					
					continue

				elif self.exp_data_type[str(i)].split(";")[1] == "end_points":

					self.temperatures[str(i)] = np.linspace(self.temperatures[str(i)][0],self.temperatures[str(i)][1],20)
					self.uncertainties[str(i)] = np.linspace(self.uncertainties[str(i)][0],self.uncertainties[str(i)][1],20)

					a = self.trans_dict[str(i)]
					y = self.uncertainties[str(i)]						
					self.uncertainties[str(i)] = y*a
					continue
				
			
			
		
		for i in self.branches.split(","):
			self.cholskyDeCorrelateMat[str(i)],self.zeta[str(i)] = self.unsrt(self.uncertainties[str(i)])
	
		
		self.nametag = self.rxn+":"+self.sub_type
		self.solution = self.cholskyDeCorrelateMat
		self.f = self.getUncertainty()
		#print(self.m_dict)
		'''
		fig = plt.figure()
		plt.xlabel('Temperature (K)')
		plt.ylabel('Uncertainity ($\sigma$)')
		plt.title( string+'{}'.format(self.name), fontsize = 10)		
		plt.plot(self.temperatures,self.uncertainties,'o',label='exp uncertainties');
		plt.ylim(0,2*max(self.uncertainties[0],self.uncertainties[-1]))	
		plt.legend()
		my_path = os.getcwd()
		plt.savefig(my_path+'/Plots/TDeptUnsrt/'+self.name+'.png')
		
		'''	
	def getDtList(self):
		self.Unsrt_dict = {}
		key = ["Tag","Solution","Class","Type","Sub_type","Branch_boolen","Branches","Pressure_limit","Common_temp","temp_list","Nominal","Exp_input_data_type","priorCovariance","Basis_vector","Uncertainties","Temperatures","unsrt_func","Data_key"]
		values = [self.tag,self.solution,self.classification,self.type,self.sub_type,self.branching,self.branches,self.pressure_limit,self.common_temp,self.temp_limit,self.m_dict,self.exp_data_type,self.cholskyDeCorrelateMat,self.zeta,self.uncertainties,self.temperatures,self.f,self.branches]
		for i,element in enumerate(key):
			self.Unsrt_dict[element] = values[i]
		return self.Unsrt_dict
	
	def getAllData(self):
		Log_string = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.classification,self.type,self.sub_type,self.branching,self.branches,self.pressure_limit,self.common_temp,self.temp_limit,self.exp_data_type,self.nametag)
		exp_unsrt_string = ""
		solver_log = "{}\n{}\n".format(self.nametag,self.sol)
		calc_cholesky = self.cholskyDeCorrelateMat
		zeta_string = "{}".format(self.zeta)
		for i in range(len(self.temperatures)):
			exp_unsrt_string += "{}\t{}\n".format(self.temperatures[i],self.uncertainties[i])
		return Log_string,exp_unsrt_string,solver_log,calc_cholesky,zeta_string
		
	def getUncertainty(self): 
		unsrt = {}
		
		for i in self.branches.split(","):
			#print(i)
			T = self.temperatures[str(i)]
			unsrt_i = np.sqrt((float(self.uncertainties[str(i)][0])*(T/T))**2)
			unsrt[str(i)] = unsrt_i
		return unsrt
	def unsrt(self,unsrt):
		L = np.array([unsrt[0]])
		zeta = 1
		return L,zeta
	
#Parent class to host the uncertainty data for all the reactions. Has instances of reaction class as members        
class uncertaintyData:
	def __init__(self,xmlPath,mechPath,thermoPath,transportPath,unsrt_type="prior"):
		self.xmlPath = xmlPath
		print(xmlPath)
		self.Reactionlist = []
		self.focList = []
		self.mList = []
		self.thermoList = []
		self.transportList = []
		self.unsrt_data = {}
		self.reactionUnsrt = {}
		self.focUnsrt = {}
		self.trdBodyUnsrt = {}
		self.thermoUnsrt = {}
		self.transportUnsrt = {}
		self.tree = ET.parse(self.xmlPath)
		self.root = self.tree.getroot()
		count_rxn = 0
		count_foc = 0
		count_m = 0
		count_thermo = 0
		count_transport = 0
		
		for child in self.root:
			if child.tag == 'reaction':
				r = reaction(child,mechPath)
				self.Reactionlist.append(r.nametag)
				self.reactionUnsrt[r.nametag] = r
				self.unsrt_data[r.nametag] = r
				count_rxn +=1
			if child.tag == 'fallOffCurve':
				foc = fallOffCurve(child,mechPath)
				self.focList.append(foc.nametag)
				self.focUnsrt[foc.nametag] = foc
				self.unsrt_data[foc.nametag] = foc
				count_foc +=1
			if child.tag == 'thermo':
				th = thermodynamic(child,thermoPath)
				self.thermoList.append(th.nametag)
				self.thermoUnsrt[th.nametag] = th
				self.unsrt_data[th.nametag] = th
				count_thermo +=1
			if child.tag == 'collisionEff':
				string = "Unsrt for third bodies\n collision efficiencies [M]:  "
				m = collision(child,string,mechPath)
				self.mList.append(m.nametag)
				self.trdBodyUnsrt[m.nametag] = m
				self.unsrt_data[m.nametag] = m
				count_m +=1
			if child.tag == 'transport':
				tr = transport(child,transportPath)
				self.transportUnsrt[tr.nametag] = tr
				self.transportList.append(tr.nametag)
				self.unsrt_data[tr.nametag] = tr
				count_transport +=1
				#print(self.transportList)
		if unsrt_type !="opt":
			print("\n\n{} Reactions are selected for optimization\n".format(count_rxn))
			print("{} Fall-off (center broadening factors) of reactions {} are selected for optimization\n\n".format(count_foc,self.focList))
			print("{} third body collision efficiency's are selected for optimization\n\n".format(count_m))
			print("{} thermo-chemical parameters are selected for optimization\n\n".format(count_thermo))
			print("{} transport parameters are selected for optimization\n\n".format(count_transport))
            
	def extract_uncertainty(self):
		#print(self.root)
		return self.unsrt_data,self.reactionUnsrt, self.focUnsrt, self.trdBodyUnsrt, self.thermoUnsrt, self.transportUnsrt, self.Reactionlist,self.focList,self.mList,self.thermoList,self.transportList
	
        
