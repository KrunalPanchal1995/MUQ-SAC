import numpy as np
from solution import Solution
from scipy.optimize import minimize 
import numpy as np
import time
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib import style
import pygad
style.use("fivethirtyeight")
import os
from scipy.optimize import rosen, differential_evolution
from scipy.optimize import NonlinearConstraint, Bounds
from scipy.optimize import shgo
from scipy.optimize import BFGS
import pickle

class OptimizationTool(object):
	def __init__(self,
		target_list=None,frequency = None):
		
		self.target_list = target_list
		self.objective = 0
		self.frequency = frequency
		self.count = 0
	def obj_func_of_selected_PRS(self,x):
		
		kappa_curve = {}
		count = 0
		for i in self.rxn_index:
			temp = []
			for j in range(len(self.T)):
				temp.append(x[count])
				count += 1
			#Kappa = np.exp(np.log(kappa_0[i]) + temp*(np.log(kappa_max[i])-np.log(kappa_0[i])))
			Kappa = self.kappa_0[i] + temp*(self.kappa_max[i]-self.kappa_0[i])
			
			kappa_curve[i] = np.asarray(Kappa).flatten()
			
		
		zeta = {}
		#print(gen)
		for i in self.rxn_index:
			zeta[i] = self.rxn_unsrt_data[i].getZeta_typeA(kappa_curve[i])
		#print(zeta)
		x_transformed = []
		string = ""
		for i in self.rxn_index:
			temp = list(zeta[i])
			for k in temp:
				string+=f"{k},"
			x_transformed.extend(temp)
		string+=f"\n"
		x_transformed = np.asarray(x_transformed)
		zeta_file = open("zeta_guess_values.txt","+a").write(string)
		"""
		Just for simplicity
		"""
		x = x_transformed
		
		obj = 0.0
		rejected_PRS = []
		rejected_PRS_index = []
		target_value = []
		target_stvd = []
		direct_target_value = []
		direct_target_stvd = []
		target_value_2 = []
		case_stvd = []
		case_systematic_error = []
		response_value = []
		response_stvd = []	
		target_weights = []
		COUNT_Tig = 0
		COUNT_Fls = 0
		COUNT_All = 0	
		COUNT_Flw = 0
		D_SET = []
		#for case in self.target_list:
		#	D_SET.append(case.d_set)
		#	for i in self.frequency:
		#		COUNT_All += int(self.frequency[case.d_set])
		#print(self.selected_prs)
		frequency = {}
		for i,case in enumerate(self.target_list):
			if self.selected_prs[i] == 1:	
				if case.target == "Tig":
					if case.d_set in frequency:
						frequency[case.d_set] += 1
					else:
						frequency[case.d_set] = 1
					COUNT_All +=1
					COUNT_Tig +=1
					#dataset_weights = (1/len(self.frequency))*float(self.frequency[case.d_set])
					#dataset_weights = (1/COUNT_All)
					val = case.calculated_target_value(x_transformed)
					#val,grad = case.evaluateResponse(x)
					response_value.append(val)
					#response_stvd.append(grad)
					target_value.append(np.log(case.observed*10))
					target_value_2.append(np.log(case.observed))
					target_stvd.append(1/(np.log(case.std_dvtn*10)))
					case_stvd.append(np.log(case.std_dvtn*10))
					case_systematic_error.append(abs(np.log(case.observed*10)-val))
					#target_weights.append(dataset_weights)				
				elif case.target == "Fls":
					if case.d_set in frequency:
						frequency[case.d_set] += 1
					else:
						frequency[case.d_set] = 1
					COUNT_All +=1
					COUNT_Fls +=1
					#dataset_weights = (1/len(self.frequency))*float(self.frequency[case.d_set])
					#dataset_weights = (1/COUNT_All)
					val = np.exp(case.calculated_target_value(x_transformed))
					response_value.append(val)
					#response_stvd.append(grad)
					target_value.append(case.observed)
					target_value_2.append(case.observed)
					target_stvd.append(1/(case.std_dvtn))
					case_stvd.append(case.std_dvtn)
					case_systematic_error.append(abs(case.observed)-val)
					#target_weights.append(dataset_weights)	
				elif case.target == "Flw":
					if case.d_set in frequency:
						frequency[case.d_set] += 1
					else:
						frequency[case.d_set] = 1
					COUNT_All +=1
					COUNT_Flw +=1
					#dataset_weights = (1/len(self.frequency))*float(self.frequency[case.d_set])
					#dataset_weights = (1/COUNT_All)
					
					val = case.calculated_target_value(x_transformed)
					response_value.append(val)
					#response_stvd.append(grad)
					target_value.append(np.log(case.observed))
					target_value_2.append(np.log(case.observed))
					target_stvd.append(1/(np.log(case.std_dvtn)+abs(np.log(case.observed)-val)))
					case_stvd.append(np.log(case.std_dvtn))
					case_systematic_error.append(abs(np.log(case.observed)-val))
					#target_weights.append(dataset_weights)	
			else:
				#rejected_PRS.append(case)
				
		self.count +=1
		#print(target_value[0])
		#print(response_value[0])
		#print(target_stvd)
		diff = np.asarray(response_value)-np.asarray(target_value)
		#multiplicating_factors = np.asarray(target_weights)*np.asarray(target_stvd)
		#multiplicating_factor = 1/COUNT_All
		#print(sum(target_weights))
		multiplicating_factors = []
		#multiplicating_factors = np.asarray(target_weights)*np.asarray(target_stvd)
		for i,case in enumerate(self.target_list):
			if self.selected_prs[i] == 1:	
				if case.target == "Tig":
					multiplicating_factors.append(0.7*(1/COUNT_Tig))
		
				elif case.target == "Fls":
					multiplicating_factors.append(0.3*(1/COUNT_Fls))
		
		multiplicating_factors= np.asarray(multiplicating_factors)	
		
		for i,dif in enumerate(diff):
			obj+= multiplicating_factors[i]*(dif)**2
			#obj+= multiplicating_factor*(dif)**2
			
		#multiplicating_factors = np.asarray(target_weights)*np.asarray(target_stvd)
		#f = diff**2*multiplicating_factors
		#print(response_value[0])
		#print(target_value[0])
		#f = np.asarray(np.dot(f,np.asarray(target_weights))).flatten()
		#df = np.asarray(response_stvd)*np.asarray(target_stvd)
		#obj = np.dot(f,f)
		#print(f"\n\tResidual={obj}")
		
		"""
		Penalty function
		"""
		"""
		for i in x:
			obj+=(3/len(x))*(0.25*i**2)
		self.objective = obj
		"""
		#print(f"obj_2 = {obj_2}")

		get_systematic_error = open("systematic_error.txt","+a").write(f"{case_stvd},{case_systematic_error}\n")
		note = open("guess_values.txt","+a").write(f"{self.count},{x}\n")
		get_target_value = open("response_values.txt","+a").write(f"\t{target_value},{response_value}\n")
		
		return obj
		
	def fitness_function_for_T_INDIPENDENT(self):
		global fitness_func_T_indi
		def fitness_func_T_indi(x,solution_idx):
			string = ""
			for i in x:
				string+=f"{i},"
			string+=f"\n"
			#x_transformed = np.asarray(x_transformed)
			zeta_file = open("zeta_guess_values.txt","+a").write(string)
			"""
			Just for simplicity
			"""
			#x = x_transformed
			
			obj = 0.0
			rejected_PRS = []
			rejected_PRS_index = []
			target_value = []
			target_stvd = []
			direct_target_value = []
			direct_target_stvd = []
			target_value_2 = []
			case_stvd = []
			case_systematic_error = []
			response_value = []
			response_stvd = []	
			target_weights = []	
			COUNT_Tig = 0
			COUNT_Fls = 0
			COUNT_All = 0	
			COUNT_Flw = 0
			frequency = {}
			#for i in self.frequency:
			#	COUNT_All += int(self.frequency[case.d_set])
			#print(self.selected_prs)
			#raise AssertionError("Stop!!")
			for i,case in enumerate(self.target_list):
				if self.selected_prs[i] == 1:	
					if case.target == "Tig":
						if case.d_set in frequency:
							frequency[case.d_set] += 1
						else:
							frequency[case.d_set] = 1
						COUNT_All +=1
						COUNT_Tig +=1
						#dataset_weights = (1/len(self.frequency))*float(self.frequency[case.d_set])
						#dataset_weights = (1/COUNT_All)
						
						val = case.calculated_target_value(x)
						#print(val)
						#val,grad = case.evaluateResponse(x)
						response_value.append(val)
						#response_stvd.append(grad)
						target_value.append(np.log(case.observed*10))
						target_value_2.append(np.log(case.observed))
						target_stvd.append(1/(np.log(case.std_dvtn*10)))
						case_stvd.append(np.log(case.std_dvtn*10))
						case_systematic_error.append(abs(np.log(case.observed*10)-val))
						#target_weights.append(dataset_weights)				
					elif case.target == "Fls":
						if case.d_set in frequency:
							frequency[case.d_set] += 1
						else:
							frequency[case.d_set] = 1
						COUNT_All +=1
						COUNT_Fls +=1
						#dataset_weights = (1/len(self.frequency))*float(self.frequency[case.d_set])
						#dataset_weights = (1/COUNT_All)
						
						val = np.exp(case.calculated_target_value(x))
						response_value.append(val)
						#response_stvd.append(grad)
						target_value.append(case.observed)
						target_value_2.append(case.observed)
						target_stvd.append(1/(case.std_dvtn))
						case_stvd.append(case.std_dvtn)
						case_systematic_error.append(abs(case.observed)-val)
						#target_weights.append(dataset_weights)	
					elif case.target == "Flw":
						if case.d_set in frequency:
							frequency[case.d_set] += 1
						else:
							frequency[case.d_set] = 1
						COUNT_All +=1
						COUNT_Flw +=1
						#dataset_weights = (1/len(self.frequency))*float(self.frequency[case.d_set])
						#dataset_weights = (1/COUNT_All)
						
						val = case.calculated_target_value(x)
						response_value.append(val)
						#response_stvd.append(grad)
						target_value.append(np.log(case.observed))
						target_value_2.append(np.log(case.observed))
						target_stvd.append(1/(np.log(case.std_dvtn)+abs(np.log(case.observed)-val)))
						case_stvd.append(np.log(case.std_dvtn))
						case_systematic_error.append(abs(np.log(case.observed)-val))
						#target_weights.append(dataset_weights)	
				
			#print(frequency)
			self.count +=1
			#print(target_value[0])
			#print(response_value[0])
			#print(target_stvd)
			diff = np.asarray(response_value)-np.asarray(target_value)
			multiplicating_factors = []
			#multiplicating_factors = np.asarray(target_weights)*np.asarray(target_stvd)
			for i,case in enumerate(self.target_list):
				if self.selected_prs[i] == 1:	
					if case.target == "Tig":
						multiplicating_factors.append(1/COUNT_Tig)
			
					elif case.target == "Fls":
						multiplicating_factors.append(0.05*(1/COUNT_Fls))
			
			multiplicating_factors= np.asarray(multiplicating_factors)		
			#Giving all datapoints equal weights
			#multiplicating_factor = 1/COUNT_All
			
			for i,dif in enumerate(diff):
				#obj+= multiplicating_factors[i]*(dif)**2
				#obj+= multiplicating_factor*(dif)**2
				obj+= multiplicating_factors[i]*(dif)**2
							
			"""
			Penalty function
			"""
			"""
			for i in x:
				obj+=(3/len(x))*(0.25*i**2)
			self.objective = obj
			"""
			#print(f"obj_2 = {obj_2}")

			get_systematic_error = open("systematic_error.txt","+a").write(f"{case_stvd},{case_systematic_error}\n")
			note = open("guess_values.txt","+a").write(f"{self.count},{x}\n")
			get_target_value = open("response_values.txt","+a").write(f"\t{target_value},{response_value}\n")
			
			#note = open("target_values.txt","+a").write(f"{self.count},{obj},{target_stvd}\n")
			fitness = 1.0 / (np.abs((obj) - 0) + 0.000001)
			record =open("samplefile.txt","+a").write(f"{self.ga_instance.generations_completed},{self.count},{self.objective},{fitness}\n")
			#print(fitness)
			return fitness
		return fitness_func_T_indi
		
	def run_optimization_with_selected_PRS(self,
					       simulator,
					       selectedPRS,
					       method="gradient-based",
					       algorithm="SLSQP",
					       initial_guess=None,
					       bounds=None,
					       bounds_array = None,
					       initial_covariance=None):
	   
		self.simulator = simulator
		self.selected_prs = selectedPRS
		self.solution = Solution(initial_guess,
					 covariance_x=initial_covariance,
					 initial_x=initial_guess,
					 initial_covariance=initial_covariance)
	 
		if "solution.save" in os.listdir():
			print("Optimization is already finished\n")
			save = open("solution.save","r").readlines()
			save_zeta = open("solution_zeta.save","r").readlines()
			optimal_parameters = []
			optimal_parameters_zeta = []
			for line in save:
				optimal_parameters.append(float(line.split("=")[1].strip("\n")))
			
			for line in save:
				optimal_parameters_zeta.append(float(line.split("=")[1].strip("\n")))
			
			
			
		else:	
			if "gradient-based" in method:
				
				self.rxn_index = self.simulator.ind		
				
				
				fitness_function = self.fitness_function_factory()
				
				
				self.rxn_index = self.simulator.ind
				self.rxn_unsrt_data = self.simulator.rxnUnsert
				cholesky_dict = {}
				self.activeParameters = {}
				
				self.kappa_0 = {}
				self.kappa_max = {}
				#T = np.array([300,1500,3500])
				self.T = np.linspace(300,2500,50)
				self.init_guess = np.zeros(len(self.T)*len(self.rxn_index))
				T = np.linspace(300,2500,50)
				theta = np.array([T/T,np.log(T),-1/T])
				#self.theta_inv = np.linalg.inv(theta.T)
				for i in self.rxn_index:
					self.activeParameters[i] = self.rxn_unsrt_data[i].activeParameters
					self.kappa_0[i] = self.rxn_unsrt_data[i].getNominal(T)
					self.kappa_max[i] = self.rxn_unsrt_data[i].getKappaMax(T)
				
				#bnds = (1/3)*np.ones(len(initial_guess))
				bounds = Bounds(list(-np.ones(len(self.init_guess))),
				
				
				
				opt_output = minimize(self.obj_func_of_selected_PRS,initial_guess,bounds=bounds,method=algorithm)
				print(opt_output)
				
				optimal_parameters = np.asarray(opt_output.x)
				"""
				Finding the optimal parameters
				"""
				print("<<<<<<<<<<<<<<<<FOUND BEST SOLUTION>>>>>>>>>>>>>>>>>>>>>\n")
				print(f"{solution}")
				kappa_curve = {}
				count = 0
				for i in rxn_index:
					temp = []
					for j in range(len(self.T)):
						temp.append(optimal_parameters[count])
						count += 1
					#Kappa = np.exp(np.log(kappa_0[i]) + temp*(np.log(kappa_max[i])-np.log(kappa_0[i])))
					Kappa = kappa_0[i] + temp*(kappa_max[i]-kappa_0[i])
					#print(Kappa)
					kappa_curve[i] = np.asarray(Kappa).flatten()
					
				
				zeta = {}
				#print(gen)
				for i in rxn_index:
					zeta[i] = rxn_unsrt_data[i].getZeta_typeA(kappa_curve[i])
				optimal_parameters_zeta = []
				for i in rxn_index:
					temp = list(zeta[i])
					optimal_parameters_zeta.extend(temp)
				optimal_parameters_zeta = np.asarray(optimal_parameters_zeta)	
			
				
				
				
				"""
				-----------------------------------------
				Transformation 1.0 of the optimal parameters:
				-----------------------------------------
				
				     X =  __ln(kappa/kappa_0)____
					   ln(kappa_max/kappa_0)
				    
				  ln(kappa) - ln(kappa_0) = X*[ln(kappa_max) - ln(kappa_0)]
				  
				  ln(kappa) = ln(kappa_0) + X*[ln(kappa_max) - ln(kappa_0)]
				     
				-----------------------------------------
				Transformation 2.0 of the optimal parameters:
				-----------------------------------------
				
				     X =  ___kappa__-__kappa_0__
					    kappa_max - kappa_0
				    
				  kappa - kappa_0 = X*[kappa_max - kappa_0]
				  
				  kappa = kappa_0 + X*[kappa_max - kappa_0]
				
				
				"""
				
		return np.asarray(optimal_parameters),np.asarray(optimal_parameters_zeta)
	
