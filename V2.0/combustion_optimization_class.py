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

	
	def _obj_fun(self,x):
		num_params = len(x)
		num_expts = len(self.target_list)
		
		f = np.empty(num_params + num_expts)
		df = np.zeros((num_params + num_expts,num_params))
		inv_covar = np.matrix(self.solution.alpha_i)
		initial_guess = self.solution.x_i
		#print(np.shape(inv_covar))
		#Set the parts of the objective function that depend on x
		f[0:num_params] = np.dot(inv_covar,(x - initial_guess))
		df[0:num_params,0:num_params] = inv_covar
		for i,case in enumerate(self.target_list):
			f_num,df_num = case.estimate(x)
			f_exp = case.observed
			w = 1/case.std_dvtn
			f[i + num_params] = (f_num - f_exp)*w
			df[i + num_params,:] = df_num*w
		print(np.shape(df))
		return f,df
	
	#def constr_R2(x):	
	#	return
	
	#def constr_R3(x):
	#	return
	
	#def constr_R4(x):
	#	return
		
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
				"""
				rejected_PRS_index.append(i)
				if case.target == "Tig":
					#response_value.append(simulator.simulate_target_value(case,x,self.iter_number))
					response_stvd.append(0)
					direct_target_value.append(np.log(case.observed*10))
					target_value_2.append(np.log(case.observed*10))
					direct_target_stvd.append(np.log(case.std_dvtn*10))
					target_weights.append(case.d_weight)
					case_stvd.append(np.log(case.std_dvtn))
					
					#target.append(case)
				elif case.target == "Fls":
					#response_value.append(case.simulate_target_value(case,x,self.iter_number))
					response_stvd.append(0)
					direct_target_value.append(np.log(case.observed))
					target_value_2.append(np.log(case.observed))
					direct_target_stvd.append(np.log(case.std_dvtn))
					case_stvd.append(np.log(case.std_dvtn))
					target_weights.append(case.d_weight)
				elif case.target == "Flw":
					#response_value.append(case.simulate_target_value(case,x,self.iter_number))
					response_stvd.append(0)
					direct_target_value.append(np.log(case.observed))
					target_value_2.append(np.log(case.observed))
					direct_target_stvd.append(np.log(case.std_dvtn))
					case_stvd.append(np.log(case.std_dvtn))
					target_weights.append(case.d_weight)
		if len(rejected_PRS) != 0:
			directResponses = self.simulator.getSimulatedValues(x,rejected_PRS_index,rejected_PRS,self.count,self.objective)
			response_value.extend(directResponses)
			direct_target_stvd = np.asarray(direct_target_stvd)+abs(np.asarray(directResponses)-np.asarray(direct_target_value))
			target_value.extend(list(direct_target_value))
			target_stvd.extend(list(1/direct_target_stvd))
			systematic_error = np.asarray(direct_target_value) - np.asarray(response_value)
			case_systematic_error.extend(list(systematic_error))	
		"""
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

	def callback(self):
		print("Generation : ", self.ga_instance.generations_completed)
		print("Fitness of the best solution :", self.ga_instance.best_solution()[1])
	
	def objective_func_T_indipendent(self,x):
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
			else:
				#rejected_PRS.append(case)
				"""
				rejected_PRS_index.append(i)
				if case.target == "Tig":
					#response_value.append(simulator.simulate_target_value(case,x,self.iter_number))
					response_stvd.append(0)
					direct_target_value.append(np.log(case.observed*10))
					target_value_2.append(np.log(case.observed*10))
					direct_target_stvd.append(np.log(case.std_dvtn*10))
					target_weights.append(case.d_weight)
					case_stvd.append(np.log(case.std_dvtn))
					
					#target.append(case)
				elif case.target == "Fls":
					#response_value.append(case.simulate_target_value(case,x,self.iter_number))
					response_stvd.append(0)
					direct_target_value.append(np.log(case.observed))
					target_value_2.append(np.log(case.observed))
					direct_target_stvd.append(np.log(case.std_dvtn))
					case_stvd.append(np.log(case.std_dvtn))
					target_weights.append(case.d_weight)
				elif case.target == "Flw":
					#response_value.append(case.simulate_target_value(case,x,self.iter_number))
					response_stvd.append(0)
					direct_target_value.append(np.log(case.observed))
					target_value_2.append(np.log(case.observed))
					direct_target_stvd.append(np.log(case.std_dvtn))
					case_stvd.append(np.log(case.std_dvtn))
					target_weights.append(case.d_weight)
		if len(rejected_PRS) != 0:
			directResponses = self.simulator.getSimulatedValues(x,rejected_PRS_index,rejected_PRS,self.count,self.objective)
			response_value.extend(directResponses)
			direct_target_stvd = np.asarray(direct_target_stvd)+abs(np.asarray(directResponses)-np.asarray(direct_target_value))
			target_value.extend(list(direct_target_value))
			target_stvd.extend(list(1/direct_target_stvd))
			systematic_error = np.asarray(direct_target_value) - np.asarray(response_value)
			case_systematic_error.extend(list(systematic_error))	
		"""
		#print(frequency)
		self.count +=1
		#print(target_value[0])
		#print(response_value[0])
		#print(target_stvd)
		diff = (np.asarray(response_value)-np.asarray(target_value))/np.asarray(case_stvd)
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
		#fitness = 1.0 / (np.abs((obj) - 0) + 0.000001)
		#record =open("samplefile.txt","+a").write(f"{self.ga_instance.generations_completed},{self.count},{self.objective},{fitness}\n")
		#print(fitness)
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
				else:
					#rejected_PRS.append(case)
					"""
					rejected_PRS_index.append(i)
					if case.target == "Tig":
						#response_value.append(simulator.simulate_target_value(case,x,self.iter_number))
						response_stvd.append(0)
						direct_target_value.append(np.log(case.observed*10))
						target_value_2.append(np.log(case.observed*10))
						direct_target_stvd.append(np.log(case.std_dvtn*10))
						target_weights.append(case.d_weight)
						case_stvd.append(np.log(case.std_dvtn))
						
						#target.append(case)
					elif case.target == "Fls":
						#response_value.append(case.simulate_target_value(case,x,self.iter_number))
						response_stvd.append(0)
						direct_target_value.append(np.log(case.observed))
						target_value_2.append(np.log(case.observed))
						direct_target_stvd.append(np.log(case.std_dvtn))
						case_stvd.append(np.log(case.std_dvtn))
						target_weights.append(case.d_weight)
					elif case.target == "Flw":
						#response_value.append(case.simulate_target_value(case,x,self.iter_number))
						response_stvd.append(0)
						direct_target_value.append(np.log(case.observed))
						target_value_2.append(np.log(case.observed))
						direct_target_stvd.append(np.log(case.std_dvtn))
						case_stvd.append(np.log(case.std_dvtn))
						target_weights.append(case.d_weight)
			if len(rejected_PRS) != 0:
				directResponses = self.simulator.getSimulatedValues(x,rejected_PRS_index,rejected_PRS,self.count,self.objective)
				response_value.extend(directResponses)
				direct_target_stvd = np.asarray(direct_target_stvd)+abs(np.asarray(directResponses)-np.asarray(direct_target_value))
				target_value.extend(list(direct_target_value))
				target_stvd.extend(list(1/direct_target_stvd))
				systematic_error = np.asarray(direct_target_value) - np.asarray(response_value)
				case_systematic_error.extend(list(systematic_error))	
			"""
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
	

	def fitness_function_factory(self):
		global fitness_func
		def fitness_func(x,solution_idx):
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
			frequency = {}
			#for i in self.frequency:
			#	COUNT_All += int(self.frequency[case.d_set])
			#print(self.selected_prs)
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
					"""
					rejected_PRS_index.append(i)
					if case.target == "Tig":
						#response_value.append(simulator.simulate_target_value(case,x,self.iter_number))
						response_stvd.append(0)
						direct_target_value.append(np.log(case.observed*10))
						target_value_2.append(np.log(case.observed*10))
						direct_target_stvd.append(np.log(case.std_dvtn*10))
						target_weights.append(case.d_weight)
						case_stvd.append(np.log(case.std_dvtn))
						
						#target.append(case)
					elif case.target == "Fls":
						#response_value.append(case.simulate_target_value(case,x,self.iter_number))
						response_stvd.append(0)
						direct_target_value.append(np.log(case.observed))
						target_value_2.append(np.log(case.observed))
						direct_target_stvd.append(np.log(case.std_dvtn))
						case_stvd.append(np.log(case.std_dvtn))
						target_weights.append(case.d_weight)
					elif case.target == "Flw":
						#response_value.append(case.simulate_target_value(case,x,self.iter_number))
						response_stvd.append(0)
						direct_target_value.append(np.log(case.observed))
						target_value_2.append(np.log(case.observed))
						direct_target_stvd.append(np.log(case.std_dvtn))
						case_stvd.append(np.log(case.std_dvtn))
						target_weights.append(case.d_weight)
			if len(rejected_PRS) != 0:
				directResponses = self.simulator.getSimulatedValues(x,rejected_PRS_index,rejected_PRS,self.count,self.objective)
				response_value.extend(directResponses)
				direct_target_stvd = np.asarray(direct_target_stvd)+abs(np.asarray(directResponses)-np.asarray(direct_target_value))
				target_value.extend(list(direct_target_value))
				target_stvd.extend(list(1/direct_target_stvd))
				systematic_error = np.asarray(direct_target_value) - np.asarray(response_value)
				case_systematic_error.extend(list(systematic_error))	
			"""
			print(frequency)
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
				
					
			#print(sum(target_weights))
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
			
			#note = open("target_values.txt","+a").write(f"{self.count},{obj},{target_stvd}\n")
			fitness = 1.0 / (np.abs((obj) - 0) + 0.000001)
			record =open("samplefile.txt","+a").write(f"{self.ga_instance.generations_completed},{self.count},{self.objective},{fitness}\n")
			#print(fitness)
			return fitness
		return fitness_func
	
	def run_optimization(self,initial_guess=None,initial_covariance=None):
	   
		self.solution = Solution(initial_guess,
					 covariance_x=initial_covariance,
					 initial_x=initial_guess,
					 initial_covariance=initial_covariance)
	 
		opt_output = spopt.root(self._obj_fun,initial_guess,method='lm',jac=True)
		
		#print (opt_output.message)
		
		optimal_parameters = np.array(opt_output.x)
		
		residuals,final_jac = self._obj_fun(optimal_parameters)
		
		#print(np.shape(final_jac))
		
		
		icov = np.dot(final_jac.T,final_jac)
		cov = np.linalg.inv(icov)
		self.solution = Solution(optimal_parameters,
								 covariance_x=cov,
								 initial_x=initial_guess,
								 initial_covariance=initial_covariance)
		
		return optimal_parameters,cov	
	
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
			
			"""
			if "GA" in method:
				rxn_index = self.simulator.ind
				rxn_unsrt_data = self.simulator.rxnUnsert
				cholesky_dict = {}
				activeParameters = {}
				
				kappa_0 = {}
				kappa_max = {}
				#T = np.array([300,1500,3500])
				self.T = np.linspace(300,2500,20)#np.array([300,1000,1300,1500,1700,2000,2500,3500])
				
				T = np.linspace(300,2500,20)
				
				theta = np.array([T/T,np.log(T),-1/T])
				#theta_inv = np.linalg.inv(theta.T)
				for i in rxn_index:
					activeParameters[i] = rxn_unsrt_data[i].activeParameters
					kappa_0[i] = rxn_unsrt_data[i].getNominal(T)
					kappa_max[i] = rxn_unsrt_data[i].getKappaMax(T)

				kappa_curve = {}
				count = 0
				for i in rxn_index:
					temp = []
					for j in range(len(self.T)):
						temp.append(optimal_parameters[count])
						count += 1
					#Kappa = np.exp(np.log(kappa_0[i]) + temp*(np.log(kappa_max[i])-np.log(kappa_0[i])))
					Kappa = kappa_0[i] + temp*(kappa_max[i]-kappa_0[i])
					
					kappa_curve[i] = np.asarray(Kappa).flatten()
					
				
				zeta = {}
				#print(gen)
				for i in rxn_index:
					zeta[i] = rxn_unsrt_data[i].getZeta_typeA(kappa_curve[i])
				#print(zeta)
				optimal_parameters_zeta = []
				for i in rxn_index:
					temp = list(zeta[i])
					optimal_parameters_zeta.extend(temp)
				optimal_parameters_zeta = np.asarray(optimal_parameters_zeta)	
				
			else:
				rxn_index = self.simulator.ind
				rxn_unsrt_data = self.simulator.rxnUnsert
				cholesky_dict = {}
				activeParameters = {}
				zeta = {}
				for i in rxn_index:
					activeParameters[i] = rxn_unsrt_data[i].activeParameters
				gen = {}
				count = 0
				for i in rxn_index:
					temp = []
					for j in range(len(activeParameters[i])):
						temp.append(optimal_parameters[count])
						count += 1
					gen[i] = temp
				#print(gen)
				for i in rxn_index:
					zeta[i] = rxn_unsrt_data[i].getZetaFromGen(gen[i])
				optimal_parameters_zeta = []
				for i in rxn_index:
					temp = list(zeta[i])
					optimal_parameters_zeta.extend(temp)
				optimal_parameters_zeta = np.asarray(optimal_parameters_zeta)
			"""
			
		else:	
			if "gradient-based" in method:
				if "A-facto" in self.simulator.design:

					self.init_guess = np.zeros(len(self.rxn_index))
					opt_output = minimize(self.obj_func_of_selected_PRS,self.init_guess,bounds=bounds,method=algorithm)
					print(opt_output)
					
					optimal_parameters = np.asarray(opt_output.x)
					optimal_parameters_zeta = optimal_parameters
				
				else:
					opt_output = minimize(self.obj_func_of_selected_PRS,initial_guess,bounds=bounds,method=algorithm)
					print(opt_output)
					
					optimal_parameters = np.asarray(opt_output.x)
			
			elif "differential_evolution" in method:
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
				bounds = Bounds(list(-np.ones(len(self.init_guess))), list(np.ones(len(self.init_guess))))
				
				"""
				Finding the constraints parameters
				
				"""
				#print(self.const[0]["kappa_min"])
				
				print("Started the optimization algorithm")
				p = pickle.dumps(self.obj_func_of_selected_PRS)
				q = pickle.loads(p)
				solution = differential_evolution(q,bounds,strategy='best2bin', 
								   maxiter=1000, popsize=100, tol=0.01, mutation=(0.5, 1), 
								   recombination=0.7, callback=None, disp=True, 
								   polish=True, init='latinhypercube', atol=0, updating='immediate', 
								   workers=50)
								   
								   #nlc_r1,nlc_r2))
								   
								   #nlc_r3,nlc_r4,nlc_r5,nlc_r6,nlc_r7,nlc_r8,nlc_r9,nlc_r10,nlc_r11,nlc_r12,nlc_r13,nlc_r14,nlc_r15,nlc_r16,nlc_r17,nlc_r18,nlc_r19,nlc_r20,nlc_r21,nlc_r22,nlc_r23,nlc_r24,nlc_r25,nlc_r26,nlc_r27,nlc_r28,nlc_r29,nlc_r30))
				
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
				
			elif "GA" in method:
				self.rxn_index = self.simulator.ind
				#print(self.simulator.design)		
				if "A-facto" in self.simulator.design:
					fitness_function = self.fitness_function_for_T_INDIPENDENT()
					self.init_guess = np.zeros(len(self.rxn_index))
					gene_space = [{'low': -1, 'high': 1} for _ in self.init_guess ]
				else:
					#fitness_function = self.fitness_function_for_T_INDIPENDENT()
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
					gene_space = [{'low': -1, 'high': 1} for _ in self.init_guess ]
					T = np.linspace(300,2500,50)
					theta = np.array([T/T,np.log(T),-1/T])
					#self.theta_inv = np.linalg.inv(theta.T)
					for i in self.rxn_index:
						self.activeParameters[i] = self.rxn_unsrt_data[i].activeParameters
						self.kappa_0[i] = self.rxn_unsrt_data[i].getNominal(T)
						self.kappa_max[i] = self.rxn_unsrt_data[i].getKappaMax(T)
				
				
				self.ga_instance = pygad.GA(num_generations=2000,
			               num_parents_mating=300,
			               fitness_func=fitness_function,
			               init_range_low=-1,
			               init_range_high=1,
			               sol_per_pop=400,
			               num_genes=len(self.init_guess),
			               crossover_type="uniform",
			               crossover_probability=0.6,
			               mutation_type="adaptive",
			               mutation_probability=(0.03, 0.008),
			               gene_type=float,
			               allow_duplicate_genes=False,
			               gene_space=gene_space,
			               keep_parents = -1,
			               save_best_solutions=True,
			               save_solutions=True,
			               stop_criteria=["reach_300"])
			               
				self.ga_instance.run()
				self.ga_instance.plot_fitness(title="GA with Adaptive Mutation", linewidth=5)
				filename = 'genetic'
				self.ga_instance.save(filename=filename)
				optimal_parameters, solution_fitness, solution_idx = self.ga_instance.best_solution()
				
				
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
				
				
				
				#optimal_parameters_zeta = optimal_parameters
				if "A-facto" in self.simulator.design:
					optimal_parameters_zeta = optimal_parameters
					print("Parameters of the best solution : {solution}".format(solution=optimal_parameters))
					print("Fitness value of the best solution = {solution_fitness}".format(solution_fitness=solution_fitness))
					print("Index of the best solution : {solution_idx}".format(solution_idx=solution_idx))
					if self.ga_instance.best_solution_generation != -1:
						print("Best fitness value reached after {best_solution_generation} generations.".format(best_solution_generation=self.ga_instance.best_solution_generation))
					
				else:
				
					rxn_index = self.simulator.ind
					rxn_unsrt_data = self.simulator.rxnUnsert
					cholesky_dict = {}
					activeParameters = {}
					
					kappa_0 = {}
					kappa_max = {}
					
					for i in rxn_index:
						activeParameters[i] = rxn_unsrt_data[i].activeParameters
						kappa_0[i] = rxn_unsrt_data[i].getNominal(T)
						kappa_max[i] = rxn_unsrt_data[i].getKappaMax(T)
						
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
				
					print("Parameters of the best solution : {solution}".format(solution=optimal_parameters))
					print("Fitness value of the best solution = {solution_fitness}".format(solution_fitness=solution_fitness))
					print("Index of the best solution : {solution_idx}".format(solution_idx=solution_idx))
					if self.ga_instance.best_solution_generation != -1:
						print("Best fitness value reached after {best_solution_generation} generations.".format(best_solution_generation=self.ga_instance.best_solution_generation))
				#residuals,final_jac = self._obj_fun(optimal_parameters)
				
				#print(np.shape(final_jac))
				
				
				#icov = np.dot(final_jac.T,final_jac)
				#cov = np.linalg.inv(icov)
				#self.solution = Solution(optimal_parameters,
				#						 covariance_x=cov,
				#						 initial_x=initial_guess,
				#						 initial_covariance=initial_covariance)
				
		return np.asarray(optimal_parameters),np.asarray(optimal_parameters_zeta)
	
		
