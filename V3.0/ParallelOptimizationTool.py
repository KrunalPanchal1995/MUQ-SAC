import numpy as np
from solution import Solution
from scipy.optimize import minimize 
from scipy import optimize as spopt
from simulation_manager2_0 import SM as simulator
import numpy as np
import time
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib import style
import pygad
import matplotlib as mpl
mpl.use('Agg')
style.use("fivethirtyeight")
import os
from scipy.optimize import rosen, differential_evolution
from scipy.optimize import NonlinearConstraint, Bounds
from scipy.optimize import shgo
from scipy.optimize import BFGS
import pickle
import functools
class EarlyStopper:
    def __init__(self, patience=30):
        self.patience = patience
        self.best_fitness = None
        self.counter = 0

    def check(self,ga_instance):
        current_best_fitness = ga_instance.best_solution()[1]
        
        if self.best_fitness is None or current_best_fitness > self.best_fitness:
            self.best_fitness = current_best_fitness
            self.counter = 0  # Reset counter
        else:
            self.counter += 1  # Increase counter if no improvement
        
        if self.counter >= self.patience:
            print(f"Stopping early: No improvement for {self.patience} generations")
            ga_instance.stop()


class OptimizationTool(object):
	def __init__(self,
		targets_LTT=None,targets_HTT=None,frequency = None):
		
		self.targets_LTT = targets_LTT
		self.targets_HTT = targets_HTT
		self.objective = 0
		self.frequency = frequency
		self.count = 0
	
	def obj_func_parallel_optimization_with_selected_prs(self,x):
		"""
		-----------------------
		Just for simplicity
		-----------------------
		# Multi-Stage Optimization (MSO) with Reaction Classification (RC)
		- self.unsrt_LTC : Uncertainty Quantification of active reactions within the class of Low Temperature Chemistry (LTC)
		- self.unsrt_HTC : Uncertainty Quantification of active reactions within the class of High Temperature Chemistry (HTC)
		- self.target_LTT: Combustion target class for Low Temperature Tragets (LTT)
		- self.target_HTT: Combustion target class for High Temperature Tragets (HTT)
		- N_LTC			 : Number of reactions for LTC		 
		- N_HTC 	     : Number of reactions for HTC
		"""
		
		N_LTC = len([rxn for rxn in self.unsrt_LTC])
		N_HTC = len([rxn for rxn in self.unsrt_HTC])
		
		####################################
		### Writing the search Iterations
		####################################
		"""
		Later these search iterations can be used as DesignMatrix for
		testing the PRS predictions of targets against the black-box simulations on the 
		optimization search iterations
		
		For More Information refer: Krunal et al. 2024
		"""
		
		string = ""
		for i in x[0:N_LTC]:
			string+=f"{i},"
		string+=f"\n"
		#x_transformed = np.asarray(x_transformed)
		zeta_file = open("optimization_search_iters_HTC.txt","+a").write(string)
		
		string = ""
		for i in x[N_LTC:]:
			string+=f"{i},"
		string+=f"\n"
		#x_transformed = np.asarray(x_transformed)
		zeta_file = open("optimization_search_iters_LTC.txt","+a").write(string)
		
		string = ""
		for i in x:
			string+=f"{i},"
		string+=f"\n"
		#x_transformed = np.asarray(x_transformed)
		zeta_file = open("optimization_search_iters_combined.txt","+a").write(string)
		
		########################################
		### Initializing the objective functions
		###
		########################################
		
		
		
		obj_LTT = 0.0
		obj_HTT = 0.0
		
		target_value_LTT = []
		target_stvd_LTT = []
		response_value_LTT = []
		target_value_HTT = []
		target_stvd_HTT = []
		response_value_HTT = []

		
		COUNT_Tig_LTT = 0
		COUNT_All_LTT = 0	
		frequency_LTT = {}
		COUNT_Tig_HTT = 0
		COUNT_All_HTT = 0	
		frequency_HTT = {}
		
		for i,case in enumerate(self.targets_LTT):
			if self.ResponseSurfacesLTT[i].selection == 1:	
				if case.target == "Tig":
					if case.d_set in frequency_LTT:
						frequency_LTT[case.d_set] += 1
					else:
						frequency_LTT[case.d_set] = 1
					COUNT_All_LTT +=1
					COUNT_Tig_LTT +=1
				
					val_LTT = self.ResponseSurfacesLTT[i].evaluate(x[0:N_LTC])#PRS predictions for LTC on LTT
					
					response_value_LTT.append(val_LTT)
					target_value_LTT.append(np.log(case.observed*10))
					target_stvd_LTT.append(1/(np.log(case.std_dvtn*10)))			
		
		for i,case in enumerate(self.targets_HTT):
			if self.ResponseSurfacesHTT[i].selection == 1:	
				if case.target == "Tig":
					if case.d_set in frequency_HTT:
						frequency_HTT[case.d_set] += 1
					else:
						frequency_HTT[case.d_set] = 1
					COUNT_All_HTT +=1
					COUNT_Tig_HTT +=1

					val_HTT = self.ResponseSurfacesHTT[i].evaluate(x[N_LTC:])

					response_value_HTT.append(val_HTT)
					target_value_HTT.append(np.log(case.observed*10))
					target_stvd_HTT.append(1/(np.log(case.std_dvtn*10)))
					
			
		self.count +=1
		diff_LTT = np.asarray(response_value_LTT)-np.asarray(target_value_LTT)
		diff_HTT = np.asarray(response_value_HTT)-np.asarray(target_value_HTT)
		multiplicating_factors_LTT = []
		multiplicating_factors_HTT = []
		#multiplicating_factors = np.asarray(target_weights)*np.asarray(target_stvd)
		for i,case in enumerate(self.targets_LTT):
			if self.ResponseSurfacesLTT[i].selection == 1:	
				if case.target == "Tig":
					multiplicating_factors_LTT.append(1/COUNT_Tig_LTT)
		
		for i,case in enumerate(self.targets_HTT):
			if self.ResponseSurfacesHTT[i].selection == 1:	
				if case.target == "Tig":
					multiplicating_factors_HTT.append(1/COUNT_Tig_HTT)
		
		multiplicating_factors_LTT= np.asarray(multiplicating_factors_LTT)	
		multiplicating_factors_HTT= np.asarray(multiplicating_factors_HTT)		
		#Giving all datapoints equal weights
		#multiplicating_factor = 1/COUNT_All
		
		for i,dif in enumerate(diff_LTT):
			#obj_LTT+= multiplicating_factors_LTT[i]*(dif*target_stvd_LTT[i])**2
			obj_LTT+= multiplicating_factors_LTT[i]*(dif)**2
		
		for i,dif in enumerate(diff_HTT):
			#obj_HTT+= multiplicating_factors_HTT[i]*(dif*target_stvd_HTT[i])**2
			obj_HTT+= multiplicating_factors_HTT[i]*(dif)**2
		
		string_target_LTT = ""
		string_response_LTT = ""
		for i,value in enumerate(response_value_LTT):
			string_response_LTT+=f"{value},"
			string_target_LTT+=f"{target_value_LTT[i]},"
		string_response_LTT+="\n"
		string_target_LTT+="\n"
		string_target_HTT = ""
		string_response_HTT = ""
		for i,value in enumerate(response_value_HTT):
			string_response_HTT+=f"{value},"
			string_target_HTT+=f"{target_value_HTT[i]},"
		string_response_HTT+="\n"
		string_target_HTT+="\n"
		
		get_target_value_LTT = open("response_values_LTT.txt","+a").write(string_response_LTT)
		get_target_value_HTT = open("response_values_HTT.txt","+a").write(string_response_HTT)
		get_opt_LTT = open("Objective_LTT.txt","+a").write(f"{obj_LTT}\n")
		get_target_value_HTT = open("target_values_LTT.txt","+a").write(string_target_LTT)
		get_target_value_HTT = open("target_values_HTT.txt","+a").write(string_target_HTT)
		get_opt_HTT = open("Objective_HTT.txt","+a").write(f"{obj_HTT}\n")
		get_opt = open("Objective.txt","+a").write(f"{obj_LTT+obj_HTT}\n")
		return obj_LTT + obj_HTT
	
	
	def run_multi_stage_parallel_optimization_with_selected_PRS(self,Unsrt_LTT,Unsrt_HTT,ResponseSurfacesLTT,ResponseSurfacesHTT,Input_data):
		"""
		-------------------------
		Multi-Stage Optimization (MSO)
		-------------------------
		- This routine is purely for MSO in Parallel
		- The study is done only for A-factor of reaction rates (ECM-2025, <Journal Paper 2025 citation>
		"""
		self.unsrt_LTC = Unsrt_LTT
		self.ResponseSurfacesLTT = ResponseSurfacesLTT
		self.unsrt_HTC = Unsrt_HTT
		self.ResponseSurfacesHTT = ResponseSurfacesHTT
		self.Input_data = Input_data
		algorithm = Input_data["Type"]["Algorithm"]
		self.rxn_index_LTT = np.asarray([rxn for rxn in Unsrt_LTT]).flatten()
		self.rxn_index_HTT = np.asarray([rxn for rxn in Unsrt_HTT]).flatten()
		
		self.init_guess = np.zeros(len(self.rxn_index_LTT)+len(self.rxn_index_HTT))	
		bounds = tuple([(-1,1) for _ in self.init_guess ])
		
		if Input_data["Stats"]["Design_of_PRS"] == "A-facto":
			opt_output = minimize(self.obj_func_parallel_optimization_with_selected_prs,
					    self.init_guess,
					    bounds=bounds,
					    method='SLSQP',  # Replace 'algorithm' with specific method if known
					    options={"maxiter": 500000}  # Replace 'maxfev' with 'maxiter'
					)
			optimal_parameters = np.asarray(opt_output.x) #Optimized parameter (x*)
			optimal_parameters_zeta = np.asarray(opt_output.x) #For A-factor optimization, zeta vector - [z_a], meaning x* is z_a
			cov = [] #Posterior covariance
		return np.asarray(optimal_parameters),np.asarray(optimal_parameters_zeta),cov
	
