import numpy as np
from solution import Solution
from scipy.optimize import minimize 
from scipy import optimize as spopt
import numpy as np
import time
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib import style
#import pygad
import matplotlib as mpl
mpl.use('Agg')
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
		
		self.count +=1
		note = open("guess_values.txt","+a").write(f"{self.count},{x}\n")
		
		kappa_curve = {}
		count = 0
		for i in self.rxn_index:
			temp = []
			for j in range(len(self.T)):
				temp.append(x[count])
				count += 1
			#trial = [temp[0],(temp[0]+temp[1])/2,temp[1]]
			Kappa = self.kappa_0[i] + temp*(self.kappa_max[i]-self.kappa_0[i])
			kappa_curve[i] = np.asarray(Kappa).flatten()
			
		
		zeta = {}
		for rxn in self.rxn_index:
			zeta[rxn] = self.unsrt[rxn].getZeta_typeA(kappa_curve[rxn])
	
		x_transformed = []
		string = ""
		for rxn in self.rxn_index:
			temp = list(zeta[rxn])
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
		diff = []
		diff_2 = []
		diff_3 = {}
		for i,case in enumerate(self.target_list):
			if self.ResponseSurfaces[i].selection == 1:	
				if case.target == "Tig":
					if case.d_set in frequency:
						frequency[case.d_set] += 1
					else:
						frequency[case.d_set] = 1
					COUNT_All +=1
					COUNT_Tig +=1
					#dataset_weights = (1/len(self.frequency))*float(self.frequency[case.d_set])
					#dataset_weights = (1/COUNT_All)
					
					val = self.ResponseSurfaces[i].evaluate(x)
					#print(val)
					#val,grad = case.evaluateResponse(x)
					f_exp = np.log(case.observed*10)
					#w = 1/(np.log(case.std_dvtn*10))
					w_ = case.std_dvtn/case.observed
					w = 1/w_
					diff.append((val-f_exp)*w)
					#diff.append((val - f_exp)/f_exp)
					diff_2.append(val - f_exp)
					diff_3[case.uniqueID] = (val-f_exp)/f_exp
					response_value.append(val)
					#response_stvd.append(grad)
					#diff.append(val - np.log(case.observed*10))
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
					
					val = np.exp(self.ResponseSurfaces[i].evaluate(x))
					response_value.append(val)
					f_exp = case.observed
					w = 1/(case.std_dvtn)
					diff.append((val - f_exp)*w)
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
					
					val = self.ResponseSurfaces[i].evaluate(x)
					response_value.append(val)
					f_exp = case.observed
					w = 1/(case.std_dvtn)
					diff.append((val - f_exp)*w)
					#response_stvd.append(grad)
					target_value.append(np.log(case.observed))
					target_value_2.append(np.log(case.observed))
					target_stvd.append(1/(np.log(case.std_dvtn)+abs(np.log(case.observed)-val)))
					case_stvd.append(np.log(case.std_dvtn))
					case_systematic_error.append(abs(np.log(case.observed)-val))
					#target_weights.append(dataset_weights)	
			
		
		diff = np.asarray(diff)
		multiplicating_factors = []
		#multiplicating_factors = np.asarray(target_weights)*np.asarray(target_stvd)
		for i,case in enumerate(self.target_list):
			if self.ResponseSurfaces[i].selection == 1:	
				if case.target == "Tig":
					multiplicating_factors.append(1/COUNT_Tig)
		
				elif case.target == "Fls":
					multiplicating_factors.append(1/COUNT_Fls)
		
		multiplicating_factors= np.asarray(multiplicating_factors)		
		#Giving all datapoints equal weights
		#multiplicating_factor = 1/COUNT_All
		
		for i,dif in enumerate(diff):
			#obj+= multiplicating_factors[i]*(dif)**2
			#obj+= multiplicating_factor*(dif)**2
			#obj+= multiplicating_factors[i]*(dif)**2
			obj+=dif**2		
		
		Diff_3 = open("Dataset_based_obj","+a").write(f"{diff_3}\n")
		get_opt = open("Objective.txt","+a").write(f"{obj}\n")
		note = open("guess_values_TRANSFORMED.txt","+a").write(f"{self.count},{x}\n")
		get_target_value = open("response_values.txt","+a").write(f"\t{target_value},{response_value}\n")
		
		
		return obj
	def plot_DATA(self,x):	
		kappa_curve = {}
		count = 0
		for i in self.rxn_index:
			temp = []
			for j in range(len(self.T)):
				temp.append(x[count])
				count += 1
			#trial = [temp[0],(temp[0]+temp[1])/2,temp[1]]
			Kappa = self.kappa_0[i] + temp*(self.kappa_max[i]-self.kappa_0[i])
			kappa_curve[i] = np.asarray(Kappa).flatten()
			
		
		zeta = {}
		for rxn in self.rxn_index:
			zeta[rxn] = self.unsrt[rxn].getZeta_typeA(kappa_curve[rxn])
	
		x_transformed = []
		string = ""
		for rxn in self.rxn_index:
			temp = list(zeta[rxn])
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
		diff = []
		diff_2 = []
		diff_3 = {}
		VALUE = []
		EXP = []
		CASE = []
		TEMPERATURE = []
		for i,case in enumerate(self.target_list):
			if self.ResponseSurfaces[i].selection == 1:	
				if case.target == "Tig":					
					t = case.temperature
					val = self.ResponseSurfaces[i].evaluate(x)
					f_exp = np.log(case.observed*10)
					VALUE.append(np.exp(val)/10)
					EXP.append(np.exp(f_exp)/10)
					CASE.append(case.dataSet_id)
					TEMPERATURE.append(t)				
				elif case.target == "Fls":
					t = case.temperature
					val = np.exp(self.ResponseSurfaces[i].evaluate(x))
					f_exp = case.observed
					VALUE.append(val)
					EXP.append(f_exp)	
					CASE.append(case.dataSet_id)
					TEMPERATURE.append(t)
		
		return VALUE,EXP,CASE,TEMPERATURE
	
	def _obj_function(self,x):
		kappa_curve = {}
		count = 0
		for i in self.rxn_index:
			temp = []
			for j in range(len(self.T)):
				temp.append(x[count])
				count += 1
			Kappa = self.kappa_0[i] + temp*(self.kappa_max[i]-self.kappa_0[i])
			kappa_curve[i] = np.asarray(Kappa).flatten()
			
		
		zeta = {}
		for rxn in self.rxn_index:
			zeta[rxn] = self.unsrt[rxn].getZeta_typeA(kappa_curve[rxn])
	
		x_transformed = []
		string = ""
		for rxn in self.rxn_index:
			temp = list(zeta[rxn])
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
		
		num_params = len(x)
		num_expts = len(self.target_list)
		f = np.empty(num_expts)
		df = np.zeros((num_expts,num_params))
		inv_covar = np.linalg.cholesky(4*np.eye(len(x)))
		initial_guess = np.zeros(len(x))
		#f[0:num_params] = np.dot(inv_covar,(x - initial_guess))
		#df[0:num_params,0:num_params] = inv_covar
		for i,case in enumerate(self.target_list):
			if self.ResponseSurfaces[i].selection == 1:	
				if case.target == "Tig":
					if case.d_set in frequency:
						frequency[case.d_set] += 1
					else:
						frequency[case.d_set] = 1
					COUNT_All +=1
					COUNT_Tig +=1
					#dataset_weights = (1/len(self.frequency))*float(self.frequency[case.d_set])
					#dataset_weights = (1/COUNT_All)
					
					val = self.ResponseSurfaces[i].evaluate(x)
					#response_value.append(val)
					f_exp = np.log(case.observed*10)
					w = 1/(np.log(case.std_dvtn*10))
					f[i] = (val - f_exp)*w
					#print(self.ResponseSurfaces[i].Jacobian(x))
					df[i,:] = np.asarray(self.ResponseSurfaces[i].Jacobian(x))*w
					#raise AssertionError("Optimization about to happen")
					#response_stvd.append(grad)
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
					
					val = np.exp(self.ResponseSurfaces[i].evaluate(x))
					#response_value.append(val)
					f_exp = case.observed
					w = 1/(case.std_dvtn)
					f[i] = (val - f_exp)*w
					df[i,:] = np.asarray(self.ResponseSurfaces[i].Jacobian(x))*w
					#response_stvd.append(grad)
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
					
					val = self.ResponseSurfaces[i].evaluate(x)
					#response_value.append(val)
					f_exp = case.observed
					w = 1/(case.std_dvtn)
					f[i] = (val - f_exp)*w
					df[i,:] = np.asarray(self.ResponseSurfaces[i].Jacobian(x))*w
					#response_stvd.append(grad)
					
					#target_weights.append(dataset_weights)	
			
		return f,df
	def _obj_func(self,x):
		string = ""
		for i in x:
			string+=f"{i},"
		string+=f"\n"
		#x_transformed = np.asarray(x_transformed)
		zeta_file = open("zeta_guess_values.txt","+a").write(string)
		"""
		Just for simplicity
		"""
		
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
		
		
		for i,case in enumerate(self.target_list):
			if self.ResponseSurfaces[i].selection == 1:	
				if case.target == "Tig":
					if case.d_set in frequency:
						frequency[case.d_set] += 1
					else:
						frequency[case.d_set] = 1
					COUNT_All +=1
					COUNT_Tig +=1
					#dataset_weights = (1/len(self.frequency))*float(self.frequency[case.d_set])
					#dataset_weights = (1/COUNT_All)
					
					val = self.ResponseSurfaces[i].evaluate(x)
					#print(val)
					#val,grad = case.evaluateResponse(x)
					response_value.append(val)
					#response_stvd.append(grad)
					target_value.append(np.log(case.observed*10))
					
					target_value_2.append(np.log(case.observed))
					target_stvd.append(1/(np.log(case.std_dvtn*10)))
					case_stvd.append(case.std_dvtn/case.observed)
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
					
					val = np.exp(self.ResponseSurfaces[i].evaluate(x))
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
					
					val = self.ResponseSurfaces[i].evaluate(x)
					response_value.append(val)
					#response_stvd.append(grad)
					target_value.append(np.log(case.observed))
					target_value_2.append(np.log(case.observed))
					target_stvd.append(1/(np.log(case.std_dvtn)+abs(np.log(case.observed)-val)))
					case_stvd.append(np.log(case.std_dvtn))
					case_systematic_error.append(abs(np.log(case.observed)-val))
					#target_weights.append(dataset_weights)	
			
		self.count +=1
		diff = np.asarray(response_value)-np.asarray(target_value)
		multiplicating_factors = []
		#multiplicating_factors = np.asarray(target_weights)*np.asarray(target_stvd)
		for i,case in enumerate(self.target_list):
			if self.ResponseSurfaces[i].selection == 1:	
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
						
		
		note = open("guess_values.txt","+a").write(f"{self.count},{x}\n")
		get_target_value = open("response_values.txt","+a").write(f"\t{target_value},{response_value}\n")
		get_opt = open("Objective.txt","+a").write(f"{obj}\n")
		return obj
	
		
	def run_optimization_with_selected_PRS(self,Unsrt_data,ResponseSurfaces,Input_data):
	   
		self.unsrt = Unsrt_data
		self.ResponseSurfaces = ResponseSurfaces
		self.Input_data = Input_data
		algorithm = Input_data["Type"]["Algorithm"]
		self.rxn_index = np.asarray([rxn for rxn in Unsrt_data]).flatten()
		self.init_guess = np.zeros(len(self.rxn_index))	
		bounds = tuple([(-1,1) for _ in self.init_guess ])
		
		if Input_data["Stats"]["Design_of_PRS"] == "A-facto":
		
			opt_output = minimize(self._obj_func,self.init_guess,bounds=bounds,method=algorithm)
			print(opt_output)
			optimal_parameters = np.asarray(opt_output.x)
			optimal_parameters_zeta = np.asarray(opt_output.x)
			cov = []
		
		else:
			self.rxn_index = []
			self.kappa_0 = {}
			self.kappa_max = {}
			
			for rxn in self.unsrt:
				self.rxn_index.append(rxn)
			self.T = np.linspace(300,2500,3)
			#self.init_guess = np.zeros(len(self.T[1:])*len(self.rxn_index))
			self.init_guess = np.zeros(len(self.T)*len(self.rxn_index))
			bounds = tuple([(-1,1) for _ in self.init_guess ])
			
			theta = np.array([self.T/self.T,np.log(self.T),-1/self.T])
			#self.theta_inv = np.linalg.inv(theta.T)
			
			for rxn in self.rxn_index:
				#self.activeParameters[rxn] = self.unsrt[rxn].activeParameters
				self.kappa_0[rxn] = self.unsrt[rxn].getNominal(self.T)
				self.kappa_max[rxn] = self.unsrt[rxn].getKappaMax(self.T)
			
			start = time.time()
			opt_output = minimize(self.obj_func_of_selected_PRS,self.init_guess,bounds=bounds,method='Powell',options={"maxfev":500000})
			stop = time.time()
			final=print(f"Time taken for optimization {stop-start}")
			#opt_output = spopt.root(self._obj_function,self.init_guess,method="lm",jac = True)
			print(opt_output)
			optimal_parameters = np.asarray(opt_output.x)
			#residuals,final_jac= self._obj_function(optimal_parameters)
			#icov = np.dot(final_jac.T,final_jac)
			#cov = np.linalg.inv(icov)
			#raise AssertionError("Optimization done!")
			"""
			Finding the optimal parameters
			"""
			print("<<<<<<<<<<<<<<<<FOUND BEST SOLUTION>>>>>>>>>>>>>>>>>>>>>\n")
			#print(f"{solution}")
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
			kappa_curve = {}
			count = 0
			
			for rxn in self.rxn_index:
				temp = []
				for j in range(len(self.T)):
					temp.append(optimal_parameters[count])
					count += 1
				#Kappa = np.exp(np.log(kappa_0[i]) + temp*(np.log(kappa_max[i])-np.log(kappa_0[i])))
				#For searching class-C type curves
				#trial = [temp[0],(temp[0]+temp[1])/2,temp[1]]
				Kappa = self.kappa_0[rxn] + temp*(self.kappa_max[rxn]-self.kappa_0[rxn])
				
				#print(Kappa)
				kappa_curve[rxn] = np.asarray(Kappa).flatten()
				
	
			zeta = {}
			#print(gen)
			for rxn in self.rxn_index:
				zeta[rxn] = self.unsrt[rxn].getZeta_typeA(kappa_curve[rxn])
			optimal_parameters_zeta = []
			for rxn in self.rxn_index:
				temp = list(zeta[rxn])
				optimal_parameters_zeta.extend(temp)
			#########################
			### For plotting purposes
			########################
			delta_n = {}
			p = {}
			V_opt = {}
			V = {}#Populating V for unshuffled portion of the design matrix
			ch = {}
			nominal = {}
			p_max = {}
			p_min = {}
			theta = {}
			Temp = {}
			d_n = {}
			for index,rxn in enumerate(self.unsrt):
				ch[rxn] = self.unsrt[rxn].cholskyDeCorrelateMat
				nominal[rxn] = self.unsrt[rxn].nominal
				p = self.unsrt[rxn].nominal
				p_max[rxn] = self.unsrt[rxn].P_max
				p_min[rxn] = self.unsrt[rxn].P_min
				theta[rxn] = self.unsrt[rxn].Theta
				Temp[rxn] = self.unsrt[rxn].temperatures
				Tp = self.unsrt[rxn].temperatures
				#Tp = np.linspace(300,2500,50)
				Theta_p = np.array([Tp/Tp,np.log(Tp),-1/(Tp)])
				#P_max = P + np.asarray(np.dot(cov,zet)).flatten();
				#P_min = P - np.asarray(np.dot(cov,zet)).flatten();
				kmax = Theta_p.T.dot(p_max[rxn])
				kmin = Theta_p.T.dot(p_min[rxn])
				ka_o = Theta_p.T.dot(nominal[rxn])
				p_zet = p + np.asarray(np.dot(ch[rxn],zeta[rxn])).flatten();
				k =  Theta_p.T.dot(p_zet)
				fig = plt.figure()
				plt.title(str(rxn))
				plt.xlabel(r"1000/T\K$^{-1}$")
				plt.ylabel(r"$log_{10}(k)$ / $s^{-1}$ or $log_{10}$(k) / $cm^{3}\,molecule^{-1}\,s^{-1}$")
				plt.plot(1/Tp,kmax,'k--',label="Uncertainty limits")
				plt.plot(1/Tp,kmin,'k--')
				plt.plot(1/Tp,ka_o,'b-',label='Prior rate constant')
				plt.plot(1/Tp,k,'r-',label='Optimized rate constant')
				d_n[rxn]=abs(p[1]-p_zet[1])
				plt.savefig("Plots/reaction_"+str(rxn)+".png",bbox_inches='tight')
			optimal_parameters_zeta = np.asarray(optimal_parameters_zeta)	
			cov = []
			print(d_n)
			
			###############
			# Printing the optimal parameters
			###############
			#Optimal parameters
			#cases
			#response surface
			VALUE,EXP,CASE,TEMPERATURE = self.plot_DATA(optimal_parameters)
			DATASET = set(CASE)
			
			for case_id in DATASET:
				file_ =open(str(case_id)+".csv","+w")
				STRING = "T(k)\tf_exp\tValue\n"
				for i,case in enumerate(CASE):
					if case ==case_id:
						STRING+=f"{TEMPERATURE[i]}\t{EXP[i]}\t{VALUE[i]}\n" 
				file_.write(STRING)
		
						
		return np.asarray(optimal_parameters),np.asarray(optimal_parameters_zeta),cov
	
