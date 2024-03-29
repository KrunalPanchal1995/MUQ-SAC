import os, shutil, re, math #default python modules
import numpy as np
import pandas as pd
from pyDOE2 import *
home_dir = os.getcwd()
import multiprocessing
import subprocess
import time
import sys
import concurrent.futures
import asyncio
import Uncertainty
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib import style
import matplotlib
#matplotlib.use('Agg')

from sklearn.datasets import make_regression
import os

from keras.models import Sequential
import tensorflow as tf
#keras = tf.compat.v1.keras
#Sequential = keras.utils.Sequence
from keras.layers import Dense
#function to create directories and the required input files in a systematic way directories are named with the reaction index
import Make_input_file

def run_sampling(sample,data,generator,length):
	A = Uncertainty.UncertaintyExtractor(data)
	a1 = generator[0]
	a2 = generator[1]
	a3 = generator[2]
	A.populateValues(a1,a2,a3)
	A.getCovariance(flag=False)
	A.getUnCorrelated(flag = False)
	zeta = A.getB2Zeta(flag=True)
	del A
	return (sample,generator,zeta,length)
	
def run_sampling_b(sample,data,generator,length):
	A = Uncertainty.UncertaintyExtractor(data)
	a1 = generator[0]
	a2 = generator[1]
	A.populateValues(a1,a2)
	A.getCovariance(flag=False)
	A.getUnCorrelated(flag = False)
	zeta = A.getB2Zeta(flag=True)
	del A
	return (sample,generator,zeta,length)
	
def run_sampling_c(sample,data,generator,length):
	A = Uncertainty.UncertaintyExtractor(data)
	a1 = generator[0]
	a2 = generator[1]
	A.populateValues(a1,a2)
	A.getCovariance(flag=False)
	A.getUnCorrelated(flag = False)
	zeta = A.getC2Zeta(flag=True)
	del A
	return (sample,generator,zeta,length)
	
class Worker():
	def __init__(self, workers):
		self.pool = multiprocessing.Pool(processes=workers)
		self.progress = []
		self.parallized_zeta = []
		self.generator = []
		self.parallel_zeta_dict = {}

	def callback(self,result):
		#print("Entered callback\n")
		self.progress.append(result[0])
		self.generator.append(result[1])
		self.parallized_zeta.append(result[2])
		sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(self.progress)/float(result[-1])*100))
		sys.stdout.flush()
	def callback_error(self,result):
    		print('error', result)
    		
	def do_unsrt(self,data,sampling_points):
		#generators = data["generators"]
		for args in range(sampling_points):
			self.pool.apply_async(run_sampling, 
				  args=(1,data,data["generators"][args],sampling_points,), 
				  callback=self.callback)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
		return self.generator,self.parallized_zeta
	def do_unsrt_b(self,data,sampling_points):
		#generators = data["generators"]
		for args in range(sampling_points):
			self.pool.apply_async(run_sampling_b, 
				  args=(1,data,data["generators_b"][args],sampling_points,), 
				  callback=self.callback,error_callback=self.callback_error)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
		return self.generator,self.parallized_zeta
	def do_unsrt_c(self,data,sampling_points):
		#generators = data["generators"]
		for args in range(sampling_points):
			self.pool.apply_async(run_sampling_c, 
				  args=(1,data,data["generators_c"][args],sampling_points,), 
				  callback=self.callback,error_callback=self.callback_error)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
		return self.generator,self.parallized_zeta

class Model(object):
	def __init__(self,x,y):
		self.xData = x[0:1,:]
		self.yData = y[0:3,:]
		#print(x[:1],y[:3])
		#print(type(x),type(y))
	def get_model(self,n_inputs, n_outputs):
		inputs = tf.keras.Input(shape=(n_inputs,))
		x = tf.keras.layers.Dense(10)(inputs)
		outputs = tf.keras.layers.Dense(n_outputs)(x)
		model = tf.keras.Model(inputs, outputs)
		# Activity regularization.
		model.add_loss(tf.abs(tf.reduce_mean(x)))

		#model = Sequential()
		
		"""
	    First layer with 20 neurons connected with input layers
		"""
		
		#model.add(Dense(32, input_dim=n_inputs, activation='relu'))
		#model.add(Dense(16, activation='relu'))
		"""
		Output layer
		"""
		
		#model.add(Dense(n_outputs))
		model.compile(loss='mae', optimizer='adam')
		return model
	
	def Train(self):
		n_inputs, n_outputs =  self.xData.shape[1],self.yData.shape[1]
		print(f"\n{n_inputs},{n_outputs}\n")
		
		#n_inputs, n_outputs = 25,75
		self.model = self.get_model(n_inputs, n_outputs)
		self.model.fit(self.xData, self.yData, verbose=0, epochs=100)
		return self.model
				
class Sample :
	def __init__(self,rxnUnsrt,reactionList,sample_length):
		self.rxnUnsert = rxnUnsrt
		self.ind = reactionList
		self.sim_ = sample_length
	
	def get_SAMAP(self,samap_executable,jpdap_executable):
		"""
		Getting the SAMAP amples
		"""
		os.mkdir("SAMAP")
		os.chdir("SAMAP")
		input_rxn_dict = {}
		for i in self.ind:
			input_dict = {}
			input_dict["samples"] = self.sim_
			input_dict["samples_skipped"] = int(0.1*self.sim_)
			input_dict["Random_seed"] = 1
			input_dict["sampling_method"] = "SOBOL"
			input_dict["sampling_distribution"] = "UNIFORM"
			input_dict["equidistant_T"] = 100
			input_dict["T_begin"] = data = self.rxnUnsert[i].temperatures[0]
			input_dict["T_end"] = self.rxnUnsert[i].temperatures[-1]
			input_dict["L"] = 0
			input_dict["len_temp_data"] = len(self.rxnUnsert[i].temperatures)
			string_unsrt_data =""
			for index,k in enumerate(self.rxnUnsert[i].temperatures):
				string_unsrt_data+=f"{k} {self.rxnUnsert[i].uncertainties[index]} \n"
			input_dict["temperature_unsrt_data"] = string_unsrt_data
			input_dict["alpha"] = np.exp(self.rxnUnsert[i].rxn_dict[0])
			input_dict["n"] = self.rxnUnsert[i].rxn_dict[1]
			input_dict["n_min"] = self.rxnUnsert[i].rxn_dict[1]-2
			input_dict["n_max"] = self.rxnUnsert[i].rxn_dict[1]+2
			input_dict["epsilon"] = self.rxnUnsert[i].rxn_dict[2]
			L = self.rxnUnsert[i].cholskyDeCorrelateMat
			#input_dict["covariance_matrix"] = str(L.dot(L.T)).strip("[]").replace("[","").replace("]","")
			input_dict["uncertainty_type"] = "2slnk"
			if len(self.rxnUnsert[i].activeParameters) == 3:
				input_dict["uncertain_parameters"] = "AnE"
			else:
				input_dict["uncertain_parameters"] = "A"
			input_rxn_dict[i] = input_dict
			string_dict = {}
			jpdap_instring = Make_input_file.create_JPDAP_input(input_dict)
			"""
			Run: JPDAP code
			"""
			
			os.mkdir(f"{self.rxnUnsert[i].rIndex}")
			os.chdir(f"{self.rxnUnsert[i].rIndex}")
			file_jpdap = open("jpdap_data_"+str(self.rxnUnsert[i].rIndex)+".txt","w").write(jpdap_instring)
			run_jpdap_string = f"""#!/bin/bash
{jpdap_executable} jpdap_data_{self.rxnUnsert[i].rIndex}.txt &> out"""
			file_print_run_jpdap = open("run_jpdap","w").write(run_jpdap_string)
			subprocess.call(["chmod","+x",'run_jpdap'])
			start_Jpdap = time.time()
			subprocess.call(["./run_jpdap"])
			stop_Jpdap = time.time()
			print(f"\n\tJPDAP code took {stop_Jpdap-start_Jpdap}s to execute\n")
			Nagy_covariance_matrix = ""
			file_name_jpdap = "jpdap_data_"+str(self.rxnUnsert[i].rIndex)+".txt_fit_minRMSD.txt"
			Nagy_covariance_matrix = open(file_name_jpdap,"r").readlines()[5:8]
			covariance_matrix = ""
			for n in Nagy_covariance_matrix:
				covariance_matrix+=str(n)
			#print(covariance_matrix)
			input_dict["covariance_matrix"] = covariance_matrix.strip("\n")
			#input_dict["covariance_matrix"] = str(L.dot(L.T)).strip("[]").replace("[","").replace("]","")
			"""
			Run: SAMAP code
			"""
			instring = Make_input_file.create_SAMAP_input(input_dict)
			string_dict[i] = instring
			file_print = open("samap_data_"+str(self.rxnUnsert[i].rIndex)+".txt","w").write(instring)
			run_string = f"""#!/bin/bash
{samap_executable} samap_data_{self.rxnUnsert[i].rIndex}.txt &> out"""
			file_print_run = open("run","w").write(run_string)						
			subprocess.call(["chmod","+x",'run'])
			subprocess.call(["./run"])
			file_name = "samap_data_"+str(self.rxnUnsert[i].rIndex)+".txt_Arrhpar.txt"
			sampling_file = open(file_name,"r").readlines()
			self.Nagy_arrhenius_samples = [] 
			for sample in sampling_file[1:]:
				self.Nagy_arrhenius_samples.append(np.asfarray(sample.strip("''").strip("\n").split()[:3],float))
			
			#print(self.Nagy_arrhenius_samples)
			Y = self.Nagy_arrhenius_samples
			data = self.rxnUnsert[i].data
			
			P = self.rxnUnsert[i].getMean()	
			print(P)
			cov = self.rxnUnsert[i].getCov()
			T_ = np.linspace(300,3000,100)
			Theta = np.array([T_/T_,np.log(T_),-1/T_])
			fig = plt.figure()
			plt.plot(1/T_,self.rxnUnsert[i].getKappaMax(T_),'k-',linewidth=0.75)
			plt.plot(1/T_,self.rxnUnsert[i].getKappaMin(T_),'k-',linewidth=0.75)
			plt.plot(1/T_,self.rxnUnsert[i].getNominal(T_),"b-",linewidth=0.75)
			for j in Y:
				#print(j)
				#print(kappa_max,np.asarray(theta.T.dot(j)).flatten())
				#print(np.asarray(theta.T.dot(j)).flatten()/kappa_max)
				#Pint = P + np.asarray(np.dot(cov,np.asarray(j))).flatten();
				Pint = np.asarray(j)
				#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
				#X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
				plt.plot(1/T_,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.10)
			
			plt.show()
			
			os.chdir("..")
		os.chdir("..")
	
	def getSamples(self):
		self.beta_=[]
		temp_a = {}
		temp_b = {}
		temp_c = {}
		gen = {}
		gen_a = {}
		gen_b = {}
		gen_c = {}
		
		xdata = {}
		ydata = {}				
		self.generator_list = []
		"""
		Min- simulation = 300
		"""
		#if self.sim_ < 300:
		#	self.sim_ = 300
		#print(self.sim_)
		"""
		Dividing the generators into three bins
		------------------------------------
		n_a = 0.33*self.sim_
		n_b = 0.33*self.sim_
		n_c = 0.34*self.sim_
		
		"""
		#n_a = int(0.1*self.sim_)
		#n_b = int(0.45*self.sim_)
		#n_c = self.sim_-n_a-n_b
		n_a = 50
		n_b = 50
		n_c = 50
		
		#n_b = int(0.5*self.sim_)
		#n_c = self.sim_-n_b
		
		
		#for i in np.eye(self.n):
		#	self.generator_list.append(i)
		#	self.generator_list.append(-1*i)
						
		"""
		Mixing three types of the extreme curves
		--------------------
		zeta-A: 1
		zeta-B: 1
		zeta-C: 1
	
		for zeta-A: one normal random parameter is required
			a matrix on n_a x n will be created
			
		for zeta-B and zeta-c: two random parameters is required
			we can arrange the two normal_random_parameter 				
			a matrix of n_b+n_c x 2*n will be created
			
		The number of reactions = n_R
		"""
		n_rxn = len(self.ind)
		self.generator_list_A = []
		self.generator_list_B = []
		self.generator_list_C = []
		
		"""
		Factorial Design
		----------------
		box-benkhen design
		bb = []
		for i in bbdesign(n_rxn)[0:int(0.3*self.sim_)]:
			bb.append(np.asarray(i))
		"""
		
				
				#For zeta-A
				
		"""
		bb_a = []
		for i in bbdesign(n_rxn):
			bb_a.append(np.asarray(i))
		self.generator_list_A.extend(bb_a)
				
				#For zeta-B
				
		bb_b = []
		for i in bbdesign(2*n_rxn):
			bb_b.append(np.asarray(i))							
		self.generator_list_B.extend(bb_b)
				
				#For zeta-C
				
		bb_c = []
		for i in bbdesign(2*n_rxn):
			bb_c.append(np.asarray(i))		
				
		self.generator_list_C.extend(bb_c)		
				
		"""
		"""
		Monte-Carlo
		
		"""
				
				#For zeta-A
		
		self.generator_list_A.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(n_a)])
		
		self.generator_list_A.extend(list(np.eye(n_rxn)))
		self.generator_list_A.extend(list(-1*np.eye(n_rxn)))
		#self.generator_list_A.extend(list(1*np.ones(n_rxn)))
		#self.generator_list_A.append(np.zeros(n_rxn))
		#self.generator_list_A.append(np.zeros(n_rxn))
		#self.generator_list_A.append(np.zeros(n_rxn))
		
				#For zeta-A
				
		self.generator_list_B.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(n_b)])	
		self.generator_list_B.extend(list(np.eye(2*n_rxn)))
		self.generator_list_B.extend(list(-1*np.eye(2*n_rxn)))
		#self.generator_list_B.append(np.zeros(2*n_rxn))
		#self.generator_list_B.append(np.zeros(2*n_rxn))
		#self.generator_list_B.append(np.zeros(2*n_rxn))
				
				#For zeta-A		
		self.generator_list_C.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(n_c)])
		self.generator_list_C.extend(list(np.eye(2*n_rxn)))
		self.generator_list_C.extend(list(-1*np.eye(2*n_rxn)))
		#self.generator_list_C.append(np.zeros(2*n_rxn))
		#self.generator_list_C.append(np.zeros(2*n_rxn))
		#self.generator_list_C.append(np.zeros(2*n_rxn))
		"""
		Uniform
		
		"""
			
		#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
		
		"""
		Latin-Hypercube
		"""
	
				
				#For zeta-A
		"""		
		lhs_a = []
		for i in lhs(n_rxn, samples=n_a):
			lhs_a.append(np.asarray(i))	
		self.generator_list_A.extend(lhs_a)
		
				
				#For zeta-A
				
		lhs_b = []
		for i in lhs(2*n_rxn, samples=n_b):
			lhs_b.append(np.asarray(i))		
		self.generator_list_B.extend(lhs_b)
				
				#For zeta-A
				
		lhs_c = []
		for i in lhs(2*n_rxn, samples=n_c):
			lhs_c.append(np.asarray(i))
		self.generator_list_C.extend(lhs_c)
		
		"""
		"""
		Monte-Carlo-Mix
		"""
		
		
		
		
		
		#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)		
		#self.generator_list.extend(bb)
		#self.generator_list.extend(list(np.eye(self.n)))
		#self.generator_list.extend(list(-1*np.eye(self.n)))
		#self.generator_list.append(np.zeros(self.n))
		#self.generator_list.append(np.zeros(self.n))
		#self.generator_list.append(np.zeros(self.n))
		#print(self.generator_list)
		#self.sim_ = len(self.generator_list)
		
		#trainers = np.asarray(trainers)
		#print(f"Generator_list {np.shape(self.generator_list)}\n")
		df = pd.DataFrame(np.asarray(self.generator_list))
		count_start_a = 0
		count_start_b = 0
		count_start_c = 0			
		for i in self.ind:
			
			len_of_Arrhenius_A = 1
			count_end_a = count_start_a+len_of_Arrhenius_A
			
			
			len_of_Arrhenius_B = 2
			count_end_b = count_start_b+len_of_Arrhenius_B
			
			
			len_of_Arrhenius_C = 2
			count_end_c = count_start_c+len_of_Arrhenius_C
			
			#print(count_start)
			#print(count_end)
			self.rxn_generators_a2 = []
			self.rxn_generators_b2 = []
			self.rxn_generators_c2 = []
			
			"""
			getting zetas from generators
			"""
			#array_a = [[1.0],[-1.0]]
			#array_b = [[1.0,-1.0],[-1.0,1.0]]
			#array_c = [[1.0,-1.0],[-1.0,1.0]]
			for j in range(len(self.generator_list_A)):
				self.rxn_generators_a2.append(self.generator_list_A[j][count_start_a:count_end_a])
			
			"""
			for j in range(len(array_a)):
				self.rxn_generators_a2.append(array_a[j])
				
			for j in range(len(array_b)):
				self.rxn_generators_b2.append(array_b[j])
				
			for j in range(len(array_c)):
				self.rxn_generators_c2.append(array_c[j])
				
			"""
			#print(self.rxn_generators_a2)
			for j in range(len(self.generator_list_B)):
			
				self.rxn_generators_b2.append(self.generator_list_B[j][count_start_b:count_end_b])
			
			for j in range(len(self.generator_list_C)):
				self.rxn_generators_c2.append(self.generator_list_C[j][count_start_c:count_end_c])
			
			self.sim_ = len(self.rxn_generators_a2)+len(self.rxn_generators_b2)+len(self.rxn_generators_c2)
			
			#self.sim_ = len(array_a)+len(array_b)+len(array_c)
			#print(self.rxn_generators)
			#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
								
			g_a = self.rxn_generators_a2
			g_b = self.rxn_generators_b2
			g_c = self.rxn_generators_c2
			#g_b = self.rxn_generators[0:n_b]
			#g_c = self.rxn_generators[n_b:n_b+n_c]
			
			
			count_start_a = count_end_a
			count_start_b = count_end_b
			count_start_c = count_end_c
			
			"""
			Getting the uncertainty data
			"""
			#fig = plt.figure(figsize=(6,8))
			#gs = fig.add_gridspec(3, 1, hspace=0.19, wspace=0.2)
			
			#(ax1),(plt),(ax3) = gs.subplots(sharex='row')
			fig = plt.figure()
						
			data = self.rxnUnsert[i].data
			Tr = self.rxnUnsert[i].temperatures
			T = np.linspace(Tr[0],Tr[-1],100)
			kappa_max = self.rxnUnsert[i].getKappaMax(T)			
			kappa_min = self.rxnUnsert[i].getKappaMin(T)
						
			Theta = np.array([T/T,np.log(T),-1/T])
			P = self.rxnUnsert[i].getMean()	
			#print(P)
			kappa_0 = Theta.T.dot(P)
			cov = self.rxnUnsert[i].getCov()
			zeta_A = self.rxnUnsert[i].zeta.x
			
			#print(os.getcwd())
			"""
			Type-A zeta samples
			"""
			generators_a = g_a
			parallel_zetas_a = []
			
			time_start = time.time()
			
			
			for k in g_a:
				Pint = P+np.asarray(np.dot(cov,k[0]*zeta_A)).flatten()
				plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"k-",linewidth=1.25)
			#	plt.plot()
				#print(i)
				#print(type(i))
				#print(k)	
				#parallel_zetas_a.append(P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten())
				parallel_zetas_a.append(k[0]*zeta_A)
								
			#print(parallel_zetas_a)
			plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"k-",linewidth=1.25,label=r"A-type samples ($\beta\zeta_{A}$)")
			plt.plot(1/T,kappa_max,"k--",linewidth=1.25,label=r"Rate uncertainty limits ($\kappa(\pm\zeta_{A}$))")
			plt.plot(1/T,kappa_min,"k--",linewidth=1.25)
			plt.plot(1/T,kappa_0,"b-",linewidth=1.25,label="Nominal rate value")
			#plt.xlim([1/2500,1/900])			
			#plt.ylim([23.3,29.67])				
			#plt.xlabel(r"1000/T \ $K^{-1}$")
			#plt.ylabel(r"$log_{10}(k)\,/\,s^{-1}\,or\,log_{10}(k)\, /\, cm^{3} \, mol^{-1} \, s^{-1}$")
			#plt.legend()
			#plt.savefig("./Special_cases/"+str(i)+"_a.pdf")
			
			#print(os.getcwd())
			#print("plot saved")
			data["generators_b"] = g_b
			callWorkForce = Worker(10)	
			generators_b,parallel_zetas_b = callWorkForce.do_unsrt_b(data,len(g_b))
			del callWorkForce
			
			data["generators_c"] = g_c
			callWorkForce = Worker(10)	
			generators_c,parallel_zetas_c = callWorkForce.do_unsrt_c(data,len(g_c))
			del callWorkForce
			
			#print(parallel_zetas_a)
			
			temp_a[i] = parallel_zetas_a
			gen_a[i] = generators_a
			temp_b[i] = parallel_zetas_b
			gen_b[i] = generators_b
			temp_c[i] = parallel_zetas_c
			gen_c[i] = generators_c
			
			time_end = time.time()
			print(f"\n\t\tTime taken = {time_end-time_start}\n")
			
			#print(parallel_zetas_b)
			Y_a = parallel_zetas_a
			Y_b = parallel_zetas_b
			Y_c = parallel_zetas_c
			T_ = np.linspace(300,3000,8)
			theta = np.array([T/T,np.log(T),-1/T])
			kappa_mean = theta.T.dot(P)
			Theta = np.array([T/T,np.log(T),-1/T])
			X = []
			Y = []
			
			
			for j in Y_a:
				Y.append(j)
				Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
				
				X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
				#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
			#fig2 = plt.figure()
			
			#fig = plt.figure()
			for j in Y_b:
				Y.append(j)
				Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
				plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"ro-",linewidth=1.25)
				X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
				#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
			plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"ro-",linewidth=1.25,label=r"B-type samples ($\zeta_{B}$)")
			#plt.plot(1/T,kappa_max,"k--",linewidth=1.25,label=r"Rate uncertainty limits ($\kappa(\pm\zeta_{A}$))")
			#plt.plot(1/T,kappa_min,"k--",linewidth=1.25)
			#plt.plot(1/T,kappa_0,"b-",linewidth=1.25,label="Nominal rate value")
			#plt.xlim([1/2500,1/900])
			#plt.ylim([23.3,29.67])				
			#plt.xlabel(r"1000/T \ $K^{-1}$")
			#plt.ylabel(r"$log_{10}(k)\,/\,s^{-1}\,or\,log_{10}(k)\, /\, cm^{3} \, mol^{-1} \, s^{-1}$")
			#plt.legend()
			#plt.savefig("./Special_cases/"+str(i)+"_b.pdf")
			
			#fig = plt.figure()
			#fig3 = plt.figure()
			for j in Y_c:
				Y.append(j)
				Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
				plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"c--",linewidth=1.25)
				X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
				#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
			plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"c--",linewidth=1.25,label=r"C-type samples ($\zeta_{C}$)")
			#plt.plot(1/T,kappa_max,"k--",linewidth=1.25,label=r"Rate uncertainty limits ($\kappa(\pm\zeta_{A}$))")
			#plt.plot(1/T,kappa_min,"k--",linewidth=1.25)
			#plt.plot(1/T,kappa_0,"b-",linewidth=1.25,label="Nominal rate value")
			#plt.xlim([1/2500,1/900])
			#plt.ylim([23.3,29.67])			
			plt.legend()
			#fig.supxlabel(r"1000/T \ $K^{-1}$", x = 0.16, y = -0.17, va = "bottom", ha = "left")
			#fig.supylabel(r"$log_{10}(k)\,/\,s^{-1}\,or\,log_{10}(k)\, /\, cm^{3} \, mol^{-1} \, s^{-1}$",x = 0.05, y = 0.53)
			#plt.show()
			plt.xlabel(r"1000/T \ $K^{-1}$")
			plt.ylabel(r"$log_{10}(k)\,/\,s^{-1}\,or\,log_{10}(k)\, /\, cm^{3} \, mol^{-1} \, s^{-1}$")
			plt.savefig("./Special_cases/"+str(i)+"_c.pdf")
			xdata[i] = X
			ydata[i] = Y
			#print(len(parallel_zetas))
			#print(len(self.rxn_generators))
		A = []
		A_gen = []
		X_dash = []
		y_dash = []
		
		for j in range(len(g_a)):
			a_row = []
			a_row_gen = []
			x_temp = []
			y_temp = []
			for i in temp_a:
				a_row.extend(list(temp_a[i][j]))
				a_row_gen.extend(list(gen_a[i][j]))
				x_temp.extend(list(xdata[i][j]))
				y_temp.extend(list(ydata[i][j]))
			A.append(np.asarray(a_row))
			A_gen.append(np.asarray(a_row_gen))
			X_dash.append(np.asarray(x_temp))
			y_dash.append(np.asarray(y_temp))
		
		
		for j in range(len(g_b)):
			a_row = []
			a_row_gen = []
			x_temp = []
			y_temp = []
			for i in temp_b:
				a_row.extend(list(temp_b[i][j]))
				a_row_gen.extend(list(gen_b[i][j]))
				x_temp.extend(list(xdata[i][j+len(g_a)]))
				y_temp.extend(list(ydata[i][j+len(g_a)]))
			A.append(np.asarray(a_row))
			A_gen.append(np.asarray(a_row_gen))
			X_dash.append(np.asarray(x_temp))
			y_dash.append(np.asarray(y_temp))
		"""
		for j in range(n_b):
			a_row = []
			a_row_gen = []
			x_temp = []
			y_temp = []
			for i in temp_b:
				a_row.extend(list(temp_b[i][j]))
				a_row_gen.extend(list(gen_b[i][j]))
				x_temp.extend(list(xdata[i][j]))
				y_temp.extend(list(ydata[i][j]))
			A.append(np.asarray(a_row))
			A_gen.append(np.asarray(a_row_gen))
			X_dash.append(np.asarray(x_temp))
			y_dash.append(np.asarray(y_temp))
		"""
		
		for j in range(len(g_c)):
			a_row = []
			a_row_gen = []
			x_temp = []
			y_temp = []
			for i in temp_c:
				a_row.extend(list(temp_c[i][j]))
				a_row_gen.extend(list(gen_c[i][j]))
				x_temp.extend(list(xdata[i][j+len(g_a)+len(g_b)]))
				y_temp.extend(list(ydata[i][j+len(g_a)+len(g_b)]))
			A.append(np.asarray(a_row))
			A_gen.append(np.asarray(a_row_gen))
			X_dash.append(np.asarray(x_temp))
			y_dash.append(np.asarray(y_temp))
		"""
		for j in range(n_c):
			a_row = []
			a_row_gen = []
			x_temp = []
			y_temp = []
			for i in temp_c:
				a_row.extend(list(temp_c[i][j]))
				a_row_gen.extend(list(gen_c[i][j]))
				x_temp.extend(list(xdata[i][j+n_b]))
				y_temp.extend(list(ydata[i][j+n_b]))
			A.append(np.asarray(a_row))
			A_gen.append(np.asarray(a_row_gen))
			X_dash.append(np.asarray(x_temp))
			y_dash.append(np.asarray(y_temp))
			
		"""
		
		#print(A)
		#self.sim_ = n_a+n_b+n_c
		#self.sim_ = n_b+n_c
		self.beta_ = np.asarray(A)
		#print(self.beta_)
		self.generators = np.asarray(A_gen)
		self.X = np.asarray(X)
		self.y = np.asarray(y_dash)
		"""
		Training the model for reduction of parameters
		"""
		
		#A = Model(self.X,self.y)
		#model = A.Train()
		#trained_model[i] = model
		#for i,ele in enumerate(self.y):
		#	print(ele,model.predict(self.X[i]))
		
		"""
		
		"""
		
		
		print(self.beta_)
		print(self.X)

	def plotSamples(self,zeta_file):
		self.beta_=[]
		temp_a = {}
		temp_b = {}
		temp_c = {}
		gen = {}
		gen_a = {}
		gen_b = {}
		gen_c = {}
		
		xdata = {}
		ydata = {}				
		self.generator_list = []
		"""
		Min- simulation = 300
		"""
		#if self.sim_ < 300:
		#	self.sim_ = 300
		#print(self.sim_)
		"""
		Dividing the generators into three bins
		------------------------------------
		n_a = 0.33*self.sim_
		n_b = 0.33*self.sim_
		n_c = 0.34*self.sim_
		
		"""
		#n_a = int(0.1*self.sim_)
		#n_b = int(0.45*self.sim_)
		#n_c = self.sim_-n_a-n_b
		n_a = 50
		n_b = 50
		n_c = 50
		
		#n_b = int(0.5*self.sim_)
		#n_c = self.sim_-n_b
		
		
		#for i in np.eye(self.n):
		#	self.generator_list.append(i)
		#	self.generator_list.append(-1*i)
						
		"""
		Mixing three types of the extreme curves
		--------------------
		zeta-A: 1
		zeta-B: 1
		zeta-C: 1
	
		for zeta-A: one normal random parameter is required
			a matrix on n_a x n will be created
			
		for zeta-B and zeta-c: two random parameters is required
			we can arrange the two normal_random_parameter 				
			a matrix of n_b+n_c x 2*n will be created
			
		The number of reactions = n_R
		"""
		n_rxn = len(self.ind)
		self.generator_list_A = []
		self.generator_list_B = []
		self.generator_list_C = []
		
		generator_file = open(zeta_file,"r").readlines() 
		self.generator_list_A = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in generator_file])	
				#For zeta-A
		
		
		df = pd.DataFrame(np.asarray(self.generator_list))
		count_start_a = 0		
		for i in self.ind:
			
			len_of_Arrhenius_A = 3
			count_end_a = count_start_a+len_of_Arrhenius_A
			
			
			#print(count_start)
			#print(count_end)
			self.rxn_generators_a2 = []
			self.rxn_generators_b2 = []
			self.rxn_generators_c2 = []
			
			"""
			getting zetas from generators
			"""
			#array_a = [[1.0],[-1.0]]
			#array_b = [[1.0,-1.0],[-1.0,1.0]]
			#array_c = [[1.0,-1.0],[-1.0,1.0]]
			for j in range(len(self.generator_list_A)):
				self.rxn_generators_a2.append(self.generator_list_A[j][count_start_a:count_end_a])
					
			
			
			#self.sim_ = len(array_a)+len(array_b)+len(array_c)
			#print(self.rxn_generators)
			#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
								
			g_a = self.rxn_generators_a2
			
			#g_b = self.rxn_generators[0:n_b]
			#g_c = self.rxn_generators[n_b:n_b+n_c]
			
			
			count_start_a = count_end_a
			
			"""
			Getting the uncertainty data
			"""
			#fig = plt.figure(figsize=(6,8))
			#gs = fig.add_gridspec(3, 1, hspace=0.19, wspace=0.2)
			
			#(ax1),(plt),(ax3) = gs.subplots(sharex='row')
			fig = plt.figure()
						
			data = self.rxnUnsert[i].data
			Tr = self.rxnUnsert[i].temperatures
			T = np.linspace(Tr[0],Tr[-1],100)
			kappa_max = self.rxnUnsert[i].getKappaMax(T)			
			kappa_min = self.rxnUnsert[i].getKappaMin(T)
						
			Theta = np.array([T/T,np.log(T),-1/T])
			P = self.rxnUnsert[i].getMean()	
			#print(P)
			kappa_0 = Theta.T.dot(P)
			cov = self.rxnUnsert[i].getCov()
			zeta_A = self.rxnUnsert[i].zeta.x
			
			#print(os.getcwd())
			"""
			Type-A zeta samples
			"""
			generators_a = g_a
			parallel_zetas_a = []
			
			time_start = time.time()
			
			
			for k in g_a:
				Pint = P+np.asarray(np.dot(cov,k)).flatten()
				plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.25)
			#	plt.plot()
				#print(i)
				#print(type(i))
				#print(k)	
				#parallel_zetas_a.append(P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten())
				#parallel_zetas_a.append(k[0]*zeta_A)
								
			#print(parallel_zetas_a)
			plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"k-",linewidth=1.25,label=r"A-type samples ($\beta\zeta_{A}$)")
			plt.plot(1/T,kappa_max,"k--",linewidth=1.25,label=r"Rate uncertainty limits ($\kappa(\pm\zeta_{A}$))")
			plt.plot(1/T,kappa_min,"k--",linewidth=1.25)
			plt.plot(1/T,kappa_0,"b-",linewidth=1.25,label="Nominal rate value")
			#plt.xlim([1/2500,1/900])			
			#plt.ylim([23.3,29.67])				
			plt.xlabel(r"1000/T \ $K^{-1}$")
			plt.ylabel(r"$log_{10}(k)\,/\,s^{-1}\,or\,log_{10}(k)\, /\, cm^{3} \, mol^{-1} \, s^{-1}$")
			plt.legend()
			plt.savefig("./Special_cases/"+str(i)+"_a.pdf")
			
			
