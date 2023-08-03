
import Uncertainty # get the uncertainty extractor
import numpy as np
import math
import multiprocessing
import subprocess
import time
import sys
import copy

class DesignMatrix(object):
	def __init__(self,UnsrtData,design,sample_length,ind):
		self.unsrt = UnsrtData
		self.sim = sample_length
		self.design = design
		self.n = ind
		self.allowed_count = 100
	def getClassA_Curves(self,n_a):
		"""
			This defination generates n_a numbers of class-A type curves 
		"""
		ClassACurve_dict = {}
		Generator_dict = {}
		for rxn in self.unsrt:
			generator = np.random.random_sample((n_a,1))
			data = self.unsrt[rxn].data
			zeta_A = self.unsrt[rxn].zeta.x
			ClassA_curves = []
			for gen in generator:
				ClassA_curves.append(gen[0]*zeta_A)
			ClassACurve_dict[rxn] = ClassA_curves
			Generator_dict[rxn] = generator	
		return ClassACurve_dict,Generator_dict
			
	def getClassB_Curves(self,n_b):
		"""
			This defination generates n_b numbers of class-B type curves
			ClassB_curves: np.asarray, size of (450,3)
		"""
		ClassBCurve_dict = {}
		Generator_dict = {}
		for rxn in self.unsrt:
			generator = np.random.random_sample((n_b,2))
			data = self.unsrt[rxn].data
			data["generators_b"] = generator
			callWorkForce = Worker(self.allowed_count)	
			ClassB_curves,generator_b = callWorkForce.do_unsrt_b(data,len(generator))
			del callWorkForce
			ClassBCurve_dict[rxn] = np.asarray(ClassB_curves)
			Generator_dict[rxn] = generator_b	
		return ClassBCurve_dict,Generator_dict
	
	def getClassC_Curves(self,n_c):
		"""
			This defination generates n_c numbers of class-C type curves 
		"""
		
		ClassC_Curve_dict = {}
		Generator_dict = {}
		for rxn in self.unsrt:
			generator = np.random.random_sample((n_c,2))
			data = self.unsrt[rxn].data
			data["generators_c"] = generator
			callWorkForce = Worker(self.allowed_count)	
			ClassC_curves,generator_c = callWorkForce.do_unsrt_c(data,len(generator))
			del callWorkForce
			ClassC_Curve_dict[rxn] = np.asarray(ClassC_curves)
			Generator_dict[rxn] = generator_c	
		return ClassC_Curve_dict,Generator_dict
	
	def generatePointOnSphere(self):
		temp = np.random.random(2)
		u1 = temp[0]
		u2 = temp[1]
		lambda_ = np.arccos(2*u1-1)-math.pi/2
		phi = 2*math.pi*u2
		x = np.cos(lambda_)*np.cos(phi)
		y = np.cos(lambda_)*np.sin(phi)
		z = np.sin(lambda_)
		#print(x,y,z)
		return x,y,z
	
	
	def getSamples(self):
		
		if self.design == 'A-facto': 
			design_matrix = list(2* np.random.random_sample((self.sim*5,self.n)) - 1)
			design_matrix.extend(list(np.eye(self.n)))
			s =""
			for row in design_matrix:
				for element in row:
					s+=f"{element},"
				s+="\n"
			ff = open('DesignMatrix.csv','w').write(s)
			return np.asarray(design_matrix)
		elif self.design == "A1+B1+C1":
			design_matrix = []
			#Take the sample_length and divide it into following
			# 		n_a:n_b:n_c = 100:450:450
			#		Random shuffle = 4000
			#       Linear combination =  5000
			#      
			#      V: Unshuffled vector (numpy-array)
			#	 V_s: Shuffled vector (
			
			n_a = int(0.1*self.sim)
			n_b = int(0.45*self.sim)
			n_c = self.sim-n_a-n_b
			
			a_curves_dict, generator_a = self.getClassA_Curves(n_a)# Returns 100 class-A Arrhenius samples
			b_curves_dict, generator_b = self.getClassB_Curves(n_b)# Returns 450 class-B Arrhenius samples
			c_curves_dict, generator_c = self.getClassC_Curves(n_c)# Returns 450 class-C Arrhenius samples
					
			V = {}#Populating V for unshuffled portion of the design matrix
			for rxn in self.unsrt:
				temp = []
				for sample_a in a_curves_dict[rxn]: 
					#print(np.shape(sample_a))
					temp.append(sample_a)
					
				for sample_b in b_curves_dict[rxn]: 
					temp.append(sample_b)
				for sample_c in c_curves_dict[rxn]: 
					temp.append(sample_c)
				V[rxn] = np.asarray(temp)
			
			
			V_s = {}#Doing random shuffling
			#Deepcopy the unshuffled samples first
			V_copy = copy.deepcopy(V)
			for rxn in self.unsrt:				
				column = []
				for i in range(4):
					np.random.shuffle(V_copy[rxn])
					column.extend(np.asarray(V_copy[rxn]))
				V_s[rxn] = np.asarray(column)	
			
				
			V_linear_comb = {}#Doing linear combination to populate the matrix
			for rxn in self.unsrt:
				temp = []
				for i in range(5000):
					zeta_a = np.array(a_curves_dict[rxn][np.random.randint(0,100)])
					zeta_b = np.array(b_curves_dict[rxn][np.random.randint(0,450)])
					zeta_c = np.array(c_curves_dict[rxn][np.random.randint(0,450)])
					x,y,z = self.generatePointOnSphere()
					####
					temp.append(x*zeta_a+y*zeta_b+z*zeta_c)
				
				V_linear_comb[rxn] = np.asarray(temp)
				
					
			for i in range(n_a+n_b+n_c):
				row = []
				for rxn in self.unsrt:
					row.extend(V[rxn][i])
					
				design_matrix.append(row)
			
			for i in range(4000):
				row = []
				for rxn in self.unsrt:
					row.extend(V_s[rxn][i])
				design_matrix.append(row)
				
			for i in range(5000):
				row = []
				for rxn in self.unsrt:
					row.extend(V_linear_comb[rxn][i])
				design_matrix.append(row)
			
			s =""
			for row in design_matrix:
				for element in row:
					s+=f"{element},"
				s+="\n"
			ff = open('DesignMatrix.csv','w').write(s)	
			
			return np.asarray(design_matrix)
	
#########################################
## Definitions to run the generation   ##
## of samples using parallel computing ##
#########################################		
	
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
		return self.parallized_zeta,self.generator
	def do_unsrt_c(self,data,sampling_points):
		#generators = data["generators"]
		for args in range(sampling_points):
			self.pool.apply_async(run_sampling_c, 
				  args=(1,data,data["generators_c"][args],sampling_points,), 
				  callback=self.callback,error_callback=self.callback_error)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
		return self.parallized_zeta,self.generator

		
