import os
import Uncertainty # get the uncertainty extractor
import numpy as np
import math
import multiprocessing
import multiprocessing as mp
from multiprocessing import shared_memory
import subprocess
import time
import sys
import copy
import matplotlib.pyplot as plt
from tqdm import tqdm
import pickle
import tempfile
sys.path.append('/shuffle.so')
import shuffle
#print(dir(shuffle))
class DesignMatrix(object):
	def __init__(self,UnsrtData,design,ind,sample_length=None):
		self.unsrt = UnsrtData
		self.sim = sample_length
		self.design = design
		#self.n = ind
		self.rxn_len = ind
		self.n = ind
		allowed_count = int(0.70*multiprocessing.cpu_count())
		self.allowed_count = allowed_count
	
	def main(self,V_, num_threads, unsrt, sim):
		shuffle.shuffle_arrays(V_, num_threads, unsrt, sim)
		return V_
	def create_memory_mapped_files(self,data_dict):
		mmap_dict = {}
		temp_dir = tempfile.mkdtemp()
		for key, value in data_dict.items():
			filename = os.path.join(temp_dir, f"{key}.dat")
			mmap = np.memmap(filename, dtype=value.dtype, mode='w+', shape=value.shape)
			mmap[:] = value[:]
			mmap_dict[key] = filename
		return mmap_dict, temp_dir
	def create_shared_memory_dict(self,data_dict):
		shm_dict = {}
		for key, value in data_dict.items():
			shm = shared_memory.SharedMemory(create=True, size=value.nbytes)
			np_array = np.ndarray(value.shape, dtype=value.dtype, buffer=shm.buf)
			np_array[:] = value[:]
			shm_dict[key] = (shm, np_array)
		return shm_dict

	def cleanup_shared_memory(shm_dict):
		for shm, _ in shm_dict.values():
			shm.close()
			shm.unlink()
		
	def shuffle_values(self,rxn, filename, shape, dtype, num_shuffles):
		values = np.memmap(filename, dtype=dtype, mode='r+', shape=shape)

		column = []
		for _ in range(num_shuffles):
			np.random.shuffle(values)
			column.append(values.copy())

		column = np.concatenate(column)
		return rxn, column
	
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
	
	def getClassC_Curves(self,n_c,generator=np.array([])):
		"""
			This defination generates n_c numbers of class-C type curves 
		"""
		
		ClassC_Curve_dict = {}
		Generator_dict = {}
		
		for rxn in self.unsrt:
			if len(generator) == 0:
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
	
	def getSA_samples(self,factor):
		design_matrix = []
		design_matrix.extend(list(float(factor)*np.eye(self.rxn_len)))
		return np.asarray(design_matrix)
	def getSA3P_samples(self,zeta_vector,flag):
		#new_zeta = []
		#if flag == "A":
		#	for i in zeta_vector:
		#		new_zeta.append([abs(i[0]),0.0,0.0])
		#elif flag == "n":
		#	for i in zeta_vector:
		#		new_zeta.append([0.0,abs(i[1]),0.0])
		#elif flag == "Ea":
		#	for i in zeta_vector:
		#		new_zeta.append([0.0,0.0,abs(i[2])])
		#else:
		#	raise AssertionError("Please give a correct flag for DesignMatrix SA3P samples!!")
		
		# Number of vectors
		num_vectors = len(zeta_vector)

		# Size of each vector (assume all vectors are the same size)
		vector_size = len(zeta_vector[0])

		# Create an empty matrix of the appropriate size (3x9 in this case)
		design_matrix = np.zeros((num_vectors, num_vectors * vector_size))

		# Place each vector on the diagonal
		for i in range(num_vectors):
		    # Insert the vector into the appropriate diagonal position
		    design_matrix[i, i*vector_size:(i+1)*vector_size] = zeta_vector[i]
			
		return np.asarray(design_matrix)
	
	def getB_TYPE_fSAC(self,n_fSAC,tag="DICT"):
		if tag == "DICT":
			N_fSAC = n_fSAC
			nominal,ch,zeta,Temp = get_VARIABLES_FOR_fSAC(self.unsrt)
			V = get_B_TYPE_from_fSAC(nominal,ch,zeta,Temp,self.unsrt,N_fSAC)
			return V
		elif tag == "MAX_N":
			n_max = {}
			t_start = time.time()
			N_fSAC = n_fSAC
			nominal,ch,zeta,Temp = get_VARIABLES_FOR_fSAC(self.unsrt)
			V = get_B_TYPE_from_fSAC(nominal,ch,zeta,Temp,self.unsrt,N_fSAC)
			for rxn in self.unsrt:
				temp = V[rxn]
				temp_n = []
				L = ch[rxn]
				for sample in temp:
					z_b = L @ sample
					n_ = float(z_b[1])
					temp_n.append(n_)
				index = temp_n.index(max(temp_n))
				n_max[rxn] = temp[index]
			t_end = time.time()
			print(f"Time taken to find the B-type zeta from fSAC method = {t_start}-{t_end}")
			return n_max
		else:
			raise AssertionError("In DesignMatrix.py, please give proper tag for 'getB_TYPE_fSAC'")
	def getNominal_samples(self):
		design_matrix = []
		design_matrix.append(list(np.zeros(self.rxn_len)))
		return	design_matrix
	
	def getSamples(self):
		print("\nStarting to generate design matrix!!\n")
		if self.design == 'A-facto': 
			#self.sim = 1 + 2* self.n+self.n*(self.n-1)/2
			design_matrix = list(2* np.random.random_sample((self.sim,self.n)) - 1)
			design_matrix.extend(list(np.eye(self.n)))
			return np.asarray(design_matrix)
		elif self.design == "A1+B1+C1":
			tic = time.time()
			design_matrix = []
			#Take the sample_length and divide it into following
			# 		n_a:n_b:n_c = 10:45:45
			#               Original samples = 3000
			#		Random shuffle = 4000
			#              Linear combination =  5000
			#      
			#      V: Unshuffled vector (numpy-array)
			#	 V_s: Shuffled vector (
			
			n_a = int(0.1*self.sim*0.2)
			n_b = int(0.45*self.sim*0.2)
			n_c = int(self.sim*0.2-n_a-n_b)
			
			if "SA_a_type_samples.pkl" not in os.listdir():
				a_curves_dict, generator_a = self.getClassA_Curves(n_a)# Returns 100 class-A Arrhenius samples
				with open('SA_a_type_samples.pkl', 'wb') as file_:
					pickle.dump(a_curves_dict, file_)
				print("\nA-type curves generated")
			else:
				# Load the object from the file
				with open('SA_a_type_samples.pkl', 'rb') as file_:
					a_curves_dict = pickle.load(file_)
				print("\nA-type curves generated")
			
			if "SA_b_type_samples.pkl" not in os.listdir():	
				b_curves_dict, generator_b = self.getClassB_Curves(n_b)# Returns 450 class-B Arrhenius samples
				with open('SA_b_type_samples.pkl', 'wb') as file_:
					pickle.dump(b_curves_dict, file_)
				print("\nB-type curves generated")
			else:
				# Load the object from the file
				with open('SA_b_type_samples.pkl', 'rb') as file_:
					b_curves_dict = pickle.load(file_)
				print("\nB-type curves generated")
			
			if "SA_c_type_samples.pkl" not in os.listdir():
				c_curves_dict, generator_c = self.getClassC_Curves(n_c)# Returns 450 class-C Arrhenius samples	
				with open('SA_c_type_samples.pkl', 'wb') as file_:
					pickle.dump(c_curves_dict, file_)
				print("\nC-type curves generated")
			else:
				# Load the object from the file
				with open('SA_c_type_samples.pkl', 'rb') as file_:
					c_curves_dict = pickle.load(file_)
				print("\nC-type curves generated")
				
			V_ = {}
			for rxn in tqdm(self.unsrt,desc="Populating V_"):
				temp = []
				for sample_a in a_curves_dict[rxn]: 
					temp.append(sample_a)
				for sample_b in b_curves_dict[rxn]: 
					temp.append(sample_b)
				for sample_c in c_curves_dict[rxn]: 
					temp.append(sample_c)
				V_[rxn] = np.asarray(temp)
			
			V_s = {}#Doing random shuffling
			V_copy = copy.deepcopy(V_)
			for rxn in tqdm(self.unsrt, desc="Doing random shuffling"):
				column = []
				# Number of shuffles
				num_shuffles = int(self.sim * 0.8)
				# Get the values to shuffle
				values = V_copy[rxn]

				# Repeat the shuffle process
				for i in range(4):
					np.random.shuffle(values)
					column.append(values.copy())

				# Flatten the list of arrays
				column = np.concatenate(column)
				V_s[rxn] = column
		
			V_opt = {}
			#SAC_3
			#f-SAC algorithm
				#1] Get P_right, P_mid and P_left
				#2] get Kappa_right, Kappa_mid and Kappa_left
				#3] Solve for B-type zeta
				#4] Check if zeta_b is inside the domain
				#5] Check if n_ is < 2.
				#6] Do it for n_b number of times
			#plot the zetas
			#print(percentage)
			#Solving for a line passing from 3 points
			N_fSAC = 100
			V = self.getB_TYPE_fSAC(N_fSAC,tag="DICT")
				
			for i in range(int(self.sim*0.2)):
				row = []
				for rxn in self.unsrt:
					row.extend(V_[rxn][i])
					
				design_matrix.append(row)
			
			
			#RANDOM SHUFFLING
			for i in range(int(self.sim*0.8)):
				row = []
				for rxn in self.unsrt:
					row.extend(V_s[rxn][i])
				design_matrix.append(row)
				
			
			for i in range(N_fSAC):
				row = []
				for rxn in self.unsrt:
					row.extend(V[rxn][i])
					
				design_matrix.append(row)
			tok = time.time()
			print("Time taken to construct Design Matrix: {}".format(tok - tic))	
			
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
def get_zeta_b(P,Theta_p,kappa_o,kappa,cov,M):
	y = kappa - kappa_o
	zeta_b = np.linalg.solve(M,y)
	P_found = P + np.asarray(np.dot(cov,zeta_b)).flatten()
	n_b = abs(P[1]-P_found[1])
	kappa_found = Theta_p.T.dot(P_found)
	return zeta_b,n_b,P_found,kappa_found
def getP_for_fSAC(P,cov,zeta):
	random = 2*np.random.rand(3)-1
	P_r = P + random[0]*np.asarray(np.dot(cov,zeta)).flatten()
	P_m = P + random[1]*(7/8)*np.asarray(np.dot(cov,zeta)).flatten()
	P_l = P + random[2]*np.asarray(np.dot(cov,zeta)).flatten()
	return P_r,P_m,P_l	

def get_VARIABLES_FOR_fSAC(unsrt):
	ch = {}
	nominal = {}
	p_max = {}
	p_min = {}
	theta = {}
	Temp ={}
	zeta = {}
	for rxn in tqdm(unsrt,desc="Populating unshuffled portion of DM"):
		ch[rxn] = unsrt[rxn].cholskyDeCorrelateMat
		nominal[rxn] = unsrt[rxn].nominal
		p_max[rxn] = unsrt[rxn].P_max
		p_min[rxn] = unsrt[rxn].P_min
		theta[rxn] = unsrt[rxn].Theta
		Temp[rxn] = unsrt[rxn].temperatures
		zeta[rxn] = unsrt[rxn].zeta.x
	return nominal,ch,zeta,Temp

def get_B_TYPE_from_fSAC(nominal,ch,zeta,Temp,unsrt,N_fSAC):
	V = {}
	d_n = {}
	delta_n = {}
	_delta_n = {}
	dict_delta_n = {}
	percentage = {}
	for rxn in tqdm(unsrt,desc="Generating fSAC samples"):
		T = np.array([Temp[rxn][0],(Temp[rxn][0] + Temp[rxn][-1])/2,Temp[rxn][-1]])
		Theta = np.array([T/T,np.log(T),-1/(T)])
		P = nominal[rxn]
		cov = ch[rxn]
		zet = zeta[rxn]
		P_max = P + np.asarray(np.dot(cov,zet)).flatten();
		Tp = Temp[rxn]
		M = Theta.T.dot(cov)
		
		Theta_p = np.array([Tp/Tp,np.log(Tp),-1/(Tp)])
		kmax = Theta_p.T.dot(P_max)
		ka_o = Theta_p.T.dot(P)
		f = abs(kmax-ka_o)
		temp = []
		outside = []
		not_selected = []
		temp_n = []
		
		while len(temp)<N_fSAC:
			P_right,P_mid,P_left = getP_for_fSAC(P,cov,zet)
			kappa_left,kappa_mid,kappa_right = getKappa_for_fSAC(P_left,P_mid,P_right,T)
			kappa_o1,kappa_o2,kappa_o3 = getKappa_for_fSAC(P,P,P,T)
			kappa = np.array([kappa_left,kappa_mid,kappa_right])
			kappa_o = np.array([kappa_o1,kappa_o2,kappa_o3])
			zeta_b,n_,P_found,kappa_found = get_zeta_b(P,Theta_p,kappa_o,kappa,cov,M)
			temp_n.append(n_)
			f_found = abs(kappa_found-ka_o)
			if max(f)<max(f_found):
				outside.append(zeta_b)
			if n_ > 0.3:	
				#plt.plot(1/Tp,kappa_found,"r-",linewidth=0.25)
				#temp.append(zeta_)#for unconstrained zeta
				not_selected.append(zeta_b)
				_delta_n[n_] = zeta_b
				
			else:
				temp.append(zeta_b)
				_delta_n[n_] = zeta_b
			#print(zeta_)
		percentage[rxn] = (len(not_selected)/len(temp))*100
		delta_n[rxn] = max(temp_n)
		
		dict_delta_n[rxn] = _delta_n[max(temp_n)]
		V[rxn] = temp
	return V
def getKappa_for_fSAC(P_l,P_m,P_r,T):
	Theta_left = np.array([T[0]/T[0],np.log(T[0]),-1/(T[0])])
	Theta_mid = np.array([T[1]/T[1],np.log(T[1]),-1/(T[1])])
	Theta_right = np.array([T[2]/T[2],np.log(T[2]),-1/(T[2])])
	K_l = Theta_left.T.dot(P_l)
	K_m = Theta_mid.T.dot(P_m)
	K_r = Theta_right.T.dot(P_r)
	return K_l,K_m,K_r

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

		
