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
	def __init__(self,UnsrtData,design,sample_length,ind=None):
		self.unsrt = UnsrtData
		self.sim = sample_length
		self.design = design
		#self.n = ind
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
			generator = np.random.random_sample((n_a,1)) #n_a: number of A-type curves required 
			data = self.unsrt[rxn].data
			zeta_max = self.unsrt[rxn].zeta.x
			ClassA_curves = []
			for gen in generator:
				ClassA_curves.append(gen[0]*zeta_max)#getting new A-type curves
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
	
	def getSA_samples(self,factor):
		design_matrix = []
		design_matrix.extend(list(float(factor)*np.eye(self.rxn_len)))
		return np.asarray(design_matrix)
	
		
	
	def getNominal_samples(self):
		design_matrix = []
		design_matrix.append(list(np.zeros(self.rxn_len)))
		return	design_matrix
	
	def getSample_partial(self,case_index,selected_params):	
		if self.design == 'A-facto': 
			if "DM_FOR_PARTIAL_PRS" not in os.listdir():
				os.mkdir("DM_FOR_PARTIAL_PRS")
				os.chdir("DM_FOR_PARTIAL_PRS")
				if f"{case_index}" not in os.listdir():
					os.mkdir(case_index)
				os.chdir("..")
			else:
				os.chdir("DM_FOR_PARTIAL_PRS")
				if f"{case_index}" not in os.listdir():
					os.mkdir(case_index)
				os.chdir("..")
			if "DesignMatrix.pkl" not in os.listdir(f"DM_FOR_PARTIAL_PRS/{case_index}"):
				#self.sim = 1 + 2* self.n+self.n*(self.n-1)/2
				design_matrix = list(2* np.random.random_sample((self.sim,self.n)) - 1)
				#design_matrix.extend(list(np.eye(self.n)))
				with open(f'DM_FOR_PARTIAL_PRS/{case_index}/DesignMatrix.pkl', 'wb') as file_:
					pickle.dump(design_matrix, file_)
			else: 
				with open(f'DM_FOR_PARTIAL_PRS/{case_index}/DesignMatrix.pkl', 'rb') as file_:
					design_matrix = pickle.load(file_)
				print("\nA-factor samples are generated")

			s =""
			p_s = "" 
			p_ss = ""
			ss = ""
			p_design_matrix = []
			p_selection_matrix = []
			selection_matrix = []
			for row in design_matrix:
				row_ = []
				temo = []
				temp = []
				for index,element in enumerate(row):
					if selected_params[index] == 1:
						s+=f"{element},"
						p_s+=f"{element},"
						p_ss+="1.0,"
						ss+="1.0,"
						row_.append(element)
						temo.append(1.0)
						temp.append(1.0)
					else:
						s+="0.0,"
						p_ss+="0.0,"
						ss+="1.0,"
						temo.append(0.0)
						temp.append(1.0)
				p_design_matrix.append(row_)
				p_selection_matrix.append(temo)
				selection_matrix.append(temp)
				s+="\n"
				p_s+="\n"
				p_ss+="\n"
				ss+="\n"

			ff = open(f'DM_FOR_PARTIAL_PRS/{case_index}/DesignMatrix.csv','w').write(s)
			ff = open(f'DM_FOR_PARTIAL_PRS/{case_index}/pDesignMatrix.csv','w').write(p_s)
			ff = open(f'DM_FOR_PARTIAL_PRS/{case_index}/pSelectionMatrix.csv','w').write(p_ss)
			ff = open(f'DM_FOR_PARTIAL_PRS/{case_index}/SelectionMatrix.csv','w').write(ss)
			
			return np.asarray(design_matrix),np.asarray(selection_matrix),np.asarray(p_design_matrix),np.asarray(p_selection_matrix)
			#return np.asarray(design_matrix)
		
		elif self.design == "A1+B1+C1":
			tic = time.time()
			design_matrix = []
						
			n_a = int(0.1*self.sim*0.2)
			n_b = int(0.45*self.sim*0.2)
			n_c = int(self.sim*0.2)-n_a-n_b
			unshuffled = int(self.sim*0.2)
			
			if "DM_FOR_PARTIAL_PRS" not in os.listdir():
				os.mkdir("DM_FOR_PARTIAL_PRS")
				os.chdir("DM_FOR_PARTIAL_PRS")
				if f"{case_index}" not in os.listdir():
					os.mkdir(case_index)
				os.chdir("..")
			else:
				os.chdir("DM_FOR_PARTIAL_PRS")
				if f"{case_index}" not in os.listdir():
					os.mkdir(case_index)
				os.chdir("..")
			
			if "a_type_samples.pkl" not in os.listdir(f"DM_FOR_PARTIAL_PRS/{case_index}"):
				#print(os.getcwd())
				a_curves_dict, generator_a = self.getClassA_Curves(n_a)# Returns 100 class-A Arrhenius samples
				with open(f'DM_FOR_PARTIAL_PRS/{case_index}/a_type_samples.pkl', 'wb') as file_:
					pickle.dump(a_curves_dict, file_)
				#print("\nA-type curves generated")
			else:
				# Load the object from the file
				with open(f'DM_FOR_PARTIAL_PRS/{case_index}/a_type_samples.pkl', 'rb') as file_:
					a_curves_dict = pickle.load(file_)
				print("\nA-type curves generated")
			
			if "b_type_samples.pkl" not in os.listdir(f"DM_FOR_PARTIAL_PRS/{case_index}"):	
				b_curves_dict, generator_b = self.getClassB_Curves(n_b)# Returns 450 class-B Arrhenius samples
				with open(f'DM_FOR_PARTIAL_PRS/{case_index}/b_type_samples.pkl', 'wb') as file_:
					pickle.dump(b_curves_dict, file_)
				#print("\nB-type curves generated")
			else:
				# Load the object from the file
				with open(f'DM_FOR_PARTIAL_PRS/{case_index}/b_type_samples.pkl', 'rb') as file_:
					b_curves_dict = pickle.load(file_)
				print("\nB-type curves generated")
			
			if "c_type_samples.pkl" not in os.listdir(f"DM_FOR_PARTIAL_PRS/{case_index}"):	
				c_curves_dict, generator_c = self.getClassC_Curves(n_c)# Returns 450 class-C Arrhenius samples	
				with open(f'DM_FOR_PARTIAL_PRS/{case_index}/c_type_samples.pkl', 'wb') as file_:
					pickle.dump(c_curves_dict, file_)
				#print("\nC-type curves generated")
			else:
				# Load the object from the file
				with open(f'DM_FOR_PARTIAL_PRS/{case_index}/c_type_samples.pkl', 'rb') as file_:
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
			
			num_shuffles = 4*len(temp)
			if "RANDOM_SHUFFLING.pkl" not in os.listdir(f"DM_FOR_PARTIAL_PRS/{case_index}"):	
				V_s = {}#Doing random shuffling
				V_copy = copy.deepcopy(V_)
				for rxn in tqdm(self.unsrt, desc="Doing random shuffling"):
					column = []
					# Number of shuffles
					#num_shuffles = int(self.sim * 0.8)
					# Get the values to shuffle
					values = V_copy[rxn]

					# Repeat the shuffle process
					for i in range(4):
						np.random.shuffle(values)
						column.append(values.copy())

					# Flatten the list of arrays
					column = np.concatenate(column)
					#num_shuffles = len(column)
					V_s[rxn] = column
								
				with open(f'DM_FOR_PARTIAL_PRS/{case_index}/RANDOM_SHUFFLING.pkl', 'wb') as file_:
					pickle.dump(V_s, file_)
				#print("\nC-type curves generated")
			else:
				
				
				with open(f'DM_FOR_PARTIAL_PRS/{case_index}/RANDOM_SHUFFLING.pkl', 'rb') as file_:
					V_s = pickle.load(file_)
				print("\nRandom shuffling Arrhenius curves generated")
			delta_n = {}
			p = {}
			V_opt = {}
			V = {}#Populating V for unshuffled portion of the design matrix
			ch = {}
			nominal = {}
			p_max = {}
			p_min = {}
			theta = {}
			Temp ={}
			zeta = {}
			for rxn in tqdm(self.unsrt,desc="Populating unshuffled portion of DM"):
				ch[rxn] = self.unsrt[rxn].cholskyDeCorrelateMat
				nominal[rxn] = self.unsrt[rxn].nominal
				p_max[rxn] = self.unsrt[rxn].P_max
				p_min[rxn] = self.unsrt[rxn].P_min
				theta[rxn] = self.unsrt[rxn].Theta
				Temp[rxn] = self.unsrt[rxn].temperatures
				zeta[rxn] = self.unsrt[rxn].zeta.x
				
			
			if "fSAC_samples.pkl" not in os.listdir(f"DM_FOR_PARTIAL_PRS/{case_index}"):	
				d_n = {}
				_delta_n = {}
				dict_delta_n = {}
				percentage = {}
				for rxn in tqdm(self.unsrt,desc="Generating fSAC samples"):
					T = np.array([Temp[rxn][0],(Temp[rxn][0] + Temp[rxn][-1])/2,Temp[rxn][-1]])	
					Theta = np.array([T/T,np.log(T),-1/(T)])
					P = nominal[rxn]
					cov = ch[rxn]
					zet = zeta[rxn]
					Tp = Temp[rxn]
					#Tp = np.linspace(300,2500,50)
					Theta_p = np.array([Tp/Tp,np.log(Tp),-1/(Tp)])
					P_max = P + np.asarray(np.dot(cov,zet)).flatten();
					P_min = P - np.asarray(np.dot(cov,zet)).flatten();
					
					d_n[rxn] = abs(P[1]-P_max[1])
					kmax = Theta_p.T.dot(P_max)
					kmin = Theta_p.T.dot(P_min)
					ka_o = Theta_p.T.dot(P)
					M = Theta.T.dot(cov)
					f = abs(kmax-ka_o)
					temp = []
					outside = []
					not_selected = []
					temp_n = []
					while len(temp)<int(self.sim*0.2):
						random = 2*np.random.rand(3)-1
						P_right = P + random[0]*np.asarray(np.dot(cov,zet)).flatten()
						P_mid = P + random[1]*(7/8)*np.asarray(np.dot(cov,zet)).flatten()
						P_left = P + random[2]*np.asarray(np.dot(cov,zet)).flatten()
						Theta_left = np.array([T[0]/T[0],np.log(T[0]),-1/(T[0])])
						Theta_mid = np.array([T[1]/T[1],np.log(T[1]),-1/(T[1])])
						Theta_right = np.array([T[2]/T[2],np.log(T[2]),-1/(T[2])])
						kappa_left = Theta_left.T.dot(P_left)
						kappa_mid = Theta_mid.T.dot(P_mid)
						kappa_right = Theta_right.T.dot(P_right)
						kappa_o1 = Theta_left.T.dot(P)
						kappa_o2= Theta_mid.T.dot(P)
						kappa_o3 = Theta_right.T.dot(P)
						kappa = np.array([kappa_left,kappa_mid,kappa_right])
						kappa_o = np.array([kappa_o1,kappa_o2,kappa_o3])
						y = kappa - kappa_o
						zeta_ = np.linalg.solve(M,y)
						func = np.asarray([(i.T.dot(cov.dot(zeta_))) for i in Theta_p.T])
						P_found = P + np.asarray(np.dot(cov,zeta_)).flatten()
						n_ = abs(P[1]-P_found[1])
						temp_n.append(n_)
						kappa_found = Theta_p.T.dot(P_found)
						f_found = abs(kappa_found-ka_o)
						if max(f)<max(f_found):
							outside.append(zeta_)
						if n_ > 0.3:	
							#plt.plot(1/Tp,kappa_found,"r-",linewidth=0.25)
							#temp.append(zeta_)#for unconstrained zeta
							not_selected.append(zeta_)
							_delta_n[n_] = zeta_
							
						else:
							temp.append(zeta_)
							_delta_n[n_] = zeta_
						#print(zeta_)
					percentage[rxn] = (len(not_selected)/len(temp))*100
					delta_n[rxn] = max(temp_n)
					
					dict_delta_n[rxn] = _delta_n[max(temp_n)]
					V[rxn] = temp
				with open(f'DM_FOR_PARTIAL_PRS/{case_index}/fSAC_samples.pkl', 'wb') as file_:
					pickle.dump(V, file_)
				#print("\nC-type curves generated")
			else:
				# Load the object from the file
				with open(f'DM_FOR_PARTIAL_PRS/{case_index}/fSAC_samples.pkl', 'rb') as file_:
					V = pickle.load(file_)
				print("\nf-SAC Arrhenius curves generated")
				
			#UNSHUFFLED SAMPLES
			#print(unshuffled)
			for i in range(unshuffled):
				row = []
				for rxn in self.unsrt:
					#print((V_[rxn]))
					row.extend(V_[rxn][i])
				design_matrix.append(row)
			
			#RANDOM SHUFFLING
			for i in range(num_shuffles):
				row = []
				for rxn in self.unsrt:
					row.extend(V_s[rxn][i])
				design_matrix.append(row)
				
			for i in range(unshuffled):
				row = []
				for rxn in self.unsrt:
					row.extend(V[rxn][i])
				design_matrix.append(row)
				
			tok = time.time()
			print("Time taken to construct Design Matrix: {}".format(tok - tic))
			
			##Creating the parameter matrix
			s =""
			p_s = "" 
			selected_string = ""
			p_selected_string = ""
			p_design_matrix = []
			p_selection_matrix = []
			selection_matrix = []
			for row in design_matrix:
				row_ = []
				temo = []
				temp = []
				for index,element in enumerate(row):
					if selected_params[index] == 1:
						s+=f"{element},"
						p_s+=f"{element},"
						selected_string+="1.0,"
						p_selected_string+="1.0,"
						row_.append(element)
						temo.append(1.0)
						temp.append(1.0)
					else:
						selected_string+="1.0,"
						p_selected_string+="0.0,"
						s+=f"{element},"
						#p_s+=f"{element},"
						#row_.append(element)
						temo.append(0.0)
						temp.append(1.0)
				p_design_matrix.append(row_)
				p_selection_matrix.append(temo)
				selection_matrix.append(temp)
				s+="\n"
				p_s+="\n"
				selected_string+="\n"
				p_selected_string+="\n"
			ff = open(f'DM_FOR_PARTIAL_PRS/{case_index}/DesignMatrix.csv','w').write(s)	
			ff = open(f'DM_FOR_PARTIAL_PRS/{case_index}/pDesignMatrix.csv','w').write(p_s)
			ff = open(f'DM_FOR_PARTIAL_PRS/{case_index}/pSelectionMatrix.csv','w').write(p_selected_string)
			ff = open(f'DM_FOR_PARTIAL_PRS/{case_index}/SelectionMatrix.csv','w').write(selected_string)
			return np.asarray(design_matrix),np.asarray(selection_matrix),np.asarray(p_design_matrix),np.asarray(p_selection_matrix)
	
	
	def getSamples(self):
		print("\nStarting to generate design matrix!!\n")
		if self.design == 'A-facto': 
			#self.sim = 1 + 2* self.n+self.n*(self.n-1)/2
			design_matrix = list(2* np.random.random_sample((self.sim,self.n)) - 1)
			design_matrix.extend(list(np.eye(self.n)))
			s =""
			for row in design_matrix:
				for element in row:
					s+=f"{element},"
				s+="\n"
			ff = open('DesignMatrix.csv','w').write(s)
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
			
			if "a_type_samples.pkl" not in os.listdir():
				a_curves_dict, generator_a = self.getClassA_Curves(n_a)# Returns 100 class-A Arrhenius samples
				with open('a_type_samples.pkl', 'wb') as file_:
					pickle.dump(a_curves_dict, file_)
				print("\nA-type curves generated")
			else:
				# Load the object from the file
				with open('a_type_samples.pkl', 'rb') as file_:
					a_curves_dict = pickle.load(file_)
				print("\nA-type curves generated")
			
			if "b_type_samples.pkl" not in os.listdir():	
				b_curves_dict, generator_b = self.getClassB_Curves(n_b)# Returns 450 class-B Arrhenius samples
				with open('b_type_samples.pkl', 'wb') as file_:
					pickle.dump(b_curves_dict, file_)
				print("\nB-type curves generated")
			else:
				# Load the object from the file
				with open('b_type_samples.pkl', 'rb') as file_:
					b_curves_dict = pickle.load(file_)
				print("\nB-type curves generated")
			
			if "c_type_samples.pkl" not in os.listdir():
				c_curves_dict, generator_c = self.getClassC_Curves(n_c)# Returns 450 class-C Arrhenius samples	
				with open('c_type_samples.pkl', 'wb') as file_:
					pickle.dump(c_curves_dict, file_)
				print("\nC-type curves generated")
			else:
				# Load the object from the file
				with open('c_type_samples.pkl', 'rb') as file_:
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
			#Deepcopy the unshuffled samples first
			#print(n_a,n_b,n_c)
			#print(len(n_a),len(n_b),len(n_c))
			V_copy = copy.deepcopy(V_)
			#for rxn in V_:
			#	print(len(V_[rxn]))
			#	print(len(V_[rxn][0]))
			#	raise AssertionError("Ho stop!")
			#shm_dict = self.create_shared_memory_dict(V_copy)
			#num_shuffles = int(self.sim * 0.8)
			#reaction_list = []
			#for rxn in self.unsrt:
			#	reaction_list.append(rxn)
			#V_s = self.main(V_copy,100,reaction_list,self.sim)
			
			#for rxn in V_s:
			#	print(len(V_s[rxn]))
			#	print(len(V_s[rxn][0]))
			#	raise AssertionError("Ho stop!")
			#raise AssertionError("Shuffling done!!")
			
			#mmap_dict, temp_dir = self.create_memory_mapped_files(V_)
			#for rxn in tqdm(self.unsrt,desc="Doing random shuffling"):				
			#	column = []
			#	for i in range(int(self.sim*0.8)):
			#		np.random.shuffle(V_copy[rxn])
			#		column.extend(np.asarray(V_copy[rxn]))
			#	V_s[rxn] = np.asarray(column)	
			
			
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
			
			
			#for rxn in V_s:
			#	print(len(V_s[rxn]))
			
			#	print(len(V_s[rxn][0]))
			#	raise AssertionError("Ho stop!")
			#num_shuffles = int(self.sim * 0.8)
			# Initialize the result dictionary
			#V_s = {}

			# Create a pool of workers
			#pool = mp.Pool(mp.cpu_count()-int(0.1*mp.cpu_count()))#Leaving 10% of CPU power
			#results = []

			# Use apply_async to distribute the work
			#for rxn in tqdm(self.unsrt, desc="Doing random shuffling"):
			#	filename = mmap_dict[rxn]
				#shm, np_array = shm_dict[rxn]
			#	result = pool.apply_async(self.shuffle_values, args=(rxn, filename, V_[rxn].shape, V_[rxn].dtype, num_shuffles))
			#	results.append(result)

			# Close the pool and wait for the work to finish
			#pool.close()
			#pool.join()

			# Collect the results
			#for result in results:
			#	rxn, column = result.get()
			#	V_s[rxn] = column
			
			
			# Cleanup shared memory
			#self.cleanup_shared_memory(shm_dict)
			
			# Cleanup memory-mapped files
			#for filename in mmap_dict.values():
			#	os.remove(filename)
			#os.rmdir(temp_dir)
			
			
			"""
			V_linear_comb = {}#Doing linear combination to populate the matrix
			for rxn in self.unsrt:
				temp = []
				for i in range(1000):
					zeta_a = np.array(a_curves_dict[rxn][np.random.randint(0,10)])
					zeta_b = np.array(b_curves_dict[rxn][np.random.randint(0,10)])
					zeta_c = np.array(c_curves_dict[rxn][np.random.randint(0,10)])
					x,y,z = self.generatePointOnSphere()
					####
					temp.append(x*zeta_a+y*zeta_b+z*zeta_c)
				
				V_linear_comb[rxn] = np.asarray(temp)
				
			"""
									
			delta_n = {}
			p = {}
			V_opt = {}
			V = {}#Populating V for unshuffled portion of the design matrix
			ch = {}
			nominal = {}
			p_max = {}
			p_min = {}
			theta = {}
			Temp ={}
			zeta = {}
			for rxn in tqdm(self.unsrt,desc="Populating unshuffled portion of DM"):
				ch[rxn] = self.unsrt[rxn].cholskyDeCorrelateMat
				nominal[rxn] = self.unsrt[rxn].nominal
				p_max[rxn] = self.unsrt[rxn].P_max
				p_min[rxn] = self.unsrt[rxn].P_min
				theta[rxn] = self.unsrt[rxn].Theta
				Temp[rxn] = self.unsrt[rxn].temperatures
				zeta[rxn] = self.unsrt[rxn].zeta.x
				#temp = []
				#for sample_a in a_curves_dict[rxn]: 
					#print(np.shape(sample_a))
				#	temp.append(sample_a)
					
				#for sample_b in b_curves_dict[rxn]: 
				#	temp.append(sample_b)
				#for sample_c in c_curves_dict[rxn]: 
				#	temp.append(sample_c)
				#V[rxn] = np.asarray(temp)
			"""
			for rxn in self.unsrt:
				temp = []
				#temp_ = []
				temp_n = []
				nom = nominal[rxn]
				cholesky = ch[rxn]
				for sample_c in c_curves_dict[rxn]: 
					k = nom + np.asarray(cholesky.dot(sample_c)).flatten()
					temp_n.append(abs(nom[1]-k[1]))
					temp.append(sample_c)
				delta_n[rxn] = max(temp_n)
				V_opt[rxn] = np.asarray(temp)
			
			print(delta_n)
			raise AssertionError("delta_n part is ongoing")
			"""
			
			"""
			for rxn in self.unsrt:
				T = (Temp[rxn][0] + Temp[rxn][-1])/2	
				Theta = np.array([T/T,np.log(T),-1/(T)])
				kappa_max_mid = Theta.T.dot(p_max[rxn])
				kappa_min_mid = Theta.T.dot(p_min[rxn])
				kappa_nominal = Theta.T.dot(nominal[rxn])
				delta = abs(kappa_max_mid - kappa_nominal)
				kappa_list = []
				for i in V[rxn]:
					p = nominal[rxn] + ch[rxn].dot(i)
					delta_mid = abs(Theta.T.dot(p) - kappa_nominal)
					if delta_mid<(6/7)*delta:
						kappa_list.append(i)
					
				percentage[rxn] = (len(kappa_list)/len(V[rxn]))*100
			"""
			
			#SAC_3
			
			#plot the zetas
			#print(percentage)
			#Solving for a line passing from 3 points
			d_n = {}
			_delta_n = {}
			dict_delta_n = {}
			percentage = {}
			for rxn in tqdm(self.unsrt,desc="Generating fSAC samples"):
				T =np.array([Temp[rxn][0],(Temp[rxn][0] + Temp[rxn][-1])/2,Temp[rxn][-1]])	
				Theta = np.array([T/T,np.log(T),-1/(T)])
				P = nominal[rxn]
				cov = ch[rxn]
				zet = zeta[rxn]
				Tp = Temp[rxn]
				#Tp = np.linspace(300,2500,50)
				Theta_p = np.array([Tp/Tp,np.log(Tp),-1/(Tp)])
				P_max = P + np.asarray(np.dot(cov,zet)).flatten();
				P_min = P - np.asarray(np.dot(cov,zet)).flatten();
				
				d_n[rxn] = abs(P[1]-P_max[1])
				kmax = Theta_p.T.dot(P_max)
				kmin = Theta_p.T.dot(P_min)
				ka_o = Theta_p.T.dot(P)
				M = Theta.T.dot(cov)
				#fig = plt.figure()
				#plt.plot(1/Tp,kmax)
				#plt.plot(1/Tp,kmin)
				#plt.plot(1/Tp,ka_o)
				f = abs(kmax-ka_o)
				temp = []
				outside = []
				not_selected = []
				temp_n = []
				while len(temp)<1000:
					random = 2*np.random.rand(3)-1
					P_right = P + random[0]*np.asarray(np.dot(cov,zet)).flatten()
					P_mid = P + random[1]*(7/8)*np.asarray(np.dot(cov,zet)).flatten()
					P_left = P + random[2]*np.asarray(np.dot(cov,zet)).flatten()
					Theta_left = np.array([T[0]/T[0],np.log(T[0]),-1/(T[0])])
					Theta_mid = np.array([T[1]/T[1],np.log(T[1]),-1/(T[1])])
					Theta_right = np.array([T[2]/T[2],np.log(T[2]),-1/(T[2])])
					kappa_left = Theta_left.T.dot(P_left)
					kappa_mid = Theta_mid.T.dot(P_mid)
					kappa_right = Theta_right.T.dot(P_right)
					kappa_o1 = Theta_left.T.dot(P)
					kappa_o2= Theta_mid.T.dot(P)
					kappa_o3 = Theta_right.T.dot(P)
					kappa = np.array([kappa_left,kappa_mid,kappa_right])
					kappa_o = np.array([kappa_o1,kappa_o2,kappa_o3])
					y = kappa - kappa_o
					zeta_ = np.linalg.solve(M,y)
					func = np.asarray([(i.T.dot(cov.dot(zeta_))) for i in Theta_p.T])
					P_found = P + np.asarray(np.dot(cov,zeta_)).flatten()
					n_ = abs(P[1]-P_found[1])
					temp_n.append(n_)
					kappa_found = Theta_p.T.dot(P_found)
					f_found = abs(kappa_found-ka_o)
					if max(f)<max(f_found):
						outside.append(zeta_)
					if n_ > 2:	
						#plt.plot(1/Tp,kappa_found,"r-",linewidth=0.25)
						#temp.append(zeta_)#for unconstrained zeta
						not_selected.append(zeta_)
						_delta_n[n_] = zeta_
						
					else:
						temp.append(zeta_)
						_delta_n[n_] = zeta_
					#print(zeta_)
				percentage[rxn] = (len(not_selected)/len(temp))*100
				delta_n[rxn] = max(temp_n)
				
				dict_delta_n[rxn] = _delta_n[max(temp_n)]
				V[rxn] = temp
				
				#print(len(temp))
				
				#plt.show()
			#print(percentage)
			#print(percentage)	
			#print("------------------------------------\nNominal!!\n\n")
			#print(nominal)
			#print("------------------------------------\nCovariance!!\n\n")
			#print(ch)
			#print("------------------------------------\nZeta!!\n\n")
			#print(zeta)
			#print(Temp)
			#print(d_n)
			#print(delta_n)
			#print(dict_delta_n)
			#raise AssertionError("New")	
			
					
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
				
			
			for i in range(1000):
				row = []
				for rxn in self.unsrt:
					row.extend(V[rxn][i])
					
				design_matrix.append(row)
			tok = time.time()
			print("Time taken to construct Design Matrix: {}".format(tok - tic))
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

		
