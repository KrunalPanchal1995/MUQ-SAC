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
from filelock import FileLock
#print(dir(shuffle))
class DesignMatrix(object):
	def __init__(self,UnsrtData,design,sample_length, ind=None):
		self.unsrt = UnsrtData
		self.sim = sample_length
		self.design = design
		#self.n = ind
		#self.rxn_len = ind
		self.n = ind
		
		allowed_count = 100#int(0.70*multiprocessing.cpu_count())
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
		
	
	def get_thermo_ConstantCurves(self,n_a):
		"""
			This defination generates n_a numbers of class-A type curves 
		"""
		Class_thermo_Curve_dict = {}
		Generator_dict = {}
		
		for index,species in enumerate(self.unsrt):
			print(species)
			species_name = self.unsrt[species].species
			generator = []
			Cp_T_mid = []
			a1_ = []
			T = self.unsrt[species].common_temp
			Theta = np.array([T/T,T,T**2,T**3,T**4])
			#if self.unsrt[species].a1 == None:
			if self.unsrt[species].a1_C is None or np.all(self.unsrt[species].a1_C)==None:
				for i in range(n_a):
					a1 = np.random.uniform(-1,1,(1,1))[0]
					#print("\t NOMINAL \t",self.unsrt[species].NominalParams  )
					P_C = self.unsrt[species].NominalParams + a1*np.asarray(np.dot(self.unsrt[species].cov,self.unsrt[species].zeta_max.x)).flatten()
					Cp_T_mid.append(Theta.dot(P_C))
					a1_.append(a1)
					generator.append([a1])
				for spec in self.unsrt:
					if self.unsrt[spec].species == species_name:
						self.unsrt[spec].a1_C = a1_
						self.unsrt[spec].Cp_T_mid_C = Cp_T_mid
						#print("\n index", index  )
						#print("\n generator",generator )
						#print("\n self.unsrt",self.unsrt[spec].a1 	)			
					else:
						continue
					
			else:
				generator = []
				species_name = self.unsrt[species].species
			
				T = self.unsrt[species].common_temp
				Theta = np.array([T/T,T,T**2,T**3,T**4])
				#print(f"{species} high")
				#print("\n"+str(species))
				a1_ = self.unsrt[species].a1_C
				Cp_T_mid = self.unsrt[species].Cp_T_mid_C
				a1_dash = []
				for i in range(len(a1_)):
					a1 = ((Cp_T_mid[i] - Theta.dot(self.unsrt[species].NominalParams))/abs(Theta.dot(np.asarray(np.dot(self.unsrt[species].cov,self.unsrt[species].zeta_max.x)).flatten())))
					a1_dash.append(np.array([a1]))
					generator.append([a1])
				self.unsrt[species].a1_C = a1_dash
				
			data = self.unsrt[species].data
			zeta_A = self.unsrt[species].zeta_max.x
			ConstantCurves = []
			for gen in generator:
				ConstantCurves.append(np.asarray(gen*zeta_A).flatten())
				
			Class_thermo_Curve_dict[species] = ConstantCurves
			Generator_dict[species] = generator
			#print("Constant Curves \n")
			#print(ConstantCurves)
				
		return Class_thermo_Curve_dict,Generator_dict
		



			
	def get_thermo_LinearCurves(self,n_a):
		"""
			This defination generates n_a numbers of linear cp type curves 
		"""
		Class_thermo_Curve_dict = {}
		Generator_dict = {}
		for index,species in enumerate(self.unsrt):
			#print(species)
			
			if self.unsrt[species].a1_L is None or np.all(self.unsrt[species].a1_L) == None:
				species_name = self.unsrt[species].species
				
			
				T = self.unsrt[species].common_temp
				Theta = np.array([T/T,T,T**2,T**3,T**4])
				a1_ = []
				a2_ = []
				a3_ = []
				Cp_T_mid = []
				generator = []
				for i in range(n_a):
					a1 = np.random.uniform(-1,1,(1,1))[0]
					a2 = np.random.uniform(a1,1,(1,1))[0]
					P_2 = self.unsrt[species].NominalParams + a2*np.asarray(np.dot(self.unsrt[species].cov,self.unsrt[species].zeta_max.x)).flatten()
					Cp_T_mid.append(Theta.dot(P_2))
					a3 = np.random.uniform(a2,1,(1,1))[0]
					a1_.append(a1)
					a2_.append(a2)
					a3_.append(a3)
					generator.append([a1,a2,a3])
				
				for spec in self.unsrt:
					if self.unsrt[spec].species == species_name:
						self.unsrt[spec].a1_L = a1_
						self.unsrt[spec].a2_L = a2_
						self.unsrt[spec].a3_L = a3_
						self.unsrt[spec].Cp_T_mid_L = Cp_T_mid
							
					
			else:
				generator = []
				species_name = self.unsrt[species].species
				
			
				T = self.unsrt[species].common_temp
				Theta = np.array([T/T,T,T**2,T**3,T**4])
				#print(f"{species} high")
				#print("\n"+str(species))
				a1_ = self.unsrt[species].a1_L
				a2_ = self.unsrt[species].a2_L
				#print("\n")
				#print(a2_)
				a2_dash = []
				a3_dash = []
				a3_ = self.unsrt[species].a3_L
				Cp_T_mid = self.unsrt[species].Cp_T_mid_L
				for i in range(len(a1_)):
					a1 = a1_[i][0]
					a2 = ((Cp_T_mid[i] - Theta.dot(self.unsrt[species].NominalParams))/abs(Theta.dot(np.asarray(np.dot(self.unsrt[species].cov,self.unsrt[species].zeta_max.x)).flatten())))
					a2_dash.append(np.array([a2]))
					a3 = np.random.uniform(a2,1,(1,1))[0]
					a3_dash.append(a3)
					
					generator.append([[a1],[a2],a3])
				self.unsrt[species].a2_L = a2_dash
				self.unsrt[species].a3_L = a3_dash
			#print(generator)
			generator = np.asarray(generator,dtype="object")
			#print(len(generator))
			data = self.unsrt[species].data
			chunk_size = 3000
			params = [generator[i:i+chunk_size] for i in range(0, len(generator), chunk_size)]
			Class_thermo_curves = []
			generator_thermo = []
			#for set_ in tqdm(params,desc="Set. progress"):
			data["generators_thermo"] = generator
			callWorkForce = Worker(self.allowed_count)	
			Class_thermo_curves_,f_c,generator_thermo_ = callWorkForce.do_unsrt_thermoLinear(data,len(generator))
			Class_thermo_curves.extend(Class_thermo_curves_)
			generator_thermo.extend(generator_thermo_)
			del callWorkForce#,generator_thermo_,Class_thermo_curves_
			Class_thermo_Curve_dict[species] = np.asarray(Class_thermo_curves)
			Generator_dict[species] = generator_thermo
			#print("\n Generating Linear Curves \n")
			#print(Class_thermo_curves)
			#raise AssertionError("line 219- DM")
		return Class_thermo_Curve_dict,Generator_dict
	

	def get_thermo_Hybrid_Constant_LinearCurves(self, n_a):
		"""
			This defination generates n_a numbers of class-C_H type thermo curves 
		"""
		Class_thermo_Curve_dict = {}
		Generator_dict = {}
		
		for index,species in enumerate(self.unsrt):
			#print(species)
			species_name = self.unsrt[species].species
			generator = []
			Cp_T_mid = []
			a1_ = []
			nominal_curve = []
			T = self.unsrt[species].common_temp
			Theta = np.array([T/T,T,T**2,T**3,T**4])
			#if self.unsrt[species].a1 == None:
			if self.unsrt[species].a1_CL is None or np.all(self.unsrt[species].a1_CL)==None:
				for i in range(n_a):
					a1 = np.random.uniform(-1,1,(1,1))[0]
					P_C = self.unsrt[species].NominalParams + a1*np.asarray(np.dot(self.unsrt[species].cov,self.unsrt[species].zeta_max.x)).flatten()
					Cp_T_mid.append(Theta.dot(P_C))
					a1_.append(a1)
					generator.append([[a1]])
					#generator.append([a1])   ### this is the probable problem
				for spec in self.unsrt:
					if self.unsrt[spec].species == species_name:
						self.unsrt[spec].a1_CL = a1_
						self.unsrt[spec].a2_CL = a1_
						self.unsrt[spec].a3_CL = a1_
						self.unsrt[spec].Cp_T_mid_CL = Cp_T_mid
						#print("\n index", index  )
						#print("\n generator",generator )
						#raise AssertionError
					else:
						continue
				data = self.unsrt[species].data
				zeta_A = self.unsrt[species].zeta_max.x
				ConstantCurves = []
				for gen in generator:
					ConstantCurves.append(np.asarray(gen*zeta_A).flatten())
				Class_thermo_Curve_dict[species] = ConstantCurves
				Generator_dict[species] = generator
			else:
				generator = []
				species_name = self.unsrt[species].species
			
				T = self.unsrt[species].common_temp
				Theta = np.array([T/T,T,T**2,T**3,T**4])
				#print(f"{species} high")
				#print("\n"+str(species))
				a1_ = self.unsrt[species].a1_CL
				a2_ = self.unsrt[species].a2_CL
				#print("\n")
				#print(a2_)
				a2_dash = []
				a3_dash = []
				a3_ = self.unsrt[species].a3_CL
				Cp_T_mid = self.unsrt[species].Cp_T_mid_CL
				for i in range(len(a1_)):
					a1 = a1_[i][0]
					a2 = ((Cp_T_mid[i] - Theta.dot(self.unsrt[species].NominalParams))/abs(Theta.dot(np.asarray(np.dot(self.unsrt[species].cov,self.unsrt[species].zeta_max.x)).flatten())))
					a2_dash.append(np.array([a2]))
					a3 = np.random.uniform(a2,1,(1,1))[0]
					a3_dash.append(a3)
					
					generator.append([[a1],[a2],a3])
					
				self.unsrt[species].a2_CL = a2_dash
				self.unsrt[species].a3_CL = a3_dash
				#print(generator)
				generator = np.asarray(generator,dtype="object")
				#print(len(generator))
				data = self.unsrt[species].data
				chunk_size = 3000
				params = [generator[i:i+chunk_size] for i in range(0, len(generator), chunk_size)]
				Class_thermo_curves = []
				generator_thermo = []
				#for set_ in tqdm(params,desc="Set. progress"):
				#	print(set_)
				#	print("Inside_linear")
				data["generators_thermo"] = generator
				callWorkForce = Worker(self.allowed_count)	
				Class_thermo_curves_,f_c,generator_thermo_ = callWorkForce.do_unsrt_thermoLinear(data,len(generator))
				Class_thermo_curves.extend(Class_thermo_curves_)
				generator_thermo.extend(generator_thermo_)
				del callWorkForce#,generator_thermo_,Class_thermo_curves_
				Class_thermo_Curve_dict[species] = np.asarray(Class_thermo_curves)
				Generator_dict[species] = generator_thermo
				#print("Generating CL curves \n")
				#print(Class_thermo_curves)
				#raise AssertionError("line 307 -DM")
		
		return Class_thermo_Curve_dict,Generator_dict
		

			
	
	
	def get_thermo_Hybrid_Constant_HigherOrderCurves(self,n_a):
		"""
			This defination generates n_a numbers of class-C_H type thermo curves 
		"""
		#print("Entered CH curve generator")
		Class_thermo_Curve_dict = {}
		Generator_dict = {}
		for index,species in enumerate(self.unsrt):
			#print(species)
			species_name = self.unsrt[species].species
			generator = []
			Cp_T_mid = []
			a1_ = []
			T = self.unsrt[species].common_temp
			Theta = np.array([T/T,T,T**2,T**3,T**4])
			#if self.unsrt[species].a1 == None:
			#print(self.unsrt[species].a1_CH)
			if self.unsrt[species].a1_CH is None or np.all(self.unsrt[species].a1_CH)==None:
				#print("Doing low curves")
				species_1 = species
				for i in range(n_a):
					a1 = np.random.uniform(-1,1,(1,1))[0]
					#print(a1)
					P_C = self.unsrt[species].NominalParams + a1*np.asarray(np.dot(self.unsrt[species].cov,self.unsrt[species].zeta_max.x)).flatten()
					Cp_T_mid.append(Theta.dot(P_C))
					a1_.append(a1)
					generator.append([a1])
				for spec in self.unsrt:
					if self.unsrt[spec].species == species_name:
						self.unsrt[spec].a1_CH = a1_
						self.unsrt[spec].a2_CH = a1_
						self.unsrt[spec].a3_CH = a1_
						self.unsrt[spec].Cp_T_mid_CH = Cp_T_mid
						#print("\n index", index  )
						#print("\n generator",generator )
						#print("\n self.unsrt",self.unsrt[spec].a1 	)	
					else:
						continue
				data = self.unsrt[species].data
				zeta_A = self.unsrt[species].zeta_max.x
				ConstantCurves = []
				for gen in generator:
					ConstantCurves.append(gen*zeta_A)
				Class_thermo_Curve_dict[species] = ConstantCurves
				#print(Class_thermo_Curve_dict)
				Generator_dict[species] = generator
			else:
				species_2 = species
				a1_ = self.unsrt[species].a1_CH
				a2_ = self.unsrt[species].a2_CH
				a2_dash = []
				a3_dash = []
				a3_ = self.unsrt[species].a3_CH
				Cp_T_mid = self.unsrt[species].Cp_T_mid_CH
				for i in range(len(a1_)):
					a1 = a1_[i][0]
					#print(a1)
					a2 = ((Cp_T_mid[i] - Theta.dot(self.unsrt[species].NominalParams))/abs(Theta.dot(np.asarray(np.dot(self.unsrt[species].cov,self.unsrt[species].zeta_max.x)).flatten())))
					a2_dash.append(a2)
					a3 = np.random.uniform(a2,1,(1,1))[0]
					a3_dash.append(a3)
					generator.append([a1,a2,a3])
				self.unsrt[species].a2_CH = a2_dash
				self.unsrt[species].a3_CH = a3_dash
			#print(generator)
				generator = np.asarray(generator,dtype="object")
				data = self.unsrt[species].data
				data["generators_thermo"] = generator
				callWorkForce = Worker(self.allowed_count)	
				Class_thermo_curves,generator_thermo = callWorkForce.do_unsrt_thermoHigherOrder(data,len(generator))
				del callWorkForce
				Class_thermo_Curve_dict[species] = np.asarray(Class_thermo_curves)
				Generator_dict[species] = generator_thermo
		return Class_thermo_Curve_dict,Generator_dict	
		
		





			
	def get_thermo_Hybrid_Linear_ConstantCurves(self,n_a):
		"""
			This defination generates n_a numbers of class-L_C type thermo curves 
		"""
		Class_thermo_Curve_dict = {}
		Generator_dict = {}
		for index,species in enumerate(self.unsrt):
			#print(species)
			
			if self.unsrt[species].a1_LC is None or np.all(self.unsrt[species].a1_LC) == None:
				species_name = self.unsrt[species].species
			
				T = self.unsrt[species].common_temp
				Theta = np.array([T/T,T,T**2,T**3,T**4])
				a1_ = []
				a2_ = []
				a3_ = []
				Cp_T_mid = []
				generator = []
				for i in range(n_a):
					a1 = np.random.uniform(-1,1,(1,1))[0]
					a2 = np.random.uniform(a1,1,(1,1))[0]
					P_2 = self.unsrt[species].NominalParams + a2*np.asarray(np.dot(self.unsrt[species].cov,self.unsrt[species].zeta_max.x)).flatten()
					Cp_T_mid.append(Theta.dot(P_2))
					a3 = np.random.uniform(a2,1,(1,1))[0]
					a1_.append(a1)
					a2_.append(a2)
					a3_.append(a3)
					generator.append([a1,a2,a3])
				
				for spec in self.unsrt:
					if self.unsrt[spec].species == species_name:
						self.unsrt[spec].a1_LC = a1_
						self.unsrt[spec].a2_LC = a2_
						self.unsrt[spec].a3_LC = a3_
						self.unsrt[spec].Cp_T_mid_LC = Cp_T_mid
									
					else:
						continue
				generator = np.asarray(generator,dtype="object")
				#print(len(generator))
				data = self.unsrt[species].data
				chunk_size = 3000
				params = [generator[i:i+chunk_size] for i in range(0, len(generator), chunk_size)]
				Class_thermo_curves = []
				generator_thermo = []
				#for set_ in tqdm(params,desc="Set. progress"):
					
				data["generators_thermo"] = generator
				callWorkForce = Worker(self.allowed_count)	
				Class_thermo_curves_,f_c,generator_thermo_ = callWorkForce.do_unsrt_thermoLinear(data,len(generator))
				Class_thermo_curves.extend(Class_thermo_curves_)
				generator_thermo.extend(generator_thermo_)
				del callWorkForce#,generator_thermo_,Class_thermo_curves_
				Class_thermo_Curve_dict[species] = np.asarray(Class_thermo_curves)
				Generator_dict[species] = generator_thermo
					
			else:
				generator = []
				species_name = self.unsrt[species].species
			
				T = self.unsrt[species].common_temp
				Theta = np.array([T/T,T,T**2,T**3,T**4])
				#print(f"{species} high")
				#print("\n"+str(species))
				a1_ = self.unsrt[species].a1_LC
				Cp_T_mid = self.unsrt[species].Cp_T_mid_LC
				a1_dash = []
				for i in range(len(a1_)):
					a1 = ((Cp_T_mid[i] - Theta.dot(self.unsrt[species].NominalParams))/abs(Theta.dot(np.asarray(np.dot(self.unsrt[species].cov,self.unsrt[species].zeta_max.x)).flatten())))
					a1_dash.append(np.array([a1]))
					generator.append([a1])
				self.unsrt[species].a1_LC = a1_dash
				
				data = self.unsrt[species].data
				zeta_A = self.unsrt[species].zeta_max.x
				ConstantCurves = []
				for gen in generator:
					ConstantCurves.append(gen*zeta_A)
				Class_thermo_Curve_dict[species] = ConstantCurves
				Generator_dict[species] = generator
				
		return Class_thermo_Curve_dict,Generator_dict
			
	def get_thermo_Hybrid_Linear_HigherOrderCurves(self,n_a):
		"""
			This defination generates n_a numbers of class-L_H type thermo curves 
		"""
		Class_thermo_Curve_dict = {}
		Generator_dict = {}
		for index,species in enumerate(self.unsrt):
			#print(species)
			
			if self.unsrt[species].a1_LH is None or np.all(self.unsrt[species].a1_LH) == None:
				species_name = self.unsrt[species].species
			
				T = self.unsrt[species].common_temp
				Theta = np.array([T/T,T,T**2,T**3,T**4])
				a1_ = []
				a2_ = []
				a3_ = []
				Cp_T_mid = []
				generator = []
				for i in range(n_a):
					a1 = np.random.uniform(-1,1,(1,1))[0]
					a2 = np.random.uniform(a1,1,(1,1))[0]
					P_2 = self.unsrt[species].NominalParams + a2*np.asarray(np.dot(self.unsrt[species].cov,self.unsrt[species].zeta_max.x)).flatten()
					Cp_T_mid.append(Theta.dot(P_2))
					a3 = np.random.uniform(a2,1,(1,1))[0]
					a1_.append(a1)
					a2_.append(a2)
					a3_.append(a3)
					generator.append([a1,a2,a3])
				
				for spec in self.unsrt:
					if self.unsrt[spec].species == species_name:
						self.unsrt[spec].a1_LH = a1_
						self.unsrt[spec].a2_LH = a2_
						self.unsrt[spec].a3_LH = a3_
						self.unsrt[spec].Cp_T_mid_LH = Cp_T_mid	
					else:
						continue
				generator = np.asarray(generator,dtype="object")
				#print(len(generator))
				data = self.unsrt[species].data
				chunk_size = 3000
				params = [generator[i:i+chunk_size] for i in range(0, len(generator), chunk_size)]
				Class_thermo_curves = []
				generator_thermo = []
				#for set_ in tqdm(params,desc="Set. progress"):
					
				data["generators_thermo"] = generator
				callWorkForce = Worker(self.allowed_count)	
				Class_thermo_curves_,f_c,generator_thermo_ = callWorkForce.do_unsrt_thermoLinear(data,len(generator))
				Class_thermo_curves.extend(Class_thermo_curves_)
				generator_thermo.extend(generator_thermo_)
				del callWorkForce#,generator_thermo_,Class_thermo_curves_
				Class_thermo_Curve_dict[species] = np.asarray(Class_thermo_curves)
				Generator_dict[species] = generator_thermo
					
			else:
				generator = []
				a1_ = self.unsrt[species].a1_LH
				a2_ = self.unsrt[species].a2_LH
				a2_dash = []
				a3_dash = []
				a3_ = self.unsrt[species].a3_H
				Cp_T_mid = self.unsrt[species].Cp_T_mid_LH
				for i in range(len(a1_)):
					a1 = a1_[i][0]
					a2 = ((Cp_T_mid[i] - Theta.dot(self.unsrt[species].NominalParams))/abs(Theta.dot(np.asarray(np.dot(self.unsrt[species].cov,self.unsrt[species].zeta_max.x)).flatten())))
					a2_dash.append(a2)
					a3 = np.random.uniform(a2,1,(1,1))[0]
					a3_dash.append(a3)
					generator.append([a1,a2,a3])
				self.unsrt[species].a2_LH = a2_dash
				self.unsrt[species].a3_LH = a3_dash
			#print(generator)
				generator = np.asarray(generator,dtype="object")
				data = self.unsrt[species].data
				data["generators_thermo"] = generator
				callWorkForce = Worker(self.allowed_count)	
				Class_thermo_curves,generator_thermo = callWorkForce.do_unsrt_thermoHigherOrder(data,len(generator))
				del callWorkForce
				Class_thermo_Curve_dict[species] = np.asarray(Class_thermo_curves)
				Generator_dict[species] = generator_thermo	
		return Class_thermo_Curve_dict,Generator_dict
	
	def get_thermo_Hybrid_HigherOrder_ConstantCurves(self,n_a):
		"""
			This defination generates n_a numbers of class-H_C type thermo curves 
		"""
		Class_thermo_Curve_dict = {}
		Generator_dict = {}
		for index,species in enumerate(self.unsrt):
			#print(species)
			species_name = self.unsrt[species].species
			generator = []
			T = self.unsrt[species].common_temp
			Theta = np.array([T/T,T,T**2,T**3,T**4])
			if self.unsrt[species].a1_HC is None or np.all(self.unsrt[species].a1_HC)==None:
				a1_ = []
				a2_ = []
				a3_ = []
				Cp_T_mid = []
				for i in range(n_a):
					a1 = np.random.uniform(-1,1,(1,1))[0]
					a2 = np.random.uniform(a1,1,(1,1))[0]
					#a1 = -0.7
					#a2 = 0.1
					P_2 = self.unsrt[species].NominalParams + a2*np.asarray(np.dot(self.unsrt[species].cov,self.unsrt[species].zeta_max.x)).flatten()
					Cp_T_mid.append(Theta.dot(P_2))
					#print(a2)
					#print(Theta.dot(P_2))
					a3 = np.random.uniform(a2,1,(1,1))[0]
					a1_.append(a1)
					a2_.append(a2)
					a3_.append(a3)
					generator.append([[a1],[a2],a3])
				for spec in self.unsrt:
					if self.unsrt[spec].species == species_name:
						self.unsrt[spec].a1_HC = a1_
						self.unsrt[spec].a2_HC = a2_
						self.unsrt[spec].a3_HC = a3_
						self.unsrt[spec].Cp_T_mid_HC = Cp_T_mid
									
					else:
						continue
				generator = np.asarray(generator,dtype="object")
				data = self.unsrt[species].data
				data["generators_thermo"] = generator
				callWorkForce = Worker(self.allowed_count)	
				Class_thermo_curves,generator_thermo = callWorkForce.do_unsrt_thermoHigherOrder(data,len(generator))
				del callWorkForce
				#print(np.asarray(Class_thermo_curves))
				#raise AssertionError("Stop")
				Class_thermo_Curve_dict[species] = np.asarray(Class_thermo_curves)
				Generator_dict[species] = generator_thermo	
					
			else:
				generator = []
				species_name = self.unsrt[species].species
			
				T = self.unsrt[species].common_temp
				Theta = np.array([T/T,T,T**2,T**3,T**4])
				#print(f"{species} high")
				#print("\n"+str(species))
				a1_ = self.unsrt[species].a1_HC
				Cp_T_mid = self.unsrt[species].Cp_T_mid_HC
				a1_dash = []
				for i in range(len(a1_)):
					a1 = ((Cp_T_mid[i] - Theta.dot(self.unsrt[species].NominalParams))/abs(Theta.dot(np.asarray(np.dot(self.unsrt[species].cov,self.unsrt[species].zeta_max.x)).flatten())))
					a1_dash.append(np.array([a1]))
					generator.append([a1])
				self.unsrt[species].a1_C = a1_dash
				
				data = self.unsrt[species].data
				zeta_A = self.unsrt[species].zeta_max.x
				ConstantCurves = []
				for gen in generator:
					ConstantCurves.append(np.asarray(gen*zeta_A).flatten())
				Class_thermo_Curve_dict[species] = ConstantCurves
				Generator_dict[species] = generator
				
		return Class_thermo_Curve_dict,Generator_dict
		
	def get_thermo_Hybrid_HigherOrder_LinearCurves(self,n_a):
		"""
			This defination generates n_a numbers of class-L_H type thermo curves 
		"""
		Class_thermo_Curve_dict = {}
		Generator_dict = {}
		for index,species in enumerate(self.unsrt):
			#print(species)
			species_name = self.unsrt[species].species
			generator = []
			T = self.unsrt[species].common_temp
			Theta = np.array([T/T,T,T**2,T**3,T**4])
			if self.unsrt[species].a1_HL is None or np.all(self.unsrt[species].a1_HL)==None:
				a1_ = []
				a2_ = []
				a3_ = []
				Cp_T_mid = []
				for i in range(n_a):
					a1 = np.random.uniform(-1,1,(1,1))[0]
					a2 = np.random.uniform(a1,1,(1,1))[0]
					#a1 = -0.7
					#a2 = 0.1
					P_2 = self.unsrt[species].NominalParams + a2*np.asarray(np.dot(self.unsrt[species].cov,self.unsrt[species].zeta_max.x)).flatten()
					Cp_T_mid.append(Theta.dot(P_2))
					#print(a2)
					#print(Theta.dot(P_2))
					a3 = np.random.uniform(a2,1,(1,1))[0]
					a1_.append(a1)
					a2_.append(a2)
					a3_.append(a3)
					generator.append([[a1],[a2],a3])
				for spec in self.unsrt:
					if self.unsrt[spec].species == species_name:
						self.unsrt[spec].a1_HL = a1_
						self.unsrt[spec].a2_HL = a2_
						self.unsrt[spec].a3_HL = a3_
						self.unsrt[spec].Cp_T_mid_HL = Cp_T_mid
					else:
						continue
				generator = np.asarray(generator,dtype="object")
				data = self.unsrt[species].data
				data["generators_thermo"] = generator
				callWorkForce = Worker(self.allowed_count)	
				Class_thermo_curves,generator_thermo = callWorkForce.do_unsrt_thermoHigherOrder(data,len(generator))
				del callWorkForce
				Class_thermo_Curve_dict[species] = np.asarray(Class_thermo_curves)
				Generator_dict[species] = generator_thermo	
					
			else:
				generator = []
				species_name = self.unsrt[species].species
			
				T = self.unsrt[species].common_temp
				Theta = np.array([T/T,T,T**2,T**3,T**4])
				#print(f"{species} high")
				#print("\n"+str(species))
				a1_ = self.unsrt[species].a1_HL
				a2_ = self.unsrt[species].a2_HL
				#print("\n")
				#print(a2_)
				a2_dash = []
				a3_dash = []
				a3_ = self.unsrt[species].a3_L
				Cp_T_mid = self.unsrt[species].Cp_T_mid_HL
				for i in range(len(a1_)):
					a1 = a1_[i][0]
					a2 = ((Cp_T_mid[i] - Theta.dot(self.unsrt[species].NominalParams))/abs(Theta.dot(np.asarray(np.dot(self.unsrt[species].cov,self.unsrt[species].zeta_max.x)).flatten())))
					a2_dash.append(np.array([a2]))
					a3 = np.random.uniform(a2,1,(1,1))[0]
					a3_dash.append(a3)
					
					generator.append([[a1],[a2],a3])
				self.unsrt[species].a2_HL = a2_dash
				self.unsrt[species].a3_HL = a3_dash
			#print(generator)
			generator = np.asarray(generator,dtype="object")
			#print(len(generator))
			data = self.unsrt[species].data
			chunk_size = 3000
			params = [generator[i:i+chunk_size] for i in range(0, len(generator), chunk_size)]
			Class_thermo_curves = []
			generator_thermo = []
			#for set_ in tqdm(params,desc="Set. progress"):
				
			data["generators_thermo"] = generator
			callWorkForce = Worker(self.allowed_count)	
			Class_thermo_curves_,f_c,generator_thermo_ = callWorkForce.do_unsrt_thermoLinear(data,len(generator))
			Class_thermo_curves.extend(Class_thermo_curves_)
			generator_thermo.extend(generator_thermo_)
			del callWorkForce#,generator_thermo_,Class_thermo_curves_
			Class_thermo_Curve_dict[species] = np.asarray(Class_thermo_curves)
			Generator_dict[species] = generator_thermo
		return Class_thermo_Curve_dict,Generator_dict
							
	
	def get_thermo_higher_order_curves(self , n_a):  ###  this creats Higher_Order cp curves 
		
		Class_thermo_Curve_dict = {}
		Generator_dict = {}
		for index,species in enumerate(self.unsrt):
			#print(species)
			species_name = self.unsrt[species].species
			generator = []
			T = self.unsrt[species].common_temp
			Theta = np.array([T/T,T,T**2,T**3,T**4]).flatten()
			if self.unsrt[species].a1_H is None or np.all(self.unsrt[species].a1_H)==None:
				# this block is creating issues with list index out of 
				a1_ = []
				a2_ = []
				a3_ = []
				Cp_T_mid = []
				for i in range(n_a):
					a1 = np.random.uniform(-1,1,(1,1))[0]
					a2 = np.random.uniform(a1,1,(1,1))[0]
					#a1 = -0.7
					#a2 = 0.1
					NominalParams = np.asarray(self.unsrt[species].NominalParams)
					P_2 = NominalParams + a2*np.asarray(np.dot(self.unsrt[species].cov,self.unsrt[species].zeta_max.x)).flatten()
					#print("NominalParams shape:", NominalParams.shape)
					#print("a2 shape:", np.array(a2).shape)
					#print("cov shape:", self.unsrt[species].cov.shape)
					#print("zeta_max.x shape:", self.unsrt[species].zeta_max.x.shape)
					#print("dot product shape:", np.dot(self.unsrt[species].cov, self.unsrt[species].zeta_max.x).shape)
					#print("flattened dot product shape:", np.asarray(np.dot(self.unsrt[species].cov, self.unsrt[species].zeta_max.x)).flatten().shape)
					#raise AssertionError("line 678")
					Cp_T_mid.append(Theta.dot(P_2))
					#print(a2)
					#print(Theta.dot(P_2))
					a3 = np.random.uniform(a2,1,(1,1))[0]
					a1_.append(a1)
					a2_.append(a2)
					a3_.append(a3)
					generator.append([[a1],[a2],[a3]])
					#print(generator)
				for spec in self.unsrt:
					if self.unsrt[spec].species == species_name:
						self.unsrt[spec].a1_H = a1_
						self.unsrt[spec].a2_H = a2_
						self.unsrt[spec].a3_H = a3_
						self.unsrt[spec].Cp_T_mid_H = Cp_T_mid
							
					else:
						continue
					
			else:
				
				a1_ = self.unsrt[species].a1_H
				a2_ = self.unsrt[species].a2_H
				a2_dash = []
				a3_dash = []
				a3_ = self.unsrt[species].a3_H
				Cp_T_mid = self.unsrt[species].Cp_T_mid_H
				for i in range(len(a1_)):
					a1 = a1_[i][0]
					a2 = ((Cp_T_mid[i] - Theta.dot(self.unsrt[species].NominalParams))/abs(Theta.dot(np.asarray(np.dot(self.unsrt[species].cov,self.unsrt[species].zeta_max.x)).flatten())))
					a2_dash.append(a2)
					a3 = np.random.uniform(a2,1,(1,1))[0]
					a3_dash.append(a3)
					generator.append([a1,a2,a3])
				self.unsrt[species].a2_H = a2_dash
				self.unsrt[species].a3_H = a3_dash
			#print("generator values for Higher samples",generator)
			#print("checkpoint 1")
			generator = np.asarray(generator,dtype="object")
			data = self.unsrt[species].data
			data["generators_thermo"] = generator
			callWorkForce = Worker(self.allowed_count)	
			Class_thermo_curves,generator_thermo = callWorkForce.do_unsrt_thermoHigherOrder(data,len(generator))
			del callWorkForce
			Class_thermo_Curve_dict[species] = np.asarray(Class_thermo_curves)
			Generator_dict[species] = generator_thermo	
		return Class_thermo_Curve_dict,Generator_dict
	
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
	
	def getSA_samples(self,flag='reaction',factor=np.log(2)):
		if flag == "reaction":
			design_matrix = []
			design_matrix.extend(list(float(factor)*np.eye(self.rxn_len)))
			return np.asarray(design_matrix)
		else:
			geo_cp = 2
			delta_h = 1.987 * 298
			delta_s = 1.987  # R
			block_list = [geo_cp,delta_h,delta_s]
			block_array = np.asarray(block_list*self.rxn_len).flatten()
		return np.asarray(np.diag(block_array))
	
	
	def getNominal_samples(self,flag):
		if flag == "reaction":
			design_matrix = []
			design_matrix.append(list(np.zeros(self.rxn_len)))
			return	design_matrix
		else:
			block_list = [0.0,0.0,0.0]
			block_array = np.asarray(block_list*self.rxn_len).flatten()
			return [block_array]
	
	def getSA3P_samples(self,zeta_vector,flag):
		new_zeta = []
		if flag == "A":
			for i in zeta_vector:
				new_zeta.append([abs(i[0]),0.0,0.0])
		elif flag == "n":
			for i in zeta_vector:
				new_zeta.append([0.0,abs(i[1]),0.0])
		elif flag == "Ea":
			for i in zeta_vector:
				new_zeta.append([0.0,0.0,abs(i[2])])
		else:
			raise AssertionError("Please give a correct flag for DesignMatrix SA3P samples!!")
		
		# Number of vectors
		num_vectors = len(new_zeta)

		# Size of each vector (assume all vectors are the same size)
		vector_size = len(new_zeta[0])

		# Create an empty matrix of the appropriate size (3x9 in this case)
		design_matrix = np.zeros((num_vectors, num_vectors * vector_size))

		# Place each vector on the diagonal
		for i in range(num_vectors):
			# Insert the vector into the appropriate diagonal position
			design_matrix[i, i*vector_size:(i+1)*vector_size] = new_zeta[i]
			
		return np.asarray(design_matrix)
	
	def get_thermo_samples(self):
		print("\nStarting to generate design matrix!!\n")
		if self.design == 'A-facto': 
			#self.sim = 1 + 2* self.n+self.n*(self.n-1)/2
			design_matrix = list(2* np.random.random_sample((self.sim,self.n)) - 1)
			design_matrix.extend(list(np.eye(self.n)))
			return np.asarray(design_matrix)
		elif self.design == "Tcube":
			tic = time.time()
			design_matrix = []
			UNSHUFFLED_DM = int(self.sim*0.2) 
			print("Unshuflled Dm", UNSHUFFLED_DM)
			if int(self.sim) > 3:            			
				n_C = int(UNSHUFFLED_DM*0.2)   			# trial-1  .1 .1 .2 .1 .1 .1 .1 .1 .1	
				n_L = int(UNSHUFFLED_DM*0.2)	 		#trial-2    .1111 .1111 .1111 .1111 .1111 .1111 .1111 .1111 .1111    # taguchi trial 1
				n_H = int(UNSHUFFLED_DM*0.1)  			# trial-3  0.0667 0.0667 0.0667 0.1333 0.1333 0.1333 0.1333 0.1333 0.1333	 # taguchi trial 2
				n_CL=int(UNSHUFFLED_DM*0.1) 	#		#trial-4    0.1 0.2 .1 .1 .1 .1 .1 .1 .1 .1
				n_CH=int(UNSHUFFLED_DM*0.1)		 	# TRIAL-5 .2 .1 .1 .1 .1 .1 .1 .1 .1 .1 .1 .1
				n_LC=int(UNSHUFFLED_DM*0.1)		# trial6 - .2 .2 .1 .1 .1 .1 .1 .05 .05
				n_LH=int(UNSHUFFLED_DM*0.1)		
				n_HC=int(UNSHUFFLED_DM*0.05)	
				n_HL=int(UNSHUFFLED_DM*0.05)		
			else:
				n_C = 50
				n_L = 50
				n_H = 50
				n_CL= 50
				n_CH= 50           ######### curves are not geiing generated..may be due to strict conditions
				n_LC= 50
				n_LH= 50
				n_HC= 50
				n_HL= 50
				#UNSHUFFLED_DM = 5
				#SHUFFLED_DM = 4*UNSHUFFLED_DM
				#self.sim =  UNSHUFFLED_DM + SHUFFLED_DM
			UNSHUFFLED_DM = n_C+ n_L + n_H + n_CL + n_CH + n_LC + n_LH + n_HC + n_HL 
			SHUFFLED_DM = 4*UNSHUFFLED_DM
			self.sim =  UNSHUFFLED_DM + SHUFFLED_DM
			print(self.sim,n_C,n_L,n_H,SHUFFLED_DM)
			
			
			if "thermo_ConstantSamples.pkl" not in os.listdir():
				curves_dict_C, generator_C = self.get_thermo_ConstantCurves(n_C)# Returns 100 class-A Arrhenius samples
				with open('thermo_ConstantSamples.pkl', 'wb') as file_:
					pickle.dump(curves_dict_C, file_)
				print("\nThermo constant curves generated")
			else:
				# Load the object from the file
				with open('thermo_ConstantSamples.pkl', 'rb') as file_:
					curves_dict_C = pickle.load(file_)
				print("\nThermo constant curves generated")
					
			
			if "thermo_LinearSamples.pkl" not in os.listdir():
				curves_dict_L, generator_L = self.get_thermo_LinearCurves(n_L)# Returns 100 class-A Arrhenius samples
				with open('thermo_LinearSamples.pkl', 'wb') as file_:
					pickle.dump(curves_dict_L, file_)
				print("\nThermo linear curves generated")
			else:
				# Load the object from the file
				with open('thermo_LinearSamples.pkl', 'rb') as file_:
					curves_dict_L = pickle.load(file_)
				print("\nThermo linear curves generated")
					
			
			if "thermo_HigherOrderSamples.pkl" not in os.listdir():
				curves_dict_H, generator_H = self.get_thermo_higher_order_curves(n_H)# Returns 100 class-A Arrhenius samples
				with open('thermo_HigherOrderSamples.pkl', 'wb') as file_:
					pickle.dump(curves_dict_H, file_)
				print("\nThermo Higher curves generated")
			else:
				# Load the object from the file
				with open('thermo_HigherOrderSamples.pkl', 'rb') as file_:
					curves_dict_H = pickle.load(file_)
				print("\nThermo Higher curves generated")
					
						
			if "thermo_Hybrid_Constant_Linear_Samples.pkl" not in os.listdir():			
				curves_dict_CL, generator_CL = self.get_thermo_Hybrid_Constant_LinearCurves(n_CL)# Returns 100 class-A Arrhenius samples
				with open('thermo_Hybrid_Constant_Linear_Samples.pkl', 'wb') as file_:
					pickle.dump(curves_dict_CL, file_)
				print("\nThermo hybrid C_L curves generated")
			else:
				# Load the object from the file
				with open('thermo_Hybrid_Constant_Linear_Samples.pkl', 'rb') as file_:
					curves_dict_CL = pickle.load(file_)
				print("\nThermo hybrid C_L curves generated")
				
						
			if "thermo_Hybrid_Constant_HigherOrder_Samples.pkl" not in os.listdir():
				curves_dict_CH, generator_CH = self.get_thermo_Hybrid_Constant_HigherOrderCurves(n_CH)
				#print(curves_dict_CH)
				with open('thermo_Hybrid_Constant_HigherOrder_Samples.pkl', 'wb') as file_:
					pickle.dump(curves_dict_CH, file_)
				print("\nThermo C_H curves generated")
			else:
				# Load the object from the file
				with open('thermo_Hybrid_Constant_HigherOrder_Samples.pkl', 'rb') as file_:
					curves_dict_CH = pickle.load(file_)
				print("\nThermo C_H curves generated")
					
					
			if "thermo_Hybrid_Linear_Constant_Samples.pkl" not in os.listdir():
				curves_dict_LC, generator_LC = self.get_thermo_Hybrid_Linear_ConstantCurves(n_LC)
				with open('thermo_Hybrid_Linear_Constant_Samples.pkl', 'wb') as file_:
					pickle.dump(curves_dict_LC, file_)
				print("\nThermo L_C curves generated")
			else:
				# Load the object from the file
				with open('thermo_Hybrid_Linear_Constant_Samples.pkl', 'rb') as file_:
					curves_dict_LC = pickle.load(file_)
				print("\nThermo L_C curves generated")

									
			if "thermo_Hybrid_Linear_HigherOrder_Samples.pkl" not in os.listdir():
				curves_dict_LH, generator_LH = self.get_thermo_Hybrid_Linear_HigherOrderCurves(n_LH)
				with open('thermo_Hybrid_Linear_HigherOrder_Samples.pkl', 'wb') as file_:
					pickle.dump(curves_dict_LH, file_)
				print("\nThermo hybrid L_H curves generated")
			else:
				# Load the object from the file
				with open('thermo_Hybrid_Linear_HigherOrder_Samples.pkl', 'rb') as file_:
					curves_dict_LH = pickle.load(file_)
				print("\nThermo hybrid L_H curves generated")
				
								
			if "thermo_Hybrid_HigherOrder_Constant_Samples.pkl" not in os.listdir():
				curves_dict_HC, generator_HC = self.get_thermo_Hybrid_HigherOrder_ConstantCurves(n_HC)
				with open('thermo_Hybrid_HigherOrder_Constant_Samples.pkl', 'wb') as file_:
					pickle.dump(curves_dict_HC, file_)
				print("\nThermo hybrid H_C curves generated")
			else:
				# Load the object from the file
				with open('thermo_Hybrid_HigherOrder_Constant_Samples.pkl', 'rb') as file_:
					curves_dict_HC = pickle.load(file_)
				print("\nThermo hybrid H_C curves generated")	
					
							
			if "thermo_Hybrid_HigherOrder_Linear_Samples.pkl" not in os.listdir():
				curves_dict_HL, generator_HL = self.get_thermo_Hybrid_HigherOrder_LinearCurves(n_HL)
				with open('thermo_Hybrid_HigherOrder_Linear_Samples.pkl', 'wb') as file_:
					pickle.dump(curves_dict_HL, file_)
				print("\nThermo hybrid H_L curves generated")
			else:
				# Load the object from the file
				with open('thermo_Hybrid_HigherOrder_Linear_Samples.pkl', 'rb') as file_:
					curves_dict_HL = pickle.load(file_)
				print("\nThermo hybrid H_L curves generated")		
			
			


			
			V_ = {}
			pairings = []

			for species in tqdm(self.unsrt, desc="Populating V_"):
				species_name = self.unsrt[species].species  # Name of the species
				if not self.unsrt[species].grouped:  # Only process ungrouped species
					pair = []

					# Find all species with the same name
					for spec in self.unsrt:
						if self.unsrt[spec].species == species_name:
							self.unsrt[spec].grouped = True  # Mark species as grouped
							pair.append(spec)
					#print(pair)
					# Ensure two entries exist for low and high T polynomials
					if len(pair) != 2:
						print(f"Warning: Incomplete pairing for species {species_name}")
						continue
					# Extract curves for the two paired species
					#print(curves_dict_C[pair[0]])
					list_c1 = curves_dict_C[pair[0]]  # First pair (NASA coefficients at condition C)
					list_c2 = curves_dict_C[pair[1]]  # Second pair
					#print(species,list_c1[0],list_c2[0])
					#print(len(list_c1))
					list_l1 = curves_dict_L[pair[0]]  
					list_l2 = curves_dict_L[pair[1]]  
					#print(species,list_l1[0],list_l2[0])
					list_h1 = curves_dict_H[pair[0]]  
					list_h2 = curves_dict_H[pair[1]]  
					#print(species,list_h1[0],list_h2[0])
					list_CL1 = curves_dict_CL[pair[0]]  
					list_CL2 = curves_dict_CL[pair[1]]  
					#print(species,list_CL1[0],list_CL2[0])
					list_CH1 = curves_dict_CH[pair[0]]  
					list_CH2 = curves_dict_CH[pair[1]]  
					#print(species,list_CH1[0],list_CH2[0])
					list_LC1 = curves_dict_LC[pair[0]]  
					list_LC2 = curves_dict_LC[pair[1]] 
					#print(species,list_LC1[0],list_LC2[0])
					list_LH1 = curves_dict_LH[pair[0]]  
					list_LH2 = curves_dict_LH[pair[1]]  
					#print(species,list_LH1[0],list_LH2[0])
					list_HC1 = curves_dict_HC[pair[0]]  
					list_HC2 = curves_dict_HC[pair[1]]  
					#print(species,list_HC1[0])#,list_HC2[0])
					list_HL1 = curves_dict_HL[pair[0]]  
					list_HL2 = curves_dict_HL[pair[1]]  
					#print(species,list_HL1[0],list_HL2[0])
					#raise AssertionError("Stop")
					# Combine low-temperature (L), high-temperature (H), and other coefficients (C)
					result = []
					result.extend([np.asarray([np.asarray(c1).flatten(),np.asarray(c2).flatten()]).flatten() for c1, c2 in zip(list_c1, list_c2)])
					result.extend([np.asarray([np.asarray(l1).flatten(),np.asarray(l2).flatten()]).flatten() for l1, l2 in zip(list_l1, list_l2)])
					result.extend([np.asarray([np.asarray(h1).flatten(),np.asarray(h2).flatten()]).flatten() for h1, h2 in zip(list_h1, list_h2)])
					result.extend([np.asarray([np.asarray(CL1).flatten(),np.asarray(CL2).flatten()]).flatten() for CL1, CL2 in zip(list_CL1, list_CL2)])
					result.extend([np.asarray([np.asarray(CH1).flatten(),np.asarray(CH2).flatten()]).flatten() for CH1, CH2 in zip(list_CH1, list_CH2)])
					result.extend([np.asarray([np.asarray(LC1).flatten(),np.asarray(LC2).flatten()]).flatten() for LC1, LC2 in zip(list_LC1, list_LC2)])
					result.extend([np.asarray([np.asarray(LH1).flatten(),np.asarray(LH2).flatten()]).flatten() for LH1, LH2 in zip(list_LH1, list_LH2)])
					result.extend([np.asarray([np.asarray(HC1).flatten(),np.asarray(HC2).flatten()]).flatten() for HC1, HC2 in zip(list_HC1, list_HC2)])
					result.extend([np.asarray([np.asarray(HL1).flatten(),np.asarray(HL2).flatten()]).flatten() for HL1, HL2 in zip(list_HL1, list_HL2)])
									
					#print(result[0])
					# Store the result in V_
					V_[species_name] = result  # Use species name as key
					pairings.append(species_name)  # Track paired species

				else:
					continue
			
			V_s = {}#Doing random shuffling
			
			V_copy = copy.deepcopy(V_)
			
			for species in tqdm(pairings, desc="Doing random shuffling"):
				column = []
				# Number of shuffles
				num_shuffles = SHUFFLED_DM
				# Get the values to shuffle
				values = V_copy[species]
				#print(len(values))
				# Repeat the shuffle process
				for i in range(4):
					np.random.shuffle(values)
					column.extend(values.copy())

				# Flatten the list of arrays
				#print(len(column))
				#column = np.concatenate(column)
				#print(column[len(list_c1)])
				V_s[species] = column
			
			SHUFFLED_DM = 4* UNSHUFFLED_DM
			
			# STEP: Check that all species in pairings have sufficient samples
			min_entries_per_species = {}
			for species in pairings:
				actual_len = len(V_[species])
				if actual_len < UNSHUFFLED_DM:
					print(f"⚠️ Warning: Species '{species}' has only {actual_len} samples. Adjusting UNSHUFFLED_DM.")
				min_entries_per_species[species] = actual_len

			# Get the new safe UNSHUFFLED_DM
			safe_DM = min(min_entries_per_species.values())
			if safe_DM < UNSHUFFLED_DM:
				print(f"\n✅ Adjusting UNSHUFFLED_DM from {UNSHUFFLED_DM} to {safe_DM} to ensure consistency.\n")
				UNSHUFFLED_DM = safe_DM
				SHUFFLED_DM = 4 * UNSHUFFLED_DM
				self.sim = UNSHUFFLED_DM + SHUFFLED_DM  # Adjust sim value too

			print('✅ Final UNSHUFFLED_DM:', UNSHUFFLED_DM)
			print('✅ Final SHUFFLED_DM  :', SHUFFLED_DM)

			# STEP: Now safely construct the design matrix
			for i in range(UNSHUFFLED_DM):
				row = []
				for species in pairings:
					try:
						row.extend(V_[species][i])
					except IndexError:
						raise IndexError(f"❌ Index {i} out of range for species '{species}' with only {len(V_[species])} samples")
				design_matrix.append(row)
			
			#RANDOM SHUFFLING
			for i in range(SHUFFLED_DM):
				row = []
				for species in pairings:
					row.extend(V_s[species][i])
				design_matrix.append(row)
			
			
			tok = time.time()
			print("Time taken to construct Design Matrix: {}".format(tok - tic))
			
			design_matrix = np.asarray(design_matrix)
			s = ""
			for row in design_matrix:
				
				row_vals = []
				for element in row:
				    if isinstance(element, (np.ndarray, list)):
				        row_vals.extend(np.ravel(element))
				    else:
				        row_vals.append(element)
				s += ",".join(f"{float(val):.16e}" for val in row_vals) + "\n"
				

			with FileLock("DesignMatrix.csv.lock"):
			
				with open('DesignMatrix.csv', 'w') as ff:
					ff.write(s)
			return np.asarray(design_matrix) 
			
			



	
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
				while len(temp)<100:
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
				
			
			for i in range(100):
				row = []
				for rxn in self.unsrt:
					row.extend(V[rxn][i])
					
				design_matrix.append(row)
			tok = time.time()
			print("Time taken to construct Design Matrix: {}".format(tok - tic))
			#s =""
			#for row in design_matrix:
			#	for element in row:
			#		s+=f"{element},"
			#	s+="\n"
			#ff = open('DesignMatrix.csv','w').write(s)	
			
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
	
def run_sampling_thermo_HigherOrder(sample,data,generator,length,row):
	A = Uncertainty.ThermoUncertaintyExtractor(data)
	a1 = generator[0]
	a2 = generator[1]
	a3 = generator[2]
	A.populateValues(a1,a2,a3)
	A.getCovariance(flag=False)
	A.getUnCorrelated(flag = False)
	#print("Uncorrelated done")
	zeta = A.get_thermo_ZetaHigherOrder(flag=True)
	del A
	return (sample,row,generator,zeta,length)

def run_sampling_thermo_LinearCurce(sample, data, generator, length, row):
	A = Uncertainty.ThermoUncertaintyExtractor(data)
	a1 = generator[0]
	a2 = generator[1]
	a3 = generator[2]
	#print(generator)
	#print("Inside thermo linear curve")
	A.populateValues(a1, a2, a3)
	A.getCovariance(flag=False)
	A.getUnCorrelated(flag=False)
	zeta, Yu = A.get_thermo_ZetaLinear(row, flag=True)
	del A
	return (sample,row, generator, zeta, Yu, length) 

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
		self.results= []
		self.parallized_zeta = []
		self.generator = []
		self.f_c = []
		self.parallel_zeta_dict = {}

	def callback(self,result):
		#print("Entered callback\n")
		self.progress.append(result[0])
		self.results.append(result)
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
	
	def do_unsrt_thermoHigherOrder(self,data,sampling_points):
		#generators = data["generators"]
		for args in range(sampling_points):
			self.pool.apply_async(run_sampling_thermo_HigherOrder, 
				  args=(1,data,data["generators_thermo"][args],sampling_points,args,), 
				  callback=self.callback,error_callback=self.callback_error)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
		
		
		self.results.sort(key=lambda x: x[1])  # Sort by the second-to-last element (row identifier)
		
		# Extract sorted values
		sorted_generators = [result[2] for result in self.results]
		sorted_zeta = [result[3] for result in self.results]
		return sorted_zeta, sorted_generators
	
	def do_unsrt_thermoLinear(self, data, sampling_points):
		for args in range(sampling_points):
			self.pool.apply_async(
				run_sampling_thermo_LinearCurce, 
				args=(1, data, data["generators_thermo"][args], sampling_points, args,), 
				callback=self.callback, 
				error_callback=self.callback_error
			)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()

		# Sort results based on the identifier (row index)
		self.results.sort(key=lambda x: x[1])  # Sort by the second-to-last element (row identifier)
		
		# Extract sorted values
		sorted_generators = [result[2] for result in self.results]
		sorted_zeta = [result[3] for result in self.results]
		sorted_fc = [result[4] for result in self.results]
		
		return sorted_zeta, sorted_fc, sorted_generators
		#return self.parallized_zeta,self.f_c,self.generator
	
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


