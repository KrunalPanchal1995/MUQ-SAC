#default python modules
import os, shutil, re, math 
import numpy as np
import pandas as pd
from pyDOE import *
home_dir = os.getcwd()
import multiprocessing
import subprocess
import time
import sys
### Program specific modules
import make_input_file 
import MechanismManipulator
import data_management
from MechanismParser import Parser

############################################
##  This module is to create directroies  ##
##  is systematic ways. Prallel computing ##
##  is done to make the process fast      ##
############################################

class Worker():
	def __init__(self, workers):
		self.pool = multiprocessing.Pool(processes=workers)
		self.progress = []
		
	def callback_run(self, result):
		self.progress.append(result[0])
		sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(self.progress)/float(result[-1])*100))
		sys.stdout.flush()
		f = open('../progress','a')
		f.write(result[0]+"/run"+"\n")
		f.close()
		
	def callback_runConvertor(self, result):
		self.progress.append(result[0])
		sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(self.progress)/float(result[-1])*100))
		sys.stdout.flush()
		#self.pool.terminate()
	def callback_create(self, result):
		self.progress.append(result[0])
		sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(self.progress)/float(result[-1])*100))
		sys.stdout.flush()
		#self.pool.terminate()

	def do_job_async_convertor(self,sim_list,string_dict,location,file_name):
		for args in sim_list:
			self.pool.apply_async(run_executable_files, 
				  args=(args,string_dict,location,file_name,len(sim_list)), 
				  callback=self.callback_runConvertor)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
	def do_job_async(self,sim_list,string_dict,location,file_name):
		for args in sim_list:
			self.pool.apply_async(run_executable_files, 
				  args=(args,string_dict,location,file_name,len(sim_list)), 
				  callback=self.callback_run)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
		
	def do_job_create_async(self,sim_list,
				sim_dict,mech_dict,thermo_dict,
				trans_dict,run_convert_dict,
				run_dict,locations):
		for args in sim_list:
			self.pool.apply_async(run_copy_files, 
				  args=(args,sim_dict,mech_dict,thermo_dict,trans_dict,run_convert_dict,run_dict,locations,len(sim_list)), 
				  callback=self.callback_create)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
	
	def do_job_map(self,locations):
		self.pool.map_async(run_generate_dir,locations)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
	
	def do_job_map_create(self,params):
		self.pool.map_async(run_map,params)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
	
	def do_job_executor(self,locations):
		with concurrent.futures.ProcessPoolExecutor() as executor:
			executor.map(run_generate_dir,locations)


class SM(object):
	def __init__(self,target_list,target_data,unsrt_data,):
		# Need to find the active parameters, perturbation of those parameters in parallel
		"""
		
		target_list   = list of combustion_target_class object 
			        for each targets
		target_data   = dictonary of the user provided information
			        in target.opt file
 		unsrt_data    = uncertainty object of all the provided 
 			        reactions
		Required Inputs:
		 1) Prior Mechanism
		 2) 
		"""
		

		self.target_list = target_list
		self.target_data = target_data
		self.case_dir = range(0,len(target_list))
		self.file_type = target_data["Inputs"]["fileType"]
		self.order = order
		self.parallel_threads = target_data["Counts"]["parallel_threads"]
		self.responseOrder = target_data["Stats"]["Order_of_PRS"]
		Prior_Mechanism = Parser(self.mech_loc).mech
		self.copy_of_mech = deepcopy(Prior_Mechanism)
	
	def getDirectoryList(self,case,ind):
		"""
		Creating each directories in parallel
		_____________________________________
		
		This module will create dictionary of the perturbed mechanism
		-------------------------------------------------------------
		1) The format of each mechanism is Yaml. The Yaml file contains thermo and transport data
		2) Only reactions are selected to be perturbed
		
		"""
		
		if os.path.isdir(os.getcwd()+"/case-"+str(case)) == True:
			os.chdir("case-"+str(case))
		else:
			os.mkdir("case-"+str(case))
			os.chdir("case-"+str(case))
		
		start = str(os.getcwd())
		yaml_dict = {}
		instring_dict = {}
		s_convert_dict = {}
		s_run_dict = {}
		extract = {}		
		dir_list = []
		run_list = {}
		sim_dict = {}
	
		for i in range(self.sim_):

			if os.path.isdir(os.getcwd()+"/"+str(i)) == True and os.path.isdir(os.getcwd()+"/"+str(i)+"/output") != True:
				shutil.rmtree(str(i))
				yaml_dict[str(i)],sim_dict[str(i)] = manipulator(self.copy_of_mech,self.unsrt,self.beta_[i]).doPerturbation()
				instring_dict[str(i)],s_convert_dict[str(i)],s_run_dict[str(i)],extract[str(i)] = make_input_file.create_input_file(case,self.target_data,self.target_list[case]) #generate input file
				dir_list.append(str(i))
				run_convert_dict[str(i)] = start+"/"+str(i)
				run_list[str(i)] = start+"/"+str(i)
				self.dir_list.append(start+"/"+str(i)+"/run")
			elif os.path.isdir(os.getcwd()+"/"+str(i)) != True:
				yaml_dict[str(i)],sim_dict[str(i)] = manipulator(self.copy_of_mech,self.unsrt,self.beta_[i]).doPerturbation()
				instring_dict[str(i)],s_convert_dict[str(i)],s_run_dict[str(i)],extract[str(i)] = make_input_file.create_input_file(case,self.target_data,self.target_list[case]) #generate input file
				dir_list.append(str(i))
				run_convert_dict[str(i)] = start+"/"+str(i)
				run_list[str(i)] = start+"/"+str(i)
				self.dir_list.append(start+"/"+str(i)+"/run")
			else:
				
				yaml_dict[str(i)],sim_dict[str(i)] = manipulator(self.copy_of_mech,self.unsrt,self.beta_[i]).doPerturbation()#Makes the perturbed mechanism
				instring_dict[str(i)],s_convert_dict[str(i)],s_run_dict[str(i)],extract[str(i)] = make_input_file.create_input_file(case,self.target_data,self.target_list[case]) #generate input file
				run_convert_dict[str(i)] = start+"/"+str(i)
				run_list[str(i)] = start+"/"+str(i)
				self.dir_list.append(start+"/"+str(i)+"/run")
				continue	
		return  yaml_dict,instring_dict,s_run_dict,dir_list,run_list,extract,sim_dict
		
	def make_dir_in_parallel(self):
		start_time = time.time()
		print("Creating Directories for simulations......\n\n\n This may take a while... Please be patient...\n\n ")
		#Parallel_jobs = multiprocessing.Pool(self.allowed_count-1)
		self.case_manipulation = {}
		allowed_count = int(self.allowed_count)
		self.dir_list = []
		self.generator_list = []
		self.generators = []
		"""
		1] Number of required simulations
		2] Generating the sampling matrix:
			A.x = b
		
		"""
			
		
		for case_index,case in enumerate(self.case_dir):
			"""
			Get number of samples
			"""
			
			self.getTotalUnknowns()
			
			"""
			Get generated samples
			"""
			self.getSamples()
			
			print(self.sim_)
			self.beta_list = []
			#Parallel_jobs.apply_async(self.makeDir, args = (case,len(self.case_dir)), callback = self.log_result)
			#mech_dict,thermo_dict,trans_dict,instring_dict,s_convert_dict,s_run_dict,dir_list,run_convert_dict,run_list,extract_list,sim_dict = self.getDirectoryList(case,case_index)
			
			yaml_dict,instring_dict,s_convert_dict,s_run_dict,dir_list,run_convert_dict,run_list,extract_list,sim_dict = self.getDirectoryList(case,case_index)
			#for i in yaml_dict:
			#	print(i)
			#print(yaml_dict['0']["reactions"][47])
			#print(yaml_dict['1']["reactions"][47])
			#progress,dict_ = self.makeDir(case,len(self.case_dir))
			start_time = time.time()
			#print(len(run_convert_dict))
			W = Worker(allowed_count)
			#W.do_job_executor(dir_list)
			W.do_job_map(dir_list)
			print("\tDirectories for case - {} is generated\n".format(case))
			
			del W
			
			V = Worker(allowed_count)
			yaml_list = []
			mech = []
			thermo = []
			trans = []
			instring = []
			convertor  = []
			#convertor_2  = []
			run = []
			locations = []
			extract = []
			for i in dir_list:
				instring.append(instring_dict[i])
				yaml_list.append(yaml_dict[i])
				#mech.append(mech_dict[i])
				#thermo.append(thermo_dict[i])
				#trans.append(trans_dict[i])
				convertor.append(s_convert_dict[i])
				#convertor_2.append(s_convert_dict_2[i])
				run.append(s_run_dict[i])
				locations.append(run_list[i])
				extract.append(extract_list[i])
			#params = list(zip(mech,thermo,trans,instring,convertor,run,locations,extract))
			params = list(zip(yaml_list,instring,convertor,run,locations,extract))
			
			#print(params[0][6])
			#V.do_job_create_async(dir_list,instring_dict,mech_dict,thermo_dict,trans_dict,s_convert_dict,s_run_dict,run_list)
			#if finished == False:
			#print(locations)
			#print(params)
			V.do_job_map_create(params)
			print("\tRequired files for case - {} is generated\n".format(case))
			
			del V
			
			file_n = []
			length = []
			for i in locations:
				file_n.append("run_convertor")
				length.append(len(locations))
			args = list(zip(locations,file_n,length))
			U = Worker(allowed_count)
			#U.do_job_async_convertor(locations,"run_convertor")
			U.do_job_async_convertor(args)
			
			del U
			
			#U2 = Worker(allowed_count)
			#U2.do_job_async_convertor(locations,"run_convertor_2")
			
			#del U2
			
			print("\tFiles converted to standard input files for case - {}\n".format(case))
			
			X = Worker(allowed_count)
			file_n = []
			length = []
			for i in locations:
				file_n.append("run")
				length.append(len(locations))
			args = list(zip(locations,file_n,length))
			#X.do_job_async(locations,"run")
			X.do_job_async(args)
			
			del X
			print("\tSimulations for case - {} is Done!!".format(case))
				
			dt = int(time.time() - start_time)
			hours = dt/3600
			minutes = (dt%3600)/60
			seconds = dt%60
			#os.system("clear")
			print("Performed {} Simulations....".format(self.sim_))
			print("Time for performing simulations : {h} hours,  {m} minutes, {s} seconds\n................................................ \n".format(h = hours, m = minutes, s =seconds))
			#print(locations)
			#del W,V,U,X
			
			self.case_manipulation[str(case)] = sim_dict
			s = ''
			string = ''
			string_norm = ''
			for k in self.beta_list:
				for l in k:
					s +='{},'.format(l)
				s +='\n'
			ff = open('../Data/Simulations/Beta_list_case-'+str(case)+'.csv','w').write(s)
			 
			for s in self.generators:
				for l in s:
					string +='{},'.format(l)
				string +='\n'
			ff = open('../Data/Simulations/Generator_list_case-'+str(case)+'.csv','w').write(string)
			
			for s in self.X:
				for l in s:
					string_norm +='{},'.format(l)
				string_norm +='\n'
			ff = open('../Data/Simulations/NormalizedBeta_list_case-'+str(case)+'.csv','w').write(string_norm)
			
			
			os.chdir('..')
			
		#Parallel_jobs.close()
		#Parallel_jobs.join()
		
		self.params_manipulation = {}
		#for j in self.selectedParams:
		#	self.case_database = {}
		#	for case in self.case_dir:
		#		temp = []
		#		for i,key in enumerate(sim_dict):
		#			temp.append(self.case_manipulation[str(case)][str(dictionary)][str(j)])
		#		self.case_database[str(case)] = temp
		#	self.params_manipulation[str(j)] = self.case_database
			
		#print(self.params_manipulation)
		
		return self.dir_list,self.params_manipulation
	
		

