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
import copy
import yaml
### Program specific modules
import make_input_file 
from MechManipulator2_0 import Manipulator as manipulator
import data_management
from MechanismParser import Parser

############################################
##  This module is to create directroies  ##
##  is systematic ways. Prallel computing ##
##  is done to make the process fast      ##
############################################
def run_generate_dir(location):
	os.mkdir(location)
	os.mkdir(location+"/output")

def run_map(params):
	location = str(params[3])
	yaml_string = yaml.dump(params[0],default_flow_style=False)
	yamlfile = open(location+"/mechanism.yaml","w").write(yaml_string)
	sim = open(location+"/cantera_.py",'w').write(params[1])
	sim = open(location+"/FlameMaster.input",'w').write(params[1])
	extract = open(location+"/extract.py",'w').write(params[4])
	#runConvertorScript = open(location+"/run_convertor",'w').write(params[2])
	runScript = open(location+"/run","w").write(params[2])
	#subprocess.call(["chmod","+x",location+"/run_convertor"])
	subprocess.call(["chmod","+x",location+"/run"])
	del yaml_string
	del location
	return (params[3])
	
def run_executable_files(args):
	os.chdir(args[0])
	subprocess.call(["./"+args[1]])
	return (args[0],args[2])

class Worker():
	def __init__(self, workers):
		self.pool = multiprocessing.Pool(processes=workers)
		self.progress = []
		
	def callback_run(self, result):
		#print(result)
		string=""
		for i in result:
			self.progress.append(i[0])
			string+=f"{i[0]}/run\n"
			total = i[1]
		#string+="\n"
		sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(self.progress)/float(total)*100))
		sys.stdout.flush()
		f = open('../progress','+a')
		f.write(string)
		f.close()
		
	def callback_create(self, result):
		self.progress.append(result[0])
		sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(self.progress)/float(result[-1])*100))
		sys.stdout.flush()
		#self.pool.terminate()
	
	def custom_error_callback(self,error):
   	 	print(f'Got an error: {error}')
   
	def callback_error(self,result):
		print('error', result)
   
	def do_job_async(self,args):
		self.pool.map_async(run_executable_files,args,callback=self.callback_run,error_callback=self.callback_error)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
	
	def do_job_map(self,locations):
		self.pool.map_async(run_generate_dir,locations)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
	
	def do_job_map_create(self,params):
		self.pool.map_async(run_map,params,error_callback=self.custom_error_callback)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
	
	def do_job_executor(self,locations):
		with concurrent.futures.ProcessPoolExecutor() as executor:
			executor.map(run_generate_dir,locations)


class SM(object):
	def __init__(self,target_list,target_data,unsrt_data,design_matrix):
		# Need to find the active parameters, perturbation of those parameters in parallel
		"""
		
		target_list   = list of combustion_target_class object 
			        for each targets
		target_data   = dictonary of the user provided information
			        in target.opt file
 		unsrt_data    = uncertainty object of all the provided 
 			        reactions
		Required Inputs:
		  
		"""
		

		self.target_list = target_list
		self.target_data = target_data
		self.unsrt = unsrt_data
		self.design_matrix = design_matrix
		self.case_dir = range(0,len(target_list))
		self.file_type = target_data["Inputs"]["fileType"]
		self.order = target_data["Stats"]["Order_of_PRS"]
		self.parallel_threads = target_data["Counts"]["parallel_threads"]
		self.responseOrder = target_data["Stats"]["Order_of_PRS"]
		self.mech_loc = target_data["Locations"]["mechanism"]
		Prior_Mechanism = Parser(self.mech_loc).mech
		self.copy_of_mech = copy.deepcopy(Prior_Mechanism)
		
	
	
	def getTotalUnknowns(self):
		if self.order == 2:
			self.n_ = 1 + 2*self.n + (self.n*(self.n-1))/2
		#Third order
		if self.order == 3:
			self.n_ = 1 + 3*self.n + (self.n*(self.n-1))/2+(self.n*(self.n-1)*(self.n-2))/6+(self.n*(self.n-1))

		#Fourth order
		if self.order == 4:		
			self.n_ = 1 + 4*self.n + (self.n*(self.n-1))/2+(self.n*(self.n-1)*(self.n-2))/6+(self.n*(self.n-1)*(self.n-2)*(self.n-3))/24+3*(self.n*(self.n-1))+ (self.n*(self.n-1)*(self.n-2))	
		#Fifth order		
		if self.order == 5:
			self.n_ = 1 + 5*self.n + (self.n*(self.n-1))/2+(self.n*(self.n-1)*(self.n-2))/6+(self.n*(self.n-1)*(self.n-2)*(self.n-3))/24+(self.n*(self.n-1)*(self.n-2)*(self.n-3)*(self.n-4))/120+5*(self.n*(self.n-1))+3*(self.n*(self.n-1)*(self.n-2))
		return self.n_
	
	
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
		run_convert_dict = {}
		dir_list = []
		run_list = {}
		sim_dict = {}
		#print(len(self.design_matrix))
		for i in range(len(self.design_matrix)):
			self.beta_ = self.design_matrix[i]
			#print(self.beta_)
			if os.path.isdir(os.getcwd()+"/"+str(i)) == True and os.path.isdir(os.getcwd()+"/"+str(i)+"/output") != True:
				shutil.rmtree(str(i))
				yaml_dict[str(i)],sim_dict[str(i)] = manipulator(self.copy_of_mech,self.unsrt,self.beta_).doPerturbation()
				#print("Generating the files\n")
				#print(yaml_dict[str(i)]["reactions"])
				#print("-----------------------------------\n\n")
				
				instring_dict[str(i)],s_convert_dict[str(i)],s_run_dict[str(i)],extract[str(i)] = make_input_file.create_input_file(case,self.target_data,self.target_list[case]) #generate input file
				dir_list.append(str(i))
				run_convert_dict[str(i)] = start+"/"+str(i)
				run_list[str(i)] = start+"/"+str(i)
				self.dir_list.append(start+"/"+str(i)+"/run")
			elif os.path.isdir(os.getcwd()+"/"+str(i)) != True:
				yaml_dict[str(i)],sim_dict[str(i)] = manipulator(self.copy_of_mech,self.unsrt,self.beta_).doPerturbation()
				instring_dict[str(i)],s_convert_dict[str(i)],s_run_dict[str(i)],extract[str(i)] = make_input_file.create_input_file(case,self.target_data,self.target_list[case]) #generate input file
				dir_list.append(str(i))
				run_convert_dict[str(i)] = start+"/"+str(i)
				run_list[str(i)] = start+"/"+str(i)
				self.dir_list.append(start+"/"+str(i)+"/run")
				#print("Generating the files\n")
				#print(yaml_dict[str(i)]["reactions"])
				#print("-----------------------------------\n\n")
				
			else:
				
				yaml_dict[str(i)],sim_dict[str(i)] = manipulator(self.copy_of_mech,self.unsrt,self.beta_).doPerturbation()#Makes the perturbed mechanism
				instring_dict[str(i)],s_convert_dict[str(i)],s_run_dict[str(i)],extract[str(i)] = make_input_file.create_input_file(case,self.target_data,self.target_list[case]) #generate input file
				#print("Generating the files\n")
				#print(yaml_dict[str(i)]["reactions"])
				#print("-----------------------------------\n\n")
				
				run_convert_dict[str(i)] = start+"/"+str(i)
				run_list[str(i)] = start+"/"+str(i)
				self.dir_list.append(start+"/"+str(i)+"/run")
				continue	
		#print(yaml_dict)
		return  yaml_dict,instring_dict,s_run_dict,dir_list,run_convert_dict,run_list,extract,sim_dict
		
	
	def make_dir_in_parallel(self):
		
		Prior_Mechanism = Parser(self.mech_loc).mech
		self.copy_of_mech = copy.deepcopy(Prior_Mechanism)
		
		start_time = time.time()
		print("Creating Directories for simulations......\n\n\n This may take a while... Please be patient...\n\n ")
		#Parallel_jobs = multiprocessing.Pool(self.allowed_count-1)
		self.case_manipulation = {}
		allowed_count = int(self.parallel_threads)
		self.dir_list = []
		self.generator_list = []
		self.generators = []
		optDir = os.getcwd()
		#self.n = 0
		#for i in self.unsrt:
		#	self.n+=len(self.unsrt[i].activeparameters)
		
		"""
		1] Number of required simulations
		2] Generating the sampling matrix:
			A.x = b
		
		"""		
		for case_index,case in enumerate(self.case_dir):
		
			self.beta_list = []
			yaml_dict,instring_dict,s_run_dict,dir_list,run_convert_dict,run_list,extract_list,sim_dict = self.getDirectoryList(case,case_index)
			
			
			
			###########################################
			##   Generated directories in parallel   ##
			##                                       ##
			###########################################
			start_time = time.time()
			W = Worker(allowed_count)
			W.do_job_map(dir_list)
			print("\tDirectories for case - {} is generated\n".format(case))
			del W
			
			###########################################
			##   Generated required files            ##
			##        in parallel                    ##
			###########################################
			
			
			V = Worker(allowed_count)
			yaml_list = []
			mech = []
			thermo = []
			trans = []
			instring = []
			#convertor  = []
			#convertor_2  = []
			run = []
			locations = []
			extract = []
			for i in range(len(self.design_matrix)):
				instring.append(instring_dict[str(i)])
				yaml_list.append(yaml_dict[str(i)])
				run.append(s_run_dict[str(i)])
				locations.append(run_list[str(i)])
				extract.append(extract_list[str(i)])
			#params = list(zip(mech,thermo,trans,instring,convertor,run,locations,extract))
			params = list(zip(yaml_list,instring,run,locations,extract))
			
			V.do_job_map_create(params)
			print("\tRequired files for case - {} is generated\n".format(case))
			del V
			
			
			###########################################
			##   Running the files                   ##
			##        in parallel                    ##
			###########################################
			
			X = Worker(allowed_count)
			file_n = []
			length = []
			for i in locations:
				file_n.append("run")
				length.append(len(locations))
			args = list(zip(locations,file_n,length))
			X.do_job_async(args)
			
			del X
			print("\tSimulations for case - {} is Done!!".format(case))
				
			dt = int(time.time() - start_time)
			hours = dt/3600
			minutes = (dt%3600)/60
			seconds = dt%60
			#os.system("clear")
			print("Performed {} Simulations....".format(len(locations)))
			print("Time for performing simulations : {h} hours,  {m} minutes, {s} seconds\n................................................ \n".format(h = hours, m = minutes, s =seconds))
			#print(locations)
			#del W,V,U,X
			
			#self.case_manipulation[str(case)] = sim_dict
			simulation_locations = open(optDir+"/locations",'+a')
			for loc in locations:
				simulation_locations.write(i+"\n")
			simulation_locations.close()
			
			os.chdir('..')
			
		#Parallel_jobs.close()
		#Parallel_jobs.join()
		
		#self.params_manipulation = {}
		#for j in self.selectedParams:
		#	self.case_database = {}
		#	for case in self.case_dir:
		#		temp = []
		#		for i,key in enumerate(sim_dict):
		#			temp.append(self.case_manipulation[str(case)][str(dictionary)][str(j)])
		#		self.case_database[str(case)] = temp
		#	self.params_manipulation[str(j)] = self.case_database
			
		#print(self.params_manipulation)
		
		#return self.dir_list,self.params_manipulation
	
		

