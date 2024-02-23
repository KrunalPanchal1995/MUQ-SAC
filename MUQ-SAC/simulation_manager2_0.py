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
import concurrent.futures
import Uncertainty
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

def run_executable_files_(args):
	os.chdir(args[0])
	subprocess.call(["./"+args[1]])
	return (args[0],args[2])
def run_executable_files(location,file_name,n):
	os.chdir(location)
	subprocess.call(["./"+file_name])
	return (location,n)
def run_generate_dir(location,total):
	os.mkdir(location)
	os.mkdir(location+"/output")
	return (location,total)

def run_map_2(params):
	location = str(params[1])
	yaml_string = yaml.dump(params[0],default_flow_style=False)
	with open(location+"/mechanism.yaml","w") as yamlfile:
		yamlfile.write(yaml_string)
	return (params[1])


def run_map(params):
	location = str(params[2])
	
	sim1 = open(location+"/cantera_.py",'w').write(params[0])
	sim2= open(location+"/FlameMaster.input",'w').write(params[0])
	extract = open(location+"/extract.py",'w').write(params[3])
	#runConvertorScript = open(location+"/run_convertor",'w').write(params[2])
	runScript = open(location+"/run","w").write(params[1])
	#subprocess.call(["chmod","+x",location+"/run_convertor"])
	subprocess.call(["chmod","+x",location+"/run"])
	
	del location
	del sim1, sim2,extract,runScript
	return (params[2])

def run_direct_map(file_dict):
	location = str(file_dict["mkdir"])
	#print(location)
	betaFile = open(location+"/beta.txt",'w').write(str(file_dict["beta"]))
	yaml_string = yaml.dump(file_dict["mechanism"],default_flow_style=False)
	with open(location+"/mechanism.yaml","w") as yamlfile:
		yamlfile.write(yaml_string)
	
	sim = open(location+"/cantera_.py",'w').write(file_dict["simulationInputString"])
	sim = open(location+"/FlameMaster.input",'w').write(file_dict["simulationInputString"])
	extract = open(location+"/extract.py",'w').write(file_dict["extractString"])
	runConvertorScript = open(location+"/run_convertor",'w').write(file_dict["file_convertor_script"])
	runScript = open(location+"/run","w").write(file_dict["run_script"])
	#print(location+"/run_convertor")
	subprocess.call(["chmod","+x",location+"/run_convertor"])
	subprocess.call(["chmod","+x",location+"/run"])
	

def update_progress(progress, total):
	sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(progress) / float(total) * 100))
	sys.stdout.flush()

def run_sampling_direct(sample,rxn,data,generator,length):
	A = Uncertainty.UncertaintyExtractor(data)
	#print(generator)
	a1 = generator[0]
	a2 = generator[1]
	#print(a1,a2,a3)
	#a1 = np.random.uniform(-1,1,1)
	#a2 = np.random.uniform(-1,1,1)
	#a3 = np.random.uniform(-1,1,1)
	#A.generator.append([a1[0],a2[0],a3[0]])
	#if a1 == float(0) and a2 == float(0) and a3 == float(0):
	#	zeta = [0,0,0]
	#else:
	A.populateValues(a1,a2)
	A.getCovariance(flag=False)
	A.getUnCorrelated(flag = False)
	zeta = A.getB2Zeta(flag=True)
	#print(data)
	#zeta =UncertaintyExtractor(data).getExtremeCurves(sample)	
	del A
	#print(zeta)
	
	return (sample,rxn,generator,zeta,length)



class Worker():
	def __init__(self, workers):
		self.pool = multiprocessing.Pool(processes=workers)
		self.pool1 = concurrent.futures.ProcessPoolExecutor(max_workers=workers)
		self.progress = []
		self.parallized_zeta = []
		self.parallel_zeta_dict = {}
		self.generator = []
	def update_progress(self, total):
        	update_progress(self.progress, total)
	def callback_runConvertor(self, result):
		self.progress.append(result[0])
		sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(self.progress)/float(result[-1])*100))
		sys.stdout.flush()
	
	def callback_direct(self,result):
		#print("Entered callback\n")
		self.progress.append(result[0])
		self.generator.append(result[2])
		self.parallized_zeta.append(result[3])
		self.parallel_zeta_dict[result[1]] = result[3]
		sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(self.progress)/float(result[-1])*100))
		sys.stdout.flush()	
	def callback_run_(self, result):
		#print(result)
		self.progress.append(result[0])
		#sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(self.progress)/float(result[-1])*100))
		#sys.stdout.flush()
	def callback_run(self, result):
		#print(result)
		self.progress.append(result[0])
		sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(self.progress)/float(result[-1])*100))
		sys.stdout.flush()
		#string=""
		#for i in result:
		#	self.progress.append(i[0])
		#	string+=f"{i[0]}/run\n"
		#	total = i[1]
		#string+="\n"
		#sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(self.progress)/float(total)*100))
		#sys.stdout.flush()
		#f = open('../progress','+a')
		#f.write(string)
		#f.close()
		
	#def callback_create(self, result):
	#	self.progress.append(result[0])
	#	sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(self.progress)/float(result[-1])*100))
	#	sys.stdout.flush()
		#self.pool.terminate()
	def callback_create(self, future):
		result = future.result()
		self.progress.append(result[0])
		self.update_progress(result[-1])
	def custom_error_callback(self,error):
   	 	print(f'Got an error: {error}')
   
	def callback_error(self,result):
		print('error', result)
   
	
	def do_job_async_unsrt_direct(self,data,generator,beta):
		#print("Entered Async\n")
		#print("Entered async\n")
		for args in data:
			self.pool.apply_async(run_sampling_direct, 
				  args=(1,args,data[args],generator[args],beta,), 
				  callback=self.callback_direct,error_callback = self.callback_error)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
		return self.parallel_zeta_dict
	def do_job_async_(self,args):
		self.pool.map_async(run_executable_files_,args,callback=self.callback_run_,error_callback=self.callback_error)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
	
	def do_job_direct_map(self,params):
		self.pool.map_async(run_direct_map,params,error_callback = self.handler)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
		
	def do_job_async(self,location,file_name):
		for args in location:
			self.pool.apply_async(run_executable_files, 
				  args=(args,file_name,len(location)), 
				  callback=self.callback_run,error_callback=self.callback_error)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()	
	def do_job_map(self, locations):
		futures = [self.pool1.submit(run_generate_dir, location,len(locations)) for location in locations]
		for future in concurrent.futures.as_completed(futures):
		    self.callback_create(future)
	#def do_job_map(self,locations):
	#	self.pool.map_async(run_generate_dir,locations)
	#	self.pool.close()
	#	self.pool.join()
	#	self.pool.terminate()
	
	#def do_job_map_create(self, params):
	 #       futures = [self.pool1.submit(run_map, param,len(params)) for param in params]
	  #      for future in concurrent.futures.as_completed(futures):
	   #     	self.callback_create(future)
	
	def handler(self,error):
		print(f'Error: {error}', flush=True)
	def do_job_map_create_2(self,params):
	#	for param in params:
	#		self.pool.apply_async(run_map,args=(param,))
		self.pool.map_async(run_map_2,params,error_callback=self.custom_error_callback)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
	def do_job_async_convertor(self,location,file_name):
		for args in location:
			self.pool.apply_async(run_executable_files, 
				  args=(args,file_name,len(location)), 
				  callback=self.callback_runConvertor,error_callback=self.callback_error)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
	def do_job_map_create(self,params):
	#	for param in params:
	#		self.pool.apply_async(run_map,args=(param,))
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
		self.allowed_count = target_data["Counts"]["parallel_threads"]
		self.fuel = "NC7H16"
	
	
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
		
	def getSimulationFiles(self,case_index,cases,beta,iter_number):
		self.simulation = "Opt"
		if "DirectSimulation" not in os.listdir():
			os.mkdir("DirectSimulation")
			os.chdir("DirectSimulation")
			
			os.mkdir("Data")
			os.mkdir(str(iter_number))
			os.chdir(str(iter_number))
			start = os.getcwd()
		else:
			os.chdir(f"DirectSimulation")
			if str(iter_number) not in os.listdir():
				os.mkdir(str(iter_number))
				os.chdir(str(iter_number))
				start = os.getcwd()
			else:
				os.chdir(str(iter_number))
				start = os.getcwd()
		dictionary_list = []
		run_convert_list = []
		mkdir_list = []
		dir_run_list = []
		output_list = []
		beta_transformed = np.ones(len(beta))
		for ind,case in enumerate(case_index):
			file_dict = {}
			file_dict["beta"] = beta
			"""
			Get the zeta's using the beta
			"""
			count = 0
			generators = {}
			data = {}
			reactionList = [i for i in self.unsrt]
			#for i,rxn in enumerate(self.unsrt):
			#	data[rxn] = self.unsrt[rxn].data
			#	activeParams = self.unsrt[rxn].activeParameters
			#	gen = []
			#	for i in range(len(activeParams)):
			#		gen.append(beta[count])
			#		count+=1
			#	generators[rxn] = np.asarray(gen)	
			#X = Worker(100)
			#zeta_dict = X.do_job_async_unsrt_direct(data,generators,len(beta))
			#del X
			#beta_ = []
			#print(zeta_dict)
			#for i in self.unsrt:
			#	beta_.extend(list(zeta_dict[i]))
			beta_ = np.asarray(beta)
			
			file_dict["mechanism"],sim = manipulator(self.copy_of_mech,self.unsrt,beta_).doPerturbation()
			file_dict["simulationInputString"],file_dict["file_convertor_script"],file_dict["run_script"],file_dict["extractString"] = make_input_file.create_input_file(case,self.target_data,self.target_list[case])	
			file_dict["run_convert"] = start+"/case-"+str(case)+"/run_convertor"
			file_dict["run_list"] = start+"/case-"+str(case)+"/run"
			file_dict["mkdir"] = start+"/case-"+str(case)
			run_convert_list.append(start+"/case-"+str(case))
			dir_run_list.append(start+"/case-"+str(case))
			mkdir_list.append("case-"+str(case))
			dictionary_list.append(file_dict)
			output_list.append(start+"/case-"+str(case)+"/output/")
			#print(file_dict)
		return dictionary_list,mkdir_list,dir_run_list,run_convert_list,output_list
	def getSimulatedValues(self,beta,case_index,cases,iter_number,objective):

		Prior_Mechanism = Parser(self.mech_loc).mech
		self.copy_of_mech = copy.deepcopy(Prior_Mechanism)
		print(f"Starting direct simulations: Iteration number {iter_number}, total cases to simulate: {len(case_index)}")
		dictionary_list,mkdir_list,dir_run_list,run_convert_list,output_list = self.getSimulationFiles(case_index,cases,beta,iter_number)
		start_time = time.time()
			
		W = Worker(int(self.allowed_count))
		#W.do_job_executor(dir_list)
		W.do_job_map(mkdir_list)
		print("\tDirectories for direct simulations are generated\n")
		
		V = Worker(int(self.allowed_count))
		#print(params[0][6])
		#V.do_job_create_async(dir_list,instring_dict,mech_dict,thermo_dict,trans_dict,s_convert_dict,s_run_dict,run_list)
		V.do_job_direct_map(dictionary_list)
		print("\tRequired files for {} iteration is generated\n".format(iter_number))
		
		U = Worker(int(self.allowed_count))
		U.do_job_async_convertor(run_convert_list,"run_convertor")
		
		print("\tFiles converted to standard input files for iteration \n".format(iter_number))
		
		X = Worker(int(self.allowed_count))
		X.do_job_async(dir_run_list,"run")
		
		print("\tSimulations for {} iteration is Done!!".format(iter_number))
			
		dt = int(time.time() - start_time)
		hours = dt/3600
		minutes = (dt%3600)/60
		seconds = dt%60
		#os.system("clear")
		print("Performed {} Simulations....".format(len(dir_run_list)))
		print("Time for performing simulations : {h} hours,  {m} minutes, {s} seconds\n................................................ \n".format(h = hours, m = minutes, s =seconds))
		del W,V,U,X
		"""
		Extracting outputs for direct simulations
		"""
		
		os.chdir('..')
		directSimulation = []
		#print(case_index)
		for i,index in enumerate(case_index):
			eta_list = data_management.extract_direct_simulation_values(index,output_list[i],self.target_list,self.fuel)
			directSimulation.extend(eta_list)
		os.chdir("..")
		#sample = open("samplefile.txt","+a")
		#sample.write(f"{iter_number},{objective}\n")
		return directSimulation
	
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
			if "case-"+str(case_index) in os.listdir():
				print("Case-{index} is generated".format(index = case_index))
			else:
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
				
				
				#V = Worker(allowed_count)
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
				params = list(zip(instring,run,locations,extract))
				params2 = list(zip(yaml_list,locations))
				tic = time.time()
				#for param in params:
				#	location = str(param[3])
				#	yaml_string = yaml.dump(param[0],default_flow_style=False)
				#	with open(location+"/mechanism.yaml","w") as yamlfile:
				#		yamlfile.write(yaml_string)
				#	sim1 = open(location+"/cantera_.py",'w').write(param[1])
				#	sim2= open(location+"/FlameMaster.input",'w').write(param[1])
				#	extract = open(location+"/extract.py",'w').write(param[4])
				#	#runConvertorScript = open(location+"/run_convertor",'w').write(params[2])
				#	runScript = open(location+"/run","w").write(param[2])
				#	#subprocess.call(["chmod","+x",location+"/run_convertor"])
				#	subprocess.call(["chmod","+x",location+"/run"])
				V = Worker(allowed_count)
				V.do_job_map_create(params)
				del V	
				#chunk_size = 500
				#chunks = [params[i:i+chunk_size] for i in range(0, len(params), chunk_size)]
				#for params in chunks:
				#	V = Worker(allowed_count)
				#	V.do_job_map_create(params)
				#	del V
				tok = time.time()
				print("\tRequired files for case - {} is generated in {} hours, {} minutes, {} seconds time\n".format(case,(tok-tic)/3600,((tok-tic)%3600)/60,(tok-tic)%60))
				tic = time.time()
				chunk_size = 250
				chunks = [params2[i:i+chunk_size] for i in range(0, len(params2), chunk_size)]
				for params2 in chunks:
					W = Worker(allowed_count)
					W.do_job_map_create_2(params2)
					del W
				tok = time.time()	
				print("\tRequired files for case - {} is generated in {} hours, {} minutes, {} seconds time\n".format(case,(tok-tic)/3600,((tok-tic)%3600)/60,(tok-tic)%60))
				
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
				X.do_job_async_(args)
				
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
	
		

