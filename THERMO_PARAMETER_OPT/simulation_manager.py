import make_input_file 
from MechManipulator import Manipulator as manipulator
#import data_management
import os, shutil, re, math 
import numpy as np
from MechanismParser import Parser
import sys, os, time
from tqdm import tqdm
import multiprocessing
import subprocess
import time
import copy
import concurrent.futures
#try:
#	import ruamel_yaml as yaml
#except ImportError:
#	from ruamel import yaml
import yaml
sys.path.append('/yamlwriter.so')  # Adjust this path to the correct build directory
#print(sys.path)
import yamlwriter

sys.path.append('/parallel_yaml_writer.so')
import parallel_yaml_writer
############################################
##  This module is to create directroies  ##
##  is systematic ways. Prallel computing ##
##  is done to make the process fast	  ##
############################################
from typing import Tuple
"""
	Runs an executable file in a specified directory with a timeout. If the executable
	exceeds the specified time, it will be terminated.

	Parameters:
	- args (Tuple[str, str]): Tuple containing the directory path and executable name.
	- total (int): Total number of tasks, used for tracking progress or other purposes.
	- timeout (int): Timeout in seconds before the process is terminated. Default is 120 seconds.

	Returns:
	- Tuple[str, int]: Directory path and total, for tracking.
	"""
def run_executable_files_with_timeout(args: Tuple[str, str], total: int, timeout: int = 1200) -> Tuple[str, int]:
	
	try:
		output_log_path = os.path.join(args[0], "output.log")
		error_log_path = os.path.join(args[0], "error.log")

		with open(output_log_path, "w") as out, open(error_log_path, "w") as err:
			subprocess.run(["./" + args[1]], cwd=args[0], stdout=out, stderr=err, timeout=timeout)
		return (args[0], total)
		
	except subprocess.TimeoutExpired:
		error_log_path = os.path.join(args[0], "error.log")
		f = open(error_log_path,"w").write(f"Process {args[1]} in directory {args[0]} exceeded {timeout} seconds and was terminated.")
		return (args[0], total)

def run_executable_files_(args,total):
	os.chdir(args[0])
	#print("Simulation started!!")
	subprocess.call(["./"+args[1]])
	#print("Simulation done!!")
	return (args[0],total)
	
def run_executable_files(location,file_name,n):
	os.chdir(location)
	subprocess.call(["./"+file_name])
	return (location,n)
def run_generate_dir(location,total):
	os.mkdir(location)
	os.mkdir(location+"/output")
	return (location,total)

def run_map_3(params,total):
	#sys.path.append('/yamlwriter.so')  # Adjust this path to the correct build directory
	#print(sys.path)
	#import yamlwriter
	location = str(params[1])
	yamlwriter.dump_to_yaml(location,f"mechanism_{params[2]}.yaml",params[0])
	#yaml_string = yaml.dump(params[1],default_flow_style=False)
	#with open(location+"/mechanism.yaml","w") as yamlfile:
	#	yamlfile.write(yaml_string)
	#	yamlfile.close()
	return (params[0],total)
def run_map_2(params,total):
	#sys.path.append('/yamlwriter.so')  # Adjust this path to the correct build directory
	#print(sys.path)
	#import yamlwriter
	location = str(params[0])
	parallel_yaml_writer.dump_to_yaml(location,"mechanism.yaml",params[1])
	#yaml_string = yaml.dump(params[1],default_flow_style=False)
	#with open(location+"/mechanism.yaml","w") as yamlfile:
	#	yamlfile.write(yaml_string)
	#	yamlfile.close()
	return (params[0],total)

def run_map_opt(params,total):
	location = str(params[2])
	sim1 = open(location+"/cantera_.py",'w').write(params[0])
	sim2= open(location+"/FlameMaster.input",'w').write(params[0])
	extract = open(location+"/extract.py",'w').write(params[3])
	runScript = open(location+"/run","w").write(params[1])
	subprocess.call(["chmod","+x",location+"/run"])
	del location#,yaml_string
	del sim1, sim2,extract,runScript
	#del params[1]
	return (params[2],total)
	
def run_map(params,total):
	location = str(params[2])
	
	sim1 = open(location+"/cantera_.py",'w').write(params[0])
	sim2= open(location+"/FlameMaster.input",'w').write(params[0])
	extract = open(location+"/extract.py",'w').write(params[3])
	#runConvertorScript = open(location+"/run_convertor",'w').write(params[2])
	runScript = open(location+"/run","w").write(params[1])
	#subprocess.call(["chmod","+x",location+"/run_convertor"])
	subprocess.call(["chmod","+x",location+"/run"])
	yaml_string = yaml.dump(params[-2],default_flow_style=False)
	with open(location+"/mechanism.yaml","w") as yamlfile:
		yamlfile.write(yaml_string)

	subprocess.call(["cp",params[-1],location])
	#subprocess.call(["cp",params[-2],location])
	del location#,yaml_string
	del sim1, sim2,extract,runScript
	#del params[1]
	return (params[2],total)

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
	a1 = generator[0]
	a2 = generator[1]
	A.populateValues(a1,a2)
	A.getCovariance(flag=False)
	A.getUnCorrelated(flag = False)
	zeta = A.getB2Zeta(flag=True)
	del A
	return (sample,rxn,generator,zeta,length)

class Worker():
	def __init__(self, workers,path_to_yamlwriter="/data2/STUDY_OF_PARALLEL_COMPUTING/build/yamlwriter.so"):
		self.pool = multiprocessing.Pool(processes=workers)
		self.pool1 = concurrent.futures.ProcessPoolExecutor(max_workers=workers)
		self.progress = []
		self.parallized_zeta = []
		self.parallel_zeta_dict = {}
		self.generator = []
		sys.path.append(path_to_yamlwriter)  
		#self.yamlwriter = import_module('yamlwriter')
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
		self.progress.append(result[0])
		sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(self.progress)/float(result[-1])*100))
		sys.stdout.flush()
	def callback_run(self, result):
		self.progress.append(result[0])
		sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(self.progress)/float(result[-1])*100))
		sys.stdout.flush()
		f = open('../progress','+a')
		f.write(result[0]+"/run"+"\n")
		f.close()
		
	def callback_create(self, future):
		result = future.result()
		self.progress.append(result[0])
		self.update_progress(result[-1])
	def custom_error_callback(self,error):
   	 	print(f'Got an error: {error}')
   
	def callback_error(self,result):
		print('error', result)
  
	def do_job_async_unsrt_direct(self,data,generator,beta):
		for args in data:
			self.pool.apply_async(run_sampling_direct, 
				  args=(1,args,data[args],generator[args],beta,), 
				  callback=self.callback_direct,error_callback = self.callback_error)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
		return self.parallel_zeta_dict
	
	
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
		for args in locations:
			self.pool.apply_async(run_generate_dir, 
			     args=(args,len(locations)),callback=self.callback_run_)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
	
	def handler(self,error):
		print(f'Error: {error}', flush=True)
	
	def do_job_map_create_3(self,params):
		for param in params:
			self.pool.apply_async(run_map_3,args=(param,len(params)),callback=self.callback_run_,error_callback=self.custom_error_callback)
		
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
	
	def do_job_map_create_2(self,params):
		for param in params:
			self.pool.apply_async(run_map_2,args=(param,len(params)),callback=self.callback_run_,error_callback=self.custom_error_callback)
		#self.pool.map_async(run_map_2,params,error_callback=self.custom_error_callback)
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
		
	def do_job_async_(self,params):
		for param in params:
			self.pool.apply_async(run_executable_files_with_timeout, 
			     args=(param,len(params)),callback=self.callback_run,error_callback=self.custom_error_callback)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
		del self.pool, params
		del self.progress, self.parallized_zeta, self.parallel_zeta_dict, self.generator
	
	def do_job_create_opt(self,params):
		for param in params:
			self.pool.apply_async(run_map_opt, 
			     args=(param,len(params)),callback=self.callback_run_,error_callback=self.custom_error_callback)
		
		self.pool.close()
		self.pool.join()
		self.pool.terminate()	
	
	def do_job_map_create(self,params):
		for param in params:
			self.pool.apply_async(run_map, 
			     args=(param,len(params)),callback=self.callback_run_,error_callback=self.custom_error_callback)
		
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
	
	def do_job_executor(self,locations):
		with concurrent.futures.ProcessPoolExecutor() as executor:
			executor.map(run_generate_dir,locations)





class SM(object):
	def __init__(self,target_list,target_data,parameter_dict,designMatrix,ParameterDictionary ={}, flag="reactions", tag="Full"):
		"""
			
		target_list   = list of combustion_target_class object 
				for each targets
		target_data   = dictonary of the user provided information
				in target.opt file
		reaction_list = selected reactions for sensitivity analysis
	 
		"""
		
		self.target_list = target_list
		self.target_data = target_data
		self.parameter_list = []
		self.parameter_index = []
		self.parameter_dict = parameter_dict
		self.param_dict = ParameterDictionary
		self.flag = flag
		if self.flag == "reactions":
			for key in parameter_dict["reaction"]:
				self.parameter_index.append(key)
				self.parameter_list.append(parameter_dict["reaction"][key])
		else:	
			#print(parameter_dict)
			for key in parameter_dict:
				self.parameter_index.extend([key+"_cp",key+"_H",key+"_S"])
				self.parameter_list.append(parameter_dict[key])				
		self.unsrt = self.parameter_list
		self.design_matrix = designMatrix
		
		#print(len(self.design_matrix))
		self.case_dir = range(0,len(target_list))
		
		if tag == "Full":
			self.design_matrix_dict = {}
			for case in self.case_dir:
				self.design_matrix_dict[case] = designMatrix
		else:
			self.design_matrix_dict = designMatrix
		
		self.file_type = target_data["Inputs"]["fileType"]
		self.order = target_data["Stats"]["Order_of_PRS"]
		self.parallel_threads = target_data["Counts"]["parallel_threads"]
		self.responseOrder = target_data["Stats"]["Order_of_PRS"]
		self.mech_loc = target_data["Locations"]["mechanism"]
		if "perturbed_mech_file" in target_data["Locations"]:
			self.pert_mech_file = open(target_data["Locations"]["perturbed_mech_file"],"r").readlines()
		else:
			self.pert_mech_file = ""
		Prior_Mechanism = Parser(self.mech_loc).mech
		self.prior_mech = copy.deepcopy(Prior_Mechanism)
		#print(self.prior_mech)
		#self.copy_of_mech = copy.deepcopy(Prior_Mechanism)
		self.allowed_count = target_data["Counts"]["parallel_threads"]
		self.fuel = target_data["Inputs"]["fuel"]
		self.build_path = target_data["Bin"]["yaml_writer"]
	
	
	def getPerturbedMechLocation(self,yaml_list,location_mech,index_list):
		params = list(zip(yaml_list,location_mech,index_list))
		W = Worker(self.allowed_count,self.build_path)
		W.do_job_map_create_3(params)
		del W
	
	
	
	def getYAML_List(self,params,selection=[]):
		yaml_list = []
		#yaml_dict = {}
		sim_dict = []
		#print(selection)
		if len(selection) != 0:
			selection_params = selection
			for i in tqdm(range(len(params)),desc="Create Perturbed YAML files"):
				beta_ = params[i]
				select = selection_params[i]
				mani = manipulator(self.prior_mech,self.unsrt,beta_,selection = select,parameter_dict=self.param_dict,flag='thermo')
				yaml,sim = mani.doPerturbation()
				yaml_list.append(yaml)		
		
		else:
			for i in tqdm(range(len(params)),desc="Create Perturbed YAML files"):
				beta_ = params[i]
				mani = manipulator(self.prior_mech,self.parameter_dict,beta_,perturbation_type = "opt",parameter_dict=self.param_dict ,flag='thermo')
				#mani = manipulator(self.prior_mech,self.unsrt,beta_)
				yaml,sim = mani.doPerturbation()
				yaml_list.append(yaml)		
			
		return yaml_list
	
	def getDirectoryList(self,case,ind,yaml_loc_dict):
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
		#yaml_dict = {}
		instring_dict = {}
		s_convert_dict = {}
		s_run_dict = {}
		extract = {}	
		run_convert_dict = {}
		dir_list = []
		run_list = {}
		#sim_dict = {}
		#memo = {}
		#print(len(self.designMatrix))
		yaml_loc = yaml_loc_dict[case]
		#if self.prs_type == "Full":
		#	design_matrix = self.design_matrix[case]
		#else:
		design_matrix = self.design_matrix_dict[case]
		#print(len(design_matrix))	
		for i in tqdm(range(len(design_matrix)),desc="Zipping all files"):
			
			if os.path.isdir(os.getcwd()+"/"+str(i)) == True and os.path.isdir(os.getcwd()+"/"+str(i)+"/output") != True:
				shutil.rmtree(str(i))
				if self.pert_mech_file == "":
					#yaml_dict[str(i)],sim_dict[str(i)] = mani.doPerturbation()				
					instring_dict[str(i)],s_convert_dict[str(i)],s_run_dict[str(i)],extract[str(i)] = make_input_file.create_input_file(case,self.target_data,self.target_list[case],mech_file = yaml_loc[i]) #generate input file
				else:
					#yaml_dict[str(i)],sim_dict[str(i)] = "",""
					instring_dict[str(i)],s_convert_dict[str(i)],s_run_dict[str(i)],extract[str(i)] = make_input_file.create_input_file(case,self.target_data,self.target_list[case],mech_file = yaml_loc[i]) #generate input file
				dir_list.append(str(i))
				run_convert_dict[str(i)] = start+"/"+str(i)
				run_list[str(i)] = start+"/"+str(i)
				self.dir_list.append(start+"/"+str(i)+"/run")
			elif os.path.isdir(os.getcwd()+"/"+str(i)) != True:
				if self.pert_mech_file == "":
					#yaml_dict[str(i)],sim_dict[str(i)] = mani.doPerturbation()				
					instring_dict[str(i)],s_convert_dict[str(i)],s_run_dict[str(i)],extract[str(i)] = make_input_file.create_input_file(case,self.target_data,self.target_list[case],mech_file = yaml_loc[i]) #generate input file
				else:
					#yaml_dict[str(i)],sim_dict[str(i)] = "",""
					instring_dict[str(i)],s_convert_dict[str(i)],s_run_dict[str(i)],extract[str(i)] = make_input_file.create_input_file(case,self.target_data,self.target_list[case],mech_file = yaml_loc[i]) #generate input file
				dir_list.append(str(i))
				run_convert_dict[str(i)] = start+"/"+str(i)
				run_list[str(i)] = start+"/"+str(i)
				self.dir_list.append(start+"/"+str(i)+"/run")
				
			else:
				
				if self.pert_mech_file == "":
					#yaml_dict[str(i)],sim_dict[str(i)] = mani.doPerturbation()				
					instring_dict[str(i)],s_convert_dict[str(i)],s_run_dict[str(i)],extract[str(i)] = make_input_file.create_input_file(case,self.target_data,self.target_list[case],mech_file = yaml_loc[i]) #generate input file
				else:
					#yaml_dict[str(i)],sim_dict[str(i)] = "",""
					instring_dict[str(i)],s_convert_dict[str(i)],s_run_dict[str(i)],extract[str(i)] = make_input_file.create_input_file(case,self.target_data,self.target_list[case],mech_file = yaml_loc[i]) #generate input file
				run_convert_dict[str(i)] = start+"/"+str(i)
				run_list[str(i)] = start+"/"+str(i)
				self.dir_list.append(start+"/"+str(i)+"/run")
				continue	
		#print(yaml_dict)
		return instring_dict,s_run_dict,dir_list,run_convert_dict,run_list,extract
	def make_dir_in_parallel(self,yaml_loc):
		
		#Prior_Mechanism = Parser(self.mech_loc).mech
		#self.copy_of_mech = copy.deepcopy(Prior_Mechanism)
		
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
		sim_location = []	
		for case_index,case in enumerate(self.case_dir):
			if "case-"+str(case_index) in os.listdir():
				print("Case-{index} is generated".format(index = case_index))
			else:
				#Old version of the code
				#yaml_dict,instring_dict,s_run_dict,dir_list,run_convert_dict,run_list,extract_list,sim_dict = self.getDirectoryList(case,case_index)
				
				instring_dict,s_run_dict,dir_list,run_convert_dict,run_list,extract_list = self.getDirectoryList(case,case_index,yaml_loc)

				
				design_matrix = self.design_matrix_dict[case]
				###########################################
				##   Generated directories in parallel   ##
				##									   ##
				###########################################
				start_time = time.time()
				W = Worker(allowed_count,self.build_path)
				
				W.do_job_map(dir_list)
				print("\n\tDirectories for case - {} is generated\n".format(case))
				del W
				
				###########################################
				##   Generated required files			##
				##		in parallel					##
				###########################################
				
				#yaml_list = []
				#mech = []
				thermo = []
				trans = []
				instring = []
				run = []
				locations = []
				extract = []
				for i in tqdm(range(len(design_matrix))):
					instring.append(instring_dict[str(i)])
					#yaml_list.append(yaml_dict[str(i)])
					run.append(s_run_dict[str(i)])
					locations.append(run_list[str(i)])
					extract.append(extract_list[str(i)])

				##############################################################
				##   Generated required files	   				 ##
				##		in parallel	 (EXCEPT YAML files)			   ##
				##############################################################
				
				params = list(zip(instring,run,locations,extract))
				tic = time.time()
				
				V = Worker(allowed_count,self.build_path)
				V.do_job_create_opt(params)
				del V
				
				print("\n\t\tImportant files other than yaml mechanism file is created\n")
				
				###########################################
				##   Running the files				   ##
				##		in parallel					##
				###########################################
				
				X = Worker(allowed_count,self.build_path)
				file_n = []
				length = []
				for i in locations:
					file_n.append("run")
				args = list(zip(locations,file_n))
				X.do_job_async_(args)
				
				del X
				print("\tSimulations for case - {} is Done!!".format(case))
					
				dt = int(time.time() - start_time)
				hours = dt/3600
				minutes = (dt%3600)/60
				seconds = dt%60
				print("Performed {} Simulations....".format(len(locations)))
				print("Time for performing simulations : {h:.5f} hours,  {m:.3f} minutes, {s:.2f} seconds\n................................................ \n".format(h = hours, m = minutes, s =seconds))
				simulation_locations = open(optDir+"/locations",'+a')
				
				for loc in locations:
					simulation_locations.write(loc+"/run\n")
					sim_location.append(loc+"/run\n")
				simulation_locations.close()
			
				os.chdir('..')
		return sim_location

