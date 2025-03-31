#default python modules
from tqdm import tqdm
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
try:
    import ruamel_yaml as yaml
except ImportError:
    from ruamel import yaml
import concurrent.futures
import Uncertainty
### Program specific modules
import make_input_file
#import make_input_file2_0_A_factor 
from MechManipulator2_0 import Manipulator as manipulator
#from MechManipulator3_0_A_factor import Manipulator as manipulator
import data_management2_0_A_factor as data_management
from MechanismParser2_0_A_factor import Parser
from importlib import import_module

sys.path.append('/yamlwriter.so')  # Adjust this path to the correct build directory
#print(sys.path)
import yamlwriter

sys.path.append('/parallel_yaml_writer.so')
import parallel_yaml_writer
#print(dir(parallel_yaml_writer))
############################################
##  This module is to create directroies  ##
##  is systematic ways. Prallel computing ##
##  is done to make the process fast      ##
############################################

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
	yamlwriter.dump_to_yaml(location,f"mechanism.yaml",params[-2])
	#yaml_string = yaml.dump(params[-2],default_flow_style=False)
	#with open(location+"/mechanism.yaml","w") as yamlfile:
	#	yamlfile.write(yaml_string)

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
			self.pool.apply_async(run_executable_files_, 
			     args=(param,len(params)),callback=self.callback_run,error_callback=self.custom_error_callback)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
	
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
	def __init__(self,target_list,target_data,unsrt_data,design_matrix,tag="Full"):
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
		self.prs_type = target_data["Stats"]["PRS_type"]
		
		self.target_list = target_list
		self.target_data = target_data
		self.unsrt = unsrt_data
		#if self.prs_type == "Full":
		#	self.design_matrix = design_matrix
		#else:
		
		self.case_dir = range(0,len(target_list))
		
		if tag == "Full":
			self.design_matrix_dict = {}
			for case in self.case_dir:
				self.design_matrix_dict[case] = design_matrix
		else:
			self.design_matrix_dict = design_matrix
		self.file_type = target_data["Inputs"]["fileType"]
		self.order = target_data["Stats"]["Order_of_PRS"]
		self.parallel_threads = target_data["Counts"]["parallel_threads"]
		self.responseOrder = target_data["Stats"]["Order_of_PRS"]
		self.mech_loc = target_data["Locations"]["mechanism"]
		self.pre_file = target_data["Locations"]["Initial_pre_file"]
		if "perturbed_mech_file" in target_data["Locations"]:
			self.pert_mech_file = open(target_data["Locations"]["perturbed_mech_file"],"r").readlines()
		else:
			self.pert_mech_file = ""
		Prior_Mechanism = Parser(self.mech_loc).mech
		self.prior_mech = Prior_Mechanism
		self.copy_of_mech = copy.deepcopy(Prior_Mechanism)
		self.allowed_count = target_data["Counts"]["parallel_threads"]
		self.fuel = "MB-C5H10O2"#Only useful when using direct simulations!!
		self.build_path = target_data["Bin"]["yaml_writer"]
		#print(type(self.build_path))
		sys.path.append(str(self.build_path))
		#sys.path.append(build_path)  
		#self.yamlwriter = import_module('yamlwriter')
		#import yamlwriter
		#yamlwriter("./","tryal.yaml",Prior_Mechanism)
		#print(self.copy_of_mech["phases"][0]["species"])
		#raise AssertionError("Stop")
	
	
	def getNominalDirectoryList(self):		
		start = str(os.getcwd())
		yaml_dict = {}
		instring_dict = {}
		s_convert_dict = {}
		s_run_dict = {}
		extract = {}	
		run_convert_dict = {}
		run_list = {}
		sim_dict = {}
		case_dir = []
		pre_file = {}
		for case in tqdm(range(len(self.target_list)),desc="Zipping all files"):
			case_dir.append(os.getcwd()+"/case-"+str(case))
			if os.path.isdir(os.getcwd()+"/"+str(case)) == True and os.path.isdir(os.getcwd()+"/"+str(case)+"/output") != True:
				shutil.rmtree(str(case))
				yaml_dict[str(case)],sim_dict[str(case)],pre_file[str(case)] = self.copy_of_mech,"",self.pre_file
				instring_dict[str(case)],s_convert_dict[str(case)],s_run_dict[str(case)],extract[str(case)] = make_input_file.create_input_file(case,self.target_data,self.target_list[case]) #generate input file
				
				run_convert_dict[str(case)] = os.getcwd()+"/case-"+str(case)
				run_list[str(case)] = os.getcwd()+"/case-"+str(case)
				#self.dir_list.append(start+"/"+str(case)+"/run")
			
			elif os.path.isdir(os.getcwd()+"/"+str(case)) != True:
				yaml_dict[str(case)],sim_dict[str(case)],pre_file[str(case)] = self.copy_of_mech,"",self.pre_file
				instring_dict[str(case)],s_convert_dict[str(case)],s_run_dict[str(case)],extract[str(case)] = make_input_file.create_input_file(case,self.target_data,self.target_list[case]) #generate input file
				run_convert_dict[str(case)] = os.getcwd()+"/case-"+str(case)
				run_list[str(case)] = os.getcwd()+"/case-"+str(case)
				#self.dir_list.append(start+"/"+str(case)+"/run")
				
			else:
				
				yaml_dict[str(case)],sim_dict[str(case)],pre_file[str(case)] = self.copy_of_mech,"",self.pre_file
				instring_dict[str(case)],s_convert_dict[str(case)],s_run_dict[str(case)],extract[str(case)] = make_input_file.create_input_file(case,self.target_data,self.target_list[case]) #generate input file
				
				run_convert_dict[str(case)] = os.getcwd()+"/case-"+str(case)
				run_list[str(case)] = os.getcwd()+"/case-"+str(case)
				#self.dir_list.append(start+"/"+str(case)+"/run")
				continue	
		#print(yaml_dict)
		return  yaml_dict,instring_dict,s_run_dict,case_dir,run_convert_dict,run_list,extract,sim_dict,pre_file
	
	
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
				mani = manipulator(self.prior_mech,self.unsrt,beta_,selection = select)
				yaml,sim = mani.doPerturbation()
				yaml_list.append(yaml)		
		
		else:
			for i in tqdm(range(len(params)),desc="Create Perturbed YAML files"):
				beta_ = params[i]
				mani = manipulator(self.prior_mech,self.unsrt,beta_)
				yaml,sim = mani.doPerturbation()
				yaml_list.append(yaml)		
			
		return yaml_list
	
	def getPerturbedMechLocation(self,yaml_list,location_mech,index_list):
		params = list(zip(yaml_list,location_mech,index_list))
		W = Worker(self.allowed_count,self.build_path)
		W.do_job_map_create_3(params)
		del W

		
	
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
		#print(len(self.design_matrix))
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
		for ind,case in enumerate(range(cases)):
			file_dict = {}
			file_dict["beta"] = beta
			"""
			Get the zeta's using the beta
			"""
			count = 0
			generators = {}
			data = {}
			reactionList = [i for i in self.unsrt]
			beta_ = np.asarray(beta)
			
			file_dict["mechanism"],sim = manipulator(self.copy_of_mech,self.unsrt,beta_).doPerturbation()
			file_dict["simulationInputString"],file_dict["file_convertor_script"],file_dict["run_script"],file_dict["extractString"] = make_input_file2_0_A_factor.create_input_file(case,self.target_data,self.target_list[case])	
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

		#Prior_Mechanism = Parser(self.mech_loc).mech
		#self.copy_of_mech = copy.deepcopy(Prior_Mechanism)
		print(f"Starting direct simulations: Iteration number {iter_number}, total cases to simulate: {cases}")
		dictionary_list,mkdir_list,dir_run_list,run_convert_list,output_list = self.getSimulationFiles(case_index,cases,beta,iter_number)
		start_time = time.time()
			
		W = Worker(int(self.allowed_count),self.build_path)
		#W.do_job_executor(dir_list)
		W.do_job_map(mkdir_list)
		print("\tDirectories for direct simulations are generated\n")
		
		V = Worker(int(self.allowed_count),self.build_path)
		#print(params[0][6])
		#V.do_job_create_async(dir_list,instring_dict,mech_dict,thermo_dict,trans_dict,s_convert_dict,s_run_dict,run_list)
		V.do_job_direct_map(dictionary_list)
		print("\tRequired files for {} iteration is generated\n".format(iter_number))
		
		U = Worker(int(self.allowed_count),self.build_path)
		U.do_job_async_convertor(run_convert_list,"run_convertor")
		
		print("\tFiles converted to standard input files for iteration \n".format(iter_number))
		
		X = Worker(int(self.allowed_count),self.build_path)
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
		for i,index in enumerate(range(cases)):
			eta_list = data_management.extract_direct_simulation_values(index,output_list[i],self.target_list,self.fuel)
			directSimulation.extend(eta_list)
		os.chdir("..")
		#sample = open("samplefile.txt","+a")
		#sample.write(f"{iter_number},{objective}\n")
		return directSimulation
	
	def do_direct_sim(self,beta,case_index,cases,iter_number,objective):
		directSimulation = self.getSimulatedValues(beta,case_index,cases,iter_number,objective)
		return np.asarray(directSimulation)
		
	def make_nominal_dir_in_parallel(self):
		#Prior_Mechanism = Parser(self.mech_loc).mech
		#self.copy_of_mech = copy.deepcopy(Prior_Mechanism)
		
		start_time = time.time()
		print("====================================\nCreating Directories for simulations......\n\n\t This may take a while... Please be patient...\n\n ")
		allowed_count = int(self.parallel_threads)
		nominalDir = os.getcwd()
		yaml_dict,instring_dict,s_run_dict,case_dir,run_convert_dict,run_list,extract_list,sim_dict,pre_file = self.getNominalDirectoryList()
		W = Worker(allowed_count,self.build_path)
		###########################################
		##   Generated directories in parallel   ##
		##                                       ##
		###########################################
		W.do_job_map(case_dir)
		print("\tDirectories for all cases is generated\n")
		del W
		
		###########################################
		##   Generated required files            ##
		##        in parallel                    ##
		###########################################
		yaml_list = []
		mech = []
		thermo = []
		trans = []
		instring = []
		run = []
		locations = []
		extract = []
		pre_file_list = []
		for i in tqdm(range(len(self.case_dir)),desc="Generating required files"):
			instring.append(instring_dict[str(i)])
			yaml_list.append(yaml_dict[str(i)])
			run.append(s_run_dict[str(i)])
			locations.append(run_list[str(i)])
			extract.append(extract_list[str(i)])
			pre_file_list.append(pre_file[str(i)])
		params = list(zip(instring,run,locations,extract,yaml_list,pre_file_list))
		tic = time.time()
		
		V = Worker(allowed_count,self.build_path)
		V.do_job_map_create(params)
		del V	
		
		tok = time.time()
		print("\tRequired files for all cases is generated in {} hours, {} minutes, {} seconds time\n".format((int(tok-tic))/3600,((int(tok-tic))%3600)/60,(int(tok-tic))%60))
		
		
		###########################################
		##   Running the files                   ##
		##        in parallel                    ##
		###########################################
		
		X = Worker(allowed_count,self.build_path)
		file_n = []
		length = []
		for i in locations:
			file_n.append("run")
		args = list(zip(locations,file_n))
		X.do_job_async_(args)
		
		del X
		print("\tSimulations for all cases is Done!!")
			
		dt = int(time.time() - start_time)
		hours = dt/3600
		minutes = (dt%3600)/60
		seconds = dt%60
		print("Performed {} Simulations....".format(len(locations)))
		print("Time for performing simulations : {h:.4f} hours,  {m:.2f} minutes, {s:.2f} seconds\n................................................ \n".format(h = hours, m = minutes, s =seconds))
		simulation_locations = open(nominalDir+"/locations",'+a')
		sim_location = []
		for loc in locations:
			simulation_locations.write(loc+"/run\n")
			sim_location.append(loc+"/run\n")
		simulation_locations.close()
		
		#os.chdir('..')
		return sim_location	
		
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
				##                                       ##
				###########################################
				start_time = time.time()
				W = Worker(allowed_count,self.build_path)
				
				W.do_job_map(dir_list)
				print("\n\tDirectories for case - {} is generated\n".format(case))
				del W
				
				###########################################
				##   Generated required files            ##
				##        in parallel                    ##
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
				##   Generated required files       			     ##
				##        in parallel     (EXCEPT YAML files)               ##
				##############################################################
				
				params = list(zip(instring,run,locations,extract))
				tic = time.time()
				
				V = Worker(allowed_count,self.build_path)
				V.do_job_create_opt(params)
				del V
				
				print("\n\t\tImportant files other than yaml mechanism file is created\n")
				
				###########################################
				##   Running the files                   ##
				##        in parallel                    ##
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

