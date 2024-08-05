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
import concurrent.futures
#try:
#    import ruamel_yaml as yaml
#except ImportError:
#    from ruamel import yaml
import yaml
############################################
##  This module is to create directroies  ##
##  is systematic ways. Prallel computing ##
##  is done to make the process fast      ##
############################################

def run_executable_files_(args,total):
	os.chdir(args[0])
	subprocess.call(["./"+args[1]])
	return (args[0],total)
def run_executable_files(location,file_name,n):
	os.chdir(location)
	subprocess.call(["./"+file_name])
	return (location,n)
def run_generate_dir(location,total):
	os.mkdir(location)
	os.mkdir(location+"/output")
	return (location,total)

def run_map_2(params,total):
	location = str(params[0])
	yaml_string = yaml.dump(params[1],default_flow_style=False)
	with open(location+"/mechanism.yaml","w") as yamlfile:
		yamlfile.write(yaml_string)
		yamlfile.close()
	return (params[0],total)


def run_map(params,total):
	location = str(params[2])
	
	sim1 = open(location+"/cantera_.py",'w').write(params[0])
	sim2= open(location+"/FlameMaster.input",'w').write(params[0])
	extract = open(location+"/extract.py",'w').write(params[3])
	perturb = open(location+"/perturb.txt",'w').write(params[4])
	#runConvertorScript = open(location+"/run_convertor",'w').write(params[2])
	runScript = open(location+"/run","w").write(params[1])
	#subprocess.call(["chmod","+x",location+"/run_convertor"])
	subprocess.call(["chmod","+x",location+"/run"])
	#yaml_string = yaml.dump(params[-1],default_flow_style=False)
	#with open(location+"/mechanism.yaml","w") as yamlfile:
	#	yamlfile.write(yaml_string)
	#	yamlfile.close()
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
	def __init__(self,target_list,target_data,reaction_list,designMatrix):
		"""
			
		target_list   = list of combustion_target_class object 
				for each targets
		target_data   = dictonary of the user provided information
				in target.opt file
		reaction_list = selected reactions for sensitivity analysis
	 
		"""
		
		self.target_list = target_list
		self.target_data = target_data
		self.reaction_list = []
		self.reaction_index = []
		for key in reaction_list["reaction"]:
			self.reaction_index.append(key)
			self.reaction_list.append(reaction_list["reaction"][key])
		
		self.unsrt = reaction_list
		self.design_matrix = designMatrix
		#print(len(self.design_matrix))
		self.case_dir = range(0,len(target_list))
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
		self.prior_mech = Prior_Mechanism
		#self.copy_of_mech = copy.deepcopy(Prior_Mechanism)
		self.allowed_count = target_data["Counts"]["parallel_threads"]
		self.fuel = target_data["Inputs"]["fuel"]
	
	
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
		
		for i in tqdm(range(len(self.design_matrix)),desc="Zipping all files"):
			dir_name = self.reaction_index[i]
			if self.pert_mech_file == "":
				self.beta_ = self.design_matrix[i]
				mani = manipulator(self.prior_mech,self.unsrt,self.beta_,perturbation_type = "SA")
			else:
				mech_file = self.pert_mech_file[i].strip("\n")
			#print(self.beta_)
			if os.path.isdir(os.getcwd()+"/"+str(dir_name)) == True and os.path.isdir(os.getcwd()+"/"+str(dir_name)+"/output") != True:
				shutil.rmtree(str(dir_name))
				if self.pert_mech_file == "":
					yaml_dict[str(dir_name)],sim_dict[str(dir_name)] = mani.doPerturbation()				
					instring_dict[str(dir_name)],s_convert_dict[str(dir_name)],s_run_dict[str(dir_name)],extract[str(dir_name)] = make_input_file.create_input_file(case,self.target_data,self.target_list[case]) #generate input file
				else:
					yaml_dict[str(dir_name)],sim_dict[str(dir_name)] = "",""
					instring_dict[str(dir_name)],s_convert_dict[str(dir_name)],s_run_dict[str(dir_name)],extract[str(dir_name)] = make_input_file.create_input_file(case,self.target_data,self.target_list[case],mech_file = mech_file) #generate input file
				dir_list.append(str(dir_name))
				run_convert_dict[str(dir_name)] = start+"/"+str(dir_name)
				run_list[str(dir_name)] = start+"/"+str(dir_name)
				self.dir_list.append(start+"/"+str(dir_name)+"/run")
			elif os.path.isdir(os.getcwd()+"/"+str(dir_name)) != True:
				if self.pert_mech_file == "":
					yaml_dict[str(dir_name)],sim_dict[str(dir_name)] = mani.doPerturbation()				
					instring_dict[str(dir_name)],s_convert_dict[str(dir_name)],s_run_dict[str(dir_name)],extract[str(dir_name)] = make_input_file.create_input_file(case,self.target_data,self.target_list[case]) #generate input file
				else:
					yaml_dict[str(dir_name)],sim_dict[str(dir_name)] = "",""
					instring_dict[str(dir_name)],s_convert_dict[str(dir_name)],s_run_dict[str(dir_name)],extract[str(dir_name)] = make_input_file.create_input_file(case,self.target_data,self.target_list[case],mech_file = mech_file) #generate input file
				dir_list.append(str(dir_name))
				run_convert_dict[str(dir_name)] = start+"/"+str(dir_name)
				run_list[str(dir_name)] = start+"/"+str(dir_name)
				self.dir_list.append(start+"/"+str(dir_name)+"/run")
				
			else:
				
				if self.pert_mech_file == "":
					yaml_dict[str(dir_name)],sim_dict[str(dir_name)] = mani.doPerturbation()				
					instring_dict[str(dir_name)],s_convert_dict[str(dir_name)],s_run_dict[str(dir_name)],extract[str(dir_name)] = make_input_file.create_input_file(case,self.target_data,self.target_list[case]) #generate input file
				else:
					yaml_dict[str(dir_name)],sim_dict[str(dir_name)] = "",""
					instring_dict[str(dir_name)],s_convert_dict[str(dir_name)],s_run_dict[str(dir_name)],extract[str(dir_name)] = make_input_file.create_input_file(case,self.target_data,self.target_list[case],mech_file = mech_file) #generate input file
				run_convert_dict[str(dir_name)] = start+"/"+str(dir_name)
				run_list[str(dir_name)] = start+"/"+str(dir_name)
				self.dir_list.append(start+"/"+str(dir_name)+"/run")
				continue	
		#print(yaml_dict)
		return  yaml_dict,instring_dict,s_run_dict,dir_list,run_convert_dict,run_list,extract,sim_dict
	
	def make_dir_in_parallel(self):
		start_time = time.time()
		print("\n\t\tCreating Directories for simulations......\n\n\n \t\tThis may take a while... Please be patient...\n\n ")
		allowed_count = int(self.parallel_threads)
		SADir = os.getcwd()
		self.dir_list = []			
		for case_index,case in enumerate(self.case_dir):
			if "case-"+str(case_index) in os.listdir():
				print("\n\t\tCase-{index} is generated".format(index = case_index))
			else:
				###############################################
				#####    Generating files for perturbation  ###
				#####	  of all reactions by x2            ###
				###############################################
				
				yaml_dict,instring_dict,s_run_dict,dir_list,run_convert_dict,run_list,extract_list,sim_dict = self.getDirectoryList(case,case_index)
						
				
				###########################################
				##   Generated directories in parallel   ##
				##                                       ##
				###########################################
				start_time = time.time()
				W = Worker(allowed_count)
				
				W.do_job_map(dir_list)
				print("\n\n\t\tDirectories for case - {} is generated\n".format(case))
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
				perturb = []
				for i in tqdm(yaml_dict):
					instring.append(instring_dict[str(i)])
					yaml_list.append(yaml_dict[str(i)])
					run.append(s_run_dict[str(i)])
					locations.append(run_list[str(i)])
					extract.append(extract_list[str(i)])
					perturb.append(sim_dict[str(i)])

				##############################################################
				##   Generated required files       			     ##
				##        in parallel     (EXCEPT YAML files)               ##
				##############################################################
				
				params = list(zip(instring,run,locations,extract,perturb))
				tic = time.time()
				
				V = Worker(allowed_count)
				V.do_job_map_create(params)
				del V
				print("\n\t\tImportant files other than yaml mechanism file is created\n")
				
				###############################################################
				##   Generated required files      (YAML files)    	      ##
				##       				                     ##
				##############################################################
				params_yaml = list(zip(locations,yaml_list))
				chunk_size = 100
				chunks = [params_yaml[i:i+chunk_size] for i in range(0, len(params_yaml), chunk_size)]
				
				tic = time.time()
				for args in chunks:
					W = Worker(allowed_count)
					W.do_job_map_create_2(args)
					del W
				
				tok = time.time()
				print("\n\t\tRequired files for case - {} is generated in {:.2f} hours, {:.2f} minutes, {:.2f} seconds time\n".format(case,(tok-tic)/3600,((tok-tic)%3600)/60,(tok-tic)%60))
				
				###########################################
				##   Running the files                   ##
				##        in parallel                    ##
				###########################################
				
				X = Worker(allowed_count)
				file_n = []
				length = []
				for i in locations:
					file_n.append("run")
				args = list(zip(locations,file_n))
				X.do_job_async_(args)
				
				del X
				print("\t\tSimulations for case - {} is Done!!".format(case))
					
				dt = int(time.time() - start_time)
				hours = dt/3600
				minutes = (dt%3600)/60
				seconds = dt%60
				print("\n\t\tPerformed {} Simulations....".format(len(locations)))
				print("\n\t\tTime for performing simulations : {h:.3f} hours,  {m:.2f} minutes, {s:.2f} seconds\n................................................ \n".format(h = hours, m = minutes, s =seconds))
				simulation_locations = open(SADir+"/locations",'+a')
				for loc in locations:
					simulation_locations.write(loc+"/run\n")
				simulation_locations.close()
			
				os.chdir('..')
		FlameMaster_Execution_location = []
		with open(SADir+"/locations") as infile:
			for line in infile:
				FlameMaster_Execution_location.append(line)
		return FlameMaster_Execution_location
