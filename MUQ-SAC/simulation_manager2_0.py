import os, shutil, re, math #default python modules
import make_input_file #program specific modules
import MechanismManipulator
import numpy as np
import pandas as pd
from pyDOE import *
home_dir = os.getcwd()
import multiprocessing
import subprocess
import time
import sys
import data_management
#function to create directories and the required input files in a systematic way directories are named with the reaction index
def run(location,n):
	return (locations,n)
	
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
		
		"""
		
		self.parameter_dict = data_management.extract_parameter_list(unsrt_data)
		self.target_list = target_list
		self.case_dir = range(0,len(target_list))
		self.file_type = target_data["Inputs"]["fileType"]
		self.order = order
		self.parallel_threads = target_data["Counts"]["parallel_threads"]
		self.responseOrder = target_data["Stats"]["Order_of_PRS"]
		

