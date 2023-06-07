import os, shutil, re, math #default python modules
import make_input_file #program specific modules
import MechanismManipulator
from MechanismParser import Parser
import numpy as np
import pandas as pd
from pyDOE2 import *
home_dir = os.getcwd()
import multiprocessing
import subprocess
import time
import sys
import concurrent.futures
import asyncio
import data_management
import Uncertainty

#import ruamel.yaml as yaml
import yaml
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib import style
import matplotlib
matplotlib.use('Agg')
style.use("fivethirtyeight")
import random
from mpire import WorkerPool
from MechManipulator2_0 import Manipulator as manipulator
#function to create directories and the required input files in a systematic way directories are named with the reaction index
from copy import deepcopy


def run_sampling(sample,data,generator,length):
	A = Uncertainty.UncertaintyExtractor(data)
	#print(generator)
	a1 = generator[0]
	a2 = generator[1]
	
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
	#zeta = A.getLinearCombinationZeta()
	zeta = A.getB2Zeta(flag=True)
	#print(data)
	#zeta =UncertaintyExtractor(data).getExtremeCurves(sample)	
	del A
	#print(zeta)
	
	return (sample,generator,zeta,length)
	
def run_sampling_c(sample,data,generator,length):
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
	#zeta = A.getLinearCombinationZeta()
	zeta = A.getC2Zeta(flag=True)
	#print(data)
	#zeta =UncertaintyExtractor(data).getExtremeCurves(sample)	
	del A
	#print(zeta)
	
	return (sample,generator,zeta,length)

def run_sampling_b(sample,data,generator,length):
	A = Uncertainty.UncertaintyExtractor(data)
	#print(generator)
	a1 = generator[0]
	a2 = generator[1]
	#print(a1,a2)
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
	#zeta = A.getLinearCombinationZeta()
	zeta = A.getB2Zeta(flag=True)
	#print(data)
	#zeta =UncertaintyExtractor(data).getExtremeCurves(sample)	
	del A
	#print(zeta)
	
	return (sample,generator,zeta,length)

def mpire_run_sampling_b(sample,data,generator,length):
	#print(d)
	#sample = d["sample"]
	#data = d['data']
	#generator = d['generator']
	#length = d['length']
	A = Uncertainty.UncertaintyExtractor(data)
	#print(generator)
	a1 = generator[0]
	a2 = generator[1]
	#print(a1,a2)
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
	#zeta = A.getLinearCombinationZeta()
	zeta = A.getB2Zeta(flag=True)
	#print(data)
	#zeta =UncertaintyExtractor(data).getExtremeCurves(sample)	
	del A
	#print(zeta)
	
	return (sample,generator,zeta,length)

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

def run_sampling_kappa_direct(sample,rxn,data,generator,length):
	A = Uncertainty.UncertaintyExtractor(data)
	A.getCovariance(flag=False)
	A.getUnCorrelated(flag = False)
	zeta = A.getZeta_typeA(generator)
	del A	
	return (sample,rxn,generator,zeta,length)

def run_generate_dir(location):
	os.mkdir(location)
	os.mkdir(location+"/output")
	
def run_executable_files(args):
	os.chdir(args[0])
	#file_paste = open("run","w")
	#file_paste.write(string)
	#file_paste.close()
	#print(location[i])	
	#subprocess.check_call([str(location[i])+"/"+file_name])
	subprocess.call(["./"+args[1]])
	return (args[0],args[2])

def run_map(params):
	location = str(params[4])
	yaml_string = yaml.dump(params[0],default_flow_style=False)
	yamlfile = open(location+"/mechanism.yaml","w").write(yaml_string)
	sim = open(location+"/cantera_.py",'w').write(params[1])
	sim = open(location+"/FlameMaster.input",'w').write(params[1])
	extract = open(location+"/extract.py",'w').write(params[5])
	runConvertorScript = open(location+"/run_convertor",'w').write(params[2])
	runScript = open(location+"/run","w").write(params[3])
	subprocess.call(["chmod","+x",location+"/run_convertor"])
	subprocess.call(["chmod","+x",location+"/run"])
	del yaml_string
	del location
	#return (params[4])

def run_direct_map(file_dict):
	location = str(file_dict["mkdir"])
	#print(location)
	betaFile = open(location+"/beta.txt",'w').write(str(file_dict["beta"]))
	#mechFile = open(location+"/mechanism.mech",'w').write(file_dict["mechanism"])
	yaml_string = yaml.dump(file_dict["mechanism"],default_flow_style=False)
	yamlfile = open(location+"/mechanism.yaml","w").write(yaml_string)
	del yaml_string
	#thermoFile = open(location+"/thermo.therm",'w').write(file_dict["thermo"])
	#transportFile = open(location+"/transport.trans",'w').write(file_dict["transport"])
	sim = open(location+"/cantera_.py",'w').write(file_dict["simulationInputString"])
	sim = open(location+"/FlameMaster.input",'w').write(file_dict["simulationInputString"])
	extract = open(location+"/extract.py",'w').write(file_dict["extractString"])
	runConvertorScript = open(location+"/run_convertor",'w').write(file_dict["file_convertor_script"])
	runScript = open(location+"/run","w").write(file_dict["run_script"])
	#print(location+"/run_convertor")
	subprocess.call(["chmod","+x",location+"/run_convertor"])
	subprocess.call(["chmod","+x",location+"/run"])
	del location
	
def run_copy_files(i,sim_dict,mech_dict,thermo_dict,trans_dict,run_convert_dict,run_dict,locations,n):
	#os.mkdir(location)
	location = locations[i]
	#print(location)
	os.chdir(location)
	NewMechanismString = mech_dict[i]
	NewTransportString = trans_dict[i]
	NewThermoString = thermo_dict[i]
	NewRunConvertorScript = run_convert_dict[i]
	NewRunScript = run_dict[i]
	NewSimScript = sim_dict[i]
	sim = open("cantera_.py",'w')
	sim.write(NewSimScript)
	sim.close()
	mechFile = open("mechanism.mech",'w')
	mechFile.write(NewMechanismString)
	mechFile.close()
	thermoFile = open("thermo.therm",'w')
	thermoFile.write(NewThermoString)
	thermoFile.close()
	transportFile = open("transport.trans",'w')
	transportFile.write(NewTransportString)
	transportFile.close()
	runConvertorScript = open("run_convertor",'w')
	runConvertorScript.write(NewRunConvertorScript)
	runConvertorScript.close()
	runScript = open("run","w")
	runScript.write(NewRunScript)
	runScript.close()
	subprocess.call(["chmod","+x","run_convertor"])
	subprocess.call(["chmod","+x","run"])
	os.mkdir("output")
	return(i,n)
	
def create_files_for_directSimulation(file_dict,n):
	#os.mkdir(location)
	location = file_dict["mkdir"]
	#print(location)
	os.chdir(location)
	NewMechanismString = file_dict["mechanism"]
	NewTransportString = file_dict["transport"]
	NewThermoString = file_dict["thermo"]
	NewRunConvertorScript = file_dict["file_convertor_script"]
	NewRunScript = file_dict["run_script"]
	NewSimScript = file_dict["simulationInputString"]
	sim = open("cantera_.py",'w')
	sim.write(NewSimScript)
	sim.close()
	mechFile = open("mechanism.mech",'w')
	mechFile.write(NewMechanismString)
	mechFile.close()
	thermoFile = open("thermo.therm",'w')
	thermoFile.write(NewThermoString)
	thermoFile.close()
	transportFile = open("transport.trans",'w')
	transportFile.write(NewTransportString)
	transportFile.close()
	runConvertorScript = open("run_convertor",'w')
	runConvertorScript.write(NewRunConvertorScript)
	runConvertorScript.close()
	runScript = open("run","w")
	runScript.write(NewRunScript)
	runScript.close()
	subprocess.call(["chmod","+x","run_convertor"])
	subprocess.call(["chmod","+x","run"])
	os.mkdir("output")
	return(dictonary,n)
	
#class Ray():
#	def __init__(self,workers):
		

class Worker():
	def __init__(self, workers):
		self.pool_mpire = WorkerPool(n_jobs=workers)
		self.pool = multiprocessing.Pool(processes=workers)
		self.progress = []
		self.parallized_zeta = []
		self.generator = []
		self.parallel_zeta_dict = {}
		
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
		
	def callback_runConvertor(self, result):
		self.progress.append(result[0])
		sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(self.progress)/float(result[-1])*100))
		sys.stdout.flush()
		#self.pool.terminate()
	def callback_create(self, result):
		self.progress.append(result[0])
		sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(self.progress)/float(result[-1])*100))
		sys.stdout.flush()
		
	def callback_create_1(self, result):
		#string=""
		#for i in result:
		#	self.progress.append(i[0])
			#string+=f"{i[0]}\n"
			#total = i[1]
		#self.progress.append(result[0])
		sys.stdout.write("\t\t\r{} is complete".format(len(result)))
		sys.stdout.flush()
		#self.pool.terminate()

	def callback(self,result):
		#print("Entered callback\n")
		self.progress.append(result[0])
		self.generator.append(result[1])
		self.parallized_zeta.append(result[2])
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
	
	def callback_kappa_direct(self,result):
		#print("Entered callback\n")
		self.progress.append(result[0])
		self.generator.append(result[2])
		self.parallized_zeta.append(result[3])
		self.parallel_zeta_dict[result[1]] = result[3]
		sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(self.progress)/float(result[-1])*100))
		sys.stdout.flush()
		
	def callback_error(self,result):
    		print('error', result)
	def do_job_async_unsrt_direct(self,data,generator,beta):
		#print("Entered Async\n")
		#print("Entered async\n")
		for args in data:
			self.pool.apply_async(run_sampling_direct, 
				  args=(1,args,data[args],generator[args],beta,), 
				  callback=self.callback_direct)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
		return self.parallel_zeta_dict
	
	def do_job_async_unsrt_kappa_direct(self,data,generator,beta):
		#print("Entered Async\n")
		#print("Entered async\n")
		for args in data:
			self.pool.apply_async(run_sampling_kappa_direct, 
				  args=(1,args,data[args],generator[args],beta,), 
				  callback=self.callback_kappa_direct,error_callback=self.callback_error)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
		return self.parallel_zeta_dict
	
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
	
	def do_unsrt_c(self,data,sampling_points):
		#generators = data["generators"]
		for args in range(sampling_points):
			self.pool.apply_async(run_sampling_c, 
				  args=(1,data,data["generators_c"][args],sampling_points,), 
				  callback=self.callback,error_callback=self.callback_error)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
		return self.generator,self.parallized_zeta
	
	
	def MPIRE_do_unsrt_b(self,data,sampling_points):
		list_compre = [({"sample": 1,"data":data,"generator":data["generators_b"][args],"length":sampling_points}) for args in range(sampling_points)]
		for result in self.pool_mpire.imap(mpire_run_sampling_b,list_compre):
			#print(result)
			self.progress.append(result[0])
			self.generator.append(result[1])
			self.parallized_zeta.append(result[2])
			sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(self.progress)/float(result[-1])*100))
			sys.stdout.flush()
		self.pool_mpire.terminate()
		return self.generator,self.parallized_zeta
				
	
	def do_unsrt_b(self,data,sampling_points):
		for args in range(sampling_points):
			self.pool.apply_async(run_sampling_b, 
				  args=(1,data,data["generators_b"][args],sampling_points,), 
				  callback=self.callback)#,error_callback=self.callback_error)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
		return self.generator,self.parallized_zeta
	
	def do_job_async_convertor(self,args):
			#self.pool.apply_async(run_executable_files, 
			#	  args=(args,file_name,len(location)), 
			#	  callback=self.callback_runConvertor,error_callback=self.callback_error)
		self.pool.map_async(run_executable_files,args,callback=self.callback_run)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
	
	def do_job_async(self,args):

			#self.pool.apply_async(run_executable_files, 
			#	  args=(arg,file_name,len(location)), 
			#	  callback=self.callback_run,error_callback=self.callback_error)
		self.pool.map_async(run_executable_files,args,callback=self.callback_run,error_callback=self.callback_error)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
		
	def do_job_create_async(self,sim_list,sim_dict,mech_dict,thermo_dict,trans_dict,run_convert_dict,run_dict,locations):
		for args in sim_list:
			self.pool.apply_async(run_copy_files, 
				  args=(args,sim_dict,mech_dict,thermo_dict,trans_dict,run_convert_dict,run_dict,locations,len(sim_list)), 
				  callback=self.callback_create,error_callback=self.callback_error)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
	
	def do_job_direct_async(dictionary_list):
		for args in dictionary_list:
			self.pool.apply_async(create_files_for_directSimulation, 
				  args=(args,len(dictionary_list)), 
				  callback=self.callback_create,error_callback=self.callback_error)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
	
	def do_job_direct_map(self,params):
		self.pool.map_async(run_direct_map,params)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
		
	def do_job_map(self,locations):
		self.pool.map_async(run_generate_dir,locations)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
	
	def do_job_map_create(self,params):
		self.pool.map_async(run_map,params,callback=self.callback_create_1)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
	
	def do_job_executor(self,locations):
		with concurrent.futures.ProcessPoolExecutor() as executor:
			executor.map(run_generate_dir,locations)

class SM():
	def __init__(self,opt_dict,old_dict,cd, order,ap, ind, foc, m, th,tr, ml,ml_type,unsrt_data, rUsrt, fUsrt, mUsrt, thUsrt, trUsrt, tl, f, gr,thfl, trfl, spl, dis, sim, ac,selectedParams,manipulationDict,activeReactions, extra_arg=None):
		self.SAMAP_INDEX = 0
		self.samap_executable = opt_dict["Bin"]["samap_executable"]
		self.PRS_type = opt_dict["Stats"]["PRS_type"]
		self.selectedParams = selectedParams
		self.activeIndexDict = extra_arg
		self.activeReactionsDict = activeReactions
		self.opt_dict = opt_dict
		self.case_dir = cd
		self.order = order
		self.activeParameters = ap 
		self.ind = ind
		self.foc = foc 
		self.mBody = m
		self.thm = th
		self.trans = tr
		self.mech_loc = ml
		self.unsrt = unsrt_data
		self.rxnUnsert = rUsrt
		self.focUnsert = fUsrt
		self.tbdUnsert = mUsrt
		self.thermoUnsert = thUsrt
		self.transUnsert = trUsrt
		self.target_list = tl 
		self.fuel = f
		self.g_reaction = gr
		self.thermo_file_location = thfl
		self.trans_file_location = trfl 
		self.s_p_location = spl 
		self.design = dis
		self.simulation = sim
		self.allowed_count = ac
		self.progress = []
		self.fileType = ml_type
		self.old_dict = old_dict
		self.manipulatorDict = manipulationDict
		#print(ml_type)
		if ml_type.strip() == "chemkin":
			self.file_specific_command = "-f chemkin"
		else:
			self.file_specific_command = ""
		#print(self.file_specific_command)
	
	def make_directories_for_simulation(self,case):
		if os.path.isdir(os.getcwd()+"/case-"+str(case)) == True:
			os.chdir("case-"+str(case))
		else:
			os.mkdir("case-"+str(case))
			os.chdir("case-"+str(case))
		
		start = os.getcwd()
		mech_dict = {}
		thermo_dict = {}
		trans_dict = {}
		instring_dict = {}
		s_convert_dict = {}
		s_run_dict = {}
		dir_list = []
		run_convert_dict = {}
		run_list = {}
		extract = {}	
		
		for case in case_dir:
			if os.path.isdir(home_dir+"/case-"+str(case)) == True:
				continue
			os.mkdir("case-"+str(case))  	# make a dir for a target
			os.chdir("case-"+str(case))		#get into that dir
			os.mkdir('0')		# make folder for first order cases
			os.chdir('0')
			os.mkdir('0')		# make folder for nominal
			os.chdir('0')
			shutil.copyfile(mech_loc,'./mechanism.mech')
			make_input_file.create_input_file(target_list[case], fuel, g_reaction, bin_file_location)
			dir_list.append(os.path.abspath('run'))
			os.mkdir("output")
			os.chdir("..")
			for i in ind:
				os.mkdir(i)
				os.chdir(i)
				os.mkdir("plus")
				os.mkdir("minus")
				os.chdir("plus")
				mechanism = open("../../0/mechanism.mech",'r')
				create_mechanism_file(mechanism,i,unsert[i],"plus")
				mechanism.close()
				make_input_file.create_input_file(target_list[case], fuel, g_reaction, bin_file_location)
				dir_list.append(os.path.abspath('run'))
				os.mkdir("output")
				os.chdir('..')
				os.chdir("minus")
				mechanism = open("../../0/mechanism.mech",'r')
				create_mechanism_file(mechanism, i, unsert[i], "minus")
				mechanism.close()
				make_input_file.create_input_file(target_list[case], fuel, g_reaction, bin_file_location)
				dir_list.append(os.path.abspath('run'))
				os.mkdir("output")
				os.chdir('../..')
			os.chdir('..')
			for i_i, i in enumerate(ind[:-1]):
				os.mkdir(i)	# make folders to do second order perturbations
				os.chdir(i)	# get into folder for second order
				os.mkdir('0')	#make folder to store mech file of  first order
				shutil.copyfile('../0/'+i+"/plus/mechanism.mech", "./0/mechanism.mech")#/plus/
				for j_i,j in enumerate(ind):
					if j_i > i_i:
						os.mkdir(j)
						os.chdir(j)
						mechanism = open("../0/mechanism.mech",'r')
						create_mechanism_file(mechanism, j, unsert[j], "plus")
						mechanism.close()
						make_input_file.create_input_file(target_list[case], fuel, g_reaction, bin_file_location)
						dir_list.append(os.path.abspath('run'))
						os.mkdir("output")
						os.chdir("..")
					
				os.chdir("..")
			os.chdir("..")
		return dir_list
	
	
	def getZetaFromKappaCurve(self,beta):
		"""
		Get the zeta's using the beta
		"""
		count = 0
		generators = {}
		data = {}
		for i,rxn in enumerate(self.ind):
			data[rxn] = self.unsrt[rxn].data
			activeParams = self.unsrt[rxn].activeParameters
			gen = []
			for i in range(len(activeParams)):
				gen.append(beta[count])
				count+=1
			generators[rxn] = np.asarray(gen)	
		X = Worker(10)
		generator,zeta_dict = X.do_job_async_unsrt_kappa_direct(data,generators,len(beta))
		del X
		beta_ = []
		for i in self.ind:
			beta_.extend(list(zeta_dict[i]))
		return beta_
	
	
	def getZetaFromGenerators(self,beta):
		"""
		Get the zeta's using the beta
		"""
		count = 0
		generators = {}
		data = {}
		for i,rxn in enumerate(self.ind):
			data[rxn] = self.unsrt[rxn].data
			activeParams = self.unsrt[rxn].activeParameters
			gen = []
			for i in range(len(activeParams)):
				gen.append(beta[count])
				count+=1
			generators[rxn] = np.asarray(gen)	
		X = Worker(10)
		generator,zeta_dict = X.do_job_async_unsrt_direct(data,generators,len(beta))
		del X
		beta_ = []
		for i in self.ind:
			beta_.extend(list(zeta_dict[i]))
		return beta_


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
			if "A-facto" in self.design:
				file_dict = {}
				file_dict["beta"] = beta
				beta_ = []
				for i in beta:
					beta_.append(i)
			else:
				file_dict = {}
				file_dict["beta"] = beta
				"""
				Get the zeta's using the beta
				"""
			"""	
				count = 0
				generators = {}
				data = {}
				reactionList = self.ind
				for i,rxn in enumerate(self.ind):
					data[rxn] = self.unsrt[rxn].data
					activeParams = self.unsrt[rxn].activeParameters
					gen = []
					for i in range(len(activeParams)):
						gen.append(beta[count])
						count+=1
					generators[rxn] = np.asarray(gen)	
				X = Worker(10)
				zeta_dict = X.do_job_async_unsrt_direct(data,generators,len(beta))
				del X
				beta_ = []
				for i in self.ind:
					beta_.extend(list(zeta_dict[i]))
				beta_ = np.asarray(beta_)
			
			"""
			beta_= beta
			#print(beta_)
			file_dict["mechanism"],sim = manipulator(self.copy_of_mech,self.unsrt,beta_,self.ind,beta_transformed).doPerturbation()
			i = 100 #samples
			
			file_dict["simulationInputString"],file_dict["file_convertor_script"],file_dict["run_script"],file_dict["extractString"] = make_input_file.create_input_file(case,i,self.opt_dict,self.old_dict,self.target_list[case], self.fuel, self.g_reaction, self.thermo_file_location, self.trans_file_location,self.s_p_location,self.file_specific_command) #generate input file
			
			#yaml_dict[str(i)],sim_dict[str(i)] = manipulator(self.copy_of_mech,self.unsrt,self.beta_[i],self.ind,beta_transformed).doPerturbation()
			#file_dict["mechanism"],file_dict["thermo"],file_dict["transport"],sim = self.MechManipulator.GeneratePerturbedMechanism(self.target_list[case],self.beta_[i],beta_transformed,reactionList,self.simulation,"False",extra_arg = self.activeIndexDict)
			#file_dict["simulationInputString"],file_dict["file_convertor_script"],file_dict["run_script"],file_dict["extractString"] = make_input_file.create_input_file(case,iter_number,self.opt_dict,self.old_dict,cases[ind], self.fuel, self.g_reaction, self.thermo_file_location, self.trans_file_location,self.s_p_location,self.file_specific_command) #generate input file
		#dir_list.append(str(i))
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
		#print(os.getcwd())
		originalMech = Parser(self.mech_loc).mech
		self.copy_of_mech = deepcopy(originalMech)#.deepcopy()
		#self.MechManipulator = MechanismManipulator.MechanismManipulator(self.mech_loc,self.fileType,self.thermo_file_location,
		#								self.trans_file_location,self.ind,self.foc,self.mBody,
		#								self.thm,self.trans,self.unsrt, self.focUnsert, 
		#								self.tbdUnsert, self.thermoUnsert, self.transUnsert,
		#								self.selectedParams,self.activeIndexDict,design_type=self.design)
		print(f"Starting direct simulations: Iteration number {iter_number}, total cases to simulate: {len(case_index)}")
		dictionary_list,mkdir_list,dir_run_list,run_convert_list,output_list = self.getSimulationFiles(case_index,cases,beta,iter_number)
		start_time = time.time()
		
		#print(mkdir_list)
		#print(dictionary_list)
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
		file_n = []
		length = []
		for i in run_convert_list:
			file_n.append("run_convertor")
			length.append(len(run_convert_list))
		args_convert = list(zip(run_convert_list,file_n,length))
		U.do_job_async_convertor(args_convert)
		
		#print("\tFiles converted to standard input files for iteration \n".format(iter_number))
		
		file_n = []
		length = []
		for i in dir_run_list:
			file_n.append("run")
			length.append(len(dir_run_list))
		args_run = list(zip(dir_run_list,file_n,length))
		
		X = Worker(int(self.allowed_count))
		X.do_job_async(args_run)
		
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
		sample = open("samplefile.txt","+a")
		sample.write(f"{iter_number},{objective}\n")
		return directSimulation
		
	def getDirectoryList(self,case,ind):
		
		if os.path.isdir(os.getcwd()+"/case-"+str(case)) == True:
			os.chdir("case-"+str(case))
			#finished = True
			#if len(os.listdir()) == self.sim_:
			#	print("Simulation for case {} is done".format(case))
			#	finished = True
				#return mech_dict,thermo_dict,trans_dict,instring_dict,s_convert_dict,s_run_dict,dir_list,run_convert_dict,run_list,extract,sim_dict
		else:
			os.mkdir("case-"+str(case))
			os.chdir("case-"+str(case))
			#finished = False
		
		start = str(os.getcwd())
		yaml_dict = {}
		mech_dict = {}
		thermo_dict = {}
		trans_dict = {}
		instring_dict = {}
		s_convert_dict = {}
		s_convert_dict_2 = {}
		s_run_dict = {}
		dir_list = []
		run_convert_dict = {}
		run_list = {}
		extract = {}
		sim_dict = {}
		
		#print(self.activeIndexDict)
		#print(self.sim_)
		#finished = False
		#print(os.getcwd())
		for i in range(self.sim_):
			if self.PRS_type == "Partial" and self.simulation != "sa":
				#print(self.n)
				#self.beta_list.append(self.beta_[i]);
				xdata = []
				beta_transformed = (self.beta_[i])
				#print(beta_transformed)
				if self.activeIndexDict is not None:
					#print("Yes !! in mechmanipulator too")
					beta_transformed = []
					count = 0
					reactionList = []
					for z in self.activeReactionsDict[case]:
						reactionList.append(z)
					for index,ele in enumerate(self.activeIndexDict[ind]):
						if self.activeIndexDict[ind][ele] == 1:
							beta_transformed.append(1)
							xdata.append(self.beta_[i][count])
							#reactionList.append(ele)
							count+=1
						else:
							beta_transformed.append(0.0)
							count+=1
							#reactionList.append(ele)
					
					"""
					# Use only for debugging:
					print(reactionList)
					print("Active Reactions\n")
					print(f"{self.activeReactionsDict[case]}")
					print(f"case is {ind}")
					print("Active Index Dictionary\n")
					print(self.activeIndexDict[ind])
					print("\nBeta Transformed\n")
					print(beta_transformed)
					print("\nBeta \n")
					print(self.beta_[0])
					print("\nXdata\n")
					print(xdata)
					"""
				else:
					reactionList = []
					for z in self.ind:
						reactionList.append(z)
				
				
				
				#print(len(xdata))
				#print(len(self.beta_[i]))
				self.beta_list.append(xdata)
			else:
				#print(self.beta_[i])
				reactionList = []
				for z in self.ind:
					reactionList.append(z)
				self.beta_list.append(self.beta_[i]);
				#beta_transformed = self.manipulatorDict["selection"]
				beta_transformed = np.ones(len(self.beta_list[0]))
				
			#print(self.beta_list[0])
			#print(beta_transformed[0])
			#raise AssertionError("Stop!!")
			if os.path.isdir(os.getcwd()+"/"+str(i)) == True and os.path.isdir(os.getcwd()+"/"+str(i)+"/output") != True:
				shutil.rmtree(str(i))
				yaml_dict[str(i)],sim_dict[str(i)] = manipulator(self.copy_of_mech,self.unsrt,self.beta_[i],self.ind,beta_transformed).doPerturbation()
				#mech_dict[str(i)],thermo_dict[str(i)],trans_dict[str(i)],sim_dict[str(i)] = self.MechManipulator.GeneratePerturbedMechanism(self.target_list[case],self.beta_[i],beta_transformed,reactionList,self.simulation,"False",extra_arg = self.activeIndexDict)
				instring_dict[str(i)],s_convert_dict[str(i)],s_run_dict[str(i)],extract[str(i)] = make_input_file.create_input_file(case,i,self.opt_dict,self.old_dict,self.target_list[case], self.fuel, self.g_reaction, self.thermo_file_location, self.trans_file_location,self.s_p_location,self.file_specific_command) #generate input file
				dir_list.append(str(i))
				run_convert_dict[str(i)] = start+"/"+str(i)
				run_list[str(i)] = start+"/"+str(i)
				self.dir_list.append(start+"/"+str(i)+"/run")
			elif os.path.isdir(os.getcwd()+"/"+str(i)) != True:
				#mech_dict[str(i)],thermo_dict[str(i)],trans_dict[str(i)],sim_dict[str(i)] = self.MechManipulator.GeneratePerturbedMechanism(self.target_list[case],self.beta_[i],beta_transformed,reactionList,self.simulation,"False",extra_arg = self.activeIndexDict)
				yaml_dict[str(i)],sim_dict[str(i)] = manipulator(self.copy_of_mech,self.unsrt,self.beta_[i],self.ind,beta_transformed).doPerturbation()
				
				instring_dict[str(i)],s_convert_dict[str(i)],s_run_dict[str(i)],extract[str(i)] = make_input_file.create_input_file(case,i,self.opt_dict,self.old_dict,self.target_list[case], self.fuel, self.g_reaction, self.thermo_file_location, self.trans_file_location,self.s_p_location,self.file_specific_command) #generate input file
				dir_list.append(str(i))
				run_convert_dict[str(i)] = start+"/"+str(i)
				run_list[str(i)] = start+"/"+str(i)
				self.dir_list.append(start+"/"+str(i)+"/run")
			else:
				#mech_dict[str(i)],thermo_dict[str(i)],trans_dict[str(i)],sim_dict[str(i)] = self.MechManipulator.GeneratePerturbedMechanism(self.target_list[case],self.beta_[i],beta_transformed,reactionList,self.simulation,"False",extra_arg = self.activeIndexDict)
				yaml_dict[str(i)],sim_dict[str(i)] = manipulator(self.copy_of_mech,self.unsrt,self.beta_[i],self.ind,beta_transformed).doPerturbation()
				
				instring_dict[str(i)],s_convert_dict[str(i)],s_run_dict[str(i)],extract[str(i)] = make_input_file.create_input_file(case,i,self.opt_dict,self.old_dict,self.target_list[case], self.fuel, self.g_reaction, self.thermo_file_location, self.trans_file_location,self.s_p_location,self.file_specific_command) #generate input file
				#dir_list.append(str(i))
				run_convert_dict[str(i)] = start+"/"+str(i)
				run_list[str(i)] = start+"/"+str(i)
				self.dir_list.append(start+"/"+str(i)+"/run")
				continue	
		#return  mech_dict,thermo_dict,trans_dict,instring_dict,s_convert_dict,s_run_dict,dir_list,run_convert_dict,run_list,extract,sim_dict#,finished
		return  yaml_dict,instring_dict,s_convert_dict,s_run_dict,dir_list,run_convert_dict,run_list,extract,sim_dict#,finished
		
	def makeDir(self,case,total_targets):
		#print(self.simulation)
		self.beta_list = []	
		if self.simulation =='sa':
			self.sim_ = 7*self.n
		# taking 2 as perturbation factor
			self.beta_=[]
			c = []
			center_point = self.n
			for i in range(self.sim_):
				temp = 2*np.random.random(self.n)-1
				c.append(temp)
			beta = np.matrix(c)
			#beta = np.vstack(beta3)
			for i in beta:
				self.beta_.append(np.array(i)[0])		
	
	
		if self.simulation =='Opt':
			#Initializing the sampling size
			self.sim_ = len(2*np.random.random(2*int(self.n_))-1)
			#Initializing the design for the response surface
			if self.design == 'LHS-maximin':  #Latin hypercube sampling (Maximin)
				self.beta_ = 2*lhs(self.n,samples = self.sim_,criterion='centermaximin')-1	
			elif self.design == 'Monte-Carlo': #Monte-carlo sampling
				self.beta_ = []
				for i in range(self.sim_):
					temp = (2*np.random.random(self.n)-1)
					self.beta_.append(temp)
			elif self.design == 'Legendre':
				self.beta_ = np.matrix(pd.read_csv(home_dir+'/outfile_ortho.csv'))
				self.beta_ = self.beta_.tolist()
				#print(self.beta_[0])
				print("Beta list (Legendre) found !! Shape of beta list {}".format(np.shape(beta_)))
	
		if self.simulation == "Original":
			self.sim_ = 1
			self.beta_ = []
			for i in range(self.sim_):
				temp = 0*(2*np.random.random(self.n)-1)
				self.beta_.append(temp)
		
		if self.simulation == "Optimized":
			self.sim_ = 1
			self.beta_ = []
			for i in range(self.sim_):
				temp = 0*(2*np.random.random(self.n)-1)
				self.beta_.append(temp)
		
		self.sim_dict = {}
		if os.path.isdir(os.getcwd()+"/case-"+str(case)) == True:
			os.chdir("case-"+str(case))
		else:
			os.mkdir("case-"+str(case))
			os.chdir("case-"+str(case))  	# make a dir for a target
			#get into that dir
		#print("Required simulations = {}".format(self.sim_))
		for i in range(self.sim_):
			if os.path.isdir(os.getcwd()+"/"+str(i)) == True and os.path.isdir(os.getcwd()+"/"+str(i)+"/output") == True: 
				continue
			elif os.path.isdir(os.getcwd()+"/"+str(i)) == True and os.path.isdir(os.getcwd()+"/"+str(i)+"/output") != True:
			
			#	print("Output is not there, re creating the folder {}".format(i))
				os.chdir(str(i))
				self.beta_list.append(self.beta_[i]);
				#print(self.beta_[i])
				#beta_rxn,beta_foc,beta_mol,beta_ther,beta_tras=
				beta_dict = self.MechManipulator.GeneratePerturbedMechanism(self.target_list[case],self.beta_[i],beta_transformed,reactionList,self.simulation,"False",extra_arg = self.activeIndexDict)
				#print(beta_dict)
				#beta_r.append(beta_rxn)
				#beta_f.append(beta_foc)
				#beta_m.append(beta_mol)
				#beta_th.append(beta_ther)
				#beta_ts.append(beta_tras)
				self.sim_dict[str(i)] = beta_dict
				#print("starting to create the input files")
				make_input_file.create_input_file(case,i,self.opt_dict,self.old_dict,self.target_list[case], self.fuel, self.g_reaction, self.thermo_file_location, self.trans_file_location,self.s_p_location,self.file_specific_command) 
					#generate input file
				self.dir_list.append(os.path.abspath('run'))
				if "output" not in os.listdir(os.getcwd()):
					os.mkdir('output')
				os.chdir('..')
				continue
		
			elif os.path.isdir(os.getcwd()+"/"+str(i)) != True: 
				os.mkdir(str(i))
				os.chdir(str(i))
				self.beta_list.append(self.beta_[i]);
				#beta_rxn,beta_foc,beta_m,beta_th,beta_ts = 
				beta_dict = self.MechManipulator.GeneratePerturbedMechanism(self.target_list[case],self.beta_[i],beta_transformed,reactionList,self.simulation,"False",extra_arg = self.activeIndexDict)
				self.sim_dict[str(i)] = beta_dict
				make_input_file.create_input_file(case,i,self.opt_dict,self.old_dict,self.target_list[case], self.fuel, self.g_reaction, self.thermo_file_location, self.trans_file_location,self.s_p_location,self.file_specific_command) #generate input file
				self.dir_list.append(os.path.abspath('run').strip("\n"))
				if "output" not in os.listdir(os.getcwd()):
					os.mkdir('output')
				os.chdir('..')
				continue
			else:
				continue
		
		return (case,total_targets),self.sim_dict

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
	
	
	def getSamples(self):
		if self.simulation =='sa':

			if "Sampling_of_PRS" not in self.opt_dict["Stats"]:
				raise AssertionError("Error in opt_dict")
			self.sim_ = 10*self.n
			if self.design == "Monte-Carlo-All":
				mean = np.zeros(self.n)
				cov = 0.33*np.eye(self.n)
				self.beta_ = list(np.random.multivariate_normal(mean,cov,self.sim_))
				self.beta_.extend(list(np.eye(self.n)))
				self.beta_.extend(list(-1*np.eye(self.n)))
				self.beta_.extend([np.zeros(self.n)])
				self.sim_ = len(self.beta_)
			
			elif self.design == 'A-facto': 
				#print(self.n)
				self.beta_=[]
				task = "complete"
				if "Data" in os.listdir():
					os.chdir("Data/Simulations")
					if f"Beta_list_case-0.csv" in os.listdir():
						#print(great)
						#print(os.getcwd())
						
						A_file = open("Beta_list_case-0.csv","r").readlines() 
						self.sim_ = len(A_file)
						print(self.sim_)
						A = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in A_file])	
						A_gen_file = open("Generator_list_case-0.csv","r").readlines() 
						A_gen = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in A_gen_file])	
						
						
						X_file = open("NormalizedBeta_list_case-0.csv","r").readlines() 
						X_dash = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in A_file])
						y_dash = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in A_file])
						self.beta_ = np.asarray(A)
						#print(self.beta_)
						self.generators = np.asarray(A_gen)
						
						self.X = np.asarray(X_dash)
						self.y = np.asarray(y_dash)
						os.chdir("../..")
					else:
						os.chdir("../..")
						task = "incomplete"
				if task != "complete":
					self.beta_ = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
					self.beta_.extend(list(np.eye(self.n)))
					#self.beta_.extend(list(-1*np.eye(self.n)))
					#self.beta_.extend(list(bbdesign(self.n)[0:int(0.3*self.sim_)]))
					#self.beta_.extend(list(np.zeros(self.n)))
					#self.beta_.extend(list(np.zeros(self.n)))
					#self.beta_.extend(list(np.zeros(self.n)))
					self.sim_ = len(self.beta_)
					self.generators = np.asarray(self.beta_)
					self.X = np.asarray(self.beta_)
					self.y = np.asarray(self.beta_)
			
			elif self.design == "B1":
				self.beta_=[]
				temp_a = {}
				temp_b = {}
				temp_c = {}
				gen = {}
				gen_a = {}
				gen_b = {}
				gen_c = {}
				
				xdata = {}
				ydata = {}				
				self.generator_list = []
				"""
				Min- simulation = 300
				"""
				#if self.sim_ < 300:
				#	self.sim_ = 300
				#print(self.sim_)
				"""
				Dividing the generators into three bins
				------------------------------------
				n_a = 0.33*self.sim_
				n_b = 0.33*self.sim_
				n_c = 0.34*self.sim_
				
				"""
				n_a = int(0.33*self.sim_)
				n_b = int(0.33*self.sim_)
				n_c = self.sim_-n_a-n_b
				
				
				#n_b = int(0.5*self.sim_)
				#n_c = self.sim_-n_b
				
				
				#for i in np.eye(self.n):
				#	self.generator_list.append(i)
				#	self.generator_list.append(-1*i)
								
				"""
				Mixing three types of the extreme curves
				--------------------
				zeta-A: 1
				zeta-B: 1
				zeta-C: 1
			
				for zeta-A: one normal random parameter is required
					a matrix on n_a x n will be created
					
				for zeta-B and zeta-c: two random parameters is required
					we can arrange the two normal_random_parameter 				
					a matrix of n_b+n_c x 2*n will be created
					
				The number of reactions = n_R
				"""
				n_rxn = len(self.ind)
				self.generator_list_A = []
				self.generator_list_B = []
				self.generator_list_C = []
				
				"""
				Factorial Design
				----------------
				box-benkhen design
				bb = []
				for i in bbdesign(n_rxn)[0:int(0.3*self.sim_)]:
					bb.append(np.asarray(i))
				"""
				"""
						
						#For zeta-A
						
				bb_a = []
				for i in bbdesign(n_rxn):
					bb_a.append(np.asarray(i))
				self.generator_list_A.extend(bb_a)
						
						#For zeta-B
						
				bb_b = []
				for i in bbdesign(2*n_rxn):
					bb_b.append(np.asarray(i))							
				self.generator_list_B.extend(bb_b)
						
						#For zeta-C
						
				bb_c = []
				for i in bbdesign(2*n_rxn):
					bb_c.append(np.asarray(i))		
						
				self.generator_list_C.extend(bb_c)		
						
				"""
				"""
				Monte-Carlo
				
				"""
				"""		
						#For zeta-A
						
				
				self.generator_list_A = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				
				
						
						#For zeta-A
						
				self.generator_list_B = list(2* np.random.random_sample((self.sim_,self.n)) - 1)		
				
						
						#For zeta-A
						
				self.generator_list_C = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				"""
				
				"""
				Uniform
				
				"""
					
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				
				"""
				Latin-Hypercube
				"""
			
						
						#For zeta-A
						
				lhs_a = []
				for i in lhs(n_rxn, samples=n_a):
					lhs_a.append(np.asarray(i))	
				self.generator_list_A.extend(lhs_a)
				
						
						#For zeta-A
						
				lhs_b = []
				for i in lhs(2*n_rxn, samples=n_b):
					lhs_b.append(np.asarray(i))		
				self.generator_list_B.extend(lhs_b)
						
						#For zeta-A
						
				lhs_c = []
				for i in lhs(2*n_rxn, samples=n_c):
					lhs_c.append(np.asarray(i))
				self.generator_list_C.extend(lhs_c)
				
				"""
				Monte-Carlo-Mix
				"""
				
				
				
				
				
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)		
				#self.generator_list.extend(bb)
				#self.generator_list.extend(list(np.eye(self.n)))
				#self.generator_list.extend(list(-1*np.eye(self.n)))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#print(self.generator_list)
				#self.sim_ = len(self.generator_list)
				
				#trainers = np.asarray(trainers)
				#print(f"Generator_list {np.shape(self.generator_list)}\n")
				df = pd.DataFrame(np.asarray(self.generator_list))
				count_start_a = 0
				count_start_b = 0
				count_start_c = 0			
				for i in self.reactions_selected:
					
					len_of_Arrhenius_A = 1
					count_end_a = count_start_a+len_of_Arrhenius_A
					
					
					len_of_Arrhenius_B = 2
					count_end_b = count_start_b+len_of_Arrhenius_B
					
					
					len_of_Arrhenius_C = 2
					count_end_c = count_start_c+len_of_Arrhenius_C
					
					#print(count_start)
					#print(count_end)
					self.rxn_generators_a2 = []
					self.rxn_generators_b2 = []
					self.rxn_generators_c2 = []
					
					"""
					getting zetas from generators
					"""
					for j in range(len(self.generator_list_A)):
						self.rxn_generators_a2.append(self.generator_list_A[j][count_start_a:count_end_a])
					#print(self.rxn_generators_a2)
					for j in range(len(self.generator_list_B)):
						self.rxn_generators_b2.append(self.generator_list_B[j][count_start_b:count_end_b])
					
					for j in range(len(self.generator_list_C)):
						self.rxn_generators_c2.append(self.generator_list_C[j][count_start_c:count_end_c])
					
					self.sim_ = len(self.rxn_generators_a2)+len(self.rxn_generators_b2)+len(self.rxn_generators_c2)
					#print(self.rxn_generators)
					#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
										
					g_a = self.rxn_generators_a2
					g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					#g_b = self.rxn_generators[0:n_b]
					#g_c = self.rxn_generators[n_b:n_b+n_c]
					
					
					count_start_a = count_end_a
					count_start_b = count_end_b
					count_start_c = count_end_c
					
					"""
					Getting the uncertainty data
					"""
					#fig = plt.figure()
					data = self.unsrt[i].data
					T = np.linspace(300,2500,8)
					kappa_max = self.unsrt[i].getKappaMax(T)			
					kappa_min = self.unsrt[i].getKappaMin(T)
					Theta = np.array([T/T,np.log(T),-1/T])
					P = self.unsrt[i].getMean()	
					cov = self.unsrt[i].getCov()
					zeta_A = self.unsrt[i].zeta.x
					
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					"""
					Type-A zeta samples
					"""
					generators_a = g_a
					parallel_zetas_a = []
					
					time_start = time.time()
					
					
					for k in g_a:
						#Pint = P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten()
					#	plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
					#	plt.plot()
						#print(i)
						#print(type(i))
						#print(k)	
						#parallel_zetas_a.append(P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten())
						parallel_zetas_a.append(k[0]*zeta_A)
										
					#print(parallel_zetas_a)
					#plt.savefig(str(i)+"_a.pdf")
					#print(os.getcwd())
					#print("plot saved")
					data["generators_b"] = g_b
					callWorkForce = Worker(self.allowed_count)	
					generators_b,parallel_zetas_b = callWorkForce.do_unsrt_b(data,len(g_b))
					del callWorkForce
					
					data["generators_c"] = g_c
					callWorkForce = Worker(self.allowed_count)	
					generators_c,parallel_zetas_c = callWorkForce.do_unsrt_c(data,len(g_c))
					del callWorkForce
					
					#print(parallel_zetas_a)
					
					temp_a[i] = parallel_zetas_a
					gen_a[i] = generators_a
					temp_b[i] = parallel_zetas_b
					gen_b[i] = generators_b
					temp_c[i] = parallel_zetas_c
					gen_c[i] = generators_c
					
					time_end = time.time()
					print(f"\n\t\tTime taken = {time_end-time_start}\n")
					
					#print(parallel_zetas_b)
					Y_a = parallel_zetas_a
					Y_b = parallel_zetas_b
					Y_c = parallel_zetas_c
					T_ = np.linspace(300,3000,8)
					theta = np.array([T/T,np.log(T),-1/T])
					kappa_mean = theta.T.dot(P)
					Theta = np.array([T/T,np.log(T),-1/T])
					X = []
					Y = []
					
					
					for j in Y_a:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					
					for j in Y_b:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_b.pdf")
					
					
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					for j in Y_c:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_c.pdf")
					
					temp_xdata = {}
					for k in range(len(Y)):
						temp_xdata[str(Y[k])] = X[k]
						
					random.shuffle(Y)
					#print(Y_)
					X_ = []
					for k in range(len(Y)):
						X_.append(temp_xdata[str(Y[k])])
						
					#print(X_[9])
					xdata[i] = X_
					ydata[i] = Y
					#print(len(parallel_zetas))
					#print(len(self.rxn_generators))
				A = []
				A_gen = []
				X_dash = []
				y_dash = []
				
				for j in range(len(g_a)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_a:
						a_row.extend(list(temp_a[i][j]))
						a_row_gen.extend(list(gen_a[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				
				
				for j in range(len(g_b)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)]))
						y_temp.extend(list(ydata[i][j+len(g_a)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				for j in range(n_b):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				
				for j in range(len(g_c)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)+len(g_b)]))
						y_temp.extend(list(ydata[i][j+len(g_a)+len(g_b)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				for j in range(n_c):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+n_b]))
						y_temp.extend(list(ydata[i][j+n_b]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
					
				"""
				
				#print(A)
				#self.sim_ = len(g_a)+len(g_b)+len(g_c)
				#self.sim_ = n_b+n_c
				self.beta_ = np.asarray(A)
				#print(self.beta_)
				self.generators = np.asarray(A_gen)
				self.X = np.asarray(X_dash)
				self.y = np.asarray(y_dash)
			
			elif self.design == "A1+C1":
				
				self.beta_=[]
				temp_a = {}
				temp_b = {}
				temp_c = {}
				gen = {}
				gen_a = {}
				gen_b = {}
				gen_c = {}
				
				xdata = {}
				ydata = {}				
				self.generator_list = []
				"""
				Min- simulation = 300
				"""
				#if self.sim_ < 300:
				#	self.sim_ = 300
				#print(self.sim_)
				"""
				Dividing the generators into three bins
				------------------------------------
				n_a = 0.33*self.sim_
				n_b = 0.33*self.sim_
				n_c = 0.34*self.sim_
				
				"""
				n_a = int(0.1*self.sim_)
				n_b = int(0.45*self.sim_)
				n_c = self.sim_-n_b-n_a
				#print(n_a)
				#print(n_b)
				#print(n_c)
				
				#n_b = int(0.5*self.sim_)
				#n_c = self.sim_-n_b
				
				
				#for i in np.eye(self.n):
				#	self.generator_list.append(i)
				#	self.generator_list.append(-1*i)
								
				"""
				Mixing three types of the extreme curves
				--------------------
				zeta-A: 1
				zeta-B: 1
				zeta-C: 1
			
				for zeta-A: one normal random parameter is required
					a matrix on n_a x n will be created
					
				for zeta-B and zeta-c: two random parameters is required
					we can arrange the two normal_random_parameter 				
					a matrix of n_b+n_c x 2*n will be created
					
				The number of reactions = n_R
				"""
				n_rxn = len(self.ind)
				self.generator_list_A = []
				self.generator_list_B = []
				self.generator_list_C = []
				
				"""
				Factorial Design
				----------------
				box-benkhen design
				bb = []
				for i in bbdesign(n_rxn)[0:int(0.3*self.sim_)]:
					bb.append(np.asarray(i))
				"""
						#For zeta-A
						
				bb_a = []
				for i in bbdesign(n_rxn)[0:int(0.4*n_a)]:
					bb_a.append(np.asarray(i))
				self.generator_list_A.extend(bb_a)
						
						#For zeta-B
						
				bb_b = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_b)]:
					bb_b.append(np.asarray(i))							
				self.generator_list_B.extend(bb_b)
						
						#For zeta-C
						
				bb_c = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_c)]:
					bb_c.append(np.asarray(i))		
						
				self.generator_list_C.extend(bb_c)		
						
				
				"""
				Monte-Carlo
				
				"""
						
						#For zeta-A
				
				self.generator_list_A.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(0.3*n_a)])
				self.generator_list_A.extend(list(np.eye(n_rxn)))
				self.generator_list_A.extend(list(-1*np.eye(n_rxn)))
				self.generator_list_A.append(np.zeros(n_rxn))
				#self.generator_list_A.append(np.zeros(n_rxn))
				#self.generator_list_A.append(np.zeros(n_rxn))
				
						#For zeta-A
						
				self.generator_list_B.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_b)])	
				self.generator_list_B.extend(list(np.eye(2*n_rxn)))
				self.generator_list_B.extend(list(-1*np.eye(2*n_rxn)))
				#self.generator_list_B.append(np.zeros(2*n_rxn))
				#self.generator_list_B.append(np.zeros(2*n_rxn))
				#self.generator_list_B.append(np.zeros(2*n_rxn))
						
						#For zeta-A		
				self.generator_list_C.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_c)])
				self.generator_list_C.extend(list(np.eye(2*n_rxn)))
				self.generator_list_C.extend(list(-1*np.eye(2*n_rxn)))
				#self.generator_list_C.append(np.zeros(2*n_rxn))
				#self.generator_list_C.append(np.zeros(2*n_rxn))
				#self.generator_list_C.append(np.zeros(2*n_rxn))
				
				"""
				Uniform
				
				"""
					
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				
				"""
				Latin-Hypercube
				"""
			
						
						#For zeta-A
						
				lhs_a = []
				for i in lhs(n_rxn, samples=int(0.3*n_a)):
					lhs_a.append(np.asarray(i))	
				self.generator_list_A.extend(lhs_a)
				
						
						#For zeta-A
						
				lhs_b = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_b)):
					lhs_b.append(np.asarray(i))		
				self.generator_list_B.extend(lhs_b)
						
						#For zeta-A
						
				lhs_c = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_c)):
					lhs_c.append(np.asarray(i))
				self.generator_list_C.extend(lhs_c)
				
				"""
				Monte-Carlo-Mix
				"""
				
				
				
				
				
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)		
				#self.generator_list.extend(bb)
				#self.generator_list.extend(list(np.eye(self.n)))
				#self.generator_list.extend(list(-1*np.eye(self.n)))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#print(self.generator_list)
				#self.sim_ = len(self.generator_list)
				
				#trainers = np.asarray(trainers)
				#print(f"Generator_list {np.shape(self.generator_list)}\n")
				df = pd.DataFrame(np.asarray(self.generator_list))
				count_start_a = 0
				count_start_b = 0
				count_start_c = 0			
				for i in self.reactions_selected:
					
					len_of_Arrhenius_A = 1
					count_end_a = count_start_a+len_of_Arrhenius_A
					
					
					len_of_Arrhenius_B = 2
					count_end_b = count_start_b+len_of_Arrhenius_B
					
					
					len_of_Arrhenius_C = 2
					count_end_c = count_start_c+len_of_Arrhenius_C
					
					#print(count_start)
					#print(count_end)
					self.rxn_generators_a2 = []
					self.rxn_generators_b2 = []
					self.rxn_generators_c2 = []
					
					"""
					getting zetas from generators
					"""
					for j in range(len(self.generator_list_A)):
						self.rxn_generators_a2.append(self.generator_list_A[j][count_start_a:count_end_a])
					#print(self.rxn_generators_a2)
					for j in range(len(self.generator_list_B)):
						self.rxn_generators_b2.append(self.generator_list_B[j][count_start_b:count_end_b])
					#print(self.rxn_generators_b2)
					
					for j in range(len(self.generator_list_C)):
						self.rxn_generators_c2.append(self.generator_list_C[j][count_start_c:count_end_c])
					#print(self.rxn_generators_c2)
					self.sim_ = len(self.rxn_generators_a2)+len(self.rxn_generators_b2)+len(self.rxn_generators_c2)
					#print(self.rxn_generators)
					#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
									
					g_a = self.rxn_generators_a2
					g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					#g_b = self.rxn_generators[0:n_b]
					#g_c = self.rxn_generators[n_b:n_b+n_c]
					#print(len(g_b))
					#print(len(g_c))	
					
					count_start_a = count_end_a
					count_start_b = count_end_b
					count_start_c = count_end_c
					
					"""
					Getting the uncertainty data
					"""
					#fig = plt.figure()
					data = self.unsrt[i].data
					T = np.linspace(300,2500,8)
					kappa_max = self.unsrt[i].getKappaMax(T)			
					kappa_min = self.unsrt[i].getKappaMin(T)
					Theta = np.array([T/T,np.log(T),-1/T])
					P = self.unsrt[i].getMean()	
					cov = self.unsrt[i].getCov()
					zeta_A = self.unsrt[i].zeta.x
					
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					"""
					Type-A zeta samples
					"""
					generators_a = g_a
					parallel_zetas_a = []
					
					time_start = time.time()
					
					
					for k in g_a:
						#Pint = P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten()
					#	plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
					#	plt.plot()
						#print(i)
						#print(type(i))
						#print(k)	
						#parallel_zetas_a.append(P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten())
						parallel_zetas_a.append(k[0]*zeta_A)
										
					#print(parallel_zetas_a)
					#plt.savefig(str(i)+"_a.pdf")
					#print(os.getcwd())
					#print("plot saved")
					data["generators_b"] = g_b
					callWorkForce = Worker(self.allowed_count)	
					generators_b,parallel_zetas_b = callWorkForce.do_unsrt_b(data,len(g_b))
					del callWorkForce
					#print("zeta b finished")
					#print(len(parallel_zetas_b))
					#print(f"{parallel_zetas_b}")
					#raise AssertionError("Stop!")
					data["generators_c"] = g_c
					callWorkForce = Worker(self.allowed_count)	
					generators_c,parallel_zetas_c = callWorkForce.do_unsrt_c(data,len(g_c))
					del callWorkForce
					#raise AssertionError("Stop!")
					#print(parallel_zetas_a)
					
					temp_a[i] = parallel_zetas_a
					gen_a[i] = generators_a
					temp_b[i] = parallel_zetas_b
					gen_b[i] = generators_b
					temp_c[i] = parallel_zetas_c
					gen_c[i] = generators_c
					
					time_end = time.time()
					print(f"\n\t\tTime taken = {time_end-time_start}\n")
					
					#print(parallel_zetas_b)
					Y_a = parallel_zetas_a
					Y_b = parallel_zetas_b
					Y_c = parallel_zetas_c
					T_ = np.linspace(300,2500,8)
					theta = np.array([T/T,np.log(T),-1/T])
					kappa_mean = theta.T.dot(P)
					Theta = np.array([T/T,np.log(T),-1/T])
					X = []
					Y = []
					
					
					for j in Y_a:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					
					for j in Y_b:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_b.pdf")
					
					
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					for j in Y_c:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_c.pdf")
					
					temp_xdata = {}
					for k in range(len(Y)):
						temp_xdata[str(Y[k])] = X[k]
						
					random.shuffle(Y)
					#print(Y_)
					X_ = []
					for k in range(len(Y)):
						X_.append(temp_xdata[str(Y[k])])
						
					#print(X_[9])
					xdata[i] = X_
					ydata[i] = Y
					#print(len(parallel_zetas))
					#print(len(self.rxn_generators))
				A = []
				A_gen = []
				X_dash = []
				y_dash = []
				
				for j in range(len(g_a)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_a:
						a_row.extend(list(temp_a[i][j]))
						a_row_gen.extend(list(gen_a[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				
				
				for j in range(len(g_b)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)]))
						y_temp.extend(list(ydata[i][j+len(g_a)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				for j in range(n_b):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				
				for j in range(len(g_c)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)+len(g_b)]))
						y_temp.extend(list(ydata[i][j+len(g_a)+len(g_b)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				for j in range(n_c):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+n_b]))
						y_temp.extend(list(ydata[i][j+n_b]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
					
				"""
				
				#print(A)
				#self.sim_ = len(g_a)+len(g_b)+len(g_c)
				#self.sim_ = n_b+n_c
				self.beta_ = np.asarray(A)
				#print(len(self.beta_))
				#print(self.beta_)
				self.generators = np.asarray(A_gen)
				self.X = np.asarray(X_dash)
				self.y = np.asarray(y_dash)
			
			elif self.design == "SAMAP":
				"""
				Getting the SAMAP amples
				"""
				os.mkdir("SAMAP_"+str(self.SAMAP_INDEX))
				os.chdir("SAMAP_"+str(self.SAMAP_INDEX))
				input_rxn_dict = {}
				samples_dict = {}
				for i in self.reactions_selected:
					os.mkdir(f"{self.unsrt[i].rIndex}")
					os.chdir(f"{self.unsrt[i].rIndex}")
					input_dict = self.unsrt[i].input_dict
					input_dict["samples"] = self.sim_
					input_dict["samples_skipped"] = int(0.1*self.sim_)
					input_dict["Random_seed"] = 1
					input_dict["sampling_method"] = "SOBOL"
					input_dict["sampling_distribution"] = "NORMAL"
					input_dict["equidistant_T"] = 100
					input_dict["T_begin"] = data = self.unsrt[i].temperatures[0]
					input_dict["T_end"] = self.unsrt[i].temperatures[-1]
					
					L = self.unsrt[i].cholskyDeCorrelateMat
					#input_dict["covariance_matrix"] = str(L.dot(L.T)).strip("[]").replace("[","").replace("]","")
					input_rxn_dict[i] = input_dict
					"""
					Run: SAMAP code
					"""
					instring = make_input_file.create_SAMAP_input(input_dict)
					#string_dict[i] = instring
					file_print = open("samap_data_"+str(self.unsrt[i].rIndex)+".txt","w").write(instring)
					run_string = f"""#!/bin/bash
		{self.samap_executable} samap_data_{self.unsrt[i].rIndex}.txt &> out"""
					file_print_run = open("run","w").write(run_string)						
					subprocess.call(["chmod","+x",'run'])
					subprocess.call(["./run"])
					file_name = "samap_data_"+str(self.unsrt[i].rIndex)+".txt_Arrhpar.txt"
					sampling_file = open(file_name,"r").readlines()
					self.Nagy_arrhenius_samples = [] 
					for sample in sampling_file[1:]:
						self.Nagy_arrhenius_samples.append(np.asfarray(sample.strip("''").strip("\n").split()[:3],float))
					
					#print(self.Nagy_arrhenius_samples)
					Y = self.Nagy_arrhenius_samples
					
					
					data = self.unsrt[i].data
					
					P = self.unsrt[i].getMean()	
					#print(P)
					#cov = self.unsrt[i].getCov()
					cov = input_dict["cov_float"]
					L_samap = np.linalg.cholesky(cov)
					T_ = np.linspace(300,2500,100)
					Theta = np.array([T_/T_,np.log(T_),-1/T_])
					
					fig = plt.figure()
					plt.plot(1/T_,self.unsrt[i].getKappaMax(T_),'k-',linewidth=0.75)
					plt.plot(1/T_,self.unsrt[i].getKappaMin(T_),'k-',linewidth=0.75)
					plt.plot(1/T_,self.unsrt[i].getNominal(T_),"b-",linewidth=0.75)
					
					temp = []
					for j in Y:
						#print(j)
						#print(kappa_max,np.asarray(theta.T.dot(j)).flatten())
						#print(np.asarray(theta.T.dot(j)).flatten()/kappa_max)
						#Pint = P + np.asarray(np.dot(cov,np.asarray(j))).flatten();
						Pint = np.asarray(j)
						ZETA = np.asarray(np.linalg.inv(L_samap).dot(Pint-P)).flatten()
						temp.append(ZETA)
							#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
						#X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						plt.plot(1/T_,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
					samples_dict[i] = temp
					plt.savefig(str(i)+"Samap.pdf",bbox_inches="tight")
				
					os.chdir("..")
					
				A = []
				A_gen = []
				X = []
				y = []
				for j in range(self.sim_):
					a_row = []
					for i in samples_dict:
						a_row.extend(list(samples_dict[i][j]))
					
					A.append(np.asarray(a_row))
					
				self.beta_ = np.asarray(A)
				#print(self.beta_)
				self.generators = []
				self.X = []
				self.y = []
				
				self.SAMAP_INDEX+=1
				os.chdir("..")
				
			elif self.design == "B1+C1":
				self.beta_=[]
				temp_a = {}
				temp_b = {}
				temp_c = {}
				gen = {}
				gen_a = {}
				gen_b = {}
				gen_c = {}
				
				xdata = {}
				ydata = {}				
				self.generator_list = []
				"""
				Min- simulation = 300
				"""
				#if self.sim_ < 300:
				#	self.sim_ = 300
				#print(self.sim_)
				"""
				Dividing the generators into three bins
				------------------------------------
				n_a = 0.33*self.sim_
				n_b = 0.33*self.sim_
				n_c = 0.34*self.sim_
				
				"""
				n_a = int(0.1*self.sim_)
				#n_b = int(0.5*self.sim_)
				#n_c = self.sim_-n_a-n_b
				
				
				n_b = int(0.4*self.sim_)
				n_c = int(0.4*self.sim_)
				
				
				#for i in np.eye(self.n):
				#	self.generator_list.append(i)
				#	self.generator_list.append(-1*i)
								
				"""
				Mixing three types of the extreme curves
				--------------------
				zeta-A: 1
				zeta-B: 1
				zeta-C: 1
			
				for zeta-A: one normal random parameter is required
					a matrix on n_a x n will be created
					
				for zeta-B and zeta-c: two random parameters is required
					we can arrange the two normal_random_parameter 				
					a matrix of n_b+n_c x 2*n will be created
					
				The number of reactions = n_R
				"""
				n_rxn = len(self.ind)
				self.generator_list_A = []
				self.generator_list_B = []
				self.generator_list_C = []
				
				"""
				Factorial Design
				----------------
				box-benkhen design
				bb = []
				for i in bbdesign(n_rxn)[0:int(0.3*self.sim_)]:
					bb.append(np.asarray(i))
				"""
						#For zeta-A
						
				bb_a = []
				for i in bbdesign(n_rxn)[0:int(0.4*n_a)]:
					bb_a.append(np.asarray(i))
				self.generator_list_A.extend(bb_a)
						
						#For zeta-B
						
				bb_b = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_b)]:
					bb_b.append(np.asarray(i))							
				self.generator_list_B.extend(bb_b)
						
						#For zeta-C
						
				bb_c = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_c)]:
					bb_c.append(np.asarray(i))		
						
				self.generator_list_C.extend(bb_c)		
						
				
				"""
				Monte-Carlo
				
				"""
						
						#For zeta-A
				
				self.generator_list_A.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(0.3*n_a)])
				self.generator_list_A.extend(list(np.eye(n_rxn)))
				self.generator_list_A.extend(list(-1*np.eye(n_rxn)))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				
						#For zeta-A
						
				self.generator_list_B.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_b)])	
				self.generator_list_B.extend(list(np.eye(2*n_rxn)))
				self.generator_list_B.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
						
						#For zeta-A		
				self.generator_list_C.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_c)])
				self.generator_list_C.extend(list(np.eye(2*n_rxn)))
				self.generator_list_C.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				
				"""
				Uniform
				
				"""
					
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				
				"""
				Latin-Hypercube
				"""
				
						#For zeta-A
						
				lhs_a = []
				for i in lhs(n_rxn, samples=int(0.3*n_a)):
					lhs_a.append(np.asarray(i))	
				self.generator_list_A.extend(lhs_a)
				
						
						#For zeta-A
						
				lhs_b = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_b)):
					lhs_b.append(np.asarray(i))		
				self.generator_list_B.extend(lhs_b)
						
						#For zeta-A
						
				lhs_c = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_c)):
					lhs_c.append(np.asarray(i))
				self.generator_list_C.extend(lhs_c)
				
				"""
				Monte-Carlo-Mix
				"""
				
				
				
				
				
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)		
				#self.generator_list.extend(bb)
				#self.generator_list.extend(list(np.eye(self.n)))
				#self.generator_list.extend(list(-1*np.eye(self.n)))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#print(self.generator_list)
				#self.sim_ = len(self.generator_list)
				
				#trainers = np.asarray(trainers)
				#print(f"Generator_list {np.shape(self.generator_list)}\n")
				df = pd.DataFrame(np.asarray(self.generator_list))
				count_start_a = 0
				count_start_b = 0
				count_start_c = 0			
				for i in self.reactions_selected:
					
					len_of_Arrhenius_A = 1
					count_end_a = count_start_a+len_of_Arrhenius_A
					
					
					len_of_Arrhenius_B = 2
					count_end_b = count_start_b+len_of_Arrhenius_B
					
					
					len_of_Arrhenius_C = 2
					count_end_c = count_start_c+len_of_Arrhenius_C
					
					#print(count_start)
					#print(count_end)
					self.rxn_generators_a2 = []
					self.rxn_generators_b2 = []
					self.rxn_generators_c2 = []
					
					"""
					getting zetas from generators
					"""
					for j in range(len(self.generator_list_A)):
						self.rxn_generators_a2.append(self.generator_list_A[j][count_start_a:count_end_a])
					#print(self.rxn_generators_a2)
					for j in range(len(self.generator_list_B)):
						self.rxn_generators_b2.append(self.generator_list_B[j][count_start_b:count_end_b])
					
					for j in range(len(self.generator_list_C)):
						self.rxn_generators_c2.append(self.generator_list_C[j][count_start_c:count_end_c])
					
					#self.sim_ = len(self.rxn_generators_a2)+len(self.rxn_generators_b2)+len(self.rxn_generators_c2)
					self.sim_ = len(self.rxn_generators_b2)+len(self.rxn_generators_c2)
					#print(self.rxn_generators)
					#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
										
					g_a = self.rxn_generators_a2
					g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					#g_b = self.rxn_generators[0:n_b]
					#g_c = self.rxn_generators[n_b:n_b+n_c]
					
					
					count_start_a = count_end_a
					count_start_b = count_end_b
					count_start_c = count_end_c
					
					"""
					Getting the uncertainty data
					"""
					#fig = plt.figure()
					data = self.unsrt[i].data
					T = np.linspace(300,2500,8)
					kappa_max = self.unsrt[i].getKappaMax(T)			
					kappa_min = self.unsrt[i].getKappaMin(T)
					Theta = np.array([T/T,np.log(T),-1/T])
					P = self.unsrt[i].getMean()	
					cov = self.unsrt[i].getCov()
					zeta_A = self.unsrt[i].zeta.x
					
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					"""
					Type-A zeta samples
					"""
					generators_a = g_a
					parallel_zetas_a = []
					
					time_start = time.time()
					
					
					for k in g_a:
						#Pint = P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten()
					#	plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
					#	plt.plot()
						#print(i)
						#print(type(i))
						#print(k)	
						#parallel_zetas_a.append(P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten())
						parallel_zetas_a.append(k[0]*zeta_A)
										
					#print(parallel_zetas_a)
					#plt.savefig(str(i)+"_a.pdf")
					#print(os.getcwd())
					#print("plot saved")
					g_a = self.rxn_generators_a2
					g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					
					generators_b = []
					parallel_zetas_b = []
					generators_c = []
					parallel_zetas_c = []
					start = 0
					for w in range(10):
						stop = int(((w+1)/10)*self.sim_)
						g_b = self.rxn_generators_b2[start:stop]
						g_c = self.rxn_generators_c2[start:stop]
						data["generators_b"] = g_b
						callWorkForce = Worker(self.allowed_count)	
						generators_b_,parallel_zetas_b_ = callWorkForce.do_unsrt_b(data,len(g_b))
						del callWorkForce
						
						
						removed = data.pop("generators_b")
						data["generators_c"] = g_c
						callWorkForce = Worker(self.allowed_count)	
						generators_c_,parallel_zetas_c_ = callWorkForce.do_unsrt_c(data,len(g_c))
						del callWorkForce
						
						generators_b.extend(generators_b_)
						parallel_zetas_b.extend(parallel_zetas_b_)
						generators_c.extend(generators_c_)
						parallel_zetas_c.extend(parallel_zetas_c_)
						removed = data.pop("generators_c")
						start = stop
					
					#print(self.sim_)
					#print(len(parallel_zetas_b)+len(parallel_zetas_c))
					g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					#print(parallel_zetas_a)
					#print(parallel_zetas_b)
					temp_a[i] = parallel_zetas_a
					gen_a[i] = generators_a
					temp_b[i] = parallel_zetas_b
					gen_b[i] = generators_b
					temp_c[i] = parallel_zetas_c
					gen_c[i] = generators_c
					
					time_end = time.time()
					print(f"\n\t\tTime taken = {time_end-time_start}\n")
					
					#print(parallel_zetas_b)
					Y_a = parallel_zetas_a
					Y_b = parallel_zetas_b
					Y_c = parallel_zetas_c
					T_ = np.linspace(300,2500,8)
					theta = np.array([T/T,np.log(T),-1/T])
					kappa_mean = theta.T.dot(P)
					Theta = np.array([T/T,np.log(T),-1/T])
					X = []
					Y = []
					
					"""
					for j in Y_a:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					"""
					for j in Y_b:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_b.pdf")
					
					
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					for j in Y_c:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_c.pdf")
					
					temp_xdata = {}
					for k in range(len(Y)):
						temp_xdata[str(Y[k])] = X[k]
						
					random.shuffle(Y)
					#print(Y_)
					X_ = []
					for k in range(len(Y)):
						X_.append(temp_xdata[str(Y[k])])
						
					#print(X_[9])
					xdata[i] = X_
					ydata[i] = Y
					#print(len(parallel_zetas))
					#print(len(self.rxn_generators))
				A = []
				A_gen = []
				X_dash = []
				y_dash = []
				"""
				for j in range(len(g_a)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_a:
						a_row.extend(list(temp_a[i][j]))
						a_row_gen.extend(list(gen_a[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				"""
				for j in range(len(g_b)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)]))
						y_temp.extend(list(ydata[i][j+len(g_a)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				print(temp_b)
				
				for j in range(len(g_b)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
					#print(a_row)
				"""
				
				for j in range(len(g_c)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)+len(g_b)]))
						y_temp.extend(list(ydata[i][j+len(g_a)+len(g_b)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				for j in range(len(g_c)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_b)]))
						y_temp.extend(list(ydata[i][j+len(g_b)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
					#print(a_row)
				
				
				#print(A)
				#self.sim_ = len(g_a)+len(g_b)+len(g_c)
				#self.sim_ = n_b+n_c
				self.beta_ = np.asarray(A)
				#print(self.beta_)
				self.generators = np.asarray(A_gen)
				self.X = np.asarray(X_dash)
				self.y = np.asarray(y_dash)
			
			
			elif self.design == "A1+B1":
				self.beta_=[]
				temp_a = {}
				temp_b = {}
				temp_c = {}
				gen = {}
				gen_a = {}
				gen_b = {}
				gen_c = {}
				
				xdata = {}
				ydata = {}				
				self.generator_list = []
				"""
				Min- simulation = 300
				"""
				#if self.sim_ < 300:
				#	self.sim_ = 300
				#print(self.sim_)
				"""
				Dividing the generators into three bins
				------------------------------------
				n_a = 0.33*self.sim_
				n_b = 0.33*self.sim_
				n_c = 0.34*self.sim_
				
				"""
				n_a = int(0.1*self.sim_)
				#n_b = int(0.5*self.sim_)
				#n_c = self.sim_-n_a-n_b
				
				
				n_b = int(0.9*self.sim_)
				n_c = self.sim_-n_b
				
				
				#for i in np.eye(self.n):
				#	self.generator_list.append(i)
				#	self.generator_list.append(-1*i)
								
				"""
				Mixing three types of the extreme curves
				--------------------
				zeta-A: 1
				zeta-B: 1
				zeta-C: 1
			
				for zeta-A: one normal random parameter is required
					a matrix on n_a x n will be created
					
				for zeta-B and zeta-c: two random parameters is required
					we can arrange the two normal_random_parameter 				
					a matrix of n_b+n_c x 2*n will be created
					
				The number of reactions = n_R
				"""
				n_rxn = len(self.ind)
				self.generator_list_A = []
				self.generator_list_B = []
				self.generator_list_C = []
				
				"""
				Factorial Design
				----------------
				box-benkhen design
				bb = []
				for i in bbdesign(n_rxn)[0:int(0.3*self.sim_)]:
					bb.append(np.asarray(i))
				"""
						#For zeta-A
						
				bb_a = []
				for i in bbdesign(n_rxn)[0:int(0.4*n_a)]:
					bb_a.append(np.asarray(i))
				self.generator_list_A.extend(bb_a)
						
						#For zeta-B
						
				bb_b = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_b)]:
					bb_b.append(np.asarray(i))							
				self.generator_list_B.extend(bb_b)
						
						#For zeta-C
						
				bb_c = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_c)]:
					bb_c.append(np.asarray(i))		
						
				self.generator_list_C.extend(bb_c)		
						
				
				"""
				Monte-Carlo
				
				"""
						
						#For zeta-A
				
				self.generator_list_A.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(0.3*n_a)])
				self.generator_list_A.extend(list(np.eye(n_rxn)))
				self.generator_list_A.extend(list(-1*np.eye(n_rxn)))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				
						#For zeta-A
						
				self.generator_list_B.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_b)])	
				self.generator_list_B.extend(list(np.eye(2*n_rxn)))
				self.generator_list_B.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
						
						#For zeta-A		
				self.generator_list_C.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_c)])
				self.generator_list_C.extend(list(np.eye(2*n_rxn)))
				self.generator_list_C.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				
				"""
				Uniform
				
				"""
					
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				
				"""
				Latin-Hypercube
				"""
			
						
						#For zeta-A
						
				lhs_a = []
				for i in lhs(n_rxn, samples=int(0.3*n_a)):
					lhs_a.append(np.asarray(i))	
				self.generator_list_A.extend(lhs_a)
				
						
						#For zeta-A
						
				lhs_b = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_b)):
					lhs_b.append(np.asarray(i))		
				self.generator_list_B.extend(lhs_b)
						
						#For zeta-A
						
				lhs_c = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_c)):
					lhs_c.append(np.asarray(i))
				self.generator_list_C.extend(lhs_c)
				
				"""
				Monte-Carlo-Mix
				"""
				
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)		
				#self.generator_list.extend(bb)
				#self.generator_list.extend(list(np.eye(self.n)))
				#self.generator_list.extend(list(-1*np.eye(self.n)))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#print(self.generator_list)
				#self.sim_ = len(self.generator_list)
				
				#trainers = np.asarray(trainers)
				#print(f"Generator_list {np.shape(self.generator_list)}\n")
				df = pd.DataFrame(np.asarray(self.generator_list))
				count_start_a = 0
				count_start_b = 0
				count_start_c = 0			
				for i in self.reactions_selected:
					
					len_of_Arrhenius_A = 1
					count_end_a = count_start_a+len_of_Arrhenius_A
					
					
					len_of_Arrhenius_B = 2
					count_end_b = count_start_b+len_of_Arrhenius_B
					
					
					len_of_Arrhenius_C = 2
					count_end_c = count_start_c+len_of_Arrhenius_C
					
					#print(count_start)
					#print(count_end)
					self.rxn_generators_a2 = []
					self.rxn_generators_b2 = []
					self.rxn_generators_c2 = []
					
					"""
					getting zetas from generators
					"""
					for j in range(len(self.generator_list_A)):
						self.rxn_generators_a2.append(self.generator_list_A[j][count_start_a:count_end_a])
					#print(self.rxn_generators_a2)
					for j in range(len(self.generator_list_B)):
						self.rxn_generators_b2.append(self.generator_list_B[j][count_start_b:count_end_b])
					
					for j in range(len(self.generator_list_C)):
						self.rxn_generators_c2.append(self.generator_list_C[j][count_start_c:count_end_c])
					
					#self.sim_ = len(self.rxn_generators_a2)+len(self.rxn_generators_b2)+len(self.rxn_generators_c2)
					self.sim_ = len(self.rxn_generators_a2)+len(self.rxn_generators_b2)
					#print(self.rxn_generators)
					#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
										
					g_a = self.rxn_generators_a2
					g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					#g_b = self.rxn_generators[0:n_b]
					#g_c = self.rxn_generators[n_b:n_b+n_c]
					
					
					count_start_a = count_end_a
					count_start_b = count_end_b
					count_start_c = count_end_c
					
					"""
					Getting the uncertainty data
					"""
					#fig = plt.figure()
					data = self.unsrt[i].data
					T = np.linspace(300,2500,8)
					kappa_max = self.unsrt[i].getKappaMax(T)			
					kappa_min = self.unsrt[i].getKappaMin(T)
					Theta = np.array([T/T,np.log(T),-1/T])
					P = self.unsrt[i].getMean()	
					cov = self.unsrt[i].getCov()
					zeta_A = self.unsrt[i].zeta.x
					
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					"""
					Type-A zeta samples
					"""
					generators_a = g_a
					parallel_zetas_a = []
					
					time_start = time.time()
					
					
					for k in g_a:
						#Pint = P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten()
					#	plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
					#	plt.plot()
						#print(i)
						#print(type(i))
						#print(k)	
						#parallel_zetas_a.append(P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten())
						parallel_zetas_a.append(k[0]*zeta_A)
										
					#print(parallel_zetas_a)
					#plt.savefig(str(i)+"_a.pdf")
					#print(os.getcwd())
					#print("plot saved")
					data["generators_b"] = g_b
					callWorkForce = Worker(self.allowed_count)	
					generators_b,parallel_zetas_b = callWorkForce.do_unsrt_b(data,len(g_b))
					del callWorkForce
					"""
					data["generators_c"] = g_c
					callWorkForce = Worker(self.allowed_count)	
					generators_c,parallel_zetas_c = callWorkForce.do_unsrt_c(data,len(g_c))
					del callWorkForce
					"""
					#print(parallel_zetas_a)
					
					temp_a[i] = parallel_zetas_a
					gen_a[i] = generators_a
					temp_b[i] = parallel_zetas_b
					gen_b[i] = generators_b
					#temp_c[i] = parallel_zetas_c
					#gen_c[i] = generators_c
					
					time_end = time.time()
					print(f"\n\t\tTime taken = {time_end-time_start}\n")
					
					#print(parallel_zetas_b)
					Y_a = parallel_zetas_a
					Y_b = parallel_zetas_b
					#Y_c = parallel_zetas_c
					T_ = np.linspace(300,2500,8)
					theta = np.array([T/T,np.log(T),-1/T])
					kappa_mean = theta.T.dot(P)
					Theta = np.array([T/T,np.log(T),-1/T])
					X = []
					Y = []
					
					
					for j in Y_a:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					
					for j in Y_b:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_b.pdf")
					
					
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					"""
					for j in Y_c:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_c.pdf")
					"""
					temp_xdata = {}
					for k in range(len(Y)):
						temp_xdata[str(Y[k])] = X[k]
						
					random.shuffle(Y)
					#print(Y_)
					X_ = []
					for k in range(len(Y)):
						X_.append(temp_xdata[str(Y[k])])
						
					#print(X_[9])
					xdata[i] = X_
					ydata[i] = Y
					#print(len(parallel_zetas))
					#print(len(self.rxn_generators))
				A = []
				A_gen = []
				X_dash = []
				y_dash = []
				
				for j in range(len(g_a)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_a:
						a_row.extend(list(temp_a[i][j]))
						a_row_gen.extend(list(gen_a[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				
				
				for j in range(len(g_b)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)]))
						y_temp.extend(list(ydata[i][j+len(g_a)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				
				
				"""
				for j in range(len(g_b)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				"""
				
				for j in range(len(g_c)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)+len(g_b)]))
						y_temp.extend(list(ydata[i][j+len(g_a)+len(g_b)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				"""
				for j in range(len(g_c)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_b)]))
						y_temp.extend(list(ydata[i][j+len(g_b)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
					
				"""
				
				#print(A)
				#self.sim_ = len(g_a)+len(g_b)+len(g_c)
				#self.sim_ = n_b+n_c
				self.beta_ = np.asarray(A)
				#print(self.beta_)
				self.generators = np.asarray(A_gen)
				self.X = np.asarray(X_dash)
				self.y = np.asarray(y_dash)
			
			elif self.design == "B1-only":
				self.beta_=[]
				temp_a = {}
				temp_b = {}
				temp_c = {}
				gen = {}
				gen_a = {}
				gen_b = {}
				gen_c = {}
				
				xdata = {}
				ydata = {}				
				self.generator_list = []
				"""
				Min- simulation = 300
				"""
				#if self.sim_ < 300:
				#	self.sim_ = 300
				#print(self.sim_)
				"""
				Dividing the generators into three bins
				------------------------------------
				n_a = 0.33*self.sim_
				n_b = 0.33*self.sim_
				n_c = 0.34*self.sim_
				
				"""
				n_a = int(0.5*self.sim_)
				#n_b = int(0.5*self.sim_)
				#n_c = self.sim_-n_a-n_b
				
				
				n_b = int(self.sim_)
				n_c = int(0.33*self.sim_)
				
				
				#for i in np.eye(self.n):
				#	self.generator_list.append(i)
				#	self.generator_list.append(-1*i)
								
				"""
				Mixing three types of the extreme curves
				--------------------
				zeta-A: 1
				zeta-B: 1
				zeta-C: 1
			
				for zeta-A: one normal random parameter is required
					a matrix on n_a x n will be created
					
				for zeta-B and zeta-c: two random parameters is required
					we can arrange the two normal_random_parameter 				
					a matrix of n_b+n_c x 2*n will be created
					
				The number of reactions = n_R
				"""
				n_rxn = len(self.ind)
				self.generator_list_A = []
				self.generator_list_B = []
				self.generator_list_C = []
				
				"""
				Factorial Design
				----------------
				box-benkhen design
				bb = []
				for i in bbdesign(n_rxn)[0:int(0.3*self.sim_)]:
					bb.append(np.asarray(i))
				"""
						#For zeta-A
						
				bb_a = []
				for i in bbdesign(n_rxn)[0:int(0.4*n_a)]:
					bb_a.append(np.asarray(i))
				self.generator_list_A.extend(bb_a)
						
						#For zeta-B
						
				bb_b = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_b)]:
					bb_b.append(np.asarray(i))							
				self.generator_list_B.extend(bb_b)
						
						#For zeta-C
						
				bb_c = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_c)]:
					bb_c.append(np.asarray(i))		
						
				self.generator_list_C.extend(bb_c)		
						
				
				"""
				Monte-Carlo
				
				"""
						
						#For zeta-A
				
				self.generator_list_A.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(0.3*n_a)])
				self.generator_list_A.extend(list(np.eye(n_rxn)))
				self.generator_list_A.extend(list(-1*np.eye(n_rxn)))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				
						#For zeta-A
						
				self.generator_list_B.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_b)])	
				self.generator_list_B.extend(list(np.eye(2*n_rxn)))
				self.generator_list_B.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
						
						#For zeta-A		
				self.generator_list_C.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_c)])
				self.generator_list_C.extend(list(np.eye(2*n_rxn)))
				self.generator_list_C.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				
				"""
				Uniform
				
				"""
					
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				
				"""
				Latin-Hypercube
				"""
			
						
						#For zeta-A
						
				lhs_a = []
				for i in lhs(n_rxn, samples=int(0.3*n_a)):
					lhs_a.append(np.asarray(i))	
				self.generator_list_A.extend(lhs_a)
				
						
						#For zeta-A
						
				lhs_b = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_b)):
					lhs_b.append(np.asarray(i))		
				self.generator_list_B.extend(lhs_b)
						
						#For zeta-A
						
				lhs_c = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_c)):
					lhs_c.append(np.asarray(i))
				self.generator_list_C.extend(lhs_c)
				
				"""
				Monte-Carlo-Mix
				"""
				
				
				df = pd.DataFrame(np.asarray(self.generator_list))
				count_start_a = 0
				count_start_b = 0
				count_start_c = 0			
				for i in self.reactions_selected:
					
					len_of_Arrhenius_A = 1
					count_end_a = count_start_a+len_of_Arrhenius_A
					
					
					len_of_Arrhenius_B = 2
					count_end_b = count_start_b+len_of_Arrhenius_B
					
					
					len_of_Arrhenius_C = 2
					count_end_c = count_start_c+len_of_Arrhenius_C
					
					#print(count_start)
					#print(count_end)
					self.rxn_generators_a2 = []
					self.rxn_generators_b2 = []
					self.rxn_generators_c2 = []
					
					"""
					getting zetas from generators
					"""
					for j in range(len(self.generator_list_A)):
						self.rxn_generators_a2.append(self.generator_list_A[j][count_start_a:count_end_a])
					#print(self.rxn_generators_a2)
					for j in range(len(self.generator_list_B)):
						self.rxn_generators_b2.append(self.generator_list_B[j][count_start_b:count_end_b])
					
					for j in range(len(self.generator_list_C)):
						self.rxn_generators_c2.append(self.generator_list_C[j][count_start_c:count_end_c])
					
					#self.sim_ = len(self.rxn_generators_a2)+len(self.rxn_generators_b2)+len(self.rxn_generators_c2)
					self.sim_ = len(self.rxn_generators_b2)
					#print(self.rxn_generators)
					#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
										
					g_a = self.rxn_generators_a2
					g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					#g_b = self.rxn_generators[0:n_b]
					#g_c = self.rxn_generators[n_b:n_b+n_c]
					
					
					count_start_a = count_end_a
					count_start_b = count_end_b
					count_start_c = count_end_c
					
					"""
					Getting the uncertainty data
					"""
					#fig = plt.figure()
					data = self.unsrt[i].data
					T = np.linspace(300,2500,8)
					kappa_max = self.unsrt[i].getKappaMax(T)			
					kappa_min = self.unsrt[i].getKappaMin(T)
					Theta = np.array([T/T,np.log(T),-1/T])
					P = self.unsrt[i].getMean()	
					cov = self.unsrt[i].getCov()
					zeta_A = self.unsrt[i].zeta.x
					
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					"""
					Type-A zeta samples
					"""
					generators_a = g_a
					parallel_zetas_a = []
					
					time_start = time.time()
					
					
					for k in g_a:
						
						parallel_zetas_a.append(k[0]*zeta_A)
										
					data["generators_b"] = g_b
					callWorkForce = Worker(self.allowed_count)	
					generators_b,parallel_zetas_b = callWorkForce.do_unsrt_b(data,len(g_b))
					del callWorkForce
					"""
					data["generators_c"] = g_c
					callWorkForce = Worker(self.allowed_count)	
					generators_c,parallel_zetas_c = callWorkForce.do_unsrt_c(data,len(g_c))
					del callWorkForce
					"""
					#print(parallel_zetas_a)
					
					temp_a[i] = parallel_zetas_a
					gen_a[i] = generators_a
					temp_b[i] = parallel_zetas_b
					gen_b[i] = generators_b
					#temp_c[i] = parallel_zetas_c
					#gen_c[i] = generators_c
					
					time_end = time.time()
					print(f"\n\t\tTime taken = {time_end-time_start}\n")
					
					#print(parallel_zetas_b)
					Y_a = parallel_zetas_a
					Y_b = parallel_zetas_b
					#Y_c = parallel_zetas_c
					T_ = np.linspace(300,2500,8)
					theta = np.array([T/T,np.log(T),-1/T])
					kappa_mean = theta.T.dot(P)
					Theta = np.array([T/T,np.log(T),-1/T])
					X = []
					Y = []
					
					for j in Y_b:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						
					
					temp_xdata = {}
					for k in range(len(Y)):
						temp_xdata[str(Y[k])] = X[k]
						
					random.shuffle(Y)
					#print(Y_)
					X_ = []
					for k in range(len(Y)):
						X_.append(temp_xdata[str(Y[k])])
						
					#print(X_[9])
					xdata[i] = X_
					ydata[i] = Y
					#print(len(parallel_zetas))
					#print(len(self.rxn_generators))
				A = []
				A_gen = []
				X_dash = []
				y_dash = []
								
				for j in range(len(g_b)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				
				
				
				self.beta_ = np.asarray(A)
				#print(self.beta_)
				self.generators = np.asarray(A_gen)
				self.X = np.asarray(X_dash)
				self.y = np.asarray(y_dash)
			
			elif self.design == "C1-only":
				self.beta_=[]
				temp_a = {}
				temp_b = {}
				temp_c = {}
				gen = {}
				gen_a = {}
				gen_b = {}
				gen_c = {}
				
				xdata = {}
				ydata = {}				
				self.generator_list = []
				"""
				Min- simulation = 300
				"""
				#if self.sim_ < 300:
				#	self.sim_ = 300
				#print(self.sim_)
				"""
				Dividing the generators into three bins
				------------------------------------
				n_a = 0.33*self.sim_
				n_b = 0.33*self.sim_
				n_c = 0.34*self.sim_
				
				"""
				n_a = int(0.5*self.sim_)
				#n_b = int(0.5*self.sim_)
				#n_c = self.sim_-n_a-n_b
				
				
				n_b = int(self.sim_)
				n_c = self.sim_
				
				
				#for i in np.eye(self.n):
				#	self.generator_list.append(i)
				#	self.generator_list.append(-1*i)
								
				"""
				Mixing three types of the extreme curves
				--------------------
				zeta-A: 1
				zeta-B: 1
				zeta-C: 1
			
				for zeta-A: one normal random parameter is required
					a matrix on n_a x n will be created
					
				for zeta-B and zeta-c: two random parameters is required
					we can arrange the two normal_random_parameter 				
					a matrix of n_b+n_c x 2*n will be created
					
				The number of reactions = n_R
				"""
				n_rxn = len(self.ind)
				self.generator_list_A = []
				self.generator_list_B = []
				self.generator_list_C = []
				
				"""
				Factorial Design
				----------------
				box-benkhen design
				bb = []
				for i in bbdesign(n_rxn)[0:int(0.3*self.sim_)]:
					bb.append(np.asarray(i))
				"""
						#For zeta-A
						
				bb_a = []
				for i in bbdesign(n_rxn)[0:int(0.4*n_a)]:
					bb_a.append(np.asarray(i))
				self.generator_list_A.extend(bb_a)
						
						#For zeta-B
						
				bb_b = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_b)]:
					bb_b.append(np.asarray(i))							
				self.generator_list_B.extend(bb_b)
						
						#For zeta-C
						
				bb_c = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_c)]:
					bb_c.append(np.asarray(i))		
						
				self.generator_list_C.extend(bb_c)		
						
				
				"""
				Monte-Carlo
				
				"""
						
						#For zeta-A
				
				self.generator_list_A.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(0.3*n_a)])
				self.generator_list_A.extend(list(np.eye(n_rxn)))
				self.generator_list_A.extend(list(-1*np.eye(n_rxn)))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				
						#For zeta-A
						
				self.generator_list_B.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_b)])	
				self.generator_list_B.extend(list(np.eye(2*n_rxn)))
				self.generator_list_B.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
						
						#For zeta-A		
				self.generator_list_C.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_c)])
				self.generator_list_C.extend(list(np.eye(2*n_rxn)))
				self.generator_list_C.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				
				"""
				Uniform
				
				"""
					
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				
				"""
				Latin-Hypercube
				"""
			
						
						#For zeta-A
						
				lhs_a = []
				for i in lhs(n_rxn, samples=int(0.3*n_a)):
					lhs_a.append(np.asarray(i))	
				self.generator_list_A.extend(lhs_a)
				
						
						#For zeta-A
						
				lhs_b = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_b)):
					lhs_b.append(np.asarray(i))		
				self.generator_list_B.extend(lhs_b)
						
						#For zeta-A
						
				lhs_c = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_c)):
					lhs_c.append(np.asarray(i))
				self.generator_list_C.extend(lhs_c)
				
				"""
				Monte-Carlo-Mix
				"""
				
				
				df = pd.DataFrame(np.asarray(self.generator_list))
				count_start_a = 0
				count_start_b = 0
				count_start_c = 0			
				for i in self.reactions_selected:
					
					len_of_Arrhenius_A = 1
					count_end_a = count_start_a+len_of_Arrhenius_A
					
					
					len_of_Arrhenius_B = 2
					count_end_b = count_start_b+len_of_Arrhenius_B
					
					
					len_of_Arrhenius_C = 2
					count_end_c = count_start_c+len_of_Arrhenius_C
					
					#print(count_start)
					#print(count_end)
					self.rxn_generators_a2 = []
					self.rxn_generators_b2 = []
					self.rxn_generators_c2 = []
					
					"""
					getting zetas from generators
					"""
					for j in range(len(self.generator_list_A)):
						self.rxn_generators_a2.append(self.generator_list_A[j][count_start_a:count_end_a])
					#print(self.rxn_generators_a2)
					for j in range(len(self.generator_list_B)):
						self.rxn_generators_b2.append(self.generator_list_B[j][count_start_b:count_end_b])
					
					for j in range(len(self.generator_list_C)):
						self.rxn_generators_c2.append(self.generator_list_C[j][count_start_c:count_end_c])
					
					#self.sim_ = len(self.rxn_generators_a2)+len(self.rxn_generators_b2)+len(self.rxn_generators_c2)
					self.sim_ = len(self.rxn_generators_c2)
					#print(self.rxn_generators)
					#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
										
					g_a = self.rxn_generators_a2
					g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					#g_b = self.rxn_generators[0:n_b]
					#g_c = self.rxn_generators[n_b:n_b+n_c]
					
					
					count_start_a = count_end_a
					count_start_b = count_end_b
					count_start_c = count_end_c
					
					"""
					Getting the uncertainty data
					"""
					#fig = plt.figure()
					data = self.unsrt[i].data
					T = np.linspace(300,2500,8)
					kappa_max = self.unsrt[i].getKappaMax(T)			
					kappa_min = self.unsrt[i].getKappaMin(T)
					Theta = np.array([T/T,np.log(T),-1/T])
					P = self.unsrt[i].getMean()	
					cov = self.unsrt[i].getCov()
					zeta_A = self.unsrt[i].zeta.x
					
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					"""
					Type-A zeta samples
					"""
					generators_a = g_a
					parallel_zetas_a = []
					
					time_start = time.time()
					
					
					for k in g_a:
						
						parallel_zetas_a.append(k[0]*zeta_A)
					generators_b = []
					parallel_zetas_b = []
					generators_c = []
					parallel_zetas_c = []
					start = 0
					for w in range(10):
						stop = int(((w+1)/10)*self.sim_)
						#g_b = self.rxn_generators_b2[start:stop]
						g_c = self.rxn_generators_c2[start:stop]
						data["generators_c"] = g_c
						callWorkForce = Worker(self.allowed_count)	
						generators_c_,parallel_zetas_c_ = callWorkForce.do_unsrt_c(data,len(g_c))
						del callWorkForce
						
						generators_c.extend(generators_c_)
						parallel_zetas_c.extend(parallel_zetas_c_)
						removed = data.pop("generators_c")
						start = stop
					
					#print(self.sim_)
					#print(len(parallel_zetas_b)+len(parallel_zetas_c))
					
					#g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					
					#print(parallel_zetas_a)
					
					temp_a[i] = parallel_zetas_a
					gen_a[i] = generators_a
					#temp_b[i] = parallel_zetas_b
					#gen_b[i] = generators_b
					temp_c[i] = parallel_zetas_c
					gen_c[i] = generators_c
					
					time_end = time.time()
					print(f"\n\t\tTime taken = {time_end-time_start}\n")
					
					#print(parallel_zetas_b)
					Y_a = parallel_zetas_a
					#Y_b = parallel_zetas_b
					Y_c = parallel_zetas_c
					T_ = np.linspace(300,2500,8)
					theta = np.array([T/T,np.log(T),-1/T])
					kappa_mean = theta.T.dot(P)
					Theta = np.array([T/T,np.log(T),-1/T])
					X = []
					Y = []
					
					for j in Y_c:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						
					
					temp_xdata = {}
					for k in range(len(Y)):
						temp_xdata[str(Y[k])] = X[k]
						
					random.shuffle(Y)
					#print(Y_)
					X_ = []
					for k in range(len(Y)):
						X_.append(temp_xdata[str(Y[k])])
						
					#print(X_[9])
					xdata[i] = X_
					ydata[i] = Y
					#print(len(parallel_zetas))
					#print(len(self.rxn_generators))
				A = []
				A_gen = []
				X_dash = []
				y_dash = []
								
				for j in range(len(g_c)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				
				
				
				self.beta_ = np.asarray(A)
				#print(self.beta_)
				self.generators = np.asarray(A_gen)
				self.X = np.asarray(X_dash)
				self.y = np.asarray(y_dash)
			
			
			
			elif self.design == "A1+B1+C1":
				self.beta_=[]
				temp_a = {}
				temp_b = {}
				temp_c = {}
				gen = {}
				gen_a = {}
				gen_b = {}
				gen_c = {}
				
				xdata = {}
				ydata = {}				
				self.generator_list = []
				"""
				Min- simulation = 300
				"""
				#if self.sim_ < 300:
				#	self.sim_ = 300
				#print(self.sim_)
				"""
				Dividing the generators into three bins
				------------------------------------
				n_a = 0.33*self.sim_
				n_b = 0.33*self.sim_
				n_c = 0.34*self.sim_
				
				"""
				n_a = int(0.1*self.sim_)
				n_b = int(0.45*self.sim_)
				n_c = self.sim_-n_a-n_b
				#print(n_a)
				#print(n_b)
				#print(n_c)
				
				#n_b = int(0.5*self.sim_)
				#n_c = self.sim_-n_b
				
				
				#for i in np.eye(self.n):
				#	self.generator_list.append(i)
				#	self.generator_list.append(-1*i)
								
				"""
				Mixing three types of the extreme curves
				--------------------
				zeta-A: 1
				zeta-B: 1
				zeta-C: 1
			
				for zeta-A: one normal random parameter is required
					a matrix on n_a x n will be created
					
				for zeta-B and zeta-c: two random parameters is required
					we can arrange the two normal_random_parameter 				
					a matrix of n_b+n_c x 2*n will be created
					
				The number of reactions = n_R
				"""
				n_rxn = len(self.ind)
				self.generator_list_A = []
				self.generator_list_B = []
				self.generator_list_C = []
				
				"""
				Factorial Design
				----------------
				box-benkhen design
				bb = []
				for i in bbdesign(n_rxn)[0:int(0.3*self.sim_)]:
					bb.append(np.asarray(i))
				"""
						#For zeta-A
						
				bb_a = []
				for i in bbdesign(n_rxn)[0:int(0.4*n_a)]:
					bb_a.append(np.asarray(i))
				self.generator_list_A.extend(bb_a)
						
						#For zeta-B
						
				bb_b = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_b)]:
					bb_b.append(np.asarray(i))							
				self.generator_list_B.extend(bb_b)
						
						#For zeta-C
						
				bb_c = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_c)]:
					bb_c.append(np.asarray(i))		
						
				self.generator_list_C.extend(bb_c)		
						
				
				"""
				Monte-Carlo
				
				"""
						
						#For zeta-A
				
				self.generator_list_A.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(0.3*n_a)])
				self.generator_list_A.extend(list(np.eye(n_rxn)))
				self.generator_list_A.extend(list(-1*np.eye(n_rxn)))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				
						#For zeta-A
						
				self.generator_list_B.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_b)])	
				self.generator_list_B.extend(list(np.eye(2*n_rxn)))
				self.generator_list_B.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
						
						#For zeta-A		
				self.generator_list_C.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_c)])
				self.generator_list_C.extend(list(np.eye(2*n_rxn)))
				self.generator_list_C.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				
				"""
				Uniform
				
				"""
					
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				
				"""
				Latin-Hypercube
				"""
			
						
						#For zeta-A
						
				lhs_a = []
				for i in lhs(n_rxn, samples=int(0.3*n_a)):
					lhs_a.append(np.asarray(i))	
				self.generator_list_A.extend(lhs_a)
				
						
						#For zeta-A
						
				lhs_b = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_b)):
					lhs_b.append(np.asarray(i))		
				self.generator_list_B.extend(lhs_b)
						
						#For zeta-A
						
				lhs_c = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_c)):
					lhs_c.append(np.asarray(i))
				self.generator_list_C.extend(lhs_c)
				
				"""
				Monte-Carlo-Mix
				"""
				
				
				
				
				
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)		
				#self.generator_list.extend(bb)
				#self.generator_list.extend(list(np.eye(self.n)))
				#self.generator_list.extend(list(-1*np.eye(self.n)))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#print(self.generator_list)
				#self.sim_ = len(self.generator_list)
				
				#trainers = np.asarray(trainers)
				#print(f"Generator_list {np.shape(self.generator_list)}\n")
				df = pd.DataFrame(np.asarray(self.generator_list))
				count_start_a = 0
				count_start_b = 0
				count_start_c = 0			
				for i in self.reactions_selected:
					
					len_of_Arrhenius_A = 1
					count_end_a = count_start_a+len_of_Arrhenius_A
					
					
					len_of_Arrhenius_B = 2
					count_end_b = count_start_b+len_of_Arrhenius_B
					
					
					len_of_Arrhenius_C = 2
					count_end_c = count_start_c+len_of_Arrhenius_C
					
					#print(count_start)
					#print(count_end)
					self.rxn_generators_a2 = []
					self.rxn_generators_b2 = []
					self.rxn_generators_c2 = []
					
					"""
					getting zetas from generators
					"""
					for j in range(len(self.generator_list_A)):
						self.rxn_generators_a2.append(self.generator_list_A[j][count_start_a:count_end_a])
					#print(self.rxn_generators_a2)
					for j in range(len(self.generator_list_B)):
						self.rxn_generators_b2.append(self.generator_list_B[j][count_start_b:count_end_b])
					#print(self.rxn_generators_b2)
					
					for j in range(len(self.generator_list_C)):
						self.rxn_generators_c2.append(self.generator_list_C[j][count_start_c:count_end_c])
					#print(self.rxn_generators_c2)
					self.sim_ = len(self.rxn_generators_a2)+len(self.rxn_generators_b2)+len(self.rxn_generators_c2)
					#print(self.rxn_generators)
					#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
									
					g_a = self.rxn_generators_a2
					g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					#g_b = self.rxn_generators[0:n_b]
					#g_c = self.rxn_generators[n_b:n_b+n_c]
					#print(len(g_b))
					#print(len(g_c))	
					
					count_start_a = count_end_a
					count_start_b = count_end_b
					count_start_c = count_end_c
					
					"""
					Getting the uncertainty data
					"""
					#fig = plt.figure()
					data = self.unsrt[i].data
					T = np.linspace(300,2500,8)
					kappa_max = self.unsrt[i].getKappaMax(T)			
					kappa_min = self.unsrt[i].getKappaMin(T)
					Theta = np.array([T/T,np.log(T),-1/T])
					P = self.unsrt[i].getMean()	
					cov = self.unsrt[i].getCov()
					zeta_A = self.unsrt[i].zeta.x
					
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					"""
					Type-A zeta samples
					"""
					generators_a = g_a
					parallel_zetas_a = []
					
					time_start = time.time()
					
					
					for k in g_a:
						#Pint = P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten()
					#	plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
					#	plt.plot()
						#print(i)
						#print(type(i))
						#print(k)	
						#parallel_zetas_a.append(P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten())
						parallel_zetas_a.append(k[0]*zeta_A)
										
					#print(parallel_zetas_a)
					#plt.savefig(str(i)+"_a.pdf")
					#print(os.getcwd())
					#print("plot saved")
					
					generators_b = []
					parallel_zetas_b = []
					generators_c = []
					parallel_zetas_c = []
					start = 0
					for w in range(10):
						stop = int(((w+1)/10)*self.sim_)
						g_b = self.rxn_generators_b2[start:stop]
						g_c = self.rxn_generators_c2[start:stop]
						data["generators_b"] = g_b
						callWorkForce = Worker(self.allowed_count)	
						generators_b_,parallel_zetas_b_ = callWorkForce.do_unsrt_b(data,len(g_b))
						del callWorkForce
						
						
						removed = data.pop("generators_b")
						data["generators_c"] = g_c
						callWorkForce = Worker(self.allowed_count)	
						generators_c_,parallel_zetas_c_ = callWorkForce.do_unsrt_c(data,len(g_c))
						del callWorkForce
						
						generators_b.extend(generators_b_)
						parallel_zetas_b.extend(parallel_zetas_b_)
						generators_c.extend(generators_c_)
						parallel_zetas_c.extend(parallel_zetas_c_)
						removed = data.pop("generators_c")
						start = stop
					
					#print(self.sim_)
					#print(len(parallel_zetas_b)+len(parallel_zetas_c))
					g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					#raise AssertionError("Stop!")
					#print(parallel_zetas_a)
					
					temp_a[i] = parallel_zetas_a
					gen_a[i] = generators_a
					temp_b[i] = parallel_zetas_b
					gen_b[i] = generators_b
					temp_c[i] = parallel_zetas_c
					gen_c[i] = generators_c
					
					time_end = time.time()
					print(f"\n\t\tTime taken = {time_end-time_start}\n")
					
					#print(parallel_zetas_b)
					Y_a = parallel_zetas_a
					Y_b = parallel_zetas_b
					Y_c = parallel_zetas_c
					T_ = np.linspace(300,2500,8)
					theta = np.array([T/T,np.log(T),-1/T])
					kappa_mean = theta.T.dot(P)
					Theta = np.array([T/T,np.log(T),-1/T])
					X = []
					Y = []
					
					
					for j in Y_a:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					
					for j in Y_b:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_b.pdf")
					
					
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					for j in Y_c:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_c.pdf")
					
					temp_xdata = {}
					for k in range(len(Y)):
						temp_xdata[str(Y[k])] = X[k]
						
					random.shuffle(Y)
					#print(Y_)
					X_ = []
					for k in range(len(Y)):
						X_.append(temp_xdata[str(Y[k])])
						
					#print(X_[9])
					xdata[i] = X_
					ydata[i] = Y
					#print(len(parallel_zetas))
					#print(len(self.rxn_generators))
				A = []
				A_gen = []
				X_dash = []
				y_dash = []
				
				for j in range(len(g_a)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_a:
						a_row.extend(list(temp_a[i][j]))
						a_row_gen.extend(list(gen_a[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				
				
				for j in range(len(g_b)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)]))
						y_temp.extend(list(ydata[i][j+len(g_a)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				for j in range(n_b):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				
				for j in range(len(g_c)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)+len(g_b)]))
						y_temp.extend(list(ydata[i][j+len(g_a)+len(g_b)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				for j in range(n_c):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+n_b]))
						y_temp.extend(list(ydata[i][j+n_b]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
					
				"""
				
				#print(A)
				#self.sim_ = len(g_a)+len(g_b)+len(g_c)
				#self.sim_ = n_b+n_c
				self.beta_ = np.asarray(A)
				#print(len(self.beta_))
				#print(self.beta_)
				self.generators = np.asarray(A_gen)
				self.X = np.asarray(X_dash)
				self.y = np.asarray(y_dash)
			
			elif self.design == "B1-factorial":
				self.beta_=[]
				temp_a = {}
				temp_b = {}
				temp_c = {}
				gen = {}
				gen_a = {}
				gen_b = {}
				gen_c = {}
				
				xdata = {}
				ydata = {}				
				self.generator_list = []
				"""
				Min- simulation = 300
				"""
				#if self.sim_ < 300:
				#	self.sim_ = 300
				#print(self.sim_)
				"""
				Dividing the generators into three bins
				------------------------------------
				n_a = 0.33*self.sim_
				n_b = 0.33*self.sim_
				n_c = 0.34*self.sim_
				
				"""
				n_a = int(0.33*self.sim_)
				n_b = int(0.33*self.sim_)
				n_c = self.sim_-n_a-n_b
				
				
				#n_b = int(0.5*self.sim_)
				#n_c = self.sim_-n_b
				
				
				#for i in np.eye(self.n):
				#	self.generator_list.append(i)
				#	self.generator_list.append(-1*i)
								
				"""
				Mixing three types of the extreme curves
				--------------------
				zeta-A: 1
				zeta-B: 1
				zeta-C: 1
			
				for zeta-A: one normal random parameter is required
					a matrix on n_a x n will be created
					
				for zeta-B and zeta-c: two random parameters is required
					we can arrange the two normal_random_parameter 				
					a matrix of n_b+n_c x 2*n will be created
					
				The number of reactions = n_R
				"""
				n_rxn = len(self.ind)
				self.generator_list_A = []
				self.generator_list_B = []
				self.generator_list_C = []
				
				"""
				Factorial Design
				----------------
				box-benkhen design
				bb = []
				for i in bbdesign(n_rxn)[0:int(0.3*self.sim_)]:
					bb.append(np.asarray(i))
				"""
						#For zeta-A
						
				bb_a = []
				for i in bbdesign(n_rxn):
					bb_a.append(np.asarray(i))
				self.generator_list_A.extend(bb_a)
				self.generator_list_A.extend(list(np.eye(n_rxn)))
				self.generator_list_A.extend(list(-1*np.eye(n_rxn)))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))		
						#For zeta-B
						
				bb_b = []
				for i in bbdesign(2*n_rxn):
					bb_b.append(np.asarray(i))							
				self.generator_list_B.extend(bb_b)
				self.generator_list_B.extend(list(np.eye(2*n_rxn)))
				self.generator_list_B.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))		
						#For zeta-C
						
				bb_c = []
				for i in bbdesign(2*n_rxn):
					bb_c.append(np.asarray(i))		
						
				self.generator_list_C.extend(bb_c)		
				self.generator_list_C.extend(list(np.eye(2*n_rxn)))
				self.generator_list_C.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))		
				
				"""
				Monte-Carlo
				
				"""
				"""		
						#For zeta-A
				
				self.generator_list_A.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(0.3*n_a)])
				
						#For zeta-A
						
				self.generator_list_B.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_a)])	
						
						#For zeta-A		
				self.generator_list_C.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_a)])
				
				"""
				"""
				Uniform
				
				"""
					
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				
				"""
				Latin-Hypercube
				"""
				"""
						
						#For zeta-A
						
				lhs_a = []
				for i in lhs(n_rxn, samples=int(0.3*n_a)):
					lhs_a.append(np.asarray(i))	
				self.generator_list_A.extend(lhs_a)
				
						
						#For zeta-A
						
				lhs_b = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_b)):
					lhs_b.append(np.asarray(i))		
				self.generator_list_B.extend(lhs_b)
						
						#For zeta-A
						
				lhs_c = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_c)):
					lhs_c.append(np.asarray(i))
				self.generator_list_C.extend(lhs_c)
				
				"""
				"""
				Monte-Carlo-Mix
				"""
				
				
				
				
				
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)		
				#self.generator_list.extend(bb)
				#self.generator_list.extend(list(np.eye(self.n)))
				#self.generator_list.extend(list(-1*np.eye(self.n)))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#print(self.generator_list)
				#self.sim_ = len(self.generator_list)
				
				#trainers = np.asarray(trainers)
				#print(f"Generator_list {np.shape(self.generator_list)}\n")
				df = pd.DataFrame(np.asarray(self.generator_list))
				count_start_a = 0
				count_start_b = 0
				count_start_c = 0			
				for i in self.reactions_selected:
					
					len_of_Arrhenius_A = 1
					count_end_a = count_start_a+len_of_Arrhenius_A
					
					
					len_of_Arrhenius_B = 2
					count_end_b = count_start_b+len_of_Arrhenius_B
					
					
					len_of_Arrhenius_C = 2
					count_end_c = count_start_c+len_of_Arrhenius_C
					
					#print(count_start)
					#print(count_end)
					self.rxn_generators_a2 = []
					self.rxn_generators_b2 = []
					self.rxn_generators_c2 = []
					
					"""
					getting zetas from generators
					"""
					for j in range(len(self.generator_list_A)):
						self.rxn_generators_a2.append(self.generator_list_A[j][count_start_a:count_end_a])
					#print(self.rxn_generators_a2)
					for j in range(len(self.generator_list_B)):
						self.rxn_generators_b2.append(self.generator_list_B[j][count_start_b:count_end_b])
					
					for j in range(len(self.generator_list_C)):
						self.rxn_generators_c2.append(self.generator_list_C[j][count_start_c:count_end_c])
					
					self.sim_ = len(self.rxn_generators_a2)+len(self.rxn_generators_b2)+len(self.rxn_generators_c2)
					#print(self.rxn_generators)
					#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
										
					g_a = self.rxn_generators_a2
					g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					#g_b = self.rxn_generators[0:n_b]
					#g_c = self.rxn_generators[n_b:n_b+n_c]
					
					
					count_start_a = count_end_a
					count_start_b = count_end_b
					count_start_c = count_end_c
					
					"""
					Getting the uncertainty data
					"""
					#fig = plt.figure()
					data = self.unsrt[i].data
					T = np.linspace(300,2500,8)
					kappa_max = self.unsrt[i].getKappaMax(T)			
					kappa_min = self.unsrt[i].getKappaMin(T)
					Theta = np.array([T/T,np.log(T),-1/T])
					P = self.unsrt[i].getMean()	
					cov = self.unsrt[i].getCov()
					zeta_A = self.unsrt[i].zeta.x
					
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					"""
					Type-A zeta samples
					"""
					generators_a = g_a
					parallel_zetas_a = []
					
					time_start = time.time()
					
					
					for k in g_a:
						#Pint = P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten()
					#	plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
					#	plt.plot()
						#print(i)
						#print(type(i))
						#print(k)	
						#parallel_zetas_a.append(P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten())
						parallel_zetas_a.append(k[0]*zeta_A)
										
					#print(parallel_zetas_a)
					#plt.savefig(str(i)+"_a.pdf")
					#print(os.getcwd())
					#print("plot saved")
					data["generators_b"] = g_b
					callWorkForce = Worker(self.allowed_count)	
					generators_b,parallel_zetas_b = callWorkForce.do_unsrt_b(data,len(g_b))
					del callWorkForce
					
					data["generators_c"] = g_c
					callWorkForce = Worker(self.allowed_count)	
					generators_c,parallel_zetas_c = callWorkForce.do_unsrt_c(data,len(g_c))
					del callWorkForce
					
					#print(parallel_zetas_a)
					
					temp_a[i] = parallel_zetas_a
					gen_a[i] = generators_a
					temp_b[i] = parallel_zetas_b
					gen_b[i] = generators_b
					temp_c[i] = parallel_zetas_c
					gen_c[i] = generators_c
					
					time_end = time.time()
					print(f"\n\t\tTime taken = {time_end-time_start}\n")
					
					#print(parallel_zetas_b)
					Y_a = parallel_zetas_a
					Y_b = parallel_zetas_b
					Y_c = parallel_zetas_c
					T_ = np.linspace(300,2500,8)
					theta = np.array([T/T,np.log(T),-1/T])
					kappa_mean = theta.T.dot(P)
					Theta = np.array([T/T,np.log(T),-1/T])
					X = []
					Y = []
					
					
					for j in Y_a:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					
					for j in Y_b:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_b.pdf")
					
					
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					for j in Y_c:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_c.pdf")
					
					temp_xdata = {}
					for k in range(len(Y)):
						temp_xdata[str(Y[k])] = X[k]
						
					random.shuffle(Y)
					#print(Y_)
					X_ = []
					for k in range(len(Y)):
						X_.append(temp_xdata[str(Y[k])])
						
					#print(X_[9])
					xdata[i] = X_
					ydata[i] = Y
					#print(len(parallel_zetas))
					#print(len(self.rxn_generators))
				A = []
				A_gen = []
				X_dash = []
				y_dash = []
				
				for j in range(len(g_a)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_a:
						a_row.extend(list(temp_a[i][j]))
						a_row_gen.extend(list(gen_a[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				
				
				for j in range(len(g_b)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)]))
						y_temp.extend(list(ydata[i][j+len(g_a)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				for j in range(n_b):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				
				for j in range(len(g_c)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)+len(g_b)]))
						y_temp.extend(list(ydata[i][j+len(g_a)+len(g_b)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				for j in range(n_c):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+n_b]))
						y_temp.extend(list(ydata[i][j+n_b]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
					
				"""
				
				#print(A)
				self.sim_ = len(g_a)+len(g_b)+len(g_c)
				#self.sim_ = n_b+n_c
				self.beta_ = np.asarray(A)
				#print(self.beta_)
				self.generators = np.asarray(A_gen)
				self.X = np.asarray(X_dash)
				self.y = np.asarray(y_dash)
						
			
			elif self.design == "B1-MC":
				self.beta_=[]
				temp_a = {}
				temp_b = {}
				temp_c = {}
				gen = {}
				gen_a = {}
				gen_b = {}
				gen_c = {}
				
				xdata = {}
				ydata = {}				
				self.generator_list = []
				"""
				Min- simulation = 300
				"""
				#if self.sim_ < 300:
				#	self.sim_ = 300
				#print(self.sim_)
				"""
				Dividing the generators into three bins
				------------------------------------
				n_a = 0.33*self.sim_
				n_b = 0.33*self.sim_
				n_c = 0.34*self.sim_
				
				"""
				n_a = int(0.33*self.sim_)
				n_b = int(0.33*self.sim_)
				n_c = self.sim_-n_a-n_b
				
				
				#n_b = int(0.5*self.sim_)
				#n_c = self.sim_-n_b
				
				
				#for i in np.eye(self.n):
				#	self.generator_list.append(i)
				#	self.generator_list.append(-1*i)
								
				"""
				Mixing three types of the extreme curves
				--------------------
				zeta-A: 1
				zeta-B: 1
				zeta-C: 1
			
				for zeta-A: one normal random parameter is required
					a matrix on n_a x n will be created
					
				for zeta-B and zeta-c: two random parameters is required
					we can arrange the two normal_random_parameter 				
					a matrix of n_b+n_c x 2*n will be created
					
				The number of reactions = n_R
				"""
				n_rxn = len(self.ind)
				self.generator_list_A = []
				self.generator_list_B = []
				self.generator_list_C = []
				
				"""
				Factorial Design
				----------------
				box-benkhen design
				bb = []
				for i in bbdesign(n_rxn)[0:int(0.3*self.sim_)]:
					bb.append(np.asarray(i))
				"""
				"""		#For zeta-A
						
				bb_a = []
				for i in bbdesign(n_rxn)[0:int(0.4*n_a)]:
					bb_a.append(np.asarray(i))
				self.generator_list_A.extend(bb_a)
						
						#For zeta-B
						
				bb_b = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_b)]:
					bb_b.append(np.asarray(i))							
				self.generator_list_B.extend(bb_b)
						
						#For zeta-C
						
				bb_c = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_c)]:
					bb_c.append(np.asarray(i))		
						
				self.generator_list_C.extend(bb_c)		
						
				"""
				"""
				Monte-Carlo
				
				"""
						
						#For zeta-A
				
				self.generator_list_A.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(0.3*n_a)])
				#self.generator_list_A.extend(list(np.eye(n_rxn)))
				#self.generator_list_A.extend(list(-1*np.eye(n_rxn)))
				#self.generator_list_A.append(np.zeros(n_rxn))
				#self.generator_list_A.append(np.zeros(n_rxn))
				#self.generator_list_A.append(np.zeros(n_rxn))
				
						#For zeta-A
						
				self.generator_list_B.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_a)])	
				#self.generator_list_B.extend(list(np.eye(2*n_rxn)))
				#self.generator_list_B.extend(list(-1*np.eye(2*n_rxn)))
				#self.generator_list_B.append(np.zeros(2*n_rxn))
				#self.generator_list_B.append(np.zeros(2*n_rxn))
				#self.generator_list_B.append(np.zeros(2*n_rxn))
						
						#For zeta-A		
				self.generator_list_C.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_a)])
				#self.generator_list_C.extend(list(np.eye(2*n_rxn)))
				#self.generator_list_C.extend(list(-1*np.eye(2*n_rxn)))
				#self.generator_list_C.append(np.zeros(2*n_rxn))
				#self.generator_list_C.append(np.zeros(2*n_rxn))
				#self.generator_list_C.append(np.zeros(2*n_rxn))
				
				"""
				Uniform
				
				"""
					
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				
				"""
				Latin-Hypercube
				"""
				"""
						
						#For zeta-A
						
				lhs_a = []
				for i in lhs(n_rxn, samples=int(0.3*n_a)):
					lhs_a.append(np.asarray(i))	
				self.generator_list_A.extend(lhs_a)
				
						
						#For zeta-A
						
				lhs_b = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_b)):
					lhs_b.append(np.asarray(i))		
				self.generator_list_B.extend(lhs_b)
						
						#For zeta-A
						
				lhs_c = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_c)):
					lhs_c.append(np.asarray(i))
				self.generator_list_C.extend(lhs_c)
				
				"""
				"""
				Monte-Carlo-Mix
				"""
				
				
				
				
				
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)		
				#self.generator_list.extend(bb)
				#self.generator_list.extend(list(np.eye(self.n)))
				#self.generator_list.extend(list(-1*np.eye(self.n)))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#print(self.generator_list)
				#self.sim_ = len(self.generator_list)
				
				#trainers = np.asarray(trainers)
				#print(f"Generator_list {np.shape(self.generator_list)}\n")
				df = pd.DataFrame(np.asarray(self.generator_list))
				count_start_a = 0
				count_start_b = 0
				count_start_c = 0			
				for i in self.reactions_selected:
					
					len_of_Arrhenius_A = 1
					count_end_a = count_start_a+len_of_Arrhenius_A
					
					
					len_of_Arrhenius_B = 2
					count_end_b = count_start_b+len_of_Arrhenius_B
					
					
					len_of_Arrhenius_C = 2
					count_end_c = count_start_c+len_of_Arrhenius_C
					
					#print(count_start)
					#print(count_end)
					self.rxn_generators_a2 = []
					self.rxn_generators_b2 = []
					self.rxn_generators_c2 = []
					
					"""
					getting zetas from generators
					"""
					for j in range(len(self.generator_list_A)):
						self.rxn_generators_a2.append(self.generator_list_A[j][count_start_a:count_end_a])
					#print(self.rxn_generators_a2)
					for j in range(len(self.generator_list_B)):
						self.rxn_generators_b2.append(self.generator_list_B[j][count_start_b:count_end_b])
					
					for j in range(len(self.generator_list_C)):
						self.rxn_generators_c2.append(self.generator_list_C[j][count_start_c:count_end_c])
					
					self.sim_ = len(self.rxn_generators_a2)+len(self.rxn_generators_b2)+len(self.rxn_generators_c2)
					#print(self.rxn_generators)
					#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
										
					g_a = self.rxn_generators_a2
					g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					#g_b = self.rxn_generators[0:n_b]
					#g_c = self.rxn_generators[n_b:n_b+n_c]
					
					
					count_start_a = count_end_a
					count_start_b = count_end_b
					count_start_c = count_end_c
					
					"""
					Getting the uncertainty data
					"""
					#fig = plt.figure()
					data = self.unsrt[i].data
					T = np.linspace(300,2500,8)
					kappa_max = self.unsrt[i].getKappaMax(T)			
					kappa_min = self.unsrt[i].getKappaMin(T)
					Theta = np.array([T/T,np.log(T),-1/T])
					P = self.unsrt[i].getMean()	
					cov = self.unsrt[i].getCov()
					zeta_A = self.unsrt[i].zeta.x
					
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					"""
					Type-A zeta samples
					"""
					generators_a = g_a
					parallel_zetas_a = []
					
					time_start = time.time()
					
					
					for k in g_a:
						#Pint = P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten()
					#	plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
					#	plt.plot()
						#print(i)
						#print(type(i))
						#print(k)	
						#parallel_zetas_a.append(P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten())
						parallel_zetas_a.append(k[0]*zeta_A)
										
					#print(parallel_zetas_a)
					#plt.savefig(str(i)+"_a.pdf")
					#print(os.getcwd())
					#print("plot saved")
					data["generators_b"] = g_b
					callWorkForce = Worker(self.allowed_count)	
					generators_b,parallel_zetas_b = callWorkForce.do_unsrt_b(data,len(g_b))
					del callWorkForce
					
					data["generators_c"] = g_c
					callWorkForce = Worker(self.allowed_count)	
					generators_c,parallel_zetas_c = callWorkForce.do_unsrt_c(data,len(g_c))
					del callWorkForce
					
					#print(parallel_zetas_a)
					
					temp_a[i] = parallel_zetas_a
					gen_a[i] = generators_a
					temp_b[i] = parallel_zetas_b
					gen_b[i] = generators_b
					temp_c[i] = parallel_zetas_c
					gen_c[i] = generators_c
					
					time_end = time.time()
					print(f"\n\t\tTime taken = {time_end-time_start}\n")
					
					#print(parallel_zetas_b)
					Y_a = parallel_zetas_a
					Y_b = parallel_zetas_b
					Y_c = parallel_zetas_c
					T_ = np.linspace(300,2500,8)
					theta = np.array([T/T,np.log(T),-1/T])
					kappa_mean = theta.T.dot(P)
					Theta = np.array([T/T,np.log(T),-1/T])
					X = []
					Y = []
					
					
					for j in Y_a:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					
					for j in Y_b:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_b.pdf")
					
					
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					for j in Y_c:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_c.pdf")
					
					temp_xdata = {}
					for k in range(len(Y)):
						temp_xdata[str(Y[k])] = X[k]
						
					random.shuffle(Y)
					#print(Y_)
					X_ = []
					for k in range(len(Y)):
						X_.append(temp_xdata[str(Y[k])])
						
					#print(X_[9])
					xdata[i] = X_
					ydata[i] = Y
					#print(len(parallel_zetas))
					#print(len(self.rxn_generators))
				A = []
				A_gen = []
				X_dash = []
				y_dash = []
				
				for j in range(len(g_a)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_a:
						a_row.extend(list(temp_a[i][j]))
						a_row_gen.extend(list(gen_a[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				
				
				for j in range(len(g_b)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)]))
						y_temp.extend(list(ydata[i][j+len(g_a)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				for j in range(n_b):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				
				for j in range(len(g_c)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)+len(g_b)]))
						y_temp.extend(list(ydata[i][j+len(g_a)+len(g_b)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				for j in range(n_c):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+n_b]))
						y_temp.extend(list(ydata[i][j+n_b]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
					
				"""
				
				#print(A)
				self.sim_ = len(g_a)+len(g_b)+len(g_c)
				#self.sim_ = n_b+n_c
				self.beta_ = np.asarray(A)
				#print(self.beta_)
				self.generators = np.asarray(A_gen)
				self.X = np.asarray(X_dash)
				self.y = np.asarray(y_dash)
			
			elif self.design == "B1-Linear":
				self.beta_=[]
				temp = {}
				gen = {}
				xdata = {}
				ydata = {}
				for i in self.reactions_selected:
					
					self.generator_list = []
					"""
					Getting linear combination of zeta
					"""
					for j in range(self.sim_):
						beta1 = 2*np.random.random_sample(1)-1
						beta2 = 2*np.random.random_sample(1)-1
						m = beta1*np.pi
						n = beta2*np.pi/2
						alpha = m + np.pi
						gamma = n + (np.pi/2)
						x = np.cos(alpha)*np.sin(gamma)
						y = np.sin(alpha)*np.sin(gamma)
						z = np.cos(gamma)
						self.generator_list.append(np.array([x[0],y[0],z[0]]))
					

					"""
					getting zetas from generators
					"""
					"""
					for j in range(self.sim_):
						self.generator_list.append(2*np.random.random_sample(3)-1)
					"""
					
					data = self.unsrt[i].data
					T = np.linspace(300,2500,8)
					kappa_max = self.unsrt[i].getKappaMax(T)			
					kappa_min = self.unsrt[i].getKappaMin(T)
					P = self.unsrt[i].getMean()	
					zetaMat = self.unsrt[i].zeta_Matrix
					parallel_zetas = []
					for k in range(self.sim_):
						zeta = np.asarray(zetaMat.T.dot(self.generator_list[k].T)).flatten()
						parallel_zetas.append(zeta)
					
					#print(parallel_zetas)
					cov = self.unsrt[i].getCov()
					data["generators"] = self.generator_list
					#callWorkForce = Worker(10)	
					#generators,parallel_zetas = callWorkForce.do_unsrt(data,int(self.sim_))
					gen[i] = self.generator_list
					#del callWorkForce
					temp[i] = parallel_zetas
					
					Y = parallel_zetas
					#print(Y)
					T_ = np.linspace(300,2500,8)
					theta = np.array([T/T,np.log(T),-1/T])
					kappa_mean = theta.T.dot(P)
					Theta = np.array([T_/T_,np.log(T_),-1/T_])
					X = []
					for j in Y:
						
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#Pint = P + np.asarray(np.dot(cov,np.asarray(j))).flatten();
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					
					temp_xdata = {}
					for k in range(len(Y)):
						temp_xdata[str(Y[k])] = X[k]
						
					random.shuffle(Y)
					#print(Y_)
					X_ = []
					for k in range(len(Y)):
						X_.append(temp_xdata[str(Y[k])])
						
					#print(X_[9])
					xdata[i] = X_
					ydata[i] = Y
					
					#print(len(parallel_zetas))
				A = []
				A_gen = []
				X = []
				y = []
				for j in range(self.sim_):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp:
						a_row.extend(list(temp[i][j]))
						a_row_gen.extend(list(gen[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X.append(np.asarray(x_temp))
					y.append(np.asarray(y_temp))
				
				self.beta_ = np.asarray(A)
				#print(self.beta_)
				self.generators = np.asarray(A_gen)
				self.X = np.asarray(X)
				self.y = np.asarray(y)
		
		elif self.simulation == "test":	
			#print(self.n)
			if self.PRS_type == "Partial":
				self.n = len(self.selectedParams)
			#print(self.n)
			if "Sampling_of_PRS" not in self.opt_dict["Stats"]:
				raise AssertionError("Error in opt_dict")
			self.sim_ = int(0.3*self.opt_dict["Stats"]["Sampling_of_PRS"]*self.n_)
			#self.sim_ = 1000		
		# taking 2 as perturbation factor
			self.beta_=[]
			#mean = np.zeros(self.n)
			#cov = 0.33*np.eye(self.n)
			#self.beta_ = list(np.random.multivariate_normal(mean,cov,self.sim_))
			if self.design == 'A-facto': 
				self.beta_=[]
				task = "complete"
				if "Data" in os.listdir():
					os.chdir("Data/Simulations")
					if f"Beta_list_case-0.csv" in os.listdir():
						#print(great)
						#print(os.getcwd())
						
						A_file = open("Beta_list_case-0.csv","r").readlines() 
						self.sim_ = len(A_file)
						print(self.sim_)
						A = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in A_file])	
						A_gen_file = open("Generator_list_case-0.csv","r").readlines() 
						A_gen = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in A_gen_file])	
						
						
						X_file = open("NormalizedBeta_list_case-0.csv","r").readlines() 
						X_dash = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in A_file])
						y_dash = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in A_file])
						self.beta_ = np.asarray(A)
						#print(self.beta_)
						self.generators = np.asarray(A_gen)
						
						self.X = np.asarray(X_dash)
						self.y = np.asarray(y_dash)
						os.chdir("../..")
					else:
						os.chdir("../..")
						task = "incomplete"
				if task != "complete":
					self.beta_ = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
					self.beta_.extend(list(np.eye(self.n)))
					#self.beta_.extend(list(-1*np.eye(self.n)))
					#self.beta_.extend(list(bbdesign(self.n)[0:int(0.3*self.sim_)]))
					#self.beta_.extend(list(np.zeros(self.n)))
					#self.beta_.extend(list(np.zeros(self.n)))
					#self.beta_.extend(list(np.zeros(self.n)))
					self.sim_ = len(self.beta_)
					#self.beta_ = np.asarray(A)
					#print(len(self.beta_))
					#print(self.beta_)
					self.generators = np.asarray(self.beta_)
					self.X = np.asarray(self.beta_)
					self.y = np.asarray(self.beta_)
			elif self.design == "Monte-Carlo-All":
				#mean = np.zeros(self.n)
				#cov = 0.33*np.eye(self.n)
				self.beta_ = list((2/3)*np.random.random_sample((self.sim_,self.n)) - (1/3))			
				
			elif self.design == "SAMAP":
				"""
				Getting the SAMAP amples
				"""
				os.mkdir("SAMAP_"+str(self.SAMAP_INDEX))
				os.chdir("SAMAP_"+str(self.SAMAP_INDEX))
				input_rxn_dict = {}
				samples_dict = {}
				for i in self.reactions_selected:
					os.mkdir(f"{self.unsrt[i].rIndex}")
					os.chdir(f"{self.unsrt[i].rIndex}")
					input_dict = self.unsrt[i].input_dict
					input_dict["samples"] = self.sim_
					input_dict["samples_skipped"] = int(0.1*self.sim_)
					input_dict["Random_seed"] = 1
					input_dict["sampling_method"] = "SOBOL"
					input_dict["sampling_distribution"] = "NORMAL"
					input_dict["equidistant_T"] = 100
					input_dict["T_begin"] = data = self.unsrt[i].temperatures[0]
					input_dict["T_end"] = self.unsrt[i].temperatures[-1]
					
					L = self.unsrt[i].cholskyDeCorrelateMat
					#input_dict["covariance_matrix"] = str(L.dot(L.T)).strip("[]").replace("[","").replace("]","")
					input_rxn_dict[i] = input_dict
					"""
					Run: SAMAP code
					"""
					instring = make_input_file.create_SAMAP_input(input_dict)
					#string_dict[i] = instring
					file_print = open("samap_data_"+str(self.unsrt[i].rIndex)+".txt","w").write(instring)
					run_string = f"""#!/bin/bash
		{self.samap_executable} samap_data_{self.unsrt[i].rIndex}.txt &> out"""
					file_print_run = open("run","w").write(run_string)						
					subprocess.call(["chmod","+x",'run'])
					subprocess.call(["./run"])
					file_name = "samap_data_"+str(self.unsrt[i].rIndex)+".txt_Arrhpar.txt"
					sampling_file = open(file_name,"r").readlines()
					self.Nagy_arrhenius_samples = [] 
					for sample in sampling_file[1:]:
						self.Nagy_arrhenius_samples.append(np.asfarray(sample.strip("''").strip("\n").split()[:3],float))
					
					#print(self.Nagy_arrhenius_samples)
					Y = self.Nagy_arrhenius_samples
					
					
					data = self.unsrt[i].data
					
					P = self.unsrt[i].getMean()	
					#print(P)
					cov = self.unsrt[i].getCov()
					T_ = np.linspace(300,2500,100)
					Theta = np.array([T_/T_,np.log(T_),-1/T_])
					"""
					fig = plt.figure()
					plt.plot(1/T_,self.unsrt[i].getKappaMax(T_),'k-',linewidth=0.75)
					plt.plot(1/T_,self.unsrt[i].getKappaMin(T_),'k-',linewidth=0.75)
					plt.plot(1/T_,self.unsrt[i].getNominal(T_),"b-",linewidth=0.75)
					"""
					temp = []
					for j in Y:
						#print(j)
						#print(kappa_max,np.asarray(theta.T.dot(j)).flatten())
						#print(np.asarray(theta.T.dot(j)).flatten()/kappa_max)
						#Pint = P + np.asarray(np.dot(cov,np.asarray(j))).flatten();
						Pint = np.asarray(j)
						ZETA = np.asarray(np.linalg.inv(cov).dot(Pint-P)).flatten()
						temp.append(ZETA)
							#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
						#X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#plt.plot(1/T_,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
					samples_dict[i] = temp
					#plt.show()
				
					os.chdir("..")
					
				A = []
				A_gen = []
				X = []
				y = []
				for j in range(self.sim_):
					a_row = []
					for i in samples_dict:
						a_row.extend(list(samples_dict[i][j]))
					
					A.append(np.asarray(a_row))
					
				self.beta_ = np.asarray(A)
				#print(self.beta_)
				self.generators = []
				self.X = []
				self.y = []
				
				self.SAMAP_INDEX+=1
				os.chdir("..")
			
			elif self.design == "B1":
				self.beta_=[]
				temp_a = {}
				temp_b = {}
				temp_c = {}
				gen = {}
				gen_a = {}
				gen_b = {}
				gen_c = {}
				
				xdata = {}
				ydata = {}				
				self.generator_list = []
				"""
				Min- simulation = 300
				"""
				#if self.sim_ < 300:
				#	self.sim_ = 300
				#print(self.sim_)
				"""
				Dividing the generators into three bins
				------------------------------------
				n_a = 0.33*self.sim_
				n_b = 0.33*self.sim_
				n_c = 0.34*self.sim_
				
				"""
				n_a = int(0.33*self.sim_)
				n_b = int(0.33*self.sim_)
				n_c = self.sim_-n_a-n_b
				
				
				#n_b = int(0.5*self.sim_)
				#n_c = self.sim_-n_b
				
				
				#for i in np.eye(self.n):
				#	self.generator_list.append(i)
				#	self.generator_list.append(-1*i)
								
				"""
				Mixing three types of the extreme curves
				--------------------
				zeta-A: 1
				zeta-B: 1
				zeta-C: 1
			
				for zeta-A: one normal random parameter is required
					a matrix on n_a x n will be created
					
				for zeta-B and zeta-c: two random parameters is required
					we can arrange the two normal_random_parameter 				
					a matrix of n_b+n_c x 2*n will be created
					
				The number of reactions = n_R
				"""
				n_rxn = len(self.ind)
				self.generator_list_A = []
				self.generator_list_B = []
				self.generator_list_C = []
				
				"""
				Factorial Design
				----------------
				box-benkhen design
				bb = []
				for i in bbdesign(n_rxn)[0:int(0.3*self.sim_)]:
					bb.append(np.asarray(i))
				"""
				"""
						
						#For zeta-A
						
				bb_a = []
				for i in bbdesign(n_rxn):
					bb_a.append(np.asarray(i))
				self.generator_list_A.extend(bb_a)
						
						#For zeta-B
						
				bb_b = []
				for i in bbdesign(2*n_rxn):
					bb_b.append(np.asarray(i))							
				self.generator_list_B.extend(bb_b)
						
						#For zeta-C
						
				bb_c = []
				for i in bbdesign(2*n_rxn):
					bb_c.append(np.asarray(i))		
						
				self.generator_list_C.extend(bb_c)		
						
				"""
				"""
				Monte-Carlo
				
				"""
				"""		
						#For zeta-A
						
				
				self.generator_list_A = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				
				
						
						#For zeta-A
						
				self.generator_list_B = list(2* np.random.random_sample((self.sim_,self.n)) - 1)		
				
						
						#For zeta-A
						
				self.generator_list_C = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				"""
				
				"""
				Uniform
				
				"""
					
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				
				"""
				Latin-Hypercube
				"""
			
						
						#For zeta-A
						
				lhs_a = []
				for i in lhs(n_rxn, samples=n_a):
					lhs_a.append(np.asarray(i))	
				self.generator_list_A.extend(lhs_a)
				
						
						#For zeta-A
						
				lhs_b = []
				for i in lhs(2*n_rxn, samples=n_b):
					lhs_b.append(np.asarray(i))		
				self.generator_list_B.extend(lhs_b)
						
						#For zeta-A
						
				lhs_c = []
				for i in lhs(2*n_rxn, samples=n_c):
					lhs_c.append(np.asarray(i))
				self.generator_list_C.extend(lhs_c)
				
				"""
				Monte-Carlo-Mix
				"""
				
				
				
				
				
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)		
				#self.generator_list.extend(bb)
				#self.generator_list.extend(list(np.eye(self.n)))
				#self.generator_list.extend(list(-1*np.eye(self.n)))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#print(self.generator_list)
				#self.sim_ = len(self.generator_list)
				
				#trainers = np.asarray(trainers)
				#print(f"Generator_list {np.shape(self.generator_list)}\n")
				df = pd.DataFrame(np.asarray(self.generator_list))
				count_start_a = 0
				count_start_b = 0
				count_start_c = 0			
				for i in self.reactions_selected:
					
					len_of_Arrhenius_A = 1
					count_end_a = count_start_a+len_of_Arrhenius_A
					
					
					len_of_Arrhenius_B = 2
					count_end_b = count_start_b+len_of_Arrhenius_B
					
					
					len_of_Arrhenius_C = 2
					count_end_c = count_start_c+len_of_Arrhenius_C
					
					#print(count_start)
					#print(count_end)
					self.rxn_generators_a2 = []
					self.rxn_generators_b2 = []
					self.rxn_generators_c2 = []
					
					"""
					getting zetas from generators
					"""
					for j in range(len(self.generator_list_A)):
						self.rxn_generators_a2.append(self.generator_list_A[j][count_start_a:count_end_a])
					#print(self.rxn_generators_a2)
					for j in range(len(self.generator_list_B)):
						self.rxn_generators_b2.append(self.generator_list_B[j][count_start_b:count_end_b])
					
					for j in range(len(self.generator_list_C)):
						self.rxn_generators_c2.append(self.generator_list_C[j][count_start_c:count_end_c])
					
					self.sim_ = len(self.rxn_generators_a2)+len(self.rxn_generators_b2)+len(self.rxn_generators_c2)
					#print(self.rxn_generators)
					#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
										
					g_a = self.rxn_generators_a2
					g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					#g_b = self.rxn_generators[0:n_b]
					#g_c = self.rxn_generators[n_b:n_b+n_c]
					
					
					count_start_a = count_end_a
					count_start_b = count_end_b
					count_start_c = count_end_c
					
					"""
					Getting the uncertainty data
					"""
					#fig = plt.figure()
					data = self.unsrt[i].data
					T = np.linspace(300,2500,8)
					kappa_max = self.unsrt[i].getKappaMax(T)			
					kappa_min = self.unsrt[i].getKappaMin(T)
					Theta = np.array([T/T,np.log(T),-1/T])
					P = self.unsrt[i].getMean()	
					cov = self.unsrt[i].getCov()
					zeta_A = self.unsrt[i].zeta.x
					
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					"""
					Type-A zeta samples
					"""
					generators_a = g_a
					parallel_zetas_a = []
					
					time_start = time.time()
					
					
					for k in g_a:
						#Pint = P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten()
					#	plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
					#	plt.plot()
						#print(i)
						#print(type(i))
						#print(k)	
						#parallel_zetas_a.append(P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten())
						parallel_zetas_a.append(k[0]*zeta_A)
										
					#print(parallel_zetas_a)
					#plt.savefig(str(i)+"_a.pdf")
					#print(os.getcwd())
					#print("plot saved")
					data["generators_b"] = g_b
					callWorkForce = Worker(self.allowed_count)	
					generators_b,parallel_zetas_b = callWorkForce.do_unsrt_b(data,len(g_b))
					del callWorkForce
					
					data["generators_c"] = g_c
					callWorkForce = Worker(self.allowed_count)	
					generators_c,parallel_zetas_c = callWorkForce.do_unsrt_c(data,len(g_c))
					del callWorkForce
					
					#print(parallel_zetas_a)
					
					temp_a[i] = parallel_zetas_a
					gen_a[i] = generators_a
					temp_b[i] = parallel_zetas_b
					gen_b[i] = generators_b
					temp_c[i] = parallel_zetas_c
					gen_c[i] = generators_c
					
					time_end = time.time()
					print(f"\n\t\tTime taken = {time_end-time_start}\n")
					
					#print(parallel_zetas_b)
					Y_a = parallel_zetas_a
					Y_b = parallel_zetas_b
					Y_c = parallel_zetas_c
					T_ = np.linspace(300,2500,8)
					theta = np.array([T/T,np.log(T),-1/T])
					kappa_mean = theta.T.dot(P)
					Theta = np.array([T/T,np.log(T),-1/T])
					X = []
					Y = []
					
					
					for j in Y_a:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					
					for j in Y_b:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_b.pdf")
					
					
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					for j in Y_c:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_c.pdf")
					
					temp_xdata = {}
					for k in range(len(Y)):
						temp_xdata[str(Y[k])] = X[k]
						
					random.shuffle(Y)
					#print(Y_)
					X_ = []
					for k in range(len(Y)):
						X_.append(temp_xdata[str(Y[k])])
						
					#print(X_[9])
					xdata[i] = X_
					ydata[i] = Y
					#print(len(parallel_zetas))
					#print(len(self.rxn_generators))
				A = []
				A_gen = []
				X_dash = []
				y_dash = []
				
				for j in range(len(g_a)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_a:
						a_row.extend(list(temp_a[i][j]))
						a_row_gen.extend(list(gen_a[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				
				
				for j in range(len(g_b)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)]))
						y_temp.extend(list(ydata[i][j+len(g_a)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				for j in range(n_b):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				
				for j in range(len(g_c)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)+len(g_b)]))
						y_temp.extend(list(ydata[i][j+len(g_a)+len(g_b)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				for j in range(n_c):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+n_b]))
						y_temp.extend(list(ydata[i][j+n_b]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
					
				"""
				
				#print(A)
				self.sim_ = len(g_a)+len(g_b)+len(g_c)
				#self.sim_ = n_b+n_c
				self.beta_ = np.asarray(A)
				#print(self.beta_)
				self.generators = np.asarray(A_gen)
				self.X = np.asarray(X_dash)
				self.y = np.asarray(y_dash)
			
			
			elif self.design == "B1+C1":
				self.beta_=[]
				temp_a = {}
				temp_b = {}
				temp_c = {}
				gen = {}
				gen_a = {}
				gen_b = {}
				gen_c = {}
				
				xdata = {}
				ydata = {}				
				self.generator_list = []
				"""
				Min- simulation = 300
				"""
				#if self.sim_ < 300:
				#	self.sim_ = 300
				#print(self.sim_)
				"""
				Dividing the generators into three bins
				------------------------------------
				n_a = 0.33*self.sim_
				n_b = 0.33*self.sim_
				n_c = 0.34*self.sim_
				
				"""
				n_a = int(0.33*self.sim_)
				#n_b = int(0.5*self.sim_)
				#n_c = self.sim_-n_a-n_b
				
				
				n_b = int(0.5*self.sim_)
				n_c = self.sim_-n_b
				
				
				#for i in np.eye(self.n):
				#	self.generator_list.append(i)
				#	self.generator_list.append(-1*i)
								
				"""
				Mixing three types of the extreme curves
				--------------------
				zeta-A: 1
				zeta-B: 1
				zeta-C: 1
			
				for zeta-A: one normal random parameter is required
					a matrix on n_a x n will be created
					
				for zeta-B and zeta-c: two random parameters is required
					we can arrange the two normal_random_parameter 				
					a matrix of n_b+n_c x 2*n will be created
					
				The number of reactions = n_R
				"""
				n_rxn = len(self.ind)
				self.generator_list_A = []
				self.generator_list_B = []
				self.generator_list_C = []
				
				"""
				Factorial Design
				----------------
				box-benkhen design
				bb = []
				for i in bbdesign(n_rxn)[0:int(0.3*self.sim_)]:
					bb.append(np.asarray(i))
				"""
						#For zeta-A
						
				bb_a = []
				for i in bbdesign(n_rxn)[0:int(0.4*n_a)]:
					bb_a.append(np.asarray(i))
				self.generator_list_A.extend(bb_a)
						
						#For zeta-B
						
				bb_b = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_b)]:
					bb_b.append(np.asarray(i))							
				self.generator_list_B.extend(bb_b)
						
						#For zeta-C
						
				bb_c = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_c)]:
					bb_c.append(np.asarray(i))		
						
				self.generator_list_C.extend(bb_c)		
						
				
				"""
				Monte-Carlo
				
				"""
						
						#For zeta-A
				
				self.generator_list_A.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(0.3*n_a)])
				self.generator_list_A.extend(list(np.eye(n_rxn)))
				self.generator_list_A.extend(list(-1*np.eye(n_rxn)))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				
						#For zeta-A
						
				self.generator_list_B.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_b)])	
				self.generator_list_B.extend(list(np.eye(2*n_rxn)))
				self.generator_list_B.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
						
						#For zeta-A		
				self.generator_list_C.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_c)])
				self.generator_list_C.extend(list(np.eye(2*n_rxn)))
				self.generator_list_C.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				
				"""
				Uniform
				
				"""
					
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				
				"""
				Latin-Hypercube
				"""
			
						
						#For zeta-A
						
				lhs_a = []
				for i in lhs(n_rxn, samples=int(0.3*n_a)):
					lhs_a.append(np.asarray(i))	
				self.generator_list_A.extend(lhs_a)
				
						
						#For zeta-A
						
				lhs_b = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_b)):
					lhs_b.append(np.asarray(i))		
				self.generator_list_B.extend(lhs_b)
						
						#For zeta-A
						
				lhs_c = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_c)):
					lhs_c.append(np.asarray(i))
				self.generator_list_C.extend(lhs_c)
				
				"""
				Monte-Carlo-Mix
				"""
				
				
				
				
				
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)		
				#self.generator_list.extend(bb)
				#self.generator_list.extend(list(np.eye(self.n)))
				#self.generator_list.extend(list(-1*np.eye(self.n)))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#print(self.generator_list)
				#self.sim_ = len(self.generator_list)
				
				#trainers = np.asarray(trainers)
				#print(f"Generator_list {np.shape(self.generator_list)}\n")
				df = pd.DataFrame(np.asarray(self.generator_list))
				count_start_a = 0
				count_start_b = 0
				count_start_c = 0			
				for i in self.reactions_selected:
					
					len_of_Arrhenius_A = 1
					count_end_a = count_start_a+len_of_Arrhenius_A
					
					
					len_of_Arrhenius_B = 2
					count_end_b = count_start_b+len_of_Arrhenius_B
					
					
					len_of_Arrhenius_C = 2
					count_end_c = count_start_c+len_of_Arrhenius_C
					
					#print(count_start)
					#print(count_end)
					self.rxn_generators_a2 = []
					self.rxn_generators_b2 = []
					self.rxn_generators_c2 = []
					
					"""
					getting zetas from generators
					"""
					for j in range(len(self.generator_list_A)):
						self.rxn_generators_a2.append(self.generator_list_A[j][count_start_a:count_end_a])
					#print(self.rxn_generators_a2)
					for j in range(len(self.generator_list_B)):
						self.rxn_generators_b2.append(self.generator_list_B[j][count_start_b:count_end_b])
					
					for j in range(len(self.generator_list_C)):
						self.rxn_generators_c2.append(self.generator_list_C[j][count_start_c:count_end_c])
					
					#self.sim_ = len(self.rxn_generators_a2)+len(self.rxn_generators_b2)+len(self.rxn_generators_c2)
					self.sim_ = len(self.rxn_generators_b2)+len(self.rxn_generators_c2)
					#print(self.rxn_generators)
					#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
										
					g_a = self.rxn_generators_a2
					g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					#g_b = self.rxn_generators[0:n_b]
					#g_c = self.rxn_generators[n_b:n_b+n_c]
					
					
					count_start_a = count_end_a
					count_start_b = count_end_b
					count_start_c = count_end_c
					
					"""
					Getting the uncertainty data
					"""
					#fig = plt.figure()
					data = self.unsrt[i].data
					T = np.linspace(300,2500,8)
					kappa_max = self.unsrt[i].getKappaMax(T)			
					kappa_min = self.unsrt[i].getKappaMin(T)
					Theta = np.array([T/T,np.log(T),-1/T])
					P = self.unsrt[i].getMean()	
					cov = self.unsrt[i].getCov()
					zeta_A = self.unsrt[i].zeta.x
					
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					"""
					Type-A zeta samples
					"""
					generators_a = g_a
					parallel_zetas_a = []
					
					time_start = time.time()
					
					
					for k in g_a:
						#Pint = P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten()
					#	plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
					#	plt.plot()
						#print(i)
						#print(type(i))
						#print(k)	
						#parallel_zetas_a.append(P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten())
						parallel_zetas_a.append(k[0]*zeta_A)
										
					#print(parallel_zetas_a)
					#plt.savefig(str(i)+"_a.pdf")
					#print(os.getcwd())
					#print("plot saved")
					g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					
					generators_b = []
					parallel_zetas_b = []
					generators_c = []
					parallel_zetas_c = []
					start = 0
					for w in range(10):
						stop = int(((w+1)/10)*self.sim_)
						g_b = self.rxn_generators_b2[start:stop]
						g_c = self.rxn_generators_c2[start:stop]
						data["generators_b"] = g_b
						callWorkForce = Worker(self.allowed_count)	
						generators_b_,parallel_zetas_b_ = callWorkForce.do_unsrt_b(data,len(g_b))
						del callWorkForce
						
						
						removed = data.pop("generators_b")
						data["generators_c"] = g_c
						callWorkForce = Worker(self.allowed_count)	
						generators_c_,parallel_zetas_c_ = callWorkForce.do_unsrt_c(data,len(g_c))
						del callWorkForce
						
						generators_b.extend(generators_b_)
						parallel_zetas_b.extend(parallel_zetas_b_)
						generators_c.extend(generators_c_)
						parallel_zetas_c.extend(parallel_zetas_c_)
						removed = data.pop("generators_c")
						start = stop
					
					#print(self.sim_)
					#print(len(parallel_zetas_b)+len(parallel_zetas_c))
					g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					
					#print(parallel_zetas_a)
					
					temp_a[i] = parallel_zetas_a
					gen_a[i] = generators_a
					temp_b[i] = parallel_zetas_b
					gen_b[i] = generators_b
					temp_c[i] = parallel_zetas_c
					gen_c[i] = generators_c
					
					time_end = time.time()
					print(f"\n\t\tTime taken = {time_end-time_start}\n")
					
					#print(parallel_zetas_b)
					Y_a = parallel_zetas_a
					Y_b = parallel_zetas_b
					Y_c = parallel_zetas_c
					T_ = np.linspace(300,2500,8)
					theta = np.array([T/T,np.log(T),-1/T])
					kappa_mean = theta.T.dot(P)
					Theta = np.array([T/T,np.log(T),-1/T])
					X = []
					Y = []
					
					"""
					for j in Y_a:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					"""
					for j in Y_b:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_b.pdf")
					
					
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					for j in Y_c:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_c.pdf")
					
					temp_xdata = {}
					for k in range(len(Y)):
						temp_xdata[str(Y[k])] = X[k]
						
					random.shuffle(Y)
					#print(Y_)
					X_ = []
					for k in range(len(Y)):
						X_.append(temp_xdata[str(Y[k])])
						
					#print(X_[9])
					xdata[i] = X_
					ydata[i] = Y
					#print(len(parallel_zetas))
					#print(len(self.rxn_generators))
				A = []
				A_gen = []
				X_dash = []
				y_dash = []
				"""
				for j in range(len(g_a)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_a:
						a_row.extend(list(temp_a[i][j]))
						a_row_gen.extend(list(gen_a[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				"""
				for j in range(len(g_b)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)]))
						y_temp.extend(list(ydata[i][j+len(g_a)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				
				
				for j in range(len(g_b)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				
				for j in range(len(g_c)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)+len(g_b)]))
						y_temp.extend(list(ydata[i][j+len(g_a)+len(g_b)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				for j in range(len(g_c)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_b)]))
						y_temp.extend(list(ydata[i][j+len(g_b)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
					
				
				
				#print(A)
				#self.sim_ = len(g_a)+len(g_b)+len(g_c)
				#self.sim_ = n_b+n_c
				self.beta_ = np.asarray(A)
				#print(self.beta_)
				self.generators = np.asarray(A_gen)
				self.X = np.asarray(X_dash)
				self.y = np.asarray(y_dash)
			
			
			elif self.design == "A1+B1":
				self.beta_=[]
				temp_a = {}
				temp_b = {}
				temp_c = {}
				gen = {}
				gen_a = {}
				gen_b = {}
				gen_c = {}
				
				xdata = {}
				ydata = {}				
				self.generator_list = []
				"""
				Min- simulation = 300
				"""
				#if self.sim_ < 300:
				#	self.sim_ = 300
				#print(self.sim_)
				"""
				Dividing the generators into three bins
				------------------------------------
				n_a = 0.33*self.sim_
				n_b = 0.33*self.sim_
				n_c = 0.34*self.sim_
				
				"""
				n_a = int(0.1*self.sim_)
				#n_b = int(0.5*self.sim_)
				#n_c = self.sim_-n_a-n_b
				
				
				n_b = int(0.9*self.sim_)
				n_c = self.sim_-n_b
				
				
				#for i in np.eye(self.n):
				#	self.generator_list.append(i)
				#	self.generator_list.append(-1*i)
								
				"""
				Mixing three types of the extreme curves
				--------------------
				zeta-A: 1
				zeta-B: 1
				zeta-C: 1
			
				for zeta-A: one normal random parameter is required
					a matrix on n_a x n will be created
					
				for zeta-B and zeta-c: two random parameters is required
					we can arrange the two normal_random_parameter 				
					a matrix of n_b+n_c x 2*n will be created
					
				The number of reactions = n_R
				"""
				n_rxn = len(self.ind)
				self.generator_list_A = []
				self.generator_list_B = []
				self.generator_list_C = []
				
				"""
				Factorial Design
				----------------
				box-benkhen design
				bb = []
				for i in bbdesign(n_rxn)[0:int(0.3*self.sim_)]:
					bb.append(np.asarray(i))
				"""
						#For zeta-A
						
				bb_a = []
				for i in bbdesign(n_rxn)[0:int(0.4*n_a)]:
					bb_a.append(np.asarray(i))
				self.generator_list_A.extend(bb_a)
						
						#For zeta-B
						
				bb_b = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_b)]:
					bb_b.append(np.asarray(i))							
				self.generator_list_B.extend(bb_b)
						
						#For zeta-C
						
				bb_c = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_c)]:
					bb_c.append(np.asarray(i))		
						
				self.generator_list_C.extend(bb_c)		
						
				
				"""
				Monte-Carlo
				
				"""
						
						#For zeta-A
				
				self.generator_list_A.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(0.3*n_a)])
				self.generator_list_A.extend(list(np.eye(n_rxn)))
				self.generator_list_A.extend(list(-1*np.eye(n_rxn)))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				
						#For zeta-A
						
				self.generator_list_B.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_b)])	
				self.generator_list_B.extend(list(np.eye(2*n_rxn)))
				self.generator_list_B.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
						
						#For zeta-A		
				self.generator_list_C.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_c)])
				self.generator_list_C.extend(list(np.eye(2*n_rxn)))
				self.generator_list_C.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				
				"""
				Uniform
				
				"""
					
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				
				"""
				Latin-Hypercube
				"""
			
						
						#For zeta-A
						
				lhs_a = []
				for i in lhs(n_rxn, samples=int(0.3*n_a)):
					lhs_a.append(np.asarray(i))	
				self.generator_list_A.extend(lhs_a)
				
						
						#For zeta-A
						
				lhs_b = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_b)):
					lhs_b.append(np.asarray(i))		
				self.generator_list_B.extend(lhs_b)
						
						#For zeta-A
						
				lhs_c = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_c)):
					lhs_c.append(np.asarray(i))
				self.generator_list_C.extend(lhs_c)
				
				"""
				Monte-Carlo-Mix
				"""
				
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)		
				#self.generator_list.extend(bb)
				#self.generator_list.extend(list(np.eye(self.n)))
				#self.generator_list.extend(list(-1*np.eye(self.n)))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#print(self.generator_list)
				#self.sim_ = len(self.generator_list)
				
				#trainers = np.asarray(trainers)
				#print(f"Generator_list {np.shape(self.generator_list)}\n")
				df = pd.DataFrame(np.asarray(self.generator_list))
				count_start_a = 0
				count_start_b = 0
				count_start_c = 0			
				for i in self.reactions_selected:
					
					len_of_Arrhenius_A = 1
					count_end_a = count_start_a+len_of_Arrhenius_A
					
					
					len_of_Arrhenius_B = 2
					count_end_b = count_start_b+len_of_Arrhenius_B
					
					
					len_of_Arrhenius_C = 2
					count_end_c = count_start_c+len_of_Arrhenius_C
					
					#print(count_start)
					#print(count_end)
					self.rxn_generators_a2 = []
					self.rxn_generators_b2 = []
					self.rxn_generators_c2 = []
					
					"""
					getting zetas from generators
					"""
					for j in range(len(self.generator_list_A)):
						self.rxn_generators_a2.append(self.generator_list_A[j][count_start_a:count_end_a])
					#print(self.rxn_generators_a2)
					for j in range(len(self.generator_list_B)):
						self.rxn_generators_b2.append(self.generator_list_B[j][count_start_b:count_end_b])
					
					for j in range(len(self.generator_list_C)):
						self.rxn_generators_c2.append(self.generator_list_C[j][count_start_c:count_end_c])
					
					#self.sim_ = len(self.rxn_generators_a2)+len(self.rxn_generators_b2)+len(self.rxn_generators_c2)
					self.sim_ = len(self.rxn_generators_a2)+len(self.rxn_generators_b2)
					#print(self.rxn_generators)
					#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
										
					g_a = self.rxn_generators_a2
					g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					#g_b = self.rxn_generators[0:n_b]
					#g_c = self.rxn_generators[n_b:n_b+n_c]
					
					
					count_start_a = count_end_a
					count_start_b = count_end_b
					count_start_c = count_end_c
					
					"""
					Getting the uncertainty data
					"""
					#fig = plt.figure()
					data = self.unsrt[i].data
					T = np.linspace(300,2500,8)
					kappa_max = self.unsrt[i].getKappaMax(T)			
					kappa_min = self.unsrt[i].getKappaMin(T)
					Theta = np.array([T/T,np.log(T),-1/T])
					P = self.unsrt[i].getMean()	
					cov = self.unsrt[i].getCov()
					zeta_A = self.unsrt[i].zeta.x
					
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					"""
					Type-A zeta samples
					"""
					generators_a = g_a
					parallel_zetas_a = []
					
					time_start = time.time()
					
					
					for k in g_a:
						parallel_zetas_a.append(k[0]*zeta_A)
										
					generators_b = []
					parallel_zetas_b = []
					generators_c = []
					parallel_zetas_c = []
					start = 0
					for w in range(10):
						stop = int(((w+1)/10)*self.sim_)
						g_b = self.rxn_generators_b2[start:stop]
						data["generators_b"] = g_b
						callWorkForce = Worker(self.allowed_count)	
						generators_b_,parallel_zetas_b_ = callWorkForce.do_unsrt_b(data,len(g_b))
						del callWorkForce
						removed = data.pop("generators_b")
						generators_b.extend(generators_b_)
						parallel_zetas_b.extend(parallel_zetas_b_)
						start = stop
					
					#print(self.sim_)
					#print(len(parallel_zetas_b)+len(parallel_zetas_c))
					g_b = self.rxn_generators_b2
					#g_c = self.rxn_generators_c2
					#print(parallel_zetas_a)
					
					temp_a[i] = parallel_zetas_a
					gen_a[i] = generators_a
					temp_b[i] = parallel_zetas_b
					gen_b[i] = generators_b
					#temp_c[i] = parallel_zetas_c
					#gen_c[i] = generators_c
					
					time_end = time.time()
					print(f"\n\t\tTime taken = {time_end-time_start}\n")
					
					#print(parallel_zetas_b)
					Y_a = parallel_zetas_a
					Y_b = parallel_zetas_b
					#Y_c = parallel_zetas_c
					T_ = np.linspace(300,2500,8)
					theta = np.array([T/T,np.log(T),-1/T])
					kappa_mean = theta.T.dot(P)
					Theta = np.array([T/T,np.log(T),-1/T])
					X = []
					Y = []
					
					
					for j in Y_a:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					
					for j in Y_b:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_b.pdf")
					
					
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					"""
					for j in Y_c:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_c.pdf")
					"""
					#xdata[i] = X
					#ydata[i] = Y

					temp_xdata = {}
					for k in range(len(Y)):
						temp_xdata[str(Y[k])] = X[k]
						
					random.shuffle(Y)
					#print(Y_)
					X_ = []
					for k in range(len(Y)):
						X_.append(temp_xdata[str(Y[k])])
						
					#print(X_[9])
					xdata[i] = X_
					ydata[i] = Y
					#print(len(parallel_zetas))
					#print(len(self.rxn_generators))
				A = []
				A_gen = []
				X_dash = []
				y_dash = []
				
				for j in range(len(g_a)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_a:
						a_row.extend(list(ydata[i][j]))
						a_row_gen.extend(list(gen_a[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				
				
				for j in range(len(g_b)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(ydata[i][j+len(g_a)]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)]))
						y_temp.extend(list(ydata[i][j+len(g_a)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				
				
				"""
				for j in range(len(g_b)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				"""
				
				for j in range(len(g_c)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)+len(g_b)]))
						y_temp.extend(list(ydata[i][j+len(g_a)+len(g_b)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				"""
				for j in range(len(g_c)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_b)]))
						y_temp.extend(list(ydata[i][j+len(g_b)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
					
				"""
				
				#print(A)
				#self.sim_ = len(g_a)+len(g_b)+len(g_c)
				#self.sim_ = n_b+n_c
				self.beta_ = np.asarray(A)
				#print(self.beta_)
				self.generators = np.asarray(A_gen)
				self.X = np.asarray(X_dash)
				self.y = np.asarray(y_dash)
			
			elif self.design == "B1-only":
				self.beta_=[]
				temp_a = {}
				temp_b = {}
				temp_c = {}
				gen = {}
				gen_a = {}
				gen_b = {}
				gen_c = {}
				
				xdata = {}
				ydata = {}				
				self.generator_list = []
				"""
				Min- simulation = 300
				"""
				#if self.sim_ < 300:
				#	self.sim_ = 300
				#print(self.sim_)
				"""
				Dividing the generators into three bins
				------------------------------------
				n_a = 0.33*self.sim_
				n_b = 0.33*self.sim_
				n_c = 0.34*self.sim_
				
				"""
				n_a = int(0.5*self.sim_)
				#n_b = int(0.5*self.sim_)
				#n_c = self.sim_-n_a-n_b
				
				
				n_b = int(self.sim_)
				n_c = int(self.sim_-0.33*n_b)
				
				
				#for i in np.eye(self.n):
				#	self.generator_list.append(i)
				#	self.generator_list.append(-1*i)
								
				"""
				Mixing three types of the extreme curves
				--------------------
				zeta-A: 1
				zeta-B: 1
				zeta-C: 1
			
				for zeta-A: one normal random parameter is required
					a matrix on n_a x n will be created
					
				for zeta-B and zeta-c: two random parameters is required
					we can arrange the two normal_random_parameter 				
					a matrix of n_b+n_c x 2*n will be created
					
				The number of reactions = n_R
				"""
				n_rxn = len(self.ind)
				self.generator_list_A = []
				self.generator_list_B = []
				self.generator_list_C = []
				
				"""
				Factorial Design
				----------------
				box-benkhen design
				bb = []
				for i in bbdesign(n_rxn)[0:int(0.3*self.sim_)]:
					bb.append(np.asarray(i))
				"""
						#For zeta-A
						
				bb_a = []
				for i in bbdesign(n_rxn)[0:int(0.4*n_a)]:
					bb_a.append(np.asarray(i))
				self.generator_list_A.extend(bb_a)
						
						#For zeta-B
						
				bb_b = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_b)]:
					bb_b.append(np.asarray(i))							
				self.generator_list_B.extend(bb_b)
						
						#For zeta-C
						
				bb_c = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_c)]:
					bb_c.append(np.asarray(i))		
						
				self.generator_list_C.extend(bb_c)		
						
				
				"""
				Monte-Carlo
				
				"""
						
						#For zeta-A
				
				self.generator_list_A.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(0.3*n_a)])
				self.generator_list_A.extend(list(np.eye(n_rxn)))
				self.generator_list_A.extend(list(-1*np.eye(n_rxn)))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				
						#For zeta-A
						
				self.generator_list_B.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_b)])	
				self.generator_list_B.extend(list(np.eye(2*n_rxn)))
				self.generator_list_B.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
						
						#For zeta-A		
				self.generator_list_C.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_c)])
				self.generator_list_C.extend(list(np.eye(2*n_rxn)))
				self.generator_list_C.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				
				"""
				Uniform
				
				"""
					
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				
				"""
				Latin-Hypercube
				"""
			
						
						#For zeta-A
						
				lhs_a = []
				for i in lhs(n_rxn, samples=int(0.3*n_a)):
					lhs_a.append(np.asarray(i))	
				self.generator_list_A.extend(lhs_a)
				
						
						#For zeta-A
						
				lhs_b = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_b)):
					lhs_b.append(np.asarray(i))		
				self.generator_list_B.extend(lhs_b)
						
						#For zeta-A
						
				lhs_c = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_c)):
					lhs_c.append(np.asarray(i))
				self.generator_list_C.extend(lhs_c)
				
				"""
				Monte-Carlo-Mix
				"""
				
				
				df = pd.DataFrame(np.asarray(self.generator_list))
				count_start_a = 0
				count_start_b = 0
				count_start_c = 0			
				for i in self.reactions_selected:
					
					len_of_Arrhenius_A = 1
					count_end_a = count_start_a+len_of_Arrhenius_A
					
					
					len_of_Arrhenius_B = 2
					count_end_b = count_start_b+len_of_Arrhenius_B
					
					
					len_of_Arrhenius_C = 2
					count_end_c = count_start_c+len_of_Arrhenius_C
					
					#print(count_start)
					#print(count_end)
					self.rxn_generators_a2 = []
					self.rxn_generators_b2 = []
					self.rxn_generators_c2 = []
					
					"""
					getting zetas from generators
					"""
					for j in range(len(self.generator_list_A)):
						self.rxn_generators_a2.append(self.generator_list_A[j][count_start_a:count_end_a])
					#print(self.rxn_generators_a2)
					for j in range(len(self.generator_list_B)):
						self.rxn_generators_b2.append(self.generator_list_B[j][count_start_b:count_end_b])
					
					for j in range(len(self.generator_list_C)):
						self.rxn_generators_c2.append(self.generator_list_C[j][count_start_c:count_end_c])
					
					#self.sim_ = len(self.rxn_generators_a2)+len(self.rxn_generators_b2)+len(self.rxn_generators_c2)
					self.sim_ = len(self.rxn_generators_b2)
					#print(self.rxn_generators)
					#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
										
					g_a = self.rxn_generators_a2
					g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					#g_b = self.rxn_generators[0:n_b]
					#g_c = self.rxn_generators[n_b:n_b+n_c]
					
					
					count_start_a = count_end_a
					count_start_b = count_end_b
					count_start_c = count_end_c
					
					"""
					Getting the uncertainty data
					"""
					#fig = plt.figure()
					data = self.unsrt[i].data
					T = np.linspace(300,2500,8)
					kappa_max = self.unsrt[i].getKappaMax(T)			
					kappa_min = self.unsrt[i].getKappaMin(T)
					Theta = np.array([T/T,np.log(T),-1/T])
					P = self.unsrt[i].getMean()	
					cov = self.unsrt[i].getCov()
					zeta_A = self.unsrt[i].zeta.x
					
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					"""
					Type-A zeta samples
					"""
					generators_a = g_a
					parallel_zetas_a = []
					
					time_start = time.time()
					
					
					for k in g_a:
						
						parallel_zetas_a.append(k[0]*zeta_A)
										
					data["generators_b"] = g_b
					callWorkForce = Worker(self.allowed_count)	
					generators_b,parallel_zetas_b = callWorkForce.do_unsrt_b(data,len(g_b))
					del callWorkForce
					"""
					data["generators_c"] = g_c
					callWorkForce = Worker(self.allowed_count)	
					generators_c,parallel_zetas_c = callWorkForce.do_unsrt_c(data,len(g_c))
					del callWorkForce
					"""
					#print(parallel_zetas_a)
					
					temp_a[i] = parallel_zetas_a
					gen_a[i] = generators_a
					temp_b[i] = parallel_zetas_b
					gen_b[i] = generators_b
					#temp_c[i] = parallel_zetas_c
					#gen_c[i] = generators_c
					
					time_end = time.time()
					print(f"\n\t\tTime taken = {time_end-time_start}\n")
					
					#print(parallel_zetas_b)
					Y_a = parallel_zetas_a
					Y_b = parallel_zetas_b
					#Y_c = parallel_zetas_c
					T_ = np.linspace(300,2500,8)
					theta = np.array([T/T,np.log(T),-1/T])
					kappa_mean = theta.T.dot(P)
					Theta = np.array([T/T,np.log(T),-1/T])
					X = []
					Y = []
					
					for j in Y_b:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						
					
					temp_xdata = {}
					for k in range(len(Y)):
						temp_xdata[str(Y[k])] = X[k]
						
					random.shuffle(Y)
					#print(Y_)
					X_ = []
					for k in range(len(Y)):
						X_.append(temp_xdata[str(Y[k])])
						
					#print(X_[9])
					xdata[i] = X_
					ydata[i] = Y
					#print(len(parallel_zetas))
					#print(len(self.rxn_generators))
				A = []
				A_gen = []
				X_dash = []
				y_dash = []
								
				for j in range(len(g_b)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				
				
				
				self.beta_ = np.asarray(A)
				#print(self.beta_)
				self.generators = np.asarray(A_gen)
				self.X = np.asarray(X_dash)
				self.y = np.asarray(y_dash)
			
			elif self.design == "C1-only":
				self.beta_=[]
				temp_a = {}
				temp_b = {}
				temp_c = {}
				gen = {}
				gen_a = {}
				gen_b = {}
				gen_c = {}
				
				xdata = {}
				ydata = {}				
				self.generator_list = []
				"""
				Min- simulation = 300
				"""
				#if self.sim_ < 300:
				#	self.sim_ = 300
				#print(self.sim_)
				"""
				Dividing the generators into three bins
				------------------------------------
				n_a = 0.33*self.sim_
				n_b = 0.33*self.sim_
				n_c = 0.34*self.sim_
				
				"""
				n_a = int(0.5*self.sim_)
				#n_b = int(0.5*self.sim_)
				#n_c = self.sim_-n_a-n_b
				
				
				n_b = int(self.sim_)
				n_c = self.sim_
				
				
				#for i in np.eye(self.n):
				#	self.generator_list.append(i)
				#	self.generator_list.append(-1*i)
								
				"""
				Mixing three types of the extreme curves
				--------------------
				zeta-A: 1
				zeta-B: 1
				zeta-C: 1
			
				for zeta-A: one normal random parameter is required
					a matrix on n_a x n will be created
					
				for zeta-B and zeta-c: two random parameters is required
					we can arrange the two normal_random_parameter 				
					a matrix of n_b+n_c x 2*n will be created
					
				The number of reactions = n_R
				"""
				n_rxn = len(self.ind)
				self.generator_list_A = []
				self.generator_list_B = []
				self.generator_list_C = []
				
				"""
				Factorial Design
				----------------
				box-benkhen design
				bb = []
				for i in bbdesign(n_rxn)[0:int(0.3*self.sim_)]:
					bb.append(np.asarray(i))
				"""
						#For zeta-A
						
				bb_a = []
				for i in bbdesign(n_rxn)[0:int(0.4*n_a)]:
					bb_a.append(np.asarray(i))
				self.generator_list_A.extend(bb_a)
						
						#For zeta-B
						
				bb_b = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_b)]:
					bb_b.append(np.asarray(i))							
				self.generator_list_B.extend(bb_b)
						
						#For zeta-C
						
				bb_c = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_c)]:
					bb_c.append(np.asarray(i))		
						
				self.generator_list_C.extend(bb_c)		
						
				
				"""
				Monte-Carlo
				
				"""
						
						#For zeta-A
				
				self.generator_list_A.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(0.3*n_a)])
				self.generator_list_A.extend(list(np.eye(n_rxn)))
				self.generator_list_A.extend(list(-1*np.eye(n_rxn)))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				
						#For zeta-A
						
				self.generator_list_B.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_b)])	
				self.generator_list_B.extend(list(np.eye(2*n_rxn)))
				self.generator_list_B.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
						
						#For zeta-A		
				self.generator_list_C.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_c)])
				self.generator_list_C.extend(list(np.eye(2*n_rxn)))
				self.generator_list_C.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				
				"""
				Uniform
				
				"""
					
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				
				"""
				Latin-Hypercube
				"""
			
						
						#For zeta-A
						
				lhs_a = []
				for i in lhs(n_rxn, samples=int(0.3*n_a)):
					lhs_a.append(np.asarray(i))	
				self.generator_list_A.extend(lhs_a)
				
						
						#For zeta-A
						
				lhs_b = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_b)):
					lhs_b.append(np.asarray(i))		
				self.generator_list_B.extend(lhs_b)
						
						#For zeta-A
						
				lhs_c = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_c)):
					lhs_c.append(np.asarray(i))
				self.generator_list_C.extend(lhs_c)
				
				"""
				Monte-Carlo-Mix
				"""
				
				
				df = pd.DataFrame(np.asarray(self.generator_list))
				count_start_a = 0
				count_start_b = 0
				count_start_c = 0			
				for i in self.reactions_selected:
					
					len_of_Arrhenius_A = 1
					count_end_a = count_start_a+len_of_Arrhenius_A
					
					
					len_of_Arrhenius_B = 2
					count_end_b = count_start_b+len_of_Arrhenius_B
					
					
					len_of_Arrhenius_C = 2
					count_end_c = count_start_c+len_of_Arrhenius_C
					
					#print(count_start)
					#print(count_end)
					self.rxn_generators_a2 = []
					self.rxn_generators_b2 = []
					self.rxn_generators_c2 = []
					
					"""
					getting zetas from generators
					"""
					for j in range(len(self.generator_list_A)):
						self.rxn_generators_a2.append(self.generator_list_A[j][count_start_a:count_end_a])
					#print(self.rxn_generators_a2)
					for j in range(len(self.generator_list_B)):
						self.rxn_generators_b2.append(self.generator_list_B[j][count_start_b:count_end_b])
					
					for j in range(len(self.generator_list_C)):
						self.rxn_generators_c2.append(self.generator_list_C[j][count_start_c:count_end_c])
					
					#self.sim_ = len(self.rxn_generators_a2)+len(self.rxn_generators_b2)+len(self.rxn_generators_c2)
					self.sim_ = len(self.rxn_generators_c2)
					#print(self.rxn_generators)
					#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
										
					g_a = self.rxn_generators_a2
					g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					#g_b = self.rxn_generators[0:n_b]
					#g_c = self.rxn_generators[n_b:n_b+n_c]
					
					
					count_start_a = count_end_a
					count_start_b = count_end_b
					count_start_c = count_end_c
					
					"""
					Getting the uncertainty data
					"""
					#fig = plt.figure()
					data = self.unsrt[i].data
					T = np.linspace(300,2500,8)
					kappa_max = self.unsrt[i].getKappaMax(T)			
					kappa_min = self.unsrt[i].getKappaMin(T)
					Theta = np.array([T/T,np.log(T),-1/T])
					P = self.unsrt[i].getMean()	
					cov = self.unsrt[i].getCov()
					zeta_A = self.unsrt[i].zeta.x
					
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					"""
					Type-A zeta samples
					"""
					generators_a = g_a
					parallel_zetas_a = []
					
					time_start = time.time()
					
					
					for k in g_a:
						
						parallel_zetas_a.append(k[0]*zeta_A)
					#generators_b = []
					#parallel_zetas_b = []
					generators_c = []
					parallel_zetas_c = []
					start = 0
					for w in range(10):
						stop = int(((w+1)/10)*self.sim_)
						g_c = self.rxn_generators_c2[start:stop]
						data["generators_c"] = g_c
						callWorkForce = Worker(self.allowed_count)	
						generators_c_,parallel_zetas_c_ = callWorkForce.do_unsrt_c(data,len(g_c))
						del callWorkForce
						generators_c.extend(generators_c_)
						parallel_zetas_c.extend(parallel_zetas_c_)
						removed = data.pop("generators_c")
						start = stop
					
					g_c = self.rxn_generators_c2
					temp_a[i] = parallel_zetas_a
					gen_a[i] = generators_a
					temp_c[i] = parallel_zetas_c
					gen_c[i] = generators_c
					
					time_end = time.time()
					print(f"\n\t\tTime taken = {time_end-time_start}\n")
					
					#print(parallel_zetas_b)
					Y_a = parallel_zetas_a
					#Y_b = parallel_zetas_b
					Y_c = parallel_zetas_c
					T_ = np.linspace(300,2500,8)
					theta = np.array([T/T,np.log(T),-1/T])
					kappa_mean = theta.T.dot(P)
					Theta = np.array([T/T,np.log(T),-1/T])
					X = []
					Y = []
					
					for j in Y_c:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						
					
					temp_xdata = {}
					for k in range(len(Y)):
						temp_xdata[str(Y[k])] = X[k]
						
					random.shuffle(Y)
					#print(Y_)
					X_ = []
					for k in range(len(Y)):
						X_.append(temp_xdata[str(Y[k])])
						
					#print(X_[9])
					xdata[i] = X_
					ydata[i] = Y
					#print(len(parallel_zetas))
					#print(len(self.rxn_generators))
				A = []
				A_gen = []
				X_dash = []
				y_dash = []
								
				for j in range(len(g_c)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(ydata[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				
				
				
				self.beta_ = np.asarray(A)
				#print(self.beta_)
				self.generators = np.asarray(A_gen)
				self.X = np.asarray(X_dash)
				self.y = np.asarray(y_dash)
			
			
			elif self.design == "A1+B1+C1":
				self.beta_=[]
				task = "complete"
				if "Data" in os.listdir():
					os.chdir("Data/Simulations")
					if "Beta_list_case-0.csv" in os.listdir():
						#print(great)
						#print(os.getcwd())
						
						A_file = open("Beta_list_case-0.csv","r").readlines() 
						self.sim_ = len(A_file)
						print(self.sim_)
						A = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in A_file])	
						A_gen_file = open("Generator_list_case-0.csv","r").readlines() 
						A_gen = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in A_gen_file])	
						
						
						X_file = open("NormalizedBeta_list_case-0.csv","r").readlines() 
						X_dash = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in A_file])
						y_dash = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in A_file])
						self.beta_ = np.asarray(A)
						#print(self.beta_)
						self.generators = np.asarray(A_gen)
						
						self.X = np.asarray(X_dash)
						self.y = np.asarray(y_dash)
						os.chdir("../..")
						#print(self.X[0])
						#raise AssertionError("Stop!!")
					else:
						os.chdir("../..")
						task = "incomplete"
				if task != "complete":
					self.beta_=[]
					temp_a = {}
					temp_b = {}
					temp_c = {}
					gen = {}
					gen_a = {}
					gen_b = {}
					gen_c = {}
					xdata = {}
					ydata = {}
					
					temp_a_ = {}
					temp_b_ = {}
					temp_c_ = {}
					gen_ = {}
					gen_a_ = {}
					gen_b_ = {}
					gen_c_ = {}
					
					xdata_ = {}
					ydata_ = {}				
					self.generator_list = []
					"""
					Min- simulation = 300
					"""
					#if self.sim_ < 300:
					#	self.sim_ = 300
					#print(self.sim_)
					"""
					Dividing the generators into three bins
					------------------------------------
					n_a = 0.33*self.sim_
					n_b = 0.33*self.sim_
					n_c = 0.34*self.sim_
					
					"""
					shuffled = int(0.8*self.sim_)
					unshuffled = self.sim_ - shuffled
					
					n_a = int(0.1*shuffled)
					n_b = int(0.45*shuffled)
					n_c = shuffled-n_a-n_b
					
					n_a_ = int(0.1*unshuffled)
					n_b_ = int(0.45*unshuffled)
					n_c_ = unshuffled-n_a_-n_b_
					
					#n_b = int(0.5*self.sim_)
					#n_c = self.sim_-n_b
					
					
					#for i in np.eye(self.n):
					#	self.generator_list.append(i)
					#	self.generator_list.append(-1*i)
									
					"""
					Mixing three types of the extreme curves
					--------------------
					zeta-A: 1
					zeta-B: 1
					zeta-C: 1
				
					for zeta-A: one normal random parameter is required
						a matrix on n_a x n will be created
						
					for zeta-B and zeta-c: two random parameters is required
						we can arrange the two normal_random_parameter 				
						a matrix of n_b+n_c x 2*n will be created
						
					The number of reactions = n_R
					"""
					
					
					n_rxn = len(self.ind)
					self.generator_list_A = []
					self.generator_list_B = []
					self.generator_list_C = []
					
					self.generator_list_A_ = []
					self.generator_list_B_ = []
					self.generator_list_C_ = []
					
					"""
					Factorial Design
					----------------
					box-benkhen design
					bb = []
					for i in bbdesign(n_rxn)[0:int(0.3*self.sim_)]:
						bb.append(np.asarray(i))
					"""
							#For zeta-A
					"""		
					bb_a = []
					for i in bbdesign(n_rxn)[0:int(0.4*n_a)]:
						bb_a.append(np.asarray(i))
					self.generator_list_A.extend(bb_a)
					
							#For zeta-B
							
					bb_b = []
					for i in bbdesign(2*n_rxn)[0:int(0.4*n_b)]:
						bb_b.append(np.asarray(i))							
					self.generator_list_B.extend(bb_b)	
					
							#For zeta-C
							
					bb_c = []
					for i in bbdesign(2*n_rxn)[0:int(0.4*n_c)]:
						bb_c.append(np.asarray(i))			
					self.generator_list_C.extend(bb_c)		
							
					"""
					"""
					Monte-Carlo
					
					"""
							
							#For zeta-A
					
					self.generator_list_A.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(n_a)])
					self.generator_list_A.extend(list(np.eye(n_rxn)))
					self.generator_list_A.extend(list(-1*np.eye(n_rxn)))
					self.generator_list_A.append(np.zeros(n_rxn))
					self.generator_list_A.append(np.zeros(n_rxn))
					self.generator_list_A.append(np.zeros(n_rxn))
					
							#For zeta-b
							
					self.generator_list_B.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(n_b)])	
					self.generator_list_B.extend(list(np.eye(2*n_rxn)))
					self.generator_list_B.extend(list(-1*np.eye(2*n_rxn)))
					self.generator_list_B.append(np.zeros(2*n_rxn))
					self.generator_list_B.append(np.zeros(2*n_rxn))
					self.generator_list_B.append(np.zeros(2*n_rxn))
							
							#For zeta-C		
					self.generator_list_C.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(n_c)])
					self.generator_list_C.extend(list(np.eye(2*n_rxn)))
					self.generator_list_C.extend(list(-1*np.eye(2*n_rxn)))
					self.generator_list_C.append(np.zeros(2*n_rxn))
					self.generator_list_C.append(np.zeros(2*n_rxn))
					self.generator_list_C.append(np.zeros(2*n_rxn))
					
					
					self.generator_list_A_.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(n_a_)])
					#self.generator_list_A_.extend(list(np.eye(n_rxn)))
					#self.generator_list_A_.extend(list(-1*np.eye(n_rxn)))
					#self.generator_list_A_.append(np.zeros(n_rxn))
					#self.generator_list_A_.append(np.zeros(n_rxn))
					#self.generator_list_A_.append(np.zeros(n_rxn))
					
							#For zeta-A
							
					self.generator_list_B_.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(n_b_)])	
					#self.generator_list_B_.extend(list(np.eye(2*n_rxn)))
					#self.generator_list_B_.extend(list(-1*np.eye(2*n_rxn)))
					#self.generator_list_B_.append(np.zeros(2*n_rxn))
					#self.generator_list_B_.append(np.zeros(2*n_rxn))
					#self.generator_list_B_.append(np.zeros(2*n_rxn))
							
							#For zeta-A		
					self.generator_list_C_.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(n_c_)])
					#self.generator_list_C_.extend(list(np.eye(2*n_rxn)))
					#self.generator_list_C_.extend(list(-1*np.eye(2*n_rxn)))
					#self.generator_list_C_.append(np.zeros(2*n_rxn))
					#self.generator_list_C_.append(np.zeros(2*n_rxn))
					#self.generator_list_C_.append(np.zeros(2*n_rxn))
					
					
					
					"""
					Uniform
					
					"""
						
					#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
					
					"""
					Latin-Hypercube
					"""
					if n_a < 10:
						n_a = 30
					if n_a_ < 10:
						n_a_ = 30
							
							#For zeta-A
					"""		
					lhs_a = []
					for i in lhs(n_rxn, samples=int(0.2*n_a)):
						lhs_a.append(np.asarray(i))	
					self.generator_list_A.extend(lhs_a)
					
							
							#For zeta-A
							
					lhs_b = []
					for i in lhs(2*n_rxn, samples=int(0.3*n_b)):
						lhs_b.append(np.asarray(i))		
					self.generator_list_B.extend(lhs_b)
							
							#For zeta-A
							
					lhs_c = []
					for i in lhs(2*n_rxn, samples=int(0.3*n_c)):
						lhs_c.append(np.asarray(i))
					self.generator_list_C.extend(lhs_c)
					
					"""
					"""
					Monte-Carlo-Mix
					"""
					#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)		
					#self.generator_list.extend(bb)
					#self.generator_list.extend(list(np.eye(self.n)))
					#self.generator_list.extend(list(-1*np.eye(self.n)))
					#self.generator_list.append(np.zeros(self.n))
					#self.generator_list.append(np.zeros(self.n))
					#self.generator_list.append(np.zeros(self.n))
					#print(self.generator_list)
					#self.sim_ = len(self.generator_list)
					
					#trainers = np.asarray(trainers)
					#print(f"Generator_list {np.shape(self.generator_list)}\n")
					df = pd.DataFrame(np.asarray(self.generator_list))
					count_start_a = 0
					count_start_b = 0
					count_start_c = 0			
					for i in self.reactions_selected:
						
						len_of_Arrhenius_A = 1
						count_end_a = count_start_a+len_of_Arrhenius_A
						
						
						len_of_Arrhenius_B = 2
						count_end_b = count_start_b+len_of_Arrhenius_B
						
						
						len_of_Arrhenius_C = 2
						count_end_c = count_start_c+len_of_Arrhenius_C
						
						#print(count_start)
						#print(count_end)
						self.rxn_generators_a2 = []
						self.rxn_generators_b2 = []
						self.rxn_generators_c2 = []
						
						self.rxn_generators_a2_ = []
						self.rxn_generators_b2_ = []
						self.rxn_generators_c2_ = []
						
						"""
						getting zetas from generators
						"""
						for j in range(len(self.generator_list_A)):
							self.rxn_generators_a2.append(self.generator_list_A[j][count_start_a:count_end_a])
						#print(self.rxn_generators_a2)
						for j in range(len(self.generator_list_B)):
							self.rxn_generators_b2.append(self.generator_list_B[j][count_start_b:count_end_b])
						
						for j in range(len(self.generator_list_C)):
							self.rxn_generators_c2.append(self.generator_list_C[j][count_start_c:count_end_c])
						
						
						for j in range(len(self.generator_list_A_)):
							self.rxn_generators_a2_.append(self.generator_list_A_[j][count_start_a:count_end_a])
						#print(self.rxn_generators_a2)
						for j in range(len(self.generator_list_B_)):
							self.rxn_generators_b2_.append(self.generator_list_B_[j][count_start_b:count_end_b])
						
						for j in range(len(self.generator_list_C_)):
							self.rxn_generators_c2_.append(self.generator_list_C_[j][count_start_c:count_end_c])
						
						
						
						
						self.sim_ = len(self.rxn_generators_a2)+len(self.rxn_generators_b2)+len(self.rxn_generators_c2)+len(self.rxn_generators_a2_)+len(self.rxn_generators_b2_)+len(self.rxn_generators_c2_)
						#print(self.rxn_generators)
						#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
											
						g_a = self.rxn_generators_a2
						g_b = self.rxn_generators_b2
						g_c = self.rxn_generators_c2
						
						g_a_ = self.rxn_generators_a2_
						g_b_ = self.rxn_generators_b2_
						g_c_ = self.rxn_generators_c2_
						
						#g_b = self.rxn_generators[0:n_b]
						#g_c = self.rxn_generators[n_b:n_b+n_c]
						
						
						count_start_a = count_end_a
						count_start_b = count_end_b
						count_start_c = count_end_c
						
						"""
						Getting the uncertainty data
						"""
						#fig = plt.figure()
						data = self.unsrt[i].data
						T = np.linspace(300,2500,8)
						kappa_max = self.unsrt[i].getKappaMax(T)			
						kappa_min = self.unsrt[i].getKappaMin(T)
						Theta = np.array([T/T,np.log(T),-1/T])
						P = self.unsrt[i].getMean()	
						cov = self.unsrt[i].getCov()
						zeta_A = self.unsrt[i].zeta.x
						
						#plt.plot(1/T,kappa_max,"k--")
						#plt.plot(1/T,kappa_min,"k--")
						"""
						Type-A zeta samples
						"""
						generators_a = g_a
						generators_a_ = g_a_
						parallel_zetas_a = []
						parallel_zetas_a_ = []
						
						time_start = time.time()
						
						
						for k in g_a:
							
							parallel_zetas_a.append(k[0]*zeta_A)
						
						for k in g_a_:
							
							parallel_zetas_a_.append(k[0]*zeta_A)
											
						
						generators_b_ = []
						parallel_zetas_b_ = []
						generators_c_ = []
						parallel_zetas_c_ = []
						
						generators_b = []
						parallel_zetas_b = []
						generators_c = []
						parallel_zetas_c = []
						start = 0
						start_ = 0
						for w in range(4):
							stop = int(((w+1)/4)*shuffled)
							stop_ = int(((w+1)/4)*unshuffled)
							g_b = self.rxn_generators_b2[start:stop]
							g_c = self.rxn_generators_c2[start:stop]
							g_b_ = self.rxn_generators_b2_[start_:stop_]
							g_c_ = self.rxn_generators_c2_[start_:stop_]
							data["generators_b"] = g_b
							callWorkForce = Worker(self.allowed_count)	
							gb,pb = callWorkForce.do_unsrt_b(data,len(g_b))
							del callWorkForce
							
							generators_b.extend(gb)
							parallel_zetas_b.extend(pb)
							
							removed = data.pop("generators_b")
							data["generators_b"] = g_b_
							callWorkForce = Worker(self.allowed_count)	
							gb,pb = callWorkForce.do_unsrt_b(data,len(g_b_))
							del callWorkForce
							
							generators_b_.extend(gb)
							parallel_zetas_b_.extend(pb)
							
							removed = data.pop("generators_b")
							data["generators_c"] = g_c
							callWorkForce = Worker(self.allowed_count)	
							gc,pc = callWorkForce.do_unsrt_c(data,len(g_c))
							del callWorkForce
							
							generators_c.extend(gc)
							parallel_zetas_c.extend(pc)
							removed = data.pop("generators_c")
							
							data["generators_c"] = g_c_
							callWorkForce = Worker(self.allowed_count)	
							gc,pc = callWorkForce.do_unsrt_c(data,len(g_c_))
							del callWorkForce
							
							generators_c_.extend(gc)
							parallel_zetas_c_.extend(pc)
							removed = data.pop("generators_c")
							
							start = stop
							start_ = stop_
						
						#print(self.sim_)
						#print(len(parallel_zetas_b)+len(parallel_zetas_c))
						g_b = self.rxn_generators_b2
						g_c = self.rxn_generators_c2
						
						#print(len(g_b))
						#print(len(parallel_zetas_b))
						
						g_b_ = self.rxn_generators_b2_
						g_c_ = self.rxn_generators_c2_
						
						
						temp_a[i] = parallel_zetas_a
						gen_a[i] = generators_a
						temp_b[i] = parallel_zetas_b
						gen_b[i] = generators_b
						temp_c[i] = parallel_zetas_c
						gen_c[i] = generators_c
						
						temp_a_[i] = parallel_zetas_a_
						gen_a_[i] = generators_a_
						temp_b_[i] = parallel_zetas_b_
						gen_b_[i] = generators_b_
						temp_c_[i] = parallel_zetas_c_
						gen_c_[i] = generators_c_
						
						time_end = time.time()
						print(f"\n\t\tTime taken = {time_end-time_start}\n")
						
						#print(parallel_zetas_b)
						Y_a = parallel_zetas_a
						Y_b = parallel_zetas_b
						Y_c = parallel_zetas_c
						
						Y_a_ = parallel_zetas_a_
						Y_b_ = parallel_zetas_b_
						Y_c_ = parallel_zetas_c_
						
						T_ = np.linspace(300,2500,8)
						theta = np.array([T/T,np.log(T),-1/T])
						kappa_mean = theta.T.dot(P)
						Theta = np.array([T/T,np.log(T),-1/T])
						X = []
						Y = []
						
						X_ = []
						Y_ = []
						
						
						for j in Y_a:
							Y.append(j)
							Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
							
							X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
							
						for j in Y_a_:
							Y_.append(j)
							Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
							
							X_.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
							
						
						for j in Y_b:
							Y.append(j)
							Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
							#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
							X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						
						for j in Y_b_:
							Y_.append(j)
							Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
							#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
							X_.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
							#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
						
						for j in Y_c:
							Y.append(j)
							Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
							#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
							X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
							#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
						#plt.savefig(str(i)+"_c.pdf")
						
						for j in Y_c_:
							Y_.append(j)
							Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
							#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
							X_.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
							#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
						#plt.savefig(str(i)+"_c.pdf")
						
						temp_xdata = {}
						for k in range(len(Y_)):
							temp_xdata[str(Y_[k])] = X_[k]
							
						random.shuffle(Y_)
						#print(Y_)
						X1_ = []
						for k in range(len(Y_)):
							X1_.append(temp_xdata[str(Y_[k])])
							
						#print(X_[9])
						xdata[i] = X
						ydata[i] = Y
						
						xdata_[i] = X1_
						ydata_[i] = Y_
						
						#print(len(X1_))
						#print(len(Y_))
						#print(len(X))
						#print(len(Y))
						print(self.sim_)
						#Was trying to shuffle the zetas but the try was not successfull
						"""
						Y_ = []
						a = int(unshuffled/3)
						b = unshuffled - 2*a
						for w in Y_a[0:a]:
							Y_.append(w)
						for w in Y_b[0:a]:
							Y_.append(w)
						for w in Y_c[0:b]:
							Y_.append(w)
						
						random.shuffle(Y_)
						ydata_[i]=Y_
						
						#ydata[i] = Y
						#print(len(Y))
						#print(len(Y_a)+len(Y_b)+len(Y_c)))
						extra = len(Y_a[0:a])+len(Y_b[0:a])+len(Y_c[0:b])
						#print(len(parallel_zetas))
						#print(len(self.rxn_generators))
						"""
					A = []
					A_gen = []
					X_dash = []
					y_dash = []
					
					for j in range(len(g_a)):
						a_row = []
						a_row_gen = []
						x_temp = []
						y_temp = []
						for i in temp_a:
							a_row.extend(list(temp_a[i][j]))
							a_row_gen.extend(list(gen_a[i][j]))
							x_temp.extend(list(xdata[i][j]))
							y_temp.extend(list(ydata[i][j]))
						A.append(np.asarray(a_row))
						A_gen.append(np.asarray(a_row_gen))
						X_dash.append(np.asarray(x_temp))
						y_dash.append(np.asarray(y_temp))
						
					for j in range(len(g_a_)):
						a_row = []
						a_row_gen = []
						x_temp = []
						y_temp = []
						for i in temp_a_:
							a_row.extend(list(ydata_[i][j]))
							a_row_gen.extend(list(gen_a_[i][j]))
							x_temp.extend(list(xdata_[i][j]))
							y_temp.extend(list(ydata_[i][j]))
						A.append(np.asarray(a_row))
						A_gen.append(np.asarray(a_row_gen))
						X_dash.append(np.asarray(x_temp))
						y_dash.append(np.asarray(y_temp))
					
					
					for j in range(len(g_b)):
						a_row = []
						a_row_gen = []
						x_temp = []
						y_temp = []
						
						for i in temp_b:
							#print(temp_b[i])
							#print(len(g_b))
							#print(len(temp_b[i]))
							a_row.extend(list(temp_b[i][j]))
							a_row_gen.extend(list(gen_b[i][j]))
							x_temp.extend(list(xdata[i][j+len(g_a)]))
							y_temp.extend(list(ydata[i][j+len(g_a)]))
						A.append(np.asarray(a_row))
						A_gen.append(np.asarray(a_row_gen))
						X_dash.append(np.asarray(x_temp))
						y_dash.append(np.asarray(y_temp))
					
					for j in range(len(g_b_)):
						a_row = []
						a_row_gen = []
						x_temp = []
						y_temp = []
						for i in temp_b_:
							a_row.extend(list(ydata_[i][j+len(g_a_)]))
							a_row_gen.extend(list(gen_b_[i][j]))
							x_temp.extend(list(xdata_[i][j+len(g_a_)]))
							y_temp.extend(list(ydata_[i][j+len(g_a_)]))
						A.append(np.asarray(a_row))
						A_gen.append(np.asarray(a_row_gen))
						X_dash.append(np.asarray(x_temp))
						y_dash.append(np.asarray(y_temp))
					"""
					for j in range(n_b):
						a_row = []
						a_row_gen = []
						x_temp = []
						y_temp = []
						for i in temp_b:
							a_row.extend(list(temp_b[i][j]))
							a_row_gen.extend(list(gen_b[i][j]))
							x_temp.extend(list(xdata[i][j]))
							y_temp.extend(list(ydata[i][j]))
						A.append(np.asarray(a_row))
						A_gen.append(np.asarray(a_row_gen))
						X_dash.append(np.asarray(x_temp))
						y_dash.append(np.asarray(y_temp))
					"""
					
					for j in range(len(g_c)):
						a_row = []
						a_row_gen = []
						x_temp = []
						y_temp = []
						for i in temp_c:
							print(len(temp_c[i]))
							print(len(g_c))
							a_row.extend(list(temp_c[i][j]))
							a_row_gen.extend(list(gen_c[i][j]))
							x_temp.extend(list(xdata[i][j+len(g_a)+len(g_b)]))
							y_temp.extend(list(ydata[i][j+len(g_a)+len(g_b)]))
						A.append(np.asarray(a_row))
						A_gen.append(np.asarray(a_row_gen))
						X_dash.append(np.asarray(x_temp))
						y_dash.append(np.asarray(y_temp))
					
					#print(len())
					for j in range(len(g_c_)):
						a_row = []
						a_row_gen = []
						x_temp = []
						y_temp = []
						for i in temp_c_:
							#print(len(ydata_[i]))
							#print(len(g_a_)+len(g_b_)+len(g_c_))
							a_row.extend(list(ydata_[i][j+len(g_a_)+len(g_b_)]))
							a_row_gen.extend(list(gen_c_[i][j]))
							x_temp.extend(list(xdata_[i][j+len(g_a_)+len(g_b_)]))
							y_temp.extend(list(ydata_[i][j+len(g_a_)+len(g_b_)]))
						A.append(np.asarray(a_row))
						A_gen.append(np.asarray(a_row_gen))
						X_dash.append(np.asarray(x_temp))
						y_dash.append(np.asarray(y_temp))
				
					
					#print(A)
					#self.sim_ = len(g_a)+len(g_b)+len(g_c)
					#self.sim_ = len(g_a)+len(g_b)+len(g_c)
					#self.sim_ = n_b+n_c
					self.beta_ = np.asarray(A)
					#print(self.beta_)
					self.generators = np.asarray(A_gen)
					self.X = np.asarray(X_dash)
					self.y = np.asarray(y_dash)
			
			elif self.design == "B1-factorial":
				self.beta_=[]
				temp_a = {}
				temp_b = {}
				temp_c = {}
				gen = {}
				gen_a = {}
				gen_b = {}
				gen_c = {}
				
				xdata = {}
				ydata = {}				
				self.generator_list = []
				"""
				Min- simulation = 300
				"""
				#if self.sim_ < 300:
				#	self.sim_ = 300
				#print(self.sim_)
				"""
				Dividing the generators into three bins
				------------------------------------
				n_a = 0.33*self.sim_
				n_b = 0.33*self.sim_
				n_c = 0.34*self.sim_
				
				"""
				n_a = int(0.33*self.sim_)
				n_b = int(0.33*self.sim_)
				n_c = self.sim_-n_a-n_b
				
				
				#n_b = int(0.5*self.sim_)
				#n_c = self.sim_-n_b
				
				
				#for i in np.eye(self.n):
				#	self.generator_list.append(i)
				#	self.generator_list.append(-1*i)
								
				"""
				Mixing three types of the extreme curves
				--------------------
				zeta-A: 1
				zeta-B: 1
				zeta-C: 1
			
				for zeta-A: one normal random parameter is required
					a matrix on n_a x n will be created
					
				for zeta-B and zeta-c: two random parameters is required
					we can arrange the two normal_random_parameter 				
					a matrix of n_b+n_c x 2*n will be created
					
				The number of reactions = n_R
				"""
				n_rxn = len(self.ind)
				self.generator_list_A = []
				self.generator_list_B = []
				self.generator_list_C = []
				
				"""
				Factorial Design
				----------------
				box-benkhen design
				bb = []
				for i in bbdesign(n_rxn)[0:int(0.3*self.sim_)]:
					bb.append(np.asarray(i))
				"""
						#For zeta-A
						
				bb_a = []
				for i in bbdesign(n_rxn):
					bb_a.append(np.asarray(i))
				self.generator_list_A.extend(bb_a)
				self.generator_list_A.extend(list(np.eye(n_rxn)))
				self.generator_list_A.extend(list(-1*np.eye(n_rxn)))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))		
						#For zeta-B
						
				bb_b = []
				for i in bbdesign(2*n_rxn):
					bb_b.append(np.asarray(i))							
				self.generator_list_B.extend(bb_b)
				self.generator_list_B.extend(list(np.eye(2*n_rxn)))
				self.generator_list_B.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))		
						#For zeta-C
						
				bb_c = []
				for i in bbdesign(2*n_rxn):
					bb_c.append(np.asarray(i))		
						
				self.generator_list_C.extend(bb_c)		
				self.generator_list_C.extend(list(np.eye(2*n_rxn)))
				self.generator_list_C.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))		
				
				"""
				Monte-Carlo
				
				"""
				"""		
						#For zeta-A
				
				self.generator_list_A.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(0.3*n_a)])
				
						#For zeta-A
						
				self.generator_list_B.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_a)])	
						
						#For zeta-A		
				self.generator_list_C.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_a)])
				
				"""
				"""
				Uniform
				
				"""
					
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				
				"""
				Latin-Hypercube
				"""
				"""
						
						#For zeta-A
						
				lhs_a = []
				for i in lhs(n_rxn, samples=int(0.3*n_a)):
					lhs_a.append(np.asarray(i))	
				self.generator_list_A.extend(lhs_a)
				
						
						#For zeta-A
						
				lhs_b = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_b)):
					lhs_b.append(np.asarray(i))		
				self.generator_list_B.extend(lhs_b)
						
						#For zeta-A
						
				lhs_c = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_c)):
					lhs_c.append(np.asarray(i))
				self.generator_list_C.extend(lhs_c)
				
				"""
				"""
				Monte-Carlo-Mix
				"""
				
				
				
				
				
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)		
				#self.generator_list.extend(bb)
				#self.generator_list.extend(list(np.eye(self.n)))
				#self.generator_list.extend(list(-1*np.eye(self.n)))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#print(self.generator_list)
				#self.sim_ = len(self.generator_list)
				
				#trainers = np.asarray(trainers)
				#print(f"Generator_list {np.shape(self.generator_list)}\n")
				df = pd.DataFrame(np.asarray(self.generator_list))
				count_start_a = 0
				count_start_b = 0
				count_start_c = 0			
				for i in self.reactions_selected:
					
					len_of_Arrhenius_A = 1
					count_end_a = count_start_a+len_of_Arrhenius_A
					
					
					len_of_Arrhenius_B = 2
					count_end_b = count_start_b+len_of_Arrhenius_B
					
					
					len_of_Arrhenius_C = 2
					count_end_c = count_start_c+len_of_Arrhenius_C
					
					#print(count_start)
					#print(count_end)
					self.rxn_generators_a2 = []
					self.rxn_generators_b2 = []
					self.rxn_generators_c2 = []
					
					"""
					getting zetas from generators
					"""
					for j in range(len(self.generator_list_A)):
						self.rxn_generators_a2.append(self.generator_list_A[j][count_start_a:count_end_a])
					#print(self.rxn_generators_a2)
					for j in range(len(self.generator_list_B)):
						self.rxn_generators_b2.append(self.generator_list_B[j][count_start_b:count_end_b])
					
					for j in range(len(self.generator_list_C)):
						self.rxn_generators_c2.append(self.generator_list_C[j][count_start_c:count_end_c])
					
					self.sim_ = len(self.rxn_generators_a2)+len(self.rxn_generators_b2)+len(self.rxn_generators_c2)
					#print(self.rxn_generators)
					#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
										
					g_a = self.rxn_generators_a2
					g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					#g_b = self.rxn_generators[0:n_b]
					#g_c = self.rxn_generators[n_b:n_b+n_c]
					
					
					count_start_a = count_end_a
					count_start_b = count_end_b
					count_start_c = count_end_c
					
					"""
					Getting the uncertainty data
					"""
					#fig = plt.figure()
					data = self.unsrt[i].data
					T = np.linspace(300,2500,8)
					kappa_max = self.unsrt[i].getKappaMax(T)			
					kappa_min = self.unsrt[i].getKappaMin(T)
					Theta = np.array([T/T,np.log(T),-1/T])
					P = self.unsrt[i].getMean()	
					cov = self.unsrt[i].getCov()
					zeta_A = self.unsrt[i].zeta.x
					
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					"""
					Type-A zeta samples
					"""
					generators_a = g_a
					parallel_zetas_a = []
					
					time_start = time.time()
					
					
					for k in g_a:
						#Pint = P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten()
					#	plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
					#	plt.plot()
						#print(i)
						#print(type(i))
						#print(k)	
						#parallel_zetas_a.append(P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten())
						parallel_zetas_a.append(k[0]*zeta_A)
										
					#print(parallel_zetas_a)
					#plt.savefig(str(i)+"_a.pdf")
					#print(os.getcwd())
					#print("plot saved")
					data["generators_b"] = g_b
					callWorkForce = Worker(self.allowed_count)	
					generators_b,parallel_zetas_b = callWorkForce.do_unsrt_b(data,len(g_b))
					del callWorkForce
					
					data["generators_c"] = g_c
					callWorkForce = Worker(self.allowed_count)	
					generators_c,parallel_zetas_c = callWorkForce.do_unsrt_c(data,len(g_c))
					del callWorkForce
					
					#print(parallel_zetas_a)
					
					temp_a[i] = parallel_zetas_a
					gen_a[i] = generators_a
					temp_b[i] = parallel_zetas_b
					gen_b[i] = generators_b
					temp_c[i] = parallel_zetas_c
					gen_c[i] = generators_c
					
					time_end = time.time()
					print(f"\n\t\tTime taken = {time_end-time_start}\n")
					
					#print(parallel_zetas_b)
					Y_a = parallel_zetas_a
					Y_b = parallel_zetas_b
					Y_c = parallel_zetas_c
					T_ = np.linspace(300,2500,8)
					theta = np.array([T/T,np.log(T),-1/T])
					kappa_mean = theta.T.dot(P)
					Theta = np.array([T/T,np.log(T),-1/T])
					X = []
					Y = []
					
					
					for j in Y_a:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					
					for j in Y_b:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_b.pdf")
					
					
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					for j in Y_c:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_c.pdf")
					
					temp_xdata = {}
					for k in range(len(Y)):
						temp_xdata[str(Y[k])] = X[k]
						
					random.shuffle(Y)
					#print(Y_)
					X_ = []
					for k in range(len(Y)):
						X_.append(temp_xdata[str(Y[k])])
						
					#print(X_[9])
					xdata[i] = X_
					ydata[i] = Y
					#print(len(parallel_zetas))
					#print(len(self.rxn_generators))
				A = []
				A_gen = []
				X_dash = []
				y_dash = []
				
				for j in range(len(g_a)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_a:
						a_row.extend(list(temp_a[i][j]))
						a_row_gen.extend(list(gen_a[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				
				
				for j in range(len(g_b)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)]))
						y_temp.extend(list(ydata[i][j+len(g_a)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				for j in range(n_b):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				
				for j in range(len(g_c)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)+len(g_b)]))
						y_temp.extend(list(ydata[i][j+len(g_a)+len(g_b)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				for j in range(n_c):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+n_b]))
						y_temp.extend(list(ydata[i][j+n_b]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
					
				"""
				
				#print(A)
				self.sim_ = len(g_a)+len(g_b)+len(g_c)
				#self.sim_ = n_b+n_c
				self.beta_ = np.asarray(A)
				#print(self.beta_)
				self.generators = np.asarray(A_gen)
				self.X = np.asarray(X_dash)
				self.y = np.asarray(y_dash)
			
			
			elif self.design == "A1+C1":
				self.beta_=[]
				temp_a = {}
				gen = {}
				gen_a = {}
				xdata = {}
				ydata = {}				
				temp_c = {}
				gen = {}
				gen_c = {}
				self.generator_list = []			
				n_a = int(0.1*self.sim_)
				n_c = int(0.9*self.sim_)
				n_rxn = len(self.ind)
				self.generator_list_A = []
				self.generator_list_C = []
				
				"""
				Factorial Design
				----------------
				box-benkhen design
				bb = []
				for i in bbdesign(n_rxn)[0:int(0.3*self.sim_)]:
					bb.append(np.asarray(i))
				"""
						#For zeta-A
						
				bb_a = []
				for i in bbdesign(n_rxn)[0:int(0.4*n_a)]:
					bb_a.append(np.asarray(i))
				self.generator_list_A.extend(bb_a)
				
				bb_c = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_c)]:
					bb_c.append(np.asarray(i))
				self.generator_list_C.extend(bb_c)
						
				"""
				Monte-Carlo
				
				"""
						#For zeta-A
				self.generator_list_A.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(0.3*n_a)])
				self.generator_list_A.extend(list(np.eye(n_rxn)))
				self.generator_list_A.extend(list(-1*np.eye(n_rxn)))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				
						#For zeta-C
						
				self.generator_list_C.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_c)])
				self.generator_list_C.extend(list(np.eye(2*n_rxn)))
				self.generator_list_C.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				
				"""
				Latin-Hypercube
				"""	
						#For zeta-A
						
				lhs_a = []
				for i in lhs(n_rxn, samples=int(0.3*n_a)):
					lhs_a.append(np.asarray(i))	
				self.generator_list_A.extend(lhs_a)
				df = pd.DataFrame(np.asarray(self.generator_list))
				count_start_a = 0		
				
				
						#For zeta-C
				
						
				lhs_c = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_c)):
					lhs_c.append(np.asarray(i))
				self.generator_list_C.extend(lhs_c)
				
				
				
				"""
				Uniform
				
				"""
				count_start_a = 0
				count_start_b = 0
				count_start_c = 0
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)				
				
				for i in self.reactions_selected:
					
										
					
					len_of_Arrhenius_A = 1
					count_end_a = count_start_a+len_of_Arrhenius_A
					
					len_of_Arrhenius_C = 2
					count_end_c = count_start_c+len_of_Arrhenius_C
					
					#print(count_start)
					#print(count_end)
					self.rxn_generators_a2 = []
					self.rxn_generators_c2 = []
					
					"""
					getting zetas from generators
					"""
					for j in range(len(self.generator_list_A)):
						self.rxn_generators_a2.append(self.generator_list_A[j][count_start_a:count_end_a])
				
					for j in range(len(self.generator_list_C)):
						self.rxn_generators_c2.append(self.generator_list_C[j][count_start_c:count_end_c])
					
					self.sim_ = len(self.rxn_generators_a2)+len(self.rxn_generators_c2)
					#print(self.rxn_generators)
					#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
										
					g_a = self.rxn_generators_a2
					g_c = self.rxn_generators_c2
					#g_b = self.rxn_generators[0:n_b]
					#g_c = self.rxn_generators[n_b:n_b+n_c]
					
					
					count_start_a = count_end_a
					count_start_c = count_end_c
					
					
					"""
					Getting the uncertainty data
					"""
					#fig = plt.figure()
					data = self.unsrt[i].data
					T = np.linspace(300,2500,8)
					kappa_max = self.unsrt[i].getKappaMax(T)			
					kappa_min = self.unsrt[i].getKappaMin(T)
					Theta = np.array([T/T,np.log(T),-1/T])
					P = self.unsrt[i].getMean()	
					cov = self.unsrt[i].getCov()
					zeta_A = self.unsrt[i].zeta.x
					"""
					Type-A zeta samples
					"""
					generators_a = g_a
					parallel_zetas_a = []
					time_start = time.time()
					for k in g_a:
						parallel_zetas_a.append(k[0]*zeta_A)
					
					generators_b = []
					parallel_zetas_b = []
					generators_c = []
					parallel_zetas_c = []
					start = 0
					for w in range(10):
						stop = int(((w+1)/10)*self.sim_)
						g_c = self.rxn_generators_c2[start:stop]
						
						data["generators_c"] = g_c
						callWorkForce = Worker(self.allowed_count)	
						generators_c_,parallel_zetas_c_ = callWorkForce.do_unsrt_c(data,len(g_c))
						del callWorkForce
						
						generators_c.extend(generators_c_)
						parallel_zetas_c.extend(parallel_zetas_c_)
						removed = data.pop("generators_c")
						start = stop
					
				
					g_c = self.rxn_generators_c2
					
					#print(parallel_zetas_a)
					
					temp_a[i] = parallel_zetas_a
					gen_a[i] = generators_a

					temp_c[i] = parallel_zetas_c
					gen_c[i] = generators_c
					
					time_end = time.time()
					print(f"\n\t\tTime taken = {time_end-time_start}\n")
					
					#print(parallel_zetas_b)
					Y_a = parallel_zetas_a
					
					Y_c = parallel_zetas_c
					T_ = np.linspace(300,2500,8)
					theta = np.array([T/T,np.log(T),-1/T])
					kappa_mean = theta.T.dot(P)
					Theta = np.array([T/T,np.log(T),-1/T])
					X = []
					Y = []
					
					for j in Y_a:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
					
					for j in Y_c:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
					
					
					temp_xdata = {}
					for k in range(len(Y)):
						temp_xdata[str(Y[k])] = X[k]
						
					random.shuffle(Y)
					#print(Y_)
					X_ = []
					for k in range(len(Y)):
						X_.append(temp_xdata[str(Y[k])])
						
					#print(X_[9])
					xdata[i] = X_
					ydata[i] = Y
					
				A = []
				A_gen = []
				X_dash = []
				y_dash = []
				for j in range(len(g_a)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_a:
						a_row.extend(list(ydata[i][j]))
						a_row_gen.extend(list(gen_a[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				
				for j in range(len(g_c)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(ydata[i][j+len(g_a)]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)]))
						y_temp.extend(list(ydata[i][j+len(g_a)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				#self.sim_ = n_a
				self.beta_ = np.asarray(A)
				self.generators = np.asarray(A_gen)
				self.X = np.asarray(X_dash)
				self.y = np.asarray(y_dash)
				
			elif self.design == "B1-MC":
				self.beta_=[]
				temp_a = {}
				temp_b = {}
				temp_c = {}
				gen = {}
				gen_a = {}
				gen_b = {}
				gen_c = {}
				
				xdata = {}
				ydata = {}				
				self.generator_list = []
				"""
				Min- simulation = 300
				"""
				#if self.sim_ < 300:
				#	self.sim_ = 300
				#print(self.sim_)
				"""
				Dividing the generators into three bins
				------------------------------------
				n_a = 0.33*self.sim_
				n_b = 0.33*self.sim_
				n_c = 0.34*self.sim_
				
				"""
				n_a = int(0.33*self.sim_)
				n_b = int(0.33*self.sim_)
				n_c = self.sim_-n_a-n_b
				
				
				#n_b = int(0.5*self.sim_)
				#n_c = self.sim_-n_b
				
				
				#for i in np.eye(self.n):
				#	self.generator_list.append(i)
				#	self.generator_list.append(-1*i)
								
				"""
				Mixing three types of the extreme curves
				--------------------
				zeta-A: 1
				zeta-B: 1
				zeta-C: 1
			
				for zeta-A: one normal random parameter is required
					a matrix on n_a x n will be created
					
				for zeta-B and zeta-c: two random parameters is required
					we can arrange the two normal_random_parameter 				
					a matrix of n_b+n_c x 2*n will be created
					
				The number of reactions = n_R
				"""
				n_rxn = len(self.ind)
				self.generator_list_A = []
				self.generator_list_B = []
				self.generator_list_C = []
				
				"""
				Factorial Design
				----------------
				box-benkhen design
				bb = []
				for i in bbdesign(n_rxn)[0:int(0.3*self.sim_)]:
					bb.append(np.asarray(i))
				"""
				"""		#For zeta-A
						
				bb_a = []
				for i in bbdesign(n_rxn)[0:int(0.4*n_a)]:
					bb_a.append(np.asarray(i))
				self.generator_list_A.extend(bb_a)
						
						#For zeta-B
						
				bb_b = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_b)]:
					bb_b.append(np.asarray(i))							
				self.generator_list_B.extend(bb_b)
						
						#For zeta-C
						
				bb_c = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_c)]:
					bb_c.append(np.asarray(i))		
						
				self.generator_list_C.extend(bb_c)		
						
				"""
				"""
				Monte-Carlo
				
				"""
						
						#For zeta-A
				
				self.generator_list_A.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(0.3*n_a)])
				#self.generator_list_A.extend(list(np.eye(n_rxn)))
				#self.generator_list_A.extend(list(-1*np.eye(n_rxn)))
				#self.generator_list_A.append(np.zeros(n_rxn))
				#self.generator_list_A.append(np.zeros(n_rxn))
				#self.generator_list_A.append(np.zeros(n_rxn))
				
						#For zeta-A
						
				self.generator_list_B.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_a)])	
				#self.generator_list_B.extend(list(np.eye(2*n_rxn)))
				#self.generator_list_B.extend(list(-1*np.eye(2*n_rxn)))
				#self.generator_list_B.append(np.zeros(2*n_rxn))
				#self.generator_list_B.append(np.zeros(2*n_rxn))
				#self.generator_list_B.append(np.zeros(2*n_rxn))
						
						#For zeta-A		
				self.generator_list_C.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_a)])
				#self.generator_list_C.extend(list(np.eye(2*n_rxn)))
				#self.generator_list_C.extend(list(-1*np.eye(2*n_rxn)))
				#self.generator_list_C.append(np.zeros(2*n_rxn))
				#self.generator_list_C.append(np.zeros(2*n_rxn))
				#self.generator_list_C.append(np.zeros(2*n_rxn))
				
				"""
				Uniform
				
				"""
					
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				
				"""
				Latin-Hypercube
				"""
				"""
						
						#For zeta-A
						
				lhs_a = []
				for i in lhs(n_rxn, samples=int(0.3*n_a)):
					lhs_a.append(np.asarray(i))	
				self.generator_list_A.extend(lhs_a)
				
						
						#For zeta-A
						
				lhs_b = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_b)):
					lhs_b.append(np.asarray(i))		
				self.generator_list_B.extend(lhs_b)
						
						#For zeta-A
						
				lhs_c = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_c)):
					lhs_c.append(np.asarray(i))
				self.generator_list_C.extend(lhs_c)
				
				"""
				"""
				Monte-Carlo-Mix
				"""
				
				
				
				
				
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)		
				#self.generator_list.extend(bb)
				#self.generator_list.extend(list(np.eye(self.n)))
				#self.generator_list.extend(list(-1*np.eye(self.n)))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#print(self.generator_list)
				#self.sim_ = len(self.generator_list)
				
				#trainers = np.asarray(trainers)
				#print(f"Generator_list {np.shape(self.generator_list)}\n")
				df = pd.DataFrame(np.asarray(self.generator_list))
				count_start_a = 0
				count_start_b = 0
				count_start_c = 0			
				for i in self.reactions_selected:
					
					len_of_Arrhenius_A = 1
					count_end_a = count_start_a+len_of_Arrhenius_A
					
					
					len_of_Arrhenius_B = 2
					count_end_b = count_start_b+len_of_Arrhenius_B
					
					
					len_of_Arrhenius_C = 2
					count_end_c = count_start_c+len_of_Arrhenius_C
					
					#print(count_start)
					#print(count_end)
					self.rxn_generators_a2 = []
					self.rxn_generators_b2 = []
					self.rxn_generators_c2 = []
					
					"""
					getting zetas from generators
					"""
					for j in range(len(self.generator_list_A)):
						self.rxn_generators_a2.append(self.generator_list_A[j][count_start_a:count_end_a])
					#print(self.rxn_generators_a2)
					for j in range(len(self.generator_list_B)):
						self.rxn_generators_b2.append(self.generator_list_B[j][count_start_b:count_end_b])
					
					for j in range(len(self.generator_list_C)):
						self.rxn_generators_c2.append(self.generator_list_C[j][count_start_c:count_end_c])
					
					self.sim_ = len(self.rxn_generators_a2)+len(self.rxn_generators_b2)+len(self.rxn_generators_c2)
					#print(self.rxn_generators)
					#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
										
					g_a = self.rxn_generators_a2
					g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					#g_b = self.rxn_generators[0:n_b]
					#g_c = self.rxn_generators[n_b:n_b+n_c]
					
					
					count_start_a = count_end_a
					count_start_b = count_end_b
					count_start_c = count_end_c
					
					"""
					Getting the uncertainty data
					"""
					#fig = plt.figure()
					data = self.unsrt[i].data
					T = np.linspace(300,2500,8)
					kappa_max = self.unsrt[i].getKappaMax(T)			
					kappa_min = self.unsrt[i].getKappaMin(T)
					Theta = np.array([T/T,np.log(T),-1/T])
					P = self.unsrt[i].getMean()	
					cov = self.unsrt[i].getCov()
					zeta_A = self.unsrt[i].zeta.x
					
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					"""
					Type-A zeta samples
					"""
					generators_a = g_a
					parallel_zetas_a = []
					
					time_start = time.time()
					
					
					for k in g_a:
						#Pint = P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten()
					#	plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
					#	plt.plot()
						#print(i)
						#print(type(i))
						#print(k)	
						#parallel_zetas_a.append(P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten())
						parallel_zetas_a.append(k[0]*zeta_A)
										
					#print(parallel_zetas_a)
					#plt.savefig(str(i)+"_a.pdf")
					#print(os.getcwd())
					#print("plot saved")
					data["generators_b"] = g_b
					callWorkForce = Worker(self.allowed_count)	
					generators_b,parallel_zetas_b = callWorkForce.do_unsrt_b(data,len(g_b))
					del callWorkForce
					
					data["generators_c"] = g_c
					callWorkForce = Worker(self.allowed_count)	
					generators_c,parallel_zetas_c = callWorkForce.do_unsrt_c(data,len(g_c))
					del callWorkForce
					
					#print(parallel_zetas_a)
					
					temp_a[i] = parallel_zetas_a
					gen_a[i] = generators_a
					temp_b[i] = parallel_zetas_b
					gen_b[i] = generators_b
					temp_c[i] = parallel_zetas_c
					gen_c[i] = generators_c
					
					time_end = time.time()
					print(f"\n\t\tTime taken = {time_end-time_start}\n")
					
					#print(parallel_zetas_b)
					Y_a = parallel_zetas_a
					Y_b = parallel_zetas_b
					Y_c = parallel_zetas_c
					T_ = np.linspace(300,2500,8)
					theta = np.array([T/T,np.log(T),-1/T])
					kappa_mean = theta.T.dot(P)
					Theta = np.array([T/T,np.log(T),-1/T])
					X = []
					Y = []
					
					
					for j in Y_a:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					
					for j in Y_b:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_b.pdf")
					
					
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					for j in Y_c:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_c.pdf")
					
					temp_xdata = {}
					for k in range(len(Y)):
						temp_xdata[str(Y[k])] = X[k]
						
					random.shuffle(Y)
					#print(Y_)
					X_ = []
					for k in range(len(Y)):
						X_.append(temp_xdata[str(Y[k])])
						
					#print(X_[9])
					xdata[i] = X_
					ydata[i] = Y
					#print(len(parallel_zetas))
					#print(len(self.rxn_generators))
				A = []
				A_gen = []
				X_dash = []
				y_dash = []
				
				for j in range(len(g_a)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_a:
						a_row.extend(list(temp_a[i][j]))
						a_row_gen.extend(list(gen_a[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				
				
				for j in range(len(g_b)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)]))
						y_temp.extend(list(ydata[i][j+len(g_a)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				for j in range(n_b):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				
				for j in range(len(g_c)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)+len(g_b)]))
						y_temp.extend(list(ydata[i][j+len(g_a)+len(g_b)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				for j in range(n_c):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+n_b]))
						y_temp.extend(list(ydata[i][j+n_b]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
					
				"""
				
				#print(A)
				self.sim_ = len(g_a)+len(g_b)+len(g_c)
				#self.sim_ = n_b+n_c
				self.beta_ = np.asarray(A)
				#print(self.beta_)
				self.generators = np.asarray(A_gen)
				self.X = np.asarray(X_dash)
				self.y = np.asarray(y_dash)
			
			elif self.design == "B1-Linear":
				self.beta_=[]
				temp = {}
				gen = {}
				xdata = {}
				ydata = {}
				for i in self.reactions_selected:
					self.generator_list = []
					"""
					Getting linear combination of zeta
					"""
					for j in range(self.sim_):
						beta1 = 2*np.random.random_sample(1)-1
						beta2 = 2*np.random.random_sample(1)-1
						m = beta1*np.pi
						n = beta2*np.pi/2
						alpha = m + np.pi
						gamma = n + (np.pi/2)
						x = np.cos(alpha)*np.sin(gamma)
						y = np.sin(alpha)*np.sin(gamma)
						z = np.cos(gamma)
						self.generator_list.append(np.array([x[0],y[0],z[0]]))
					
					
					"""
					getting zetas from generators
					"""
					#for j in range(self.sim_):
					#	self.generator_list.append(2*np.random.random_sample(3)-1)
					data = self.unsrt[i].data
					T = np.linspace(300,2500,8)
					kappa_max = self.unsrt[i].getKappaMax(T)			
					kappa_min = self.unsrt[i].getKappaMin(T)
					P = self.unsrt[i].getMean()	
					zetaMat = self.unsrt[i].zeta_Matrix
					
					parallel_zetas = []
					for k in range(self.sim_):
						zeta = np.asarray(zetaMat.T.dot(self.generator_list[k].T)).flatten()
						parallel_zetas.append(zeta)
			
					
					cov = self.unsrt[i].getCov()
					#data["generators"] = self.generator_list
					#callWorkForce = Worker(10)	
					#generators,parallel_zetas = callWorkForce.do_unsrt(data,int(self.sim_))
					
					gen[i] = self.generator_list
					#del callWorkForce
					temp[i] = parallel_zetas
					
					Y = parallel_zetas
					T_ = np.linspace(300,2500,8)
					theta = np.array([T/T,np.log(T),-1/T])
					kappa_mean = theta.T.dot(P)
					Theta = np.array([T_/T_,np.log(T_),-1/T_])
					X = []
					for j in Y:
						#print(i)
						#print(kappa_max,np.asarray(theta.T.dot(j)).flatten())
						#print(np.asarray(theta.T.dot(j)).flatten()/kappa_max)
						#Pint = P + np.asarray(np.dot(cov,np.asarray(j))).flatten();
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					temp_xdata = {}
					for k in range(len(Y)):
						temp_xdata[str(Y[k])] = X[k]
						
					random.shuffle(Y)
					#print(Y_)
					X_ = []
					for k in range(len(Y)):
						X_.append(temp_xdata[str(Y[k])])
						
					#print(X_[9])
					xdata[i] = X_
					ydata[i] = Y
					
					#print(len(parallel_zetas))
				A = []
				A_gen = []
				X = []
				y = []
				for j in range(self.sim_):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp:
						a_row.extend(list(temp[i][j]))
						a_row_gen.extend(list(gen[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X.append(np.asarray(x_temp))
					y.append(np.asarray(y_temp))
				self.beta_ = np.asarray(A)
				self.generators = np.asarray(A_gen)
				self.X = np.asarray(X)
				self.y = np.asarray(y)
		elif self.simulation =='Opt':
			#Initializing the sampling size
			if "Sampling_of_PRS" not in self.opt_dict["Stats"]:
				self.sim_ = len(2*np.random.random(2*int(self.n_))-1)
			else:
				self.sim_ = len(2*np.random.random(int(0.7*self.opt_dict["Stats"]["Sampling_of_PRS"]*self.n_))-1)
			#print(self.sim_)
			#Initializing the design for the response surface
			#print(self.n)
			if self.PRS_type == "Partial":
				self.n = len(self.selectedParams)
			#print(self.n)
			if self.design == "Full-Factorial":
				self.beta_ = ff2n(self.n)
			elif self.design == "Plackett-Burman":
				self.beta_ = pbdesign(self.n)
			elif self.design == "Box-Benkhen":
				self.beta_ = bbdesign(self.n)
				self.sim_ = len(self.beta_)
			elif self.design == "center-composite-rotational":
				self.beta_ = ccdesign(self.n, center=(0, 1), alpha='r',face='cci')
			elif self.design == "center-composite-orthogonal-ccc":
				self.beta_ = ccdesign(self.n, center=(0, 1), alpha='o', face='ccc')
			elif self.design == "center-composite-orthogonal-cci":
				self.beta_ = ccdesign(self.n, center=(0, 1), alpha='o', face='cci')
			elif self.design == "center-composite-orthogonal-ccf":
				self.beta_ = ccdesign(self.n, center=(0, 1), alpha='o', face='ccf')
			elif self.design == 'LHS-maximin':  #Latin hypercube sampling (Maximin)
				self.beta_ = list(2*lhs(self.n,samples = self.sim_,criterion='centermaximin')-1)
				self.beta_.extend(list(np.eye(self.n)))
				self.beta_.extend(list(-1*np.eye(self.n)))
				self.beta_.extend(list(np.zeros(self.n)))
				self.sim_ = len(self.beta_)
			elif self.design == 'SAB': #Monte-carlo sampling
				self.beta_ = []
				self.beta_.extend(list(np.zeros(self.n)))
				self.beta_.extend(list(0.5*np.eye(self.n)))
				self.beta_.extend(list(-0.5*np.eye(self.n)))
				#mean = np.zeros(self.n)
				#cov = 0.33*np.eye(self.n)
				#self.beta_ = list(np.random.multivariate_normal(mean,cov,self.sim_))
				#self.beta_.extend(list(np.eye(self.n)))
				#self.beta_.extend(list(-1*np.eye(self.n)))
				#self.beta_.extend(list(np.zeros(self.n)))
				self.sim_ = len(self.beta_)
			elif self.design == 'Monte-Carlo': #Monte-carlo sampling
				#mean = np.zeros(self.n)
				#cov = 0.33*np.eye(self.n)
				#self.beta_ = list(np.random.multivariate_normal(mean,cov,self.sim_))
				self.beta_ = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				self.beta_.extend(list(np.eye(self.n)))
				self.beta_.extend(list(-1*np.eye(self.n)))
				self.beta_.extend(list(bbdesign(self.n)[0:int(0.3*self.sim_)]))
				self.beta_.extend(list(np.zeros(self.n)))
				self.beta_.extend(list(np.zeros(self.n)))
				self.beta_.extend(list(np.zeros(self.n)))
				self.sim_ = len(self.beta_)
			elif self.design == 'A-facto': 
				self.beta_=[]
				task = "complete"
				if "Data" in os.listdir():
					os.chdir("Data/Simulations")
					if f"Beta_list_case-0.csv" in os.listdir():
						#print(great)
						#print(os.getcwd())
						
						A_file = open("Beta_list_case-0.csv","r").readlines() 
						self.sim_ = len(A_file)
						print(self.sim_)
						A = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in A_file])	
						A_gen_file = open("Generator_list_case-0.csv","r").readlines() 
						A_gen = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in A_gen_file])	
						
						
						X_file = open("NormalizedBeta_list_case-0.csv","r").readlines() 
						X_dash = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in A_file])
						y_dash = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in A_file])
						self.beta_ = np.asarray(A)
						#print(self.beta_)
						self.generators = np.asarray(A_gen)
						
						self.X = np.asarray(X_dash)
						self.y = np.asarray(y_dash)
						os.chdir("../..")
					else:
						os.chdir("../..")
						task = "incomplete"
				if task != "complete":
					self.beta_ = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
					#self.beta_.extend(list(np.eye(self.n)))
					#self.beta_.extend(list(-1*np.eye(self.n)))
					#self.beta_.extend(list(bbdesign(self.n)[0:int(0.3*self.sim_)]))
					#self.beta_.extend(list(np.zeros(self.n)))
					#self.beta_.extend(list(np.zeros(self.n)))
					#self.beta_.extend(list(np.zeros(self.n)))
					self.sim_ = len(self.beta_)
					self.generators = np.asarray(self.beta_)
					self.X = np.asarray(self.beta_)
					self.y = np.asarray(self.beta_)
			
			elif self.design == 'Monte-Carlo-All': #Monte-carlo sampling
				#mean = np.zeros(self.n)
				#cov = 0.33*np.eye(self.n)
				#self.beta_ = list(np.random.multivariate_normal(mean,cov,self.sim_))
				self.beta_ = list(2*(1/3)*np.random.random_sample((self.sim_,self.n)) - (1/3))
				self.beta_.extend(list((1/3)*np.eye(self.n)))
				self.beta_.extend(list((-1/3)*np.eye(self.n)))
				#self.beta_.extend(list(bbdesign(self.n)[0:int(0.3*self.sim_)]))
				self.beta_.extend([np.zeros(self.n)])
				self.beta_.extend([np.zeros(self.n)])
				self.beta_.extend([np.zeros(self.n)])
				self.sim_ = len(self.beta_)
			elif self.design == "Mix-Monte-Carlo":
				self.beta_ = list(bbdesign(self.n))
				self.beta_.extend(list(np.eye(self.n)))
				self.beta_.extend(list(-1*np.eye(self.n)))
				mean = np.zeros(self.n)
				cov = 0.33*np.eye(self.n)
				g = list(np.random.multivariate_normal(mean,cov,1000))
				self.beta_.extend(g)
				self.sim_ = len(self.beta_)
			elif self.design == 'Legendre':
				self.beta_ = np.matrix(pd.read_csv(home_dir+'/outfile_ortho.csv'))
				self.beta_ = self.beta_.tolist()
				#print(self.beta_[0])
				print("Beta list (Legendre) found !! Shape of beta list {}".format(np.shape(beta_)))
			
			elif self.design == "B1-factorial":
				self.beta_=[]
				temp_a = {}
				temp_b = {}
				temp_c = {}
				gen = {}
				gen_a = {}
				gen_b = {}
				gen_c = {}
				
				xdata = {}
				ydata = {}				
				self.generator_list = []
				"""
				Min- simulation = 300
				"""
				#if self.sim_ < 300:
				#	self.sim_ = 300
				#print(self.sim_)
				"""
				Dividing the generators into three bins
				------------------------------------
				n_a = 0.33*self.sim_
				n_b = 0.33*self.sim_
				n_c = 0.34*self.sim_
				
				"""
				n_a = int(0.33*self.sim_)
				n_b = int(0.33*self.sim_)
				n_c = self.sim_-n_a-n_b
				
				
				#n_b = int(0.5*self.sim_)
				#n_c = self.sim_-n_b
				
				
				#for i in np.eye(self.n):
				#	self.generator_list.append(i)
				#	self.generator_list.append(-1*i)
								
				"""
				Mixing three types of the extreme curves
				--------------------
				zeta-A: 1
				zeta-B: 1
				zeta-C: 1
			
				for zeta-A: one normal random parameter is required
					a matrix on n_a x n will be created
					
				for zeta-B and zeta-c: two random parameters is required
					we can arrange the two normal_random_parameter 				
					a matrix of n_b+n_c x 2*n will be created
					
				The number of reactions = n_R
				"""
				n_rxn = len(self.ind)
				self.generator_list_A = []
				self.generator_list_B = []
				self.generator_list_C = []
				
				"""
				Factorial Design
				----------------
				box-benkhen design
				bb = []
				for i in bbdesign(n_rxn)[0:int(0.3*self.sim_)]:
					bb.append(np.asarray(i))
				"""
						#For zeta-A
						
				bb_a = []
				for i in bbdesign(n_rxn):
					bb_a.append(np.asarray(i))
				self.generator_list_A.extend(bb_a)
				self.generator_list_A.extend(list(np.eye(n_rxn)))
				self.generator_list_A.extend(list(-1*np.eye(n_rxn)))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))		
						#For zeta-B
						
				bb_b = []
				for i in bbdesign(2*n_rxn):
					bb_b.append(np.asarray(i))							
				self.generator_list_B.extend(bb_b)
				self.generator_list_B.extend(list(np.eye(2*n_rxn)))
				self.generator_list_B.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))		
						#For zeta-C
						
				bb_c = []
				for i in bbdesign(2*n_rxn):
					bb_c.append(np.asarray(i))		
						
				self.generator_list_C.extend(bb_c)		
				self.generator_list_C.extend(list(np.eye(2*n_rxn)))
				self.generator_list_C.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))		
				
				"""
				Monte-Carlo
				
				"""
				"""		
						#For zeta-A
				
				self.generator_list_A.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(0.3*n_a)])
				
						#For zeta-A
						
				self.generator_list_B.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_a)])	
						
						#For zeta-A		
				self.generator_list_C.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_a)])
				
				"""
				"""
				Uniform
				
				"""
					
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				
				"""
				Latin-Hypercube
				"""
				"""
						
						#For zeta-A
						
				lhs_a = []
				for i in lhs(n_rxn, samples=int(0.3*n_a)):
					lhs_a.append(np.asarray(i))	
				self.generator_list_A.extend(lhs_a)
				
						
						#For zeta-A
						
				lhs_b = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_b)):
					lhs_b.append(np.asarray(i))		
				self.generator_list_B.extend(lhs_b)
						
						#For zeta-A
						
				lhs_c = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_c)):
					lhs_c.append(np.asarray(i))
				self.generator_list_C.extend(lhs_c)
				
				"""
				"""
				Monte-Carlo-Mix
				"""
				
				
				
				
				
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)		
				#self.generator_list.extend(bb)
				#self.generator_list.extend(list(np.eye(self.n)))
				#self.generator_list.extend(list(-1*np.eye(self.n)))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#print(self.generator_list)
				#self.sim_ = len(self.generator_list)
				
				#trainers = np.asarray(trainers)
				#print(f"Generator_list {np.shape(self.generator_list)}\n")
				df = pd.DataFrame(np.asarray(self.generator_list))
				count_start_a = 0
				count_start_b = 0
				count_start_c = 0			
				for i in self.reactions_selected:
					
					len_of_Arrhenius_A = 1
					count_end_a = count_start_a+len_of_Arrhenius_A
					
					
					len_of_Arrhenius_B = 2
					count_end_b = count_start_b+len_of_Arrhenius_B
					
					
					len_of_Arrhenius_C = 2
					count_end_c = count_start_c+len_of_Arrhenius_C
					
					#print(count_start)
					#print(count_end)
					self.rxn_generators_a2 = []
					self.rxn_generators_b2 = []
					self.rxn_generators_c2 = []
					
					"""
					getting zetas from generators
					"""
					for j in range(len(self.generator_list_A)):
						self.rxn_generators_a2.append(self.generator_list_A[j][count_start_a:count_end_a])
					#print(self.rxn_generators_a2)
					for j in range(len(self.generator_list_B)):
						self.rxn_generators_b2.append(self.generator_list_B[j][count_start_b:count_end_b])
					
					for j in range(len(self.generator_list_C)):
						self.rxn_generators_c2.append(self.generator_list_C[j][count_start_c:count_end_c])
					
					self.sim_ = len(self.rxn_generators_a2)+len(self.rxn_generators_b2)+len(self.rxn_generators_c2)
					#print(self.rxn_generators)
					#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
										
					g_a = self.rxn_generators_a2
					g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					#g_b = self.rxn_generators[0:n_b]
					#g_c = self.rxn_generators[n_b:n_b+n_c]
					
					
					count_start_a = count_end_a
					count_start_b = count_end_b
					count_start_c = count_end_c
					
					"""
					Getting the uncertainty data
					"""
					#fig = plt.figure()
					data = self.unsrt[i].data
					T = np.linspace(300,2500,8)
					kappa_max = self.unsrt[i].getKappaMax(T)			
					kappa_min = self.unsrt[i].getKappaMin(T)
					Theta = np.array([T/T,np.log(T),-1/T])
					P = self.unsrt[i].getMean()	
					cov = self.unsrt[i].getCov()
					zeta_A = self.unsrt[i].zeta.x
					
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					"""
					Type-A zeta samples
					"""
					generators_a = g_a
					parallel_zetas_a = []
					
					time_start = time.time()
					
					
					for k in g_a:
						#Pint = P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten()
					#	plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
					#	plt.plot()
						#print(i)
						#print(type(i))
						#print(k)	
						#parallel_zetas_a.append(P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten())
						parallel_zetas_a.append(k[0]*zeta_A)
										
					#print(parallel_zetas_a)
					#plt.savefig(str(i)+"_a.pdf")
					#print(os.getcwd())
					#print("plot saved")
					data["generators_b"] = g_b
					callWorkForce = Worker(self.allowed_count)	
					generators_b,parallel_zetas_b = callWorkForce.do_unsrt_b(data,len(g_b))
					del callWorkForce
					
					data["generators_c"] = g_c
					callWorkForce = Worker(self.allowed_count)	
					generators_c,parallel_zetas_c = callWorkForce.do_unsrt_c(data,len(g_c))
					del callWorkForce
					
					#print(parallel_zetas_a)
					
					temp_a[i] = parallel_zetas_a
					gen_a[i] = generators_a
					temp_b[i] = parallel_zetas_b
					gen_b[i] = generators_b
					temp_c[i] = parallel_zetas_c
					gen_c[i] = generators_c
					
					time_end = time.time()
					print(f"\n\t\tTime taken = {time_end-time_start}\n")
					
					#print(parallel_zetas_b)
					Y_a = parallel_zetas_a
					Y_b = parallel_zetas_b
					Y_c = parallel_zetas_c
					T_ = np.linspace(300,2500,8)
					theta = np.array([T/T,np.log(T),-1/T])
					kappa_mean = theta.T.dot(P)
					Theta = np.array([T/T,np.log(T),-1/T])
					X = []
					Y = []
					
					
					for j in Y_a:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					
					for j in Y_b:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_b.pdf")
					
					
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					for j in Y_c:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_c.pdf")
					
					temp_xdata = {}
					for k in range(len(Y)):
						temp_xdata[str(Y[k])] = X[k]
						
					random.shuffle(Y)
					#print(Y_)
					X_ = []
					for k in range(len(Y)):
						X_.append(temp_xdata[str(Y[k])])
						
					#print(X_[9])
					xdata[i] = X_
					ydata[i] = Y
					#print(len(parallel_zetas))
					#print(len(self.rxn_generators))
				A = []
				A_gen = []
				X_dash = []
				y_dash = []
				
				for j in range(len(g_a)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_a:
						a_row.extend(list(temp_a[i][j]))
						a_row_gen.extend(list(gen_a[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				
				
				for j in range(len(g_b)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)]))
						y_temp.extend(list(ydata[i][j+len(g_a)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				for j in range(n_b):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				
				for j in range(len(g_c)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)+len(g_b)]))
						y_temp.extend(list(ydata[i][j+len(g_a)+len(g_b)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				for j in range(n_c):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+n_b]))
						y_temp.extend(list(ydata[i][j+n_b]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
					
				"""
				
				#print(A)
				self.sim_ = len(g_a)+len(g_b)+len(g_c)
				#self.sim_ = n_b+n_c
				self.beta_ = np.asarray(A)
				#print(self.beta_)
				self.generators = np.asarray(A_gen)
				self.X = np.asarray(X_dash)
				self.y = np.asarray(y_dash)
			
			elif self.design == "SAMAP":
				"""
				Getting the SAMAP amples
				"""
				os.mkdir("SAMAP_"+str(self.SAMAP_INDEX))
				os.chdir("SAMAP_"+str(self.SAMAP_INDEX))
				input_rxn_dict = {}
				samples_dict = {}
				for i in self.reactions_selected:
					os.mkdir(f"{self.unsrt[i].rIndex}")
					os.chdir(f"{self.unsrt[i].rIndex}")
					input_dict = self.unsrt[i].input_dict
					input_dict["samples"] = self.sim_
					input_dict["samples_skipped"] = int(0.1*self.sim_)
					input_dict["Random_seed"] = 1
					input_dict["sampling_method"] = "SOBOL"
					input_dict["sampling_distribution"] = "NORMAL"
					input_dict["equidistant_T"] = 100
					input_dict["T_begin"] = data = self.unsrt[i].temperatures[0]
					input_dict["T_end"] = self.unsrt[i].temperatures[-1]
					
					L = self.unsrt[i].cholskyDeCorrelateMat
					#input_dict["covariance_matrix"] = str(L.dot(L.T)).strip("[]").replace("[","").replace("]","")
					input_rxn_dict[i] = input_dict
					"""
					Run: SAMAP code
					"""
					instring = make_input_file.create_SAMAP_input(input_dict)
					#string_dict[i] = instring
					file_print = open("samap_data_"+str(self.unsrt[i].rIndex)+".txt","w").write(instring)
					run_string = f"""#!/bin/bash
		{self.samap_executable} samap_data_{self.unsrt[i].rIndex}.txt &> out"""
					file_print_run = open("run","w").write(run_string)						
					subprocess.call(["chmod","+x",'run'])
					subprocess.call(["./run"])
					file_name = "samap_data_"+str(self.unsrt[i].rIndex)+".txt_Arrhpar.txt"
					sampling_file = open(file_name,"r").readlines()
					self.Nagy_arrhenius_samples = [] 
					for sample in sampling_file[1:]:
						self.Nagy_arrhenius_samples.append(np.asfarray(sample.strip("''").strip("\n").split()[:3],float))
					
					#print(self.Nagy_arrhenius_samples)
					Y = self.Nagy_arrhenius_samples
					
					
					data = self.unsrt[i].data
					
					P = self.unsrt[i].getMean()	
					#print(P)
					cov = self.unsrt[i].getCov()
					T_ = np.linspace(300,2500,100)
					Theta = np.array([T_/T_,np.log(T_),-1/T_])
					"""
					fig = plt.figure()
					plt.plot(1/T_,self.unsrt[i].getKappaMax(T_),'k-',linewidth=0.75)
					plt.plot(1/T_,self.unsrt[i].getKappaMin(T_),'k-',linewidth=0.75)
					plt.plot(1/T_,self.unsrt[i].getNominal(T_),"b-",linewidth=0.75)
					"""
					temp = []
					for j in Y:
						#print(j)
						#print(kappa_max,np.asarray(theta.T.dot(j)).flatten())
						#print(np.asarray(theta.T.dot(j)).flatten()/kappa_max)
						#Pint = P + np.asarray(np.dot(cov,np.asarray(j))).flatten();
						Pint = np.asarray(j)
						ZETA = np.asarray(np.linalg.inv(cov).dot(Pint-P)).flatten()
						temp.append(ZETA)
							#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
						#X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#plt.plot(1/T_,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
					samples_dict[i] = temp
					#plt.show()
				
					os.chdir("..")
					
				A = []
				A_gen = []
				X = []
				y = []
				for j in range(self.sim_):
					a_row = []
					for i in samples_dict:
						a_row.extend(list(samples_dict[i][j]))
					
					A.append(np.asarray(a_row))
					
				self.beta_ = np.asarray(A)
				#print(self.beta_)
				self.generators = []
				self.X = []
				self.y = []
				
				self.SAMAP_INDEX+=1
				os.chdir("..")
			
			elif self.design == "A1+C1":
				self.beta_=[]
				temp_a = {}
				temp_b = {}
				temp_c = {}
				gen = {}
				gen_a = {}
				gen_b = {}
				gen_c = {}
				xdata = {}
				ydata = {}
				
				temp_a_ = {}
				temp_b_ = {}
				temp_c_ = {}
				gen_ = {}
				gen_a_ = {}
				gen_b_ = {}
				gen_c_ = {}
				
				xdata_ = {}
				ydata_ = {}				
				self.generator_list = []
				"""
				Min- simulation = 300
				"""
				#if self.sim_ < 300:
				#	self.sim_ = 300
				#print(self.sim_)
				"""
				Dividing the generators into three bins
				------------------------------------
				n_a = 0.33*self.sim_
				n_b = 0.33*self.sim_
				n_c = 0.34*self.sim_
				
				"""
				shuffled = int(0.8*self.sim_)
				unshuffled = self.sim_ - shuffled
				
				n_a = int(0.1*shuffled)
				n_b = int(0.45*shuffled)
				n_c = shuffled-n_a-n_b
				
				n_a_ = int(0.1*unshuffled)
				n_b_ = int(0.45*unshuffled)
				n_c_ = unshuffled-n_a_-n_b_
				
				#n_b = int(0.5*self.sim_)
				#n_c = self.sim_-n_b
				
				
				#for i in np.eye(self.n):
				#	self.generator_list.append(i)
				#	self.generator_list.append(-1*i)
								
				"""
				Mixing three types of the extreme curves
				--------------------
				zeta-A: 1
				zeta-B: 1
				zeta-C: 1
			
				for zeta-A: one normal random parameter is required
					a matrix on n_a x n will be created
					
				for zeta-B and zeta-c: two random parameters is required
					we can arrange the two normal_random_parameter 				
					a matrix of n_b+n_c x 2*n will be created
					
				The number of reactions = n_R
				"""
				
				
				n_rxn = len(self.ind)
				self.generator_list_A = []
				self.generator_list_B = []
				self.generator_list_C = []
				
				self.generator_list_A_ = []
				self.generator_list_B_ = []
				self.generator_list_C_ = []
				
				"""
				Factorial Design
				----------------
				box-benkhen design
				bb = []
				for i in bbdesign(n_rxn)[0:int(0.3*self.sim_)]:
					bb.append(np.asarray(i))
				"""
				"""		#For zeta-A
					
				bb_a = []
				for i in bbdesign(n_rxn)[0:int(0.4*n_a)]:
					bb_a.append(np.asarray(i))
				self.generator_list_A.extend(bb_a)
				
						#For zeta-B
						
				bb_b = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_b)]:
					bb_b.append(np.asarray(i))							
				self.generator_list_B.extend(bb_b)	
				
						#For zeta-C
						
				bb_c = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_c)]:
					bb_c.append(np.asarray(i))			
				self.generator_list_C.extend(bb_c)		
				
				"""		
					
				bb_a = []
				for i in bbdesign(n_rxn):
					bb_a.append(np.asarray(i))
				self.generator_list_A.extend(bb_a)
				
						#For zeta-B
						
				bb_b = []
				for i in bbdesign(2*n_rxn):
					bb_b.append(np.asarray(i))							
				self.generator_list_B.extend(bb_b)	
				
						#For zeta-C
						
				bb_c = []
				for i in bbdesign(2*n_rxn):
					bb_c.append(np.asarray(i))			
				self.generator_list_C.extend(bb_c)		
						
				
				"""
				Shuffled cases
				"""
					
				bb_a_ = []
				for i in bbdesign(n_rxn):
					bb_a_.append(np.asarray(i))
				self.generator_list_A_.extend(bb_a_)
				
						#For zeta-B
						
				bb_b_ = []
				for i in bbdesign(2*n_rxn):
					bb_b_.append(np.asarray(i))							
				self.generator_list_B_.extend(bb_b_)	
				
						#For zeta-C
						
				bb_c_ = []
				for i in bbdesign(2*n_rxn):
					bb_c_.append(np.asarray(i))			
				self.generator_list_C_.extend(bb_c_)		
				
				shuffled = len(self.generator_list_A)+len(self.generator_list_C)
				unshuffled = len(self.generator_list_A_)+len(self.generator_list_C_)
				
				
				"""
				Monte-Carlo
				
				"""
				"""		
						#For zeta-A
				
				self.generator_list_A.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(0.4*n_a)])
				self.generator_list_A.extend(list(np.eye(n_rxn)))
				self.generator_list_A.extend(list(-1*np.eye(n_rxn)))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				
						#For zeta-b
						
				self.generator_list_B.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.4*n_b)])	
				self.generator_list_B.extend(list(np.eye(2*n_rxn)))
				self.generator_list_B.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
						
						#For zeta-C		
				self.generator_list_C.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.4*n_c)])
				self.generator_list_C.extend(list(np.eye(2*n_rxn)))
				self.generator_list_C.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				
				
				self.generator_list_A_.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(n_a_)])
				#self.generator_list_A_.extend(list(np.eye(n_rxn)))
				#self.generator_list_A_.extend(list(-1*np.eye(n_rxn)))
				#self.generator_list_A_.append(np.zeros(n_rxn))
				#self.generator_list_A_.append(np.zeros(n_rxn))
				#self.generator_list_A_.append(np.zeros(n_rxn))
				
						#For zeta-A
						
				self.generator_list_B_.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(n_b_)])	
				#self.generator_list_B_.extend(list(np.eye(2*n_rxn)))
				#self.generator_list_B_.extend(list(-1*np.eye(2*n_rxn)))
				#self.generator_list_B_.append(np.zeros(2*n_rxn))
				#self.generator_list_B_.append(np.zeros(2*n_rxn))
				#self.generator_list_B_.append(np.zeros(2*n_rxn))
						
						#For zeta-A		
				self.generator_list_C_.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(n_c_)])
				#self.generator_list_C_.extend(list(np.eye(2*n_rxn)))
				#self.generator_list_C_.extend(list(-1*np.eye(2*n_rxn)))
				#self.generator_list_C_.append(np.zeros(2*n_rxn))
				#self.generator_list_C_.append(np.zeros(2*n_rxn))
				#self.generator_list_C_.append(np.zeros(2*n_rxn))
				
				
				"""
				"""
				For testing the distribution type
				"""
				"""		
						#For zeta-A
				
				self.generator_list_A.extend(list(2* np.random.random_sample((int(0.33*shuffled),n_rxn)) - 1))
				self.generator_list_A.extend(list(np.eye(n_rxn)))
				self.generator_list_A.extend(list(-1*np.eye(n_rxn)))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				
						#For zeta-b
						
				self.generator_list_B.extend(list(2* np.random.random_sample((int(0.33*shuffled),2*n_rxn)) - 1))	
				self.generator_list_B.extend(list(np.eye(2*n_rxn)))
				self.generator_list_B.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
						
						#For zeta-C		
				self.generator_list_C.extend(list(2* np.random.random_sample((int(0.33*shuffled),2*n_rxn)) - 1))
				self.generator_list_C.extend(list(np.eye(2*n_rxn)))
				self.generator_list_C.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				
				
				self.generator_list_A_.extend(list(2* np.random.random_sample((int(0.33*unshuffled),n_rxn)) - 1))
				#self.generator_list_A_.extend(list(np.eye(n_rxn)))
				#self.generator_list_A_.extend(list(-1*np.eye(n_rxn)))
				#self.generator_list_A_.append(np.zeros(n_rxn))
				#self.generator_list_A_.append(np.zeros(n_rxn))
				#self.generator_list_A_.append(np.zeros(n_rxn))
				
						#For zeta-A
						
				self.generator_list_B_.extend(list(2* np.random.random_sample((int(0.33*unshuffled),2*n_rxn)) - 1))	
				#self.generator_list_B_.extend(list(np.eye(2*n_rxn)))
				#self.generator_list_B_.extend(list(-1*np.eye(2*n_rxn)))
				#self.generator_list_B_.append(np.zeros(2*n_rxn))
				#self.generator_list_B_.append(np.zeros(2*n_rxn))
				#self.generator_list_B_.append(np.zeros(2*n_rxn))
						
						#For zeta-A		
				self.generator_list_C_.extend(list(2* np.random.random_sample((int(0.33*unshuffled),2*n_rxn)) - 1))
				#self.generator_list_C_.extend(list(np.eye(2*n_rxn)))
				#self.generator_list_C_.extend(list(-1*np.eye(2*n_rxn)))
				#self.generator_list_C_.append(np.zeros(2*n_rxn))
				#self.generator_list_C_.append(np.zeros(2*n_rxn))
				#self.generator_list_C_.append(np.zeros(2*n_rxn))
				"""
				"""
				Uniform
				
				"""
					
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				
				"""
				Latin-Hypercube
				"""
						#For zeta-A
				"""		
				lhs_a = []
				for i in lhs(n_rxn, samples=int(0.3*n_a)):
					lhs_a.append(np.asarray(i))	
				self.generator_list_A.extend(lhs_a)
				
						
						#For zeta-A
						
				lhs_b = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_b)):
					lhs_b.append(np.asarray(i))		
				self.generator_list_B.extend(lhs_b)
						
						#For zeta-A
						
				lhs_c = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_c)):
					lhs_c.append(np.asarray(i))
				self.generator_list_C.extend(lhs_c)
				
				"""
				"""		
				lhs_a = []
				for i in lhs(n_rxn, samples=int(0.33*self.sim_)):
					lhs_a.append(np.asarray(i))	
				self.generator_list_A.extend(lhs_a)
				
						
						#For zeta-A
						
				lhs_b = []
				for i in lhs(2*n_rxn, samples=int(0.33*self.sim_)):
					lhs_b.append(np.asarray(i))		
				self.generator_list_B.extend(lhs_b)
						
						#For zeta-A
						
				lhs_c = []
				for i in lhs(2*n_rxn, samples=int(0.33*self.sim_)):
					lhs_c.append(np.asarray(i))
				self.generator_list_C.extend(lhs_c)
				
				"""
				"""
				lhs_a_ = []
				for i in lhs(n_rxn, samples=int(0.33*unshuffled)):
					lhs_a_.append(np.asarray(i))	
				self.generator_list_A_.extend(lhs_a_)
				
						
						#For zeta-A
						
				lhs_b_ = []
				for i in lhs(2*n_rxn, samples=int(0.33*unshuffled)):
					lhs_b_.append(np.asarray(i))		
				self.generator_list_B_.extend(lhs_b_)
						
						#For zeta-A
						
				lhs_c_ = []
				for i in lhs(2*n_rxn, samples=int(0.33*unshuffled)):
					lhs_c_.append(np.asarray(i))
				self.generator_list_C_.extend(lhs_c_)
				
				"""
				"""
				Monte-Carlo-Mix
				"""
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)		
				#self.generator_list.extend(bb)
				#self.generator_list.extend(list(np.eye(self.n)))
				#self.generator_list.extend(list(-1*np.eye(self.n)))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#print(self.generator_list)
				#self.sim_ = len(self.generator_list)
				
				#trainers = np.asarray(trainers)
				#print(f"Generator_list {np.shape(self.generator_list)}\n")
				df = pd.DataFrame(np.asarray(self.generator_list))
				count_start_a = 0
				count_start_b = 0
				count_start_c = 0			
				for i in self.reactions_selected:
					
					len_of_Arrhenius_A = 1
					count_end_a = count_start_a+len_of_Arrhenius_A
					
					
					len_of_Arrhenius_B = 2
					count_end_b = count_start_b+len_of_Arrhenius_B
					
					
					len_of_Arrhenius_C = 2
					count_end_c = count_start_c+len_of_Arrhenius_C
					
					#print(count_start)
					#print(count_end)
					self.rxn_generators_a2 = []
					self.rxn_generators_b2 = []
					self.rxn_generators_c2 = []
					
					self.rxn_generators_a2_ = []
					self.rxn_generators_b2_ = []
					self.rxn_generators_c2_ = []
					
					"""
					getting zetas from generators
					"""
					for j in range(len(self.generator_list_A)):
						self.rxn_generators_a2.append(self.generator_list_A[j][count_start_a:count_end_a])
					#print(self.rxn_generators_a2)
					for j in range(len(self.generator_list_B)):
						self.rxn_generators_b2.append(self.generator_list_B[j][count_start_b:count_end_b])
					
					for j in range(len(self.generator_list_C)):
						self.rxn_generators_c2.append(self.generator_list_C[j][count_start_c:count_end_c])
					
					
					for j in range(len(self.generator_list_A_)):
						self.rxn_generators_a2_.append(self.generator_list_A_[j][count_start_a:count_end_a])
					#print(self.rxn_generators_a2)
					for j in range(len(self.generator_list_B_)):
						self.rxn_generators_b2_.append(self.generator_list_B_[j][count_start_b:count_end_b])
					
					for j in range(len(self.generator_list_C_)):
						self.rxn_generators_c2_.append(self.generator_list_C_[j][count_start_c:count_end_c])
					
					
					
					
					self.sim_ = len(self.rxn_generators_a2)+len(self.rxn_generators_c2)+len(self.rxn_generators_a2_)+len(self.rxn_generators_c2_)
					#print(self.rxn_generators)
					#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
										
					g_a = self.rxn_generators_a2
					g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					
					g_a_ = self.rxn_generators_a2_
					g_b_ = self.rxn_generators_b2_
					g_c_ = self.rxn_generators_c2_
					
					#g_b = self.rxn_generators[0:n_b]
					#g_c = self.rxn_generators[n_b:n_b+n_c]
					
					
					count_start_a = count_end_a
					count_start_b = count_end_b
					count_start_c = count_end_c
					
					"""
					Getting the uncertainty data
					"""
					#fig = plt.figure()
					data = self.unsrt[i].data
					T = np.linspace(300,2500,8)
					kappa_max = self.unsrt[i].getKappaMax(T)			
					kappa_min = self.unsrt[i].getKappaMin(T)
					Theta = np.array([T/T,np.log(T),-1/T])
					P = self.unsrt[i].getMean()	
					cov = self.unsrt[i].getCov()
					zeta_A = self.unsrt[i].zeta.x
					
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					"""
					Type-A zeta samples
					"""
					generators_a = g_a
					generators_a_ = g_a_
					parallel_zetas_a = []
					parallel_zetas_a_ = []
					
					time_start = time.time()
					
					
					for k in g_a:
						
						parallel_zetas_a.append(k[0]*zeta_A)
					
					for k in g_a_:
						
						parallel_zetas_a_.append(k[0]*zeta_A)
										
					
					generators_b_ = []
					parallel_zetas_b_ = []
					generators_c_ = []
					parallel_zetas_c_ = []
					
					generators_b = []
					parallel_zetas_b = []
					generators_c = []
					parallel_zetas_c = []
					start = 0
					start_ = 0
					for w in range(8):
						stop = int(((w+1)/8)*shuffled)
						stop_ = int(((w+1)/8)*unshuffled)
						#g_b = self.rxn_generators_b2[start:stop]
						g_c = self.rxn_generators_c2[start:stop]
						#g_b_ = self.rxn_generators_b2_[start_:stop_]
						g_c_ = self.rxn_generators_c2_[start_:stop_]
						
						
						data["generators_c"] = g_c
						callWorkForce = Worker(self.allowed_count)	
						gc,pc = callWorkForce.do_unsrt_c(data,len(g_c))
						del callWorkForce
						
						generators_c.extend(gc)
						parallel_zetas_c.extend(pc)
						removed = data.pop("generators_c")
						
						data["generators_c"] = g_c_
						callWorkForce = Worker(self.allowed_count)	
						gc,pc = callWorkForce.do_unsrt_c(data,len(g_c_))
						del callWorkForce
						
						generators_c_.extend(gc)
						parallel_zetas_c_.extend(pc)
						removed = data.pop("generators_c")
						
						start = stop
						start_ = stop_
					
					#print(self.sim_)
					#print(len(parallel_zetas_b)+len(parallel_zetas_c))
					#g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					
					#g_b_ = self.rxn_generators_b2_
					g_c_ = self.rxn_generators_c2_
					
					
					temp_a[i] = parallel_zetas_a
					gen_a[i] = generators_a
					#temp_b[i] = parallel_zetas_b
					#gen_b[i] = generators_b
					temp_c[i] = parallel_zetas_c
					gen_c[i] = generators_c
					
					temp_a_[i] = parallel_zetas_a_
					gen_a_[i] = generators_a_
					#temp_b_[i] = parallel_zetas_b_
					#gen_b_[i] = generators_b_
					temp_c_[i] = parallel_zetas_c_
					gen_c_[i] = generators_c_
					
					time_end = time.time()
					print(f"\n\t\tTime taken = {time_end-time_start}\n")
					
					#print(parallel_zetas_b)
					Y_a = parallel_zetas_a
					#Y_b = parallel_zetas_b
					Y_c = parallel_zetas_c
					
					Y_a_ = parallel_zetas_a_
					#Y_b_ = parallel_zetas_b_
					Y_c_ = parallel_zetas_c_
					
					
					
					T_ = np.linspace(300,2500,8)
					theta = np.array([T/T,np.log(T),-1/T])
					kappa_mean = theta.T.dot(P)
					Theta = np.array([T/T,np.log(T),-1/T])
					X = []
					Y = []
					
					X_ = []
					Y_ = []
					
					
					for j in Y_a:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						
					for j in Y_a_:
						Y_.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						
						X_.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						
					
					
					for j in Y_c:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_c.pdf")
					
					for j in Y_c_:
						Y_.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X_.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_c.pdf")
					
					temp_xdata = {}
					for k in range(len(Y_)):
						temp_xdata[str(Y_[k])] = X_[k]
						
					random.shuffle(Y_)
					#print(Y_)
					X1_ = []
					for k in range(len(Y_)):
						X1_.append(temp_xdata[str(Y_[k])])
						
					#print(X_[9])
					xdata[i] = X
					ydata[i] = Y
					
					xdata_[i] = X1_
					ydata_[i] = Y_
					
					#print(len(X1_))
					#print(len(Y_))
					#print(len(X))
					#print(len(Y))
					print(self.sim_)
					#Was trying to shuffle the zetas but the try was not successfull
					
				A = []
				A_gen = []
				X_dash = []
				y_dash = []
				
				for j in range(len(g_a)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_a:
						a_row.extend(list(temp_a[i][j]))
						a_row_gen.extend(list(gen_a[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
					
				for j in range(len(g_a_)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_a_:
						a_row.extend(list(ydata_[i][j]))
						a_row_gen.extend(list(gen_a_[i][j]))
						x_temp.extend(list(xdata_[i][j]))
						y_temp.extend(list(ydata_[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				
				for j in range(len(g_c)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)]))
						y_temp.extend(list(ydata[i][j+len(g_a)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				
				#print(len())
				for j in range(len(g_c_)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c_:
						#print(len(ydata_[i]))
						#print(len(g_a_)+len(g_b_)+len(g_c_))
						a_row.extend(list(ydata_[i][j+len(g_a_)]))
						a_row_gen.extend(list(gen_c_[i][j]))
						x_temp.extend(list(xdata_[i][j+len(g_a_)]))
						y_temp.extend(list(ydata_[i][j+len(g_a_)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				
				
			
				self.beta_ = np.asarray(A)
				#print(self.beta_)
				self.generators = np.asarray(A_gen)
				self.X = np.asarray(X_dash)
				self.y = np.asarray(y_dash)
				
			elif self.design == "B1":
				self.beta_=[]
				temp_a = {}
				temp_b = {}
				temp_c = {}
				gen = {}
				gen_a = {}
				gen_b = {}
				gen_c = {}
				
				xdata = {}
				ydata = {}				
				self.generator_list = []
				"""
				Min- simulation = 300
				"""
				#if self.sim_ < 300:
				#	self.sim_ = 300
				#print(self.sim_)
				"""
				Dividing the generators into three bins
				------------------------------------
				n_a = 0.33*self.sim_
				n_b = 0.33*self.sim_
				n_c = 0.34*self.sim_
				
				"""
				n_a = int(0.33*self.sim_)
				n_b = int(0.33*self.sim_)
				n_c = self.sim_-n_a-n_b
				
				
				#n_b = int(0.5*self.sim_)
				#n_c = self.sim_-n_b
				
				
				#for i in np.eye(self.n):
				#	self.generator_list.append(i)
				#	self.generator_list.append(-1*i)
								
				"""
				Mixing three types of the extreme curves
				--------------------
				zeta-A: 1
				zeta-B: 1
				zeta-C: 1
			
				for zeta-A: one normal random parameter is required
					a matrix on n_a x n will be created
					
				for zeta-B and zeta-c: two random parameters is required
					we can arrange the two normal_random_parameter 				
					a matrix of n_b+n_c x 2*n will be created
					
				The number of reactions = n_R
				"""
				n_rxn = len(self.ind)
				self.generator_list_A = []
				self.generator_list_B = []
				self.generator_list_C = []
				
				"""
				Factorial Design
				----------------
				box-benkhen design
				bb = []
				for i in bbdesign(n_rxn)[0:int(0.3*self.sim_)]:
					bb.append(np.asarray(i))
				"""
				"""
						
						#For zeta-A
						
				bb_a = []
				for i in bbdesign(n_rxn):
					bb_a.append(np.asarray(i))
				self.generator_list_A.extend(bb_a)
						
						#For zeta-B
						
				bb_b = []
				for i in bbdesign(2*n_rxn):
					bb_b.append(np.asarray(i))							
				self.generator_list_B.extend(bb_b)
						
						#For zeta-C
						
				bb_c = []
				for i in bbdesign(2*n_rxn):
					bb_c.append(np.asarray(i))		
						
				self.generator_list_C.extend(bb_c)		
						
				"""
				"""
				Monte-Carlo
				
				"""
				"""		
						#For zeta-A
						
				
				self.generator_list_A = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				
				
						
						#For zeta-A
						
				self.generator_list_B = list(2* np.random.random_sample((self.sim_,self.n)) - 1)		
				
						
						#For zeta-A
						
				self.generator_list_C = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				"""
				
				"""
				Uniform
				
				"""
					
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				
				"""
				Latin-Hypercube
				"""
			
						
						#For zeta-A
						
				lhs_a = []
				for i in lhs(n_rxn, samples=n_a):
					lhs_a.append(np.asarray(i))	
				self.generator_list_A.extend(lhs_a)
				
						
						#For zeta-A
						
				lhs_b = []
				for i in lhs(2*n_rxn, samples=n_b):
					lhs_b.append(np.asarray(i))		
				self.generator_list_B.extend(lhs_b)
						
						#For zeta-A
						
				lhs_c = []
				for i in lhs(2*n_rxn, samples=n_c):
					lhs_c.append(np.asarray(i))
				self.generator_list_C.extend(lhs_c)
				
				"""
				Monte-Carlo-Mix
				"""
				
				
				
				
				
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)		
				#self.generator_list.extend(bb)
				#self.generator_list.extend(list(np.eye(self.n)))
				#self.generator_list.extend(list(-1*np.eye(self.n)))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#print(self.generator_list)
				#self.sim_ = len(self.generator_list)
				
				#trainers = np.asarray(trainers)
				#print(f"Generator_list {np.shape(self.generator_list)}\n")
				df = pd.DataFrame(np.asarray(self.generator_list))
				count_start_a = 0
				count_start_b = 0
				count_start_c = 0			
				for i in self.reactions_selected:
					
					len_of_Arrhenius_A = 1
					count_end_a = count_start_a+len_of_Arrhenius_A
					
					
					len_of_Arrhenius_B = 2
					count_end_b = count_start_b+len_of_Arrhenius_B
					
					
					len_of_Arrhenius_C = 2
					count_end_c = count_start_c+len_of_Arrhenius_C
					
					#print(count_start)
					#print(count_end)
					self.rxn_generators_a2 = []
					self.rxn_generators_b2 = []
					self.rxn_generators_c2 = []
					
					"""
					getting zetas from generators
					"""
					for j in range(len(self.generator_list_A)):
						self.rxn_generators_a2.append(self.generator_list_A[j][count_start_a:count_end_a])
					#print(self.rxn_generators_a2)
					for j in range(len(self.generator_list_B)):
						self.rxn_generators_b2.append(self.generator_list_B[j][count_start_b:count_end_b])
					
					for j in range(len(self.generator_list_C)):
						self.rxn_generators_c2.append(self.generator_list_C[j][count_start_c:count_end_c])
					
					self.sim_ = len(self.rxn_generators_a2)+len(self.rxn_generators_b2)+len(self.rxn_generators_c2)
					#print(self.rxn_generators)
					#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
										
					g_a = self.rxn_generators_a2
					g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					#g_b = self.rxn_generators[0:n_b]
					#g_c = self.rxn_generators[n_b:n_b+n_c]
					
					
					count_start_a = count_end_a
					count_start_b = count_end_b
					count_start_c = count_end_c
					
					"""
					Getting the uncertainty data
					"""
					#fig = plt.figure()
					data = self.unsrt[i].data
					T = np.linspace(300,2500,8)
					kappa_max = self.unsrt[i].getKappaMax(T)			
					kappa_min = self.unsrt[i].getKappaMin(T)
					Theta = np.array([T/T,np.log(T),-1/T])
					P = self.unsrt[i].getMean()	
					cov = self.unsrt[i].getCov()
					zeta_A = self.unsrt[i].zeta.x
					
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					"""
					Type-A zeta samples
					"""
					generators_a = g_a
					parallel_zetas_a = []
					
					time_start = time.time()
					
					
					for k in g_a:
						#Pint = P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten()
					#	plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
					#	plt.plot()
						#print(i)
						#print(type(i))
						#print(k)	
						#parallel_zetas_a.append(P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten())
						parallel_zetas_a.append(k[0]*zeta_A)
										
					#print(parallel_zetas_a)
					#plt.savefig(str(i)+"_a.pdf")
					#print(os.getcwd())
					#print("plot saved")
					data["generators_b"] = g_b
					callWorkForce = Worker(self.allowed_count)	
					generators_b,parallel_zetas_b = callWorkForce.do_unsrt_b(data,len(g_b))
					del callWorkForce
					
					data["generators_c"] = g_c
					callWorkForce = Worker(self.allowed_count)	
					generators_c,parallel_zetas_c = callWorkForce.do_unsrt_c(data,len(g_c))
					del callWorkForce
					
					#print(parallel_zetas_a)
					
					temp_a[i] = parallel_zetas_a
					gen_a[i] = generators_a
					temp_b[i] = parallel_zetas_b
					gen_b[i] = generators_b
					temp_c[i] = parallel_zetas_c
					gen_c[i] = generators_c
					
					time_end = time.time()
					print(f"\n\t\tTime taken = {time_end-time_start}\n")
					
					#print(parallel_zetas_b)
					Y_a = parallel_zetas_a
					Y_b = parallel_zetas_b
					Y_c = parallel_zetas_c
					T_ = np.linspace(300,2500,8)
					theta = np.array([T/T,np.log(T),-1/T])
					kappa_mean = theta.T.dot(P)
					Theta = np.array([T/T,np.log(T),-1/T])
					X = []
					Y = []
					
					
					for j in Y_a:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					
					for j in Y_b:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_b.pdf")
					
					
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					for j in Y_c:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_c.pdf")
					
					temp_xdata = {}
					for k in range(len(Y)):
						temp_xdata[str(Y[k])] = X[k]
						
					random.shuffle(Y)
					#print(Y_)
					X_ = []
					for k in range(len(Y)):
						X_.append(temp_xdata[str(Y[k])])
						
					#print(X_[9])
					xdata[i] = X_
					ydata[i] = Y
					#print(len(parallel_zetas))
					#print(len(self.rxn_generators))
				A = []
				A_gen = []
				X_dash = []
				y_dash = []
				
				for j in range(len(g_a)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_a:
						a_row.extend(list(temp_a[i][j]))
						a_row_gen.extend(list(gen_a[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				
				
				for j in range(len(g_b)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)]))
						y_temp.extend(list(ydata[i][j+len(g_a)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				for j in range(n_b):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				
				for j in range(len(g_c)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)+len(g_b)]))
						y_temp.extend(list(ydata[i][j+len(g_a)+len(g_b)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				for j in range(n_c):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+n_b]))
						y_temp.extend(list(ydata[i][j+n_b]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
					
				"""
				
				#print(A)
				self.sim_ = len(g_a)+len(g_b)+len(g_c)
				#self.sim_ = n_b+n_c
				self.beta_ = np.asarray(A)
				#print(self.beta_)
				self.generators = np.asarray(A_gen)
				self.X = np.asarray(X_dash)
				self.y = np.asarray(y_dash)
			
			elif self.design == "B1+C1":
				self.beta_=[]
				temp_a = {}
				temp_b = {}
				temp_c = {}
				gen = {}
				gen_a = {}
				gen_b = {}
				gen_c = {}
				xdata = {}
				ydata = {}
				
				temp_a_ = {}
				temp_b_ = {}
				temp_c_ = {}
				gen_ = {}
				gen_a_ = {}
				gen_b_ = {}
				gen_c_ = {}
				
				xdata_ = {}
				ydata_ = {}				
				self.generator_list = []
				"""
				Min- simulation = 300
				"""
				#if self.sim_ < 300:
				#	self.sim_ = 300
				#print(self.sim_)
				"""
				Dividing the generators into three bins
				------------------------------------
				n_a = 0.33*self.sim_
				n_b = 0.33*self.sim_
				n_c = 0.34*self.sim_
				
				"""
				shuffled = int(0.8*self.sim_)
				unshuffled = self.sim_ - shuffled
				
				n_a = int(0.1*shuffled)
				n_b = int(0.45*shuffled)
				n_c = shuffled-n_a-n_b
				
				n_a_ = int(0.1*unshuffled)
				n_b_ = int(0.45*unshuffled)
				n_c_ = unshuffled-n_a_-n_b_
				
				#n_b = int(0.5*self.sim_)
				#n_c = self.sim_-n_b
				
				
				#for i in np.eye(self.n):
				#	self.generator_list.append(i)
				#	self.generator_list.append(-1*i)
								
				"""
				Mixing three types of the extreme curves
				--------------------
				zeta-A: 1
				zeta-B: 1
				zeta-C: 1
			
				for zeta-A: one normal random parameter is required
					a matrix on n_a x n will be created
					
				for zeta-B and zeta-c: two random parameters is required
					we can arrange the two normal_random_parameter 				
					a matrix of n_b+n_c x 2*n will be created
					
				The number of reactions = n_R
				"""
				
				
				n_rxn = len(self.ind)
				self.generator_list_A = []
				self.generator_list_B = []
				self.generator_list_C = []
				
				self.generator_list_A_ = []
				self.generator_list_B_ = []
				self.generator_list_C_ = []
				
				"""
				Factorial Design
				----------------
				box-benkhen design
				bb = []
				for i in bbdesign(n_rxn)[0:int(0.3*self.sim_)]:
					bb.append(np.asarray(i))
				"""
						#For zeta-A
					
				bb_a = []
				for i in bbdesign(n_rxn)[0:int(0.4*n_a)]:
					bb_a.append(np.asarray(i))
				self.generator_list_A.extend(bb_a)
				
						#For zeta-B
						
				bb_b = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_b)]:
					bb_b.append(np.asarray(i))							
				self.generator_list_B.extend(bb_b)	
				
						#For zeta-C
						
				bb_c = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_c)]:
					bb_c.append(np.asarray(i))			
				self.generator_list_C.extend(bb_c)		
				
				"""		
					
				bb_a = []
				for i in bbdesign(n_rxn):
					bb_a.append(np.asarray(i))
				self.generator_list_A.extend(bb_a)
				
						#For zeta-B
						
				bb_b = []
				for i in bbdesign(2*n_rxn):
					bb_b.append(np.asarray(i))							
				self.generator_list_B.extend(bb_b)	
				
						#For zeta-C
						
				bb_c = []
				for i in bbdesign(2*n_rxn):
					bb_c.append(np.asarray(i))			
				self.generator_list_C.extend(bb_c)		
				"""		
				
				"""
				Shuffled cases
				"""
				"""	
				bb_a_ = []
				for i in bbdesign(n_rxn):
					bb_a_.append(np.asarray(i))
				self.generator_list_A_.extend(bb_a_)
				
						#For zeta-B
						
				bb_b_ = []
				for i in bbdesign(2*n_rxn):
					bb_b_.append(np.asarray(i))							
				self.generator_list_B_.extend(bb_b_)	
				
						#For zeta-C
						
				bb_c_ = []
				for i in bbdesign(2*n_rxn):
					bb_c_.append(np.asarray(i))			
				self.generator_list_C_.extend(bb_c_)		
				
				shuffled = len(self.generator_list_A)+len(self.generator_list_B)+len(self.generator_list_C)
				unshuffled = len(self.generator_list_A_)+len(self.generator_list_B_)+len(self.generator_list_C_)
				
				"""
				"""
				Monte-Carlo
				
				"""
						
						#For zeta-A
				
				self.generator_list_A.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(0.4*n_a)])
				self.generator_list_A.extend(list(np.eye(n_rxn)))
				self.generator_list_A.extend(list(-1*np.eye(n_rxn)))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				
						#For zeta-b
						
				self.generator_list_B.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.4*n_b)])	
				self.generator_list_B.extend(list(np.eye(2*n_rxn)))
				self.generator_list_B.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
						
						#For zeta-C		
				self.generator_list_C.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.4*n_c)])
				self.generator_list_C.extend(list(np.eye(2*n_rxn)))
				self.generator_list_C.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				
				
				self.generator_list_A_.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(n_a_)])
				#self.generator_list_A_.extend(list(np.eye(n_rxn)))
				#self.generator_list_A_.extend(list(-1*np.eye(n_rxn)))
				#self.generator_list_A_.append(np.zeros(n_rxn))
				#self.generator_list_A_.append(np.zeros(n_rxn))
				#self.generator_list_A_.append(np.zeros(n_rxn))
				
						#For zeta-A
						
				self.generator_list_B_.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(n_b_)])	
				#self.generator_list_B_.extend(list(np.eye(2*n_rxn)))
				#self.generator_list_B_.extend(list(-1*np.eye(2*n_rxn)))
				#self.generator_list_B_.append(np.zeros(2*n_rxn))
				#self.generator_list_B_.append(np.zeros(2*n_rxn))
				#self.generator_list_B_.append(np.zeros(2*n_rxn))
						
						#For zeta-A		
				self.generator_list_C_.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(n_c_)])
				#self.generator_list_C_.extend(list(np.eye(2*n_rxn)))
				#self.generator_list_C_.extend(list(-1*np.eye(2*n_rxn)))
				#self.generator_list_C_.append(np.zeros(2*n_rxn))
				#self.generator_list_C_.append(np.zeros(2*n_rxn))
				#self.generator_list_C_.append(np.zeros(2*n_rxn))
				
				
				
				"""
				For testing the distribution type
				"""
				"""		
						#For zeta-A
				
				self.generator_list_A.extend(list(2* np.random.random_sample((int(0.33*shuffled),n_rxn)) - 1))
				self.generator_list_A.extend(list(np.eye(n_rxn)))
				self.generator_list_A.extend(list(-1*np.eye(n_rxn)))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				
						#For zeta-b
						
				self.generator_list_B.extend(list(2* np.random.random_sample((int(0.33*shuffled),2*n_rxn)) - 1))	
				self.generator_list_B.extend(list(np.eye(2*n_rxn)))
				self.generator_list_B.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
						
						#For zeta-C		
				self.generator_list_C.extend(list(2* np.random.random_sample((int(0.33*shuffled),2*n_rxn)) - 1))
				self.generator_list_C.extend(list(np.eye(2*n_rxn)))
				self.generator_list_C.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				
				
				self.generator_list_A_.extend(list(2* np.random.random_sample((int(0.33*unshuffled),n_rxn)) - 1))
				#self.generator_list_A_.extend(list(np.eye(n_rxn)))
				#self.generator_list_A_.extend(list(-1*np.eye(n_rxn)))
				#self.generator_list_A_.append(np.zeros(n_rxn))
				#self.generator_list_A_.append(np.zeros(n_rxn))
				#self.generator_list_A_.append(np.zeros(n_rxn))
				
						#For zeta-A
						
				self.generator_list_B_.extend(list(2* np.random.random_sample((int(0.33*unshuffled),2*n_rxn)) - 1))	
				#self.generator_list_B_.extend(list(np.eye(2*n_rxn)))
				#self.generator_list_B_.extend(list(-1*np.eye(2*n_rxn)))
				#self.generator_list_B_.append(np.zeros(2*n_rxn))
				#self.generator_list_B_.append(np.zeros(2*n_rxn))
				#self.generator_list_B_.append(np.zeros(2*n_rxn))
						
						#For zeta-A		
				self.generator_list_C_.extend(list(2* np.random.random_sample((int(0.33*unshuffled),2*n_rxn)) - 1))
				#self.generator_list_C_.extend(list(np.eye(2*n_rxn)))
				#self.generator_list_C_.extend(list(-1*np.eye(2*n_rxn)))
				#self.generator_list_C_.append(np.zeros(2*n_rxn))
				#self.generator_list_C_.append(np.zeros(2*n_rxn))
				#self.generator_list_C_.append(np.zeros(2*n_rxn))
				"""
				"""
				Uniform
				
				"""
					
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				
				"""
				Latin-Hypercube
				"""
						#For zeta-A
						
				lhs_a = []
				for i in lhs(n_rxn, samples=int(0.3*n_a)):
					lhs_a.append(np.asarray(i))	
				self.generator_list_A.extend(lhs_a)
				
						
						#For zeta-A
						
				lhs_b = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_b)):
					lhs_b.append(np.asarray(i))		
				self.generator_list_B.extend(lhs_b)
						
						#For zeta-A
						
				lhs_c = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_c)):
					lhs_c.append(np.asarray(i))
				self.generator_list_C.extend(lhs_c)
				
				
				"""		
				lhs_a = []
				for i in lhs(n_rxn, samples=int(0.33*self.sim_)):
					lhs_a.append(np.asarray(i))	
				self.generator_list_A.extend(lhs_a)
				
						
						#For zeta-A
						
				lhs_b = []
				for i in lhs(2*n_rxn, samples=int(0.33*self.sim_)):
					lhs_b.append(np.asarray(i))		
				self.generator_list_B.extend(lhs_b)
						
						#For zeta-A
						
				lhs_c = []
				for i in lhs(2*n_rxn, samples=int(0.33*self.sim_)):
					lhs_c.append(np.asarray(i))
				self.generator_list_C.extend(lhs_c)
				
				"""
				"""
				lhs_a_ = []
				for i in lhs(n_rxn, samples=int(0.33*unshuffled)):
					lhs_a_.append(np.asarray(i))	
				self.generator_list_A_.extend(lhs_a_)
				
						
						#For zeta-A
						
				lhs_b_ = []
				for i in lhs(2*n_rxn, samples=int(0.33*unshuffled)):
					lhs_b_.append(np.asarray(i))		
				self.generator_list_B_.extend(lhs_b_)
						
						#For zeta-A
						
				lhs_c_ = []
				for i in lhs(2*n_rxn, samples=int(0.33*unshuffled)):
					lhs_c_.append(np.asarray(i))
				self.generator_list_C_.extend(lhs_c_)
				
				"""
				"""
				Monte-Carlo-Mix
				"""
				df = pd.DataFrame(np.asarray(self.generator_list))
				count_start_a = 0
				count_start_b = 0
				count_start_c = 0			
				for i in self.reactions_selected:
					
					#len_of_Arrhenius_A = 1
					#count_end_a = count_start_a+len_of_Arrhenius_A
					
					
					len_of_Arrhenius_B = 2
					count_end_b = count_start_b+len_of_Arrhenius_B
					
					
					len_of_Arrhenius_C = 2
					count_end_c = count_start_c+len_of_Arrhenius_C
					
					#print(count_start)
					#print(count_end)
					#self.rxn_generators_a2 = []
					self.rxn_generators_b2 = []
					self.rxn_generators_c2 = []
					
					#self.rxn_generators_a2_ = []
					self.rxn_generators_b2_ = []
					self.rxn_generators_c2_ = []
					
					"""
					getting zetas from generators
					"""
					#for j in range(len(self.generator_list_A)):
					#	self.rxn_generators_a2.append(self.generator_list_A[j][count_start_a:count_end_a])
					#print(self.rxn_generators_a2)
					for j in range(len(self.generator_list_B)):
						self.rxn_generators_b2.append(self.generator_list_B[j][count_start_b:count_end_b])
					
					for j in range(len(self.generator_list_C)):
						self.rxn_generators_c2.append(self.generator_list_C[j][count_start_c:count_end_c])
					
					
					#for j in range(len(self.generator_list_A_)):
					#	self.rxn_generators_a2_.append(self.generator_list_A_[j][count_start_a:count_end_a])
					#print(self.rxn_generators_a2)
					for j in range(len(self.generator_list_B_)):
						self.rxn_generators_b2_.append(self.generator_list_B_[j][count_start_b:count_end_b])
					
					for j in range(len(self.generator_list_C_)):
						self.rxn_generators_c2_.append(self.generator_list_C_[j][count_start_c:count_end_c])
					
					
					
					
					self.sim_ =len(self.rxn_generators_b2)+len(self.rxn_generators_c2)+len(self.rxn_generators_b2_)+len(self.rxn_generators_c2_)
					#print(self.rxn_generators)
					#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
										
					#g_a = self.rxn_generators_a2
					g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					
					#g_a_ = self.rxn_generators_a2_
					g_b_ = self.rxn_generators_b2_
					g_c_ = self.rxn_generators_c2_
					#g_b = self.rxn_generators[0:n_b]
					#g_c = self.rxn_generators[n_b:n_b+n_c]
					
					
					#count_start_a = count_end_a
					count_start_b = count_end_b
					count_start_c = count_end_c
					
					"""
					Getting the uncertainty data
					"""
					#fig = plt.figure()
					data = self.unsrt[i].data
					T = np.linspace(300,2500,8)
					kappa_max = self.unsrt[i].getKappaMax(T)			
					kappa_min = self.unsrt[i].getKappaMin(T)
					Theta = np.array([T/T,np.log(T),-1/T])
					P = self.unsrt[i].getMean()	
					cov = self.unsrt[i].getCov()
					zeta_A = self.unsrt[i].zeta.x
					
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					"""
					Type-A zeta samples
					"""
					#generators_a = g_a
					parallel_zetas_a = []
					
					time_start = time.time()
					
					
					#for k in g_a:
						#Pint = P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten()
					#	plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
					#	plt.plot()
						#print(i)
						#print(type(i))
						#print(k)	
						#parallel_zetas_a.append(P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten())
					#	parallel_zetas_a.append(k[0]*zeta_A)
										
					
					
					
					
					parallel_zetas_a = []
					parallel_zetas_a_ = []
					
					time_start = time.time()
					
					
					#for k in g_a:
						
					#	parallel_zetas_a.append(k[0]*zeta_A)
					
					#for k in g_a_:
						
					#	parallel_zetas_a_.append(k[0]*zeta_A)
										
					
					generators_b_ = []
					parallel_zetas_b_ = []
					generators_c_ = []
					parallel_zetas_c_ = []
					
					generators_b = []
					parallel_zetas_b = []
					generators_c = []
					parallel_zetas_c = []
					start = 0
					start_ = 0
					for w in range(8):
						stop = int(((w+1)/8)*shuffled)
						stop_ = int(((w+1)/8)*unshuffled)
						g_b = self.rxn_generators_b2[start:stop]
						g_c = self.rxn_generators_c2[start:stop]
						g_b_ = self.rxn_generators_b2_[start_:stop_]
						g_c_ = self.rxn_generators_c2_[start_:stop_]
						data["generators_b"] = g_b
						callWorkForce = Worker(self.allowed_count)	
						gb,pb = callWorkForce.do_unsrt_b(data,len(g_b))
						del callWorkForce
						
						generators_b.extend(gb)
						parallel_zetas_b.extend(pb)
						
						removed = data.pop("generators_b")
						data["generators_b"] = g_b_
						callWorkForce = Worker(self.allowed_count)	
						gb,pb = callWorkForce.do_unsrt_b(data,len(g_b_))
						del callWorkForce
						
						generators_b_.extend(gb)
						parallel_zetas_b_.extend(pb)
						
						removed = data.pop("generators_b")
						data["generators_c"] = g_c
						callWorkForce = Worker(self.allowed_count)	
						gc,pc = callWorkForce.do_unsrt_c(data,len(g_c))
						del callWorkForce
						
						generators_c.extend(gc)
						parallel_zetas_c.extend(pc)
						removed = data.pop("generators_c")
						
						data["generators_c"] = g_c_
						callWorkForce = Worker(self.allowed_count)	
						gc,pc = callWorkForce.do_unsrt_c(data,len(g_c_))
						del callWorkForce
						
						generators_c_.extend(gc)
						parallel_zetas_c_.extend(pc)
						removed = data.pop("generators_c")
						
						start = stop
						start_ = stop_
					
					#print(self.sim_)
					#print(len(parallel_zetas_b)+len(parallel_zetas_c))
					g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					
					g_b_ = self.rxn_generators_b2_
					g_c_ = self.rxn_generators_c2_
					
					
					#temp_a[i] = parallel_zetas_a
					#gen_a[i] = generators_a
					temp_b[i] = parallel_zetas_b
					gen_b[i] = generators_b
					temp_c[i] = parallel_zetas_c
					gen_c[i] = generators_c
					
					#temp_a_[i] = parallel_zetas_a_
					#gen_a_[i] = generators_a_
					temp_b_[i] = parallel_zetas_b_
					gen_b_[i] = generators_b_
					temp_c_[i] = parallel_zetas_c_
					gen_c_[i] = generators_c_
					
					time_end = time.time()
					print(f"\n\t\tTime taken = {time_end-time_start}\n")
					
					#print(parallel_zetas_b)
					#Y_a = parallel_zetas_a
					Y_b = parallel_zetas_b
					Y_c = parallel_zetas_c
					
					#Y_a_ = parallel_zetas_a_
					Y_b_ = parallel_zetas_b_
					Y_c_ = parallel_zetas_c_
					T_ = np.linspace(300,2500,8)
					theta = np.array([T/T,np.log(T),-1/T])
					kappa_mean = theta.T.dot(P)
					Theta = np.array([T/T,np.log(T),-1/T])
					X = []
					Y = []
					
					X_ = []
					Y_ = []
					
					
						
					
					for j in Y_b:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
					
					for j in Y_b_:
						Y_.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X_.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					
					for j in Y_c:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_c.pdf")
					
					for j in Y_c_:
						Y_.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X_.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_c.pdf")
					
					temp_xdata = {}
					for k in range(len(Y_)):
						temp_xdata[str(Y_[k])] = X_[k]
						
					random.shuffle(Y_)
					#print(Y_)
					X1_ = []
					for k in range(len(Y_)):
						X1_.append(temp_xdata[str(Y_[k])])
						
					#print(X_[9])
					xdata[i] = X
					ydata[i] = Y
					
					xdata_[i] = X1_
					ydata_[i] = Y_
					
					#print(len(X1_))
					#print(len(Y_))
					#print(len(X))
					#print(len(Y))
					print(self.sim_)
					#Was trying to shuffle the zetas but the try was not successfull
					"""
					Y_ = []
					a = int(unshuffled/3)
					b = unshuffled - 2*a
					for w in Y_a[0:a]:
						Y_.append(w)
					for w in Y_b[0:a]:
						Y_.append(w)
					for w in Y_c[0:b]:
						Y_.append(w)
					
					random.shuffle(Y_)
					ydata_[i]=Y_
					
					#ydata[i] = Y
					#print(len(Y))
					#print(len(Y_a)+len(Y_b)+len(Y_c)))
					extra = len(Y_a[0:a])+len(Y_b[0:a])+len(Y_c[0:b])
					#print(len(parallel_zetas))
					#print(len(self.rxn_generators))
					"""
				A = []
				A_gen = []
				X_dash = []
				y_dash = []
				
				
				
				for j in range(len(g_b)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				
				for j in range(len(g_b_)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b_:
						a_row.extend(list(ydata_[i][j]))
						a_row_gen.extend(list(gen_b_[i][j]))
						x_temp.extend(list(xdata_[i][j]))
						y_temp.extend(list(ydata_[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				for j in range(n_b):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				
				for j in range(len(g_c)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_b)]))
						y_temp.extend(list(ydata[i][j+len(g_b)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				
				#print(len())
				for j in range(len(g_c_)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c_:
						#print(len(ydata_[i]))
						#print(len(g_a_)+len(g_b_)+len(g_c_))
						a_row.extend(list(ydata_[i][j+len(g_b_)]))
						a_row_gen.extend(list(gen_c_[i][j]))
						x_temp.extend(list(xdata_[i][j+len(g_b_)]))
						y_temp.extend(list(ydata_[i][j+len(g_b_)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				
				
				
				"""
				for j in range(extra):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_a:
						a_row.extend(list(ydata_[i][j]))
						y_temp.extend(list(ydata_[i][j]))
					A.append(np.asarray(a_row))
					y_dash.append(np.asarray(y_temp))
				"""
				"""
				for j in range(n_c):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+n_b]))
						y_temp.extend(list(ydata[i][j+n_b]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
					
				"""
				
				#print(A)
				#self.sim_ = len(g_a)+len(g_b)+len(g_c)
				#self.sim_ = n_b+n_c
				self.beta_ = np.asarray(A)
				#print(self.beta_)
				self.generators = np.asarray(A_gen)
				self.X = np.asarray(X_dash)
				self.y = np.asarray(y_dash)
			
			
			elif self.design == "A1+B1":
				self.beta_=[]
				if "Data" in os.listdir():
					os.chdir("Data/Simulations")
					if "Beta_list_case-0.csv" in os.listdir():
						#print(great)
						#print(os.getcwd())
						
						A_file = open("Beta_list_case-0.csv","r").readlines() 
						self.sim_ = len(A_file)
						print(self.sim_)
						A = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in A_file])	
						A_gen_file = open("Generator_list_case-0.csv","r").readlines() 
						A_gen = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in A_gen_file])	
						
						
						X_file = open("NormalizedBeta_list_case-0.csv","r").readlines() 
						X_dash = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in A_file])
						y_dash = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in A_file])
						self.beta_ = np.asarray(A)
						#print(self.beta_)
						self.generators = np.asarray(A_gen)
						
						self.X = np.asarray(X_dash)
						self.y = np.asarray(y_dash)
						os.chdir("../..")
					else:
						print("awe")
				else:
					print(os.getcwd())
					temp_a = {}
					temp_b = {}
					temp_c = {}
					gen = {}
					gen_a = {}
					gen_b = {}
					gen_c = {}
					xdata = {}
					ydata = {}
					
					temp_a_ = {}
					temp_b_ = {}
					temp_c_ = {}
					gen_ = {}
					gen_a_ = {}
					gen_b_ = {}
					gen_c_ = {}
					
					xdata_ = {}
					ydata_ = {}				
					self.generator_list = []
					"""
					Min- simulation = 300
					"""
					#if self.sim_ < 300:
					#	self.sim_ = 300
					#print(self.sim_)
					"""
					Dividing the generators into three bins
					------------------------------------
					n_a = 0.33*self.sim_
					n_b = 0.33*self.sim_
					n_c = 0.34*self.sim_
					
					"""
					#shuffled = int(0.8*self.sim_)
					#unshuffled = self.sim_ - shuffled
					
					#shuffled = int(0.2*self.sim_)
					#unshuffled = self.sim_ - shuffled
					
					shuffled = 100
					unshuffled = self.sim_
					
					n_a = int(0.1*shuffled)
					n_b = int(0.9*shuffled)
					n_c = shuffled-n_a-n_b
					
					n_a_ = int(0.1*unshuffled)
					n_b_ = int(0.9*unshuffled)
					n_c_ = unshuffled-n_a_-n_b_
					
					#n_b = int(0.5*self.sim_)
					#n_c = self.sim_-n_b
					
					
					#for i in np.eye(self.n):
					#	self.generator_list.append(i)
					#	self.generator_list.append(-1*i)
									
					"""
					Mixing three types of the extreme curves
					--------------------
					zeta-A: 1
					zeta-B: 1
					zeta-C: 1
				
					for zeta-A: one normal random parameter is required
						a matrix on n_a x n will be created
						
					for zeta-B and zeta-c: two random parameters is required
						we can arrange the two normal_random_parameter 				
						a matrix of n_b+n_c x 2*n will be created
						
					The number of reactions = n_R
					"""
					
					
					n_rxn = len(self.ind)
					self.generator_list_A = []
					self.generator_list_B = []
					self.generator_list_C = []
					
					self.generator_list_A_ = []
					self.generator_list_B_ = []
					self.generator_list_C_ = []
					
					"""
					Factorial Design
					----------------
					box-benkhen design
					bb = []
					for i in bbdesign(n_rxn)[0:int(0.3*self.sim_)]:
						bb.append(np.asarray(i))
					"""
							#For zeta-A
					"""	
					bb_a = []
					for i in bbdesign(n_rxn)[0:int(0.4*n_a)]:
						bb_a.append(np.asarray(i))
					self.generator_list_A.extend(bb_a)
					
							#For zeta-B
							
					bb_b = []
					for i in bbdesign(2*n_rxn)[0:int(0.4*n_b)]:
						bb_b.append(np.asarray(i))							
					self.generator_list_B.extend(bb_b)	
					
							#For zeta-C
							
					bb_c = []
					for i in bbdesign(2*n_rxn)[0:int(0.4*n_c)]:
						bb_c.append(np.asarray(i))			
					self.generator_list_C.extend(bb_c)		
					"""
							
					"""	
					bb_a = []
					for i in bbdesign(n_rxn):
						bb_a.append(np.asarray(i))
					self.generator_list_A.extend(bb_a)
					
							#For zeta-B
							
					bb_b = []
					for i in bbdesign(2*n_rxn):
						bb_b.append(np.asarray(i))							
					self.generator_list_B.extend(bb_b)	
					
							#For zeta-C
							
					bb_c = []
					for i in bbdesign(2*n_rxn):
						bb_c.append(np.asarray(i))			
					self.generator_list_C.extend(bb_c)		
							
					"""
					"""
					Shuffled cases
					"""
					"""	
					bb_a_ = []
					for i in bbdesign(n_rxn):
						bb_a_.append(np.asarray(i))
					self.generator_list_A_.extend(bb_a_)
					
							#For zeta-B
							
					bb_b_ = []
					for i in bbdesign(2*n_rxn):
						bb_b_.append(np.asarray(i))							
					self.generator_list_B_.extend(bb_b_)	
					
							#For zeta-C
							
					bb_c_ = []
					for i in bbdesign(2*n_rxn):
						bb_c_.append(np.asarray(i))			
					self.generator_list_C_.extend(bb_c_)		
					
					shuffled = len(self.generator_list_A)+len(self.generator_list_B)#+len(self.generator_list_C)
					unshuffled = len(self.generator_list_A_)+len(self.generator_list_B_)#+len(self.generator_list_C_)
					
					"""
					"""
					Monte-Carlo
					
					"""
							
							#For zeta-A
					
					self.generator_list_A.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(n_a)])
					self.generator_list_A.extend(list(np.eye(n_rxn)))
					self.generator_list_A.extend(list(-1*np.eye(n_rxn)))
					self.generator_list_A.append(np.zeros(n_rxn))
					self.generator_list_A.append(np.zeros(n_rxn))
					self.generator_list_A.append(np.zeros(n_rxn))
					
							#For zeta-b
							
					self.generator_list_B.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(n_b)])	
					self.generator_list_B.extend(list(np.eye(2*n_rxn)))
					self.generator_list_B.extend(list(-1*np.eye(2*n_rxn)))
					self.generator_list_B.append(np.zeros(2*n_rxn))
					self.generator_list_B.append(np.zeros(2*n_rxn))
					self.generator_list_B.append(np.zeros(2*n_rxn))
							
							#For zeta-C		
					self.generator_list_C.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(n_c)])
					self.generator_list_C.extend(list(np.eye(2*n_rxn)))
					self.generator_list_C.extend(list(-1*np.eye(2*n_rxn)))
					self.generator_list_C.append(np.zeros(2*n_rxn))
					self.generator_list_C.append(np.zeros(2*n_rxn))
					self.generator_list_C.append(np.zeros(2*n_rxn))
					
					
					self.generator_list_A_.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(n_a_)])
					#self.generator_list_A_.extend(list(np.eye(n_rxn)))
					#self.generator_list_A_.extend(list(-1*np.eye(n_rxn)))
					#self.generator_list_A_.append(np.zeros(n_rxn))
					#self.generator_list_A_.append(np.zeros(n_rxn))
					#self.generator_list_A_.append(np.zeros(n_rxn))
					
							#For zeta-A
							
					self.generator_list_B_.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(n_b_)])	
					#self.generator_list_B_.extend(list(np.eye(2*n_rxn)))
					#self.generator_list_B_.extend(list(-1*np.eye(2*n_rxn)))
					#self.generator_list_B_.append(np.zeros(2*n_rxn))
					#self.generator_list_B_.append(np.zeros(2*n_rxn))
					#self.generator_list_B_.append(np.zeros(2*n_rxn))
							
							#For zeta-A		
					self.generator_list_C_.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(n_c_)])
					#self.generator_list_C_.extend(list(np.eye(2*n_rxn)))
					#self.generator_list_C_.extend(list(-1*np.eye(2*n_rxn)))
					#self.generator_list_C_.append(np.zeros(2*n_rxn))
					#self.generator_list_C_.append(np.zeros(2*n_rxn))
					#self.generator_list_C_.append(np.zeros(2*n_rxn))
					
					
					
					"""
					For testing the distribution type
					"""
					"""		
							#For zeta-A
					
					self.generator_list_A.extend(list(2* np.random.random_sample((int(0.33*shuffled),n_rxn)) - 1))
					self.generator_list_A.extend(list(np.eye(n_rxn)))
					self.generator_list_A.extend(list(-1*np.eye(n_rxn)))
					self.generator_list_A.append(np.zeros(n_rxn))
					self.generator_list_A.append(np.zeros(n_rxn))
					self.generator_list_A.append(np.zeros(n_rxn))
					
							#For zeta-b
							
					self.generator_list_B.extend(list(2* np.random.random_sample((int(0.33*shuffled),2*n_rxn)) - 1))	
					self.generator_list_B.extend(list(np.eye(2*n_rxn)))
					self.generator_list_B.extend(list(-1*np.eye(2*n_rxn)))
					self.generator_list_B.append(np.zeros(2*n_rxn))
					self.generator_list_B.append(np.zeros(2*n_rxn))
					self.generator_list_B.append(np.zeros(2*n_rxn))
							
							#For zeta-C		
					self.generator_list_C.extend(list(2* np.random.random_sample((int(0.33*shuffled),2*n_rxn)) - 1))
					self.generator_list_C.extend(list(np.eye(2*n_rxn)))
					self.generator_list_C.extend(list(-1*np.eye(2*n_rxn)))
					self.generator_list_C.append(np.zeros(2*n_rxn))
					self.generator_list_C.append(np.zeros(2*n_rxn))
					self.generator_list_C.append(np.zeros(2*n_rxn))
					
					
					self.generator_list_A_.extend(list(2* np.random.random_sample((int(0.33*unshuffled),n_rxn)) - 1))
					#self.generator_list_A_.extend(list(np.eye(n_rxn)))
					#self.generator_list_A_.extend(list(-1*np.eye(n_rxn)))
					#self.generator_list_A_.append(np.zeros(n_rxn))
					#self.generator_list_A_.append(np.zeros(n_rxn))
					#self.generator_list_A_.append(np.zeros(n_rxn))
					
							#For zeta-A
							
					self.generator_list_B_.extend(list(2* np.random.random_sample((int(0.33*unshuffled),2*n_rxn)) - 1))	
					#self.generator_list_B_.extend(list(np.eye(2*n_rxn)))
					#self.generator_list_B_.extend(list(-1*np.eye(2*n_rxn)))
					#self.generator_list_B_.append(np.zeros(2*n_rxn))
					#self.generator_list_B_.append(np.zeros(2*n_rxn))
					#self.generator_list_B_.append(np.zeros(2*n_rxn))
							
							#For zeta-A		
					self.generator_list_C_.extend(list(2* np.random.random_sample((int(0.33*unshuffled),2*n_rxn)) - 1))
					#self.generator_list_C_.extend(list(np.eye(2*n_rxn)))
					#self.generator_list_C_.extend(list(-1*np.eye(2*n_rxn)))
					#self.generator_list_C_.append(np.zeros(2*n_rxn))
					#self.generator_list_C_.append(np.zeros(2*n_rxn))
					#self.generator_list_C_.append(np.zeros(2*n_rxn))
					"""
					"""
					Uniform
					
					"""
						
					#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
					
					"""
					Latin-Hypercube
					"""
							#For zeta-A
					"""		
					lhs_a = []
					for i in lhs(n_rxn, samples=int(0.3*n_a)):
						lhs_a.append(np.asarray(i))	
					self.generator_list_A.extend(lhs_a)
					
							
							#For zeta-A
							
					lhs_b = []
					for i in lhs(2*n_rxn, samples=int(0.3*n_b)):
						lhs_b.append(np.asarray(i))		
					self.generator_list_B.extend(lhs_b)
							
							#For zeta-A
							
					lhs_c = []
					for i in lhs(2*n_rxn, samples=int(0.3*n_c)):
						lhs_c.append(np.asarray(i))
					self.generator_list_C.extend(lhs_c)
					
					"""
					"""		
					lhs_a = []
					for i in lhs(n_rxn, samples=int(0.33*self.sim_)):
						lhs_a.append(np.asarray(i))	
					self.generator_list_A.extend(lhs_a)
					
							
							#For zeta-A
							
					lhs_b = []
					for i in lhs(2*n_rxn, samples=int(0.33*self.sim_)):
						lhs_b.append(np.asarray(i))		
					self.generator_list_B.extend(lhs_b)
							
							#For zeta-A
							
					lhs_c = []
					for i in lhs(2*n_rxn, samples=int(0.33*self.sim_)):
						lhs_c.append(np.asarray(i))
					self.generator_list_C.extend(lhs_c)
					
					"""
					"""
					lhs_a_ = []
					for i in lhs(n_rxn, samples=int(0.33*unshuffled)):
						lhs_a_.append(np.asarray(i))	
					self.generator_list_A_.extend(lhs_a_)
					
							
							#For zeta-A
							
					lhs_b_ = []
					for i in lhs(2*n_rxn, samples=int(0.33*unshuffled)):
						lhs_b_.append(np.asarray(i))		
					self.generator_list_B_.extend(lhs_b_)
							
							#For zeta-A
							
					lhs_c_ = []
					for i in lhs(2*n_rxn, samples=int(0.33*unshuffled)):
						lhs_c_.append(np.asarray(i))
					self.generator_list_C_.extend(lhs_c_)
					
					"""
					"""
					Monte-Carlo-Mix
					"""
					#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)		
					#self.generator_list.extend(bb)
					#self.generator_list.extend(list(np.eye(self.n)))
					#self.generator_list.extend(list(-1*np.eye(self.n)))
					#self.generator_list.append(np.zeros(self.n))
					#self.generator_list.append(np.zeros(self.n))
					#self.generator_list.append(np.zeros(self.n))
					#print(self.generator_list)
					#self.sim_ = len(self.generator_list)
					
					#trainers = np.asarray(trainers)
					#print(f"Generator_list {np.shape(self.generator_list)}\n")
					df = pd.DataFrame(np.asarray(self.generator_list))
					count_start_a = 0
					count_start_b = 0
					count_start_c = 0			
					for i in self.reactions_selected:
						
						len_of_Arrhenius_A = 1
						count_end_a = count_start_a+len_of_Arrhenius_A
						
						
						len_of_Arrhenius_B = 2
						count_end_b = count_start_b+len_of_Arrhenius_B
						
						
						len_of_Arrhenius_C = 2
						count_end_c = count_start_c+len_of_Arrhenius_C
						
						#print(count_start)
						#print(count_end)
						self.rxn_generators_a2 = []
						self.rxn_generators_b2 = []
						self.rxn_generators_c2 = []
						
						self.rxn_generators_a2_ = []
						self.rxn_generators_b2_ = []
						self.rxn_generators_c2_ = []
						
						"""
						getting zetas from generators
						"""
						for j in range(len(self.generator_list_A)):
							self.rxn_generators_a2.append(self.generator_list_A[j][count_start_a:count_end_a])
						#print(self.rxn_generators_a2)
						for j in range(len(self.generator_list_B)):
							self.rxn_generators_b2.append(self.generator_list_B[j][count_start_b:count_end_b])
						
						for j in range(len(self.generator_list_C)):
							self.rxn_generators_c2.append(self.generator_list_C[j][count_start_c:count_end_c])
						
						
						for j in range(len(self.generator_list_A_)):
							self.rxn_generators_a2_.append(self.generator_list_A_[j][count_start_a:count_end_a])
						#print(self.rxn_generators_a2)
						for j in range(len(self.generator_list_B_)):
							self.rxn_generators_b2_.append(self.generator_list_B_[j][count_start_b:count_end_b])
						
						for j in range(len(self.generator_list_C_)):
							self.rxn_generators_c2_.append(self.generator_list_C_[j][count_start_c:count_end_c])
						
						
						
						
						self.sim_ = len(self.rxn_generators_a2)+len(self.rxn_generators_b2)+len(self.rxn_generators_a2_)+len(self.rxn_generators_b2_)
						#print(self.rxn_generators)
						#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
											
						g_a = self.rxn_generators_a2
						g_b = self.rxn_generators_b2
						g_c = self.rxn_generators_c2
						
						g_a_ = self.rxn_generators_a2_
						g_b_ = self.rxn_generators_b2_
						g_c_ = self.rxn_generators_c2_
						
						#g_b = self.rxn_generators[0:n_b]
						#g_c = self.rxn_generators[n_b:n_b+n_c]
						
						
						count_start_a = count_end_a
						count_start_b = count_end_b
						count_start_c = count_end_c
						
						"""
						Getting the uncertainty data
						"""
						#fig = plt.figure()
						data = self.unsrt[i].data
						T = np.linspace(300,2500,8)
						kappa_max = self.unsrt[i].getKappaMax(T)			
						kappa_min = self.unsrt[i].getKappaMin(T)
						Theta = np.array([T/T,np.log(T),-1/T])
						P = self.unsrt[i].getMean()	
						cov = self.unsrt[i].getCov()
						zeta_A = self.unsrt[i].zeta.x
						
						#plt.plot(1/T,kappa_max,"k--")
						#plt.plot(1/T,kappa_min,"k--")
						"""
						Type-A zeta samples
						"""
						generators_a = g_a
						generators_a_ = g_a_
						parallel_zetas_a = []
						parallel_zetas_a_ = []
						
						time_start = time.time()
						
						
						for k in g_a:
							
							parallel_zetas_a.append(k[0]*zeta_A)
						
						for k in g_a_:
							
							parallel_zetas_a_.append(k[0]*zeta_A)
											
						
						generators_b_ = []
						parallel_zetas_b_ = []
						generators_c_ = []
						parallel_zetas_c_ = []
						
						generators_b = []
						parallel_zetas_b = []
						generators_c = []
						parallel_zetas_c = []
						start = 0
						start_ = 0
						for w in range(8):
							stop = int(((w+1)/8)*shuffled)
							stop_ = int(((w+1)/8)*unshuffled)
							g_b = self.rxn_generators_b2[start:stop]
							g_c = self.rxn_generators_c2[start:stop]
							g_b_ = self.rxn_generators_b2_[start_:stop_]
							g_c_ = self.rxn_generators_c2_[start_:stop_]
							data["generators_b"] = g_b
							callWorkForce = Worker(self.allowed_count)	
							gb,pb = callWorkForce.do_unsrt_b(data,len(g_b))
							del callWorkForce
							
							generators_b.extend(gb)
							parallel_zetas_b.extend(pb)
							
							removed = data.pop("generators_b")
							data["generators_b"] = g_b_
							callWorkForce = Worker(self.allowed_count)	
							gb_,pb_ = callWorkForce.do_unsrt_b(data,len(g_b_))
							del callWorkForce
							
							generators_b_.extend(gb_)
							parallel_zetas_b_.extend(pb_)
							
							removed = data.pop("generators_b")
							
							
						g_b = self.rxn_generators_b2
						
						g_b_ = self.rxn_generators_b2_
						
						
						temp_a[i] = parallel_zetas_a
						gen_a[i] = generators_a
						temp_b[i] = parallel_zetas_b
						gen_b[i] = generators_b
						
						temp_a_[i] = parallel_zetas_a_
						gen_a_[i] = generators_a_
						temp_b_[i] = parallel_zetas_b_
						gen_b_[i] = generators_b_
						
						time_end = time.time()
						print(f"\n\t\tTime taken = {time_end-time_start}\n")
						
						Y_a = parallel_zetas_a
						Y_b = parallel_zetas_b
						
						Y_a_ = parallel_zetas_a_
						Y_b_ = parallel_zetas_b_
						
						
						T_ = np.linspace(300,2500,8)
						theta = np.array([T/T,np.log(T),-1/T])
						kappa_mean = theta.T.dot(P)
						Theta = np.array([T/T,np.log(T),-1/T])
						X = []
						Y = []
						
						X_ = []
						Y_ = []
						
						
						for j in Y_a:
							Y.append(j)
							Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
							
							X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
							
						for j in Y_a_:
							Y_.append(j)
							Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
							
							X_.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
							
						
						for j in Y_b:
							Y.append(j)
							Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
							#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
							X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						
						for j in Y_b_:
							Y_.append(j)
							Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
							#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
							X_.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
							#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
						
						temp_xdata = {}
						for k in range(len(Y_)):
							temp_xdata[str(Y_[k])] = X_[k]
							
						random.shuffle(Y_)
						#print(Y_)
						X1_ = []
						for k in range(len(Y_)):
							X1_.append(temp_xdata[str(Y_[k])])
							
						#print(X_[9])
						xdata[i] = X
						ydata[i] = Y
						
						xdata_[i] = X1_
						ydata_[i] = Y_
						
						
						print(self.sim_)
						
					A = []
					A_gen = []
					X_dash = []
					y_dash = []
					
					for j in range(len(g_a)):
						a_row = []
						a_row_gen = []
						x_temp = []
						y_temp = []
						for i in temp_a:
							a_row.extend(list(temp_a[i][j]))
							a_row_gen.extend(list(gen_a[i][j]))
							x_temp.extend(list(xdata[i][j]))
							y_temp.extend(list(ydata[i][j]))
						A.append(np.asarray(a_row))
						A_gen.append(np.asarray(a_row_gen))
						X_dash.append(np.asarray(x_temp))
						y_dash.append(np.asarray(y_temp))
						
					for j in range(len(g_a_)):
						a_row = []
						a_row_gen = []
						x_temp = []
						y_temp = []
						for i in temp_a_:
							a_row.extend(list(ydata_[i][j]))
							a_row_gen.extend(list(gen_a_[i][j]))
							x_temp.extend(list(xdata_[i][j]))
							y_temp.extend(list(ydata_[i][j]))
						A.append(np.asarray(a_row))
						A_gen.append(np.asarray(a_row_gen))
						X_dash.append(np.asarray(x_temp))
						y_dash.append(np.asarray(y_temp))
					
					
					for j in range(len(g_b)):
						a_row = []
						a_row_gen = []
						x_temp = []
						y_temp = []
						for i in temp_b:
							a_row.extend(list(temp_b[i][j]))
							a_row_gen.extend(list(gen_b[i][j]))
							x_temp.extend(list(xdata[i][j+len(g_a)]))
							y_temp.extend(list(ydata[i][j+len(g_a)]))
						A.append(np.asarray(a_row))
						A_gen.append(np.asarray(a_row_gen))
						X_dash.append(np.asarray(x_temp))
						y_dash.append(np.asarray(y_temp))
					
					for j in range(len(g_b_)):
						a_row = []
						a_row_gen = []
						x_temp = []
						y_temp = []
						for i in temp_b_:
							a_row.extend(list(ydata_[i][j+len(g_a_)]))
							a_row_gen.extend(list(gen_b_[i][j]))
							x_temp.extend(list(xdata_[i][j+len(g_a_)]))
							y_temp.extend(list(ydata_[i][j+len(g_a_)]))
						A.append(np.asarray(a_row))
						A_gen.append(np.asarray(a_row_gen))
						X_dash.append(np.asarray(x_temp))
						y_dash.append(np.asarray(y_temp))
					"""
					for j in range(n_b):
						a_row = []
						a_row_gen = []
						x_temp = []
						y_temp = []
						for i in temp_b:
							a_row.extend(list(temp_b[i][j]))
							a_row_gen.extend(list(gen_b[i][j]))
							x_temp.extend(list(xdata[i][j]))
							y_temp.extend(list(ydata[i][j]))
						A.append(np.asarray(a_row))
						A_gen.append(np.asarray(a_row_gen))
						X_dash.append(np.asarray(x_temp))
						y_dash.append(np.asarray(y_temp))
					"""
					
					self.beta_ = np.asarray(A)
					#print(self.beta_)
					self.generators = np.asarray(A_gen)
					self.X = np.asarray(X_dash)
					self.y = np.asarray(y_dash)
			
			elif self.design == "B1-only":
				self.beta_=[]
				temp_a = {}
				temp_b = {}
				temp_c = {}
				gen = {}
				gen_a = {}
				gen_b = {}
				gen_c = {}
				
				xdata = {}
				ydata = {}				
				self.generator_list = []
				"""
				Min- simulation = 300
				"""
				#if self.sim_ < 300:
				#	self.sim_ = 300
				#print(self.sim_)
				"""
				Dividing the generators into three bins
				------------------------------------
				n_a = 0.33*self.sim_
				n_b = 0.33*self.sim_
				n_c = 0.34*self.sim_
				
				"""
				n_a = int(0.5*self.sim_)
				#n_b = int(0.5*self.sim_)
				#n_c = self.sim_-n_a-n_b
				
				
				n_b = int(self.sim_)
				n_c = int(self.sim_-0.33*n_b)
				
				
				#for i in np.eye(self.n):
				#	self.generator_list.append(i)
				#	self.generator_list.append(-1*i)
								
				"""
				Mixing three types of the extreme curves
				--------------------
				zeta-A: 1
				zeta-B: 1
				zeta-C: 1
			
				for zeta-A: one normal random parameter is required
					a matrix on n_a x n will be created
					
				for zeta-B and zeta-c: two random parameters is required
					we can arrange the two normal_random_parameter 				
					a matrix of n_b+n_c x 2*n will be created
					
				The number of reactions = n_R
				"""
				n_rxn = len(self.ind)
				self.generator_list_A = []
				self.generator_list_B = []
				self.generator_list_C = []
				
				"""
				Factorial Design
				----------------
				box-benkhen design
				bb = []
				for i in bbdesign(n_rxn)[0:int(0.3*self.sim_)]:
					bb.append(np.asarray(i))
				"""
						#For zeta-A
						
				bb_a = []
				for i in bbdesign(n_rxn)[0:int(0.4*n_a)]:
					bb_a.append(np.asarray(i))
				self.generator_list_A.extend(bb_a)
						
						#For zeta-B
						
				bb_b = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_b)]:
					bb_b.append(np.asarray(i))							
				self.generator_list_B.extend(bb_b)
						
						#For zeta-C
						
				bb_c = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_c)]:
					bb_c.append(np.asarray(i))		
						
				self.generator_list_C.extend(bb_c)		
						
				
				"""
				Monte-Carlo
				
				"""
						
						#For zeta-A
				
				self.generator_list_A.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(0.3*n_a)])
				self.generator_list_A.extend(list(np.eye(n_rxn)))
				self.generator_list_A.extend(list(-1*np.eye(n_rxn)))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				
						#For zeta-A
						
				self.generator_list_B.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_b)])	
				self.generator_list_B.extend(list(np.eye(2*n_rxn)))
				self.generator_list_B.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
						
						#For zeta-A		
				self.generator_list_C.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_c)])
				self.generator_list_C.extend(list(np.eye(2*n_rxn)))
				self.generator_list_C.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				
				"""
				Uniform
				
				"""
					
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				
				"""
				Latin-Hypercube
				"""
			
						
						#For zeta-A
						
				lhs_a = []
				for i in lhs(n_rxn, samples=int(0.3*n_a)):
					lhs_a.append(np.asarray(i))	
				self.generator_list_A.extend(lhs_a)
				
						
						#For zeta-A
						
				lhs_b = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_b)):
					lhs_b.append(np.asarray(i))		
				self.generator_list_B.extend(lhs_b)
						
						#For zeta-A
						
				lhs_c = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_c)):
					lhs_c.append(np.asarray(i))
				self.generator_list_C.extend(lhs_c)
				
				"""
				Monte-Carlo-Mix
				"""
				
				
				df = pd.DataFrame(np.asarray(self.generator_list))
				count_start_a = 0
				count_start_b = 0
				count_start_c = 0			
				for i in self.reactions_selected:
					
					len_of_Arrhenius_A = 1
					count_end_a = count_start_a+len_of_Arrhenius_A
					
					
					len_of_Arrhenius_B = 2
					count_end_b = count_start_b+len_of_Arrhenius_B
					
					
					len_of_Arrhenius_C = 2
					count_end_c = count_start_c+len_of_Arrhenius_C
					
					#print(count_start)
					#print(count_end)
					self.rxn_generators_a2 = []
					self.rxn_generators_b2 = []
					self.rxn_generators_c2 = []
					
					"""
					getting zetas from generators
					"""
					for j in range(len(self.generator_list_A)):
						self.rxn_generators_a2.append(self.generator_list_A[j][count_start_a:count_end_a])
					#print(self.rxn_generators_a2)
					for j in range(len(self.generator_list_B)):
						self.rxn_generators_b2.append(self.generator_list_B[j][count_start_b:count_end_b])
					
					for j in range(len(self.generator_list_C)):
						self.rxn_generators_c2.append(self.generator_list_C[j][count_start_c:count_end_c])
					
					#self.sim_ = len(self.rxn_generators_a2)+len(self.rxn_generators_b2)+len(self.rxn_generators_c2)
					self.sim_ = len(self.rxn_generators_b2)
					#print(self.rxn_generators)
					#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
										
					g_a = self.rxn_generators_a2
					g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					#g_b = self.rxn_generators[0:n_b]
					#g_c = self.rxn_generators[n_b:n_b+n_c]
					
					
					count_start_a = count_end_a
					count_start_b = count_end_b
					count_start_c = count_end_c
					
					"""
					Getting the uncertainty data
					"""
					#fig = plt.figure()
					data = self.unsrt[i].data
					T = np.linspace(300,2500,8)
					kappa_max = self.unsrt[i].getKappaMax(T)			
					kappa_min = self.unsrt[i].getKappaMin(T)
					Theta = np.array([T/T,np.log(T),-1/T])
					P = self.unsrt[i].getMean()	
					cov = self.unsrt[i].getCov()
					zeta_A = self.unsrt[i].zeta.x
					
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					"""
					Type-A zeta samples
					"""
					generators_a = g_a
					parallel_zetas_a = []
					
					time_start = time.time()
					
					
					for k in g_a:
						
						parallel_zetas_a.append(k[0]*zeta_A)
										
					data["generators_b"] = g_b
					callWorkForce = Worker(self.allowed_count)	
					generators_b,parallel_zetas_b = callWorkForce.do_unsrt_b(data,len(g_b))
					del callWorkForce
					"""
					data["generators_c"] = g_c
					callWorkForce = Worker(self.allowed_count)	
					generators_c,parallel_zetas_c = callWorkForce.do_unsrt_c(data,len(g_c))
					del callWorkForce
					"""
					#print(parallel_zetas_a)
					
					temp_a[i] = parallel_zetas_a
					gen_a[i] = generators_a
					temp_b[i] = parallel_zetas_b
					gen_b[i] = generators_b
					#temp_c[i] = parallel_zetas_c
					#gen_c[i] = generators_c
					
					time_end = time.time()
					print(f"\n\t\tTime taken = {time_end-time_start}\n")
					
					#print(parallel_zetas_b)
					Y_a = parallel_zetas_a
					Y_b = parallel_zetas_b
					#Y_c = parallel_zetas_c
					T_ = np.linspace(300,2500,8)
					theta = np.array([T/T,np.log(T),-1/T])
					kappa_mean = theta.T.dot(P)
					Theta = np.array([T/T,np.log(T),-1/T])
					X = []
					Y = []
					
					for j in Y_b:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						
					
					temp_xdata = {}
					for k in range(len(Y)):
						temp_xdata[str(Y[k])] = X[k]
						
					random.shuffle(Y)
					#print(Y_)
					X_ = []
					for k in range(len(Y)):
						X_.append(temp_xdata[str(Y[k])])
						
					#print(X_[9])
					xdata[i] = X_
					ydata[i] = Y
					#print(len(parallel_zetas))
					#print(len(self.rxn_generators))
				A = []
				A_gen = []
				X_dash = []
				y_dash = []
								
				for j in range(len(g_b)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				
				
				
				self.beta_ = np.asarray(A)
				#print(self.beta_)
				self.generators = np.asarray(A_gen)
				self.X = np.asarray(X_dash)
				self.y = np.asarray(y_dash)
			
			elif self.design == "C1-only":
				self.beta_=[]
				temp_a = {}
				temp_b = {}
				temp_c = {}
				gen = {}
				gen_a = {}
				gen_b = {}
				gen_c = {}
				
				xdata = {}
				ydata = {}				
				self.generator_list = []
				"""
				Min- simulation = 300
				"""
				#if self.sim_ < 300:
				#	self.sim_ = 300
				#print(self.sim_)
				"""
				Dividing the generators into three bins
				------------------------------------
				n_a = 0.33*self.sim_
				n_b = 0.33*self.sim_
				n_c = 0.34*self.sim_
				
				"""
				n_a = int(0.5*self.sim_)
				#n_b = int(0.5*self.sim_)
				#n_c = self.sim_-n_a-n_b
				
				
				n_b = int(self.sim_)
				n_c = self.sim_
				
				
				#for i in np.eye(self.n):
				#	self.generator_list.append(i)
				#	self.generator_list.append(-1*i)
								
				"""
				Mixing three types of the extreme curves
				--------------------
				zeta-A: 1
				zeta-B: 1
				zeta-C: 1
			
				for zeta-A: one normal random parameter is required
					a matrix on n_a x n will be created
					
				for zeta-B and zeta-c: two random parameters is required
					we can arrange the two normal_random_parameter 				
					a matrix of n_b+n_c x 2*n will be created
					
				The number of reactions = n_R
				"""
				n_rxn = len(self.ind)
				self.generator_list_A = []
				self.generator_list_B = []
				self.generator_list_C = []
				
				"""
				Factorial Design
				----------------
				box-benkhen design
				bb = []
				for i in bbdesign(n_rxn)[0:int(0.3*self.sim_)]:
					bb.append(np.asarray(i))
				"""
						#For zeta-A
						
				bb_a = []
				for i in bbdesign(n_rxn)[0:int(0.4*n_a)]:
					bb_a.append(np.asarray(i))
				self.generator_list_A.extend(bb_a)
						
						#For zeta-B
						
				bb_b = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_b)]:
					bb_b.append(np.asarray(i))							
				self.generator_list_B.extend(bb_b)
						
						#For zeta-C
						
				bb_c = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_c)]:
					bb_c.append(np.asarray(i))		
						
				self.generator_list_C.extend(bb_c)		
						
				
				"""
				Monte-Carlo
				
				"""
						
						#For zeta-A
				
				self.generator_list_A.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(0.3*n_a)])
				self.generator_list_A.extend(list(np.eye(n_rxn)))
				self.generator_list_A.extend(list(-1*np.eye(n_rxn)))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				self.generator_list_A.append(np.zeros(n_rxn))
				
						#For zeta-A
						
				self.generator_list_B.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_b)])	
				self.generator_list_B.extend(list(np.eye(2*n_rxn)))
				self.generator_list_B.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
				self.generator_list_B.append(np.zeros(2*n_rxn))
						
						#For zeta-A		
				self.generator_list_C.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_c)])
				self.generator_list_C.extend(list(np.eye(2*n_rxn)))
				self.generator_list_C.extend(list(-1*np.eye(2*n_rxn)))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				self.generator_list_C.append(np.zeros(2*n_rxn))
				
				"""
				Uniform
				
				"""
					
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				
				"""
				Latin-Hypercube
				"""
			
						
						#For zeta-A
						
				lhs_a = []
				for i in lhs(n_rxn, samples=int(0.3*n_a)):
					lhs_a.append(np.asarray(i))	
				self.generator_list_A.extend(lhs_a)
				
						
						#For zeta-A
						
				lhs_b = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_b)):
					lhs_b.append(np.asarray(i))		
				self.generator_list_B.extend(lhs_b)
						
						#For zeta-A
						
				lhs_c = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_c)):
					lhs_c.append(np.asarray(i))
				self.generator_list_C.extend(lhs_c)
				
				"""
				Monte-Carlo-Mix
				"""
				
				
				df = pd.DataFrame(np.asarray(self.generator_list))
				count_start_a = 0
				count_start_b = 0
				count_start_c = 0			
				for i in self.reactions_selected:
					
					len_of_Arrhenius_A = 1
					count_end_a = count_start_a+len_of_Arrhenius_A
					
					
					len_of_Arrhenius_B = 2
					count_end_b = count_start_b+len_of_Arrhenius_B
					
					
					len_of_Arrhenius_C = 2
					count_end_c = count_start_c+len_of_Arrhenius_C
					
					#print(count_start)
					#print(count_end)
					self.rxn_generators_a2 = []
					self.rxn_generators_b2 = []
					self.rxn_generators_c2 = []
					
					"""
					getting zetas from generators
					"""
					for j in range(len(self.generator_list_A)):
						self.rxn_generators_a2.append(self.generator_list_A[j][count_start_a:count_end_a])
					#print(self.rxn_generators_a2)
					for j in range(len(self.generator_list_B)):
						self.rxn_generators_b2.append(self.generator_list_B[j][count_start_b:count_end_b])
					
					for j in range(len(self.generator_list_C)):
						self.rxn_generators_c2.append(self.generator_list_C[j][count_start_c:count_end_c])
					
					#self.sim_ = len(self.rxn_generators_a2)+len(self.rxn_generators_b2)+len(self.rxn_generators_c2)
					self.sim_ = len(self.rxn_generators_c2)
					#print(self.rxn_generators)
					#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
										
					g_a = self.rxn_generators_a2
					g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					#g_b = self.rxn_generators[0:n_b]
					#g_c = self.rxn_generators[n_b:n_b+n_c]
					
					
					count_start_a = count_end_a
					count_start_b = count_end_b
					count_start_c = count_end_c
					
					"""
					Getting the uncertainty data
					"""
					#fig = plt.figure()
					data = self.unsrt[i].data
					T = np.linspace(300,2500,8)
					kappa_max = self.unsrt[i].getKappaMax(T)			
					kappa_min = self.unsrt[i].getKappaMin(T)
					Theta = np.array([T/T,np.log(T),-1/T])
					P = self.unsrt[i].getMean()	
					cov = self.unsrt[i].getCov()
					zeta_A = self.unsrt[i].zeta.x
					
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					"""
					Type-A zeta samples
					"""
					generators_a = g_a
					parallel_zetas_a = []
					
					time_start = time.time()
					
					
					for k in g_a:
						
						parallel_zetas_a.append(k[0]*zeta_A)
					
					generators_c = []
					parallel_zetas_c = []
					start = 0
					for w in range(10):
						stop = int(((w+1)/10)*self.sim_)
						g_c = self.rxn_generators_c2[start:stop]
						data["generators_c"] = g_c
						callWorkForce = Worker(self.allowed_count)	
						generators_c_,parallel_zetas_c_ = callWorkForce.do_unsrt_c(data,len(g_c))
						del callWorkForce
						generators_c.extend(generators_c_)
						parallel_zetas_c.extend(parallel_zetas_c_)
						removed = data.pop("generators_c")
						start = stop
					
					#print(self.sim_)
					#print(len(parallel_zetas_b)+len(parallel_zetas_c))
					#g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					
					temp_a[i] = parallel_zetas_a
					gen_a[i] = generators_a
					#temp_b[i] = parallel_zetas_b
					#gen_b[i] = generators_b
					temp_c[i] = parallel_zetas_c
					gen_c[i] = generators_c
					
					time_end = time.time()
					print(f"\n\t\tTime taken = {time_end-time_start}\n")
					
					#print(parallel_zetas_b)
					Y_a = parallel_zetas_a
					#Y_b = parallel_zetas_b
					Y_c = parallel_zetas_c
					T_ = np.linspace(300,2500,8)
					theta = np.array([T/T,np.log(T),-1/T])
					kappa_mean = theta.T.dot(P)
					Theta = np.array([T/T,np.log(T),-1/T])
					X = []
					Y = []
					
					for j in Y_c:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						
					
					temp_xdata = {}
					for k in range(len(Y)):
						temp_xdata[str(Y[k])] = X[k]
						
					random.shuffle(Y)
					#print(Y_)
					X_ = []
					for k in range(len(Y)):
						X_.append(temp_xdata[str(Y[k])])
						
					#print(X_[9])
					xdata[i] = X_
					ydata[i] = Y
				A = []
				A_gen = []
				X_dash = []
				y_dash = []
								
				for j in range(len(g_c)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(ydata[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				self.beta_ = np.asarray(A)
				#print(self.beta_)
				self.generators = np.asarray(A_gen)
				self.X = np.asarray(X_dash)
				self.y = np.asarray(y_dash)
			
			
			elif self.design == "A1+B1+C1":
				self.beta_=[]
				if "Data" in os.listdir():
					#print(os.getcwd())
					os.chdir("Data/Simulations")
					if "Beta_list_case-0.csv" in os.listdir():
						#print(great)
						#print(os.getcwd())
						
						A_file = open("Beta_list_case-0.csv","r").readlines() 
						self.sim_ = len(A_file)
						print(self.sim_)
						#raise AssertionError("Stop\n")
						A = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in A_file])	
						A_gen_file = open("Generator_list_case-0.csv","r").readlines() 
						A_gen = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in A_gen_file])	
						
						
						X_file = open("NormalizedBeta_list_case-0.csv","r").readlines() 
						X_dash = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in A_file])
						y_dash = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in A_file])
						self.beta_ = np.asarray(A)
						#print(self.beta_)
						self.generators = np.asarray(A_gen)
						
						self.X = np.asarray(X_dash)
						self.y = np.asarray(y_dash)
						os.chdir("../..")
				else:
				
					self.beta_=[]
					temp_a = {}
					temp_b = {}
					temp_c = {}
					gen = {}
					gen_a = {}
					gen_b = {}
					gen_c = {}
					xdata = {}
					ydata = {}
					
					temp_a_ = {}
					temp_b_ = {}
					temp_c_ = {}
					gen_ = {}
					gen_a_ = {}
					gen_b_ = {}
					gen_c_ = {}
					
					xdata_ = {}
					ydata_ = {}				
					self.generator_list = []
					"""
					Min- simulation = 300
					"""
					#if self.sim_ < 300:
					#	self.sim_ = 300
					#print(self.sim_)
					"""
					Dividing the generators into three bins
					------------------------------------
					n_a = 0.33*self.sim_
					n_b = 0.33*self.sim_
					n_c = 0.34*self.sim_
					
					"""
					shuffled = int(0.8*self.sim_)
					unshuffled = self.sim_ - shuffled
					
					n_a = int(0.1*shuffled)
					n_b = int(0.45*shuffled)
					n_c = shuffled-n_a-n_b
					
					n_a_ = int(0.1*unshuffled)
					n_b_ = int(0.45*unshuffled)
					n_c_ = unshuffled-n_a_-n_b_
					
					#n_b = int(0.5*self.sim_)
					#n_c = self.sim_-n_b
					
					
					#for i in np.eye(self.n):
					#	self.generator_list.append(i)
					#	self.generator_list.append(-1*i)
									
					"""
					Mixing three types of the extreme curves
					--------------------
					zeta-A: 1
					zeta-B: 1
					zeta-C: 1
				
					for zeta-A: one normal random parameter is required
						a matrix on n_a x n will be created
						
					for zeta-B and zeta-c: two random parameters is required
						we can arrange the two normal_random_parameter 				
						a matrix of n_b+n_c x 2*n will be created
						
					The number of reactions = n_R
					"""
					
					
					n_rxn = len(self.ind)
					self.generator_list_A = []
					self.generator_list_B = []
					self.generator_list_C = []
					
					self.generator_list_A_ = []
					self.generator_list_B_ = []
					self.generator_list_C_ = []
					
					"""
					Factorial Design
					----------------
					box-benkhen design
					bb = []
					for i in bbdesign(n_rxn)[0:int(0.3*self.sim_)]:
						bb.append(np.asarray(i))
					"""
							#For zeta-A
						
					bb_a = []
					for i in bbdesign(n_rxn)[0:int(0.4*n_a)]:
						bb_a.append(np.asarray(i))
					self.generator_list_A.extend(bb_a)
					
							#For zeta-B
							
					bb_b = []
					for i in bbdesign(2*n_rxn)[0:int(0.4*n_b)]:
						bb_b.append(np.asarray(i))							
					self.generator_list_B.extend(bb_b)	
					
							#For zeta-C
							
					bb_c = []
					for i in bbdesign(2*n_rxn)[0:int(0.4*n_c)]:
						bb_c.append(np.asarray(i))			
					self.generator_list_C.extend(bb_c)		
					
							
					"""
					bb_a = []
					for i in bbdesign(n_rxn):
						bb_a.append(np.asarray(i))
					self.generator_list_A.extend(bb_a)
					
							#For zeta-B
							
					bb_b = []
					for i in bbdesign(2*n_rxn):
						bb_b.append(np.asarray(i))							
					self.generator_list_B.extend(bb_b)	
					
							#For zeta-C
							
					bb_c = []
					for i in bbdesign(2*n_rxn):
						bb_c.append(np.asarray(i))			
					self.generator_list_C.extend(bb_c)		
							
					"""
					"""
					Shuffled cases
					"""
					"""	
					bb_a_ = []
					for i in bbdesign(n_rxn):
						bb_a_.append(np.asarray(i))
					self.generator_list_A_.extend(bb_a_)
					
							#For zeta-B
							
					bb_b_ = []
					for i in bbdesign(2*n_rxn):
						bb_b_.append(np.asarray(i))							
					self.generator_list_B_.extend(bb_b_)	
					
							#For zeta-C
							
					bb_c_ = []
					for i in bbdesign(2*n_rxn):
						bb_c_.append(np.asarray(i))			
					self.generator_list_C_.extend(bb_c_)		
					
					shuffled = len(self.generator_list_A)+len(self.generator_list_B)+len(self.generator_list_C)
					unshuffled = len(self.generator_list_A_)+len(self.generator_list_B_)+len(self.generator_list_C_)
					"""
					
					"""
					Monte-Carlo
					
					"""
							
							#For zeta-A
					
					self.generator_list_A.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(0.4*n_a)])
					self.generator_list_A.extend(list(np.eye(n_rxn)))
					self.generator_list_A.extend(list(-1*np.eye(n_rxn)))
					self.generator_list_A.append(np.zeros(n_rxn))
					self.generator_list_A.append(np.zeros(n_rxn))
					self.generator_list_A.append(np.zeros(n_rxn))
					
							#For zeta-b
							
					self.generator_list_B.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.4*n_b)])	
					self.generator_list_B.extend(list(np.eye(2*n_rxn)))
					self.generator_list_B.extend(list(-1*np.eye(2*n_rxn)))
					self.generator_list_B.append(np.zeros(2*n_rxn))
					self.generator_list_B.append(np.zeros(2*n_rxn))
					self.generator_list_B.append(np.zeros(2*n_rxn))
							
							#For zeta-C		
					self.generator_list_C.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.4*n_c)])
					self.generator_list_C.extend(list(np.eye(2*n_rxn)))
					self.generator_list_C.extend(list(-1*np.eye(2*n_rxn)))
					self.generator_list_C.append(np.zeros(2*n_rxn))
					self.generator_list_C.append(np.zeros(2*n_rxn))
					self.generator_list_C.append(np.zeros(2*n_rxn))
					
					
					self.generator_list_A_.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(n_a_)])
					#self.generator_list_A_.extend(list(np.eye(n_rxn)))
					#self.generator_list_A_.extend(list(-1*np.eye(n_rxn)))
					#self.generator_list_A_.append(np.zeros(n_rxn))
					#self.generator_list_A_.append(np.zeros(n_rxn))
					#self.generator_list_A_.append(np.zeros(n_rxn))
					
							#For zeta-A
							
					self.generator_list_B_.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(n_b_)])	
					#self.generator_list_B_.extend(list(np.eye(2*n_rxn)))
					#self.generator_list_B_.extend(list(-1*np.eye(2*n_rxn)))
					#self.generator_list_B_.append(np.zeros(2*n_rxn))
					#self.generator_list_B_.append(np.zeros(2*n_rxn))
					#self.generator_list_B_.append(np.zeros(2*n_rxn))
							
							#For zeta-A		
					self.generator_list_C_.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(n_c_)])
					#self.generator_list_C_.extend(list(np.eye(2*n_rxn)))
					#self.generator_list_C_.extend(list(-1*np.eye(2*n_rxn)))
					#self.generator_list_C_.append(np.zeros(2*n_rxn))
					#self.generator_list_C_.append(np.zeros(2*n_rxn))
					#self.generator_list_C_.append(np.zeros(2*n_rxn))
					
					
					
					"""
					For testing the distribution type
					"""
					"""		
							#For zeta-A
					
					self.generator_list_A.extend(list(2* np.random.random_sample((int(0.33*shuffled),n_rxn)) - 1))
					self.generator_list_A.extend(list(np.eye(n_rxn)))
					self.generator_list_A.extend(list(-1*np.eye(n_rxn)))
					self.generator_list_A.append(np.zeros(n_rxn))
					self.generator_list_A.append(np.zeros(n_rxn))
					self.generator_list_A.append(np.zeros(n_rxn))
					
							#For zeta-b
							
					self.generator_list_B.extend(list(2* np.random.random_sample((int(0.33*shuffled),2*n_rxn)) - 1))	
					self.generator_list_B.extend(list(np.eye(2*n_rxn)))
					self.generator_list_B.extend(list(-1*np.eye(2*n_rxn)))
					self.generator_list_B.append(np.zeros(2*n_rxn))
					self.generator_list_B.append(np.zeros(2*n_rxn))
					self.generator_list_B.append(np.zeros(2*n_rxn))
							
							#For zeta-C		
					self.generator_list_C.extend(list(2* np.random.random_sample((int(0.33*shuffled),2*n_rxn)) - 1))
					self.generator_list_C.extend(list(np.eye(2*n_rxn)))
					self.generator_list_C.extend(list(-1*np.eye(2*n_rxn)))
					self.generator_list_C.append(np.zeros(2*n_rxn))
					self.generator_list_C.append(np.zeros(2*n_rxn))
					self.generator_list_C.append(np.zeros(2*n_rxn))
					
					
					self.generator_list_A_.extend(list(2* np.random.random_sample((int(0.33*unshuffled),n_rxn)) - 1))
					#self.generator_list_A_.extend(list(np.eye(n_rxn)))
					#self.generator_list_A_.extend(list(-1*np.eye(n_rxn)))
					#self.generator_list_A_.append(np.zeros(n_rxn))
					#self.generator_list_A_.append(np.zeros(n_rxn))
					#self.generator_list_A_.append(np.zeros(n_rxn))
					
							#For zeta-A
							
					self.generator_list_B_.extend(list(2* np.random.random_sample((int(0.33*unshuffled),2*n_rxn)) - 1))	
					#self.generator_list_B_.extend(list(np.eye(2*n_rxn)))
					#self.generator_list_B_.extend(list(-1*np.eye(2*n_rxn)))
					#self.generator_list_B_.append(np.zeros(2*n_rxn))
					#self.generator_list_B_.append(np.zeros(2*n_rxn))
					#self.generator_list_B_.append(np.zeros(2*n_rxn))
							
							#For zeta-A		
					self.generator_list_C_.extend(list(2* np.random.random_sample((int(0.33*unshuffled),2*n_rxn)) - 1))
					#self.generator_list_C_.extend(list(np.eye(2*n_rxn)))
					#self.generator_list_C_.extend(list(-1*np.eye(2*n_rxn)))
					#self.generator_list_C_.append(np.zeros(2*n_rxn))
					#self.generator_list_C_.append(np.zeros(2*n_rxn))
					#self.generator_list_C_.append(np.zeros(2*n_rxn))
					"""
					"""
					Uniform
					
					"""
						
					#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
					
					"""
					Latin-Hypercube
					"""
							#For zeta-A
							
					lhs_a = []
					for i in lhs(n_rxn, samples=int(0.3*n_a)):
						lhs_a.append(np.asarray(i))	
					self.generator_list_A.extend(lhs_a)
					
							
							#For zeta-A
							
					lhs_b = []
					for i in lhs(2*n_rxn, samples=int(0.3*n_b)):
						lhs_b.append(np.asarray(i))		
					self.generator_list_B.extend(lhs_b)
							
							#For zeta-A
							
					lhs_c = []
					for i in lhs(2*n_rxn, samples=int(0.3*n_c)):
						lhs_c.append(np.asarray(i))
					self.generator_list_C.extend(lhs_c)
					
					
					"""		
					lhs_a = []
					for i in lhs(n_rxn, samples=int(0.33*self.sim_)):
						lhs_a.append(np.asarray(i))	
					self.generator_list_A.extend(lhs_a)
					
							
							#For zeta-A
							
					lhs_b = []
					for i in lhs(2*n_rxn, samples=int(0.33*self.sim_)):
						lhs_b.append(np.asarray(i))		
					self.generator_list_B.extend(lhs_b)
							
							#For zeta-A
							
					lhs_c = []
					for i in lhs(2*n_rxn, samples=int(0.33*self.sim_)):
						lhs_c.append(np.asarray(i))
					self.generator_list_C.extend(lhs_c)
					
					"""
					"""
					lhs_a_ = []
					for i in lhs(n_rxn, samples=int(0.33*unshuffled)):
						lhs_a_.append(np.asarray(i))	
					self.generator_list_A_.extend(lhs_a_)
					
							
							#For zeta-A
							
					lhs_b_ = []
					for i in lhs(2*n_rxn, samples=int(0.33*unshuffled)):
						lhs_b_.append(np.asarray(i))		
					self.generator_list_B_.extend(lhs_b_)
							
							#For zeta-A
							
					lhs_c_ = []
					for i in lhs(2*n_rxn, samples=int(0.33*unshuffled)):
						lhs_c_.append(np.asarray(i))
					self.generator_list_C_.extend(lhs_c_)
					
					"""
					"""
					Monte-Carlo-Mix
					"""
					#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)		
					#self.generator_list.extend(bb)
					#self.generator_list.extend(list(np.eye(self.n)))
					#self.generator_list.extend(list(-1*np.eye(self.n)))
					#self.generator_list.append(np.zeros(self.n))
					#self.generator_list.append(np.zeros(self.n))
					#self.generator_list.append(np.zeros(self.n))
					#print(self.generator_list)
					#self.sim_ = len(self.generator_list)
					
					#trainers = np.asarray(trainers)
					#print(f"Generator_list {np.shape(self.generator_list)}\n")
					df = pd.DataFrame(np.asarray(self.generator_list))
					count_start_a = 0
					count_start_b = 0
					count_start_c = 0			
					for i in self.reactions_selected:
						
						len_of_Arrhenius_A = 1
						count_end_a = count_start_a+len_of_Arrhenius_A
						
						
						len_of_Arrhenius_B = 2
						count_end_b = count_start_b+len_of_Arrhenius_B
						
						
						len_of_Arrhenius_C = 2
						count_end_c = count_start_c+len_of_Arrhenius_C
						
						#print(count_start)
						#print(count_end)
						self.rxn_generators_a2 = []
						self.rxn_generators_b2 = []
						self.rxn_generators_c2 = []
						
						self.rxn_generators_a2_ = []
						self.rxn_generators_b2_ = []
						self.rxn_generators_c2_ = []
						
						"""
						getting zetas from generators
						"""
						for j in range(len(self.generator_list_A)):
							self.rxn_generators_a2.append(self.generator_list_A[j][count_start_a:count_end_a])
						#print(self.rxn_generators_a2)
						for j in range(len(self.generator_list_B)):
							self.rxn_generators_b2.append(self.generator_list_B[j][count_start_b:count_end_b])
						
						for j in range(len(self.generator_list_C)):
							self.rxn_generators_c2.append(self.generator_list_C[j][count_start_c:count_end_c])
						
						
						for j in range(len(self.generator_list_A_)):
							self.rxn_generators_a2_.append(self.generator_list_A_[j][count_start_a:count_end_a])
						#print(self.rxn_generators_a2)
						for j in range(len(self.generator_list_B_)):
							self.rxn_generators_b2_.append(self.generator_list_B_[j][count_start_b:count_end_b])
						
						for j in range(len(self.generator_list_C_)):
							self.rxn_generators_c2_.append(self.generator_list_C_[j][count_start_c:count_end_c])
						
						
						
						
						self.sim_ = len(self.rxn_generators_a2)+len(self.rxn_generators_b2)+len(self.rxn_generators_c2)+len(self.rxn_generators_a2_)+len(self.rxn_generators_b2_)+len(self.rxn_generators_c2_)
						#print(self.rxn_generators)
						#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
											
						g_a = self.rxn_generators_a2
						g_b = self.rxn_generators_b2
						g_c = self.rxn_generators_c2
						
						g_a_ = self.rxn_generators_a2_
						g_b_ = self.rxn_generators_b2_
						g_c_ = self.rxn_generators_c2_
						
						#g_b = self.rxn_generators[0:n_b]
						#g_c = self.rxn_generators[n_b:n_b+n_c]
						
						
						count_start_a = count_end_a
						count_start_b = count_end_b
						count_start_c = count_end_c
						
						"""
						Getting the uncertainty data
						"""
						#fig = plt.figure()
						data = self.unsrt[i].data
						T = np.linspace(300,2500,8)
						kappa_max = self.unsrt[i].getKappaMax(T)			
						kappa_min = self.unsrt[i].getKappaMin(T)
						Theta = np.array([T/T,np.log(T),-1/T])
						P = self.unsrt[i].getMean()	
						cov = self.unsrt[i].getCov()
						zeta_A = self.unsrt[i].zeta.x
						
						#plt.plot(1/T,kappa_max,"k--")
						#plt.plot(1/T,kappa_min,"k--")
						"""
						Type-A zeta samples
						"""
						generators_a = g_a
						generators_a_ = g_a_
						parallel_zetas_a = []
						parallel_zetas_a_ = []
						
						time_start = time.time()
						
						
						for k in g_a:
							
							parallel_zetas_a.append(k[0]*zeta_A)
						
						for k in g_a_:
							
							parallel_zetas_a_.append(k[0]*zeta_A)
											
						
						generators_b_ = []
						parallel_zetas_b_ = []
						generators_c_ = []
						parallel_zetas_c_ = []
						
						generators_b = []
						parallel_zetas_b = []
						generators_c = []
						parallel_zetas_c = []
						start = 0
						start_ = 0
						for w in range(8):
							stop = int(((w+1)/8)*shuffled)
							stop_ = int(((w+1)/8)*unshuffled)
							g_b = self.rxn_generators_b2[start:stop]
							g_c = self.rxn_generators_c2[start:stop]
							g_b_ = self.rxn_generators_b2_[start_:stop_]
							g_c_ = self.rxn_generators_c2_[start_:stop_]
							data["generators_b"] = g_b
							callWorkForce = Worker(self.allowed_count)	
							gb,pb = callWorkForce.do_unsrt_b(data,len(g_b))
							del callWorkForce
							
							generators_b.extend(gb)
							parallel_zetas_b.extend(pb)
							
							removed = data.pop("generators_b")
							data["generators_b"] = g_b_
							callWorkForce = Worker(self.allowed_count)	
							gb,pb = callWorkForce.do_unsrt_b(data,len(g_b_))
							del callWorkForce
							
							generators_b_.extend(gb)
							parallel_zetas_b_.extend(pb)
							
							removed = data.pop("generators_b")
							data["generators_c"] = g_c
							callWorkForce = Worker(self.allowed_count)	
							gc,pc = callWorkForce.do_unsrt_c(data,len(g_c))
							del callWorkForce
							
							generators_c.extend(gc)
							parallel_zetas_c.extend(pc)
							removed = data.pop("generators_c")
							
							data["generators_c"] = g_c_
							callWorkForce = Worker(self.allowed_count)	
							gc,pc = callWorkForce.do_unsrt_c(data,len(g_c_))
							del callWorkForce
							
							generators_c_.extend(gc)
							parallel_zetas_c_.extend(pc)
							removed = data.pop("generators_c")
							
							start = stop
							start_ = stop_
						
						#print(self.sim_)
						#print(len(parallel_zetas_b)+len(parallel_zetas_c))
						g_b = self.rxn_generators_b2
						g_c = self.rxn_generators_c2
						
						g_b_ = self.rxn_generators_b2_
						g_c_ = self.rxn_generators_c2_
						
						
						temp_a[i] = parallel_zetas_a
						gen_a[i] = generators_a
						temp_b[i] = parallel_zetas_b
						gen_b[i] = generators_b
						temp_c[i] = parallel_zetas_c
						gen_c[i] = generators_c
						
						temp_a_[i] = parallel_zetas_a_
						gen_a_[i] = generators_a_
						temp_b_[i] = parallel_zetas_b_
						gen_b_[i] = generators_b_
						temp_c_[i] = parallel_zetas_c_
						gen_c_[i] = generators_c_
						
						time_end = time.time()
						print(f"\n\t\tTime taken = {time_end-time_start}\n")
						
						#print(parallel_zetas_b)
						Y_a = parallel_zetas_a
						Y_b = parallel_zetas_b
						Y_c = parallel_zetas_c
						
						Y_a_ = parallel_zetas_a_
						Y_b_ = parallel_zetas_b_
						Y_c_ = parallel_zetas_c_
						
						#print(self.sim_)
						#print(len(Y_a))
						#print(len(Y_b))
						#print(len(Y_c))
						#print(len(Y_a_))
						#print(len(Y_b_))
						#print(len(Y_c_))
						
						
						T_ = np.linspace(300,2500,8)
						theta = np.array([T/T,np.log(T),-1/T])
						kappa_mean = theta.T.dot(P)
						Theta = np.array([T/T,np.log(T),-1/T])
						X = []
						Y = []
						
						X_ = []
						Y_ = []
						
						
						for j in Y_a:
							Y.append(j)
							Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
							
							X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
							
						for j in Y_a_:
							Y_.append(j)
							Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
							
							X_.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
							
						
						for j in Y_b:
							Y.append(j)
							Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
							#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
							X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						
						for j in Y_b_:
							Y_.append(j)
							Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
							#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
							X_.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
							#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
						
						for j in Y_c:
							Y.append(j)
							Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
							#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
							X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
							#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
						#plt.savefig(str(i)+"_c.pdf")
						
						for j in Y_c_:
							Y_.append(j)
							Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
							#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
							X_.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
							#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
						#plt.savefig(str(i)+"_c.pdf")
						
						temp_xdata = {}
						for k in range(len(Y_)):
							temp_xdata[str(Y_[k])] = X_[k]
							
						random.shuffle(Y_)
						#print(Y_)
						X1_ = []
						for k in range(len(Y_)):
							X1_.append(temp_xdata[str(Y_[k])])
							
						#print(X_[9])
						xdata[i] = X
						ydata[i] = Y
						
						xdata_[i] = X1_
						ydata_[i] = Y_
						
						#print(len(X1_))
						#print(len(Y_))
						#print(len(X))
						#print(len(Y))
						print(self.sim_)
						#Was trying to shuffle the zetas but the try was not successfull
						"""
						Y_ = []
						a = int(unshuffled/3)
						b = unshuffled - 2*a
						for w in Y_a[0:a]:
							Y_.append(w)
						for w in Y_b[0:a]:
							Y_.append(w)
						for w in Y_c[0:b]:
							Y_.append(w)
						
						random.shuffle(Y_)
						ydata_[i]=Y_
						
						#ydata[i] = Y
						#print(len(Y))
						#print(len(Y_a)+len(Y_b)+len(Y_c)))
						extra = len(Y_a[0:a])+len(Y_b[0:a])+len(Y_c[0:b])
						#print(len(parallel_zetas))
						#print(len(self.rxn_generators))
						"""
					A = []
					A_gen = []
					X_dash = []
					y_dash = []
					
					for j in range(len(g_a)):
						a_row = []
						a_row_gen = []
						x_temp = []
						y_temp = []
						for i in temp_a:
							a_row.extend(list(temp_a[i][j]))
							a_row_gen.extend(list(gen_a[i][j]))
							x_temp.extend(list(xdata[i][j]))
							y_temp.extend(list(ydata[i][j]))
						A.append(np.asarray(a_row))
						A_gen.append(np.asarray(a_row_gen))
						X_dash.append(np.asarray(x_temp))
						y_dash.append(np.asarray(y_temp))
						
					for j in range(len(g_a_)):
						a_row = []
						a_row_gen = []
						x_temp = []
						y_temp = []
						for i in temp_a_:
							a_row.extend(list(ydata_[i][j]))
							a_row_gen.extend(list(gen_a_[i][j]))
							x_temp.extend(list(xdata_[i][j]))
							y_temp.extend(list(ydata_[i][j]))
						A.append(np.asarray(a_row))
						A_gen.append(np.asarray(a_row_gen))
						X_dash.append(np.asarray(x_temp))
						y_dash.append(np.asarray(y_temp))
					
					
					for j in range(len(g_b)):
						a_row = []
						a_row_gen = []
						x_temp = []
						y_temp = []
						for i in temp_b:
							a_row.extend(list(temp_b[i][j]))
							a_row_gen.extend(list(gen_b[i][j]))
							x_temp.extend(list(xdata[i][j+len(g_a)]))
							y_temp.extend(list(ydata[i][j+len(g_a)]))
						A.append(np.asarray(a_row))
						A_gen.append(np.asarray(a_row_gen))
						X_dash.append(np.asarray(x_temp))
						y_dash.append(np.asarray(y_temp))
					
					for j in range(len(g_b_)):
						a_row = []
						a_row_gen = []
						x_temp = []
						y_temp = []
						for i in temp_b_:
							a_row.extend(list(ydata_[i][j+len(g_a_)]))
							a_row_gen.extend(list(gen_b_[i][j]))
							x_temp.extend(list(xdata_[i][j+len(g_a_)]))
							y_temp.extend(list(ydata_[i][j+len(g_a_)]))
						A.append(np.asarray(a_row))
						A_gen.append(np.asarray(a_row_gen))
						X_dash.append(np.asarray(x_temp))
						y_dash.append(np.asarray(y_temp))
					"""
					for j in range(n_b):
						a_row = []
						a_row_gen = []
						x_temp = []
						y_temp = []
						for i in temp_b:
							a_row.extend(list(temp_b[i][j]))
							a_row_gen.extend(list(gen_b[i][j]))
							x_temp.extend(list(xdata[i][j]))
							y_temp.extend(list(ydata[i][j]))
						A.append(np.asarray(a_row))
						A_gen.append(np.asarray(a_row_gen))
						X_dash.append(np.asarray(x_temp))
						y_dash.append(np.asarray(y_temp))
					"""
					
					for j in range(len(g_c)):
						a_row = []
						a_row_gen = []
						x_temp = []
						y_temp = []
						for i in temp_c:
							a_row.extend(list(temp_c[i][j]))
							a_row_gen.extend(list(gen_c[i][j]))
							x_temp.extend(list(xdata[i][j+len(g_a)+len(g_b)]))
							y_temp.extend(list(ydata[i][j+len(g_a)+len(g_b)]))
						A.append(np.asarray(a_row))
						A_gen.append(np.asarray(a_row_gen))
						X_dash.append(np.asarray(x_temp))
						y_dash.append(np.asarray(y_temp))
					
					#print(len())
					for j in range(len(g_c_)):
						a_row = []
						a_row_gen = []
						x_temp = []
						y_temp = []
						for i in temp_c_:
							#print(len(ydata_[i]))
							#print(len(g_a_)+len(g_b_)+len(g_c_))
							a_row.extend(list(ydata_[i][j+len(g_a_)+len(g_b_)]))
							a_row_gen.extend(list(gen_c_[i][j]))
							x_temp.extend(list(xdata_[i][j+len(g_a_)+len(g_b_)]))
							y_temp.extend(list(ydata_[i][j+len(g_a_)+len(g_b_)]))
						A.append(np.asarray(a_row))
						A_gen.append(np.asarray(a_row_gen))
						X_dash.append(np.asarray(x_temp))
						y_dash.append(np.asarray(y_temp))
					
					
					
					"""
					for j in range(extra):
						a_row = []
						a_row_gen = []
						x_temp = []
						y_temp = []
						for i in temp_a:
							a_row.extend(list(ydata_[i][j]))
							y_temp.extend(list(ydata_[i][j]))
						A.append(np.asarray(a_row))
						y_dash.append(np.asarray(y_temp))
					"""
					"""
					for j in range(n_c):
						a_row = []
						a_row_gen = []
						x_temp = []
						y_temp = []
						for i in temp_c:
							a_row.extend(list(temp_c[i][j]))
							a_row_gen.extend(list(gen_c[i][j]))
							x_temp.extend(list(xdata[i][j+n_b]))
							y_temp.extend(list(ydata[i][j+n_b]))
						A.append(np.asarray(a_row))
						A_gen.append(np.asarray(a_row_gen))
						X_dash.append(np.asarray(x_temp))
						y_dash.append(np.asarray(y_temp))
						
					"""
					
					#print(A)
					#self.sim_ = len(g_a)+len(g_b)+len(g_c)
					#self.sim_ = n_b+n_c
					self.beta_ = np.asarray(A)
					#print(self.beta_)
					self.generators = np.asarray(A_gen)
					self.X = np.asarray(X_dash)
					self.y = np.asarray(y_dash)
			
			elif self.design == "B1-MC":
				self.beta_=[]
				temp_a = {}
				temp_b = {}
				temp_c = {}
				gen = {}
				gen_a = {}
				gen_b = {}
				gen_c = {}
				
				xdata = {}
				ydata = {}				
				self.generator_list = []
				"""
				Min- simulation = 300
				"""
				#if self.sim_ < 300:
				#	self.sim_ = 300
				#print(self.sim_)
				"""
				Dividing the generators into three bins
				------------------------------------
				n_a = 0.33*self.sim_
				n_b = 0.33*self.sim_
				n_c = 0.34*self.sim_
				
				"""
				n_a = int(0.33*self.sim_)
				n_b = int(0.33*self.sim_)
				n_c = self.sim_-n_a-n_b
				
				
				#n_b = int(0.5*self.sim_)
				#n_c = self.sim_-n_b
				
				
				#for i in np.eye(self.n):
				#	self.generator_list.append(i)
				#	self.generator_list.append(-1*i)
								
				"""
				Mixing three types of the extreme curves
				--------------------
				zeta-A: 1
				zeta-B: 1
				zeta-C: 1
			
				for zeta-A: one normal random parameter is required
					a matrix on n_a x n will be created
					
				for zeta-B and zeta-c: two random parameters is required
					we can arrange the two normal_random_parameter 				
					a matrix of n_b+n_c x 2*n will be created
					
				The number of reactions = n_R
				"""
				n_rxn = len(self.ind)
				self.generator_list_A = []
				self.generator_list_B = []
				self.generator_list_C = []
				
				"""
				Factorial Design
				----------------
				box-benkhen design
				bb = []
				for i in bbdesign(n_rxn)[0:int(0.3*self.sim_)]:
					bb.append(np.asarray(i))
				"""
				"""		#For zeta-A
						
				bb_a = []
				for i in bbdesign(n_rxn)[0:int(0.4*n_a)]:
					bb_a.append(np.asarray(i))
				self.generator_list_A.extend(bb_a)
						
						#For zeta-B
						
				bb_b = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_b)]:
					bb_b.append(np.asarray(i))							
				self.generator_list_B.extend(bb_b)
						
						#For zeta-C
						
				bb_c = []
				for i in bbdesign(2*n_rxn)[0:int(0.4*n_c)]:
					bb_c.append(np.asarray(i))		
						
				self.generator_list_C.extend(bb_c)		
						
				"""
				"""
				Monte-Carlo
				
				"""
						
						#For zeta-A
				
				self.generator_list_A.extend(list(2* np.random.random_sample((self.sim_,n_rxn)) - 1)[0:int(0.3*n_a)])
				#self.generator_list_A.extend(list(np.eye(n_rxn)))
				#self.generator_list_A.extend(list(-1*np.eye(n_rxn)))
				#self.generator_list_A.append(np.zeros(n_rxn))
				#self.generator_list_A.append(np.zeros(n_rxn))
				#self.generator_list_A.append(np.zeros(n_rxn))
				
						#For zeta-A
						
				self.generator_list_B.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_a)])	
				#self.generator_list_B.extend(list(np.eye(2*n_rxn)))
				#self.generator_list_B.extend(list(-1*np.eye(2*n_rxn)))
				#self.generator_list_B.append(np.zeros(2*n_rxn))
				#self.generator_list_B.append(np.zeros(2*n_rxn))
				#self.generator_list_B.append(np.zeros(2*n_rxn))
						
						#For zeta-A		
				self.generator_list_C.extend(list(2* np.random.random_sample((self.sim_,2*n_rxn)) - 1)[0:int(0.3*n_a)])
				#self.generator_list_C.extend(list(np.eye(2*n_rxn)))
				#self.generator_list_C.extend(list(-1*np.eye(2*n_rxn)))
				#self.generator_list_C.append(np.zeros(2*n_rxn))
				#self.generator_list_C.append(np.zeros(2*n_rxn))
				#self.generator_list_C.append(np.zeros(2*n_rxn))
				
				"""
				Uniform
				
				"""
					
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				
				"""
				Latin-Hypercube
				"""
				"""
						
						#For zeta-A
						
				lhs_a = []
				for i in lhs(n_rxn, samples=int(0.3*n_a)):
					lhs_a.append(np.asarray(i))	
				self.generator_list_A.extend(lhs_a)
				
						
						#For zeta-A
						
				lhs_b = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_b)):
					lhs_b.append(np.asarray(i))		
				self.generator_list_B.extend(lhs_b)
						
						#For zeta-A
						
				lhs_c = []
				for i in lhs(2*n_rxn, samples=int(0.3*n_c)):
					lhs_c.append(np.asarray(i))
				self.generator_list_C.extend(lhs_c)
				
				"""
				"""
				Monte-Carlo-Mix
				"""
				
				
				
				
				
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)		
				#self.generator_list.extend(bb)
				#self.generator_list.extend(list(np.eye(self.n)))
				#self.generator_list.extend(list(-1*np.eye(self.n)))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#print(self.generator_list)
				#self.sim_ = len(self.generator_list)
				
				#trainers = np.asarray(trainers)
				#print(f"Generator_list {np.shape(self.generator_list)}\n")
				df = pd.DataFrame(np.asarray(self.generator_list))
				count_start_a = 0
				count_start_b = 0
				count_start_c = 0			
				for i in self.reactions_selected:
					
					len_of_Arrhenius_A = 1
					count_end_a = count_start_a+len_of_Arrhenius_A
					
					
					len_of_Arrhenius_B = 2
					count_end_b = count_start_b+len_of_Arrhenius_B
					
					
					len_of_Arrhenius_C = 2
					count_end_c = count_start_c+len_of_Arrhenius_C
					
					#print(count_start)
					#print(count_end)
					self.rxn_generators_a2 = []
					self.rxn_generators_b2 = []
					self.rxn_generators_c2 = []
					
					"""
					getting zetas from generators
					"""
					for j in range(len(self.generator_list_A)):
						self.rxn_generators_a2.append(self.generator_list_A[j][count_start_a:count_end_a])
					#print(self.rxn_generators_a2)
					for j in range(len(self.generator_list_B)):
						self.rxn_generators_b2.append(self.generator_list_B[j][count_start_b:count_end_b])
					
					for j in range(len(self.generator_list_C)):
						self.rxn_generators_c2.append(self.generator_list_C[j][count_start_c:count_end_c])
					
					self.sim_ = len(self.rxn_generators_a2)+len(self.rxn_generators_b2)+len(self.rxn_generators_c2)
					#print(self.rxn_generators)
					#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
										
					g_a = self.rxn_generators_a2
					g_b = self.rxn_generators_b2
					g_c = self.rxn_generators_c2
					#g_b = self.rxn_generators[0:n_b]
					#g_c = self.rxn_generators[n_b:n_b+n_c]
					
					
					count_start_a = count_end_a
					count_start_b = count_end_b
					count_start_c = count_end_c
					
					"""
					Getting the uncertainty data
					"""
					#fig = plt.figure()
					data = self.unsrt[i].data
					T = np.linspace(300,2500,8)
					kappa_max = self.unsrt[i].getKappaMax(T)			
					kappa_min = self.unsrt[i].getKappaMin(T)
					Theta = np.array([T/T,np.log(T),-1/T])
					P = self.unsrt[i].getMean()	
					cov = self.unsrt[i].getCov()
					zeta_A = self.unsrt[i].zeta.x
					
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					"""
					Type-A zeta samples
					"""
					generators_a = g_a
					parallel_zetas_a = []
					
					time_start = time.time()
					
					
					for k in g_a:
						#Pint = P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten()
					#	plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
					#	plt.plot()
						#print(i)
						#print(type(i))
						#print(k)	
						#parallel_zetas_a.append(P+np.asarray(np.dot(cov,np.array(k)*zeta_A)).flatten())
						parallel_zetas_a.append(k[0]*zeta_A)
										
					#print(parallel_zetas_a)
					#plt.savefig(str(i)+"_a.pdf")
					#print(os.getcwd())
					#print("plot saved")
					data["generators_b"] = g_b
					callWorkForce = Worker(self.allowed_count)	
					generators_b,parallel_zetas_b = callWorkForce.do_unsrt_b(data,len(g_b))
					del callWorkForce
					
					data["generators_c"] = g_c
					callWorkForce = Worker(self.allowed_count)	
					generators_c,parallel_zetas_c = callWorkForce.do_unsrt_c(data,len(g_c))
					del callWorkForce
					
					#print(parallel_zetas_a)
					
					temp_a[i] = parallel_zetas_a
					gen_a[i] = generators_a
					temp_b[i] = parallel_zetas_b
					gen_b[i] = generators_b
					temp_c[i] = parallel_zetas_c
					gen_c[i] = generators_c
					
					time_end = time.time()
					print(f"\n\t\tTime taken = {time_end-time_start}\n")
					
					#print(parallel_zetas_b)
					Y_a = parallel_zetas_a
					Y_b = parallel_zetas_b
					Y_c = parallel_zetas_c
					T_ = np.linspace(300,2500,8)
					theta = np.array([T/T,np.log(T),-1/T])
					kappa_mean = theta.T.dot(P)
					Theta = np.array([T/T,np.log(T),-1/T])
					X = []
					Y = []
					
					
					for j in Y_a:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					
					for j in Y_b:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_b.pdf")
					
					
					#fig = plt.figure()
					#plt.plot(1/T,kappa_max,"k--")
					#plt.plot(1/T,kappa_min,"k--")
					for j in Y_c:
						Y.append(j)
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						#plt.plot(1/T,np.asarray(Theta.T.dot(Pint)).flatten(),"r--",linewidth=0.15)
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					#plt.savefig(str(i)+"_c.pdf")
					
					temp_xdata = {}
					for k in range(len(Y)):
						temp_xdata[str(Y[k])] = X[k]
						
					random.shuffle(Y)
					#print(Y_)
					X_ = []
					for k in range(len(Y)):
						X_.append(temp_xdata[str(Y[k])])
						
					#print(X_[9])
					xdata[i] = X_
					ydata[i] = Y
					#print(len(parallel_zetas))
					#print(len(self.rxn_generators))
				A = []
				A_gen = []
				X_dash = []
				y_dash = []
				
				for j in range(len(g_a)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_a:
						a_row.extend(list(temp_a[i][j]))
						a_row_gen.extend(list(gen_a[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				
				
				for j in range(len(g_b)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)]))
						y_temp.extend(list(ydata[i][j+len(g_a)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				for j in range(n_b):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_b:
						a_row.extend(list(temp_b[i][j]))
						a_row_gen.extend(list(gen_b[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				
				for j in range(len(g_c)):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+len(g_a)+len(g_b)]))
						y_temp.extend(list(ydata[i][j+len(g_a)+len(g_b)]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
				"""
				for j in range(n_c):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp_c:
						a_row.extend(list(temp_c[i][j]))
						a_row_gen.extend(list(gen_c[i][j]))
						x_temp.extend(list(xdata[i][j+n_b]))
						y_temp.extend(list(ydata[i][j+n_b]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X_dash.append(np.asarray(x_temp))
					y_dash.append(np.asarray(y_temp))
					
				"""
				
				#print(A)
				self.sim_ = len(g_a)+len(g_b)+len(g_c)
				#self.sim_ = n_b+n_c
				self.beta_ = np.asarray(A)
				#print(self.beta_)
				self.generators = np.asarray(A_gen)
				self.X = np.asarray(X_dash)
				self.y = np.asarray(y_dash)
			
			elif self.design == "B1-Linear":
				self.beta_=[]
				temp = {}
				gen = {}
				xdata = {}
				ydata = {}				
				self.generator_list = []
				#self.sim_ = int((25/9)*self.sim_)
				#self.generator_list = list(2* np.random.random_sample((self.sim_,self.n)) - 1)
				#for i in np.eye(self.n):
				#	self.generator_list.append(i)
				#	self.generator_list.append(-1*i)
								
				"""
				box-benkhen design
				"""
				bb = []
				for i in bbdesign(self.n)[0:int(0.3*self.sim_)]:
					bb.append(np.asarray(i))
					
				
				#self.generator_list.extend(bb)
				#self.generator_list.extend(list(np.eye(self.n)))
				#self.generator_list.extend(list(-1*np.eye(self.n)))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#self.generator_list.append(np.zeros(self.n))
				#print(self.generator_list)
				#self.sim_ = len(self.generator_list)
				
				#trainers = []
				#for i in self.generator_list:
				
				#	trainers.append(np.asarray(i))
				#trainers = np.asarray(trainers)
				#print(f"Generator_list {np.shape(self.generator_list)}\n")
				df = pd.DataFrame(np.asarray(self.generator_list))
				count_start = 0			
				for i in self.reactions_selected:
					len_of_Arrhenius = len(self.unsrt[i].activeParameters)
					count_end = count_start+len_of_Arrhenius
					bb = bbdesign(len_of_Arrhenius)
					#print(count_start)
					#print(count_end)
					self.rxn_generators = []
					
					"""
					Getting linear combination of zeta
					"""
					
					for j in range(self.sim_):
						beta1 = 2*np.random.random_sample(1)-1
						beta2 = 2*np.random.random_sample(1)-1
						m = beta1*np.pi
						n = beta2*np.pi/2
						alpha = m + np.pi
						gamma = n + (np.pi/2)
						x = np.cos(alpha)*np.sin(gamma)
						y = np.sin(alpha)*np.sin(gamma)
						z = np.cos(gamma)
						self.rxn_generators.append(np.array([x[0],y[0],z[0]]))
					
					
					"""
					getting zetas from generators
					"""
					#for j in range(len(self.generator_list)):
					#	self.rxn_generators.append(self.generator_list[j][count_start:count_end])
					"""
					for k in bb:
						self.rxn_generators.append(np.asarray(k))
					self.rxn_generators.append(np.array([1.,1.,1.]))
					self.rxn_generators.append(np.array([-1.,-1.,-1.]))
					self.rxn_generators.append(np.array([1.,0.,0.]))
					self.rxn_generators.append(np.array([0.,1.,0.]))
					self.rxn_generators.append(np.array([0.,0.,1.]))
					self.rxn_generators.append(np.array([0.,-1.,0.]))
					self.rxn_generators.append(np.array([-1.,0.,0.]))
					self.rxn_generators.append(np.array([0.,0.,-1.]))
					#self.rxn_generators.append(np.array([0.,0.,0.]))
					#self.rxn_generators.append(np.array([0.,0.,0.]))
					#self.rxn_generators.append(np.array([0.,0.,0.]))			
					self.sim_ = len(self.rxn_generators)
					#print(self.rxn_generators)
					#print(f"len of rxn_generator {np.shape(self.rxn_generators)}\n")
					"""
					count_start = count_end
					
					data = self.unsrt[i].data
					T = np.linspace(300,2500,8)
					kappa_max = self.unsrt[i].getKappaMax(T)			
					kappa_min = self.unsrt[i].getKappaMin(T)
					P = self.unsrt[i].getMean()	
					zetaMat = self.unsrt[i].zeta_Matrix
					
					parallel_zetas = []
					for k in range(self.sim_):
						zeta = np.asarray(zetaMat.T.dot(self.rxn_generators[k].T)).flatten()
						parallel_zetas.append(zeta)
					
					cov = self.unsrt[i].getCov()
					#data["generators"] = self.rxn_generators
					#callWorkForce = Worker(10)	
					#generators,parallel_zetas = callWorkForce.do_unsrt(data,self.sim_)
					#del callWorkForce
					temp[i] = parallel_zetas
					gen[i] = self.rxn_generators
					Y = parallel_zetas
					T_ = np.linspace(300,2500,8)
					theta = np.array([T/T,np.log(T),-1/T])
					kappa_mean = theta.T.dot(P)
					Theta = np.array([T_/T_,np.log(T_),-1/T_])
					X = []
					for j in Y:
						#print(i)
						#print(kappa_max,np.asarray(theta.T.dot(j)).flatten())
						#print(np.asarray(theta.T.dot(j)).flatten()/kappa_max)
						#Pint = P + np.asarray(np.dot(cov,np.asarray(j))).flatten();
						Pint = P + np.asarray(np.dot(cov,np.asarray(j).T)).flatten();
						X.append((np.asarray(theta.T.dot(Pint)).flatten()-kappa_mean)/(kappa_max-kappa_mean))
						#X.append(np.log(np.asarray(theta.T.dot(Pint)).flatten()/kappa_mean)/np.log(kappa_max/kappa_mean))
					temp_xdata = {}
					for k in range(len(Y)):
						temp_xdata[str(Y[k])] = X[k]
						
					random.shuffle(Y)
					#print(Y_)
					X_ = []
					for k in range(len(Y)):
						X_.append(temp_xdata[str(Y[k])])
						
					#print(X_[9])
					xdata[i] = X_
					ydata[i] = Y
					#print(len(parallel_zetas))
					#print(len(self.rxn_generators))
				A = []
				A_gen = []
				X = []
				y = []
				for j in range(self.sim_):
					a_row = []
					a_row_gen = []
					x_temp = []
					y_temp = []
					for i in temp:
						a_row.extend(list(temp[i][j]))
						a_row_gen.extend(list(gen[i][j]))
						x_temp.extend(list(xdata[i][j]))
						y_temp.extend(list(ydata[i][j]))
					A.append(np.asarray(a_row))
					A_gen.append(np.asarray(a_row_gen))
					X.append(np.asarray(x_temp))
					y.append(np.asarray(y_temp))
				self.beta_ = np.asarray(A)
				#print(self.beta_)
				self.generators = np.asarray(A_gen)
				self.X = np.asarray(X)
				self.y = np.asarray(y)
		
		elif self.simulation == "Original":
			self.sim_ = 1
			self.beta_ = []
			for i in range(self.sim_):
				temp = 0*(2*np.random.random(self.n)-1)
				self.beta_.append(temp)
			self.X = ""
		elif self.simulation == "Optimized":
			self.sim_ = 1
			self.beta_ = []
			for i in range(self.sim_):
				temp = 0*(2*np.random.random(self.n)-1)
				self.beta_.append(temp)
			self.X = ""
	def make_dir_in_parallel(self):
		#self.MechManipulator = MechanismManipulator.MechanismManipulator(self.mech_loc,self.fileType,self.thermo_file_location,self.trans_file_location,self.ind,self.foc,self.mBody,self.thm,self.trans,self.unsrt, self.focUnsert, self.tbdUnsert, self.thermoUnsert, self.transUnsert,self.selectedParams,self.activeIndexDict,design_type=self.design)
		originalMech = Parser(self.mech_loc).mech
		self.copy_of_mech = deepcopy(originalMech)#.deepcopy()
		
		#self.manipulator(copy_of_mech,)		
		self.R_list=self.ind
		self.foc_list = self.foc
		self.molecule = []
		for i in self.mBody:
			self.molecule.append(self.tbdUnsert[i].branches.split(","))
		
		self.thermo_list = self.thm
		self.transport_list = self.trans
		count_Arrhenius_params = 0
		for i in self.ind:
			if "all" in self.unsrt[i].perturbation_type:
				count_Arrhenius_params+=3
			elif "factor" in self.unsrt[i].perturbation_type:
				count_Arrhenius_params+=1
			else:
				raise AssertionError("Incorrect perturbation type")
				
		
		self.n = len(self.selectedParams)
		#self.n = len(self.ind)*8
		#print(f"selected params are:{self.selectedParams}\n\tlength is {len(self.selectedParams)}")
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
		if self.PRS_type == "Full":
			self.reactions_selected = self.ind
			self.getTotalUnknowns()
			self.getSamples()						
		
		elif self.PRS_type == "Partial" and self.simulation == "sa":
			self.reactions_selected = self.ind
			self.getTotalUnknowns()
			self.getSamples()		
		
		for case_index,case in enumerate(self.case_dir):
			if self.PRS_type == "Partial" and self.simulation != "sa":
				if self.activeIndexDict is not None:
					count_params = 0
					for i in self.activeIndexDict[case_index]:
						if self.activeIndexDict[case_index][i] == 1:
							count_params +=1
					self.n = count_params
					self.reactions_selected = self.activeReactionsDict[case]
					print(f"Reactions active are {self.reactions_selected}")
					print(f"Active params are: {self.n}\n")
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

	def make_directories_for_simulation(self):
		self.MechManipulator = MechanismManipulator.MechanismManipulator(self.mech_loc,self.fileType,self.thermo_file_location,self.trans_file_location,self.ind,self.foc,self.mBody,self.thm,self.trans,self.unsrt, self.focUnsert, self.tbdUnsert, self.thermoUnsert, self.transUnsert,self.selectedParams,self.activeIndexDict,design_type=self.design)
		self.R_list=self.ind
		self.foc_list = self.foc
		self.molecule = []
		for i in self.mBody:
			self.molecule.append(self.tbdUnsert[i].branches.split(","))
		
		self.thermo_list = self.thm
		self.transport_list = self.trans
		
		"""
		count_Arrhenius_params = 0
		for i in self.ind:
			if "all" in self.unsrt[i].perturbation_type:
				count_Arrhenius_params+=3
			elif "factor" in self.unsrt[i].perturbation_type:
				count_Arrhenius_params+=1
			else:
				raise AssertionError("Incorrect perturbation type")		
		self.n = count_Arrhenius_params+len(self.foc)+len(np.asarray(self.molecule).flatten())+7*len(self.thm)+2*len(self.trans)
		"""
		self.n = len(self.manipulatorDict["activeParameters"])
		
		
		
		self.dir_list = []	
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

		#Initializing sampling size
	
		count = 0
		start_time = time.time()
		print("Creating Directories for simulations......\n\n\n This may take a while... Please be patient...\n\n ")
		#Parallel_jobs = multiprocessing.Pool(self.allowed_count-1)
		self.case_manipulation = {}
		
		for case in self.case_dir:
			#Parallel_jobs.apply_async(self.makeDir, args = (case,len(self.case_dir)), callback = self.log_result)
			progress,dict_ = self.makeDir(case,len(self.case_dir))
			self.case_manipulation[str(case)] = dict_
			s = ''
			for k in self.beta_list:
				for l in k:
					s +='{},'.format(l)
				s +='\n'
			ff = open('../Data/Simulations/Beta_list_case-'+str(case)+'.csv','w')
			ff.write(s)
			ff.close(); 
			os.chdir('..')
		#Parallel_jobs.close()
		#Parallel_jobs.join()
		
		self.params_manipulation = {}
		for j in self.activeParameters:
			self.case_database = {}
			self.param_case_manipulation = {}
			for case in self.case_dir:
				temp = []
				for i,dictionary in enumerate(self.sim_dict):
					temp.append(self.case_manipulation[str(case)][str(dictionary)][str(j)])
				self.case_database[str(case)] = temp
			self.params_manipulation[str(j)] = self.case_database
			
		#print(self.params_manipulation)
		
		return self.dir_list,self.params_manipulation
