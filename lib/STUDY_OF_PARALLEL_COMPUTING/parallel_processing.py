import os,time,sys
import multiprocessing
import concurrent.futures
import copy
import marshal
import subprocess
from tqdm import tqdm
try:
    import ruamel_yaml as yaml
except ImportError:
    from ruamel import yaml

#### Read
with open("mechanism.yaml","r") as my_file:
	#print(my_file)
	yaml_files = yaml.safe_load(my_file)

def make_input_files():
	cantera_file = """#!/usr/bin/python
from __future__ import division
from __future__ import print_function
import pandas as pd
import numpy as np
import time
import cantera as ct
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

def find_nearest(array, value):
    array = np.asarray(np.abs(np.log(array)))
    value = np.abs(np.log(value))
    idx = (np.abs(array - value)).argmin()
    return idx,array[idx]

def fast_nearest_interp(xi, x, y):
    # Shift x points to centers
    spacing = np.diff(x) / 2
    x = x + np.hstack([spacing, spacing[-1]])
    # Append the last point in y twice for ease of use
    y = np.hstack([y, y[-1]])
    return y[np.searchsorted(x, xi)]

gas = ct.Solution('mechanism.yaml')

reactorTemperature = 1293.30970887598 #Kelvin
reactorPressure = 1368225.0

gas.TPX = reactorTemperature, reactorPressure,{'NC7H16': 0.017857, 'O2': 0.19643, 'N2': 0.785713}
r = ct.IdealGasReactor(contents=gas, name="Batch Reactor")
reactorNetwork = ct.ReactorNet([r])

stateVariableNames = [r.component_name(item) for item in range(r.n_vars)]
timeHistory = pd.DataFrame(columns=stateVariableNames)
#index_a = stateVariableNames.index("CH")

def ignitionDelay(df,pList, species,cond="max",specified_value ="None;None",exp_conc = "None"):
	if cond == "max":
		tau = df[species].idxmax()
	elif cond == "onset" and species != "p":
		time = df.index.to_numpy()
		conc = (df[species]).to_numpy()
		dtime = np.diff(df.index.to_numpy())
		dconc = np.diff((df[species]).to_numpy())
		slope = dconc/dtime
		intercept = int(np.diff(slope).argmax())
		tau = df.index.to_numpy()[intercept]
	elif cond == "onset" and species == "p":
		time = df.index.to_numpy()
		conc = np.asarray(pList)
		dtime = np.diff(df.index.to_numpy())
		dconc = np.diff(conc)
		slope = dconc/dtime
		intercept = int(np.diff(slope).argmax())
		tau = df.index.to_numpy()[intercept]
	elif cond == "dt-max" and species == "p":
		time = np.diff(df.index.to_numpy())
		conc = np.diff(np.asarray(pList))
		slope = conc/time
		index = int(slope.argmax())
		tau = df.index.to_numpy()[index]
	elif cond == "dt-max" and species != "p":
		time = np.diff(df.index.to_numpy())
		conc = np.diff(df[species].to_numpy())
		slope = conc/time
		index = int(slope.argmax())
		tau = df.index.to_numpy()[index]
	elif cond == "specific":
		if specified_value.split(";")[0] == None:
			raise Assertionerror("Input required for specified_value in ignition delay")
		else:
			target = float(specified_value.split(";")[0])
			unit = specified_value.split(";")[1]
			if unit == "molecule":
				avogadro = 6.02214E+23
				target = target/avogadro
				molecular_wt = gas.atomic_weight(species)
				target = molecular_wt*target	
				index,value = find_nearest(df[species],target)
				tau = df.index.to_numpy()[index[0]]	
			
			elif unit == "molecule/cm3":
				unit_conversion = 10**(6)
				time = df.index.to_numpy()
				conc = df[species].to_numpy()
				avogadro = 6.02214E+23
				target = (target/avogadro)*unit_conversion
				if exp_conc != "":
					exp_conc = (exp_conc/avogadro)*unit_conversion
					conc = conc*exp_conc
			
				f = CubicSpline(time,conc, bc_type='natural')
				time_new = np.arange(min(time),max(time),10**(-8))			
				conc_new = f(time_new)
				tau = fast_nearest_interp(target,conc_new,time_new)
				
			elif unit == "mole/cm3":
				time = df.index.to_numpy()
				conc = df[species].to_numpy()
				if exp_conc != "":
					conc = conc*exp_conc
				
				f = CubicSpline(time,conc, bc_type='natural')
				time_new = np.arange(min(time),max(time),10**(-8))
				conc_new = f(time_new)
				tau = fast_nearest_interp(target,conc_new,time_new)
	return tau

t0 = time.time()
estimatedIgnitionDelayTime = 2.0
t = 0
pressureList = []
counter = 1;
while(t < estimatedIgnitionDelayTime):
    t = reactorNetwork.step()
    #print(t)
    if (counter%1 == 0):
        #if reactorNetwork.get_state().astype('float64')[index_a] <0:
         #   raise AssertionError("Numerical error!!")        	
        timeHistory.loc[t] = reactorNetwork.get_state().astype('float64')
        pressureList.append(gas.P)
    counter+=1

tau = ignitionDelay(timeHistory,pressureList,"CH","max","None",)
# Toc
t1 = time.time()
tau_file = open("output/tau.out",'w')
tau_file.write("#T(K)    tau(us)"+"\\n"+"{}    {}".format(reactorTemperature,tau*(10**6)))
tau_file.close()
#timeHistory.to_csv("time_history.csv")
	"""
	extract_file = """
	"""
	flame_master_file = """#!/usr/bin/python
from __future__ import division
from __future__ import print_function
import pandas as pd
import numpy as np
import time
import cantera as ct
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

def find_nearest(array, value):
    array = np.asarray(np.abs(np.log(array)))
    value = np.abs(np.log(value))
    idx = (np.abs(array - value)).argmin()
    return idx,array[idx]

def fast_nearest_interp(xi, x, y):
    # Shift x points to centers
    spacing = np.diff(x) / 2
    x = x + np.hstack([spacing, spacing[-1]])
    # Append the last point in y twice for ease of use
    y = np.hstack([y, y[-1]])
    return y[np.searchsorted(x, xi)]

gas = ct.Solution('mechanism.yaml')

reactorTemperature = 1293.30970887598 #Kelvin
reactorPressure = 1368225.0

gas.TPX = reactorTemperature, reactorPressure,{'NC7H16': 0.017857, 'O2': 0.19643, 'N2': 0.785713}
r = ct.IdealGasReactor(contents=gas, name="Batch Reactor")
reactorNetwork = ct.ReactorNet([r])

stateVariableNames = [r.component_name(item) for item in range(r.n_vars)]
timeHistory = pd.DataFrame(columns=stateVariableNames)
#index_a = stateVariableNames.index("CH")

def ignitionDelay(df,pList, species,cond="max",specified_value ="None;None",exp_conc = "None"):
	if cond == "max":
		tau = df[species].idxmax()
	elif cond == "onset" and species != "p":
		time = df.index.to_numpy()
		conc = (df[species]).to_numpy()
		dtime = np.diff(df.index.to_numpy())
		dconc = np.diff((df[species]).to_numpy())
		slope = dconc/dtime
		intercept = int(np.diff(slope).argmax())
		tau = df.index.to_numpy()[intercept]
	elif cond == "onset" and species == "p":
		time = df.index.to_numpy()
		conc = np.asarray(pList)
		dtime = np.diff(df.index.to_numpy())
		dconc = np.diff(conc)
		slope = dconc/dtime
		intercept = int(np.diff(slope).argmax())
		tau = df.index.to_numpy()[intercept]
	elif cond == "dt-max" and species == "p":
		time = np.diff(df.index.to_numpy())
		conc = np.diff(np.asarray(pList))
		slope = conc/time
		index = int(slope.argmax())
		tau = df.index.to_numpy()[index]
	elif cond == "dt-max" and species != "p":
		time = np.diff(df.index.to_numpy())
		conc = np.diff(df[species].to_numpy())
		slope = conc/time
		index = int(slope.argmax())
		tau = df.index.to_numpy()[index]
	elif cond == "specific":
		if specified_value.split(";")[0] == None:
			raise Assertionerror("Input required for specified_value in ignition delay")
		else:
			target = float(specified_value.split(";")[0])
			unit = specified_value.split(";")[1]
			if unit == "molecule":
				avogadro = 6.02214E+23
				target = target/avogadro
				molecular_wt = gas.atomic_weight(species)
				target = molecular_wt*target	
				index,value = find_nearest(df[species],target)
				tau = df.index.to_numpy()[index[0]]	
			
			elif unit == "molecule/cm3":
				unit_conversion = 10**(6)
				time = df.index.to_numpy()
				conc = df[species].to_numpy()
				avogadro = 6.02214E+23
				target = (target/avogadro)*unit_conversion
				if exp_conc != "":
					exp_conc = (exp_conc/avogadro)*unit_conversion
					conc = conc*exp_conc
			
				f = CubicSpline(time,conc, bc_type='natural')
				time_new = np.arange(min(time),max(time),10**(-8))			
				conc_new = f(time_new)
				tau = fast_nearest_interp(target,conc_new,time_new)
				
			elif unit == "mole/cm3":
				time = df.index.to_numpy()
				conc = df[species].to_numpy()
				if exp_conc != "":
					conc = conc*exp_conc
				
				f = CubicSpline(time,conc, bc_type='natural')
				time_new = np.arange(min(time),max(time),10**(-8))
				conc_new = f(time_new)
				tau = fast_nearest_interp(target,conc_new,time_new)
	return tau

t0 = time.time()
estimatedIgnitionDelayTime = 2.0
t = 0
pressureList = []
counter = 1;
while(t < estimatedIgnitionDelayTime):
    t = reactorNetwork.step()
    #print(t)
    if (counter%1 == 0):
        #if reactorNetwork.get_state().astype('float64')[index_a] <0:
         #   raise AssertionError("Numerical error!!")        	
        timeHistory.loc[t] = reactorNetwork.get_state().astype('float64')
        pressureList.append(gas.P)
    counter+=1

tau = ignitionDelay(timeHistory,pressureList,"CH","max","None",)
# Toc
t1 = time.time()
tau_file = open("output/tau.out",'w')
tau_file.write("#T(K)    tau(us)"+"\\n"+"{}    {}".format(reactorTemperature,tau*(10**6)))
tau_file.close()
#timeHistory.to_csv("time_history.csv")
"""
	run_file = """#!/bin/bash
python3 cantera_.py &> solve
	"""
	return cantera_file,extract_file,flame_master_file,run_file

def callback_run(self, result):
	#print(result)
	self.progress.append(result[0])
	sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(self.progress)/float(result[-1])*100))
	sys.stdout.flush()

def run_generate_dir(location,total):
	os.mkdir(location)
	os.mkdir(location+"/output")
	return (location,total)
def run_executable_files_(args,total):
	os.chdir(args[0])
	subprocess.call(["./"+args[1]])
	return (args[0],total)
	
def run_map(params,total):
	location = str(params[2])
	sim1 = open(location+"/cantera_.py",'w').write(params[0])
	sim2= open(location+"/FlameMaster.input",'w').write(params[0])
	extract = open(location+"/extract.py",'w').write(params[3])
	#runConvertorScript = open(location+"/run_convertor",'w').write(params[2])
	runScript = open(location+"/run","w").write(params[1])
	#subprocess.call(["chmod","+x",location+"/run_convertor"])
	subprocess.call(["chmod","+x",location+"/run"])
	yaml_string = yaml.dump(params[-1],default_flow_style=False)
	with open(location+"/mechanism.yaml","+w") as yamlfile:
		yamlfile.write(yaml_string)
		yamlfile.close()
	del location
	del sim1, sim2,extract,runScript
	return (params[2],total)
		
class Worker():
	def __init__(self,workers):
		self.pool1 = concurrent.futures.ProcessPoolExecutor(max_workers=workers)
		self.pool = multiprocessing.Pool(processes=workers)
		self.progress = []
	def callback_error(self,result):
		print('error', result)
	def custom_error_callback(self,error):
   	 	print(f'Got an error: {error}')
   	 	
	def callback_run(self, result):
		#print(type(result[-1]))
		self.progress.append(result[0])
		
		sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(self.progress)/float(result[-1])*100))
		sys.stdout.flush()
	#def test_async(self,location_id):
	#	for                                 
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
			     args=(args,len(locations)),callback=self.callback_run)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
	
	def do_job_map_create(self,params):
		for param in params:
			self.pool.apply_async(run_map, 
			     args=(param,len(params)),callback=self.callback_run,error_callback=self.custom_error_callback)
	     	#self.pool.map_async(run_map,params,callback=self.callback_run,error_callback=self.custom_error_callback)
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

cantera_file = {}
exact_file = {}
flame_master_file = {}
dir_list = {}
yaml_mechanism = {}
run_file = {}
dir_list1 = []
os.mkdir("Test")
os.chdir("Test")
my_location = os.getcwd()

def deepcopy_v2(original):
    copied_dict = {}
    for key,value in original.items():
        if isinstance(value, dict):
            copied_dict[key] = deepcopy_v2(value)
        elif isinstance(value, list):
            copied_dict[key] = [deepcopy_v2(item) if isinstance(item, dict) else item for item in value]
        elif isinstance(value, tuple):
            copied_dict[key] = tuple(deepcopy_v2(item) if isinstance(item, dict) else item for item in value)
        elif isinstance(value, set):
            copied_dict[key] = {deepcopy_v2(item) if isinstance(item, dict) else item for item in value}
        else:
            copied_dict[key] = value
    return copied_dict

##############################################
###### Make directonary for input files ######
##############################################
for i in tqdm(range(1000)):
	memo = {}
	dir_list1.append(my_location+"/"+str(i))
	dir_list[str(i)] = my_location+"/"+str(i)
	yaml_mechanism[str(i)] = copy.deepcopy(yaml_files)
	#yaml_mechanism[str(i)] = yaml_files
	#yaml_mechanism[str(i)] = deepcopy_v2(yaml_files)
	#yaml_mechanism[str(i)] = marshal.loads(marshal.dumps(yaml_files))	
	cantera_file[str(i)],exact_file[str(i)],flame_master_file[str(i)],run_file[str(i)] = make_input_files()
#for i in yaml_mechanism:
#	print(id(i),sys. getsizeof(i))
#raise AssertionError("Stop")
###############################################
###### Make dir in parallel            ########
###############################################
#for i in tqdm(range(36000)):
#	os.mkdir(str(i))
#	os.chdir(str(i))
#	os.mkdir("output")
#	os.chdir("..")
	
start_dir = time.time()
workers = 100
W = Worker(workers)
W.do_job_map(dir_list1)
stop_dir = time.time()
print(f"\nTime taken for generating dir is {stop_dir - start_dir}")
del W
###############################################
###### Make input files in parallel    ########
###############################################

yaml_list = []
mech = []
thermo = []
trans = []
instring = []
run = []
extract = []

for i in tqdm(range(1000)):
	instring.append(cantera_file[str(i)])
	yaml_list.append(yaml_mechanism[str(i)])
	run.append(run_file[str(i)])
	extract.append(exact_file[str(i)])


params = list(zip(instring,run,dir_list1,extract,yaml_list))
#print(f"size of zip file is {sys.getsizeof(params)}")
#raise AssertionError("Stop!")
start_input = time.time()
V = Worker(100)
V.do_job_map_create(params)
stop_input = time.time()
print(f"\nTime taken for generating dir is {stop_input - start_input}")
del V	
start_exu = time.time()
X = Worker(150)
file_n = []
length = []
for i in tqdm(range(1000)):
	file_n.append("run")
	#length.append(len(dir_list1))
args = list(zip(dir_list1,file_n))
X.do_job_async_(args)
stop_exu = time.time()
print(f"\nTime taken for generating dir is {stop_exu - start_exu}")
del X
