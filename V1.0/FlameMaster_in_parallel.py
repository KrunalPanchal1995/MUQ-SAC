import multiprocessing
from multiprocessing import Pool
import concurrent.futures
import subprocess
import time
import os
import sys

#This function requires a list containing possible locations where FlameMaster can be executed. Function moves the current directory to those from the list and execute a run commnad which is 
#a shell script to run ScanMan and FlameMaster.

def run_script(location, n):
	os.chdir(location[:-3])
	subprocess.call(["./run_convertor"])
	subprocess.call(["./run"])
	return (location, n)

def run_script_map(location):
	os.chdir(location[:-3])
	subprocess.call(location)
	#print(os.getcwd())
	#subprocess.call(location)

#	print("\rExecuted FlameMaster at {}".format(location[location.index("case"):-4]))

#def run_script(location):
#	os.chdir(location[:-3])
#	subprocess.call(location)
	#return (location, n)	
	
#Create a bash script containing the execution calls for FlameMaster and ScanMan. run file is made executable to be called to execute.
#Callback function of the command apply_async function (part of python multiprocessing module
#prints the progress of simulation and writes the execution location to the progress file	
#progress = []
#def log_result(result):
#	global progress
#	progress.append(result[0])
#	sys.stdout.write("\r{:06.2f}% of simulation is complete".format(len(progress)/float(result[1])*100))
#	sys.stdout.flush()
#	Parallel_jobs.terminate()
	#l.acquire()
	#f = open('progress','a')
	#f.write(result[0]+"\n")
	#f.close()
	#l.release()

class Worker():
	def __init__(self, workers,locations):
		self.pool = multiprocessing.Pool(processes=workers)
		self.locations = locations
		self.progress = []
		
	def callback(self, result):
		self.progress.append(result[0])
		sys.stdout.write("\r{:06.2f}% of simulation is complete".format(len(self.progress)/float(result[1])*100))
		sys.stdout.flush()
		#self.pool.terminate()
		#l.acquire()
		f = open('progress','a')
		f.write(result[0]+"\n")
		f.close()
		#l.release()

	def do_job_async(self):
		for args in self.locations:
			self.pool.apply_async(run_script, 
                                  args=(args,len(self.locations)), 
                                  callback=self.callback)
		self.pool.close()
		self.pool.join()
	
	def do_job_map(self):
		self.pool.map_async(run_script_map, self.locations)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()
	
	def do_job_executor(self):
		with concurrent.futures.ProcessPoolExecutor() as executor:
			executor.map(run_script_map,self.locations)
#performs FlameMaster simulations using multiprocessing module. location is given as input and function iterate through all the locations and execute run script.. finally time for simulations is printed.
def run_FlameMaster_parallel(locations,allowed_count, thermo_loc, trans_loc,s_p_location,fsc):
	home_dir = os.getcwd()
	count = 0
	t_list = []
	start_time = time.time()
	allowed_count = int(allowed_count)
	#print("Setting up the permissions for simulation......\n\n\n This may take a while... Please be patient...\n\n ")
		
	#for i in locations:
	#	make_run_file(i, thermo_loc, trans_loc,fsc)
	
	print("Executing the code......\n\n\n This may take a while... Please be patient...\n\n ")	
		#Safety measure when user gives a parallel thread count greater than the computer can handle
	if allowed_count > multiprocessing.cpu_count():
		allowed_count = int(multiprocessing.cpu_count()/2)
		
	"""
	#	Threading using concurrent.futures
	with concurrent.futures.ThreadPoolExecutor() as executor:
		results = executor.map(run_script,locations)
		for result in results:
			print(result)	
	
	"""
	"""
	#	Threading using concurrent.futures
	#	(Threading is useful for io bound processes)
	
	with concurrent.futures.ThreadPoolExecutor() as executor:
		results = [executor.submit(run_script,loc) for loc in locations]
		
		for f in concurrent.futures.as_completed(results):
			print(f.result())	
	
	"""
	"""
	#	Multiprocessing using concurrent.futures
	
	with concurrent.futures.ProcessPoolExecutor() as executor:
		executor.map(run_script,locations)
		#executor.shutdown(wait=True)
		#for result in results:
		#	print(result)	
	
	"""
	"""
	#	Multiprocessing using concurrent.futures
	#	(Threading is useful for io bound processes)
	
	with concurrent.futures.ProcessPoolExecutor() as executor:
		results = [executor.submit(run_script,loc) for loc in locations]
		
		for f in concurrent.futures.as_completed(results):
			print(f.result())	
	
	"""
	
	w = Worker(allowed_count,locations)
	w.do_job_async()	
	#w.do_job_map()
	#w.do_job_executor()
	"""
	#	Multiprocessing using Processes
	#	Use pickel
	processes = []
	for loc in locations:
		p = multiprocessing.Process(target=run_script,args=[loc])
		p.start()
		processes.append(p)
		
	for process in processes:
		process.join()
	
	"""			
	dt = int(time.time() - start_time)
	hours = dt/3600
	minutes = (dt%3600)/60
	seconds = dt%60
	#os.system("clear")
	print("Performed {} Simulations....".format(len(locations)))
	print("Time for performing simulations : {h} hours,  {m} minutes, {s} seconds\n................................................ ".format(h = hours, m = minutes, s =seconds))
	os.chdir(home_dir)
