import multiprocessing
import subprocess
import time
import os
import sys
#This function requires a list containing possible locations where FlameMaster can be executed. Function moves the current directory to those from the list and execute a run commnad which is 
#a shell script to run ScanMan and FlameMaster.

def run_script(location, n):
	os.chdir(location[:-3])
	subprocess.call(location)
	return (location, n)
	
	
#Create a bash script containing the execution calls for FlameMaster and ScanMan. run file is made executable to be called to execute.
def make_run_file(location,  thermo_loc, trans_loc,fsc):
	'''
	s = """#!/bin/bash
ScanMan -i {location}mechanism.mech -t {thermo_l} -m {trans_l} {fsc} -o {location}mechanism.pre -3sr &> {location}Scan
FlameMan -i {location}FlameMaster.input -r {location}mechanism.pre -o {location}output &> {location}Flame
""".format(location = location[:-3], thermo_l = thermo_loc, trans_l = trans_loc)
	#f = open(location,'w')
	#f.write(s)
	f.close()
	'''
	#print(location)
	subprocess.call(["chmod","+rwx",location])


#Callback function of the command apply_async function (part of python multiprocessing module
#prints the progress of simulation and writes the execution location to the progress file	
progress = []
def log_result(result):
	global progress
	progress.append(result[0])
	sys.stdout.write("\r{:06.2f}% of simulation is complete".format(len(progress)/float(result[1])*100))
	sys.stdout.flush()

	#l.acquire()
	f = open('progress','a')
	f.write(result[0]+"\n")
	f.close()
	#l.release()


#performs FlameMaster simulations using multiprocessing module. location is given as input and function iterate through all the locations and execute run script.. finally time for simulations is printed.
def run_FlameMaster_parallel(locations,allowed_count, thermo_loc, trans_loc,s_p_location,fsc):
	home_dir = os.getcwd()
	count = 0
	t_list = []
	start_time = time.time()
	
	#print("Setting up the permissions for simulation......\n\n\n This may take a while... Please be patient...\n\n ")
		
	#for i in locations:
	#	make_run_file(i, thermo_loc, trans_loc,fsc)
	
	print("Executing the code......\n\n\n This may take a while... Please be patient...\n\n ")	
		#Safety measure when user gives a parallel thread count greater than the computer can handle
	if allowed_count > multiprocessing.cpu_count():
		allowed_count = multiprocessing.cpu_count()/2
		
		
	Parallel_jobs = multiprocessing.Pool(allowed_count-1)
	for i in locations:
		Parallel_jobs.apply_async(run_script, args = (i,len(locations)), callback = log_result)
		
	Parallel_jobs.close()
	Parallel_jobs.join()

		
	dt = int(time.time() - start_time)
	hours = dt/3600
	minutes = (dt%3600)/60
	seconds = dt%60
	#os.system("clear")
	print("Performed {} Simulations....".format(len(locations)))
	print("Time for performing simulations : {h} hours,  {m} minutes, {s} seconds\n................................................ ".format(h = hours, m = minutes, s =seconds))
	os.chdir(home_dir)
