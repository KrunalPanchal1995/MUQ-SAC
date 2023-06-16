import re, os, math
import time
import numpy as np
import shutil
import subprocess
import subprocess as sub
import threading
import combustion_target_class


class RunCmd(threading.Thread):
    def __init__(self, cmd, timeout):
        threading.Thread.__init__(self)
        self.cmd = cmd
        self.timeout = timeout

    def run(self):
        self.p = sub.Popen(self.cmd)
        self.p.wait()

    def Run(self):
        self.start()
        self.join(self.timeout)

        if self.is_alive():
            self.p.terminate()      #use self.p.kill() if process needs a kill -9
            self.join()

#Extract_reactinons pre exponential factor from mechanism file using sorted numbers. dummy function not called in main code. replaced with a regular expression.
def extract_reaction_coeff(mech):
	data = mech.readlines()
	reaction_set = []
	reaction = ''
	pre_exp_factor = {}
	for n,line in enumerate(data):
		if('->' in line):
			reaction = line
			i = 0
			while(True):
				if('}' not in reaction):
					reaction = reaction + data[n+i+1]
				else:
					break
				i +=1
			reaction_set.append(reaction)

	for r in reaction_set:
		param = r.split()
		index = param[0][:-1]
		for a in param:
			if a=="a" or a== "A":
				pre_exp_factor[index] = float(param[param.index(a)+2]) 
	
	return pre_exp_factor	
		
#extract uncertainity data for each reaction from user created uncertainty file. It should be a two column file with first column containing index and second one with uncertainties
def extract_index_and_uncertainty(uns):
	data = uns.readlines()
	uncertainty_a = {}
	uncertainty_n = {}
	uncertainty_e = {}
	index = []
	for i in data:
		if "#" in i:
			i = i[:i.index('#')]
		seg = i.split()
		if len(seg) == 4:
			uncertainty_a[seg[0]] = float(seg[1]);
			uncertainty_n[seg[0]] = float(seg[2]);
			uncertainty_e[seg[0]] = float(seg[3]);
			index.append(seg[0])
		else:
			continue
	return uncertainty_a, uncertainty_n, uncertainty_e, index

#Once the simulations are completed, this function goes through all the locations where FlameMaster simulation is performed and read the output value based on the type of target.
#and print as a text file containg four columns. First three contain location information and last column contain  target value.
def generate_target_value_tables(locations, t_list, case, fuel):
	#print(locations)
	list_fuel = []
	if "dict" in str(type(fuel)):
		for i in fuel:
			list_fuel.append(fuel[i])
		fuel = list_fuel[0]
	
	#print(fuel)
	data_loc = []
	for location in locations:
		#print(location)
		list_loc = location.split("/")
		for i in list_loc:
			if i == "case-"+str(case).strip():
				data_loc.append(location.strip("\n"))

	
	data =''
	failed_sim = ""
	ETA = []
	for i in data_loc:
		pathList = i.split("/")
		file_loc = open("./eta_file_location.txt","+a")
		start = pathList.index("case-{}".format(case))
		#print(start)
		#print(i[:-4]+"output/")
		#print(t_list[case])
		#print(data_loc.index(i))
		#print(i)
		
		eta,file_path = extract_output(t_list[case],fuel,i[:-3]+"output/",data_loc.index(i))
		if "N/A" in str(eta):
			#print(eta)
			#print(file_path)
			file_loc.write(file_path+"\n")
			folderName = pathList[start+1]
			#print(folderName)
			#print(eta)
			failed_sim += "{}\t{}\n".format(folderName , eta)
			file_loc.close()
			
		else:	
			#print(eta)
			#print(file_path)
			file_loc.write(file_path+"\n")
			folderName = pathList[start+1]
			#print(folderName)
			#print(eta)
			data += "{}\t{}\n".format(folderName , eta)
			file_loc.close()
			ETA.append(eta)
	return data,failed_sim,ETA

def extract_direct_simulation_values(case,loc,target_list,fuel):
	data =''
	eta_list = []
	failed_list = []
	pathList = loc.strip("\n").split("/")
	#print(pathList)
	#start = pathList.index("case-{}".format(case))
	file_loc = open("./eta_file_location.txt","w")
	eta,file_path = extract_output(target_list[case],fuel,loc,case)
	if "N/A" in str(eta):
		file_loc.write(file_path+"\n")
		folderName = pathList[-3]
		data += "{}\t{}\n".format(folderName , eta)
		file_loc.close()
		failed_list.append(float(eta))	
	else:
		file_loc.write(file_path+"\n")
		folderName = pathList[-3]
		data += "{}\t{}\n".format(folderName , eta)
		file_loc.close()
		eta_list.append(float(eta))
	return eta_list

#function extracts output from the given text file based on the called location where it currently is. Called in the previous function.		
def extract_output(case,fuel,path,index):
	
	eta = string = None
	#print(case.target)

	if case.target.strip() == "Tig":
		if "cantera" in case.add["solver"]:
			if "tau.out" in os.listdir(path):
				out_file = open(path+"tau.out",'r').readlines()
				string = path +"tau.out"
				#print(out_file)
				line = out_file[1].split()
				#print(len(line))
				#print(line)
				if len(line) == 2:
					eta = np.log(float(line[1])*10)	#us/micro seconds
				else:
					eta = np.log(100*10000)
				
			else:
				eta = "N/A"
				string = path
		else:
			#print("Tig is the target")
			out_file = open(path+"tau.out",'r').readlines()
			string = path +"tau.out"
			#out_file = open(path + fuel+"_IgniDelTimes.dout",'r').readlines()
			#string = path + fuel+"_IgniDelTimes.dout"
			
			line = out_file[0].split("	")
			
			#print(line)
			if len(line) == 3:
				#eta = np.log(float(line[2])*1000)
				print(line)
				eta = np.log(float(line[1])*1000000)#us /micro seconds
			else:
				eta = np.log(100*10000)
			#return eta,string
	elif case.target.strip() == "Fls":
		if "cantera" in case.add["solver"]:
			out_file = open(path+"Su.out",'r').readlines()
			string = path +"Su.out"
			line = out_file[1].split()
			#print(line)
			if len(line) == 2:
				eta = np.log(float(line[1]))	#cm/seconds
			else:
				eta = np.log(200)
		else:
			#start = os.getcwd()
			#print(path)
			os.chdir(path)
			start = path
			outfile =""
			flist = os.listdir()
			for i in flist:
				if fuel in i:
					if i.endswith("{}".format(int(case.temperature))):
						outfile =open(i,'r').readlines()
						string = path+i
						for line in outfile:
							if "burningVelocity" in line:
								#print(line.split()[2])
								#eta = np.log(float(line.split()[2]))
								eta = float(line.split()[2])
								#os.chdir("..")
								#return eta,string
					else:
						if i.endswith("noC"):
							#os.remove(i)
							continue			
			if outfile == "":
				#print("{sim_index} in {case_index} is not converged. Changing the start profile and re-calculating the flame speed as outfile is {outfile}".format(sim_index = index, case_index = case.case_index,outfile = outfile))
				#print("\n Current location is {}".format(os.getcwd()))
				start_profiles = open("../../eta_file_location.txt","r").readlines()			
				count = 1
				os.chdir(path)
				os.chdir("..")
				start = os.getcwd() 
				while outfile == "":
					print(start_profiles)
					for start_profile in start_profiles:
						os.chdir(start)
						print("{}: Using profile location of {}".format(count,start_profile))
						print(os.getcwd())
						#print(os.listdir())
						FM_input = open("FlameMaster.input","r").readlines()
						temp = open("New.txt","+a")
						for line in FM_input:
							if "StartProfilesFile" in line.split():
								string = "StartProfilesFile is "+start_profile
								temp.write(line.replace(line,string))
							else:
								temp.write(line)
						temp.close()
						print("\n \t Creating FlameMaster input file with new start profile \n")
						#print(os.listdir())
						os.remove("FlameMaster.input")
						os.rename("New.txt","FlameMaster.input")
						#run FlameMaster locally
						#print(os.listdir())
						kill = int(180)
						print("\n \t Running the FlameMaster locally...\n \t Kill Time out of {} sec".format(kill))
						RunCmd(["./run"],kill).Run()
						 
						count +=1
						os.chdir("output")
						flist = os.listdir()
						for i in flist:
							if fuel in i:
								if i.endswith("{}".format(int(case.temperature))):
									outfile =open(i,'r').readlines()
									string = path+i
									for line in outfile:
										if "burningVelocity" in line:
											#print(line.split()[2])
											#eta = np.log(float(line.split()[2]))
											eta = float(line.split()[2])
											string = path+i
											os.chdir(start)
											#return eta,string	
								else:
									if i.endswith("noC"):
										#os.remove(i)
										print("Case not converged yet!!!")
										continue
									else:
										continue	
				'''
				for i in flist: 
					if i.endswith("noC"):
						outfile =open(i,'r').readlines()
						for line in outfile:
							if "burningVelocity" in line:
								print(line.split()[2])
								eta = line.split()[2]
								os.chdir(start)
								return eta
					'''
		
			#exit() 
		#for i in outfile:
		#	if "burningVelocity" in i:
		#		eta = i.split()[2]
		#return eta,string
	elif case.target.strip() == "Flf":
		if "cantera" in case.add["solver"]:
			out_file = open(path+"result.dout",'r').readlines()
			string = path +"result.dout"
			line = out_file[0].split("	")
			if len(line) == 2:
				eta = float(line[1])*100	#%mole
			else:
				eta = 100
		else:
			#print("Tig is the target")
			out_file = open(path+"result.dout",'r').readlines()
			string = path +"result.dout"
			line = out_file[0].split()
			#print(line)
			if len(line) == 2:
				#eta = np.log(float(line[2])*1000)
				eta = float(line[1])*100 #%mole
			else:
				eta = 100
	
	elif case.target.strip() == "Flw":
		if "cantera" in case.add["solver"]:
			if "slope" in case.add["flw_method"]:
				out_file = open(path+"rate.csv",'r').readlines()
				string = path +"rate.csv"
				line = out_file[0].split()
				if len(line) == 3:
					#eta = np.log(float(line[2])*1000)
					eta = float(line[1]) #ppm/ms
				else:
					eta = 100
			else:	
				#print("Tig is the target")
				out_file = open(path+"rate.csv",'r').readlines()
				string = path +"time.csv"
				line = out_file[0].split()
				#print(line)
				if len(line) == 3:
					#eta = np.log(float(line[2])*1000)
					eta = float(line[1]) #ms
				else:
					eta = 100
		else:
			if "slope" in case.add["flw_method"]:
				out_file = open(path+"rate.csv",'r').readlines()
				string = path +"rate.csv"
				line = out_file[0].split()
				if len(line) == 3:
					#eta = np.log(float(line[2])*1000)
					eta = float(line[1]) #ppm/ms
				else:
					eta = 100
			else:	
				#print("Tig is the target")
				out_file = open(path+"rate.csv",'r').readlines()
				string = path +"time.csv"
				line = out_file[0].split()
				#print(line)
				if len(line) == 3:
					#eta = np.log(float(line[2])*1000)
					eta = float(line[1]) #ms
				else:
					eta = 100
#			start = os.getcwd()
#			os.chdir(path)
#			output_file = [i for i in os.listdir() if i.startswith("X1_")]
#			print(output_file)
#			outfile = open(path+"/"+output_file[0],"r").readlines()
#			string = path+output_file[0]
#			header_line = outfile[1]
#			header = header_line.split()
#			if "X-"+case.species in header:
#				ind = header.index('X-'+case.species)
#			else:
#				print("Invalid species name: {} in input".format(case.species))
#			eta = float(outfile[-1].split()[ind])*100
			#os.chdir(start)
		
	#print(eta)
	#print(string)
	return eta,string			
				

#Generate an optimized mechanism based on the optimized vector. 
def generate_optimized_mechanism(mech_file_location, reaction_index, unsrt_data, opt_x):
	mech_file = open(mech_file_location,'r')
	mech = new_mech = mech_file.read()	
	for i in reaction_index:
		w = re.compile(r'\n{}:.*?->.*?a\s*?=.*?(\S*?E\S\d*?\s).*?\}}'.format(i), re.DOTALL | re.IGNORECASE) #Group(1) will extract pre exp factor
		match = re.search(w, mech)
		if match == None:
			print ("Unable to perturb reaction rate... Exiting")
			exit()
		reaction = match.group()
		pre_exp = match.group(1)
		pre_exp_p = "{:.4E} ".format(float(pre_exp)*math.exp(opt_x[reaction_index.index(i)]*math.log(unsrt_data[i])))
		new_reaction = reaction.replace(pre_exp, pre_exp_p)
		new_mech = new_mech.replace(reaction, new_reaction)
		#Perturb backward reaction
		if i[-1] == 'f':
			kb = i[:-1]+'b'
			w = re.compile(r'\n{}:.*?->.*?a\s*?=.*?(\S*?E\S\d*?\s).*?\}}'.format(kb), re.DOTALL | re.IGNORECASE) #Group(1) will extract pre exp factor
			match = re.search(w, mech)
			if match != None:
				reaction = match.group()
				pre_exp = match.group(1)
				pre_exp_p = "{:.4E} ".format(float(pre_exp)*math.exp(opt_x[reaction_index.index(i)]*math.log(unsrt_data[i])))
				new_reaction = reaction.replace(pre_exp, pre_exp_p)
				new_mech = new_mech.replace(reaction, new_reaction)
	
	mech_name = mech_file_location.split('/')[-1]
	opt_mech_name = mech_name.replace(".mech","_optimized.mech")		
	ff = open(opt_mech_name,'w')
	ff.write(new_mech)
	ff.close()	

def getTestingData(sensDir,case):
	#print(os.getcwd())
	#print(f"{sensDir}")
	y_file = open(sensDir+"/Data/Simulations/sim_data_case-"+str(case)+".lst","r").readlines()
	#print(y_file)
	#x_file = open(sensDir+"/Data/Simulations/Generator_list_case-"+str(case)+".csv","r").readlines()
	x_file = open(sensDir+"/Data/Simulations/Beta_list_case-"+str(case)+".csv","r").readlines()
	#x_file = open(sensDir+"/Data/Simulations/NormalizedBeta_list_case-"+str(case)+".csv","r").readlines()
	x_data = np.asarray([np.asarray(np.float_(i[:-1].strip(",").split(","))) for i in x_file])
	y_data = np.asarray([float(i.strip("\n").split()[1]) for i in y_file])
		
	return  x_data, y_data
#create log files containing errors response surface values and optimized values					
def make_log_files(case, reaction_index, opt_x):
	co_eff = open("response_surface_co_efficients.out",'w')
	co_eff.write("#Response surface coefficients\n")
	for i in case:
		co_eff.write("{}\t".format(i.case_index))
		for j in i.co_efficients:
			co_eff.write("{}\t".format(j))
		co_eff.write("\n")
	co_eff.close()
	
	vect = open("Optimized_vector",'w')
	vect.write("Normalized Pre Exponential Factors\n")
	for i, j in enumerate(reaction_index):
		vect.write("{}\t{}\n".format(j, opt_x[i]))
	vect.close()
	
	co_eff = open("results_and_errors.out",'w')
	co_eff.write("Unoptimized\tOptimized\t\Experiment\tOld_error\tNew_error\n")
	for i in case:
		opt_eta = i.calculated_target_value(opt_x)
		old_error = abs(i.observed - i.calculated)
		new_error = abs(i.observed - opt_eta)
		co_eff.write("{}\t{}\t{}\t{}\t{}\n".format(i.calculated, opt_eta, i.observed, old_error, new_error))
		
	co_eff.close()
	
#Checks the contents of two files created using simulations. One file containing all the locations where simulation
#should be performed and another containing the locations where it is performed. The differene is returned and used for 
#performing remaining simulations. Helps to resume simulations in case of power failure. 		
def find_missing_location(initial, progress):
	#print(len(initial))
	#print(len(progress))
	missing_locations = []
	if len(initial) == len(progress):
		print("All simulations are completed\n")
	else:	
		for i in initial:
			if i not in progress:
				missing_locations.append(i[:-1])
		print("Found all missing locations\n")
	#print(len(missing_locations))
	return missing_locations
	
