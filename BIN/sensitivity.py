import os,json
import data_management,FlameMaster_in_parallel
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import simulation_manager

def analysis(optInputs,iFile,case_dir,rps_order,activeParameters, reaction_index,fallOffCurve_index, thirdBody_index, thermo_index, transport_index, mech_file_location, fileType, rxnUnsrt_data, focUnsrt_data,tbdUnsrt_data, thermoUnsrt_data, transportUnsrt_data, target_list, fuel, global_reaction, thermo_file_location, trans_file_location,startProfile_location,design_type,parallel_threads,file_specific_input,rIndex,unsrt_data,manipulation_dict,sim_dict,response_surface_dict,selectedParams,manipulationDict):
	if "SA" not in os.listdir():
		os.mkdir("SA")
		os.chdir("SA")
		sensDir = os.getcwd()
		os.mkdir("Data")
		os.chdir("Data")
		os.mkdir("NumericalAnalysis")
		os.mkdir("ResponseSurface")
		os.mkdir("Simulations")
		os.chdir("..")
		os.mkdir("Plots")
		os.chdir("Plots")
		os.mkdir("SingularValues")
		os.chdir("..")
	else:
		os.chdir("SA")
		
		sensDir = os.getcwd()

	#################################################################################
	########           SIMULATIONS              #####################################
	########              FOR                   #####################################
	########        SENSITIVITY ANALYSIS        #####################################
	#################################################################################

	sim_type = 'sa'
	if os.path.isfile("progress") == False:
		#print(type(optInputs))
		activeReactions = []
		#FlameMaster_Execution_location,sa_manipulation_dict = simulation_manager.SM(optInputs,iFile,case_dir,rps_order,activeParameters, reaction_index,fallOffCurve_index, thirdBody_index, thermo_index, transport_index, mech_file_location, fileType, rxnUnsrt_data, focUnsrt_data,tbdUnsrt_data, thermoUnsrt_data, transportUnsrt_data, target_list, fuel, global_reaction, thermo_file_location, trans_file_location,startProfile_location,design_type,sim_type,parallel_threads).make_directories_for_simulation()	
		FlameMaster_Execution_location,sa_manipulation_dict = simulation_manager.SM(optInputs,iFile,case_dir,rps_order,activeParameters, reaction_index,fallOffCurve_index, thirdBody_index, thermo_index, transport_index, mech_file_location, fileType,unsrt_data, rxnUnsrt_data, focUnsrt_data,tbdUnsrt_data, thermoUnsrt_data, transportUnsrt_data, target_list, fuel, global_reaction, thermo_file_location, trans_file_location,startProfile_location,design_type,sim_type,parallel_threads,selectedParams,manipulationDict,activeReactions).make_dir_in_parallel()
		#print(FlameMaster_Execution_location)
		locations = open(sensDir+"/locations",'+a')
		mani_list = open(sensDir+"/manipulation_list",'w')
		mani_key  = open(sensDir+"/manipulation_key",'w')
		for i in sa_manipulation_dict:
			mani_key.write(i+'\n')
		mani_list.write(str(sa_manipulation_dict))
		for i in FlameMaster_Execution_location:
			locations.write(i+'\n')
		locations.close()
		mani_key.close()
		mani_list.close()
		#FlameMaster_in_parallel.run_FlameMaster_parallel(FlameMaster_Execution_location, parallel_threads, thermo_file_location, trans_file_location,startProfile_location,file_specific_input)
	else:
		print("Progress file detected\n")
		sa_manipulation_dict = {}
		progress = open(sensDir+"/progress",'r').readlines()
		#manipulation_list = open(sensDir+"/manipulation_list",'r').read().replace("\'","\"")
		#sa_manipulation_dict  = json.loads(manipulation_list)
		print("Reading the progress file \n")
		FlameMaster_Execution_location = []
		with open(sensDir+"/locations") as infile:
    			for line in infile:
        			FlameMaster_Execution_location.append(line)
		#FlameMaster_Execution_location = open(sensDir+"/locations",'r').readlines()
		missing_location = data_management.find_missing_location(FlameMaster_Execution_location, progress)
		if len(missing_location) != 0:
			FlameMaster_in_parallel.run_FlameMaster_parallel(missing_location, parallel_threads, thermo_file_location, trans_file_location,startProfile_location,file_specific_input)


	manipulation_dict[sim_type] = sa_manipulation_dict
	sim_dict = {}
	temp_sim = {}	
	for case in case_dir:
		os.chdir("case-"+str(case))
		data_sheet,failed_sim = data_management.generate_target_value_tables(FlameMaster_Execution_location, target_list, case, fuel,sim_type)
		temp_sim[str(case)] = data_sheet
		f = open('../Data/Simulations/sim_data_case-'+str(case)+'.lst','w').write(data_sheet)
		g = open('../Data/Simulations/failed_sim_data_case-'+str(case)+'.lst','w').write(failed_sim)
		#f.close()
		os.chdir(sensDir)
	sim_dict[sim_type] = temp_sim

	temp_rs = {}	
	response_surface_dict = {}
	active_PRS_dict = {}
	active_reaction_per_case = {}
	M = 3.0/np.log(10)
	for k,case in enumerate(target_list):
		if design_type == "A-facto":	
			T = case.temperature
			theta = np.array([1,np.log(T),-1/T])
			bounds,init_guess = case.create_response_surface(activeParameters,unsrt_data,1)
			temp_rs[str(k)] = case.resCoef
			s_a = []
			I_a = []
			count = 1
			coefficientDict = {} 
			for i in range(len(case.resCoef[1:])):
				coefficientDict[selectedParams[i]] = case.resCoef[i]
				
			count = 1
			sensitivity_dictionary = {}
			impactFactorDict = {}
			"""
			Map the active parameters of each reaction to the reaction
			"""
			unsrt_params = []
			rxn_params = {}
			factor_per_params = {}
			perturbation_factor = {}
			#print(reaction_index)
			for rxn in unsrt_data:
				#print(rxn)
				temp = list(unsrt_data[rxn].activeParameters)
				perturbation_factor[rxn] = unsrt_data[rxn].perturb_factor[0]
				mat = unsrt_data[rxn].perturb_factor
				#print(mat)
				for i,j in enumerate(temp):
					rxn_params[j] = rxn
					factor_per_params[j] = mat
					unsrt_params.append(mat[i])
					
			#print(unsrt_params)
			"""
			Map the sensitive parameters to the reaction
			"""
			for i in unsrt_data:
				#print(i)
				unsrt_fact = perturbation_factor[i]
				#print(unsrt_fact)
				param_length = len(unsrt_data[i].activeParameters)
				
				s_a.append(case.resCoef[count])
				I_a.append(case.resCoef[count]*unsrt_fact)
				sensitivity_dictionary[i] = np.array([case.resCoef[count]])
				impactFactorDict[i] = np.array([case.resCoef[count]])
				count+=1
					
				
			s_a = np.asarray(s_a)/case.resCoef[0]
						
			"""
			Finding the Impact factor
			"""
			I_a = np.asarray(I_a)/case.resCoef[0]
			fig = plt.figure()
			
			"""
			Sorting the impact factor
			-------------------------
			making dict and list for both the sensitivity analysis and Impact factor
			"""
			
			data = {}			
			Idata = {}
			ticks = []
			Iticks = []
			for ind,i in enumerate(s_a):
				data[str(i)] = rIndex[ind]
				ticks.append(ind)
			for ind,i in enumerate(I_a):
				Idata[str(i)] = rIndex[ind]
				Iticks.append(ind)
			
			"""
			Sorting the impact factor
			-------------------------
			Sorting process
			"""
			
			sort_rlist = []			
			sort_Irlist = []
			for ind,i in enumerate(sorted(s_a,key=abs)):
				sort_rlist.append(data[str(i)])
			sort_alist = sorted(s_a,key=abs)
				
			for ind,i in enumerate(sorted(I_a,key=abs)):
				sort_Irlist.append(Idata[str(i)])	
			sort_Ialist = sorted(I_a,key=abs)
			
			
			"""
			Plotting the sensitivites and Impact factors
			---------------------------
			"""
			
			y_pos = range(0,len(s_a))
			#print(y_pos)
			plt.barh(y_pos,sorted(s_a,key =abs), alpha=0.51)
			plt.yticks(y_pos, sort_rlist)
			plt.xlabel(r'normalized sensitivities $S_i = \frac{\partial ln(S_{u}^{o})}{\partial x_i}}$')
			#plt.title('Sensitivity Analysis using Response surface method')
			plt.savefig('./Plots/sensitivity_A_'+str(target_list.index(case))+'.png',bbox_inches="tight")
			plt.close()
						
			y_pos = range(0,len(I_a))
			#print(y_pos)
			plt.barh(y_pos,sorted(I_a,key =abs), alpha=0.51)
			plt.yticks(y_pos, sort_Irlist)
			plt.xlabel(r'normalized optimization potential $I_i = f(T)* \frac{\partial ln(S_{u}^{o})}{\partial x_i}$')
			#plt.title('Sensitivity Analysis using Response surface method')
			plt.savefig('./Plots/impact_A_'+str(target_list.index(case))+'.png',bbox_inches="tight")
			plt.close()
									
			
			sorted_Params = []
			
			Params = (np.asarray(case.resCoef[1:])*unsrt_params)/float(case.resCoef[0])
			#Params = (np.asarray(case.resCoef[1:]))/float(case.resCoef[0])
			
			Max_param = max(abs(np.asarray(Params)))
			
			
			"""
			Selection criteria for the reactions in sparcity effect
			"""
			
			criteria = 0.05*Max_param
			#criteria = 0.75*Max_param
			activeindex = {}
			count = 0
			
			"""
			Need to resuffle the rxn whose params are active
			rather than the individual parameters, we take the 
			rxn as active or not
			
			If rxn is active, all the rxn parameters are flagged as
			active
			==================================================
			Map the sensitive params to the active parameters
			"""
			for i,ele in enumerate(Params):
				if abs(ele)>criteria:
					activeindex[selectedParams[i]] = int(1)
					count+=1
				else:
					activeindex[selectedParams[i]] = int(0)
					
			#print(f"\tActiveParams before suffling \n \t\t{activeindex}\n")
			"""
			Map the sensitive parameters to the reactions
			"""
			active_rxn_index = {}
			for i in activeindex:
				if activeindex[i] == 1:
					active_rxn_index[rxn_params[i]] = 1
			
			active_rxn_dict = {}
			#active_rxn_list = []
			for i in unsrt_data:
				if i in active_rxn_index:
					active_rxn_dict[i] = 1
					#active_rxn_list.append(i)
					continue
				else:
					active_rxn_dict[i] = 0
					active_rxn_index[i] = 0
			#print(active_rxn_dict)
			
			#print(active_rxn_index)	
			final_active_index = {}
			count = 0
			for ele in active_rxn_dict:
				if active_rxn_dict[ele] == 1:
					for params in unsrt_data[ele].activeParameters:
						final_active_index[params] = 1
						count+=1
				else:
					for params in unsrt_data[ele].activeParameters:
						final_active_index[params] = 0
			
			
			string_sensitivity = ""
			for i in active_rxn_dict:
				string_sensitivity+=f"{i}\t{active_rxn_dict[i]}\n"
			
			file_sensitivity = open("sensitivity_"+str(target_list.index(case))+".csv","w").write(string_sensitivity)
			#print(f"\tActiveParams after suffling \n \t\t{final_active_index}\n")	
			#print(active_rxn_index)
			#print(count)	
			#print(k,count)
			active_PRS_dict[k] = final_active_index
			active_reaction_per_case[k] = active_rxn_dict
			
		else:
			
			T = case.temperature
			theta = np.array([1,np.log(T),-1/T])
			bounds,init_guess = case.create_response_surface(activeParameters,unsrt_data,1)
			temp_rs[str(k)] = case.resCoef
			s_a = []
			s_n = []
			s_ea = []
			I_a = []
			I_n = []
			I_ea = []
			count = 1
			#print(len(case.resCoef[1:]))
			coefficientDict = {} 
			for i in range(len(case.resCoef[1:])):
				coefficientDict[selectedParams[i]] = case.resCoef[i]
				
			count = 1
			sensitivity_dictionary = {}
			impactFactorDict = {}
			"""
			Map the active parameters of each reaction to the reaction
			"""
			unsrt_params = []
			rxn_params = {}
			cholesky_per_params = {}
			choleskyMatrix = {}
			for rxn in unsrt_data:
				temp = list(rxnUnsrt_data[rxn].activeParameters)
				choleskyMatrix[rxn] = rxnUnsrt_data[rxn].cholskyDeCorrelateMat
				mat = rxnUnsrt_data[rxn].cholskyDeCorrelateMat
				for i in temp:
					rxn_params[i] = rxn
					cholesky_per_params[i] = mat
					unsrt_params.append(M*np.linalg.norm(mat.T.dot(theta)))
					
			"""
			Map the sensitive parameters to the reaction
			"""
			for i in unsrt_data:
				#print(i)
				unsrt_fact = M*np.linalg.norm(choleskyMatrix[i].T.dot(theta))
				param_length = len(rxnUnsrt_data[i].activeParameters)
				#print(rxnUnsrt_data[i].activeParameters)
				#print(param_length)
				if param_length == 1:
					s_a.append(case.resCoef[count])
					I_a.append(case.resCoef[count]*unsrt_fact)
					s_n.append(0)
					s_ea.append(0)
					I_n.append(0)
					I_ea.append(0)
					sensitivity_dictionary[i] = np.array([case.resCoef[count],0.0,0.0])
					impactFactorDict[i] = np.array([case.resCoef[count]*unsrt_fact,0.0,0.0])
					count+=1
					
				else:
					s_a.append(case.resCoef[count])
					I_a.append(case.resCoef[count]*unsrt_fact)
					s_n.append(case.resCoef[count+1])
					s_ea.append(case.resCoef[count+2])
					I_n.append(case.resCoef[count+1]*unsrt_fact)
					I_ea.append(case.resCoef[count+2]*unsrt_fact)
					sensitivity_dictionary[i] = np.array([case.resCoef[count],case.resCoef[count+1],case.resCoef[count+2]])
					impactFactorDict[i] = np.array([case.resCoef[count]*unsrt_fact,case.resCoef[count+1]*unsrt_fact,case.resCoef[count+2]*unsrt_fact])
					count+=3
					
				
			s_a = np.asarray(s_a)/case.resCoef[0]
			s_n = np.asarray(s_n)/case.resCoef[0]
			s_ea = np.asarray(s_ea)/case.resCoef[0]
			
			"""
			Finding the Impact factor
			"""
			
			I_a = np.asarray(I_a)/case.resCoef[0]
			I_n = np.asarray(I_n)/case.resCoef[0]
			I_ea = np.asarray(I_ea)/case.resCoef[0]
			fig = plt.figure()
			
			"""
			Sorting the impact factor
			-------------------------
			making dict and list for both the sensitivity analysis and Impact factor
			"""
			
			data = {}
			data_n = {}
			data_ea = {}
			
			Idata = {}
			Idata_n = {}
			Idata_ea = {}
			ticks = []
			Iticks = []
			for ind,i in enumerate(s_a):
				data[str(i)] = rIndex[ind]
				data_n[rIndex[ind]] = s_n[ind]
				data_ea[rIndex[ind]] = s_ea[ind]
				ticks.append(ind)
			for ind,i in enumerate(I_a):
				Idata[str(i)] = rIndex[ind]
				Idata_n[rIndex[ind]] = I_n[ind]
				Idata_ea[rIndex[ind]] = I_ea[ind]
				Iticks.append(ind)
			
			"""
			Sorting the impact factor
			-------------------------
			Sorting process
			"""
			
			sort_rlist = []
			sort_nlist = []
			sort_ealist = []
			
			sort_Irlist = []
			sort_Inlist = []
			sort_Iealist = []
			
			for ind,i in enumerate(sorted(s_a,key=abs)):
				sort_rlist.append(data[str(i)])
				sort_nlist.append(data_n[data[str(i)]])
				sort_ealist.append(data_ea[data[str(i)]])
			sort_alist = sorted(s_a,key=abs)
				
			for ind,i in enumerate(sorted(I_a,key=abs)):
				sort_Irlist.append(Idata[str(i)])
				sort_Inlist.append(Idata_n[Idata[str(i)]])
				sort_Iealist.append(Idata_ea[Idata[str(i)]])	
			sort_Ialist = sorted(I_a,key=abs)
			
			
			"""
			Plotting the sensitivites and Impact factors
			---------------------------
			"""
			
			y_pos = range(0,len(s_a))
			#print(y_pos)
			plt.barh(y_pos,sorted(s_a,key =abs), alpha=0.51)
			plt.yticks(y_pos, sort_rlist)
			plt.xlabel(r'normalized sensitivities $S_i = \frac{\partial ln(S_{u}^{o})}{\partial x_i}}$')
			#plt.title('Sensitivity Analysis using Response surface method')
			plt.savefig('./Plots/sensitivity_A_'+str(target_list.index(case))+'.png',bbox_inches="tight")
			plt.close()
			
			fig = plt.figure()
			y_pos = range(0,len(s_n))
			plt.barh(y_pos, sort_nlist, alpha=0.5)
			plt.yticks(y_pos, sort_rlist)
			plt.xlabel('Sensitivity of temperature exponent (n)')
			#plt.title('Sensitivity Analysis using Response surface method')
			plt.savefig('./Plots/sensitivity_n_'+str(target_list.index(case))+'.png')
			plt.close()
			
			fig = plt.figure()
			y_pos = range(0,len(s_ea))
			plt.barh(y_pos, sort_ealist,alpha=0.5)
			plt.yticks(y_pos, sort_rlist)
			plt.xlabel('Sensitivity of Activation Energy (Ea)')
			#plt.title('Sensitivity Analysis using Response surface method')
			plt.savefig('./Plots/sensitivity_Ea_'+str(target_list.index(case))+'.png')
			plt.close()
			
			
			y_pos = range(0,len(I_a))
			#print(y_pos)
			plt.barh(y_pos,sorted(I_a,key =abs), alpha=0.51)
			plt.yticks(y_pos, sort_Irlist)
			plt.xlabel(r'normalized optimization potential $I_i = f(T)* \frac{\partial ln(S_{u}^{o})}{\partial x_i}$')
			#plt.title('Sensitivity Analysis using Response surface method')
			plt.savefig('./Plots/impact_A_'+str(target_list.index(case))+'.png',bbox_inches="tight")
			plt.close()
			
			fig = plt.figure()
			y_pos = range(0,len(I_n))
			plt.barh(y_pos, sort_Inlist, alpha=0.5)
			plt.yticks(y_pos, sort_Irlist)
			plt.xlabel('Impact factor for temperature exponent (n)')
			#plt.title('Sensitivity Analysis using Response surface method')
			plt.savefig('./Plots/impact_n_'+str(target_list.index(case))+'.png')
			plt.close()
			
			fig = plt.figure()
			y_pos = range(0,len(I_ea))
			plt.barh(y_pos, sort_Iealist,alpha=0.5)
			plt.yticks(y_pos, sort_Irlist)
			plt.xlabel('Impact factor for Activation Energy (Ea)')
			#plt.title('Sensitivity Analysis using Response surface method')
			plt.savefig('./Plots/impact_Ea_'+str(target_list.index(case))+'.png')
			plt.close()
			
			"""
			Plotting the sensitivity in same fig
			"""
			
			#print(sort_rlist)
			fake_data = pd.DataFrame({"index": list(sort_rlist), 0: sort_alist , 1: sort_nlist, 2: np.asarray(sort_ealist)*10})
			fake_data.set_index("index",drop=False)
			fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(8, 8), frameon=False)
			
			fake_data[0].plot.barh(ax=ax1)
			fake_data[1].plot.barh(ax=ax2)
			fake_data[2].plot.barh(ax=ax3)
			
			ax1.set_yticks(ticks,sort_rlist)
			ax1.set_xlabel(r'$\partial ln(\eta) / \partial \zeta_{\alpha}$')
			ax2.set_xlabel(r'$\partial ln(\eta) / \partial \zeta_{n}$')
			ax3.set_xlabel(r'$\partial ln(\eta) / \partial \zeta_{\epsilon} (\times 10^{-1})$')
			fig.savefig('./Plots/sensitivity'+str(target_list.index(case))+'.png',bbox_inches="tight")
			plt.close()
			"""
			Plotting the Impact factor in same plot
			"""
			
			fake_data = pd.DataFrame({"index": list(sort_Irlist), 0: sort_Ialist , 1: sort_Inlist, 2: np.asarray(sort_Iealist)*10})
			fake_data.set_index("index",drop=False)
			fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(8, 8), frameon=False)
			
			fake_data[0].plot.barh(ax=ax1)
			fake_data[1].plot.barh(ax=ax2)
			fake_data[2].plot.barh(ax=ax3)
			
			ax1.set_yticks(ticks,sort_rlist)
			ax1.set_xlabel(r'$f(T)*\frac{\partial ln(\eta) }{\partial \alpha}$')
			ax2.set_xlabel(r'$f(T)*\frac{\partial ln(\eta)}{\partial n}$')
			ax3.set_xlabel(r'$f(T)*\frac{\partial ln(\eta)}{\partial \epsilon} (\times 10^{-1})$')
			fig.savefig('./Plots/Impact'+str(target_list.index(case))+'.png',bbox_inches="tight")
			plt.close()
			
			
			
			sorted_Params = []
			
			Params = (np.asarray(case.resCoef[1:])*unsrt_params)/float(case.resCoef[0])
			#Params = (np.asarray(case.resCoef[1:]))/float(case.resCoef[0])
			
			Max_param = max(abs(np.asarray(Params)))
			
			
			"""
			Selection criteria for the reactions in sparcity effect
			"""
			
			criteria = 0.05*Max_param
			#criteria = 0.75*Max_param
			activeindex = {}
			count = 0
			
			"""
			Need to resuffle the rxn whose params are active
			rather than the individual parameters, we take the 
			rxn as active or not
			
			If rxn is active, all the rxn parameters are flagged as
			active
			==================================================
			Map the sensitive params to the active parameters
			"""
			for i,ele in enumerate(Params):
				if abs(ele)>criteria:
					activeindex[selectedParams[i]] = int(1)
					count+=1
				else:
					activeindex[selectedParams[i]] = int(0)
					
			#print(f"\tActiveParams before suffling \n \t\t{activeindex}\n")
			"""
			Map the sensitive parameters to the reactions
			"""
			active_rxn_index = {}
			for i in activeindex:
				if activeindex[i] == 1:
					active_rxn_index[rxn_params[i]] = 1
			
			active_rxn_dict = {}
			#active_rxn_list = []
			for i in unsrt_data:
				if i in active_rxn_index:
					active_rxn_dict[i] = 1
					#active_rxn_list.append(i)
					continue
				else:
					active_rxn_dict[i] = 0
					active_rxn_index[i] = 0
			#print(active_rxn_dict)
			
			#print(active_rxn_index)	
			final_active_index = {}
			count = 0
			for ele in active_rxn_dict:
				if active_rxn_dict[ele] == 1:
					for params in rxnUnsrt_data[ele].activeParameters:
						final_active_index[params] = 1
						count+=1
				else:
					for params in rxnUnsrt_data[ele].activeParameters:
						final_active_index[params] = 0
			
			
			string_sensitivity = ""
			for i in active_rxn_dict:
				string_sensitivity+=f"{i}\t{active_rxn_dict[i]}\n"
			
			file_sensitivity = open("sensitivity_"+str(target_list.index(case))+".csv","w").write(string_sensitivity)
			#print(f"\tActiveParams after suffling \n \t\t{final_active_index}\n")	
			#print(active_rxn_index)
			#print(count)	
			#print(k,count)
			active_PRS_dict[k] = final_active_index
			active_reaction_per_case[k] = active_rxn_dict
	response_surface_dict[sim_type]= temp_rs
	os.chdir("..")
	return sensDir,manipulation_dict,sim_dict,response_surface_dict,coefficientDict,active_PRS_dict,bounds,init_guess,active_reaction_per_case
