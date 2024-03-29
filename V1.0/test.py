import os,json
import data_management,FlameMaster_in_parallel
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import simulation_manager

def simulations(optInputs,iFile,case_dir,rps_order,activeParameters, reaction_index,fallOffCurve_index, thirdBody_index, thermo_index, transport_index, mech_file_location, fileType, rxnUnsrt_data, focUnsrt_data,tbdUnsrt_data, thermoUnsrt_data, transportUnsrt_data, target_list, fuel, global_reaction, thermo_file_location, trans_file_location,startProfile_location,design_type,parallel_threads,file_specific_input,rIndex,unsrt_data,manipulation_dict,sim_dict,response_surface_dict,selectedParams,manipulationDict,activeIndexDict,activeReactions):

	if "test" not in os.listdir():
		os.mkdir("test")
		os.chdir("test")
		testDir = os.getcwd()
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
		os.chdir("test")
		testDir = os.getcwd()


	print("Generate Directories and mechanism files for testing set\n")
	sim_type = "test"
	if os.path.isfile("progress") == False:
		#FlameMaster_Execution_location , opt_manipulation_dict  = simulation_manager.SM(optInputs,iFile,case_dir,rps_order,activeParameters, reaction_index,fallOffCurve_index,thirdBody_index,thermo_index,transport_index, mech_file_location,fileType, rxnUnsrt_data, focUnsrt_data, tbdUnsrt_data, thermoUnsrt_data, transportUnsrt_data, target_list, fuel, global_reaction, thermo_file_location, trans_file_location,startProfile_location,design_type,sim_type,parallel_threads).make_directories_for_simulation()
		FlameMaster_Execution_location , test_manipulation_dict  = simulation_manager.SM(optInputs,iFile,case_dir,rps_order,activeParameters, reaction_index,fallOffCurve_index,thirdBody_index,thermo_index,transport_index, mech_file_location,fileType,unsrt_data,rxnUnsrt_data, focUnsrt_data, tbdUnsrt_data, thermoUnsrt_data, transportUnsrt_data, target_list, fuel, global_reaction, thermo_file_location, trans_file_location,startProfile_location,design_type,sim_type,parallel_threads,selectedParams,manipulationDict,activeReactions,extra_arg = activeIndexDict).make_dir_in_parallel()
		locations = open(testDir+"/locations",'w')
		mani_list = open(testDir+"/manipulation_list",'w')
		mani_key  = open(testDir+"/manipulation_key",'w')
		for i in test_manipulation_dict:
			mani_key.write(i+'\n')
		mani_list.write(str(test_manipulation_dict))
		mani_key.close()
		mani_list.close()
		for i in FlameMaster_Execution_location:
			locations.write(i+"\n")
		locations.close()
		#FlameMaster_in_parallel.run_FlameMaster_parallel(FlameMaster_Execution_location, parallel_threads, thermo_file_location, trans_file_location,startProfile_location,file_specific_input)

	else:
		print("Progress file detected")
		test_manipulation_dict = {}
		#manipulation_list = open(testDir+"/manipulation_list",'r').read().replace("\'","\"")
		#test_manipulation_dict = json.loads(manipulation_list)
		#print(manipulation_list)
		#manipulation_list = filter_list(manipulation_list)
		#manipulation_key = open(testDir+"/manipulation_key",'r').readlines()
		#for i,key in enumerate(manipulation_key):
		#	test_manipulation_dict[str(key).strip("\n")] = dict(manipulation_list[i].replace('"',"").strip('\n').strip('{').strip('}'))
		#progress = open(testDir+"/progress",'r').readlines()
		FlameMaster_Execution_location = []
		with open(testDir+"/locations") as infile:
    			for line in infile:
        			FlameMaster_Execution_location.append(line)
		
		#FlameMaster_Execution_location = open(testDir+"/locations",'r').readlines()
		#missing_location = data_management.find_missing_location(FlameMaster_Execution_location, progress)
		#if len(missing_location) !=0:
		#	FlameMaster_in_parallel.run_FlameMaster_parallel(missing_location, parallel_threads, thermo_file_location, trans_file_location,startProfile_location,file_specific_input)

	manipulation_dict[sim_type] = test_manipulation_dict
	temp_sim_test = {}
	for case in case_dir:	
		os.chdir("case-"+str(case))
		#f = open('../Data/Simulations/sim_data_case-'+str(case)+'.lst','w')
		data_sheet,failed_sim = data_management.generate_target_value_tables(FlameMaster_Execution_location, target_list, case, fuel,sim_type)
		temp_sim_test[str(case)] = data_sheet
		f = open('../Data/Simulations/sim_data_case-'+str(case)+'.lst','w').write(data_sheet)
		g = open('../Data/Simulations/failed_sim_data_case-'+str(case)+'.lst','w').write(failed_sim)
		#f.write(data_sheet)
		#f.close()
		os.chdir(testDir)		
	sim_dict[sim_type] = temp_sim_test
	os.chdir("..")
	print("----------------\nTest Case is over\n-------------------\n")
	return testDir,manipulation_dict,sim_dict,response_surface_dict
