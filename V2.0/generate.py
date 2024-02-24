import os, json,time
import data_management
import numpy as np
import scipy as sp
import scipy.stats as stats
import matplotlib.pyplot as plt

                    
def statistics(error):
	h = list(error)
	h.sort()
	hmean = np.mean(h)
	hstd = np.std(h)
	return h,hmean,hstd
def response_surface(sensDir,testDir,optDir,optInputs,case_dir,response_surface_dict,target_list,activeParameters,unsrt_data,stats_,selectedParams,activeIndexDict):
	case_max_error = {}
	case_mean_test_error = {}
	case_max_train_error = {}
	case_mean_train_error = {}
	temp_rs_opt = {}
	string_error = "#case/sim	Traning 	Testing\n"
	string_mean = "#case/sim	Traning 	Testing\n"
	for k,case in enumerate(target_list):
		start = time.time()
		bounds,init_guess = case.create_response_surface(activeParameters,unsrt_data,stats_["Order_of_PRS"],selectedParams,activeIndexDict,PRS_type=optInputs["Stats"]["PRS_type"])
		#print(f"Bounds = {bounds}")
		#New = case.create_response_surface(activeParameters,unsrt_data,3)
		print(f"\n\t\tTime taken to generate the response surface is {time.time()-start}")
		yData_training = case.ydata
		
		#yData_training_new = New.ydata
		
		Sim_value_training = case.resFramWrk
		#Sim_value_training_new = New.resFramWrk
		
		error_training = [abs(yData_training[i] - Sim_value_training[i]) for i in range(len(yData_training))]
		#error_training_new = [abs(yData_training_new[i] - Sim_value_training_new[i]) for i in range(len(yData_training_new))]
		
		temp_rs_opt[str(k)] = case.resCoef
		print("Testing the response surface for case-{}\n".format(target_list.index(case)))
		xData_testing, yData_testing = data_management.getTestingData(testDir,case_dir[target_list.index(case)])
		#print(f"xTest = {xData_testing}\n\ny_data = {yData_testing}")
		actual_error_testing = []
		error_testing = []
		error_testing_relative = []
		c = []
		Sim_value_testing = []
		res_Surface_Testing = open("./Data/Res_testing_case-"+str(target_list.index(case))+".txt","+w")
		string = ""
		#print(yData_testing)
		error_residual = []
		#,stats_["Order_of_PRS"]
		for index,data in enumerate(yData_testing):
			c.append(np.log(case.std_dvtn/1000))
			actual_error_testing.append(((case.evaluate(np.asarray(xData_testing[index]),stats_["Order_of_PRS"])-data)/data)*100)
			temp_relative = abs((case.evaluate(np.asarray(xData_testing[index]),stats_["Order_of_PRS"])-data)/data)*100
			actual_value = case.evaluate(np.asarray(xData_testing[index]),stats_["Order_of_PRS"])
			Sim_value_testing.append(actual_value)
			error_testing.append(abs((case.evaluate(np.asarray(xData_testing[index]),stats_["Order_of_PRS"])-data)))
			error_testing_relative.append((abs((case.evaluate(np.asarray(xData_testing[index]),stats_["Order_of_PRS"])-data))/data)*100)
			#error_residual.append((abs((case.evaluate(np.asarray(xData_testing[index]),stats_["Order_of_PRS"])-data))))
			string +="{},{},{}\n".format(data,actual_value,temp_relative)
		#print(Sim_value_testing)
		#print(np.exp(np.asarray(yData)))
		#print(c)
		c = np.asarray(c)
		x_ln = np.linspace(0.001,5,len(yData_testing))
		LN_x = np.log(np.linspace(0.001,5,len(yData_testing)))
		LN_domain = np.asarray(yData_testing)
		
		#print(yData_testing)
		res_Surface_Testing.write(string)
		MaxError_testing = max(error_testing)
		MeanError_testing = sum(error_testing)/len(error_testing)
		
		MaxError_testing_relative = max(error_testing_relative)
		MeanError_testing_relative = sum(error_testing_relative)/len(error_testing_relative)
		
		#MaxError_testing_relative = max(error_testing_relative)
		#MeanError_testing_relative = sum(error_testing_relative)/len(error_testing_relative)
		
		"""
		record the max percentage error for each cases
		"""
		case_max_error[k] = max(actual_error_testing)
		case_mean_test_error[k] = MeanError_testing_relative
		
		#MaxError_testing_new = max(error_testing_new)
		#MeanError_testing_new = sum(error_testing_new)/len(error_testing_new)
		"""
		Generating error, training samples and testing samples distribution plot
		"""
		h_train,hmean_train,hstd_train = statistics(list(yData_training))
		pdf_train = stats.norm.pdf(h_train, hmean_train, hstd_train)
		
		#print(hmean_train)
		h_test,hmean_test,hstd_test = statistics(list(yData_testing))
		pdf_test = stats.norm.pdf(h_test, hmean_test, hstd_test)
		
		
		"""
		Plotting the experimental pdf for the target
		"""
		if "Fls" in case.target:
			data_exp = np.random.normal(np.log(case.observed),np.log(case.std_dvtn),1000)
		elif "Tig" in case.target:
			data_exp = np.random.normal(np.log(case.observed*10),np.log(case.std_dvtn*10),1000)
		
		h_exp,h_exp_mean,h_exp_std = statistics(data_exp)
		
		pdf_exp = stats.norm.pdf(h_exp,h_exp_mean,h_exp_std)
		
		fig = plt.figure()
		plt.plot(h_train, pdf_train,"b-",label="Training pdf")
		plt.plot(h_test, pdf_test,"r-",label="Testing pdf")
		plt.axvline(x = hmean_train, color = 'b', label = 'training mean')
		plt.axvline(x = hmean_test, color = 'r', label = 'testing mean')		
		#plt.plot(h_exp,pdf_exp,label="Experimental pdf")
		#plt.fill_between(h_exp,pdf_exp,color='grey')
		#plt.hist(h_test,bins=15)
		#plt.plot(h_train, hmean_train,label="Training pdf")
		#plt.plot(h_test, hmean_test,label="Testing pdf")
		plt.legend()
		plt.ylabel("Probability distribution function (PDF)")
		
		if "Fls" in case.target:	
			plt.xlabel(r"Laminar burning velocities (ln{$Su_0$})")
		elif "Tig" in case.target:
			plt.xlabel(r"Ignition delay times (ln{$\tau$})")
			
		plt.savefig('Plots/Dist_training_testing_case_'+str(target_list.index(case))+'.png',bbox_inches="tight")
		string_error += "{}/{}\t{}\t{}\n".format(case_dir[target_list.index(case)],error_testing.index(MaxError_testing),case.MaxError*100,MaxError_testing)
		string_mean += "{}/{}\t{}\t{}\n".format(case_dir[target_list.index(case)],"N/A",case.MeanError*100,MeanError_testing)
		
		plt.xlim(min(h_train)*0.9,max(h_train)*1.1)	
		
		h,hmean,hstd = statistics(list(actual_error_testing))
		pdf = stats.norm.pdf(h, hmean, hstd)
		#print(max(pdf))
		#y = np.arange(0,max(pdf),100)
		#x = hmean*np.ones(len(y))
		
		fig = plt.figure()
		plt.plot(h, pdf,label="pdf")
		plt.axvline(x = hmean, color = 'b', label = 'error mean')
		#plt.plot(x,y,"--",label= "Mean")
		plt.legend()
		plt.savefig('Plots/Dist_error_case_'+str(target_list.index(case))+'.png',bbox_inches="tight")
		
		#string_error += "{}/{}\t{}\t{}\n".format(case_dir[target_list.index(case)],error_testing.index(MaxError_testing),case.MaxError*100,MaxError_testing)
		#string_mean += "{}/{}\t{}\t{}\n".format(case_dir[target_list.index(case)],"N/A",case.MeanError*100,MeanError_testing)
		
		
		
		fig = plt.figure()
		ax = fig.add_subplot()
		plt.xlabel("Response Surface estimation")
		plt.ylabel("FlameMaster Simulation")
		#plt.title("Parity plot: \nshowing error distribution in Response Surface\n testing and traning")
		if "Tig" in case.target:
			plt.xlabel(r"Black box simulations ($log(\tau)$)")
			plt.ylabel(r"Surrogate model estimation ($log(\tau)$)")
		elif "Fls" in case.target:
			plt.xlabel(r"Black box simulations ($Su^o$)")
			plt.ylabel(r"Surrogate model estimation ($Su^o$)")
		else:
			plt.xlabel(r"Black box simulations ($log(\tau)$)")
			plt.ylabel(r"Surrogate model estimation ($log(\tau)$)")
	#	plt.plot(np.asarray(yData_testing),np.asarray(Sim_value_testing),"b1",ms=4,label="scatter")
		alpha = 0.3
		plt.scatter(np.asarray(yData_testing), np.asarray(Sim_value_testing), color="none", edgecolor="blue")
			#plt.fill_between(np.asarray(yData_testing), np.asarray(Sim_value_testing)-np.percentile(np.asarray(Sim_value_testing), 0.05),np.asarray(Sim_value_testing)+np.percentile(np.asarray(Sim_value_testing), 0.05), color="k", alpha = alpha)
		#plt.errorbar(np.asarray(yData_testing),np.asarray(Sim_value),yerr=c, linestyle="None")
		#plt.text(0.95,0.01,'max error = {}\nmean error = {}'.format(MaxError,MeanError), fontsize = 5)
		ax.text(0.95, 0.01, f'Max error = {MaxError_testing}, mean error = {MeanError_testing}',
       		verticalalignment='bottom', 
       		horizontalalignment='right',
        		transform=ax.transAxes,
        		color='green', fontsize=8)
		
		x = np.linspace(0,500,1000)
		plt.xlim(min(np.asarray(Sim_value_testing))*0.98,max(np.asarray(Sim_value_testing))*1.02)
		plt.ylim(min(np.asarray(Sim_value_testing))*0.98,max(np.asarray(Sim_value_testing))*1.02)
		plt.plot(x,x,"-",label="parity line")
		plt.legend()
		plt.savefig('Plots/Parity_plot_case_'+str(target_list.index(case))+'_testing.png',bbox_inches="tight")
		
		fig = plt.figure()
		ax = fig.add_subplot()
		plt.xlabel("Response Surface estimation")
		plt.ylabel("FlameMaster Simulation")
		#plt.title("Parity plot: \nshowing error distribution in Response Surface\n testing and traning")
		if "Tig" in case.target:
			plt.xlabel(r"Black box simulations / $ln(\tau)$")
			plt.ylabel(r"Surrogate model estimation / $ln(\tau)$")
			#plt.plot(np.asarray(yData_training),np.asarray(Sim_value_training),"b1",ms=4,label="scatter")
		elif "Fls" in case.target:
			plt.xlabel(r"Black box simulations / $ln(Su^o)$")
			plt.ylabel(r"Surrogate model estimation / $ln(Su^o)$")
			#plt.plot(np.exp(np.asarray(yData_training)),np.exp(np.asarray(Sim_value_training)),"b1",ms=4,label="scatter")
		else:
			plt.xlabel(r"Black box simulations / $ln(\tau)$")
			plt.ylabel(r"Surrogate model estimation / $ln(\tau)$")
			#plt.plot(np.asarray(yData_training),np.asarray(Sim_value_training),"b1",ms=4,label="scatter")
	#	plt.plot(np.asarray(yData_training),np.asarray(Sim_value_training),"bs",ms=6,label="scatter")
		
		plt.plot(np.asarray(yData_testing),np.asarray(Sim_value_testing),"k.",ms=8,label=f"Testing (max error = {MaxError_testing_relative:.3f}%, mean error = {MeanError_testing_relative:.3f}%)")
		#plt.scatter(np.asarray(yData_testing), np.asarray(Sim_value_testing), color="none", edgecolor="black")
		plt.scatter(np.asarray(yData_training), np.asarray(Sim_value_training), color="none", edgecolor="green",label=f"Training (max error = {case.MaxError_relative:.3f}%, mean error = {case.MeanError_relative:.3f}%)")
		
		case_max_train_error[k] = case.MaxError_relative
		case_mean_train_error[k] = case.MeanError_relative
		
		#plt.errorbar(np.asarray(yData),np.asarray(Sim_value),yerr=c, linestyle="None")
		#plt.text(0.95,0.01,'max error = {}\nmean error = {}'.format(MaxError,MeanError), fontsize = 5)
		#ax.text(0.95, 0.01, f'Max error = {case.MaxError}, mean error = {case.MeanError}',
       		#verticalalignment='bottom', 
       		#horizontalalignment='right',
        		#transform=ax.transAxes,
        		#color='green', fontsize=8)
		
		x = np.linspace(0,500,1000)
		plt.xlim(min(np.asarray(yData_training))*0.98,max(np.asarray(yData_training))*1.02)
		plt.ylim(min(np.asarray(Sim_value_training))*0.98,max(np.asarray(Sim_value_training))*1.02)
		plt.plot(x,x,"-",label="parity line")
		plt.legend(loc="upper left")
		plt.savefig('Plots/Parity_plot_case_'+str(target_list.index(case))+'_training.png',bbox_inches="tight")
		
		
		
		fig = plt.figure()
		ax = fig.add_subplot()
		plt.xlabel("Response Surface estimation")
		plt.ylabel("FlameMaster Simulation")
		#plt.title("Parity plot: \nshowing error distribution in Response Surface\n testing and traning")
		if "Tig" in case.target:
			plt.xlabel(r"Black box simulations ($log(\tau)$)")
			plt.ylabel(r"Surrogate model estimation ($log(\tau)$)")
			#plt.plot(np.asarray(yData_training),np.asarray(Sim_value_training),"b1",ms=4,label="scatter")
		elif "Fls" in case.target:
			plt.xlabel(r"Black box simulations ($Su^o$)")
			plt.ylabel(r"Surrogate model estimation ($Su^o$)")
			#plt.plot(np.exp(np.asarray(yData_training)),np.exp(np.asarray(Sim_value_training)),"b1",ms=4,label="scatter")
		else:
			plt.xlabel(r"Black box simulations ($log(\tau)$)")
			plt.ylabel(r"Surrogate model estimation ($log(\tau)$)")
			#plt.plot(np.asarray(yData_training),np.asarray(Sim_value_training),"b1",ms=4,label="scatter")
	#	plt.plot(np.asarray(yData_training),np.asarray(Sim_value_training),"bs",ms=6,label="scatter")
		plt.scatter(np.asarray(yData_training), np.asarray(Sim_value_training), color="none", edgecolor="blue")
		plt.scatter(np.asarray(yData_training), np.asarray(Sim_value_training), color="none", edgecolor="blue")
		#plt.errorbar(np.asarray(yData),np.asarray(Sim_value),yerr=c, linestyle="None")
		#plt.text(0.95,0.01,'max error = {}\nmean error = {}'.format(MaxError,MeanError), fontsize = 5)
		ax.text(0.95, 0.01, f'Max error = {case.MaxError}, mean error = {case.MeanError}',
       		verticalalignment='bottom', 
       		horizontalalignment='right',
        		transform=ax.transAxes,
        		color='green', fontsize=8)
		
		x = np.linspace(0,500,1000)
		plt.xlim(min(np.asarray(yData_training))*0.98,max(np.asarray(yData_training))*1.02)
		plt.ylim(min(np.asarray(Sim_value_training))*0.98,max(np.asarray(Sim_value_training))*1.02)
		plt.plot(x,x,"-",label="parity line")
		plt.legend()
		plt.savefig('Plots/Parity_production_plot_case_'+str(target_list.index(case))+'_training.png',bbox_inches="tight")
			
		"""
		fig = plt.figure()
		plt.xlabel("Response Surface curve")
		plt.ylabel("FlameMaster Simulation")
		plt.title("Parity plot: \nshowing error distribution in Response Surface\n testing and traning")
		plt.xlabel("FlameMaster Simulations")
		plt.ylabel("Response Surface estimation")
		plt.plot(x_ln,LN_x,"--")
		plt.plot(x_ln,LN_domain,"-")
		
		plt.legend()
		plt.savefig('Plots/LN_plot_case_'+str(target_list.index(case))+'.png')
		
		"""
		"""
		LN_data = open("LN_file.csv","+w")
		string = ""
		for ind,y in enumerate(yData):
			string += "{}\t{}\n".format(y,Sim_value[ind])
		LN_data.write(string)
		LN_data.close()
		"""
	
	
	
	
	response_surface_dict["Opt"] = temp_rs_opt
	error_file = open("./Max_error.txt","w")
	error_file.write(string_error)	
	error_file.close()
	mean_file = open("./Mean_error.txt","w")
	mean_file.write(string_mean)	
	mean_file.close()
	"""
	generate a dict for cases that will be optimized using response surface and cases that will be run using direct simulations
	"""
	#max_testing_error_criteria = 2 #percentage
	max_testing_error_criteria = 10 #percentage
	
	case_prs_selection_dict = {}
	#print(case_max_error)
	count = 0
	for i in case_max_error:
		if case_max_error[i] > max_testing_error_criteria:
			case_prs_selection_dict[i] = 0
			count+=1
		else:
			case_prs_selection_dict[i] = 1
	
	#print(os.getcwd())
	if "Rejected_PRS" in os.listdir():
		fil_prs = open("Rejected_PRS","r").readlines()
		for i in fil_prs:
			case_prs_selection_dict[int(i)] = 0
	#print(case_prs_selection_dict)
	#raise AssertionError("Stop!")
	#print(k,count)
	os.chdir("..")
	
	print_selection = ""
	print_errors = "#Case\tMax_error\tMean_error\tMax_training\tMean_training\n"
	
	for i,ele in enumerate(case_prs_selection_dict):
		print_selection+=f"{i}\t{case_prs_selection_dict[ele]}\n"
	file_selection = open("selected_PRS.csv","w+").write(print_selection)
	
	
	for i,ele in enumerate(case_max_error):
		print_selection+=f"{i}\t{case_max_error[ele]}\t{case_mean_test_error}\t{case_max_train_error}\t{case_mean_train_error}\n"
	file_error = open("error_prs_PRS.csv","w+").write(print_selection)
	
	
	
	return case_prs_selection_dict,bounds,init_guess

