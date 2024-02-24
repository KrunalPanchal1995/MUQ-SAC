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
def response_surface(sensDir,testDir,optDir,case_dir,response_surface_dict,target_list,activeParameters,unsrt_data):
	temp_rs_opt = {}
	string_error = "#case/sim	Traning 	Testing\n"
	string_mean = "#case/sim	Traning 	Testing\n"
	for k,case in enumerate(target_list):
		start = time.time()
		case.create_response_surface(activeParameters,unsrt_data,2)
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
		xData_testing, yData_testing = data_management.getTestingData(sensDir,case_dir[target_list.index(case)])
		actual_error_testing = []
		error_testing = []
		c = []
		Sim_value_testing = []
		res_Surface_Testing = open("./Data/Res_testing_case-"+str(target_list.index(case))+".txt","+w")
		string = ""
		for index,data in enumerate(yData_testing):
			c.append(np.log(case.std_dvtn/1000))
			actual_error_testing.append(((case.calculated_target_value(np.asarray(xData_testing[index]))-data)/data)*100)
			temp_relative = abs((case.calculated_target_value(np.asarray(xData_testing[index]))-data)/data)*100
			actual_value = case.calculated_target_value(np.asarray(xData_testing[index]))
			Sim_value_testing.append(actual_value)
			error_testing.append(abs((case.calculated_target_value(np.asarray(xData_testing[index]))-data)))
			string +="{},{},{}\n".format(data,actual_value,temp_relative)
			
		#print(np.exp(np.asarray(yData)))
		#print(c)
		c = np.asarray(c)
		LN_x = np.linspace(0.001,5,len(yData_testing))
		LN_domain = np.log(np.asarray(yData_testing))
		#print(yData_testing)
		res_Surface_Testing.write(string)
		MaxError_testing = max(error_testing)
		MeanError_testing = sum(error_testing)/len(error_testing)
		#MaxError_testing_new = max(error_testing_new)
		#MeanError_testing_new = sum(error_testing_new)/len(error_testing_new)
		"""
		Generating error, training samples and testing samples distribution plot
		"""
		h_train,hmean_train,hstd_train = statistics(list(yData_training))
		pdf_train = stats.norm.pdf(h_train, hmean_train, hstd_train)
		
		print(hmean_train)
		h_test,hmean_test,hstd_test = statistics(list(yData_testing))
		pdf_test = stats.norm.pdf(h_test, hmean_test, hstd_test)
		
		
		fig = plt.figure()
		plt.plot(h_train, pdf_train,label="Training pdf")
		plt.plot(h_test, pdf_test,label="Testing pdf")
		plt.legend()
		plt.savefig('Plots/Dist_training_testing_case_'+str(target_list.index(case))+'.png',bbox_inches="tight")
		string_error += "{}/{}\t{}\t{}\n".format(case_dir[target_list.index(case)],error_testing.index(MaxError_testing),case.MaxError*100,MaxError_testing)
		string_mean += "{}/{}\t{}\t{}\n".format(case_dir[target_list.index(case)],"N/A",case.MeanError*100,MeanError_testing)
			
		
		h,hmean,hstd = statistics(list(actual_error_testing))
		pdf = stats.norm.pdf(h, hmean, hstd)
		fig = plt.figure()
		plt.plot(h, pdf,label="pdf")
		plt.legend()
		plt.savefig('Plots/Dist_error_case_'+str(target_list.index(case))+'.png',bbox_inches="tight")
		
		#string_error += "{}/{}\t{}\t{}\n".format(case_dir[target_list.index(case)],error_testing.index(MaxError_testing),case.MaxError*100,MaxError_testing)
		#string_mean += "{}/{}\t{}\t{}\n".format(case_dir[target_list.index(case)],"N/A",case.MeanError*100,MeanError_testing)
		
		"""
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
		
	
		plt.plot(np.asarray(yData_training_new),np.asarray(Sim_value_training_new),"b1",ms=4,label="scatter")
		#plt.errorbar(np.asarray(yData_testing),np.asarray(Sim_value),yerr=c, linestyle="None")
		#plt.text(0.95,0.01,'max error = {}\nmean error = {}'.format(MaxError,MeanError), fontsize = 5)
		ax.text(0.95, 0.01, f'Max error = {MaxError_training_new}, mean error = {MeanError_training_new}',
       		verticalalignment='bottom', 
       		horizontalalignment='right',
        		transform=ax.transAxes,
        		color='green', fontsize=8)
		
		x = np.linspace(0,500,1000)
		plt.xlim(min(np.asarray(yData_training_new))*0.98,max(np.asarray(Sim_value_training_new))*1.02)
		plt.ylim(min(np.asarray(Sim_value_training_new))*0.98,max(np.asarray(Sim_value_training_new))*1.02)
		plt.plot(x,x,"-",label="parity line")
		plt.legend()
		plt.savefig('Plots/Parity_plot_case_'+str(target_list.index(case))+'_training_3rd_order.png',bbox_inches="tight")
		
		"""
		
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
		plt.plot(np.asarray(yData_testing),np.asarray(Sim_value_testing),"b1",ms=4,label="scatter")
		#plt.errorbar(np.asarray(yData_testing),np.asarray(Sim_value),yerr=c, linestyle="None")
		#plt.text(0.95,0.01,'max error = {}\nmean error = {}'.format(MaxError,MeanError), fontsize = 5)
		ax.text(0.95, 0.01, f'Max error = {MaxError_testing}, mean error = {MeanError_testing}',
       		verticalalignment='bottom', 
       		horizontalalignment='right',
        		transform=ax.transAxes,
        		color='green', fontsize=8)
		
		x = np.linspace(0,500,1000)
		plt.xlim(min(np.asarray(yData_testing))*0.98,max(np.asarray(Sim_value_testing))*1.02)
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
		plt.plot(np.asarray(yData_training),np.asarray(Sim_value_training),"b1",ms=4,label="scatter")
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
		plt.savefig('Plots/Parity_plot_case_'+str(target_list.index(case))+'_training.png',bbox_inches="tight")
			
		#fig = plt.figure()
		#plt.xlabel("Response Surface curve")
		#plt.ylabel("FlameMaster Simulation")
		#plt.title("Parity plot: \nshowing error distribution in Response Surface\n testing and traning")
		#plt.xlabel("FlameMaster Simulations")
		#plt.ylabel("Response Surface estimation")
		#plt.plot(LN_x,LN_domain,"-")
		#plt.plot(np.asarray(yData),LN_domain,"-")
		#plt.plot(np.asarray(yData),np.asarray(Sim_value),"go",label="Simulations")
		#if "LN_file.csv" in os.listdir("Plots"):
			#file_LN = open("Plots/LN_file.csv","r").readlines()
			#x_ln = []
			#y_ln = []
			#for i in file_LN:
				#data = i.split("	")
				#x_ln.append(float(data[0].strip()))
				#y_ln.append(float(data[1].strip()))
			#plt.plot(x_ln,y_ln,":")
		#plt.plot()
		#plt.plot(np.asarray(yData),np.asarray(Sim_value),"go",label="Simulations")
		#x = np.linspace(-50,50,1000)
		#plt.xlim(np.min(np.asarray(yData))*0.999,np.max(np.asarray(yData))*1.001)
		#plt.ylim(np.min(np.asarray(Sim_value))*0.999,np.max(np.asarray(Sim_value))*1.001)
		#plt.plot(x,x,"-",label="Parity Line")
		
		#plt.legend()
		#plt.savefig('Plots/LN_plot_case_'+str(target_list.index(case))+'.png')
		
		
		#LN_data = open("LN_file.csv","+w")
		#string = ""
		#for ind,y in enumerate(yData):
		#	string += "{}\t{}\n".format(y,Sim_value[ind])
		#LN_data.write(string)
		#LN_data.close()
		
	response_surface_dict["Opt"] = temp_rs_opt
	error_file = open("./Max_error.txt","w")
	error_file.write(string_error)	
	error_file.close()
	mean_file = open("./Mean_error.txt","w")
	mean_file.write(string_mean)	
	mean_file.close()
	os.chdir("..")

