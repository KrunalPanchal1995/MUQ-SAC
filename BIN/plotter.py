import matplotlib as mpl
import matplotlib.pyplot as plt ; plt.rcdefaults()
import matplotlib.font_manager as font_manager
import numpy as np
import scipy as sp
import Uncertainty as uncertainty
from collections import OrderedDict
import re,os
import MechanismManipulator
import Input_file_reader as IFR
import glob
import sys
import warnings

if not sys.warnoptions:
    warnings.simplefilter("ignore")

site = ('Python', 'C++', 'Java', 'Perl', 'Scala', 'Lisp', 'จาวาสคลิป')
usage = [10,8,6,4,2,1, 2]
mpl.use('Agg')
#path = os.path.join(mpl.rcParams["datapath"], "fonts/ttf/DejaVuSans.ttf")
#prop = font_manager.FontProperties(fname=path)
#plt.rcParams['font.family'] = prop.get_name()


def plot_errors(target_list,opt_x):
	x = np.arange(len(target_list))
	y = []
	z = []
	cases = []
	for i in target_list:
		y.append(abs(i.calculated - i.observed)/i.observed*100)
		z.append(abs(i.calculated_target_value(opt_x) - i.observed)/i.observed*100)
		cases.append(i.case_index)
		
	fig = plt.figure()
	plt.title("Errors before and after optimization")
	plt.xlabel("Target cases")
	plt.ylabel("Errors (%)")
	sub1 = plt.subplot()

	plt.subplots_adjust(left=0.15, bottom=0.2, right=0.95, top=0.8, wspace=0, hspace=0)
	plt.xticks(x,cases,rotation='vertical')
	sub1.bar(x-0.1, y, width=0.2,color='b',align='center',label="Old error")
	sub1.bar(x+0.1, z, width=0.2,color='r', align='center',label="New Error")
	plt.legend(loc='upper left', shadow=True)
	plt.rc('axes', unicode_minus=False)
	plt.savefig("Error_plot.eps",format="eps", figsize=(20,10), dpi = 1200)
	plt.close()

def plot_vector(reaction_index, opt_x):
	x = np.arange(len(reaction_index))
	fig = plt.figure()
	plt.title("Normalized Optimum Vector")
	plt.xlabel("Reaction index")
	plt.ylabel("Perturbed rate constant (Normalized)")
	sub1 = plt.subplot()
	
	plt.subplots_adjust(left=0.15, bottom=0.2, right=0.95, top=0.8, wspace=0, hspace=0)
	plt.xticks(x,reaction_index,rotation='vertical')
	plt.xlim(-1,len(reaction_index))
	sub1.bar(x, opt_x,width=0.2, color='b', align='center')
	plt.savefig("Optimum_vector.eps",format="eps", figsize=(20,10), dpi = 1200)
	plt.rc('axes', unicode_minus=False)
	plt.close()

def extract_block_diag(A,M,k=0):
	"""https://stackoverflow.com/questions/10831417/extracting-diagonal-blocks-from-a-numpy-array"""
	"""Extracts blocks of size M from the kth diagonal
	of square matrix A, whose size must be a multiple of M."""

	# Check that the matrix can be block divided
	if A.shape[0] != A.shape[1] or A.shape[0] % M != 0:
		print(A)
		raise ValueError('Matrix must be square and a multiple of block size')

	# Assign indices for offset from main diagonal
	if abs(k) > M - 1:
		raise StandardError('kth diagonal does not exist in matrix')
	elif k > 0:
		ro = 0
		co = abs(k)*M 
	elif k < 0:
		ro = abs(k)*M
		co = 0
	else:
		ro = 0
		co = 0

	blocks = np.array([A[i+ro:i+ro+M,i+co:i+co+M] for i in range(0,len(A)-abs(k)*M,M)])
	return blocks

##To plot the posterior curves

def plot(parameters,parameter_list,case_dir,target_sets):
	os.chdir("./Plots")
	#print(os.listdir())
	if "reaction" not in os.listdir():
		os.mkdir("reaction")
	if "thermo" not in os.listdir():
		os.mkdir("thermo")
	if "transport" not in os.listdir():
		os.mkdir("transport")
	if "collisionEff" not in os.listdir():
		os.mkdir("collisionEff")
	if "fallOffCurve" not in os.listdir():
		os.mkdir("fallOffCurve")
	
	for i in os.listdir():
		os.chdir(i)
		if "uncert_func" not in os.listdir():
			os.mkdir("uncert_func")
			os.mkdir("uncert_limits")
			os.mkdir("optimized")
			os.chdir("..")
		else:
			os.chdir("..")
			continue
	
	os.chdir("..")
	if "BranchRatio" not in os.listdir():
		os.mkdir("BranchRatio")
#	os.chdir("thermo")
#	os.mkdir("uncert_func")
#	os.mkdir("uncert_limits")
#	os.mkdir("optimized")
#	os.chdir("..")
#	
#	os.chdir("transport")
#	os.mkdir("uncert_func")
#	os.mkdir("uncert_limits")
#	os.mkdir("optimized")
#	os.chdir("../..")
	for index,params in enumerate(parameter_list):
		plot_uncertFunc(params)
		plot_uncertLimits(params)
		plot_optimized(params)

	#for case,dataset in enumerate(target_sets):
	#	plot_optResults(dataset)

def plot_optResults(index,add,dataSet,opt_tag="on"):
	#if "plotx" in add:
		#print(add["plotx"])
	
	if dataSet.dataType:
		fig = plt.figure()
		model_response = []
		model_nominal = []
		if "Tig" in dataSet.dataType[0].strip():
			string = "1000/T\K^-1\tobs\tstd\tNominal\tOpt\n"
			plt.xlabel('1000/T ($K^{-1}$)')
			plt.ylabel(r'Ignition Delay, $\tau$ / $\mu$s')
			plt.yscale('log')
			#plt.yticks(np.geomspace(min(np.exp(np.asarray(dataSet.posterior_model_response))/10),max(np.exp(np.asarray(dataSet.posterior_model_response))/10),5))
			#plt.ylim(0,10**5)
			x = 1000.0/np.asarray(dataSet.Temperature)
			plt.plot(x,np.exp(np.asarray(dataSet.nominalSimulations))/10,'r-.',linewidth=0.8,label = "Prior estimation")
			model_nominal.extend(np.exp(np.asarray(dataSet.nominalSimulations))/10)
			if opt_tag == "on":
				plt.plot(x,np.exp(np.asarray(dataSet.posterior_model_response))/10,'k-',label = "Current work")
				model_response.extend(list(np.exp(np.asarray(dataSet.posterior_model_response))/10))
				plt.plot(x,np.exp(np.asarray(dataSet.prior_model_response)),'r-',linewidth=0.8,label = "Prior estimation")
				error_prior = np.exp(np.asarray(dataSet.prior_model_unsrt))
				#error_posterior = np.exp(np.asarray(dataSet.posterior_model_unsrt))
				#print(error_posterior)
				#plt.fill_between(x,np.exp(np.asarray(dataSet.posterior_model_response))+ error_posterior, np.exp(np.asarray(dataSet.posterior_model_response))-error_posterior,label="Posterior uncertainty")
		if "Fls" in dataSet.dataType[0].strip():
			#print(dataSet.fls_id)
			if "plotx" in add:
				string = "P[atm]\tobs\tstd\tNominal\tOpt\n"
				if add["plotx"] == 'P':
					plt.xlabel('Pressure / atm')
					x = np.asarray(dataSet.Pressure)/101325
					#plt.title(r"$\phi$ = {}".format(dataSet.fls_x[0]))
					#print(x)
			else:
				string = "Phi\tobs\tstd\tNominal\tOpt\n"
				plt.xlabel('{}'.format(dataSet.fls_id[0]))
				x = np.asarray(dataSet.fls_x)
				#plt.title("Pressure = {} atm".format(dataSet.Pressure[0]/101325))
			plt.plot(x,np.exp(np.asarray(dataSet.nominalSimulations)),'r-.',linewidth=0.8,label = "Prior estimation")
			model_nominal.extend(np.exp(np.asarray(dataSet.nominalSimulations)))
			if opt_tag == "on":
				
				plt.plot(x,np.exp(np.asarray(dataSet.posterior_model_response)),'k-',label = "Current work")
				plt.plot(x,np.exp(np.asarray(dataSet.prior_model_response)),'r-',linewidth=0.8,label = "Prior estimation")
				model_response.extend(list(np.exp(np.asarray(dataSet.posterior_model_response))))
				error_prior = np.exp(np.asarray(dataSet.prior_model_unsrt))
				#error_posterior = np.exp(np.asarray(dataSet.posterior_model_unsrt))
				#print(error_posterior)
				#plt.fill_between(x,np.exp(np.asarray(dataSet.posterior_model_response))+ error_posterior, np.exp(np.asarray(dataSet.posterior_model_response))-error_posterior,label="Posterior uncertainty")
			plt.ylabel(r'Laminar Burning Velocities, $S_u$ / $cm\,s^{-1}$')
			
		if "Flw" in dataSet.dataType[0].strip():
			plt.xlabel('1000/T ($\[K^{-1}])')
			plt.ylabel(r'Species concentration [$ppm\,ms^{-1}$]')
			x = 1000.0/np.asarray(dataSet.Temperature)
			if opt_tag == "on":
				plt.plot(x,np.exp(np.asarray(dataSet.posterior_model_response)),'k-',label = "Current work")
				model_response.extend(list(np.exp(np.asarray(dataSet.posterior_model_response))))
				#plt.plot(x,np.exp(np.asarray(dataSet.prior_model_response)),'r-',linewidth=0.8,label = "Prior estimation")
				error_prior = np.exp(np.asarray(dataSet.prior_model_unsrt))
				#error_posterior = np.exp(np.asarray(dataSet.posterior_model_unsrt))
				#print(error_posterior)
				#plt.fill_between(x,np.exp(np.asarray(dataSet.posterior_model_response))+ error_posterior, np.exp(np.asarray(dataSet.posterior_model_response))-error_posterior,label="Posterior uncertainty")
		if "Flf" in dataSet.dataType[0].strip():
			plt.xlabel('1000/T ($\[K^{-1}])')
			plt.ylabel(r'Species concentration $X_{%(species)s,%(cond)s}$ [$ mol$]'%{"species":dataSet.addendum[0]["flf_target"],"cond": dataSet.addendum[0]["flf_cond"]})
			plt.ylim(0,np.max(np.asarray(dataSet.observed))*1.02)
			x = 1000.0/np.asarray(dataSet.Temperature)
			if opt_tag == "on":
				plt.plot(x,np.exp(np.asarray(dataSet.posterior_model_response)),'k-',label = "Current work")
				model_response.extend(list(np.exp(np.asarray(dataSet.posterior_model_response))))
				#plt.plot(x,np.exp(np.asarray(dataSet.prior_model_response)),'r-',linewidth=0.8,label = "Prior estimation")
				error_prior = np.exp(np.asarray(dataSet.prior_model_unsrt))
				#error_posterior = np.exp(np.asarray(dataSet.posterior_model_unsrt))
				#print(error_posterior)
				#plt.fill_between(x,np.exp(np.asarray(dataSet.posterior_model_response))+ error_posterior, np.exp(np.asarray(dataSet.posterior_model_response))-error_posterior,label="Posterior uncertainty")
			
		
		
		home = os.getcwd()
		os.chdir("./Plots")
		dir_list = os.listdir()
		if "Dataset" not in dir_list:
			os.mkdir("Dataset")
			os.chdir("Dataset")
		else:
			os.chdir("Dataset")
		
		plt.errorbar(x,(np.asarray(dataSet.observed)),yerr = (np.asarray(dataSet.Unsrt_fct)),fmt='.',ecolor='green',markersize=8,capsize=2,label = "Observed")
		print(len(x))
		for i in range(len(x)):
			string+=f"{x[i]}\t{dataSet.observed[i]}\t{dataSet.Unsrt_fct[i]}\t{model_nominal[i]}\n"
		plotter_file = open(f"{str(index)}"+".csv",'w').write(string)
			
		#print(error_prior)
		#plt.fill_between(x,np.exp(np.asarray(dataSet.prior_model_response))+ error_prior, np.exp(np.asarray(dataSet.nominalSimulations))-error_prior,label="Prior uncertainty")
		#if opt_tag == "on":
		#	plt.plot(x,np.exp(np.asarray(dataSet.posterior_model_response))/10,'k-.',label = "Current work")
			#plt.plot(x,np.exp(np.asarray(dataSet.prior_model_response)),'r-',linewidth=0.8,label = "Prior estimation")
		#	error_prior = np.exp(np.asarray(dataSet.prior_model_unsrt))
			#error_posterior = np.exp(np.asarray(dataSet.posterior_model_unsrt))
			#print(error_posterior)
			#plt.fill_between(x,np.exp(np.asarray(dataSet.posterior_model_response))+ error_posterior, np.exp(np.asarray(dataSet.posterior_model_response))-error_posterior,label="Posterior uncertainty")
		
		plt.legend()
		plt.savefig(str(index)+'.pdf',bbox_inches="tight")
		plt.close()
		os.chdir(home)
		
		
def plot_uncertFunc(variable):
	if variable.dataKey.strip() != "":
		for i in variable.dataKey.split(","):
			string = str(i)
			T = np.asarray(variable.temperatures[str(i)])
			u = np.asarray(variable.priorUnsrt[str(i)])
			u_k = np.asarray(variable.unsrtBasisFunction[str(i)])
			u_z = np.asarray(variable.unsrt_function[str(i)])
			fig = plt.figure()
			plt.yscale('log')
			plt.title('Uncertainty function for {}\t\n\t\t:{}]'.format(variable.tag,variable.name))
			plt.plot((1000/T),u,'o',label='Uncertainties (Baulch (2002))')
			plt.plot((1000/T),u_k,"-",label='Uncertainty function ($\sigma_{\kappa}$)')
			plt.plot((1000/T),u_z,"--",label='Uncertainty function (change of base ;$\sigma_{\zeta}$)')
			plt.rcParams['font.sans-serif'] = ['sans-serif']
			plt.rcParams['axes.unicode_minus'] = False
			plt.savefig("./Plots/"+str(variable.tag)+"/uncert_func/"+str(variable.name)+"_"+string+".pdf")
			plt.legend()
			plt.close()
	else:
		T = np.asarray(variable.temperatures)
		u = np.asarray(variable.priorUnsrt)
		u_k = np.asarray(variable.unsrtBasisFunction)
		u_z = np.asarray(variable.unsrt_function)
		fig = plt.figure()
		plt.yscale('log')
		plt.title('Uncertainty function for {}\t\n\t\t:{}]'.format(variable.tag,variable.name))
		plt.plot((1000/T),u,'o',label='Uncertainties (Baulch (2002))')
		plt.plot((1000/T),u_k,"-",label='Uncertainty function ($\sigma_{\kappa}$)')
		plt.plot((1000/T),u_z,"--",label='Uncertainty function (change of base ;$\sigma_{\zeta}$)')
		
		plt.rcParams['font.sans-serif'] = ['sans-serif']
		plt.rcParams['axes.unicode_minus'] = False
		plt.savefig("./Plots/"+str(variable.tag)+"/uncert_func/"+str(variable.name)+".pdf",bbbox_inches="tight")
		plt.legend()
		plt.close()

def plot_uncertLimits(variable):
	#print(variable.tag)
	if variable.dataKey.strip() != "":	
		for i in variable.dataKey.split(","):
			T = np.asarray(variable.temperatures[str(i)])
			N = np.asarray(Func(variable.tag,variable.nominal[str(i)],T,str(i)))
			#print(N)
			#print(variable.unsrt_function[str(i)])
			N_upper = N + variable.unsrt_function[str(i)]
			N_lower = N - variable.unsrt_function[str(i)]
			N_basisUpper = np.asarray(Func(variable.tag,variable.upperBasisLimit[str(i)],T,str(i)))
			N_basisLower = np.asarray(Func(variable.tag,variable.lowerBasisLimit[str(i)],T,str(i)))
			fig = plt.figure()
			plt.yscale('log')
			plt.title('Uncertainty band for {}\t\n\t\t:{}]'.format(variable.tag,variable.name))
			plt.plot((1000/T),N,"b-",label='Nominal')
			plt.plot((1000/T),N_upper,"k-",label='Upper uncertainty limit ($\sigma_{\kappa}$)',linewidth=0.75)
			plt.plot((1000/T),N_lower,"k-",label='Lower uncertainty limit ($\sigma_{\kappa}$)',linewidth=0.75)
			#plt.plot((1/T),N_basisLower,"r--",label='Upper uncertainty limit ($\sigma_{\zeta}$)')
			#plt.plot((1/T),N_basisUpper,"r--",label='Lower uncertainty limit ($\sigma_{\zeta}$)')
			
			plt.rcParams['font.sans-serif'] = ['sans-serif']
			plt.rcParams['axes.unicode_minus'] = False
			plt.savefig("./Plots/"+str(variable.tag)+"/uncert_limits/"+str(variable.name)+"_"+str(i)+".pdf")
			plt.legend()			
			plt.close()
	
	else:
	
		
		T = np.asarray(variable.temperatures)
		N = np.asarray(Func(variable.tag,variable.nominal,T))
		
		#print(N)
		#print(variable.unsrt_function)
		N_upper = N + variable.unsrt_function
		N_lower = N - variable.unsrt_function
		N_basisUpper = np.asarray(Func(variable.tag,variable.upperBasisLimit,T))
		N_basisLower = np.asarray(Func(variable.tag,variable.lowerBasisLimit,T))
		fig = plt.figure()
		plt.yscale('log')
		plt.title('Uncertainty band for {}\t\n\t\t:{}]'.format(variable.tag,variable.name))
		plt.plot((1000/T),N,"b-",label='Nominal')
		plt.plot((1000/T),N_upper,"k-",label='Upper uncertainty limit ($\sigma_{\kappa}$)',linewidth=0.75)
		plt.plot((1000/T),N_lower,"k-",label='Lower uncertainty limit ($\sigma_{\kappa}$)',linewidth=0.75)
		#plt.plot((1/T),N_basisLower,"r--",label='Upper uncertainty limit ($\sigma_{\zeta}$)')
		#plt.plot((1/T),N_basisUpper,"r--",label='Lower uncertainty limit ($\sigma_{\zeta}$)')
		plt.rcParams['font.sans-serif'] = ['sans-serif']
		plt.rcParams['axes.unicode_minus'] = False
		plt.savefig("./Plots/"+str(variable.tag)+"/uncert_limits/"+str(variable.name)+".pdf")
		plt.legend()
		plt.close()

def plot_optimized(variable):

	if variable.dataKey.strip() != "":	
		for i in variable.dataKey.split(","):
			T = np.asarray(variable.temperatures[str(i)])
			N = np.asarray(Func(variable.tag,variable.nominal[str(i)],T,str(i)))
			N_opt = np.asarray(Func(variable.tag,variable.optPerturbation[str(i)],T,str(i)))
			N_basisUpper = np.asarray(Func(variable.tag,variable.upperBasisLimit[str(i)],T,str(i)))
			N_basisLower = np.asarray(Func(variable.tag,variable.lowerBasisLimit[str(i)],T,str(i)))
			fig = plt.figure()
			plt.yscale('log')
			plt.title('Optimized parameters for {}\t\n\t\t:{}]'.format(variable.tag,variable.name))
			plt.plot((1000/T),N,"b-",label='Prior rate')
			plt.plot((1000/T),N_opt,"r-",label="Optimized rate")
			#plt.plot((1/T),N_basisLower,"k--",label='Upper uncertainty limit ($\sigma_{\zeta}$)')
			#plt.plot((1/T),N_basisUpper,"k--",label='Lower uncertainty limit ($\sigma_{\zeta}$)')
			plt.rcParams['font.sans-serif'] = ['sans-serif']
			plt.rcParams['axes.unicode_minus'] = False		
			plt.savefig("./Plots/"+str(variable.tag)+"/optimized/"+str(variable.name)+"_"+str(i)+".pdf")
			plt.legend()
			plt.close()
	
	else:
		#print("Entered")
		T = np.asarray(variable.temperatures)
		N = np.log10(np.asarray(Func(variable.tag,variable.nominal,T)))
		
		#N_opt = np.asarray(Func(variable.tag,variable.optPerturbation,T))
		N_opt = np.log10(np.asarray(Func(variable.tag,variable.optArrhenius,T)))
		N_basisUpper = np.log10(np.asarray(Func(variable.tag,variable.upperBasisLimit,T)))
		N_basisLower = np.log10(np.asarray(Func(variable.tag,variable.lowerBasisLimit,T)))
		N_opt_basisUpper = np.log10(np.asarray(Func(variable.tag,variable.optUpperBasisLimit,T)))
		N_opt_basisLower = np.log10(np.asarray(Func(variable.tag,variable.optLowerBasisLimit,T)))
		fig = plt.figure()
		#plt.yscale('log')
		#plt.yticks(np.geomspace(min(N_opt),max(N_opt),4).round())
		plt.title('Optimized parameters for {}\t\n\t\t:{}]'.format(variable.tag,variable.name))
		plt.plot((1000/T),N,"b-",label='Prior rate')
		plt.plot((1000/T),N_basisLower,"k-",label='Lower uncertainty limit, $\kappa_{max}(\sigma_{\zeta})$',linewidth=0.75)
		plt.plot((1000/T),N_basisUpper,"k-",label='Upper uncertainty limit, $\kappa_{min}(\sigma_{\zeta})$',linewidth=0.75)
		plt.plot((1000/T),N_opt,"r-",label="Optimized rate")
		#plt.plot((1/T),N_opt_basisLower,"k--",label='Opt upper uncertainty limit ($\sigma_{\zeta^*}$)')
		#plt.plot((1/T),N_opt_basisUpper,"k--",label='Opt lower uncertainty limit ($\sigma_{\zeta^*}$)')
		plt.xlabel(r"1000/T\K$^{-1}$")
		plt.ylabel(r"$log_{10}(k)$ / $s^{-1}$ or $log_{10}$(k) / $cm^{3}\,molecule^{-1}\,s^{-1}$")
		plt.legend()
		plt.rcParams['font.sans-serif'] = ['sans-serif']
		plt.rcParams['axes.unicode_minus'] = False
		plt.savefig("./Plots/"+str(variable.tag)+"/optimized/"+str(variable.name)+".pdf",bbox_inches="tight")
		plt.close()
		plt.legend()
		x = 1000/T
		string = "1000/T\K^-1\t Prior \t unsrt_max \t unsrt_min \tOpt\n"
		for i in range(len(x)):
			string+=f"{x[i]}\t{N[i]}\t{N_basisUpper[i]}\t{N_basisLower[i]}\t{N_opt[i]}\n"
		plotter_file = open("./Plots/"+str(variable.tag)+"/optimized/"+str(variable.name)+".csv",'w').write(string)		
				


def Func(tag,vector,T, string = ""):
	if tag == "reaction":
		function = np.exp(float(vector[0]))*T**float(vector[1])*np.exp(-float(vector[2])/(T))
		return function
	if tag == "fallOffCurve":
		function = float(vector)*(T/T)
		return function
	if tag == "thermo":
		if string == "Hcp":
			function = float(vector[0])*(T/T) + float(vector[1])*T + float(vector[2])*T**2  +float(vector[3])*T**3 + float(vector[4])*T**4
			return function
		if string == "h":
			function = float(vector[0])*(T/T)
			return function
		if string == "e":
			function = float(vector[0])*(T/T)
			return function
	if tag == "transport":
		#print(vector)
		function = float(vector)*(T/T)
		return function
	if tag == "collisionEff":
		#print(vector)
		function = float(vector)*(T/T)
		return function

def plot_datasets():		
		'''		
		#print(k_pmax)
		fig = plt.figure()
		#plt.title('{}: Zeta = [{} {} {}]'.format(r_,zeta.x[0],zeta.x[1],zeta.x[2]))
		plt.plot((1/T),k,label='Nominal rate constant (k)')
		plt.plot((1/T),k_max,label='Upper prior TDeptUnsrt limit (k+f)')
		plt.plot((1/T),k_min,label='Lower prior TDeptUnsrt limit (k-f)')
		plt.plot((1/T),kappa_max,'--',label='Cholesky_kappa_max')
		plt.plot((1/T),kappa_min,'--',label ='Cholesky_kappa_min')
		plt.legend()
		plt.savefig(os.getcwd()+'/Plots/Kappa_vs_T/'+r_+'.png')
		plt.close()
		
		fig = plt.figure()
		plt.plot((1/T),k,"-b",label='Nominal rate constant (k)')
		plt.plot((1/T),k_max,"--b",label='Upper prior TDeptUnsrt limit (k+f)')
		plt.plot((1/T),k_min,"--b",label='Lower prior TDeptUnsrt limit (k-f)')
		plt.plot((1/T),k_opt,'-r',label='optimized rate constant')
		plt.plot((1/T),k_pmax,'--r',label='Upper posterior TdeptUnsrt limit')
		plt.plot((1/T),k_pmin,'--r',label ='Lower posterior TdeptUnsrt limit')
		plt.legend()
		plt.savefig(os.getcwd()+'/Plots/Posterior_unsrt/'+r_+'.png')
		plt.close()
		
		plot_k = open("./Plots/PlotData/Nominal_"+r_+".txt","+a")
		plot_k_max = open("./Plots/PlotData/Prior_Max"+r_+".txt","+a")
		plot_k_min = open("./Plots/PlotData/Prior_Min"+r_+".txt","+a")
		plot_kappa_max = open("./Plots/PlotData/Prior_Max_zeta"+r_+".txt","+a")
		plot_kappa_min = open("./Plots/PlotData/Prior_Min_zeta_"+r_+".txt","+a")
		printTDeptUnsrt = open("./Plots/PlotData/TDeptUnsrt_"+r_+".txt","+a")
		
		temp_k = ""
		temp_k_max = ""
		temp_k_min = ""
		temp_kappa_max = ""
		temp_kappa_min = ""
		temp_Tdeptunsrt = ""
		
		for index,i in enumerate(T):
			temp_k += "{},{}\n".format(1/i,k[index])
			temp_k_max += "{},{}\n".format(1/i,k_max[index])
			temp_k_min += "{},{}\n".format(1/i,k_min[index])
			temp_kappa_max += "{},{}\n".format(1/i,kappa_max[index])
			temp_kappa_min +="{},{}\n".format(1/i,kappa_min[index])
			temp_Tdeptunsrt += "{},{}\n".format(1/i,TDeptUnsrt[index])	
			
		plot_k.write(temp_k) 
		plot_k_max.write(temp_k_max) 
		plot_k_min.write(temp_k_min ) 
		plot_kappa_max.write(temp_kappa_max) 
		plot_kappa_min.write(temp_kappa_min) 
		printTDeptUnsrt.write(temp_Tdeptunsrt) 
		
		plot_k.close()
		plot_k_max.close()
		plot_k_min.close()
		plot_kappa_max.close()
		plot_kappa_min.close()
		printTDeptUnsrt.close()
		
		#####for fallOffCurve ####
	for index,foc in enumerate(fallOffCurve_index):
		foc_ = foc.split("-")[0]
		choleskyMatrix = unsrt_data[foc].getCholeskyMat()
		
		#print(choleskyMatrix)
		zeta = unsrt_data[foc].getZeta()
		#print(zeta)
		T = unsrt_data[foc].getTemperature()
		TDeptUnsrt = unsrt_data[foc].getUncertainty(T)
		nominal = Input_file_reader.MechParsing(mech_file_location).getFocData(foc)		
		beta_opt = np.asarray(Beta_Dict[foc])
		
		
		alpha =np.exp(-T[0]/nominal[1]) 
		sigma = (1-nominal[0])*np.exp(-T[0]/nominal[1])
		gamma = nominal[2]*np.exp(-T[0]/nominal[3])
		delta = nominal[4]*np.exp(-nominal[5]/T[0])
		
		
		nominal_ = np.array([alpha,sigma,gamma,delta])
		
		parameter = alpha*(T/T)-sigma*(T/T)+gamma*(T/T)+delta*(T/T)
		pMin = np.array([alpha,sigma,gamma,delta]) - np.asarray(np.dot(choleskyMatrix,zeta)).flatten()
		pMax =  np.array([alpha,sigma,gamma,delta]) + np.asarray(np.dot(choleskyMatrix,zeta)).flatten()
		
		
		pExp_max = TDeptUnsrt+parameter
		pExp_min = parameter-TDeptUnsrt
		
		
		zeta_star = beta_opt*zeta
		
		zeta_pmax = zeta_star+np.dot(Posterior_cov[1][index],np.asarray([1,1,1,1]))
		zeta_pmin = zeta_star-np.dot(Posterior_cov[1][index],np.asarray([1,1,1,1]))	
		p_opt = nominal_  + np.asarray(np.dot(choleskyMatrix,zeta_star)).flatten()		
		parameter_opt = p_opt[0]*(T/T)-p_opt[1]*(T/T)+p_opt[2]*(T/T)+p_opt[3]*(T/T)
		
		p_opt_max =   np.array([alpha,sigma,gamma,delta])  + np.asarray(np.dot(choleskyMatrix,zeta_pmax)).flatten()
		p_opt_min = np.array([alpha,sigma,gamma,delta])   + np.asarray(np.dot(choleskyMatrix,zeta_pmin)).flatten()
		
		
		
		parameter_max =pMax[0]*(T/T)-pMax[1]*(T/T)+pMax[2]*(T/T)+pMax[3]*(T/T)
		parameter_min = pMin[0]*(T/T)-pMin[1]*(T/T)+pMin[2]*(T/T)+pMin[3]*(T/T)
		
		
		parameter_opt_max = p_opt_max[0]*(T/T)-p_opt_max[1]*(T/T)+p_opt_max[2]*(T/T)+p_opt_max[3]*(T/T)
		parameter_opt_min = p_opt_min[0]*(T/T)-p_opt_min[1]*(T/T)+p_opt_min[2]*(T/T)+p_opt_min[3]*(T/T) 	
		fig = plt.figure()
		
		plt.plot((1/T),parameter,label='Nominal {}'.format(foc))
		plt.plot((1/T),pExp_max,label='Upper exp unsrt')
		plt.plot((1/T),pExp_min,label='Lower exp unsrt')
		plt.plot((1/T),parameter_max,'--',label='Prior_Max')
		plt.plot((1/T),parameter_min,'--',label ='Prior_Min')
		plt.legend()
		plt.savefig(os.getcwd()+'/Plots/Kappa_vs_T/'+foc+'.png')
		plt.close()
		
		fig = plt.figure()
		plt.plot((1/T),parameter,"-b",label='Nominal rate constant (k)')
		plt.plot((1/T),parameter_max,"--b",label='Upper prior TDeptUnsrt limit (k+f)')
		plt.plot((1/T),parameter_min,"--b",label='Lower prior TDeptUnsrt limit (k-f)')
		plt.plot((1/T),parameter_opt,'-r',label='opt {}'.format(foc))
		plt.plot((1/T),parameter_opt_max,'--r',label='Upper posterior unsrt limit')
		plt.plot((1/T),parameter_opt_min,'--r',label ='Lower posterior unsrt limit')
		plt.legend()
		plt.savefig(os.getcwd()+'/Plots/Posterior_unsrt/'+foc+'.png')
		plt.close()
		plot_k = open("./Plots/PlotData/Nominal_"+foc+".txt","+a")
		plot_k_max = open("./Plots/PlotData/Prior_Exp_Max"+foc+".txt","+a")
		plot_k_min = open("./Plots/PlotData/Prior_Exp_Min"+foc+".txt","+a")
		plot_kappa_max = open("./Plots/PlotData/Prior_Max_zeta"+foc+".txt","+a")
		plot_kappa_min = open("./Plots/PlotData/Prior_Min_zeta_"+foc+".txt","+a")
		plot_p_opt = open("./Plots/PlotData/Optimal_"+foc+".txt","+a")
		plot_p_opt_max = open("./Plots/PlotData/posterior_Max_zeta"+foc+".txt","+a")
		plot_p_opt_min = open("./Plots/PlotData/posterior_Min_zeta_"+foc+".txt","+a")
		printTDeptUnsrt = open("./Plots/PlotData/TDeptUnsrt_"+foc+".txt","+a")
		
		temp_k = ""
		temp_k_max = ""
		temp_k_min = ""
		temp_kappa_max = ""
		temp_kappa_min = ""
		temp_p_opt = ""
		temp_p_opt_max = ""
		temp_p_opt_min = ""
		temp_Tdeptunsrt = ""
		
		for index,i in enumerate(T):
			temp_k += "{},{}\n".format(1/i,parameter[index])
			temp_k_max += "{},{}\n".format(1/i,pExp_max[index])
			temp_k_min += "{},{}\n".format(1/i,pExp_min[index])
			temp_kappa_max += "{},{}\n".format(1/i,parameter_max[index])
			temp_kappa_min +="{},{}\n".format(1/i,parameter_min[index])
			temp_p_opt += "{},{}\n".format(1/i,parameter_opt[index])
			temp_p_opt_max += "{},{}\n".format(1/i,parameter_opt_max[index])
			temp_p_opt_min += "{},{}\n".format(1/i,parameter_opt_min[index])
			temp_Tdeptunsrt += "{},{}\n".format(1/i,TDeptUnsrt[index])	
			
		plot_k.write(temp_k) 
		plot_k_max.write(temp_k_max) 
		plot_k_min.write(temp_k_min ) 
		plot_kappa_max.write(temp_kappa_max) 
		plot_kappa_min.write(temp_kappa_min) 
		plot_p_opt.write(temp_p_opt)
		plot_p_opt_max.write(temp_p_opt_max) 
		plot_p_opt_min.write(temp_p_opt_min) 
		printTDeptUnsrt.write(temp_Tdeptunsrt) 
		
		plot_k.close()
		plot_k_max.close()
		plot_k_min.close()
		plot_kappa_max.close()
		plot_kappa_min.close()
		plot_p_opt.close()
		plot_p_opt_max.close() 
		plot_p_opt_min.close() 
		printTDeptUnsrt.close()
	
	
		
		####for thermo ########
	for index,hcp in enumerate(heatCap_index):
		hcp_ = hcp.split("-")[0]
		hcp_type = hcp.split("-")[2]
		if hcp_type == "low":
			choleskyMatrix = unsrt_data[hcp].getCholeskyMat()
			zeta = unsrt_data[hcp].getZeta()
			T = unsrt_data[hcp].getTemperature()
			TDeptUnsrt = unsrt_data[hcp].getUncertainty(T)
			nominal = Input_file_reader.ThermoParsing(thermo_file_location).getHCP_low(hcp)
			beta_opt = np.asarray(Beta_Dict[hcp])
			tag = "low"		
		elif hcp_type == "high":
			choleskyMatrix = unsrt_data[hcp].getCholeskyMat()
			zeta = unsrt_data[hcp].getZeta()
			T = unsrt_data[hcp].getTemperature()
			TDeptUnsrt = unsrt_data[hcp].getUncertainty(T)
			nominal = Input_file_reader.ThermoParsing(thermo_file_location).getHCP_high(hcp)
			beta_opt = np.asarray(Beta_Dict[hcp])
			tag = "high"	
		else:
			print("Invalid Input for Heat Capacities")
		
		pMin = nominal  - np.asarray(np.dot(choleskyMatrix,zeta)).flatten()
		pMax = nominal  + np.asarray(np.dot(choleskyMatrix,zeta)).flatten()
			
		parameter = nominal[0]*(T/T)+nominal[1]*T+nominal[2]*T**2+nominal[3]*T**3+nominal[4]*T**4
		
		pExp_max = TDeptUnsrt+parameter
		pExp_min = parameter-TDeptUnsrt
		
		
		zeta_star = beta_opt*zeta
		zeta_pmax = zeta_star+np.dot(Posterior_cov[2][index],np.asarray([1,1,1,1,1]))
		zeta_pmin = zeta_star-np.dot(Posterior_cov[2][index],np.asarray([1,1,1,1,1]))
		p_opt = nominal  + np.asarray(np.dot(choleskyMatrix,zeta_star)).flatten()
		parameter_opt = p_opt[0]*(T/T)+p_opt[1]*T+p_opt[2]*T**2+p_opt[3]*T**3+p_opt[4]*T**4 
		p_opt_max =  nominal + np.asarray(np.dot(choleskyMatrix,zeta_pmax)).flatten()
		p_opt_min = nominal  + np.asarray(np.dot(choleskyMatrix,zeta_pmin)).flatten()
		
		
		
		
		parameter_max =pMax[0]*(T/T)+pMax[1]*T+pMax[2]*T**2+pMax[3]*T**3+pMax[4]*T**4
		parameter_min = pMin[0]*(T/T)+pMin[1]*T+pMin[2]*T**2+pMin[3]*T**3+pMin[4]*T**4
		
		
		parameter_opt_max = p_opt_max[0]*(T/T)+p_opt_max[1]*T+p_opt_max[2]*T**2+p_opt_max[3]*T**3+p_opt_max[4]*T**4 
		parameter_opt_min = p_opt_min[0]*(T/T)+p_opt_min[1]*T+p_opt_min[2]*T**2+p_opt_min[3]*T**3+p_opt_min[4]*T**4 	
		fig = plt.figure()
		
		plt.plot((1/T),parameter,label='Nominal {}-{}'.format(hcp_,tag))
		plt.plot((1/T),pExp_max,label='Upper exp unsrt')
		plt.plot((1/T),pExp_min,label='Lower exp unsrt')
		plt.plot((1/T),parameter_max,'--',label='Prior_Max')
		plt.plot((1/T),parameter_min,'--',label ='Prior_Min')
		plt.legend()
		plt.savefig(os.getcwd()+'/Plots/Kappa_vs_T/'+hcp_+":"+tag+'.png')
		plt.close()
		
		fig = plt.figure()
		plt.plot((1/T),parameter,"-b",label='Nominal rate constant (k)')
		plt.plot((1/T),parameter_max,"--b",label='Upper prior TDeptUnsrt limit (k+f)')
		plt.plot((1/T),parameter_min,"--b",label='Lower prior TDeptUnsrt limit (k-f)')
		plt.plot((1/T),parameter_opt,'-r',label='opt {}-{}'.format(hcp_,tag))
		plt.plot((1/T),parameter_opt_max,'--r',label='Upper posterior unsrt limit')
		plt.plot((1/T),parameter_opt_min,'--r',label ='Lower posterior unsrt limit')
		plt.legend()
		plt.savefig(os.getcwd()+'/Plots/Posterior_unsrt/'+hcp_+":"+tag+'.png')
		plt.close()
	
		
		plot_k = open("./Plots/PlotData/Nominal_"+hcp_+":"+tag+".txt","+a")
		plot_k_max = open("./Plots/PlotData/Prior_Exp_Max"+hcp_+":"+tag+".txt","+a")
		plot_k_min = open("./Plots/PlotData/Prior_Exp_Min"+hcp_+":"+tag+".txt","+a")
		plot_kappa_max = open("./Plots/PlotData/Prior_Max_zeta"+hcp_+":"+tag+".txt","+a")
		plot_kappa_min = open("./Plots/PlotData/Prior_Min_zeta_"+hcp_+":"+tag+".txt","+a")
		plot_p_opt = open("./Plots/PlotData/Optimal_"+hcp_+":"+tag+".txt","+a")
		plot_p_opt_max = open("./Plots/PlotData/posterior_Max_zeta"+hcp_+":"+tag+".txt","+a")
		plot_p_opt_min = open("./Plots/PlotData/posterior_Min_zeta_"+hcp_+":"+tag+".txt","+a")
		printTDeptUnsrt = open("./Plots/PlotData/TDeptUnsrt_"+hcp_+":"+tag+".txt","+a")
		
		temp_k = ""
		temp_k_max = ""
		temp_k_min = ""
		temp_kappa_max = ""
		temp_kappa_min = ""
		temp_p_opt = ""
		temp_p_opt_max = ""
		temp_p_opt_min = ""
		temp_Tdeptunsrt = ""
		
		for index,i in enumerate(T):
			temp_k += "{},{}\n".format(1/i,parameter[index])
			temp_k_max += "{},{}\n".format(1/i,pExp_max[index])
			temp_k_min += "{},{}\n".format(1/i,pExp_min[index])
			temp_kappa_max += "{},{}\n".format(1/i,parameter_max[index])
			temp_kappa_min +="{},{}\n".format(1/i,parameter_min[index])
			temp_p_opt += "{},{}\n".format(1/i,parameter_opt[index])
			temp_p_opt_max += "{},{}\n".format(1/i,parameter_opt_max[index])
			temp_p_opt_min += "{},{}\n".format(1/i,parameter_opt_min[index])
			temp_Tdeptunsrt += "{},{}\n".format(1/i,TDeptUnsrt[index])	
			
		plot_k.write(temp_k) 
		plot_k_max.write(temp_k_max) 
		plot_k_min.write(temp_k_min ) 
		plot_kappa_max.write(temp_kappa_max) 
		plot_kappa_min.write(temp_kappa_min) 
		plot_p_opt.write(temp_p_opt)
		plot_p_opt_max.write(temp_p_opt_max) 
		plot_p_opt_min.write(temp_p_opt_min) 
		printTDeptUnsrt.write(temp_Tdeptunsrt) 
		
		plot_k.close()
		plot_k_max.close()
		plot_k_min.close()
		plot_kappa_max.close()
		plot_kappa_min.close()
		plot_p_opt.close()
		plot_p_opt_max.close() 
		plot_p_opt_min.close() 
		printTDeptUnsrt.close()
		
		
	####f for third Body collision ####
	
	for index,m in enumerate(thirdBody_index):
		m_ = th.split("-")[0]
		m_species = th.split("-")[1] 		
		nominal = float(Input_file_reader.MechParsing(mech_file_location).getThirdBodyCollisionEff(m,m_species))	
		
		choleskyMatrix = unsrt_data[m].getCholeskyMat()
		zeta = unsrt_data[m].getZeta()
		T = np.asarray(unsrt_data[m].getTemperature())
		TDeptUnsrt = unsrt_data[m].getUncertainty()
		
		beta_opt = np.asarray(Beta_Dict[m])
				
		
		pMin = nominal  - np.asarray(np.dot(choleskyMatrix,zeta)).flatten()
		pMax = nominal  + np.asarray(np.dot(choleskyMatrix,zeta)).flatten()
			
		parameter = nominal*(T/T)
		pExp_max = TDeptUnsrt+parameter
		pExp_min = parameter-TDeptUnsrt
		
		
		zeta_star = beta_opt*zeta
		zeta_pmax = zeta_star+np.dot(Posterior_cov[3][index],np.asarray([1]))
		zeta_pmin = zeta_star-np.dot(Posterior_cov[3][index],np.asarray([1]))
		p_opt =  nominal + np.asarray(np.dot(choleskyMatrix,zeta_star)).flatten()
		parameter_opt = p_opt*(T/T)
		p_opt_max =  nominal + np.asarray(np.dot(choleskyMatrix,zeta_pmax)).flatten()
		p_opt_min = nominal  + np.asarray(np.dot(choleskyMatrix,zeta_pmin)).flatten()
		
		parameter_max = pMax*(T/T)
		parameter_min = pMin*(T/T)
		parameter_opt_max = p_opt_max*(T/T)
		parameter_opt_min = p_opt_min*(T/T)	
		fig = plt.figure()
		
		plt.plot((1/T),parameter,label='Nominal {}-{}'.format(m_,m_species))
		plt.plot((1/T),pExp_max,label='Upper exp unsrt')
		plt.plot((1/T),pExp_min,label='Lower exp unsrt')
		plt.plot((1/T),parameter_max,'--',label='Prior_Max')
		plt.plot((1/T),parameter_min,'--',label ='Prior_Min')
		plt.legend()
		plt.savefig(os.getcwd()+'/Plots/Kappa_vs_T/'+m_+":"+m_species+'.png')
		plt.close()
		
		fig = plt.figure()
		plt.plot((1/T),parameter,"-b",label='Nominal rate constant (k)')
		plt.plot((1/T),parameter_max,"--b",label='Upper prior TDeptUnsrt limit (k+f)')
		plt.plot((1/T),parameter_min,"--b",label='Lower prior TDeptUnsrt limit (k-f)')
		plt.plot((1/T),parameter_opt,'-r',label='opt {}-{}'.format(m_,m_species))
		plt.plot((1/T),parameter_opt_max,'--r',label='Upper posterior unsrt limit')
		plt.plot((1/T),parameter_opt_min,'--r',label ='Lower posterior unsrt limit')
		plt.legend()
		plt.savefig(os.getcwd()+'/Plots/Posterior_unsrt/'+m_+":"+m_species+'.png')
		plt.close()
		
		
		plot_k = open("./Plots/PlotData/Nominal_"+m_+":"+m_species+".txt","+a")
		plot_k_max = open("./Plots/PlotData/Prior_Exp_Max"+m_+":"+m_species+".txt","+a")
		plot_k_min = open("./Plots/PlotData/Prior_Exp_Min"+m_+":"+m_species+".txt","+a")
		plot_kappa_max = open("./Plots/PlotData/Prior_Max_zeta"+m_+":"+m_species+".txt","+a")
		plot_kappa_min = open("./Plots/PlotData/Prior_Min_zeta_"+m_+":"+m_species+".txt","+a")
		plot_p_opt = open("./Plots/PlotData/Optimal_"+m_+":"+m_species+".txt","+a")
		plot_p_opt_max = open("./Plots/PlotData/posterior_Max_zeta"+m_+":"+m_species+".txt","+a")
		plot_p_opt_min = open("./Plots/PlotData/posterior_Min_zeta_"+m_+":"+m_species+".txt","+a")
		printTDeptUnsrt = open("./Plots/PlotData/TDeptUnsrt_"+m_+":"+m_species+".txt","+a")
		
		temp_k = ""
		temp_k_max = ""
		temp_k_min = ""
		temp_kappa_max = ""
		temp_kappa_min = ""
		temp_p_opt = ""
		temp_p_opt_max = ""
		temp_p_opt_min = ""
		temp_Tdeptunsrt = ""
		
		for index,i in enumerate(T):
			temp_k += "{},{}\n".format(1/i,parameter[index])
			temp_k_max += "{},{}\n".format(1/i,pExp_max[index])
			temp_k_min += "{},{}\n".format(1/i,pExp_min[index])
			temp_kappa_max += "{},{}\n".format(1/i,parameter_max[index])
			temp_kappa_min +="{},{}\n".format(1/i,parameter_min[index])
			temp_p_opt += "{},{}\n".format(1/i,parameter_opt[index])
			temp_p_opt_max += "{},{}\n".format(1/i,parameter_opt_max[index])
			temp_p_opt_min += "{},{}\n".format(1/i,parameter_opt_min[index])
			temp_Tdeptunsrt += "{},{}\n".format(1/i,TDeptUnsrt[index])	
			
		plot_k.write(temp_k) 
		plot_k_max.write(temp_k_max) 
		plot_k_min.write(temp_k_min ) 
		plot_kappa_max.write(temp_kappa_max) 
		plot_kappa_min.write(temp_kappa_min) 
		plot_p_opt.write(temp_p_opt)
		plot_p_opt_max.write(temp_p_opt_max) 
		plot_p_opt_min.write(temp_p_opt_min) 
		printTDeptUnsrt.write(temp_Tdeptunsrt) 
		
		plot_k.close()
		plot_k_max.close()
		plot_k_min.close()
		plot_kappa_max.close()
		plot_kappa_min.close()
		plot_p_opt.close()
		plot_p_opt_max.close() 
		plot_p_opt_min.close() 
		printTDeptUnsrt.close()
	
	######## Thermo #########
	
	for index,th in enumerate(thermo_index):
		th_ = th.split("-")[0]
		th_type = th.split("-")[1]
		if th_type == "low":
			continue
		
		elif th_type == "high":
			continue
		
		else:
			print("Invalid Input for ENTHALPY")
		
		
		
		if th_classification == "h":
			if th_type == "low":
				choleskyMatrix = unsrt_data[th].getCholeskyMat()
				zeta = unsrt_data[th].getZeta()
				T = unsrt_data[th].getTemperature()
				TDeptUnsrt = unsrt_data[th].getUncertainty()
				nominal = Input_file_reader.ThermoParsing(thermo_file_location).getEnthalpy_low(th)
				beta_opt = np.asarray(Beta_Dict[th])
				tag = "low"
				pClass = "Enthalpy"	
			elif th_type == "high":
				choleskyMatrix = unsrt_data[th].getCholeskyMat()
				zeta = unsrt_data[th].getZeta()
				T = unsrt_data[th].getTemperature()
				TDeptUnsrt = unsrt_data[th].getUncertainty()
				nominal = Input_file_reader.ThermoParsing(thermo_file_location).getEnthalpy_high(th)
				beta_opt = np.asarray(Beta_Dict[th])
				tag = "high"	
				pClass = "Enthalpy"
			else:
				print("Invalid Input for ENTHALPY")
		
		
		
		elif th_classification == "e":
			if th_type == "low":
				choleskyMatrix = unsrt_data[th].getCholeskyMat()
				zeta = unsrt_data[th].getZeta()
				T = unsrt_data[th].getTemperature()
				TDeptUnsrt = unsrt_data[th].getUncertainty()
				nominal = Input_file_reader.ThermoParsing(thermo_file_location).getEntropy_low(th)
				beta_opt = np.asarray(Beta_Dict[th])
				tag = "low"
				pClass = "Entropy"				
			elif th_type == "high":
				choleskyMatrix = unsrt_data[th].getCholeskyMat()
				zeta = unsrt_data[th].getZeta()
				T = unsrt_data[th].getTemperature()
				TDeptUnsrt = unsrt_data[th].getUncertainty()
				nominal = Input_file_reader.ThermoParsing(thermo_file_location).getEntropy_high(th)
				beta_opt = np.asarray(Beta_Dict[th])
				tag = "high"
				pClass = "Entropy"	
			else:
				print("Invalid Input for ENTROPY")		
		else:		
			print("Invalid Input for thermo")
			
		
		pMin = nominal  - np.asarray(np.dot(choleskyMatrix,zeta)).flatten()
		pMax = nominal  + np.asarray(np.dot(choleskyMatrix,zeta)).flatten()
			
		parameter = nominal*(T/T)
		pExp_max = TDeptUnsrt+parameter
		pExp_min = parameter-TDeptUnsrt
		
		
		zeta_star = beta_opt*zeta
		zeta_pmax = zeta_star+np.dot(Posterior_cov[4][index],np.asarray([1]))
		zeta_pmin = zeta_star-np.dot(Posterior_cov[4][index],np.asarray([1]))
		p_opt =  nominal + np.asarray(np.dot(choleskyMatrix,zeta_star)).flatten()
		parameter_opt = p_opt*(T/T)
		p_opt_max =  nominal + np.asarray(np.dot(choleskyMatrix,zeta_pmax)).flatten()
		p_opt_min = nominal  + np.asarray(np.dot(choleskyMatrix,zeta_pmin)).flatten()
		
		parameter_max = pMax*(T/T)
		parameter_min = pMin*(T/T)
		parameter_opt_max = p_opt_max*(T/T)
		parameter_opt_min = p_opt_min*(T/T)	
		fig = plt.figure()
		
		plt.plot((1/T),parameter,label='Nominal {}-{}'.format(pClass,tag))
		plt.plot((1/T),pExp_max,label='Upper exp unsrt')
		plt.plot((1/T),pExp_min,label='Lower exp unsrt')
		plt.plot((1/T),parameter_max,'--',label='Prior_Max')
		plt.plot((1/T),parameter_min,'--',label ='Prior_Min')
		plt.legend()
		plt.savefig(os.getcwd()+'/Plots/Kappa_vs_T/'+th_+":"+tag+'.png')
		plt.close()
		
		fig = plt.figure()
		plt.plot((1/T),parameter,"-b",label='Nominal rate constant (k)')
		plt.plot((1/T),parameter_max,"--b",label='Upper prior TDeptUnsrt limit (k+f)')
		plt.plot((1/T),parameter_min,"--b",label='Lower prior TDeptUnsrt limit (k-f)')
		plt.plot((1/T),parameter_opt,'-r',label='opt {}-{}'.format(pClass,tag))
		plt.plot((1/T),parameter_opt_max,'--r',label='Upper posterior unsrt limit')
		plt.plot((1/T),parameter_opt_min,'--r',label ='Lower posterior unsrt limit')
		plt.legend()
		plt.savefig(os.getcwd()+'/Plots/Posterior_unsrt/'+th_+":"+tag+'.png')
		plt.close()
		
		plot_k = open("./Plots/PlotData/Nominal_"+th_+":"+tag+".txt","+a")
		plot_k_max = open("./Plots/PlotData/Prior_Exp_Max"+th_+":"+tag+".txt","+a")
		plot_k_min = open("./Plots/PlotData/Prior_Exp_Min"+th_+":"+tag+".txt","+a")
		plot_kappa_max = open("./Plots/PlotData/Prior_Max_zeta"+th_+":"+tag+".txt","+a")
		plot_kappa_min = open("./Plots/PlotData/Prior_Min_zeta_"+th_+":"+tag+".txt","+a")
		plot_p_opt = open("./Plots/PlotData/Optimal_"+th_+":"+tag+".txt","+a")
		plot_p_opt_max = open("./Plots/PlotData/posterior_Max_zeta"+th_+":"+tag+".txt","+a")
		plot_p_opt_min = open("./Plots/PlotData/posterior_Min_zeta_"+th_+":"+tag+".txt","+a")
		printTDeptUnsrt = open("./Plots/PlotData/TDeptUnsrt_"+th_+":"+tag+".txt","+a")
		
		temp_k = ""
		temp_k_max = ""
		temp_k_min = ""
		temp_kappa_max = ""
		temp_kappa_min = ""
		temp_p_opt = ""
		temp_p_opt_max = ""
		temp_p_opt_min = ""
		temp_Tdeptunsrt = ""
		
		for index,i in enumerate(T):
			temp_k += "{},{}\n".format(1/i,parameter[index])
			temp_k_max += "{},{}\n".format(1/i,pExp_max[index])
			temp_k_min += "{},{}\n".format(1/i,pExp_min[index])
			temp_kappa_max += "{},{}\n".format(1/i,parameter_max[index])
			temp_kappa_min +="{},{}\n".format(1/i,parameter_min[index])
			temp_p_opt += "{},{}\n".format(1/i,parameter_opt[index])
			temp_p_opt_max += "{},{}\n".format(1/i,parameter_opt_max[index])
			temp_p_opt_min += "{},{}\n".format(1/i,parameter_opt_min[index])
			temp_Tdeptunsrt += "{},{}\n".format(1/i,TDeptUnsrt[index])	
			
		plot_k.write(temp_k) 
		plot_k_max.write(temp_k_max) 
		plot_k_min.write(temp_k_min ) 
		plot_kappa_max.write(temp_kappa_max) 
		plot_kappa_min.write(temp_kappa_min) 
		plot_p_opt.write(temp_p_opt)
		plot_p_opt_max.write(temp_p_opt_max) 
		plot_p_opt_min.write(temp_p_opt_min) 
		printTDeptUnsrt.write(temp_Tdeptunsrt) 
		
		plot_k.close()
		plot_k_max.close()
		plot_k_min.close()
		plot_kappa_max.close()
		plot_kappa_min.close()
		plot_p_opt.close()
		plot_p_opt_max.close() 
		plot_p_opt_min.close() 
		printTDeptUnsrt.close()
	

	
		
	for index,ts in enumerate(transport_index):
		ts_ = ts.split("-")[0]
		ts_type = ts.split("-")[1]
		if ts_type == "LJe":
			choleskyMatrix = unsrt_data[ts].getCholeskyMat()
			zeta = unsrt_data[ts].getZeta()
			T = unsrt_data[ts].getTemperature()
			TDeptUnsrt = unsrt_data[ts].getUncertainty()
			nominal =Input_file_reader.TransportParsing(transport_file_location).getTransportData(ts)
			norm = float(nominal[1])
			beta_opt = np.asarray(Beta_Dict[ts])
			tag="Lennard-Jonnes potential well depth"
		elif ts_type == "LJs":
			
			choleskyMatrix = unsrt_data[ts].getCholeskyMat()
			zeta = unsrt_data[ts].getZeta()
			T = unsrt_data[ts].getTemperature()
			TDeptUnsrt = unsrt_data[ts].getUncertainty()
			nominal = Input_file_reader.TransportParsing(transport_file_location).getTransportData(ts)
			norm = float(nominal[2])
			beta_opt = np.asarray(Beta_Dict[ts])
			tag="Hard sphere Dia"	
		######## Third Body collision efficiencies ########
		
		pMin = np.array([norm])  - np.array(np.dot(choleskyMatrix,zeta)).flatten()
		pMax = np.array([norm])  + np.array(np.dot(choleskyMatrix,zeta)).flatten()
			
		parameter = norm*(T/T)
		pExp_max = TDeptUnsrt+parameter
		pExp_min = parameter-TDeptUnsrt
		
		
		zeta_star = beta_opt*zeta
		zeta_pmax = zeta_star+np.dot(Posterior_cov[5][index],np.asarray([1]))
		zeta_pmin = zeta_star-np.dot(Posterior_cov[5][index],np.asarray([1]))
		p_opt =  np.array([norm]) + np.asarray(np.dot(choleskyMatrix,zeta_star)).flatten()
		parameter_opt = p_opt*(T/T)
		p_opt_max =  np.array([norm]) + np.asarray(np.dot(choleskyMatrix,zeta_pmax)).flatten()
		p_opt_min = np.array([norm])  + np.asarray(np.dot(choleskyMatrix,zeta_pmin)).flatten()
		parameter_max = pMax*(T/T)
		parameter_min = pMin*(T/T)
		parameter_opt_max = p_opt_max*(T/T)
		parameter_opt_min = p_opt_min*(T/T)	
		fig = plt.figure()
		
		plt.plot((1/T),parameter,label='Nominal {}-{}'.format(pClass,tag))
		plt.plot((1/T),pExp_max,label='Upper exp unsrt')
		plt.plot((1/T),pExp_min,label='Lower exp unsrt')
		plt.plot((1/T),parameter_max,'--',label='Prior_Max')
		plt.plot((1/T),parameter_min,'--',label ='Prior_Min')
		plt.legend()
		plt.savefig(os.getcwd()+'/Plots/Kappa_vs_T/'+ts_+":"+tag+'.png')
		plt.close()
		
		fig = plt.figure()
		plt.plot((1/T),parameter,"-b",label='Nominal rate constant (k)')
		plt.plot((1/T),parameter_max,"--b",label='Upper prior TDeptUnsrt limit (k+f)')
		plt.plot((1/T),parameter_min,"--b",label='Lower prior TDeptUnsrt limit (k-f)')
		plt.plot((1/T),parameter_opt,'-r',label='opt {}-{}'.format(pClass,tag))
		plt.plot((1/T),parameter_opt_max,'--r',label='Upper posterior unsrt limit')
		plt.plot((1/T),parameter_opt_min,'--r',label ='Lower posterior unsrt limit')
		plt.legend()
		plt.savefig(os.getcwd()+'/Plots/Posterior_unsrt/'+ts_+":"+tag+'.png')
		plt.close()
		
		plot_k = open("./Plots/PlotData/Nominal_"+ts_+":"+tag+".txt","+a")
		plot_k_max = open("./Plots/PlotData/Prior_Exp_Max"+ts_+":"+tag+".txt","+a")
		plot_k_min = open("./Plots/PlotData/Prior_Exp_Min"+ts_+":"+tag+".txt","+a")
		plot_kappa_max = open("./Plots/PlotData/Prior_Max_zeta"+ts_+":"+tag+".txt","+a")
		plot_kappa_min = open("./Plots/PlotData/Prior_Min_zeta_"+ts_+":"+tag+".txt","+a")
		plot_p_opt = open("./Plots/PlotData/Optimal_"+ts_+":"+tag+".txt","+a")
		plot_p_opt_max = open("./Plots/PlotData/posterior_Max_zeta"+ts_+":"+tag+".txt","+a")
		plot_p_opt_min = open("./Plots/PlotData/posterior_Min_zeta_"+ts_+":"+tag+".txt","+a")
		printTDeptUnsrt = open("./Plots/PlotData/TDeptUnsrt_"+ts_+":"+tag+".txt","+a")
		
		temp_k = ""
		temp_k_max = ""
		temp_k_min = ""
		temp_kappa_max = ""
		temp_kappa_min = ""
		temp_p_opt = ""
		temp_p_opt_max = ""
		temp_p_opt_min = ""
		temp_Tdeptunsrt = ""
		
		for index,i in enumerate(T):
			temp_k += "{},{}\n".format(1/i,parameter[index])
			temp_k_max += "{},{}\n".format(1/i,pExp_max[index])
			temp_k_min += "{},{}\n".format(1/i,pExp_min[index])
			temp_kappa_max += "{},{}\n".format(1/i,parameter_max[index])
			temp_kappa_min +="{},{}\n".format(1/i,parameter_min[index])
			temp_p_opt += "{},{}\n".format(1/i,parameter_opt[index])
			temp_p_opt_max += "{},{}\n".format(1/i,parameter_opt_max[index])
			temp_p_opt_min += "{},{}\n".format(1/i,parameter_opt_min[index])
			temp_Tdeptunsrt += "{},{}\n".format(1/i,TDeptUnsrt[index])	
			
		plot_k.write(temp_k) 
		plot_k_max.write(temp_k_max) 
		plot_k_min.write(temp_k_min ) 
		plot_kappa_max.write(temp_kappa_max) 
		plot_kappa_min.write(temp_kappa_min) 
		plot_p_opt.write(temp_p_opt)
		plot_p_opt_max.write(temp_p_opt_max) 
		plot_p_opt_min.write(temp_p_opt_min) 
		printTDeptUnsrt.write(temp_Tdeptunsrt) 
		
		plot_k.close()
		plot_k_max.close()
		plot_k_min.close()
		plot_kappa_max.close()
		plot_kappa_min.close()
		plot_p_opt.close()
		plot_p_opt_max.close() 
		plot_p_opt_min.close() 
		printTDeptUnsrt.close()
		
def plot_graphs(target_list,init_guess,Posterior_cov,opt):
	case_dir = range(0,len(target_list))
	DataSet_Tig = []
	DataSet_Fsl = []
	DataSet_Scp = []

	for case in case_dir:
		if target_list[case].target == "Tig":
			DataSet_Tig.append(target_list[case].d_set)
		if target_list[case].target == "Fsl":
			DataSet_Fsl.append(target_list[case].d_set)
		if target_list[case].target == "Scp":
			DataSet_Scp.append(target_list[case].d_set)

	DataSet_Tig = list(OrderedDict.fromkeys(DataSet_Tig))
	DataSet_Fsl = list(OrderedDict.fromkeys(DataSet_Fsl))
	DataSet_Scp = list(OrderedDict.fromkeys(DataSet_Scp))


	for ind in DataSet_Tig:
		temp_T = []
		temp_Tig = []
		temp_P = []
		temp_sigma_exp = []
		temp_error = []
		temp_phi = []
		temp_opt_sim = []
		temp_init_sim = []
		prior_unsrt_res = []
		posterior_unsrt_res = []
		for case in target_list:
			if case.target == "Tig":
				if case.d_set == ind:
					temp_T.append(case.temperature)
					temp_Tig.append(case.observed)
					temp_sigma_exp.append(case.std_dvtn)
					temp_opt_sim.append(case.calculated_target_value(opt))
					temp_init_sim.append(case.calculated_target_value(init_guess))
					temp_P.append(case.pressure)
					temp_phi.append(case.phi)
					#1/9 if specified
					prior_cov = (1)*np.identity(len(opt))
					prior_unsrt_res.append(case.model_response_uncertainty(init_guess,prior_cov,2))
					posterior_unsrt_res.append(case.model_response_uncertainty(opt,Posterior_cov,2)) #Error in micro seconds
		#print(prior_unsrt_res)
		#print(posterior_unsrt_res)
		fig = plt.figure()
		plt.xlabel('Temperature (1000/K)')
		plt.ylabel('Ignition Delay ($ tau$) [$\mu$ s]')
		plt.yscale('log')
		plt.title('P =  {} and Phi = {}'.format(temp_P[0],temp_phi[0]), fontsize = 10)
		plt.plot(np.asarray(temp_T),np.asarray(temp_Tig),'o',label = "Observed")
		plt.errorbar(np.asarray(temp_T),np.asarray(temp_Tig),yerr = np.asarray(temp_sigma_exp),fmt='o',label = "exp std_dvtn")
		plt.plot(np.asarray(temp_T),(np.asarray(temp_init_sim)),'--',label = "Li et al. (2007)")
		plt.plot(np.asarray(temp_T),(np.asarray(temp_opt_sim)),'-',label = "Optimized mechanism")
		plt.legend()
		plt.savefig('Plots/Tig/DataSet_Tig_'+str(ind)+'.png')
		plt.close()
		fig = plt.figure()
		plt.xlabel('Temperature (1000/K)')
		plt.ylabel('Ignition Delay ($ tau$) [$\mu$ s]')
		plt.yscale("log")
		plt.title('Uncertainty Plot [P =  {} and Phi = {}]'.format(temp_P[0],temp_phi[0]), fontsize = 10)
		plt.plot(np.asarray(temp_T),np.asarray(temp_Tig),'o',label = "Observerd")
		plt.errorbar(np.asarray(temp_T),np.asarray(temp_Tig),yerr = np.asarray(temp_sigma_exp),fmt='o',label = "exp std_dvtn")
		plt.fill_between(np.asarray(temp_T), (np.asarray(temp_init_sim)-np.asarray(prior_unsrt_res)),(np.asarray(temp_init_sim)+np.asarray(prior_unsrt_res)),facecolor='r',alpha=0.5)
		plt.plot(np.asarray(temp_T),np.asarray(temp_init_sim),'-r',label = "Li et al. (2007)") #Ignition delay in micro seconds
		plt.fill_between(np.asarray(temp_T), np.asarray(temp_opt_sim)-np.asarray(posterior_unsrt_res), np.asarray(temp_opt_sim)+np.asarray(posterior_unsrt_res),facecolor='b',alpha=0.5)
		plt.plot(np.asarray(temp_T),(np.asarray(temp_opt_sim)),'-b',label = "Optimized mechanism") #Ignition delay in micro seconds
		plt.legend()
		plt.savefig('Plots/Tig/Unsrt_Tig_'+str(ind)+'.png')
		plt.close()
	for ind in DataSet_Fsl:
		temp_T = []
		temp_fsl = []
		temp_P = []
		temp_sigma_exp = []
		temp_error = []
		temp_phi = []
		temp_opt_sim = []
		temp_init_sim = []
		prior_unsrt_res = []
		posterior_unsrt_res = []
		for case in target_list:
			if case.target == "Fsl":
				if case.d_set == ind:
					temp_T.append(case.temperature)
					temp_fsl.append(case.observed)
					temp_opt_sim.append(case.calculated_target_value(opt))
					temp_sigma_exp.append(case.std_dvtn)
					temp_init_sim.append(case.calculated_target_value(init_guess))
					temp_P.append(case.pressure)
					temp_phi.append(case.phi)
					prior_cov = np.identity(len(opt))
					prior_unsrt_res.append(case.model_response_uncertainty(init_guess,prior_cov,2))
					posterior_unsrt_res.append(case.model_response_uncertainty(opt,Posterior_cov,2)) #Error in micro seconds
		fig = plt.figure()
		plt.xlabel('Equivalence Ratio')
		plt.ylabel('Laminar Burning Velocity ($\S_u /cm s^{-1}$)')
		plt.title('T =  {} and P = {}'.format(temp_T[0],temp_P[0]), fontsize = 10)
		plt.plot(np.asarray(temp_phi),np.asarray(temp_fsl),'.',label = "Observerd")
		plt.errorbar(np.asarray(temp_T),np.asarray(temp_Tig),yerr = np.asarray(temp_sigma_exp),label = "exp std_dvtn")
		plt.plot(np.asarray(temp_phi),np.exp(np.asarray(temp_init_sim)),'--',label = "Li et al. (2007)")
		plt.plot(np.asarray(temp_phi),np.exp(np.asarray(temp_opt_sim)),'-',label = "optimized mechanism")
		plt.legend()
		plt.savefig('Plots/Fsl/DataSet_Fsl'+str(ind)+'.png')
		plt.close()
	
		fig = plt.figure()
		plt.xlabel('Equivalence ratio')
		plt.ylabel('Laminar Flame Speed $cm s^{-1}$')
		plt.title('Uncertainty Plot [P =  {} and Phi = {}]'.format(temp_P[0],temp_phi[0]), fontsize = 10)
		plt.plot(np.asarray(temp_phi),np.asarray(temp_fsl),'.',label = "Observerd")
		plt.errorbar(np.asarray(temp_phi),np.asarray(temp_Tig),yerr = np.asarray(temp_sigma_exp),label = "exp std_dvtn")
		plt.fill_between(np.asarray(temp_phi), (np.asarray(temp_init_sim)-np.asarray(prior_unsrt_res)),(np.asarray(temp_init_sim)+np.asarray(prior_unsrt_res)),facecolor='r',alpha=0.5)
		plt.plot(np.asarray(temp_phi),np.asarray(temp_init_sim),'-r',label = "optimized mechanism") #Ignition delay in micro seconds
		plt.fill_between(np.asarray(temp_phi), np.asarray(temp_opt_sim)-np.asarray(posterior_unsrt_res), np.asarray(temp_opt_sim)+np.asarray(posterior_unsrt_res),facecolor='b',alpha=0.5)
		plt.plot(np.asarray(temp_phi),(np.asarray(temp_opt_sim)),'-b',label = "optimized mechanism") #Ignition delay in micro seconds
		plt.legend()
		plt.savefig('Plots/Fsl/Unsrt_Fsl_'+str(ind)+'.png')
		plt.close()
	
	for ind in DataSet_Scp:
		temp_T = []
		temp_Tig = []
		temp_P = []
		temp_time = []
		temp_sigma_exp = []
		temp_error = []
		temp_species = []
		temp_phi = []
		temp_opt_sim = []
		temp_init_sim = []
		prior_unsrt_res = []
		posterior_unsrt_res = []
		for case in target_list:
			if case.target == "Scp":
				if case.d_set == ind:
					temp_T.append(case.temperature)
					temp_Tig.append(case.observed)
					temp_sigma_exp.append(case.std_dvtn)
					temp_opt_sim.append(case.calculated_target_value(opt))
					temp_init_sim.append(case.calculated_target_value(init_guess))
					temp_species.append(case.species)
					temp_P.append(case.pressure)
					temp_phi.append(case.phi)
					temp_time.append(case.Tend)
					prior_cov = np.identity(len(opt))
					
					prior_unsrt_res.append(case.model_response_uncertainty(init_guess,prior_cov,2))
					posterior_unsrt_res.append(case.model_response_uncertainty(opt,Posterior_cov,2)) #Error in micro seconds
		#print(temp_unsrt_res)			
		fig = plt.figure()
		plt.xlabel('Time [ms]')
		plt.ylabel('Species mole fraction (%)')
		plt.title('Species = {}, P =  {} and Phi = {}'.format(temp_species[0],temp_P[0],temp_phi[0]), fontsize = 10)
		plt.plot(np.asarray(temp_phi),np.asarray(temp_Tig),'.',label = "Observerd")
		plt.errorbar(np.asarray(temp_T),np.asarray(temp_Tig),yerr = np.asarray(temp_sigma_exp),label = "exp std_dvtn")
		plt.plot(np.asarray(temp_time),np.asarray(temp_init_sim),'--',label = "Li et al. (2007)")
		plt.plot(np.asarray(temp_time),np.asarray(temp_opt_sim),'-',label = "optimized mechanism")
		plt.legend()
		plt.savefig('Plots/Scp/DataSet_Scp'+str(ind)+'.png')
		plt.close()
		fig = plt.figure()
		plt.xlabel('Time [ms]')
		plt.ylabel('Species mole fraction (%)')
		#plt.xlim(0,temp_time[-1])
		#plt.ylim(0,1.1)
		plt.title('Uncertainty Plot [P =  {} and Phi = {}]'.format(temp_P[0],temp_phi[0]), fontsize = 10)
		plt.plot(np.asarray(temp_phi),np.asarray(temp_Tig),'.',label = "Observerd")
		plt.errorbar(np.asarray(temp_T),np.asarray(temp_Tig),yerr = np.asarray(temp_sigma_exp),label = "exp std_dvtn")
		plt.fill_between(np.asarray(temp_T), (np.asarray(temp_init_sim)-np.asarray(prior_unsrt_res)),(np.asarray(temp_init_sim)+np.asarray(prior_unsrt_res)),facecolor='r',alpha=0.5)
		plt.plot(np.asarray(temp_T),np.asarray(temp_init_sim),'-r',label = "optimized mechanism") #Ignition delay in micro seconds
		plt.fill_between(np.asarray(temp_phi), np.asarray(temp_opt_sim)-np.asarray(posterior_unsrt_res), np.asarray(temp_opt_sim)+np.asarray(posterior_unsrt_res),facecolor='b',alpha=0.5)
		plt.plot(np.asarray(temp_T),(np.asarray(temp_opt_sim)),'-b',label = "optimized mechanism") #Ignition delay in micro seconds
		plt.legend()
		plt.savefig('Plots/Scp/Unsrt_Scp_'+str(ind)+'.png')
		plt.close()
def plot_BranchingRatio(unoptimized,optimized,rxnList,unsrt):
	ratioDic = MechanismManipulator.MechanismManipulator(unoptimized,rxnList,unsrt).BranchingDict
	#print(ratioDic)
'''
		
