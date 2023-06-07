import numpy as np
import scipy as sp
from scipy.optimize import minimize
import os, sys, re, threading, subprocess, time
from sklearn.preprocessing import PolynomialFeatures
import matplotlib.pyplot as plt
#program specific modules
import uncertainty


### KEY WORDS #######
mechanism = "mechanism"
targets_count = "targets_count"
targets = "targets:"
home_dir = os.getcwd()
fuel_name = "fuel"
global_reaction_eqn = "global_reaction"
parallel_thread_count = "parallel_threads"
plus = "plus"
minus = "minus"
unsrt = "uncertainty_data"
thermo_file = "thermo_file"
trans_file = "trans_file"

if len(sys.argv) > 1:
	input_file = open(sys.argv[1],'r')
	lines = input_file.readlines()
	print("Input file found\n")
else:
	print("Please enter a valid input file name as arguement. \n For details of preparing the input file, please see the UserManual\n\nProgram exiting")
	exit()
#!!!!!!! GET MECHANISM FILE , number of targets  from the input file !!!!!!!!!
for line in lines:
	######check for comments!!!!!!!!!!
	if "#" in line:
		line = line[:line.index('#')]
		
	word = line.split()
	if mechanism in word:
		mech_file_location = os.path.abspath(word[word.index(mechanism) + 2])
			
	if thermo_file in word:
		thermo_file_location = os.path.abspath(word[word.index(thermo_file) + 2])
	
	if trans_file in word:
		trans_file_location = os.path.abspath(word[word.index(trans_file) + 2])
		
	if targets_count in word:
		no_of_targets = int(word[word.index(targets_count) + 2])
		
	if fuel_name in word:
		fuel = word[word.index(fuel_name) + 2]
		
	if global_reaction_eqn in word:
		global_reaction = word[word.index(global_reaction_eqn) + 2:]
		
	if parallel_thread_count in word:
		parallel_threads = int(word[word.index(parallel_thread_count) + 2])+1
		
	if unsrt in word:
		unsrt_location = os.path.abspath(word[word.index(unsrt) + 2])
		
	if targets in word:
		start_line = lines.index(line) + 1
		stop_line = start_line + no_of_targets

UncertDataSet = uncertainty.uncertaintyData(unsrt_location);
unsrt_data, reaction_index = UncertDataSet.extract_uncertainty();

f = open(mech_file_location).readlines()
mechLine = ""
PerturbingReactions = {}
for r in reaction_index:
	for i in f:
		if r+":" == i.split(" ")[0]: #logic to extract reaction index from the string
			if  "}" in i:
				mechLine += i
				print(mechLine)
				PerturbingReactions[i.split(" ")[0]] = mechLine
				mechLine = ""
			else:
				j = f.index(i)
				while "}" not in mechLine:
					mechLine += f[j]
					j = j+1
				print(mechLine)					
				PerturbingReactions[i.split(" ")[0]] = mechLine
				mechLine = ""       	

def Get_cholesky(index,uncert):
	choleskyMatrix = uncert[index].getCholeskyMat()
	return choleskyMatrix
def Objective (z):
	f = (TDepUnsrt*np.log(10)-(L11*z[0]+(L21*z[0]+L22*z[1])*np.log(T)-(L31*z[0]+L32*z[1]+L33*z[2])*(1/T)))
	#print(np.dot(np.transpose(theta),np.dot(L1,eta)))
	obj = np.dot(f,f)
	return obj

#self.uncert.getMatrix(reaction))#cholesky matrix should be numpy array
def Get_TDeptUnsrt(index,uncert):
	T1, T2 = uncert[index].getTempRange()
	N = 1000;	
	T = sp.linspace(T1,T2,N);
	zeta = uncert[index].getZeta()
	T_test = sp.linspace(500,2500,30);
	print("The upper and lower temp are:\n{}\t{}\n\n".format(T1,T2))
	f = uncert[index].getUncertainty(T)		
	return f,T,T_test,zeta

ff = open(os.getcwd()+'/zeta_bounds.txt','a+')
data = ''


print("Reaction index found");
print("################ Reaction list #################")


for r_ in reaction_index:
	k = r_+':';
	C = Get_cholesky(r_,unsrt_data)	#self.uncert.getMatrix(reaction))#cholesky matrix should be numpy array
	print("\nReaction is {}".format(PerturbingReactions[k]))
	M = 3/(np.log(10.0));		
	w = re.compile(r'.*:?.*->.*A\s*?=\s*?(\S*E\S*)?\s*n\s*?=\s*(\S*)?\s*E\s*?=\s*(\S*)?\s*.*\}}?'.format(), re.DOTALL | re.IGNORECASE)
	match1 = re.search(w,PerturbingReactions[k] )
	aString = match1.group(1)
	a = float(aString)
	alpha = np.log(a)
	match1 = re.search(w, PerturbingReactions[k])
	nString = match1.group(2)
	n = float(nString)
	match1 = re.search(w,PerturbingReactions[k])
	eString = match1.group(3)
	e = float(eString)
	epsilon = e/8.314E-3				

	TDepUnsrt,T,T_test,optZeta = Get_TDeptUnsrt(r_,unsrt_data)	

	#optZeta = uncert[index].getZeta()
	
	logArrheniusMin = np.array([alpha,n,epsilon]  - np.dot(C,optZeta)).flatten()
	logArrheniusMax = np.array([alpha,n,epsilon]  + np.dot(C,optZeta)).flatten()
	ArrheniusMax = [logArrheniusMax[0], logArrheniusMax[1], logArrheniusMax[2]]
	ArrheniusMin = [logArrheniusMin[0], logArrheniusMin[1], logArrheniusMin[2]]


	k = alpha + n*np.log(T)-epsilon*(1/T)
	k_max = TDepUnsrt*np.log(10)+k
	k_min = k-TDepUnsrt*np.log(10)
	kappa_min = ArrheniusMin[0] + ArrheniusMin[1]*np.log(T)-ArrheniusMin[2]*(1/T)
	kappa_max = ArrheniusMax[0] + ArrheniusMax[1]*np.log(T)-ArrheniusMax[2]*(1/T)	
	
	fig = plt.figure()
	#plt.title('{}: Zeta = [{} {} {}]'.format(r_,zeta.x[0],zeta.x[1],zeta.x[2]))
	plt.plot((1/T),k,label='Nominal rate constant (k)')
	plt.plot((1/T),k_max,label='Upper TDeptUnsrt limit (k+f)')
	plt.plot((1/T),k_min,label='Lower TDeptUnsrt limit (k-f)')
	plt.plot((1/T),kappa_max,'--',label='Cholesky_kappa_max')
	plt.plot((1/T),kappa_min,'--',label ='Cholesky_kappa_min')
	plt.legend()
	plt.savefig(os.getcwd()+'/Plots/Kappa_vs_T/'+r_+'.png')
	print("Cholesky_mat for this (rxn {}) is:\n {}\n\n".format(r_,C));
	#print("Nominal (rxn {}) = {}\n".format(r_,Arrhenius_o));
	print("Max (rxn {}) = {}\n".format(r_,ArrheniusMax));
	print("Mix (rxn {}) = {}\n\n".format(r_,ArrheniusMin));
	print("######################")
	print("The k-value to select the uncertainity is as follows:")
	for i in range(0,len(T_test)):
		a = np.exp(alpha);
		print("{}\t{:4E}\n".format(T_test[i],a*(T_test[i])**(n)*np.exp(-epsilon/T_test[i])))
	print("######################")	
	#print("Nominal, Min and Max rate of reactions for rxn {} are as follows:\n{}\t{}\t{}\n".format(r_,kappa,kappa_min,kappa_max));
print("###########################################################################")
#for i in range(len(unsrt_data)):
#print(unsrt_data[[:-1]].get_choleskyDecomposedMatrix());
#for r in reaction_index:	
#	PltArrCurve.(mech_file_location,r,reaction_index);
