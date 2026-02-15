import os
import Uncertainty # get the uncertainty extractor
import numpy as np
import math
import multiprocessing
import multiprocessing as mp
from multiprocessing import shared_memory
import subprocess
import time
import sys
import copy
import matplotlib.pyplot as plt
from tqdm import tqdm
import pickle
import tempfile
sys.path.append('/shuffle.so')
import shuffle
from scipy.stats.qmc import Sobol

class DesignMatrix(object):
	def __init__(self,UnsrtData,design,sample_length,ind=None):
		self.unsrt = UnsrtData
		self.sim = sample_length
		self.design = design
		#self.n = ind
		self.n = ind
		allowed_count = int(0.70*multiprocessing.cpu_count())
		self.allowed_count = allowed_count
		
	def getClassA_Curves(self,n_a):
	    """
		    This defination generates n_a numbers of class-A type curves 
	    """
	    sobol_sampler = Sobol(d=1, scramble=True)
	    ClassACurve_dict = {}
	    Generator_dict = {}
	    for rxn in self.unsrt:
		    generator = sobol_sampler.random(n=n_a) #n_a: number of A-type curves required 
		    generator = 2 * generator - 1
		    data = self.unsrt[rxn].data
		    zeta_max = self.unsrt[rxn].zeta.x
		    ClassA_curves = []
		    for gen in generator:
			    ClassA_curves.append(gen[0]*zeta_max)#getting new A-type curves
		    ClassACurve_dict[rxn] = ClassA_curves
		    Generator_dict[rxn] = generator	
	    return ClassACurve_dict,Generator_dict
	    
	def getSamples(self):
		    print("\nStarting to generate design matrix!!\n")
		    if self.design == 'A-facto': 
			    #self.sim = 1 + 2* self.n+self.n*(self.n-1)/2
			    sobol_sampler = Sobol(d=self.n, scramble=True)
			    sobol_samples = sobol_sampler.random(n=self.sim)
			    design_matrix = list(2 * sobol_samples - 1)
			    design_matrix.extend(list(np.eye(self.n)))
			    s =""
			    for row in design_matrix:
				    for element in row:
					    s+=f"{element},"
				    s+="\n"
			    ff = open('DesignMatrix.csv','w').write(s)
			    return np.asarray(design_matrix)
			    
			    if "a_type_samples.pkl" not in os.listdir():
				    a_curves_dict, generator_a = self.getClassA_Curves(n_a)# Returns 100 class-A Arrhenius samples
				    with open('a_type_samples.pkl', 'wb') as file_:
					    pickle.dump(a_curves_dict, file_)
				    print("\nA-type curves generated")
			    else:
				    # Load the object from the file
				    with open('a_type_samples.pkl', 'rb') as file_:
					    a_curves_dict = pickle.load(file_)
				    print("\nA-type curves generated")
				    
			    V_ = {}
			    for rxn in tqdm(self.unsrt,desc="Populating V_"):
				    temp = []
				    for sample_a in a_curves_dict[rxn]: 
					    temp.append(sample_a)
				    for sample_b in b_curves_dict[rxn]: 
					    temp.append(sample_b)
				    for sample_c in c_curves_dict[rxn]: 
					    temp.append(sample_c)
				    V_[rxn] = np.asarray(temp)
			    
			    V_s = {}
			    #Doing random shuffling
			    #Deepcopy the unshuffled samples first
			    #print(n_a,n_b,n_c)
			    #print(len(n_a),len(n_b),len(n_c))
			    V_copy = copy.deepcopy(V_)
			    for rxn in tqdm(self.unsrt, desc="Doing random shuffling"):
				    column = []
				    # Number of shuffles
				    num_shuffles = int(self.sim * 0.8)
				    # Get the values to shuffle
				    values = V_copy[rxn]

				    # Repeat the shuffle process
				    for i in range(4):
					    np.random.shuffle(values)
					    column.append(values.copy())

				    # Flatten the list of arrays
				    column = np.concatenate(column)
				    V_s[rxn] = column
									    
			    delta_n = {}
			    p = {}
			    V_opt = {}
			    V = {}    #Populating V for unshuffled portion of the design matrix
			    ch = {}
			    nominal = {}
			    p_max = {}
			    p_min = {}
			    theta = {}
			    Temp ={}
			    zeta = {}
			    for rxn in tqdm(self.unsrt,desc="Populating unshuffled portion of DM"):
				    ch[rxn] = self.unsrt[rxn].cholskyDeCorrelateMat
				    nominal[rxn] = self.unsrt[rxn].nominal
				    p_max[rxn] = self.unsrt[rxn].P_max
				    p_min[rxn] = self.unsrt[rxn].P_min
				    theta[rxn] = self.unsrt[rxn].Theta
				    Temp[rxn] = self.unsrt[rxn].temperatures
				    zeta[rxn] = self.unsrt[rxn].zeta.x
				    
			    #SAC_3		
			    #plot the zetas
			    #print(percentage)
			    #Solving for a line passing from 3 points
			    d_n = {}
			    _delta_n = {}
			    dict_delta_n = {}
			    percentage = {}
			    for rxn in tqdm(self.unsrt,desc="Generating fSAC samples"):
				    T =np.array([Temp[rxn][0],(Temp[rxn][0] + Temp[rxn][-1])/2,Temp[rxn][-1]])	
				    Theta = np.array([T/T,np.log(T),-1/(T)])
				    P = nominal[rxn]
				    cov = ch[rxn]
				    zet = zeta[rxn]
				    Tp = Temp[rxn]
				    #Tp = np.linspace(300,2500,50)
				    Theta_p = np.array([Tp/Tp,np.log(Tp),-1/(Tp)])
				    P_max = P + np.asarray(np.dot(cov,zet)).flatten();
				    P_min = P - np.asarray(np.dot(cov,zet)).flatten();
				    
				    d_n[rxn] = abs(P[1]-P_max[1])
				    kmax = Theta_p.T.dot(P_max)
				    kmin = Theta_p.T.dot(P_min)
				    ka_o = Theta_p.T.dot(P)
				    M = Theta.T.dot(cov)
				    #fig = plt.figure()
				    #plt.plot(1/Tp,kmax)
				    #plt.plot(1/Tp,kmin)
				    #plt.plot(1/Tp,ka_o)
				    f = abs(kmax-ka_o)
				    temp = []
				    outside = []
				    not_selected = []
				    temp_n = []
				    while len(temp)<1000:
					    random = 2*np.random.rand(3)-1
					    P_right = P + random[0]*np.asarray(np.dot(cov,zet)).flatten()
					    P_mid = P + random[1]*(7/8)*np.asarray(np.dot(cov,zet)).flatten()
					    P_left = P + random[2]*np.asarray(np.dot(cov,zet)).flatten()
					    Theta_left = np.array([T[0]/T[0],np.log(T[0]),-1/(T[0])])
					    Theta_mid = np.array([T[1]/T[1],np.log(T[1]),-1/(T[1])])
					    Theta_right = np.array([T[2]/T[2],np.log(T[2]),-1/(T[2])])
					    kappa_left = Theta_left.T.dot(P_left)
					    kappa_mid = Theta_mid.T.dot(P_mid)
					    kappa_right = Theta_right.T.dot(P_right)
					    kappa_o1 = Theta_left.T.dot(P)
					    kappa_o2= Theta_mid.T.dot(P)
					    kappa_o3 = Theta_right.T.dot(P)
					    kappa = np.array([kappa_left,kappa_mid,kappa_right])
					    kappa_o = np.array([kappa_o1,kappa_o2,kappa_o3])
					    y = kappa - kappa_o
					    zeta_ = np.linalg.solve(M,y)
					    func = np.asarray([(i.T.dot(cov.dot(zeta_))) for i in Theta_p.T])
					    P_found = P + np.asarray(np.dot(cov,zeta_)).flatten()
					    n_ = abs(P[1]-P_found[1])
					    temp_n.append(n_)
					    kappa_found = Theta_p.T.dot(P_found)
					    f_found = abs(kappa_found-ka_o)
					    if max(f)<max(f_found):
						    outside.append(zeta_)
					    if n_ > 2:	
						    #plt.plot(1/Tp,kappa_found,"r-",linewidth=0.25)
						    #temp.append(zeta_)#for unconstrained zeta
						    not_selected.append(zeta_)
						    _delta_n[n_] = zeta_
						    
					    else:
						    temp.append(zeta_)
						    _delta_n[n_] = zeta_
					    #print(zeta_)
				    percentage[rxn] = (len(not_selected)/len(temp))*100
				    delta_n[rxn] = max(temp_n)
				    
				    dict_delta_n[rxn] = _delta_n[max(temp_n)]
				    V[rxn] = temp
					    
			    for i in range(int(self.sim*0.2)):
				    row = []
				    for rxn in self.unsrt:
					    row.extend(V_[rxn][i])
					    
				    design_matrix.append(row)
			    
			    
			    #RANDOM SHUFFLING
			    for i in range(int(self.sim*0.8)):
				    row = []
				    for rxn in self.unsrt:
					    row.extend(V_s[rxn][i])
				    design_matrix.append(row)
				    		    
			    for i in range(1000):
				    row = []
				    for rxn in self.unsrt:
					    row.extend(V[rxn][i])
					    
				    design_matrix.append(row)
			    tok = time.time()
			    print("Time taken to construct Design Matrix: {}".format(tok - tic))
			    s =""
			    for row in design_matrix:
				    for element in row:
					    s+=f"{element},"
				    s+="\n"
			    ff = open('DesignMatrix.csv','w').write(s)	
			    
			    return np.asarray(design_matrix)
