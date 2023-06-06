import xml.etree.ElementTree as ET
import scipy as sp
import numpy as np
from scipy.optimize import minimize
from scipy.optimize import curve_fit
import os,re
import Input_file_reader as IFR
import matplotlib.pyplot as plt
from scipy.linalg import block_diag
from scipy.interpolate import CubicSpline
import multiprocessing
import sys,time
import make_input_file
from scipy.optimize import shgo
from scipy.optimize import BFGS
from MechanismParser import Parser
import subprocess
#############################################################
###	   Uncertainty for arrhenius parameters		 ######
###	   of elementary reactions					  ######
#############################################################
def run(sample,data,length):
	A = UncertaintyExtractor(data)
	a1 =2*np.random.random_sample(1)-1
	a2 = 2*np.random.random_sample(1)-1
	a3 = 2*np.random.random_sample(1)-1
	print(a1,a2,a3)
	#a1 = np.random.uniform(-1,1,1)
	#a2 = np.random.uniform(-1,1,1)
	#a3 = np.random.uniform(-1,1,1)
	#A.generator.append([a1[0],a2[0],a3[0]])
	A.populateValues(a1,a2,a3)
	A.getCovariance(flag=False)
	A.getUnCorrelated(flag = False)
	zeta = A.getB2Zeta(flag=True)
	#print(data)
	#zeta =UncertaintyExtractor(data).getExtremeCurves(sample)	
	del A
	print(zeta)
	
	return (sample,zeta,length)

class workers(object):
	def __init__(self,workers):
		#print("Initialized\n")
		self.pool = multiprocessing.Pool(processes=workers)
		self.progress = []
		self.parallized_zeta = []
	def callback(self,result):
		#print("Entered callback\n")
		self.progress.append(result[0])
		self.parallized_zeta.append(result[1])
		sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(self.progress)/float(result[-1])*100))
		sys.stdout.flush()
		
	def do_job_async(self,data,sampling_points):
		#print("Entered Async\n")
		#print("Entered async\n")
		for args in range(sampling_points):
			x = self.pool.apply_async(run, 
				  args=(1,data,sampling_points), 
				  callback=self.callback)
		self.pool.close()
		self.pool.join()
		print(x.get())
		self.pool.terminate()
		return self.parallized_zeta	

class UncertaintyExtractor(object):
	def __init__(self,data):
		self.data = data
		self.temperatures = data["temperatures"]
		self.uncertainties = data["uncertainties"]
		#cs = CubicSpline(self.temperatures,self.uncertainties)
		#self.temperatures = np.arange(self.temperatures[0],self.temperatures[-1],25)
		#self.uncertainties = cs(self.temperatures)
		self.ArrheniusParams = data["Arrhenius"]
		self.Theta = np.array([self.temperatures/self.temperatures,np.log(self.temperatures),-1/(self.temperatures)])
		self.M = 3.0/np.log(10.0)
		self.guess = np.array([-10.0,-10.0,0.5,200.0,10.0,5.0,1,1])
		self.guess_z = np.array([1,1,1])
		self.guess_z2 = np.array([10,10,100,200])
		self.parallized_zeta = []
		
		self.generator = []
		self.samples = []
		self.kleft_fact = None
		self.kright_fact = None
		self.kmiddle_fact = None
	
	def callback(self,result):
		self.progress.append(result[0])
		self.parallized_zeta.append(result[1])
		sys.stdout.write("\t\t\r{:06.2f}% is complete".format(len(self.progress)/float(result[-1])*100))
		sys.stdout.flush()
	def getPool(self,workers):
		#print("Entered pool")
		self.pool = multiprocessing.Pool(processes=workers)
	def do_job_async(self,sampling_points):
		#print("Entered Async\n")
		
		for args in range(sampling_points):
			self.pool.apply_async(run, 
				  args=(1,self.data,sampling_points), 
				  callback=self.callback)
		self.pool.close()
		self.pool.join()
		self.pool.terminate()	
		
	def getUncertFunc(self,L):
		func = [self.M*np.linalg.norm(np.dot(L.T,i)) for i in self.Theta.T]
		return np.asarray(func)
	
	def getZetaUnsrtKappaFunc(self,L,z):
		func = [(i.T.dot(L.dot(z))) for i in self.theta_for_kappa.T]
		return np.asarray(func)
	
	
	def getZetaUnsrtFunc(self,L,z):
		func = [(i.T.dot(L.dot(z))) for i in self.Theta.T]
		return np.asarray(func)	

	def obj_func(self,guess):
		z = guess
		cov = np.array([[z[0],0,0],[z[1],z[2],0],[z[3],z[4],z[5]]]);#cholesky lower triangular matric		
		f = (self.uncertainties - self.getUncertFunc(cov))/(self.uncertainties/self.M)
		obj = np.dot(f,f)
		return obj

	def const_func(self,guess):
		self.z = guess
		f = np.zeros(len(self.temperatures))
		f = self.uncertainties - self.getUncertFunc()
		return np.amin(f)

	def obj_func_zeta(self,guess):
		M = self.M
		T = self.temperatures
		cov = self.L
		Theta = np.array([T/T,np.log(T),-1/(T)])
		f = (self.uncertainties-self.getZetaUnsrtFunc(cov,guess))
		obj = np.dot(f,f)
		return obj
	
	def obj_func_zeta_b2(self,guess):
		M = self.M
		T = self.temperatures
		cov = self.L
		Theta = np.array([T/T,np.log(T),-1/(T)])
		f = (self.uncertainties-self.getZetaUnsrtFunc(cov,guess[0:-1]))
		obj = np.dot(f,f)
		return obj
		
	def obj_func_zeta_c2(self,guess):
		M = self.M
		T = self.temperatures
		cov = self.L
		Theta = np.array([T/T,np.log(T),-1/(T)])
		f = (self.Yu-self.getZetaUnsrtFunc(cov,guess))
		obj = np.dot(f,f)
		return obj
	
	def const_func_zeta_1(self,z):
		M = self.M
		T = self.temperatures
		cov = self.L
		Theta = np.array([T/T,np.log(T),-1/(T)])
		normTheta = np.sqrt((T/T)**2 + (np.log(T))**2 + (1/T)**2)
		unsrtFunc = M*np.linalg.norm(np.dot(cov.T,Theta))
		uncorrFunc = np.linalg.norm(np.dot(cov,guess))
		QLT = np.asarray(np.dot(Theta.T,np.dot(cov,guess))).flatten()
		f = (self.uncertainties - self.getZetaUnsrtFunc(cov,z))
		#print(np.dot(np.transpose(theta),np.dot(L1,eta)))
		obj = np.dot(f,f)
		return np.amin(f)

	def const_func_zeta_2(self,z):
		M = self.M
		T = self.temperatures[1:-1]
		cov = self.L
		Theta = np.array([T/T,np.log(T),-1/(T)])
		QtLZ = np.asarray([(i.dot(cov.dot(z))) for i in Theta.T])
		f = (self.uncertainties[1:-1]-QtLZ)
		return np.amin(f)
	
	def const_1_typeB_Zeta(self,z):
		M = self.M
		#T = self.temperatures
		P = self.ArrheniusParams
		T = self.temperatures[0]
		cov = self.L
		Pmin = self.P_min
		
		Theta = np.array([T/T,np.log(T),-1/(T)])
		k_min = Theta.dot(Pmin)
		QtLZ = (Theta.T.dot(cov.dot(z)))
		f = Theta.dot(P)-QtLZ
		return k_min - f
	
	def const_2_typeB_Zeta(self,z):
		M = self.M
		#T = self.temperatures
		T = self.const_T
		cov = self.L
		Pmax = self.P_max
		P = self.ArrheniusParams
		Theta = np.array([T/T,np.log(T),-1/(T)])
		k_max = Theta.dot(Pmax)
		QtLZ = (Theta.T.dot(cov.dot(z)))
		f = (Theta.dot(P)-QtLZ)
		return k_max - f
	
	def const_3_typeB_Zeta(self,z):
		M = self.M
		#T = self.temperatures
		T = self.temperatures[-1]
		cov = self.L
		Pmin = self.P_min
		P = self.ArrheniusParams
		Theta = np.array([T/T,np.log(T),-1/(T)])
		k_min = Theta.dot(Pmin)
		QtLZ = (Theta.T.dot(cov.dot(z)))  
		f = (Theta.dot(P)-QtLZ)
		return k_min - f
	
	def const_2_typeC_Zeta(self,z):
		M = self.M
		T = self.temperatures[-1]
		cov = self.L
		Pmax = self.P_max
		P = self.ArrheniusParams
		Theta = np.array([T/T,np.log(T),-1/(T)])
		k_max = Theta.dot(Pmax)
		QtLZ = (Theta.T.dot(cov.dot(z)))
		f = (Theta.dot(P)-QtLZ)
		return k_max - f
	"""
	def cons_derivative_b2(self,z):
		
		if self.kright_fact < 0 and self.kleft_fact < 0:
			T = z[-1]
			guess = z[0:-1]
			P = self.ArrheniusParams
			cov = self.L
			P_max = P + abs(self.kmiddle_fact)*np.asarray(np.dot(cov,self.zeta.x)).flatten()
			theta = np.array([0,1/T,1/T**2])
			dk_dt = theta.T.dot(P_max)
			dko_dt = theta.T.dot(P)
			dQtLZ_dt =(theta.T.dot(cov.dot(z[0:-1])))
			obj = dk_dt - dko_dt - dQtLZ_dt
		
		elif self.kright_fact >0 and self.kleft_fact>0:
			T = z[-1]
			guess = z[0:-1]
			P = self.ArrheniusParams
			cov = self.L
			P_min = P - abs(self.kmiddle_fact)*np.asarray(np.dot(cov,self.zeta.x)).flatten()
			theta = np.array([0,1/T,1/T**2])
			dk_dt = theta.T.dot(P_min)
			dko_dt = theta.T.dot(P)
			dQtLZ_dt =(theta.T.dot(cov.dot(z[0:-1])))
			obj = dk_dt - dko_dt - dQtLZ_dt

		elif self.kright_fact>0 and self.kleft_fact<0:
			
			T = z[-1]
			guess = z[0:-1]
			P = self.ArrheniusParams
			cov = self.L
			P_min = P + abs(self.kmiddle_fact)*np.asarray(np.dot(cov,self.zeta.x)).flatten()
			theta = np.array([0,1/T,1/T**2])
			dk_dt = theta.T.dot(P_min)
			dko_dt = theta.T.dot(P)
			dQtLZ_dt =(theta.T.dot(cov.dot(z[0:-1])))
			obj = dk_dt - dko_dt - dQtLZ_dt
			
			#obj = 0.0
		elif self.kright_fact<0 and self.kleft_fact>0:
			
			T = z[-1]
			guess = z[0:-1]
			P = self.ArrheniusParams
			cov = self.L
			P_max = P - abs(self.kmiddle_fact)*np.asarray(np.dot(cov,self.zeta.x)).flatten()
			theta = np.array([0,1/T,1/T**2])
			dk_dt = theta.T.dot(P_max)
			dko_dt = theta.T.dot(P)
			dQtLZ_dt =(theta.T.dot(cov.dot(z[0:-1])))
			obj = dk_dt - dko_dt - dQtLZ_dt
			
			#obj = 0.0
			
		else:
			
			fact = 0
			T = z[-1]
			guess = z[0:-1]
			P = self.ArrheniusParams
			cov = self.L
			P_max = P - fact*np.asarray(np.dot(cov,self.zeta.x)).flatten()
			theta = np.array([0,1/T,1/T**2])
			dk_dt = theta.T.dot(P_max)
			dko_dt = theta.T.dot(P)
			dQtLZ_dt =(theta.T.dot(cov.dot(z[0:-1])))
			obj = dk_dt - dko_dt - dQtLZ_dt
		return obj
	
	def cons_derivative_b2(self,z):
		
		if self.kright_fact < 0:
			T = z[-1]
			guess = z[0:-1]
			P = self.ArrheniusParams
			cov = self.L
			P_max = P + np.asarray(np.dot(cov,self.zeta.x)).flatten()
			theta = np.array([0,1/T,1/T**2])
			dk_dt = theta.T.dot(P_max)
			dko_dt = theta.T.dot(P)
			dQtLZ_dt =(theta.T.dot(cov.dot(z[0:-1])))
			obj = dk_dt - dko_dt - dQtLZ_dt
		
		elif self.kright_fact >0:
			T = z[-1]
			guess = z[0:-1]
			P = self.ArrheniusParams
			cov = self.L
			P_min = P - np.asarray(np.dot(cov,self.zeta.x)).flatten()
			theta = np.array([0,1/T,1/T**2])
			dk_dt = theta.T.dot(P_min)
			dko_dt = theta.T.dot(P)
			dQtLZ_dt =(theta.T.dot(cov.dot(z[0:-1])))
			obj = dk_dt - dko_dt + dQtLZ_dt
		
		else:
			T = z[-1]
			guess = z[0:-1]
			P = self.ArrheniusParams
			cov = self.L
			P_max = P + np.asarray(np.dot(cov,self.zeta.x)).flatten()
			theta = np.array([0,1/T,1/T**2])
			dk_dt = theta.T.dot(P_max)
			dko_dt = theta.T.dot(P)
			dQtLZ_dt =(theta.T.dot(cov.dot(z[0:-1])))
			obj = dk_dt - dko_dt - dQtLZ_dt
		return obj
	
	def const_2_typeB2_Zeta(self,z):
		
		if self.kright_fact>0:
			M = self.M
			T = z[-1]
			cov = self.L
			P = self.ArrheniusParams
			Pmin = P - np.asarray(np.dot(self.cov,self.zeta.x)).flatten()
			Theta = np.array([1,np.log(T),-1/(T)])
			k_min = Theta.dot(Pmin)
			QtLZ = (Theta.T.dot(cov.dot(z[0:-1])))
			f = (Theta.dot(P)+QtLZ)
			obj = abs(k_min - f)
		
		elif self.kright_fact<0:	
			M = self.M
			T = z[-1]
			cov = self.L
			P = self.ArrheniusParams
			Pmax = P + np.asarray(np.dot(self.cov,self.zeta.x)).flatten()
			Theta = np.array([1,np.log(T),-1/(T)])
			k_max = Theta.dot(Pmax)
			QtLZ = (Theta.T.dot(cov.dot(z[0:-1])))
			f = (Theta.dot(P)+QtLZ)
			obj = abs(k_max - f)
		
		
		else:
			M = self.M
			T = z[-1]
			cov = self.L
			P = self.ArrheniusParams
			Pmin = P - np.asarray(np.dot(self.cov,self.zeta.x)).flatten()
			Theta = np.array([1,np.log(T),-1/(T)])
			k_min = Theta.dot(Pmin)
			QtLZ = (Theta.T.dot(cov.dot(z[0:-1])))
			f = (Theta.dot(P)+QtLZ)
			obj = abs(k_min - f)
		return obj
		
		"""
	
	def cons_derivative_b2(self,z):
		if self.kright_fact < 0 and self.kleft_fact<0:
		    T = z[-1]
		    guess = z[0:-1]
		    P = self.ArrheniusParams
		    cov = self.L
		    P_max = P + np.asarray(np.dot(cov,self.kmiddle_fact*self.zeta.x)).flatten()
		    theta = np.array([0*(T/T),1/T,1/T**2])
		    dk_dt = theta.T.dot(P_max)
		    dko_dt = theta.T.dot(P)
		    dQtLZ_dt =(theta.T.dot(cov.dot(z[0:-1])))
		    obj = dk_dt - dko_dt - dQtLZ_dt
		elif self.kright_fact >0 and self.kleft_fact>0:
		    T = z[-1]
		    guess = z[0:-1]
		    P = self.ArrheniusParams
		    cov = self.L
		    P_min = P - np.asarray(np.dot(cov,self.kmiddle_fact*self.zeta.x)).flatten()
		    theta = np.array([0,1/T,1/T**2])
		    dk_dt = theta.T.dot(P_min)
		    dko_dt = theta.T.dot(P)
		    dQtLZ_dt =(theta.T.dot(cov.dot(z[0:-1])))
		    obj = dk_dt - dko_dt - dQtLZ_dt
		    
		elif self.kright_fact>0 and self.kleft_fact<0:
		    T = z[-1]
		    guess = z[0:-1]
		    P = self.ArrheniusParams
		    cov = self.L
		    P_min = P + self.kmiddle_fact*np.asarray(np.dot(cov,self.zeta.x)).flatten()
		    theta = np.array([0,1/T,1/T**2])
		    dk_dt = theta.T.dot(P_min)
		    dko_dt = theta.T.dot(P)
		    dQtLZ_dt =(theta.T.dot(cov.dot(z[0:-1])))
		    obj = dk_dt - dko_dt - dQtLZ_dt
		elif self.kright_fact<0 and self.kleft_fact>0:
		    T = z[-1]
		    guess = z[0:-1]
		    P = self.ArrheniusParams
		    cov = self.L
		    P_max = P - self.kmiddle_fact*np.asarray(np.dot(cov,self.zeta.x)).flatten()
		    theta = np.array([0,1/T,1/T**2])
		    dk_dt = theta.T.dot(P_max)
		    dko_dt = theta.T.dot(P)
		    dQtLZ_dt =(theta.T.dot(cov.dot(z[0:-1])))
		    obj = dk_dt - dko_dt - dQtLZ_dt
		else:
		    T = z[-1]
		    guess = z[0:-1]
		    P = self.ArrheniusParams
		    cov = self.L
		    P_max = P + self.kmiddle_fact*np.asarray(np.dot(cov,self.zeta.x)).flatten()
		    theta = np.array([0,1/T,1/T**2])
		    dk_dt = theta.T.dot(P_max)
		    dko_dt = theta.T.dot(P)
		    dQtLZ_dt =(theta.T.dot(cov.dot(z[0:-1])))
		    obj = dk_dt - dko_dt - dQtLZ_dt
		return obj
	    
	def const_2_typeB2_Zeta(self,z):
		if self.kleft_fact >0 and self.kright_fact>0:
		    M = self.M
		    T = z[-1]
		    cov = self.L
		    P = self.ArrheniusParams
		    Pmin = P - self.kmiddle_fact*np.asarray(np.dot(self.cov,self.zeta.x)).flatten()
		    Theta = np.array([T/T,np.log(T),-1/(T)]).astype(float)
		    k_min = Theta.dot(Pmin)
		    QtLZ = (Theta.T.dot(cov.dot(z[0:-1])))
		    f = (Theta.dot(P)+QtLZ)
		    obj = k_min - f
		elif self.kleft_fact <0 and self.kright_fact<0:
		    M = self.M
		    T = z[-1]
		    cov = self.L
		    P = self.ArrheniusParams
		    Pmax = P + self.kmiddle_fact*np.asarray(np.dot(self.cov,self.zeta.x)).flatten()
		    Theta = np.array([1,np.log(T),-1/(T)])
		    k_max = Theta.dot(Pmax)
		    QtLZ = (Theta.T.dot(cov.dot(z[0:-1])))
		    f = (Theta.dot(P)+QtLZ)
		    obj = k_max - f
		elif self.kleft_fact <0 and self.kright_fact>0:
		    M = self.M
		    T = z[-1]
		    cov = self.L
		    P = self.ArrheniusParams
		    Pmax = P + self.kmiddle_fact*np.asarray(np.dot(self.cov,self.zeta.x)).flatten()
		    Theta = np.array([1,np.log(T),-1/(T)])
		    k_max = Theta.dot(Pmax)
		    QtLZ = (Theta.T.dot(cov.dot(z[0:-1])))
		    f = (Theta.dot(P)+QtLZ)
		    obj = k_max - f
		    """
		    M = self.M
		    T = self.temperatures[1:-1]
		    cov = self.L
		    Theta = np.array([T/T,np.log(T),-1/(T)])
		    QtLZ = np.asarray([(i.dot(cov.dot(z[0:-1]))) for i in Theta.T])
		    f = (self.uncertainties[1:-1]-QtLZ)
		    obj = np.amin(f)
		    """
		elif self.kleft_fact >0 and self.kright_fact<0:
		    M = self.M
		    T = z[-1]
		    cov = self.L
		    P = self.ArrheniusParams
		    Pmin = P - self.kmiddle_fact*np.asarray(np.dot(self.cov,self.zeta.x)).flatten()
		    Theta = np.array([1,np.log(T),-1/(T)])
		    k_min = Theta.dot(Pmin)
		    QtLZ = (Theta.T.dot(cov.dot(z[0:-1])))
		    f = (Theta.dot(P)+QtLZ)
		    obj = k_min - f
		    """
		    M = self.M
		    T = self.temperatures[1:-1]
		    cov = self.L
		    Theta = np.array([T/T,np.log(T),-1/(T)])
		    QtLZ = np.asarray([(i.dot(cov.dot(z[0:-1]))) for i in Theta.T])
		    f = (self.uncertainties[1:-1]-QtLZ)
		    obj = np.amin(f)
		    """
		else:
		    M = self.M
		    T = z[-1]
		    cov = self.L
		    P = self.ArrheniusParams
		    Pmin = P - self.kmiddle_fact*np.asarray(np.dot(self.cov,self.zeta.x)).flatten()
		    Theta = np.array([1,np.log(T),-1/(T)])
		    k_min = Theta.dot(Pmin)
		    QtLZ = (Theta.T.dot(cov.dot(z[0:-1])))
		    f = (Theta.dot(P)+QtLZ)
		    obj = k_min - f
		return obj

	def const_1_typeB2_Zeta(self,z):
		if self.kleft_fact <0:
			M = self.M
			P = self.ArrheniusParams
			T = self.temperatures[0]
			cov = self.L
			P_left = P - abs(self.kleft_fact)*np.asarray(np.dot(self.cov,self.zeta.x)).flatten()
			Theta = np.array([T/T,np.log(T),-1/(T)])
			k_left = Theta.dot(P_left)
			QtLZ = (Theta.T.dot(cov.dot(z[0:-1])))
			f = Theta.dot(P)+QtLZ
		
		elif self.kleft_fact>0:
			M = self.M
			P = self.ArrheniusParams
			T = self.temperatures[0]
			cov = self.L
			P_left = P + abs(self.kleft_fact)*np.asarray(np.dot(self.cov,self.zeta.x)).flatten()
			Theta = np.array([T/T,np.log(T),-1/(T)])
			k_left = Theta.dot(P_left)
			QtLZ = (Theta.T.dot(cov.dot(z[0:-1])))
			f = Theta.dot(P)+QtLZ
		else:
			M = self.M
			P = self.ArrheniusParams
			#T = self.temperatures
			T = self.temperatures[0]
			cov = self.L
			P_right = P 
			Theta = np.array([T/T,np.log(T),-1/(T)])
			k_left = Theta.dot(P_right)
			QtLZ = (Theta.T.dot(cov.dot(z[0:-1])))
			f = Theta.dot(P)+QtLZ
		return k_left - f
	
	def const_1_typeC2_Zeta(self,z):
		if self.kleft_fact <0:
			M = self.M
			P = self.ArrheniusParams
			T = self.temperatures[0]
			cov = self.L
			P_left = P - abs(self.kleft_fact)*np.asarray(np.dot(self.cov,self.zeta.x)).flatten()
			Theta = np.array([T/T,np.log(T),-1/(T)])
			k_left = Theta.dot(P_left)
			QtLZ = (Theta.T.dot(cov.dot(z[0:-1])))
			f = Theta.dot(P)+QtLZ
		elif self.kleft_fact>0:
			M = self.M
			P = self.ArrheniusParams
			T = self.temperatures[0]
			cov = self.L
			P_left = P + abs(self.kleft_fact)*np.asarray(np.dot(self.cov,self.zeta.x)).flatten()
			Theta = np.array([T/T,np.log(T),-1/(T)])
			k_left = Theta.dot(P_left)
			QtLZ = (Theta.T.dot(cov.dot(z[0:-1])))
			f = Theta.dot(P)+QtLZ
		else:
			M = self.M
			P = self.ArrheniusParams
			#T = self.temperatures
			T = self.temperatures[0]
			cov = self.L
			P_right = P 
			Theta = np.array([T/T,np.log(T),-1/(T)])
			k_left = Theta.dot(P_right)
			QtLZ = (Theta.T.dot(cov.dot(z[0:-1])))
			f = Theta.dot(P)+QtLZ
		return k_left - f
	
	def populateValues(self,a1,a2):
		self.kleft_fact = a1
		self.kright_fact = a2
		self.kmiddle_fact = 1.0
		#print(f"In populate{a1},{a2}\n")
	def const_3_typeB2_Zeta(self,z):
		if self.kright_fact <0:
			M = self.M
			P = self.ArrheniusParams
			T = self.temperatures[-1]
			cov = self.L
			P_right = P - abs(self.kright_fact)*np.asarray(np.dot(self.cov,self.zeta.x)).flatten()
			Theta = np.array([1,np.log(T),-1/(T)])
			k_right = Theta.dot(P_right)
			QtLZ = (Theta.T.dot(cov.dot(z[0:-1])))
			f = Theta.dot(P)+QtLZ
		elif self.kright_fact >0:
			M = self.M
			P = self.ArrheniusParams
			T = self.temperatures[-1]
			cov = self.L
			P_right = P + abs(self.kright_fact)*np.asarray(np.dot(self.cov,self.zeta.x)).flatten()
			Theta = np.array([1,np.log(T),-1/(T)])
			k_right = Theta.dot(P_right)
			QtLZ = (Theta.T.dot(cov.dot(z[0:-1])))
			f = Theta.dot(P)+QtLZ
		else:
			M = self.M
			P = self.ArrheniusParams
			#T = self.temperatures
			T = self.temperatures[0]
			cov = self.L
			P_right = P 
			Theta = np.array([T/T,np.log(T),-1/(T)])
			k_right = Theta.dot(P_right)
			QtLZ = (Theta.T.dot(cov.dot(z[0:-1])))
			f = Theta.dot(P)+QtLZ
		return k_right - f
	
	def const_3_typeC2_Zeta(self,z):
		if self.kright_fact <0:
			M = self.M
			P = self.ArrheniusParams
			T = self.temperatures[-1]
			cov = self.L
			P_right = P - abs(self.kright_fact)*np.asarray(np.dot(self.cov,self.zeta.x)).flatten()
			Theta = np.array([1,np.log(T),-1/(T)])
			k_right = Theta.dot(P_right)
			QtLZ = (Theta.T.dot(cov.dot(z[0:-1])))
			f = Theta.dot(P)+QtLZ
		elif self.kright_fact >0:
			M = self.M
			P = self.ArrheniusParams
			T = self.temperatures[-1]
			cov = self.L
			P_right = P + abs(self.kright_fact)*np.asarray(np.dot(self.cov,self.zeta.x)).flatten()
			Theta = np.array([1,np.log(T),-1/(T)])
			k_right = Theta.dot(P_right)
			QtLZ = (Theta.T.dot(cov.dot(z[0:-1])))
			f = Theta.dot(P)+QtLZ
		else:
			M = self.M
			P = self.ArrheniusParams
			#T = self.temperatures
			T = self.temperatures[0]
			cov = self.L
			P_right = P 
			Theta = np.array([T/T,np.log(T),-1/(T)])
			k_right = Theta.dot(P_right)
			QtLZ = (Theta.T.dot(cov.dot(z[0:-1])))
			f = Theta.dot(P)+QtLZ
		return k_right - f
	
	
	def getCovariance(self,flag = False):
		if flag == True:
			constraints = {'type': 'ineq', 'fun': self.const_func }
			self.const = [constraints]
			START_MUQ = time.time()
			self.solution = minimize(self.obj_func,self.guess,constraints=self.cons)
			STOP_MUQ = time.time()
			#print(f"\nMUQ-SAC method takes {STOP_MUQ-START_MUQ}\n")
		else:
			START_MUQ = time.time()
			self.solution = minimize(self.obj_func,self.guess,method="SLSQP")
			STOP_MUQ = time.time()
			#print(f"\nMUQ-SAC method takes {STOP_MUQ-START_MUQ}\n")
			
		#print(self.solution.x)
		self.L = np.array([[self.solution.x[0],0,0],[self.solution.x[1],self.solution.x[2],0],[self.solution.x[3],self.solution.x[4],self.solution.x[5]]]);#cholesky lower triangular matric
		cov1 = self.L
		cov2 = np.dot(self.L,self.L.T)
		#print(cov2)
		#print(np.exp(self.ArrheniusParams[0]),self.ArrheniusParams[1],self.ArrheniusParams[2])
		#print(self.temperatures[0],self.temperatures[-1])
		D,Q = np.linalg.eigh(cov2)
		self.A = Q.dot(sp.linalg.sqrtm(np.diag(D)))
		self.cov = self.L
	
	def getUnCorrelated(self,flag = False):
		self.getCovariance(flag = False)
		if flag == True:
			#con1 = {'type': 'ineq', 'fun': self.const_func_zeta_1}
			#con2 = {'type': 'ineq', 'fun': self.const_func_zeta_2}
			con3 = {'type': 'eq', 'fun': self.const_nmax}
			con4 = {'type': 'eq', 'fun': self.const_nmin}
			self.const_zeta = [con3,con4]
			zeta = minimize(self.obj_func_zeta,self.guess_z,constraints=self.const_zeta)
		else:
			zeta = minimize(self.obj_func_zeta,self.guess_z,method="Nelder-Mead")

		self.zeta = zeta
		P = self.ArrheniusParams
		self.P_max = P + np.asarray(np.dot(self.L,self.zeta.x)).flatten();
		self.P_min = P - np.asarray(np.dot(self.L,self.zeta.x)).flatten();
		self.kmax = self.Theta.T.dot(self.P_max)
		self.kmin = self.Theta.T.dot(self.P_min)
		self.kappa = self.Theta.T.dot(P)
		return zeta
	
	def getConstrainedUnsrtZeta(self,flag=False):
		if flag == True:
			con1 = {'type': 'eq', 'fun': self.const_1_typeB_Zeta}
			con2 = {'type': 'eq', 'fun': self.const_2_typeB_Zeta}
			con3 = {'type': 'eq', 'fun': self.const_3_typeB_Zeta}
			self.const_zeta = [con1,con2,con3]
			zeta_list = []
			obj_val = []
			alpha = []
			n = []
			epsilon = []
			
			for i,T in enumerate(self.temperatures):
				if T<self.temperatures[-1]-300 and T> self.temperatures[0]+300:
					self.const_T = T
					zeta = minimize(self.obj_func_zeta,self.guess_z,method="SLSQP",constraints=self.const_zeta)
					#zeta = minimize(self.obj_curved_zeta,self.guess_z,method="Nelder-Mead")
					zeta_list.append(zeta.x)
					alpha.append(zeta.x[0])
					n.append(zeta.x[1])
					epsilon.append(zeta.x[2])
					obj_val.append(abs(self.obj_func_zeta(zeta.x))+abs(self.const_1_typeB_Zeta(zeta.x))+abs(self.const_2_typeB_Zeta(zeta.x))+abs(self.const_3_typeB_Zeta(zeta.x)))
			#print(self.rIndex)
			#print(obj_val)
			alpha_square = [i**2 for i in alpha]
			n_square = [i**2 for i in n]
			epsilon_square = [i**2 for i in epsilon]
			index = obj_val.index(min(obj_val))
			
			return [zeta_list[index]],index,np.array([alpha[alpha_square.index(max(alpha_square))],n[n_square.index(max(n_square))],epsilon[epsilon_square.index(max(epsilon_square))]])
		else:
			con1 = {'type': 'eq', 'fun': self.const_1_typeB_Zeta}
			con2 = {'type': 'eq', 'fun': self.const_2_typeC_Zeta}
			con3 = {'type': 'ineq', 'fun': self.const_func_zeta_2}
			self.const_zeta = [con1,con2]
			zeta = minimize(self.obj_func_zeta,self.guess_z,constraints=self.const_zeta)
			return zeta
		#return zeta_list
	
	
	def getC2Zeta(self,flag):
		self.getUnCorrelated(flag=False)
		#self.kleft_fact = 0.1
		#self.kright_fact = -0.5
		if flag == True:
			if self.kleft_fact <0:
				M = self.M
				P = self.ArrheniusParams
				T = self.temperatures[0]
				cov = self.L
				P_left = P - abs(self.kleft_fact)*np.asarray(np.dot(self.cov,self.zeta.x)).flatten()
				Theta = np.array([T/T,np.log(T),-1/(T)])
				k0_left = Theta.dot(P)
				k_left = Theta.dot(P_left)
				
				
			elif self.kleft_fact>0:
				M = self.M
				P = self.ArrheniusParams
				T = self.temperatures[0]
				cov = self.L
				P_left = P + abs(self.kleft_fact)*np.asarray(np.dot(self.cov,self.zeta.x)).flatten()
				Theta = np.array([T/T,np.log(T),-1/(T)])
				k0_left = Theta.dot(P)
				k_left = Theta.dot(P_left)
				
				
			else:
				M = self.M
				P = self.ArrheniusParams
				#T = self.temperatures
				T = self.temperatures[0]
				cov = self.L
				P_right = P 
				Theta = np.array([T/T,np.log(T),-1/(T)])
				k0_left = Theta.dot(P)
				k_left = Theta.dot(P_right)
				
			
			if self.kright_fact <0:
				M = self.M
				P = self.ArrheniusParams
				T = self.temperatures[-1]
				cov = self.L
				P_right = P - abs(self.kright_fact)*np.asarray(np.dot(self.cov,self.zeta.x)).flatten()
				Theta = np.array([1,np.log(T),-1/(T)])
				k0_right = Theta.dot(P)
				k_right = Theta.dot(P_right)
				
				
			elif self.kright_fact >0:
				M = self.M
				P = self.ArrheniusParams
				T = self.temperatures[-1]
				cov = self.L
				P_right = P + abs(self.kright_fact)*np.asarray(np.dot(self.cov,self.zeta.x)).flatten()
				Theta = np.array([1,np.log(T),-1/(T)])
				k0_right = Theta.dot(P)
				k_right = Theta.dot(P_right)
				
				
			else:
				M = self.M
				P = self.ArrheniusParams
				#T = self.temperatures
				T = self.temperatures[0]
				
				cov = self.L
				P_right = P 
				Theta = np.array([T/T,np.log(T),-1/(T)])
				k0_right = Theta.dot(P)
				k_right = Theta.dot(P_right)
				
			
			"""
			Find slope
			------
			To find that we need the f(T) from the kappa_left and kappa_right
			"""
			T2 = self.temperatures[-1]
			FT2 = k_left-k0_left
			T1 = self.temperatures[0]
			FT1 = k_right-k0_right
			
			slope = (FT2 - FT1)/(T2-T1)
			constant = FT2-slope*(T2)
			
			self.Yu = slope*(self.temperatures) + constant
								
			#con1 = {'type': 'eq', 'fun': self.const_1_typeC2_Zeta}
			#con3 = {'type': 'eq', 'fun': self.const_3_typeC2_Zeta}
			#self.const_zeta = [con1,con3]
			#bnds = ((float("-inf"),float("inf")),(float("-inf"),float("inf")),(float("-inf"),float("inf")),(200,3500))
			#zeta = minimize(self.obj_func_zeta_b2,self.guess_z2,method="SLSQP",constraints=self.const_zeta,bounds=bnds)
			zeta = minimize(self.obj_func_zeta_c2,self.guess_z,method="SLSQP")
			#zeta = shgo(self.obj_func_zeta_b2,bnds,constraints=self.const_zeta)
			#zeta = minimize(self.obj_func_zeta_c2,self.guess_z2,method="trust-constr",constraints=self.const_zeta)
			#print(zeta)
		else:
			con5 = {'type': 'ineq', 'fun': self.cons_T}
			self.const_zeta = [con5]
			bnds = ((float("-inf"),float("inf")),(float("-inf"),float("inf")),(float("-inf"),float("inf")),(200,3500))
			zeta = minimize(self.obj_func_zeta_b2,self.guess_z2,method="SLSQP",constraints=self.const_zeta,bounds=bnds)
		return [zeta.x[0],zeta.x[1],zeta.x[2]]
	
	def getB2Zeta(self,flag): 
		self.getUnCorrelated(flag=False)
		if flag == True:
			con1 = {'type': 'eq', 'fun': self.const_1_typeB2_Zeta}
			con2 = {'type': 'eq', 'fun': self.const_2_typeB2_Zeta}
			con3 = {'type': 'eq', 'fun': self.const_3_typeB2_Zeta}
			con4 = {'type': 'eq', 'fun': self.cons_derivative_b2}
			self.const_zeta = [con1,con2,con3,con4]
			bnds = ((float("-inf"),float("inf")),(float("-inf"),float("inf")),(float("-inf"),float("inf")),(200,3500))
			#zeta = minimize(self.obj_func_zeta_b2,self.guess_z2,method="SLSQP",constraints=self.const_zeta,bounds=bnds)
			zeta = shgo(self.obj_func_zeta_b2,bnds,constraints=self.const_zeta)
			#zeta = minimize(self.obj_func_zeta_b2,self.guess_z2,method="trust-constr",constraints=self.const_zeta,bounds=bnds)#,options={'xtol': 1e-05, 'gtol': 1e-05, 'barrier_tol': 1e-05})
			#print(f"{zeta.success}\t{self.kright_fact}\t{self.kleft_fact}\n")
			#if zeta.success == False:
			#	print(zeta)
			#raise AssertionError("stop")
		else:
			con5 = {'type': 'ineq', 'fun': self.cons_T}
			self.const_zeta = [con5]
			bnds = ((float("-inf"),float("inf")),(float("-inf"),float("inf")),(float("-inf"),float("inf")),(200,3500))
			zeta = minimize(self.obj_func_zeta_b2,self.guess_z2,method="SLSQP",constraints=self.const_zeta,bounds=bnds)
		return [zeta.x[0:-1][0],zeta.x[0:-1][1],zeta.x[0:-1][2]]
		
	def obj_get_kappa(self,guess):
		cov = self.L
		opt_kappa = self.kappa
		QtLZ = self.kappa_0 + self.theta_for_kappa.T.dot(cov.dot(guess))
		f = opt_kappa - QtLZ
		obj = np.dot(f,f)
		return obj
		
	def getZeta_typeA(self,kappa):
		T = np.linspace(300,2500,50)
		self.theta_for_kappa = np.array([T/T,np.log(T),-1/T])
		P = self.ArrheniusParams
		#print(P)
		self.kappa_0 = self.theta_for_kappa.T.dot(P)
		#print(self.kappa_0)
		#self.kappa = kappa
		y = kappa - self.kappa_0
		A = self.theta_for_kappa.T
		Q,R = np.linalg.qr(A)
		y_dash = Q.T.dot(y)
		x = np.linalg.solve(R,y_dash.T)
		zeta = np.linalg.solve(self.L,x.T)
		
		#print(self.kappa)
		#K = self.kappa - self.kappa_0
		#theta_inv = np.linalg.inv(self.theta_for_kappa.T)
		#zeta = np.linalg.inv(self.L).dot(theta_inv.dot(K))
		#zeta = minimize(self.obj_get_kappa,self.guess_z,method="Nelder-Mead")
		return zeta
	
	def getZetaFromGen(self,generator):
		P = self.ArrheniusParams
		self.cov = self.getCovariance()
		self.zeta = self.getUnCorrelated(flag=False)
		self.unsrtFunc = self.getUncertFunc(self.cov)
		self.zetaUnsrt = self.getZetaUnsrtFunc(self.cov,self.zeta.x)
		zeta = self.zeta.x
		self.P_max = P + np.asarray(np.dot(self.cov,zeta)).flatten();
		self.P_min = P - np.asarray(np.dot(self.cov,zeta)).flatten();
		self.kmax = self.Theta.T.dot(self.P_max)
		self.kmin = self.Theta.T.dot(self.P_min)
		self.kappa = self.Theta.T.dot(P)
		#self.zeta_curved_type_B,index,zeta_lim = self.getConstrainedUnsrtZeta(flag=True)
		#self.zeta_curved_type_C = self.getConstrainedUnsrtZeta(flag=False)
		#self.zeta_curved_type_B1 = self.getConstrainedUnsrtZeta_typeB1(flag=True)
		self.generator = []
		self.samples = []
		
		self.kleft_fact = generator[0]
		self.kright_fact = generator[1]
		#self.kmiddle_fact = generator[2]
		zeta_B2 = self.getB2Zeta(flag=True)
		return zeta_B2
		
	def getExtremeCurves(self,tag,zeta_type,sample_points):
		P = self.ArrheniusParams
		self.cov = self.getCovariance()
		self.zeta = self.getUnCorrelated(flag=False)
		self.unsrtFunc = self.getUncertFunc(self.cov)
		self.zetaUnsrt = self.getZetaUnsrtFunc(self.cov,self.zeta.x)
		zeta = self.zeta.x
		self.P_max = P + np.asarray(np.dot(self.cov,zeta)).flatten();
		self.P_min = P - np.asarray(np.dot(self.cov,zeta)).flatten();
		self.kmax = self.Theta.T.dot(self.P_max)
		self.kmin = self.Theta.T.dot(self.P_min)
		self.kappa = self.Theta.T.dot(P)
		#self.zeta_curved_type_B,index,zeta_lim = self.getConstrainedUnsrtZeta(flag=True)
		#self.zeta_curved_type_C = self.getConstrainedUnsrtZeta(flag=False)
		#self.zeta_curved_type_B1 = self.getConstrainedUnsrtZeta_typeB1(flag=True)
		self.generator = []
		self.samples = []
		for i in range(int(sample_points)):
		    a1 =2*np.random.random_sample(1)-1
		    a2 = 2*np.random.random_sample(1)-1
		    #a3 = 2*np.random.random_sample(1)-1
		    #a3 = [1.0]
		    #a1 =np.random.uniform(-1,1,1)
		    #a2 = np.random.uniform(-1,1,1)
		    #a3 = np.random.uniform(-1,1,1)
		    
		    generator.append([a1[0],a2[0]])
		    self.kleft_fact = a1[0]
		    self.kright_fact = a2[0]
		    self.kmiddle_fact = abs(a1[0])
		    zeta_B2 = self.getB2Zeta(flag=True)
		    self.samples.append(zeta_B2)
		return self.zeta_B2  
	
	def getExtremeCurves_fast(self,sample_points):
		X = workers(100)
		zeta_list = X.do_job_async(self.data,sample_points)
		del X
		#print(self.parallized_zeta)
		return zeta_list
		
	def getUncorreationMatrix(self,tag):
		T = self.temperatures
		P = self.ArrheniusParams
		self.cov = self.getCovariance()
		self.zeta = self.getUnCorrelated(flag=False)
		self.unsrtFunc = self.getUncertFunc(self.cov)
		self.zetaUnsrt = self.getZetaUnsrtFunc(self.cov,self.zeta.x)
		self.P_max = P + np.asarray(np.dot(self.cov,self.zeta.x)).flatten();
		self.P_min = P - np.asarray(np.dot(self.cov,self.zeta.x)).flatten();
		self.kmax = self.Theta.T.dot(self.P_max)
		self.kmin = self.Theta.T.dot(self.P_min)
		self.kappa = self.Theta.T.dot(P)
		self.zeta_curved_type_B,index,zeta_lim = self.getConstrainedUnsrtZeta(flag=True)
		self.zeta_curved_type_C = self.getConstrainedUnsrtZeta(flag=False)
		zeta = np.array([[self.zeta.x[0],self.zeta.x[1],self.zeta.x[2]],[-self.zeta_curved_type_B[0][0],-self.zeta_curved_type_B[0][1],-self.zeta_curved_type_B[0][2]],[self.zeta_curved_type_C.x[0],self.zeta_curved_type_C.x[1],self.zeta_curved_type_C.x[2]]])
		self.zeta_matrix = np.matrix(zeta)
		"""
		fig, axs = plt.subplots(2, 1, figsize=(15,20))

		for zeta in self.zeta_curved_type_A:
			temp_unsrtFunc = np.asarray([self.M*np.dot(i.T,np.dot(self.cov,zeta)) for i in self.Theta.T])
			axs[0].plot(self.temperatures,temp_unsrtFunc,'k--')
			axs[0].plot(self.temperatures,-temp_unsrtFunc,'k--')
			
		axs[0].set_xlabel('Temperature (K)')
		axs[0].set_ylabel('Uncertainity ($f$)')
		temp_unsrtFunc = np.asarray([np.dot(i.T,np.dot(self.cov,self.zeta.x)) for i in self.Theta.T])
		temp_unsrtFunc_C = np.asarray([np.dot(i.T,np.dot(self.cov,self.zeta_curved_type_C.x)) for i in self.Theta.T])

		axs[0].plot(self.temperatures,temp_unsrtFunc_C,'y--')
		axs[0].plot(self.temperatures,-temp_unsrtFunc_C,'y--')
		axs[0].plot(self.temperatures,self.unsrtFunc,'r--',label='present study (MUQ)')
		axs[0].plot(self.temperatures,-self.unsrtFunc,'r--')
		axs[0].plot(self.temperatures,self.zetaUnsrt,'b-',label='present study (zeta) (MUQ)')
		axs[0].plot(self.temperatures,-self.zetaUnsrt,'b-')
		axs[0].set_ylim(-2*max(self.unsrtFunc),2*max(self.unsrtFunc))
		axs[0].plot(self.temperatures,self.uncertainties,'go',label='Exp. data')
		plt.savefig(f"{tag}.pdf",bbox_inches="tight")
		"""
		
		return self.zeta_matrix,P,self.P_max,self.P_min,self.cov
						
class reaction(UncertaintyExtractor):
	def __init__(self, Element,mechPath,binary_files):
		#self.samap_executable = binary_files["samap_executable"]
		#self.jpdap_executable = binary_files["jpdap_executable"]
		self.rxn = self.classification = self.type = self.sub_type = self.exp_data_type = self.temperatures = self.uncertainties = self.branching  = self.branches = self.pressure_limit = self.common_temp = self.temp_limit = None
		
		DATA = Parser(mechPath).mech
		RXN_LIST = Parser(mechPath).rxnList
		self.tag = Element.tag
		self.rxn = str(Element.attrib["rxn"])
		self.rIndex = str(Element.attrib["no"])
		#self.rxn_dict = IFR.MechParsing(mechPath).getKappa(self.rxn)
		#print(self.rIndex,IFR.MechParsing(mechPath).getArrhenius(self.rxn))
		if self.rxn in RXN_LIST:
			self.index = RXN_LIST.index(self.rxn)
		else:
			raise AssertionError(f"Rxn {self.rxn} not in the mechanism. Kindly check the uncertainty file that you have submitted !!\n")
		
		for item in Element:
			if item.tag == "class":
				self.classification = item.text
			if item.tag == "type":
				self.type = item.text
			if item.tag == "perturbation_type":
				self.perturbation_type = item.text
			if item.tag == "perturbation_factor":
				self.perturbation_factor = float(item.text)
			if item.tag == "sub_type":
				self.sub_type = item.attrib["name"]
				for subitem in item:
					if subitem.tag == "multiple":
						self.multiple = subitem.text
					if subitem.tag == "branching":
						self.branching = subitem.text
					if subitem.tag == "branches":
						self.branches = subitem.text
					if subitem.tag == "pressure_limit":
						self.pressure_limit = subitem.text
					if subitem.tag == "common_temp":
						self.common_temp = subitem.text
					if subitem.tag == "temp_limit":
						self.temp_limit = subitem.text
			
			if item.tag == "data_type":
				self.exp_data_type = item.text
			if item.tag == "file":
				self.exp_data_file = item.text
			if item.tag == "temp":
				#print(item.text)
				self.temperatures = np.asarray([float(i) for i in item.text.split(",")])
			if item.tag == "unsrt":
				self.uncertainties = np.asarray([float(i) for i in item.text.split(",")])
			
		if self.exp_data_type.split(";")[0] == "constant":
			if self.exp_data_type.split(";")[1] == "array":	
				self.temperatures = self.temperatures
				self.uncertainties = self.uncertainties
				
			elif self.exp_data_type.split(";")[1] == "end_points":
				self.temperatures = np.linspace(self.temperatures[0],self.temperatures[1],200)
				self.uncertainties = np.linspace(self.uncertainties[0],self.uncertainties[1],200)
		elif self.exp_data_type.split(";")[0] == "file":
			unsrt_file = open(str(self.exp_data_file),"r").readlines()
			unsrtData = [np.asfarray(i.strip("\n").strip("''").split(","),float) for i in unsrt_file]
			self.temperatures = np.asarray([i[0] for i in unsrtData])
			self.uncertainties = np.asarray([i[1] for i in unsrtData])
		
		if len(self.temperatures) != len(self.uncertainties):
			print("Error in unsrt data for {}".format(self.rxn))
	
		
		if self.type == "pressure_dependent" and self.pressure_limit.strip() != "":
			if self.pressure_limit == "High":
				self.rxn_Details = DATA["reactions"][self.index]
				self.rxn_dict = self.rxn_Details["high-P-rate-constant"]
				self.nominal = [np.log(self.rxn_dict["A"]),self.rxn_dict["b"],self.rxn_dict["Ea"]/1.987]
			else:
				self.rxn_Details = DATA["reactions"][self.index]
				self.rxn_dict = self.rxn_Details["low-P-rate-constant"]
				self.nominal = [np.log(self.rxn_dict["A"]),self.rxn_dict["b"],self.rxn_dict["Ea"]/1.987]
			self.nametag = self.rxn+":"+self.pressure_limit
			
		elif self.type == "pressure_independent" and self.sub_type == "duplicate":
			if self.branches.strip() == "A":
				self.index = self.index
			else:
				self.index = self.index+1
			self.rxn_Details = DATA["reactions"][self.index]
			self.nametag = self.rxn+":"+self.branches
			self.rxn_dict = self.rxn_Details["rate-constant"]
			self.nominal = [np.log(self.rxn_dict["A"]),self.rxn_dict["b"],self.rxn_dict["Ea"]/1.987]
		else:
			self.rxn_Details = DATA["reactions"][self.index]
			self.nametag = self.rxn
			self.rxn_dict = self.rxn_Details["rate-constant"]
			self.nominal = [np.log(self.rxn_dict["A"]),self.rxn_dict["b"],self.rxn_dict["Ea"]/1.987]
		
		if self.branching == "True":
			self.branches = self.branches.strip('"').split(",")
			self.branches = [int(i)-1 for i in self.branches]
			#print(self.branches)
		#print(self.nametag)
		#print(self.classification)
		
		
		
		data = {}
		data["temperatures"] = self.temperatures
		data["uncertainties"] = self.uncertainties
		data["Arrhenius"] = self.nominal
		
		super().__init__(data)
		self.zeta_Matrix,self.P,self.P_max,self.P_min,self.cov = self.getUncorreationMatrix(self.rIndex)
		self.solution = self.zeta
		self.cholskyDeCorrelateMat = self.L
		#print(self.rIndex,self.L.dot(self.L.T))
		self.activeParameters = [self.rIndex+'_A',self.rIndex+'_n',self.rIndex+'_Ea']
		self.perturb_factor = self.zeta.x
		self.selection = [1.0,1.0,1.0]
		#print(self.zeta.x)
		#print(self.L)
		"""
		if "JPDAP" not in os.listdir():
			os.mkdir("JPDAP")
			os.chdir("JPDAP")
			print(f"{self.rIndex}")
			start = time.time()
			self.input_dict = self.getJPDAP()
			stop = time.time()
			print(f"{stop-start}")
		else:	
			os.chdir("JPDAP")
			if str(self.rIndex) in os.listdir():
				print(f"{self.rIndex}")
				print("Uncertainty_analysis is done!!")
				os.chdir(str(self.rIndex))
				self.input_dict = self.readJPDAP()
			else:
				print(f"{self.rIndex}")
				start = time.time()
				self.input_dict = self.getJPDAP()
				stop = time.time()
				print(f"{stop-start}")
			
		os.chdir("..")
		"""
		#print(f"{self.rIndex}")
		#print(f"{self.cholskyDeCorrelateMat}")
		#print(f"{self.zeta.x}")
		if "factor" in self.perturbation_type:
			#self.perturb_factor =  [min(self.uncertainties),0,0]
			#self.solution = 1.0
			#self.getUncorreationMatrix = np.array([[1.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]])
			#self.cholskyDeCorrelateMat = np.array([[1.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]])
			#self.zeta_Matrix = 1
			#self.activeParameters = [self.rIndex+'_A',self.rIndex+'_n',self.rIndex+'_Ea']
			#self.selection = [1.0,0.0,0.0]
			self.perturb_factor =  [min(self.uncertainties)]
			self.solution = 1.0
			#self.getUncorreationMatrix = np.array([[1.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]])
			self.cholskyDeCorrelateMat = np.array([1.0])
			self.zeta_Matrix = 1
			self.activeParameters = [self.rIndex+'_A']
			self.selection = [1.0,0.0,0.0]
		#print(f"{self.rIndex}= {self.zeta_Matrix}")
		"""
		file_unsrt = open("Reaction_detail_nominal.csv","a+")
		string_rates = f"{self.rIndex},"
		for i in self.rxn_dict:
			string_rates+=f"{i},"
		string_rates+="\n"
		file_unsrt.write(string_rates)
		file_unsrt.close()
		
		file_mat = open("cholesky.csv","a+")
		string_cholesky = f"{self.rIndex},"
		for i in list(self.cholskyDeCorrelateMat):
			for j in i:
				string_cholesky+=f"{j},"
		string_cholesky+="\n"
		file_mat.write(string_cholesky)
		file_mat.close()
		
		file_zeta = open("rxn_zeta_data.csv","a+")
		string_zeta = f"{self.rIndex},"
		for i in self.zeta.x:
			string_zeta+=f"{i},"
		string_zeta+="\n"
		file_zeta.write(string_zeta)
		file_zeta.close()
		"""	
#public function to get the uncertainity values for discrete temperatures
	def getDtList(self):
		self.Unsrt_dict = {}
		key = ["Tag","Solution","Class","Type","Perturbation_type","Sub_type","Branch_boolen","Branches","Pressure_limit","Common_temp","temp_list","Nominal","Exp_input_data_type","priorCovariance","Perturb_factor","Basis_vector","Uncertainties","Temperatures","unsrt_func","Data_key"]
		values = [self.tag,self.solution,self.classification,self.type,self.perturbation_type,self.sub_type,self.branching,self.branches,self.pressure_limit,self.common_temp,self.temp_limit,self.rxn_dict,self.exp_data_type,self.cholskyDeCorrelateMat,self.perturb_factor,self.zeta.x,self.uncertainties,self.temperatures,self.unsrtFunc,""]
		for i,element in enumerate(key):
			self.Unsrt_dict[element] = values[i]
		return self.Unsrt_dict
	
	def getKappaMax(self,T):
		T = np.asarray(T)
		theta = np.array([T/T,np.log(T),-1/T])
		return np.asarray(theta.T.dot(self.P_max)).flatten()
	
	def readJPDAP(self):
		input_dict = {}
		#input_dict["samples"] = self.sim_
		#input_dict["samples_skipped"] = int(0.1*self.sim_)
		#input_dict["Random_seed"] = 1
		#input_dict["sampling_method"] = "SOBOL"
		#input_dict["sampling_distribution"] = "NORMAL"
		#input_dict["equidistant_T"] = 100
		#input_dict["T_begin"] = data = self.rxnUnsert[i].temperatures[0]
		#input_dict["T_end"] = self.rxnUnsert[i].temperatures[-1]
		input_dict["L"] = 0
		input_dict["len_temp_data"] = len(self.temperatures)
		string_unsrt_data =""
		for index,k in enumerate(self.temperatures):
			string_unsrt_data+=f"{k} {self.uncertainties[index]} \n"
		input_dict["temperature_unsrt_data"] = string_unsrt_data
		input_dict["alpha"] = np.exp(self.rxn_dict[0])
		input_dict["n"] = self.rxn_dict[1]
		input_dict["n_min"] = self.rxn_dict[1]-2
		input_dict["n_max"] = self.rxn_dict[1]+2
		input_dict["epsilon"] = self.rxn_dict[2]
		L = self.cholskyDeCorrelateMat
		#input_dict["covariance_matrix"] = str(L.dot(L.T)).strip("[]").replace("[","").replace("]","")
		input_dict["uncertainty_type"] = "2slnk"
		if len(self.activeParameters) == 3:
			input_dict["uncertain_parameters"] = "AnE"
		else:
			input_dict["uncertain_parameters"] = "A"
		Nagy_covariance_matrix = ""
		file_name_jpdap = "jpdap_data_"+str(self.rIndex)+".txt_fit_minRMSD.txt"
		Nagy_covariance_matrix = open(file_name_jpdap,"r").readlines()[5:8]
		covariance_matrix = ""
		cov_float = []
		for n in Nagy_covariance_matrix:
			covariance_matrix+=str(n)
			cov_float.append(np.asfarray(n.strip("''").strip("\n").split(),float))
		#print(covariance_matrix)
		#print(cov_float)
		#print(type(np.asarray(cov_float)))
		os.chdir("..")
		input_dict["cov_float"] = np.asarray(cov_float)
		input_dict["covariance_matrix"] = covariance_matrix.strip("\n")
		return input_dict
	
	def getJPDAP(self):
		input_dict = {}
		#input_dict["samples"] = self.sim_
		#input_dict["samples_skipped"] = int(0.1*self.sim_)
		#input_dict["Random_seed"] = 1
		#input_dict["sampling_method"] = "SOBOL"
		#input_dict["sampling_distribution"] = "NORMAL"
		#input_dict["equidistant_T"] = 100
		#input_dict["T_begin"] = data = self.rxnUnsert[i].temperatures[0]
		#input_dict["T_end"] = self.rxnUnsert[i].temperatures[-1]
		input_dict["L"] = 0
		input_dict["len_temp_data"] = len(self.temperatures)
		string_unsrt_data =""
		for index,k in enumerate(self.temperatures):
			string_unsrt_data+=f"{k} {self.uncertainties[index]} \n"
		input_dict["temperature_unsrt_data"] = string_unsrt_data
		input_dict["alpha"] = np.exp(self.rxn_dict[0])
		input_dict["n"] = self.rxn_dict[1]
		input_dict["n_min"] = self.rxn_dict[1]-2
		input_dict["n_max"] = self.rxn_dict[1]+2
		input_dict["epsilon"] = self.rxn_dict[2]
		L = self.cholskyDeCorrelateMat
		#input_dict["covariance_matrix"] = str(L.dot(L.T)).strip("[]").replace("[","").replace("]","")
		input_dict["uncertainty_type"] = "2slnk"
		if len(self.activeParameters) == 3:
			input_dict["uncertain_parameters"] = "AnE"
		else:
			input_dict["uncertain_parameters"] = "A"
		#input_rxn_dict[i] = input_dict
		string_dict = {}
		jpdap_instring = make_input_file.create_JPDAP_input(input_dict)
		"""
		Run: JPDAP code
		"""
		os.mkdir(f"{self.rIndex}")
		os.chdir(f"{self.rIndex}")
		file_jpdap = open("jpdap_data_"+str(self.rIndex)+".txt","w").write(jpdap_instring)
		run_jpdap_string = f"""#!/bin/bash
{self.jpdap_executable} jpdap_data_{self.rIndex}.txt &> out"""
		file_print_run_jpdap = open("run_jpdap","w").write(run_jpdap_string)
		subprocess.call(["chmod","+x",'run_jpdap'])
		start_Jpdap = time.time()
		subprocess.call(["./run_jpdap"])
		stop_Jpdap = time.time()
		print(f"\n\tJPDAP code took {stop_Jpdap-start_Jpdap}s to execute\n")
		Nagy_covariance_matrix = ""
		file_name_jpdap = "jpdap_data_"+str(self.rIndex)+".txt_fit_minRMSD.txt"
		Nagy_covariance_matrix = open(file_name_jpdap,"r").readlines()[5:8]
		covariance_matrix = ""
		for n in Nagy_covariance_matrix:
			covariance_matrix+=str(n)
		#print(covariance_matrix)
		os.chdir("..")
		input_dict["covariance_matrix"] = covariance_matrix.strip("\n")
		return input_dict
	
	def getKappaMin(self,T):
		T = np.asarray(T)
		theta = np.array([T/T,np.log(T),-1/T])
		return np.asarray(theta.T.dot(self.P_min)).flatten()
	
	def getMean(self):
		return self.P
	
	def getNominal(self,T):
		T = np.asarray(T)
		theta = np.array([T/T,np.log(T),-1/T])
		return np.asarray(theta.T.dot(self.P)).flatten()
	
	def getCov(self):
		return self.cov
	
	def getAllData(self):
		if len(self.branches.split(","))==1:
			b1 = self.branches.split(",")[0]
			b2 = ""
			b3 = ""
		elif len(self.branches.split(",")) >1:
			b1 = self.branches.split(",")[0]
			b2 = self.branches.split(",")[1]
			if len(self.branches.split(",")) >2:
				b3 = self.branches.split(",")[3]
			else:
				b3 = ""
		else:
			b1 = ""
			b2 = ""
			b3 = ""
		exp_data_type = self.exp_data_type.split(";")[0]
		exp_format = self.exp_data_type.split(";")[1]
		Log_string = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.classification,self.type,self.sub_type,self.branching,b1,b2,b3,self.pressure_limit,self.common_temp,self.temp_limit,exp_data_type,exp_format,self.nametag)
		exp_unsrt_string = ""
		solver_log = "######################\n{}######################\n\t\t{}\n\n".format(self.nametag,self.solution)
		calc_cholesky = self.cholskyDeCorrelateMat
		zeta_string = "{}".format(self.zeta)
		for i in range(len(self.temperatures)):
			exp_unsrt_string += "{}\t{}\n".format(self.temperatures[i],self.uncertainties[i])
		return Log_string,exp_unsrt_string,solver_log,calc_cholesky,zeta_string
 
#public function to get the temperature range for the uncertainity of a perticular reaction
	def getTempRange(self):
		return self.temperatures[0], self.temperatures[-1]
	
	def getTemperature(self):
		return self.temperatures
	
	def getRxnType(self):
		return self.type,self.branching,self.branches

#public function to get the zeta values for a perticulat reaction
	def getData(self):	
		return self.zeta.x
	
	def zetaValues(self):
		return self.zeta.x
#public function to get the cholesky decomposed matrix for normalization of variables
		
	def getCholeskyMat(self):
		return self.cholskyDeCorrelateMat


#############################################################
###	   center broadening factors for showing		######
###		  decomposition and recombination of		###### 
###		  pressure dependent reactions			  ######
#############################################################
class PLOG_Interconnectedness(UncertaintyExtractor):
	def __init__(self,plog_object_list,index,count):
		"""
		Creates the class for the PLOG reactions
		
		"""
		
		parent = plog_object_list[0]
		parent2 = plog_object_list[1]
		self.rxn = parent.rxn
		self.classification = parent.classification
		self.tag = parent.tag
		self.rIndex = str(parent.rIndex.split(":")[0])+":"+str(index)
		#print(self.rxn)
		fraction = int(count+1)
		alpha = float(int(index)/int(fraction))
		#print(self.rIndex)
		self.index = parent.index
		#print(self.index)
		self.rxn_Details = parent.rxn_Details
		#self.perturbation_factor = parent.perturbation_factor
		self.perturbation_type = parent.perturbation_type
		
		self.type = parent.type
		self.sub_type = parent.sub_type
		self.exp_data_type = "Interpolation"
		self.temperatures = parent.temperatures
		self.uncertainties = parent.uncertainties
		self.branching  = parent.branching
		self.branches = parent.branches
		self.pressure_limit = "PLOG_"+str(index)
		self.common_temp = parent.common_temp
		self.temp_limit = parent.temp_limit
		self.rxn_dict = self.rxn_Details["rate-constants"][index]
		self.nominal = [np.log(self.rxn_dict["A"]),self.rxn_dict["b"],self.rxn_dict["Ea"]/1.987]
		
		
		
		if self.type == "pressure_dependent" and self.pressure_limit.strip() != "" and self.classification == "PLOG":
			#print(self.pressure_limit)
			self.nametag = str(self.rxn)+":"+str(self.pressure_limit)
		
		elif self.type == "pressure_dependent" and self.pressure_limit.strip() != "" and self.classification == "PLOG-Duplicate":
			#print(self.pressure_limit)
			self.nametag = str(self.rxn)+":"+str(self.branches)+"-"+str(self.pressure_limit)
		
		elif self.type == "pressure_independent" and self.sub_type == "duplicate":
			self.nametag = str(self.rxn)+":"+str(self.branches)
		else:
			self.nametag = self.rxn
		
		#print(self.nametag)
		#print(self.classification)
		"""
		Interpolation for uncertainty 
		"""
		
		
		data = {}
		data["temperatures"] = self.temperatures
		data["uncertainties"] = self.uncertainties
		data["Arrhenius"] = self.nominal
		
		LOW = None
		HIGH = None
		if "High" in parent.rIndex:
			HIGH = parent.uncertainties
			LOW = parent2.uncertainties
		else:
			HIGH = parent2.uncertainties
			LOW = parent.uncertainties
		#print(alpha)
		interpolation = alpha*LOW + (1-alpha)*HIGH
		#print(LOW)
		#print(HIGH)
		#print(interpolation)
		super().__init__(data)
		self.zeta_Matrix,self.P,self.P_max,self.P_min,self.cov = self.getUncorreationMatrix(self.rIndex)
		self.solution = self.zeta.x
		self.cholskyDeCorrelateMat = self.L
		self.activeParameters = [self.rIndex+'_A',self.rIndex+'_n',self.rIndex+'_Ea']
		self.selection = [1.0,1.0,1.0]
		self.perturb_factor = self.zeta.x
		
		"""
		For perturbing only A-factor
		"""
		if "factor" in self.perturbation_type:
			#self.perturb_factor =  [min(self.uncertainties),0,0]
			#self.solution = 1.0
			#self.getUncorreationMatrix = np.array([[1.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]])
			#self.cholskyDeCorrelateMat = np.array([[1.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]])
			#self.zeta_Matrix = 1
			#self.activeParameters = [self.rIndex+'_A',self.rIndex+'_n',self.rIndex+'_Ea']
			#self.selection = [1.0,0.0,0.0]
			self.perturb_factor =  [min(self.uncertainties)]
			self.solution = 1.0
			#self.getUncorreationMatrix = np.array([[1.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]])
			self.cholskyDeCorrelateMat = np.array([1.0])
			self.zeta_Matrix = 1
			self.activeParameters = [self.rIndex+'_A']
			self.selection = [1.0,0.0,0.0]
		
	def getKappaMin(self,T):
		T = np.asarray(T)
		theta = np.array([T/T,np.log(T),-1/T])
		return np.asarray(theta.T.dot(self.P_min)).flatten()
	
	def getMean(self):
		return self.P
	
	def getNominal(self,T):
		T = np.asarray(T)
		theta = np.array([T/T,np.log(T),-1/T])
		return np.asarray(theta.T.dot(self.P)).flatten()
	
	def getCov(self):
		return self.cholskyDeCorrelateMat
	
	def getKappaMax(self,T):
		T = np.asarray(T)
		theta = np.array([T/T,np.log(T),-1/T])
		return np.asarray(theta.T.dot(self.P_max)).flatten()
		
class PLOG(UncertaintyExtractor):
		
	def __init__(self, Element,mechPath,binary_files):
		"""
		for both types of PLOG
		PLOG 
		
		PLOG-DUPLICATE
		
		
		"""
		self.rxn = self.classification = self.type = self.sub_type = self.exp_data_type = self.temperatures = self.uncertainties = self.branching  = self.branches = self.pressure_limit = self.common_temp = self.temp_limit = None
		#super().__init__(Element,mechPath,binary_files)
		DATA = Parser(mechPath).mech
		RXN_LIST = Parser(mechPath).rxnList
		
		self.tag = Element.tag
		self.rxn = str(Element.attrib["rxn"])
		#print(self.rxn)
		self.rIndex = str(Element.attrib["no"])
		
		if self.rxn in RXN_LIST:
			self.index = RXN_LIST.index(self.rxn)
		else:
			raise AssertionError("Rxn not in the mechanism. Kindly check the uncertainty file that you have submitted !!\n")
		
		#self.rxn_Details = DATA["reactions"][self.index]
		#print(self.rxn_Details)
		
		#self.rxn_dict = DATA["Reactions"]["rate-constants"]		
		for item in Element:
			if item.tag == "class":
				self.classification = item.text
			if item.tag == "type":
				self.type = item.text
			if item.tag == "perturbation_type":
				self.perturbation_type = item.text
			if item.tag == "perturbation_factor":
				self.perturbation_factor = float(item.text)
			if item.tag == "sub_type":
				self.sub_type = item.attrib["name"]
				for subitem in item:
					if subitem.tag == "multiple":
						self.branching = subitem.text
					if subitem.tag == "branches":
						self.branches = subitem.text
					if subitem.tag == "pressure_limit":
						self.pressure_limit = subitem.text.strip()
					if subitem.tag == "common_temp":
						self.common_temp = subitem.text
					if subitem.tag == "temp_limit":
						self.temp_limit = subitem.text
			
			if item.tag == "data_type":
				self.exp_data_type = item.text
			if item.tag == "file":
				self.exp_data_file = item.text
			if item.tag == "temp":
				#print(item.text)
				self.temperatures = np.asarray([float(i) for i in item.text.split(",")])
			if item.tag == "unsrt":
				self.uncertainties = np.asarray([float(i) for i in item.text.split(",")])
			
		if self.exp_data_type.split(";")[0] == "constant":
			if self.exp_data_type.split(";")[1] == "array":	
				self.temperatures = self.temperatures
				self.uncertainties = self.uncertainties
				
			elif self.exp_data_type.split(";")[1] == "end_points":
				self.temperatures = np.linspace(self.temperatures[0],self.temperatures[1],200)
				self.uncertainties = np.linspace(self.uncertainties[0],self.uncertainties[1],200)
		elif self.exp_data_type.split(";")[0] == "file":
			unsrt_file = open(str(self.exp_data_file),"r").readlines()
			unsrtData = [np.asfarray(i.strip("\n").strip("''").split(","),float) for i in unsrt_file]
			self.temperatures = np.asarray([i[0] for i in unsrtData])
			self.uncertainties = np.asarray([i[1] for i in unsrtData])
		
		if len(self.temperatures) != len(self.uncertainties):
			print("Error in unsrt data for {}".format(self.rxn))
	
		if self.type == "pressure_dependent" and self.pressure_limit.strip() != "" and self.classification == "PLOG":
			#print(self.pressure_limit)
			self.nametag = str(self.rxn)+":"+str(self.pressure_limit)
		
		elif self.type == "pressure_dependent" and self.pressure_limit.strip() != "" and self.classification == "PLOG-Duplicate":
			#print(self.pressure_limit)
			self.nametag = str(self.rxn)+":"+str(self.branches)+"-"+str(self.pressure_limit)
		
		elif self.type == "pressure_independent" and self.sub_type == "duplicate":
			self.nametag = str(self.rxn)+":"+str(self.branches)
		else:
			self.nametag = self.rxn
		
		#print(self.nametag)
		#print(self.classification)
		if self.classification == "PLOG-Duplicate":
			if self.branches == "A":
				self.index = self.index
			else:
				self.index = self.index+1
			
			self.rxn_Details = DATA["reactions"][self.index]
			if self.pressure_limit == "Low":
				self.rxn_dict = DATA["reactions"][self.index]["rate-constants"][0]
			elif self.pressure_limit == "High":
				self.rxn_dict = DATA["reactions"][self.index]["rate-constants"][-1]
			else:
				raise AssertionError(f"Please give a valid input for pressure_limit in the uncertainty file for PLOG!! The reaction is question is \n{self.rxn_details}\n")
			#print(self.rxn_dict)
			self.nominal = [np.log(self.rxn_dict["A"]),self.rxn_dict["b"],self.rxn_dict["Ea"]/1.987]
			
			
		else:
			self.rxn_Details = DATA["reactions"][self.index]
			if self.pressure_limit == "Low":
				self.rxn_dict = DATA["reactions"][self.index]["rate-constants"][0]
			elif self.pressure_limit == "High":
				self.rxn_dict = DATA["reactions"][self.index]["rate-constants"][-1]
			else:
				raise AssertionError(f"Please give a valid input for pressure_limit in the uncertainty file for PLOG!! The reaction is question is \n{self.rxn_details}\n")
				
			#print(self.rxn_dict)
			self.nominal = [np.log(self.rxn_dict["A"]),self.rxn_dict["b"],self.rxn_dict["Ea"]/1.987]
		
		
		data = {}
		data["temperatures"] = self.temperatures
		data["uncertainties"] = self.uncertainties
		data["Arrhenius"] = self.nominal
		
		"""
		For perturbing all the three reactions
		"""
		super().__init__(data)
		self.zeta_Matrix,self.P,self.P_max,self.P_min,self.cov = self.getUncorreationMatrix(self.rIndex)
		self.solution = self.zeta
		self.cholskyDeCorrelateMat = self.L
		self.activeParameters = [self.rIndex+'_A',self.rIndex+'_n',self.rIndex+'_Ea']
		self.perturb_factor = self.zeta.x
		self.selection = [1.0,1.0,1.0]
		"""
		For perturbing only A-factor
		"""
		if "factor" in self.perturbation_type:
			#self.perturb_factor =  [min(self.uncertainties),0,0]
			#self.solution = 1.0
			#self.getUncorreationMatrix = np.array([[1.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]])
			#self.cholskyDeCorrelateMat = np.array([[1.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]])
			#self.zeta_Matrix = 1
			#self.activeParameters = [self.rIndex+'_A',self.rIndex+'_n',self.rIndex+'_Ea']
			#self.selection = [1.0,0.0,0.0]
			self.perturb_factor =  [min(self.uncertainties)]
			self.solution = 1.0
			#self.getUncorreationMatrix = np.array([[1.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]])
			self.cholskyDeCorrelateMat = np.array([1.0])
			self.zeta_Matrix = 1
			self.activeParameters = [self.rIndex+'_A']
			self.selection = [1.0,0.0,0.0]
		#print(f"{self.rIndex}= {self.zeta_Matrix}")
		
		"""
		file_unsrt = open("Reaction_detail_nominal.csv","a+")
		string_rates = f"{self.rIndex},"
		for i in self.rxn_dict:
			string_rates+=f"{i},"
		string_rates+="\n"
		file_unsrt.write(string_rates)
		file_unsrt.close()
		
		file_mat = open("cholesky.csv","a+")
		string_cholesky = f"{self.rIndex},"
		for i in list(self.cholskyDeCorrelateMat):
			for j in i:
				string_cholesky+=f"{j},"
		string_cholesky+="\n"
		file_mat.write(string_cholesky)
		file_mat.close()
		
		file_zeta = open("rxn_zeta_data.csv","a+")
		string_zeta = f"{self.rIndex},"
		for i in self.zeta.x:
			string_zeta+=f"{i},"
		string_zeta+="\n"
		file_zeta.write(string_zeta)
		file_zeta.close()
		"""
	def getKappaMin(self,T):
		T = np.asarray(T)
		theta = np.array([T/T,np.log(T),-1/T])
		return np.asarray(theta.T.dot(self.P_min)).flatten()
	
	def getMean(self):
		return self.P
	
	def getNominal(self,T):
		T = np.asarray(T)
		theta = np.array([T/T,np.log(T),-1/T])
		return np.asarray(theta.T.dot(self.P)).flatten()
	
	def getCov(self):
		return self.cholskyDeCorrelateMat
	
	def getKappaMax(self,T):
		T = np.asarray(T)
		theta = np.array([T/T,np.log(T),-1/T])
		return np.asarray(theta.T.dot(self.P_max)).flatten()
		
class fallOffCurve:
	def __init__(self, Element,mechPath):
		self.rxn = self.classification = self.type = self.sub_type = self.exp_data_type = self.temperatures = self.uncertainties = None
		
		self.rxn = Element.attrib["rxn"]
		self.tag = Element.tag	
		self.foc_dict = IFR.MechParsing(mechPath).getFocData(self.rxn)[0]
		
		
		#print(self.rxn_dict)
		for item in Element:
			if item.tag == "class":
				self.classification = item.text
			if item.tag == "type":
				self.type = item.text
			if item.tag == "sub_type":
				self.sub_type = item.attrib["name"]
				for subitem in item:
					if subitem.tag == "multiple":
						self.branching = subitem.text
					if subitem.tag == "branches":
						self.branches = subitem.text
					if subitem.tag == "pressure_limit":
						self.pressure_limit = subitem.text
					if subitem.tag == "common_temp":
						self.common_temp = subitem.text
					if subitem.tag == "temp_limit":
						self.temp_limit = subitem.text
			
			if item.tag == "data_type":
				self.exp_data_type = item.text
			if item.tag == "temp":
				#print(item.text)
				self.temperatures = np.asarray([float(i) for i in item.text.split(",")])
			if item.tag == "unsrt":
				self.uncertainties = np.asarray([float(i) for i in item.text.split(",")])
			
		if self.exp_data_type.split(";")[0] == "constant":
			if self.exp_data_type.split(";")[1] == "array":	
				self.temperatures = self.temperatures
				self.uncertainties = self.uncertainties
				
			elif self.exp_data_type.split(";")[1] == "end_points":
				self.temperatures = np.linspace(self.temperatures[0],self.temperatures[1],20)
				self.uncertainties = np.linspace(self.uncertainties[0],self.uncertainties[1],20)
		
		
#		for item in Element:
#			if item.tag == "class":
#				self.classification = item.text
#			if item.tag == "type":
#				self.type = item.text
#			if item.tag == "sub_type":
#				self.sub_type = item.attrib["name"]
#				for subitem in item:
#					if subitem.tag == "branch":
#						self.branching = subitem.text
#					if subitem.tag == "branches":
#						self.branches = subitem.text
#					if subitem.tag == "pressure_limit":
#						self.pressure_limit = subitem.text
#					if subitem.tag == "common_temp":
#						self.common_temp = subitem.text
#					if subitem.tag == "temp_limit":
#						self.temp_limit = subitem.text
#			
#			if item.tag == "data_type":
#				self.exp_data_type = item.text
#			if item.tag == "temp":
#				self.temperatures = np.asarray([float(i) for i in item.text.split(",")])
#			if item.tag == "unsrt":
#				self.uncertainties = np.asarray([float(i) for i in item.text.split(",")])
#			
#		if self.exp_data_type.split(";")[0] == "constant":
#			if self.exp_data_type.split(";")[1] == "array":
#				self.temperatures = self.temperatures
#				self.uncertainties = self.uncertainties
#			elif self.exp_data_type.split(";")[1] == "end_points":
#				self.temperatures = np.linspace(self.temperatures[0],self.temperatures[1],20)
#				self.uncertainties = np.linspace(self.uncertainties[0],self.uncertainties[1],20)
		
		self.nametag = self.rxn+":"+self.sub_type

		L11 = self.uncertainties[0]/3
		L = np.array([L11]) #Cholesky_lower_triangular_matrix
		self.cholskyDeCorrelateMat = L
		#print(L)
		self.zeta = 1
		
		self.T = self.temperatures
		self.f = self.getUncertainty(self.T)		
		self.solution = self.cholskyDeCorrelateMat
		
		
#public function to get the uncertainity values for discrete temperatures
	def getAllData(self):
		Log_string = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.classification,self.type,self.sub_type,self.branching,self.branches,self.pressure_limit,self.common_temp,self.temp_limit,self.exp_data_type,self.nametag)
		exp_unsrt_string = ""
		solver_log = "{}\n{}\n".format(self.nametag,self.solution)
		calc_cholesky = self.cholskyDeCorrelateMat
		zeta_string = "{}".format(self.zeta)
		for i in range(len(self.temperatures)):
			exp_unsrt_string += "{}\t{}\n".format(self.temperatures[i],self.uncertainties[i])
		return Log_string,exp_unsrt_string,solver_log,calc_cholesky,zeta_string

	def getDtList(self):
		self.Unsrt_dict = {}
		key = ["Tag","Solution","Class","Type","Sub_type","Branch_boolen","Branches","Pressure_limit","Common_temp","temp_list","Nominal","Exp_input_data_type","priorCovariance","Basis_vector","Uncertainties","Temperatures","unsrt_func","Data_key"]
		values = [self.tag,self.solution,self.classification,self.type,self.sub_type,self.branching,self.branches,self.pressure_limit,self.common_temp,self.temp_limit,self.foc_dict,self.exp_data_type,self.cholskyDeCorrelateMat,self.zeta,self.uncertainties,self.temperatures,self.f,""]
		for i,element in enumerate(key):
			self.Unsrt_dict[element] = values[i]
		return self.Unsrt_dict
	
	def getUncertainty(self,T):
		L11 = self.uncertainties[0]/3
		
		Foc_unsrt = 3*np.sqrt((L11*(T/T))**2)
		return Foc_unsrt
 
#public function to get the temperature range for the uncertainity of a perticular reaction
	def getTempRange(self):
		return self.temperatures[0], self.temperatures[-1]
	def getTemperature(self):
		return self.temperatures
	def getRxnType(self):
		
		return self.type,self.IsBranching,self.branchRxn
	def fit_zeta(self,T,z):
		L11 = self.uncertainties[0]/3
		func =z*L11*(T/T)
		
		return func

	def getZeta(self):
		#printZeta = "{}\t{}\t{}\n".format(self.zeta.x[0],self.zeta.x[1],self.zeta.x[2])
		#printL ="{}".format(self.cholskyDeCorrelateMat)
		#printL+="\n"
		#fileZeta = open('../zetaValues.txt','a')
		#fileL = open('../Cholesky.txt','a')
		#fileZeta.write(printZeta)
		#fileL.write(printL)
		#fileZeta.close()
		#fileL.close()
		#print("\n{}\n".format(self.zeta.x));
		return self.zeta
	
	def zetaValues(self):
		return self.zeta
#public function to get the cholesky decomposed matrix for normalization of variables
	
	def getCholeskyMat(self):
		return self.cholskyDeCorrelateMat

	def getKappaMax(self,T):
		T = np.asarray(T)
		theta = np.array([T/T,np.log(T),-1/T])
		return np.asarray(theta.T.dot(self.P_max)).flatten()

#############################################################
###	   Uncertainty for heat capacities			  ######
###	   of kinetic species						   ######
#############################################################


class thermodynamic:
	def __init__(self, Element,thermo_loc):
		self.species = self.classification = self.type = self.sub_type =  self.branching = self.branches  = self.pressure_limit  = self. common_temp = self.temp_limit= None
		self.tag = Element.tag
		IFRT = IFR.ThermoParsing(thermo_loc)
				
		self.exp_data_type = {}
		self.temperatures = {}
		self.uncertainties = {}
		self.cholskyDeCorrelateMat = {}
		self.zeta = {}
		self.species = Element.attrib["species"]
		self.nominal = {}
		
		
		for item in Element:
			#print(item.tag)
			if item.tag == "class":
				self.classification = item.text
				continue
			if item.tag == "type":
				self.type = item.text
				continue
			if item.tag == "sub_type":
				self.sub_type = item.attrib["name"]
				
				for subitem in item:
					if "multiple" in subitem.tag:
						self.branching = str(subitem.text)
						continue
					if "branches" in subitem.tag:
						self.branches = str(subitem.text)
						continue
						
					if "pressure_limit" in subitem.tag:
						self.pressure_limit = str(subitem.text)
						continue
					if "common_temp" in subitem.tag:
						self.common_temp = str(subitem.text)
						continue
					if  "temp_limit" in subitem.tag :
						self.temp_limit = str(subitem.text)
						#print(self.temp_limit)
						if self.temp_limit == "Low":
		
							self.thermo_dict = IFR.ThermoParsing(thermo_loc).getThermoLow(self.species)
						else:
							self.thermo_dict = IFR.ThermoParsing(thermo_loc).getThermoHigh(self.species)
						continue
				continue
				
		for i in self.branches.split(","):
			
			for item in Element:
				if item.tag == "data_type_"+str(i):
					self.exp_data_type[str(i)] = item.text
					continue
					
				if item.tag == "temp_"+str(i):
					#print(item.text)
					self.temperatures[str(i)] = np.asarray([float(j) for j in item.text.split(",")])
					continue
						
				if item.tag == "unsrt_"+str(i):
					self.uncertainties[str(i)] = np.asarray([float(i) for i in item.text.split(",")])
					continue
			
		for i in self.exp_data_type:
			if self.exp_data_type[str(i)].split(";")[0] ==  "constant":
				if self.exp_data_type[str(i)].split(";")[1] == "array":
					continue
				elif self.exp_data_type[str(i)].split(";")[1] == "end_points":
					self.temperatures[str(i)] = np.linspace(self.temperatures[str(i)][0],self.temperatures[str(i)][1],20)
					self.uncertainties[str(i)] = np.linspace(self.uncertainties[str(i)][0],self.uncertainties[str(i)][1],20)
					continue

			elif self.exp_data_type[str(i)].split(";")[0] == "percentage":
				if self.exp_data_type[str(i)].split(";")[1] == "array":
					a = self.thermo_dict[str(i)]
					func = IFRT.function(str(i),a,self.temperatures[str(i)])
					y = self.uncertainties[str(i)]
					#print(y)						
					self.uncertainties[str(i)] = np.asarray(np.dot(y,func)).flatten()
					continue

				elif self.exp_data_type[str(i)].split(";")[1] == "end_points":

					self.temperatures[str(i)] = np.linspace(self.temperatures[str(i)][0],self.temperatures[str(i)][1],20)
					self.uncertainties[str(i)] = np.linspace(self.uncertainties[str(i)][0],self.uncertainties[str(i)][1],20)

					a = self.thermo_dict[str(i)]
					func = IFRT.function(str(i),a,self.temperatures[str(i)])
					#print(str(i),a,self.temperatures[str(i)])
					#print(func)
					y = self.uncertainties[str(i)]						
					#print(y*func)
					self.uncertainties[str(i)] = np.asarray((y*func)).flatten()
					continue
			
		
			
		self.doUnsrtAnalysis()
		for i in self.branches.split(","):
			self.cholskyDeCorrelateMat[str(i)],self.zeta[str(i)] = self.unsrt(str(i))
		self.nametag = self.species+":"+self.temp_limit	
			
		self.corelate_block = block_diag(*(self.cholskyDeCorrelateMat["Hcp"],self.cholskyDeCorrelateMat["h"],self.cholskyDeCorrelateMat["e"]))
		
		self.f = self.getUncertainty(self.temperatures)
		
		
		
		
	def getAllData(self):
		b1 = self.branches.split(",")[0]
		b2 = self.branches.split(",")[1]
		if len(self.branches.split(",")) >2:
			b3 = self.branches.split(",")[3]
		else:
			b3 = ""
		exp_data_type = self.exp_data_type["Hcp"].split(";")[0]
		exp_format = self.exp_data_type["Hcp"].split(";")[1]
		Log_string = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.classification,self.type,self.sub_type,self.branching,b1,b2,b3,self.pressure_limit,self.common_temp,self.temp_limit,exp_data_type,exp_format,self.nametag)
		exp_unsrt_string = ""
		solver_log = "######################\n{}######################\n\t\t{}\n\n".format(self.nametag,self.solution)
		calc_cholesky = self.corelate_block
		zeta_string = "{}\t{}\t{}\n".format(self.zeta["Hcp"],self.zeta["h"],self.zeta["e"])
		#print(self.temperatures)
		#print(self.uncertainties)
		for i in range(len(self.temperatures["Hcp"])):
			exp_unsrt_string += "{}\t{}\t{}\t{}\t{}\t{}\n".format(self.temperatures["Hcp"][i],self.uncertainties["Hcp"][i],self.temperatures["h"][i],self.uncertainties["h"][i],self.temperatures["e"][i],self.uncertainties["e"][i])
		string_2 = "temp\tHcp\ttemp\th\ttemp\te\t\n"
		file_unsrt = open("./Data/"+self.nametag+"_usrtData.log","w")
		file_unsrt.write(string_2+"\n"+exp_unsrt_string)
		file_unsrt.close()
		return Log_string,exp_unsrt_string,solver_log,calc_cholesky,zeta_string
	
	
	def getDtList(self):
		self.Unsrt_dict = {}
		key = ["Tag","Solution","Class","Type","Sub_type","Branch_boolen","Branches","Pressure_limit","Common_temp","temp_list","Nominal","Exp_input_data_type","priorCovariance","Basis_vector","Uncertainties","Temperatures","unsrt_func","Data_key"]
		values = [self.tag,self.solution,self.classification,self.type,self.sub_type,self.branching,self.branches,self.pressure_limit,self.common_temp,self.temp_limit,self.thermo_dict,self.exp_data_type,self.cholskyDeCorrelateMat,self.zeta,self.uncertainties,self.temperatures,self.f,"Hcp,h,e"]		
		
		for i,element in enumerate(key):
			self.Unsrt_dict[element] = values[i]
		return self.Unsrt_dict
		
	def doUnsrtAnalysis(self):
		T = self.temperatures[self.branches.split(",")[0]]
		guess = 0.001*np.ones(17)
		self.solution = minimize(self.uncertPriorObjective,guess,method="Nelder-Mead",options={'maxiter': 100000, 'maxfev': 100000, 'disp': False, 'return_all': False, 'initial_simplex': None, 'xatol': 1E-05, 'fatol': 1E-05, 'adaptive': True})
		Tscale = 5000
		#print(self.solution)
		L11 = self.solution.x[0]
		L21 = self.solution.x[1]/Tscale
		L22 = self.solution.x[2]/Tscale
		L31 = self.solution.x[3]/Tscale**2
		L32 = self.solution.x[4]/Tscale**2
		L33 = self.solution.x[5]/Tscale**2
		L41 = self.solution.x[6]/Tscale**3
		L42 = self.solution.x[7]/Tscale**3
		L43 = self.solution.x[8]/Tscale**3
		L44 = self.solution.x[9]/Tscale**3
		L51 = self.solution.x[10]/Tscale**4
		L52 = self.solution.x[11]/Tscale**4
		L53 = self.solution.x[12]/Tscale**4
		L54 = self.solution.x[13]/Tscale**4
		L55 = self.solution.x[14]/Tscale**4
		L66 = self.solution.x[15]
		L77 = self.solution.x[16]
		self.Lcp = np.array([[L11,L21,L31,L41,L51],[0,L22,L32,L42,L52],[0,0,L33,L43,L53],[0,0,0,L44,L54],[0,0,0,0,L55]])
		self.LH = np.array([L66])
		self.LS = np.array([L77])
		
		self.cholskyDeCorrelateMat_cp = np.matrix(self.Lcp.T)
		self.cholskyDeCorrelateMat_H  = np.matrix(self.LH.T)
		self.cholskyDeCorrelateMat_S  = np.matrix(self.LS.T)
		
		
		theta_cp = np.array([T/T,T,T**2,T**3,T**4])
		theta_H = np.array([T/T])
		theta_S = np.array([T/T])
		
		#Find zeta values
		guess_zeta = 0.01*np.array([1,1,1,1,1,1,1])
		self.zeta = minimize(self.obj_zeta,guess_zeta)
		self.zeta_cp = self.zeta.x[0:5]
		self.zeta_h = self.zeta.x[5]
		self.zeta_s = self.zeta.x[6]
	
	def unsrt(self,index):
		if index == "Hcp":
			L = self.Lcp
			z = self.zeta_cp
		if index == "h":
			L = self.LH
			z = self.zeta_h
		if index == "e":
			L = self.LS
			z = self.zeta_s
		return L,z
		
		
	def func_4(self,T,L11,L12,L22,L13,L23,L33,L14,L24,L34,L44,L15,L25,L35,L45,L55):
		unsrt = 9*((L55*T**4)**2+(L44*T**3+L55*T**4)**2+(L33*T**2+L34*T**3+L35*T**4)**2+(L22*T+L23*T**2+L24*T**3+L25*T**4)**2+(L11+L12*T+L13*T**2+L14*T**3+L15*T**4)**2)
		return unsrt
	def func_zeta(self,T,L0,L1,L2,L3,L4):
		z = np.array([L0,L1,L2,L3,L4])
		L11 = self.solution[0]
		L12 = self.solution[1]
		L22 = self.solution[2]
		L13 = self.solution[3]
		L23 = self.solution[4]
		L33 = self.solution[5]
		L14 = self.solution[6]
		L24 = self.solution[7]
		L34 = self.solution[8]
		L44 = self.solution[9]
		L15 = self.solution[10]
		L25 = self.solution[11]
		L35 = self.solution[12]
		L45 = self.solution[13]
		L55 = self.solution[14]
		fdiff = ((z[0]*L11)+(z[0]*L12+z[1]*L22)*T+(z[0]*L13+z[1]*L23+z[2]*L33)*T**2+(z[0]*L14+z[1]*L24+z[2]*L34+z[3]*L44)*(T**3)+(L15*z[0]+L25*z[1]+L35*z[2]+L45*z[3]+L55*z[4])*T**4)	
		return fdiff

	def getUncertainty(self,T): 
		unsrt = {}
		Tscale = 5000
		T = T["Hcp"]
		L11 = self.solution.x[0]
		L21 = self.solution.x[1]/Tscale
		L22 = self.solution.x[2]/Tscale
		L31 = self.solution.x[3]/Tscale**2
		L32 = self.solution.x[4]/Tscale**2
		L33 = self.solution.x[5]/Tscale**2
		L41 = self.solution.x[6]/Tscale**3
		L42 = self.solution.x[7]/Tscale**3
		L43 = self.solution.x[8]/Tscale**3
		L44 = self.solution.x[9]/Tscale**3
		L51 = self.solution.x[10]/Tscale**4
		L52 = self.solution.x[11]/Tscale**4
		L53 = self.solution.x[12]/Tscale**4
		L54 = self.solution.x[13]/Tscale**4
		L55 = self.solution.x[14]/Tscale**4
		L66 = self.solution.x[15]
		L77 = self.solution.x[16]
		# Truncation at 3 sigma
		unsrt_cp = 3*np.sqrt((L55*T**4)**2+(L44*T**3+L55*T**4)**2+(L33*T**2+L43*T**3+L53*T**4)**2+(L22*T+L32*T**2+L42*T**3+L52*T**4)**2+(L11+L21*T+L31*T**2+L41*T**3+L51*T**4)**2)
		unsrt_H = 3*np.sqrt((L66*(T/T))**2)
		unsrt_S = 3*np.sqrt((L77*(T/T))**2)
		unsrt["Hcp"] = unsrt_cp
		unsrt["h"] = unsrt_H
		unsrt["e"] = unsrt_S
		return unsrt
	
	def uncertPriorObjective (self,guess):
		Tscale = 5000
		R = 8.314
		Z= guess
		#z = np.array([Z[15],Z[16],Z[17],Z[18],Z[19]])
		L11 =Z[0]
		L21 =Z[1]
		L22 =Z[2]
		L31 =Z[3]
		L32 =Z[4]
		L33 =Z[5]
		L41 =Z[6]
		L42 =Z[7]
		L43 =Z[8]
		L44 =Z[9]
		L51 =Z[10]
		L52 =Z[11]
		L53 =Z[12]
		L54 =Z[13]
		L55 =Z[14]
		L66=Z[15]
		L77=Z[16]
		Lcp = np.array([[L11,L21,L31,L41,L51],[0,L22,L32,L42,L52],[0,0,L33,L43,L53],[0,0,0,L44,L54],[0,0,0,0,L55]])
		Lh = np.array([L66])
		Ls = np.array([L77])

		if "h" in self.sub_type:
			Y_h = self.uncertainties["h"]
			thetaH = np.array([T/T])
			sigma_H = 9*(np.dot(Lh,thetaH))**2

		if "e" in self.sub_type:
			Y_s = self.uncertainties["e"]
			thetaS = np.array([T/T])
			sigma_S = 9*(np.dot(Ls,thetaS))**2

		if "Hcp" in self.sub_type:
			Y_cp = self.uncertainties["Hcp"]
			thetaCP = np.array([T/T,T,T**2,T**3,T**4])
			unsrt = 9*((L55*T**4)**2+(L44*T**3+L55*T**4)**2+(L33*T**2+L43*T**3+L53*T**4)**2+(L22*T+L32*T**2+L42*T**3+L52*T**4)**2+(L11+L21*T+L31*T**2+L41*T**3+L51*T**4)**2)

		T = self.temperatures[self.branches.split(",")[0]]/Tscale


		unsrt = 9*((L55*T**4)**2+(L44*T**3+L55*T**4)**2+(L33*T**2+L43*T**3+L53*T**4)**2+(L22*T+L32*T**2+L42*T**3+L52*T**4)**2+(L11+L21*T+L31*T**2+L41*T**3+L51*T**4)**2)


		if "Hcp" in self.sub_type:
			residual_cp = (Y_cp-np.sqrt(unsrt))/(Y_cp/3)
		else:
			residual_cp = np.array([0])
		if "h" in self.sub_type:
			residual_h =(Y_h-np.sqrt(sigma_H))/(Y_h/3)
		else:
			residual_h =  np.array([0])
		if "e" in self.sub_type:
			residual_s =(Y_s-np.sqrt(sigma_S))/(Y_s/3)
		else:
			residual_s = np.array([0])
	
		obj = np.dot(residual_cp.T,residual_cp)+np.dot(residual_h.T,residual_h)+np.dot(residual_s.T,residual_s)
		return obj
	
	def obj_zeta(self,guess):
		T = self.temperatures[self.branches.split(",")[0]]
		#T = np.linspace(1000,5000,100)
		Tscale =5000
		z = np.ones(7)
		L11 = self.solution.x[0]
		L21 = self.solution.x[1]/Tscale
		L22 = self.solution.x[2]/Tscale
		L31 = self.solution.x[3]/Tscale**2
		L32 = self.solution.x[4]/Tscale**2
		L33 = self.solution.x[5]/Tscale**2
		L41 = self.solution.x[6]/Tscale**3
		L42 = self.solution.x[7]/Tscale**3
		L43 = self.solution.x[8]/Tscale**3
		L44 = self.solution.x[9]/Tscale**3
		L51 = self.solution.x[10]/Tscale**4
		L52 = self.solution.x[11]/Tscale**4
		L53 = self.solution.x[12]/Tscale**4
		L54 = self.solution.x[13]/Tscale**4
		L55 = self.solution.x[14]/Tscale**4
		L66 = self.solution.x[15]
		L77 = self.solution.x[16]
		z[0] = guess[0]
		z[1] = guess[1]
		z[2] = guess[2]
		z[3] = guess[3]
		z[4] = guess[4]
		z[5] = guess[5]
		z[6] = guess[6]
		
		zetaFunc_cp = ((z[0]*L11)*(T/T)+(z[0]*L21+z[1]*L22)*T+(z[0]*L31+z[1]*L32+z[2]*L33)*T**2+(z[0]*L41+z[1]*L42+z[2]*L43+z[3]*L44)*(T**3)+(L51*z[0]+L52*z[1]+L53*z[2]+L54*z[3]+L55*z[4])*T**4)
		zetaFunc_H = (T/T)*z[5]*L66
		zetaFunc_S = (T/T)*z[6]*L77
		
		
		if "Hcp" in self.sub_type:
			residual_cp =self.uncertainties["Hcp"]-zetaFunc_cp
		else:
			residual_cp = 0
		if "h" in self.sub_type:
			residual_H = self.H_uncertainties["h"]-zetaFunc_H
		else:
			residual_H = 0
		if "e" in self.sub_type:
			residual_S = self.S_uncertainties["e"]-zetaFunc_S
		else:
			residual_S = 0	
		
		obj = np.dot(residual_cp,residual_cp)+np.dot(residual_H,residual_H)+np.dot(residual_S,residual_S)
		return obj
		
#############################################################
###	   Uncertainty analysis for fallOffCurves,	  ######
###	   enthalpy, entropy and transport			  ######
###	   properties for kinetic species			   ######
#############################################################		
#Temperature independent uncetainties
#zeta = 3
#sigma_p = L11
	
class transport:
	def __init__(self, Element,transport_loc):
		self.species = self.classification = self.type = self.sub_type  = self. branching = self.branches  = self.pressure_limit  = self. common_temp = self.temp_limit =  None
		self.tag = Element.tag
		self.exp_data_type = {}
		self.temperatures = {}
		self.uncertainties = {}
		self.cholskyDeCorrelateMat = {}
		self.zeta = {} 
		self.species = Element.attrib["species"]
		
		self.trans_dict = IFR.TransportParsing(transport_loc).getTransportData(self.species)
		
		for item in Element:
			#print(item.tag)
			if item.tag == "class":
				self.classification = item.text
				continue
			if item.tag == "type":
				self.type = item.text
				continue
			if item.tag == "sub_type":
				self.sub_type = item.attrib["name"]
				
				for subitem in item:
					if "multiple" in subitem.tag:
						self.branching = str(subitem.text)
						continue
					if "branches" in subitem.tag:
						self.branches = str(subitem.text)
						continue
						
					if "pressure_limit" in subitem.tag:
						self.pressure_limit = str(subitem.text)
						continue
					if "common_temp" in subitem.tag:
						self.common_temp = str(subitem.text)
						continue
					if  "temp_limit" in subitem.tag :
						self.temp_limit = str(subitem.text)
						
						continue
				continue
					
		for i in self.branches.split(","):
			
			for item in Element:
				if item.tag == "data_type_"+str(i):

					self.exp_data_type[str(i)] = item.text

					continue
					
				if item.tag == "temp_"+str(i):
					#print(item.text)
					self.temperatures[str(i)] = np.asarray([float(j) for j in item.text.split(",")])
					continue
						
				if item.tag == "unsrt_"+str(i):
					self.uncertainties[str(i)] = np.asarray([float(i) for i in item.text.split(",")])
					continue
		
		
		for i in self.exp_data_type:
			
			if self.exp_data_type[str(i)].split(";")[0] ==  "constant":
				
				if self.exp_data_type[str(i)].split(";")[1] == "array":
				
					continue
				
				elif self.exp_data_type[str(i)].split(";")[1] == "end_points":
				
					self.temperatures[str(i)] = np.linspace(self.temperatures[str(i)][0],self.temperatures[str(i)][1],20)
				
					self.uncertainties[str(i)] = np.linspace(self.uncertainties[str(i)][0],self.uncertainties[str(i)][1],20)
				
					continue

			elif self.exp_data_type[str(i)].split(";")[0] == "percentage":
				
				if self.exp_data_type[str(i)].split(";")[1] == "array":
				
					a = self.trans_dict[str(i)]
					y = self.uncertainties[str(i)]						
					self.uncertainties[str(i)] = y*a
					
					continue

				elif self.exp_data_type[str(i)].split(";")[1] == "end_points":

					self.temperatures[str(i)] = np.linspace(self.temperatures[str(i)][0],self.temperatures[str(i)][1],20)
					self.uncertainties[str(i)] = np.linspace(self.uncertainties[str(i)][0],self.uncertainties[str(i)][1],20)

					a = self.trans_dict[str(i)]
					y = self.uncertainties[str(i)]						
					self.uncertainties[str(i)] = y*a
					continue
				
			
			
		
		for i in self.branches.split(","):
			self.cholskyDeCorrelateMat[str(i)],self.zeta[str(i)] = self.unsrt(str(i))
			
			
		self.nametag = self.species	
		self.f = self.getUncertainty(self.temperatures)
		'''
		fig = plt.figure()
		plt.xlabel('Temperature (K)')
		plt.ylabel('Uncertainity ($\sigma$)')
		plt.title( string+'{}'.format(self.name), fontsize = 10)		
		plt.plot(self.temperatures,self.uncertainties,'o',label='exp uncertainties');
		plt.ylim(0,2*max(self.uncertainties[0],self.uncertainties[-1]))	
		plt.legend()
		my_path = os.getcwd()
		plt.savefig(my_path+'/Plots/TDeptUnsrt/'+self.name+'.png')
		
		'''	
		self.solution = self.cholskyDeCorrelateMat
		self.corelate_block = block_diag(*(self.cholskyDeCorrelateMat["LJe"],self.cholskyDeCorrelateMat["LJs"]))
	
	
	
	def getDtList(self):
		self.Unsrt_dict = {}
		key = ["Tag","Solution","Class","Type","Sub_type","Branch_boolen","Branches","Pressure_limit","Common_temp","temp_list","Nominal","Exp_input_data_type","priorCovariance","Basis_vector","Uncertainties","Temperatures","unsrt_func","Data_key"]
		values = [self.tag,self.solution,self.classification,self.type,self.sub_type,self.branching,self.branches,self.pressure_limit,self.common_temp,self.temp_limit,self.trans_dict,self.exp_data_type,self.cholskyDeCorrelateMat,self.zeta,self.uncertainties,self.temperatures,self.f,"LJe,LJs"]
		
		for i,element in enumerate(key):
			self.Unsrt_dict[element] = values[i]
		return self.Unsrt_dict
		
	def getAllData(self):
		b1 = self.branches.split(",")[0]
		b2 = self.branches.split(",")[1]
		if len(self.branches.split(",")) >2:
			b3 = self.branches.split(",")[3]
		else:
			b3 = ""
		exp_data_type = self.exp_data_type["LJe"].split(";")[0]
		exp_format = self.exp_data_type["LJe"].split(";")[1]
		Log_string = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.classification,self.type,self.sub_type,self.branching,b1,b2,b3,self.pressure_limit,self.common_temp,self.temp_limit,exp_data_type,exp_format,self.nametag)
		exp_unsrt_string = ""
		solver_log = "######################\n{}######################\n\t\t{}\n\n".format(self.nametag,self.solution)
		calc_cholesky = self.corelate_block
		zeta_string = "{}\t{}\n".format(self.zeta["LJe"],self.zeta["LJs"])
		#print(self.temperatures)
		#print(self.uncertainties)
		for i in range(len(self.temperatures["LJe"])):
			exp_unsrt_string += "{}\t{}\t{}\t{}\n".format(self.temperatures["LJe"][i],self.uncertainties["LJe"][i],self.temperatures["LJs"][i],self.uncertainties["LJs"][i])
		string_2 = "temp\tLJe\ttemp\tLJs\n"
		file_unsrt = open("./Data/"+self.nametag+"_usrtData.log","w")
		file_unsrt.write(string_2+"\n"+exp_unsrt_string)
		file_unsrt.close()
		return Log_string,exp_unsrt_string,solver_log,calc_cholesky,zeta_string
	
	
	def getUncertainty(self,T): 
		unsrt = {}
		T1 = T["LJe"]
		T2 = T["LJs"]
		unsrt_LJe = np.sqrt((float(self.uncertainties["LJe"][0])*(T1/T1))**2)
		unsrt_LJs = np.sqrt((float(self.uncertainties["LJs"][0])*(T2/T2))**2)
		unsrt["LJe"] = unsrt_LJe
		unsrt["LJs"] = unsrt_LJs
		return unsrt
	
	def unsrt(self,index):
		if index == "LJe":
			L = self.uncertainties[index][0]
			z = 1
		if index == "LJs":
			L = self.uncertainties[index][0]
			z = 1
		#print(L)
		return L,z		
	

class collision:
	def __init__(self, Element,string,mechanism_loc):
		self.rxn = self.classification = self.type = self.sub_type = self. branching = self.branches  = self.pressure_limit  = self. common_temp = self.temp_limit= None
		
		self.tag = Element.tag
		self.cholskyDeCorrelateMat = {}
		self.zeta = {}
		self.exp_data_type = {}
		self.temperatures = {} 
		self.uncertainties = {}
		self.nominal = {}
		self.rxn = Element.attrib["rxn"]
		
		
		for item in Element:
			#print(item.tag)
			if item.tag == "class":
				self.classification = item.text
				continue
			if item.tag == "type":
				self.type = item.text
				continue
			if item.tag == "sub_type":
				self.sub_type = item.attrib["name"]
				
				for subitem in item:
					if "multiple" in subitem.tag:
						self.branching = str(subitem.text)
						continue
					if "branches" in subitem.tag:
						self.branches = str(subitem.text)
						self.m_dict = IFR.MechParsing(mechanism_loc).getThirdBodyCollisionEff(self.rxn,self.branches)
						continue
						
					if "pressure_limit" in subitem.tag:
						self.pressure_limit = str(subitem.text)
						continue
					if "common_temp" in subitem.tag:
						self.common_temp = str(subitem.text)
						continue
					if  "temp_limit" in subitem.tag :
						self.temp_limit = str(subitem.text)
						
						continue
				continue
					
		for i in self.branches.split(","):
			
			for item in Element:
				if item.tag == "data_type_"+str(i):

					self.exp_data_type[str(i)] = item.text

					continue
					
				if item.tag == "temp_"+str(i):
					#print(item.text)
					self.temperatures[str(i)] = np.asarray([float(j) for j in item.text.split(",")])
					continue
						
				if item.tag == "unsrt_"+str(i):
					self.uncertainties[str(i)] = np.asarray([float(i) for i in item.text.split(",")])
					continue
		
		
		for i in self.exp_data_type:
			
			if self.exp_data_type[str(i)].split(";")[0] ==  "constant":
				
				if self.exp_data_type[str(i)].split(";")[1] == "array":
				
					continue
				
				elif self.exp_data_type[str(i)].split(";")[1] == "end_points":
				
					self.temperatures[str(i)] = np.linspace(self.temperatures[str(i)][0],self.temperatures[str(i)][1],20)
				
					self.uncertainties[str(i)] = np.linspace(self.uncertainties[str(i)][0],self.uncertainties[str(i)][1],20)
				
					continue

			elif self.exp_data_type[str(i)].split(";")[0] == "percentage":
				
				if self.exp_data_type[str(i)].split(";")[1] == "array":
				
					a = self.trans_dict[str(i)]
					y = self.uncertainties[str(i)]						
					self.uncertainties[str(i)] = y*a
					
					continue

				elif self.exp_data_type[str(i)].split(";")[1] == "end_points":

					self.temperatures[str(i)] = np.linspace(self.temperatures[str(i)][0],self.temperatures[str(i)][1],20)
					self.uncertainties[str(i)] = np.linspace(self.uncertainties[str(i)][0],self.uncertainties[str(i)][1],20)

					a = self.trans_dict[str(i)]
					y = self.uncertainties[str(i)]						
					self.uncertainties[str(i)] = y*a
					continue
				
			
			
		
		for i in self.branches.split(","):
			self.cholskyDeCorrelateMat[str(i)],self.zeta[str(i)] = self.unsrt(self.uncertainties[str(i)])
	
		
		self.nametag = self.rxn+":"+self.sub_type
		self.solution = self.cholskyDeCorrelateMat
		self.f = self.getUncertainty()
		#print(self.m_dict)
		'''
		fig = plt.figure()
		plt.xlabel('Temperature (K)')
		plt.ylabel('Uncertainity ($\sigma$)')
		plt.title( string+'{}'.format(self.name), fontsize = 10)		
		plt.plot(self.temperatures,self.uncertainties,'o',label='exp uncertainties');
		plt.ylim(0,2*max(self.uncertainties[0],self.uncertainties[-1]))	
		plt.legend()
		my_path = os.getcwd()
		plt.savefig(my_path+'/Plots/TDeptUnsrt/'+self.name+'.png')
		
		'''	
	def getDtList(self):
		self.Unsrt_dict = {}
		key = ["Tag","Solution","Class","Type","Sub_type","Branch_boolen","Branches","Pressure_limit","Common_temp","temp_list","Nominal","Exp_input_data_type","priorCovariance","Basis_vector","Uncertainties","Temperatures","unsrt_func","Data_key"]
		values = [self.tag,self.solution,self.classification,self.type,self.sub_type,self.branching,self.branches,self.pressure_limit,self.common_temp,self.temp_limit,self.m_dict,self.exp_data_type,self.cholskyDeCorrelateMat,self.zeta,self.uncertainties,self.temperatures,self.f,self.branches]
		for i,element in enumerate(key):
			self.Unsrt_dict[element] = values[i]
		return self.Unsrt_dict
	
	def getAllData(self):
		Log_string = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.classification,self.type,self.sub_type,self.branching,self.branches,self.pressure_limit,self.common_temp,self.temp_limit,self.exp_data_type,self.nametag)
		exp_unsrt_string = ""
		solver_log = "{}\n{}\n".format(self.nametag,self.sol)
		calc_cholesky = self.cholskyDeCorrelateMat
		zeta_string = "{}".format(self.zeta)
		for i in range(len(self.temperatures)):
			exp_unsrt_string += "{}\t{}\n".format(self.temperatures[i],self.uncertainties[i])
		return Log_string,exp_unsrt_string,solver_log,calc_cholesky,zeta_string
		
	def getUncertainty(self): 
		unsrt = {}
		
		for i in self.branches.split(","):
			#print(i)
			T = self.temperatures[str(i)]
			unsrt_i = np.sqrt((float(self.uncertainties[str(i)][0])*(T/T))**2)
			unsrt[str(i)] = unsrt_i
		return unsrt
	def unsrt(self,unsrt):
		L = np.array([unsrt[0]])
		zeta = 1
		return L,zeta
	
#Parent class to host the uncertainty data for all the reactions. Has instances of reaction class as members		
class uncertaintyData:
	def __init__(self,pathDictionary,binary_files,unsrt_type=None):
	
		self.xmlPath = pathDictionary["uncertainty_data"]
		self.mechPath = pathDictionary["mechanism"]
		self.thermoPath = pathDictionary["thermo_file"]
		self.transportPath = pathDictionary["trans_file"]
		
		self.Reactionlist = []
		self.interConnectedlist = []
		self.PlogRxnlist = []
		self.focList = []
		self.mList = []
		self.thermoList = []
		self.transportList = []
		self.PlogRxnIndex = {}
		self.unsrt_data = {}
		self.reactionUnsrt = {}
		self.focUnsrt = {}
		self.plogUnsrt = {}
		self.Plog_index = {}
		self.trdBodyUnsrt = {}
		self.interconnectedRxn = {}
		#self.Plog_additional_Rxn = {}
		self.thermoUnsrt = {}
		self.transportUnsrt = {}
		self.tree = ET.parse(self.xmlPath)
		self.root = self.tree.getroot()
		count_rxn = 0
		count_foc = 0
		count_m = 0
		count_thermo = 0
		count_transport = 0
		
		for child in self.root:
			if child.tag == 'reaction':
				r = reaction(child,self.mechPath,binary_files)
				self.Reactionlist.append(r.nametag)
				self.reactionUnsrt[r.nametag] = r
				self.unsrt_data[r.nametag] = r
				count_rxn +=1
			if child.tag == 'PLOG':
				p = PLOG(child,self.mechPath,binary_files)
				self.PlogRxnlist.append(p.rIndex)
				self.PlogRxnIndex[p.rIndex] = p.nametag
				self.Plog_index[p.rIndex] = p.index
				self.plogUnsrt[p.nametag] = p
				self.unsrt_data[p.nametag] = p
				self.Reactionlist.append(p.nametag)
				count_rxn +=1
			
			if child.tag == "PLOG-Interconnectedness":
				plog_object_list = []
				list_of_PLOG_rxn = None
				#print(self.PlogRxnIndex)
				#print(self.Plog_index)
				for item in child:
					if item.tag == "InterConnectedRxns":				
						list_of_PLOG_rxn = item.text.split(",")
						
					if item.tag == "RxnCount":
						count = int(item.text)
				#print(list_of_PLOG_rxn)
				for i in list_of_PLOG_rxn:
					if i in self.PlogRxnIndex:
						#print(i)
						plog_object_list.append(self.plogUnsrt[self.PlogRxnIndex[i]])
					else:
						raise AssertionError(f"Invalid connected reactions are identified. Please check {item.text}")
				for i in range(count):
					index = i+1
					
					q = PLOG_Interconnectedness(plog_object_list,index,count)
					self.interconnectedRxn[q.nametag] = q
					self.interConnectedlist.append(q.nametag)
					self.unsrt_data[q.nametag] = q
					count_rxn +=1
					
			if child.tag == 'fallOffCurve':
				foc = fallOffCurve(child,self.mechPath)
				self.focList.append(foc.nametag)
				self.focUnsrt[foc.nametag] = foc
				self.unsrt_data[foc.nametag] = foc
				count_foc +=1
			if child.tag == 'thermo':
				th = thermodynamic(child,self.thermoPath)
				self.thermoList.append(th.nametag)
				self.thermoUnsrt[th.nametag] = th
				self.unsrt_data[th.nametag] = th
				count_thermo +=1
			if child.tag == 'collisionEff':
				string = "Unsrt for third bodies\n collision efficiencies [M]:  "
				m = collision(child,string,self.mechPath)
				self.mList.append(m.nametag)
				self.trdBodyUnsrt[m.nametag] = m
				self.unsrt_data[m.nametag] = m
				count_m +=1
			if child.tag == 'transport':
				tr = transport(child,self.transportPath)
				self.transportUnsrt[tr.nametag] = tr
				self.transportList.append(tr.nametag)
				self.unsrt_data[tr.nametag] = tr
				count_transport +=1
				#print(self.transportList)
		if unsrt_type !="opt":
			print("\n\n{} Reactions are selected for optimization\n".format(count_rxn))
			print("{} Fall-off (center broadening factors) of reactions {} are selected for optimization\n\n".format(count_foc,self.focList))
			print("{} third body collision efficiency's are selected for optimization\n\n".format(count_m))
			print("{} thermo-chemical parameters are selected for optimization\n\n".format(count_thermo))
			print("{} transport parameters are selected for optimization\n\n".format(count_transport))
			
	def extract_uncertainty(self):
		#print(self.root)
		return self.unsrt_data,self.reactionUnsrt,self.plogUnsrt,self.interconnectedRxn, self.focUnsrt, self.trdBodyUnsrt, self.thermoUnsrt, self.transportUnsrt, self.Reactionlist,self.PlogRxnlist,self.interConnectedlist,self.focList,self.mList,self.thermoList,self.transportList
	
		
