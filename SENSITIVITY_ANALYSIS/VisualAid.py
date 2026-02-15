import numpy as np
import os
import matplotlib.pyplot as plt

class ArrheniusPlotter(object):
	def __init__(self,unsrt_object,reaction):
		self.unsrt_data = unsrt_object
		self.rxn = reaction
		self.M = 3.0/np.log(10.0)
	
	def getNominalParams(self):
		Nom = self.unsrt_data[self.rxn].nominal
		return Nom
	
	def getCholeskyCovariance(self):
		return self.unsrt_data[self.rxn].cov
	
	def getZetaMax(self):
		return self.unsrt_data[self.rxn].zeta.x

	def getTemperatures(self):
		return self.unsrt_data[self.rxn].temperatures
	
	def getTheta(self):
		T = self.getTemperatures()
		Theta = np.array([T/T,np.log(T),-1/T])
		return Theta																						
	
	def getUncertFunc(self):
		L = self.getCholeskyCovariance()
		Theta =   self.getTheta()
		func = [self.M*np.linalg.norm(np.dot(L.T,i)) for i in Theta.T]
		return np.asarray(func)
	
	def getZetaUnsrtFunc(self):
		L = self.getCholeskyCovariance()
		Theta =   self.getTheta()
		z = self.getZetaMax()
		func = [(i.T.dot(L.dot(z))) for i in Theta.T]
		return np.asarray(func)
	
	def getPerturbed_A_curve(self):
		L = self.getCholeskyCovariance()
		Theta =   self.getTheta()
		z = np.array([1,0,0])
		func = [(i.T.dot(L.dot(z))) for i in Theta.T]
		return np.asarray(func)
	
	def getPerturbed_n_curve(self):
		L = self.getCholeskyCovariance()
		Theta =   self.getTheta()
		z = np.array([0,1,0])
		func = [(i.T.dot(L.dot(z))) for i in Theta.T]
		return np.asarray(func)
	
	def getPerturbed_Ea_curve(self):
		L = self.getCholeskyCovariance()
		Theta =   self.getTheta()
		z = np.array([0,0,100])
		func = [(i.T.dot(L.dot(z))) for i in Theta.T]
		return np.asarray(func)
	
	def getNominalCurve(self):
		P = self.getNominalParams()
		Theta =   self.getTheta()
		func =  [(i.T.dot(P)) for i in Theta.T]
		return np.asarray(func)
	
	def plot_uncertainty_limits(self,location="Plots"):
		self.UQ_plot_loc = location
		os.makedirs(location,exist_ok = True)
		fig = plt.figure()
		T = self.getTemperatures()
		Kappa_o = self.getNominalCurve()
		Kappa_max = self.getZetaUnsrtFunc()
		UQ_limit = self.getUncertFunc()
		plt.plot(1/T,Kappa_o,"b-",label="Nominal Curve")
		plt.plot(1/T,Kappa_o + Kappa_max,"r-",label=r"Arrhenius Curve (f($\zeta$))")
		plt.plot(1/T,Kappa_o-Kappa_max,"r-")
		plt.plot(1/T,Kappa_o+UQ_limit,"k--",label=r"Uncertainty Limits")
		plt.plot(1/T,Kappa_o-UQ_limit,"k--")
		plt.xlabel("Temperatures (1/K)")
		plt.ylabel(r"Rate Coefficient $(\kappa)$")
		plt.legend()
		plt.savefig(location+f"/{self.rxn}.pdf",bbox_inches="tight")
	
	def plot_perturbed_Arrhenius_parameters(self,location="Plots"):
		self.UQ_plot_loc = location
		os.makedirs(location,exist_ok = True)
		fig = plt.figure()
		T = self.getTemperatures()
		Kappa_o = self.getNominalCurve()
		Kappa_max = self.getZetaUnsrtFunc()
		UQ_limit = self.getUncertFunc()
		Z_a = self.getPerturbed_A_curve()
		Z_n = self.getPerturbed_n_curve()
		Z_e = self.getPerturbed_Ea_curve()
		plt.plot(1/T,Kappa_o,"b-",label="Nominal Curve")
		plt.plot(1/T,Kappa_o + Kappa_max,"r-",label=r"Arrhenius Curve (f($\zeta$))")
		plt.plot(1/T,Kappa_o-Kappa_max,"r-")
		plt.plot(1/T,Kappa_o+UQ_limit,"k--",label=r"Uncertainty Limits")
		plt.plot(1/T,Kappa_o-UQ_limit,"k--")
		plt.plot(1/T,Kappa_o+Z_a,"b--",label="Perturbing A-factor")
		plt.plot(1/T,Kappa_o+Z_n,"c--",label="Perturbing n parameter")
		plt.plot(1/T,Kappa_o+Z_e,"y--",label="Perturbing Ea")
		plt.xlabel("Temperatures (1/K)")
		plt.ylabel(r"Rate Coefficient $(\kappa)$")
		plt.legend()
		plt.savefig(location+f"/{self.rxn}.pdf",bbox_inches="tight")		
