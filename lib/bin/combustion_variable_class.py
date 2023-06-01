import numpy as np
import os
import Input_file_reader as read
from dataclasses import dataclass, field

class combustion_variable():
	def __init__(self,variable,variable_list):
		#Name:
		#Class:
		
		self.tag = self.name = self.solution = self.classification = self.type = self.subType = self.dataType = self.nominal = self.theta = self.priorUnsrt = self.posterior = self.dataRange = self.unsrtData = self.optimized = self.pressure_limit = self.common_temp = self.temp_limit  = self.branch_bool = self.branches = self.opt_branching_ratio = self.init_branching_ratio = self.prior_covariance = self.posterior_covariance = self.basisVector = self.dataKey = None
		
		self.sensitivity = {}
		self.name = variable
		self.tag = variable_list["unsrtDatabase"]["Tag"]
		self.solution = variable_list["unsrtDatabase"]["Solution"]
		self.classification = variable_list["unsrtDatabase"]["Class"]
		self.type = variable_list["unsrtDatabase"]["Type"]
		self.subType = variable_list["unsrtDatabase"]["Sub_type"]
		self.branch_bool = variable_list["unsrtDatabase"]["Branch_boolen"]
		self.branches = variable_list["unsrtDatabase"]["Branches"]
		self.pressure_limit = variable_list["unsrtDatabase"]["Pressure_limit"]
		self.common_temp = variable_list["unsrtDatabase"]["Common_temp"]
		self.temp_limit = variable_list["unsrtDatabase"]["temp_list"]
		self.dataType = variable_list["unsrtDatabase"]["Exp_input_data_type"]
		self.priorCovariance = variable_list["unsrtDatabase"]["priorCovariance"]
		self.nominal = variable_list["unsrtDatabase"]["Nominal"]
		self.basisVector= variable_list["unsrtDatabase"]["Basis_vector"]
		self.priorUnsrt = variable_list["unsrtDatabase"]["Uncertainties"]
		self.temperatures = variable_list["unsrtDatabase"]["Temperatures"]
		self.unsrt_function = variable_list["unsrtDatabase"]["unsrt_func"]
		self.dataKey = variable_list["unsrtDatabase"]["Data_key"]
		self.sa_Beta = variable_list["sensitivityManipulations"]	
		self.opt_Beta = variable_list["optManipulations"]
		self.originalBeta = variable_list["originalManipulations"]
		self.opt = variable_list["optimizedParameter"]
		self.posteriorCovariance = variable_list["posteriorCovariance"]
		
		if self.tag == "reaction":
			T = np.asarray(self.temperatures).flatten()
			self.theta = np.array([T/T,np.log(T),-1/T])
		if self.tag == "fallOffCurve":
			T = np.asarray(self.temperatures)
			self.theta = np.array([T/T])
		if self.tag == "thermo":
			self.theta = {}
			for i in self.dataKey.split(","):
				if i.strip() == "Hcp":
					T = np.asarray(self.temperatures[str(i)])
					self.theta[str(i)] = np.array([T/T,T,T**2,T**3,T**4])
				
				else:
					T = np.asarray(self.temperatures[str(i)])
					self.theta[str(i)] = np.array([T/T])
		if self.tag == "collisionEff":
			self.theta = {}
			for i in self.dataKey.split(","):
				T = np.asarray(self.temperatures[str(i)])
				self.theta[str(i)] = np.array([T/T])
		
		if self.tag == "transport":
			self.theta = {}
			for i in self.dataKey.split(","):
				T = np.asarray(self.temperatures[str(i)])
				self.theta[str(i)] = np.array([T/T])
		

		if self.dataKey.strip() == "":		
			#print(np.array([self.nominal]).flatten())
			#print(np.asarray(np.dot(self.priorCovariance,self.basisVector)).flatten())
			
			arrheniusShiftUp = np.asarray(self.nominal) + np.asarray(np.dot(self.priorCovariance,self.basisVector)).flatten()
			arrheniusShiftDown = np.asarray(self.nominal) - np.array(np.dot(self.priorCovariance,self.basisVector)).flatten()
		
			optShift = np.asarray(self.nominal) + np.array(np.dot(self.priorCovariance,np.asarray(self.opt)*self.basisVector)).flatten()
			
			self.upperBasisLimit = arrheniusShiftUp
			self.lowerBasisLimit = arrheniusShiftDown
			
			self.unsrtBasisFunction = np.asarray(np.dot(self.theta.T,np.array(np.dot(self.priorCovariance,self.basisVector)).flatten())).flatten()
			self.optPerturbation = optShift
			self.cases_sa = {}
			self.cases_opt = {}
			for k in self.sa_Beta:
				self.cases_sa[str(k)] = []
				self.cases_opt[str(k)] = []
				
				for index,m in enumerate(self.sa_Beta[str(k)]):
					
					betaShift = np.asarray(self.nominal) + np.array(np.dot(self.priorCovariance,np.asarray(m)*self.basisVector)).flatten()
					opt_betaShift = np.asarray(self.nominal) + np.array(np.dot(self.priorCovariance,np.asarray(self.opt_Beta[str(k)][index])*self.basisVector)).flatten()
					
					self.cases_sa[str(k)].append(betaShift)
					self.cases_opt[str(k)].append(opt_betaShift)
					
			self.sensitivity_perturbation = self.cases_sa
			self.optPRS_perturbation = self.cases_opt
			#self.postUpperUnsrtLimit = 
			#self.postLowerUnsrtLimit = 
		else:
				
			self.upperBasisLimit = {}
			self.lowerBasisLimit = {}
			self.unsrtBasisFunction = {}
			
			self.optPerturbation  = {}
			self.optPRS_perturbation = {}
			self.sensitivity_perturbation = {}
			for i in self.dataKey.split(","):
#				print(str(i))
#				print(self.tag)
#				print(np.asarray(self.nominal))
#				print(self.basisVector[str(i)])
#				print(self.priorCovariance[str(i)])
#				print(np.array(np.dot(self.priorCovariance[str(i)],self.basisVector[str(i)])).flatten())
				#print()
				arrheniusShiftUp = np.asarray(self.nominal[str(i)]) + np.array(np.dot(self.priorCovariance[str(i)],self.basisVector[str(i)])).flatten()
				arrheniusShiftDown = np.asarray(self.nominal[str(i)]) - np.array(np.dot(self.priorCovariance[str(i)],self.basisVector[str(i)])).flatten()
				
				#print(self.theta[str(i)].T)
				#print(np.dot(self.priorCovariance[str(i)],self.basisVector[str(i)]))
				sigmaZeta = np.asarray(np.dot(self.theta[str(i)].T,np.dot(self.priorCovariance[str(i)],self.basisVector[str(i)]))).flatten()
				#print(sigmaZeta)
				#sigmaZeta = np.array(np.dot(self.theta[str(i)].T,np.asarray(np.dot(self.priorCovariance[str(i)],self.basisVector[str(i)])).flatten())).flatten()
				
				#print(self.opt)
				optShift = np.asarray(self.nominal[str(i)]) + np.array(np.dot(self.priorCovariance[str(i)],np.asarray(self.opt[str(i)])*self.basisVector[str(i)])).flatten() 
				
				
				
				self.upperBasisLimit[str(i)] = arrheniusShiftUp
				self.lowerBasisLimit[str(i)] = arrheniusShiftDown
				self.unsrtBasisFunction[str(i)] = sigmaZeta
				self.optPerturbation[str(i)] = optShift
				
				self.cases_sa = {}
				self.cases_opt = {}
				for k in self.sa_Beta:
					self.cases_sa[str(k)] = []
					self.cases_opt[str(k)] = []
					for index,m in enumerate(self.sa_Beta[str(k)]):
						#print(m)
						betaShift_opt = np.asarray(self.nominal[str(i)]) + np.array(np.dot(self.priorCovariance[str(i)],np.asarray(m[str(i)])*self.basisVector[str(i)])).flatten() 
						betaShift = np.asarray(self.nominal[str(i)]) + np.array(np.dot(self.priorCovariance[str(i)],np.asarray(self.opt_Beta[str(k)][index][str(i)])*self.basisVector[str(i)])).flatten() 
						self.cases_sa[str(k)].append(betaShift)
				self.sensitivity_perturbation[str(i)] = self.cases_sa
				self.optPRS_perturbation[str(i)] = self.cases_opt
