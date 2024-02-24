import os
import numpy as np
from copy import deepcopy
###Tasks
# Take a copy of the mechanism in yaml format and manipulate it
# After manipulating, convert it back to ck format
####################################################3
###

class Manipulator:
	def __init__(self,copy_of_mech,unsrt_object,perturbation):
		
		self.mechanism = copy_of_mech
		self.unsrt = unsrt_object
		self.perturbation = perturbation
		self.rxn_list = [rxn for rxn in self.unsrt]
		
	def getRxnPerturbationDict(self):
		perturb = {}
		count = 0
		for rxn in self.rxn_list:
			len_active_parameters = len(self.unsrt[rxn].activeParameters)
			temp = []			
			for i in range(len_active_parameters):
				temp.append(self.perturbation[count])
				count+=1
			perturb[rxn] = np.asarray(temp)		
		return perturb
	
	
	def getRxnType(self):
		rxn_type = {}
		for rxn in self.rxn_list:
			rxn_type[rxn] = self.unsrt[rxn].classification
		return rxn_type
		
	def ElementaryPerturbation(self,rxn,beta,mechanism):
		
		"""
		This function perturbes the Elementary reactions:
			Inputs:
				1) Rxn object
				2) Rxn perturbation factor (beta)
					- For A-factor type perturbation the L (choleskyDeCorrelateMat) is a factor
					- For all Arrhenius parameter perturbation, L is a matrix storing the uncertainty data
			Output:
				perturbed reaction using the following formula
				
				$p = p_0 + L\zeta$		
		
		"""
		
		
		index = self.unsrt[rxn].index# To search the reaction from Yaml dictionary
		cov = self.unsrt[rxn].cholskyDeCorrelateMat #The L matrix storing the uncertainty data
		zeta = self.unsrt[rxn].solution
		perturbation_factor = self.unsrt[rxn].perturb_factor
		p0 = self.unsrt[rxn].nominal
		unsrt_perturbation = np.asarray(cov.dot(beta)).flatten()
		convertor = np.asarray(self.unsrt[rxn].selection)
		p = p0+convertor*unsrt_perturbation
		
		"""
		The new perturbed reaction replaces the prior Arrhenius parameters 
		"""
		reaction_details = mechanism["reactions"][index]["rate-constant"]
		reaction_details["A"] = float(np.exp(p[0]))
		reaction_details["b"] = float(p[1])
		reaction_details["Ea"] = float(p[2]*1.987)
		mechanism["reactions"][index]["rate-constant"] = deepcopy(reaction_details)
		
		return mechanism
				
	def PlogPerturbation(self,rxn,beta,mechanism):
		if "PLOG" in rxn:
			rxn_split_index = int(rxn.split(":")[1].split("_")[1])
			
		index = self.unsrt[rxn].index
		cov = self.unsrt[rxn].cholskyDeCorrelateMat
		zeta = self.unsrt[rxn].solution
		p0 = self.unsrt[rxn].nominal
		perturbation_factor = self.unsrt[rxn].perturb_factor
		#selection = np.asarray(self.unsrt[rxn].selection)
		unsrt_perturbation = np.asarray(cov.dot(beta.dot(perturbation_factor))).flatten()
		convertor = np.asarray(self.unsrt[rxn].selection)
		
		p = p0+convertor*unsrt_perturbation
		pressure_limit = self.unsrt[rxn].pressure_limit
		if pressure_limit == "High":
			pressure_index = -1
		elif pressure_limit == "Low":
			pressure_index = 0
		else:
			pressure_index = rxn_split_index
		
		reaction_details = mechanism["reactions"][index]["rate-constants"][pressure_index]
		#print(f"==================\n\tBefore\t{reaction_details}\n")
		reaction_details["A"] = float(np.exp(p[0]))
		reaction_details["b"] = float(p[1])
		reaction_details["Ea"] = float(p[2]*1.987)
		#print(reaction_details)
		mechanism["reactions"][index]["rate-constants"][pressure_index] = deepcopy(reaction_details)
		after = mechanism["reactions"][index]["rate-constants"][pressure_index]
		#print(f"\tAfter\t{after}\n")
		return mechanism
		
	def BranchingReactions(self,rxn,beta,mechanism):
		indexes = []
		rxn_index = self.unsrt[rxn].index
		indexes.append(rxn_index)
		branches = self.unsrt[rxn].branches
		indexes.extend(branches)
		cov = self.unsrt[rxn].cholskyDeCorrelateMat
		zeta = self.unsrt[rxn].solution
		perturbation_factor = self.unsrt[rxn].perturb_factor
		#p0 = self.unsrt[rxn].nominal
		#unsrt_perturbation = np.asarray(cov.dot(selection*beta.dot(perturbation_factor))).flatten()
		convertor = np.asarray(self.unsrt[rxn].selection)
		#p = p0+convertor*unsrt_perturbation
		#print("==================================\nPerturbation Started!!\n")
		#for index in indexes:
			#reaction_details = mechanism["reactions"][index]["rate-constant"]
			#print(reaction_details)
			#p0 = np.array([reaction_details["A"],reaction_details["b"],reaction_details["Ea"]])
			#print(p0)
		
		for index in indexes:
		#	print(index)
			reaction_details = mechanism["reactions"][index]["rate-constant"]
			p0 = np.asarray([float(np.log(reaction_details["A"])), float(reaction_details["b"]), float(reaction_details["Ea"]/1.987)])
		#	print(p0)
		#	print(type(p0))
			unsrt_perturbation = np.asarray(cov.dot(beta.dot(perturbation_factor))).flatten()
		#	print(unsrt_perturbation)
		#	print(convertor*unsrt_perturbation)
		#	print(type(unsrt_perturbation))
		#	print(type(convertor))
		#	print(type(convertor*unsrt_perturbation))
			p = p0+np.asarray(convertor*unsrt_perturbation)
		#	print(p)
		#	print(type(p))
			reaction_details["A"] = float(np.exp(p[0]))
			reaction_details["b"] = float(p[1])
			reaction_details["Ea"] = float(p[2]*1.987)
			mechanism["reactions"][index]["rate-constant"] = deepcopy(reaction_details)
		
		#for index in indexes:
		#	reaction_details = mechanism["reactions"][index]["rate-constant"]
			#p0 = np.asrray([reaction_details["A"],reaction_details["b"],reaction_details["Ea"]])
		#	print(reaction_details)
		
		#print("Perturbation Done!!")
		
		return mechanism

	def ThirdBodyPerturbation(self,rxn,beta,mechanism):
		pass
		
	def TroePerturbation(self,rxn,beta,mechanism):
		
		P_limit = self.unsrt[rxn].pressure_limit
		index = self.unsrt[rxn].index
		cov = self.unsrt[rxn].cholskyDeCorrelateMat
		zeta = self.unsrt[rxn].solution
		perturbation_factor = self.unsrt[rxn].perturb_factor
		p0 = self.unsrt[rxn].nominal
		unsrt_perturbation = np.asarray(cov.dot(beta)).flatten()
		convertor = np.asarray(self.unsrt[rxn].selection)
		
		p = p0+convertor*unsrt_perturbation
		
		if P_limit == "High":
			reaction_details = mechanism["reactions"][index]["high-P-rate-constant"]
		else:
			reaction_details = mechanism["reactions"][index]["low-P-rate-constant"]
		reaction_details["A"] = float(np.exp(p[0]))
		reaction_details["b"] = float(p[1])
		reaction_details["Ea"] = float(p[2]*1.987)
		if P_limit == "High":
			mechanism["reactions"][index]["low-P-rate-constant"] = deepcopy(reaction_details)
		else:
			mechanism["reactions"][index]["low-P-rate-constant"] = deepcopy(reaction_details)
		
		return mechanism
		
		
	def doPerturbation(self):
		rxn_type = self.getRxnType()
		perturb = self.getRxnPerturbationDict()
		mechanism = self.mechanism
		
		for rxn in self.rxn_list:
			beta = np.asarray(perturb[rxn])
			type_of_rxn = rxn_type[rxn]
			
			if type_of_rxn == "Elementary":
				index = self.unsrt[rxn].index
				new_mechanism = self.ElementaryPerturbation(rxn,beta,mechanism)
				mechanism = new_mechanism
			
			if type_of_rxn == "PLOG":
				index = self.unsrt[rxn].index
				new_mechanism = self.PlogPerturbation(rxn,beta,mechanism)
				mechanism = new_mechanism
				
			if type_of_rxn == "PLOG-Duplicate":
				index = self.unsrt[rxn].index
				new_mechanism = self.PlogPerturbation(rxn,beta,mechanism)
				mechanism = new_mechanism
				
			if type_of_rxn == "Duplicate":
				index = self.unsrt[rxn].index
				new_mechanism = self.ElementaryPerturbation(rxn,beta,mechanism)
				mechanism = new_mechanism
				
			if type_of_rxn == "Falloff":
				index = self.unsrt[rxn].index
				new_mechanism = self.TroePerturbation(rxn,beta,mechanism)
				mechanism = new_mechanism
				
			if type_of_rxn == "ThirdBody":
				index = self.unsrt[rxn].index
				new_mechanism = self.ElementaryPerturbation(rxn,beta,mechanism)
				mechanism = new_mechanism
				
			if type_of_rxn == "BranchingRxn":
				index = self.unsrt[rxn].index
				new_mechanism = self.BranchingReactions(rxn,beta,mechanism)
				mechanism = new_mechanism
		return mechanism,perturb
	
	
