import os
import numpy as np
from copy import deepcopy
###Tasks
# Take a copy of the mechanism in yaml format and manipulate it
# After manipulating, convert it back to ck format
####################################################3
###

class Manipulator:
	def __init__(self,copy_of_mech,unsrt_object,perturbation,rxn_list,parameter_selection):
		self.mechanism = deepcopy(copy_of_mech)
		#self.initial_mech = copy_of_mech
		self.unsrt = unsrt_object
		self.perturbation = perturbation
		self.rxn_list = rxn_list
		self.selection = parameter_selection
		
	def getRxnPerturbationDict(self):
		perturb = {}
		selection = {}
		#print(len(self.selection))
		#print(self.selection)
		#print(len(self.perturbation))
		count = 0
		print(self.rxn_list)
		#print(len(self.rxn_list))
		for rxn in self.rxn_list:
			
			len_active_parameters = len(self.unsrt[rxn].activeParameters)
			#print(len_active_parameters)
			temp = []
			temp_1 = []
			#print(len_active_parameters)
			#print(count)
			for i in range(len_active_parameters):
				temp.append(self.perturbation[count])
				temp_1.append(self.selection[count])
				count+=1
			#print(count)
			perturb[rxn] = np.asarray(temp)
			selection[rxn] = np.asarray(temp_1)
		
		return perturb,selection	
	
	def getRxnType(self):
		rxn_type = {}
		#print(self.rxn_list)
		#for i in self.unsrt:
		#	print(i)
		for rxn in self.rxn_list:
			rxn_type[rxn] = self.unsrt[rxn].classification
		return rxn_type
		
	def ElementaryPerturbation(self,rxn,beta,selection,mechanism):
		index = self.unsrt[rxn].index
		cov = self.unsrt[rxn].cholskyDeCorrelateMat
		zeta = self.unsrt[rxn].solution
		perturbation_factor = self.unsrt[rxn].perturb_factor
		#perturbation_factor = np.array([1,1,1])#this is kept to see how the search iteration data is coming
		p0 = self.unsrt[rxn].nominal
		#print(p0)
		#unsrt_perturbation = np.asarray(cov.dot(selection*beta.dot(perturbation_factor))).flatten()
		unsrt_perturbation = np.asarray(cov.dot(beta)).flatten()
		convertor = np.asarray(self.unsrt[rxn].selection)
		#print(f"Perturbation factor = {perturbation_factor}")
		#print("\n")
		#print(beta)
		#print("\n")
		#print(f"Convertor = {convertor}")
		#print("\n")
		#print(unsrt_perturbation)
		p = p0+convertor*unsrt_perturbation
		#print(p)
		reaction_details = mechanism["reactions"][index]["rate-constant"]
		reaction_details["A"] = float(np.exp(p[0]))
		reaction_details["b"] = float(p[1])
		reaction_details["Ea"] = float(p[2]*1.987)
		mechanism["reactions"][index]["rate-constant"] = deepcopy(reaction_details)
		after = mechanism["reactions"][index]["rate-constant"]
		#print(f"{rxn}\t{cov}\t{p}\t{beta}\t{after}\n")
		return mechanism
		
	def DuplicatePerturbation(self,rxn,beta,selection,mechanism):
		index = self.unsrt[rxn].index
		cov = self.unsrt[rxn].cholskyDeCorrelateMat
		zeta = self.unsrt[rxn].solution
		perturbation_factor = self.unsrt[rxn].perturb_factor
		#perturbation_factor = np.array([1,1,1])
		p0 = self.unsrt[rxn].nominal
		#selection = np.asarray(self.unsrt[rxn].selection)
		#unsrt_perturbation = np.asarray(cov.dot(selection*beta.dot(perturbation_factor))).flatten()
		unsrt_perturbation = np.asarray(cov.dot(beta)).flatten()
		convertor = np.asarray(self.unsrt[rxn].selection)
		
		p = p0+convertor*unsrt_perturbation
		reaction_details = mechanism["reactions"][index]["rate-constant"]
		#print(reaction_details)
		#raise AssertionError("STOP!!")
		reaction_details["A"] = float(np.exp(p[0]))
		reaction_details["b"] = float(p[1])
		reaction_details["Ea"] = float(p[2]*1.987)
		mechanism["reactions"][index]["rate-constant"] = deepcopy(reaction_details)
		return mechanism
		
	def PlogPerturbation(self,rxn,beta,selection,mechanism):
		if "PLOG" in rxn:
			rxn_split_index = int(rxn.split(":")[1].split("_")[1])
			
		index = self.unsrt[rxn].index
		cov = self.unsrt[rxn].cholskyDeCorrelateMat
		zeta = self.unsrt[rxn].solution
		p0 = self.unsrt[rxn].nominal
		perturbation_factor = self.unsrt[rxn].perturb_factor
		#selection = np.asarray(self.unsrt[rxn].selection)
		unsrt_perturbation = np.asarray(cov.dot(selection*beta.dot(perturbation_factor))).flatten()
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
		
	def BranchingReactions(self,rxn,beta,selection,mechanism):
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
			unsrt_perturbation = np.asarray(cov.dot(selection*beta.dot(perturbation_factor))).flatten()
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
		
	def TroePerturbation(self,rxn,beta,selection,mechanism):
		P_limit = self.unsrt[rxn].pressure_limit
		index = self.unsrt[rxn].index
		cov = self.unsrt[rxn].cholskyDeCorrelateMat
		zeta = self.unsrt[rxn].solution
		perturbation_factor = self.unsrt[rxn].perturb_factor
		#perturbation_factor = np.array([1,1,1])
		p0 = self.unsrt[rxn].nominal
		#selection = np.asarray(self.unsrt[rxn].selection)
		#unsrt_perturbation = np.asarray(cov.dot(selection*beta.dot(perturbation_factor))).flatten()
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
		
		#print(mechanism["reactions"][index])
		#raise AssertionError("TROE check")
		return mechanism
		
		
	def doPerturbation(self):
		rxn_type = self.getRxnType()
		perturb,selection = self.getRxnPerturbationDict()
		mechanism = self.mechanism
		#print("\n==========Before====\n")
		#print(mechanism["reactions"][47])
		#print(perturb)
		#raise AssertionError("Stop!!")
		
		for rxn in self.rxn_list:
			beta = np.asarray(perturb[rxn])
			type_of_rxn = rxn_type[rxn]
			
			if type_of_rxn == "Elementary":
				index = self.unsrt[rxn].index
				#a = mechanism["reactions"][index]
				#print(f"Before\t{a}\n")
				#print(index)
				new_mechanism = self.ElementaryPerturbation(rxn,beta,selection[rxn],mechanism)
				mechanism = new_mechanism
				#mechanism["reactions"][index]
				#print(f"After\t{b}\n")
				#print(a==b)
			if type_of_rxn == "PLOG":
				index = self.unsrt[rxn].index
				new_mechanism = self.PlogPerturbation(rxn,beta,selection[rxn],mechanism)
				mechanism = new_mechanism
				
			if type_of_rxn == "PLOG-Duplicate":
				index = self.unsrt[rxn].index
				new_mechanism = self.PlogPerturbation(rxn,beta,selection[rxn],mechanism)
				mechanism = new_mechanism
				
			if type_of_rxn == "Duplicate":
				index = self.unsrt[rxn].index
				new_mechanism = self.DuplicatePerturbation(rxn,beta,selection[rxn],mechanism)
				mechanism = new_mechanism
				
			if type_of_rxn == "Falloff":
				index = self.unsrt[rxn].index
				new_mechanism = self.TroePerturbation(rxn,beta,selection[rxn],mechanism)
				mechanism = new_mechanism
				
			if type_of_rxn == "ThirdBody":
				index = self.unsrt[rxn].index
				new_mechanism = self.ElementaryPerturbation(rxn,beta,selection[rxn],mechanism)
				mechanism = new_mechanism
				
			if type_of_rxn == "BranchingRxn":
				index = self.unsrt[rxn].index
				new_mechanism = self.BranchingReactions(rxn,beta,selection[rxn],mechanism)
				mechanism = new_mechanism
		#string = yaml.dump(new_mechanism,default_flow_style=False)
		#print(mechanism == self.mechanism)
		#print(mechanism["reactions"][47])
		#raise AssertionError("Elementary rxn")
		return mechanism,perturb
	
	
