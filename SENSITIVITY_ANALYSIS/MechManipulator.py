import os
import numpy as np
from copy import deepcopy
from MechanismParser import Parser
###Tasks
# Take a copy of the mechanism in yaml format and manipulate it
# After manipulating, convert it back to ck format
####################################################3
###

class Manipulator:
	def __init__(self,copy_of_mech,unsrt_object,perturbation,perturbation_type = "opt"):
		#self.mechanism = copy_of_mech
		self.mechanism = deepcopy(copy_of_mech)
		self.unsrt = unsrt_object
		self.rxn_type = unsrt_object["type"]
		self.rxn_data = unsrt_object["data"]
		self.perturbation = perturbation
		self.rxn_list = [self.unsrt["reaction"][index] for index in self.unsrt["reaction"]]
		#print(self.rxn_list)
		self.rxn_dict = unsrt_object["reaction"]
		self.perturbation_type = perturbation_type
	def getRxnDetails(self):
		rxn_dict = {}
		rxn_data = self.mechanism["reactions"]
		for rxn in self.rxn_list:
			new_rxn_data = {}
			temp = []
			index_ = []
			for index,data in enumerate(rxn_data):
				if rxn == data["equation"]:
					temp.append(data)
					index_.append(index)
			new_rxn_data["temp"] = temp
			new_rxn_data["index"] = index_
			rxn_dict[rxn] = new_rxn_data
		return rxn_dict
	
	def getRxnType(self):
		rxn_type = {}
		rxn_data = self.mechanism["reactions"]
		for rxn in self.rxn_list:
			for data in rxn_data:
				if rxn in data["equation"]:
					if "type" in data:
						if data["type"] == "three-body":
							rxn_type[data["equation"]] = "ThirdBody"
							#print("Rxn is three-body")
						elif data["type"] == "falloff":
							rxn_type[data["equation"]] = "Falloff"
							#print("Rxn is fall-off")
						elif data["type"] == "pressure-dependent-Arrhenius" and "duplicate" not in data:
							rxn_type[data["equation"]] = "PLOG"
							#print("Rxn is Plog")
						elif data["type"] == "pressure-dependent-Arrhenius" and "duplicate" in data:
							rxn_type[data["equation"]] = "PLOG-Duplicate"
							break
							#print("Rxn is Plog duplicate")
					elif "duplicate" in data:
						rxn_type[data["equation"]] = "Duplicate"
						break
						#print("Rxn is duplicate")
					
					else:
						rxn_type[data["equation"]] = "Elementary"
						#print("Rxn is elementary")
		return rxn_type
	
	
	def ElementaryPerturbation(self,index,beta,mechanism):
		#index = rxn_object["index"][0]	
		perturbation_factor = beta
		#print(perturbation_factor)	
		"""
		The new perturbed reaction replaces the prior Arrhenius parameters 
		"""
		#print("###############################")
		#print("\n\t\tBefore")
		#print(f"\n{mechanism['reactions'][index]}")
		#print(mechanism["reactions"][index]["rate-constant"])
		#print(index)
		reaction_details = mechanism["reactions"][index]["rate-constant"]
		pre_exponential_factor = np.log(float(reaction_details["A"]))
		reaction_details["A"] = float(np.exp(pre_exponential_factor+perturbation_factor))
		mechanism["reactions"][index]["rate-constant"] = deepcopy(reaction_details)
		#print("\n\t\tAfter")
		#print(mechanism["reactions"][index]["rate-constant"])
		#print("###############################")
		#raise AssertionError("Stop")
		return mechanism
		
	def PlogPerturbation(self,index,beta,mechanism):
		#index = rxn_object["index"][0]	
		if beta == 0:
			perturbation_factor = 1.0
		else:
			perturbation_factor = np.exp(beta)
			
		#print("###############################")
		#print("\n\t\tBefore")
		#print(f"\n{mechanism['reactions'][index]}")
		#print(rxn_object)
		reaction_details = mechanism["reactions"][index]["rate-constants"]
		#print(reaction_details)
		new_rxn_details = []
		#print(mechanism["reactions"][index])
		for rxn in reaction_details:
			temp = {}
			temp["P"] = rxn["P"]
			temp["A"] = float(float(rxn["A"])*perturbation_factor)
			temp["b"] = rxn["b"]
			temp["Ea"] = rxn["Ea"]
			new_rxn_details.append(temp)
		mechanism["reactions"][index]["rate-constants"] = deepcopy(new_rxn_details)
		#print("\n\t\tAfter")
		#print(mechanism["reactions"][index]["rate-constants"])
		#print("###############################")
		#raise AssertionError("Stop")
		return mechanism
	
	def DupPlogPerturbation(self,rxn_object,beta,mechanism):
		rxn_object_a = {}
		rxn_object_a["temp"] = []
		rxn_object_a["index"] = []
		rxn_object_a["temp"].append(rxn_object["temp"][0])
		rxn_object_a["index"].append(int(rxn_object["index"][0])-1)
		rxn_object_b = {}
		rxn_object_b["temp"] = []
		rxn_object_b["index"] = []
		rxn_object_b["temp"].append(rxn_object["temp"][1])
		rxn_object_b["index"].append(int(rxn_object["index"][1])-1)

		new_mechanism = self.PlogPerturbation(rxn_object_a,beta,mechanism)
		new_mechanism = self.PlogPerturbation(rxn_object_b,beta,new_mechanism)
		return new_mechanism
	
	def DupElementaryPerturbation(self,rxn_object,beta,mechanism):	
		rxn_object_a = {}
		rxn_object_a["temp"] = []
		rxn_object_a["index"] = []
		rxn_object_a["temp"].append(rxn_object["temp"][0])
		rxn_object_a["index"].append(int(rxn_object["index"][0])-1)
		rxn_object_b = {}
		rxn_object_b["temp"] = []
		rxn_object_b["index"] = []
		rxn_object_b["temp"].append(rxn_object["temp"][1])
		rxn_object_b["index"].append(int(rxn_object["index"][1])-1)
		
		new_mechanism = self.ElementaryPerturbation(rxn_object_a,beta,mechanism)
		new_mechanism = self.ElementaryPerturbation(rxn_object_a,beta,new_mechanism)
		return new_mechanism
		
	def TroePerturbation(self,index,beta,mechanism):
		#index = rxn_object["index"][0]	
		perturbation_factor = beta	
		#print("###############################")
		#print("\n\t\tBefore")
		#print(f"\n{mechanism['reactions'][index]}")
		#print(mechanism["reactions"][index])
		reaction_details_low = mechanism["reactions"][index]["low-P-rate-constant"]
		#reaction_details_high = mechanism["reactions"][index]["high-P-rate-constant"]
		pre_exponential_factor_low = np.log(float(reaction_details_low["A"]))
		#pre_exponential_factor_high = np.log(float(reaction_details_high["A"]))
		
		reaction_details_low["A"] = float(np.exp(pre_exponential_factor_low+perturbation_factor))
		#reaction_details_high["A"] = float(np.exp(pre_exponential_factor_high+perturbation_factor))
		
		mechanism["reactions"][index]["low-P-rate-constant"] = deepcopy(reaction_details_low)
		#mechanism["reactions"][index]["high-P-rate-constant"] = deepcopy(reaction_details_high)
		#print("\n\t\tAfter")
		#print(mechanism["reactions"][index])
		#print("###############################")
		#raise AssertionError("Stop")
		return mechanism
	
	def doPerturbation(self):
		rxn_type = self.rxn_type
		rxn_data = self.rxn_data
		mechanism = self.mechanism
		rxn_dict = self.rxn_dict
		string = ""
		for index in rxn_type:
			#print(index,rxn_type[index],rxn_dict[index])
			string+=f"{index}\t{rxn_type[index]}\t{rxn_dict[index]}\n"
		file_open =  open("rxn_type.txt","w").write(string)
		
		"""
		Perturbations happens for each reactions
			- Reaction index is identified
			- Duplicate reactions is treated as elementary reactions
		"""
		perturb = ""
		for i,index in enumerate(rxn_dict):
			#print(self.perturbation)
			rxn = rxn_dict[index]
			index_ = index-1
			beta = np.asarray(self.perturbation[i])
			if float(abs(beta)) > 0.0:
				perturb = f"{rxn}\t{beta}"
				type_of_rxn = rxn_type[index]
				data = rxn_data[rxn]
				
				if type_of_rxn == "Elementary":
					new_mechanism = self.ElementaryPerturbation(index_,beta,mechanism)
					mechanism = new_mechanism
					
				if type_of_rxn == "PLOG":
					new_mechanism = self.PlogPerturbation(index_,beta,mechanism)
					mechanism = new_mechanism
					
				if type_of_rxn == "PLOG-Duplicate":
					new_mechanism = self.DupPlogPerturbation(data,beta,mechanism)
					mechanism = new_mechanism
					
				if type_of_rxn == "Duplicate":
					new_mechanism = self.DupElementaryPerturbation(data,beta,mechanism)
					mechanism = new_mechanism
					
				if type_of_rxn == "Falloff":
					new_mechanism = self.TroePerturbation(index_,beta,mechanism)
					mechanism = new_mechanism
					
				if type_of_rxn == "ThirdBody":
					new_mechanism = self.ElementaryPerturbation(index_,beta,mechanism)
					mechanism = new_mechanism
			
		return mechanism,perturb
			
