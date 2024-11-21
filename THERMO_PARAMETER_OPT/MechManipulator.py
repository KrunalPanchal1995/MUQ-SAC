import os
import numpy as np
from copy import deepcopy
from MechanismParser import Parser

class Manipulator:
	def __init__(self, copy_of_mech, unsrt_object, perturbation, perturbation_type="opt", parameter_dict=None,flag="reaction"):
		self.mechanism = deepcopy(copy_of_mech)
		self.unsrt = unsrt_object
		self.perturbation_type = perturbation_type
		self.parameter_dict = parameter_dict
		self.flag = flag
		if flag == "reaction":
			self.rxn_type = parameter_dict["type"]
			self.rxn_data = parameter_dict["data"]
			self.rxn_dict = parameter_dict["reaction"]
			self.rxn_list = [parameter_dict["reaction"][index] for index in parameter_dict["reaction"]]
		self.perturbation = perturbation

	import numpy as np

	def calculate_enthalpy_and_entropy(self,T, a1, a2, a3, a4, a5, a6, a7, R=8.314):
		"""
		Calculate enthalpy and entropy at a given temperature T using NASA polynomial coefficients.
		R is in cal/(mol K).
		"""
		# Calculate enthalpy (H) at temperature T
		H_RT = (
		a1 + a2 * T / 2 + a3 * T**2 / 3 +
		a4 * T**3 / 4 + a5 * T**4 / 5 + a6 / T
		)
		H = H_RT * R * T  # Returns enthalpy H in cal/mol

		# Calculate entropy (S) at temperature T
		S_R = (
		a1 * np.log(T) + a2 * T + a3 * T**2 / 2 +
		a4 * T**3 / 3 + a5 * T**4 / 4 + a7
		)
		S = S_R * R  # Returns entropy S in 

		return H, S

	def adjust_a6_a7_for_perturbation(self,coefficients):
		"""
		Adjust a6 and a7 to ensure H_mod(298 K) and S_mod(298 K) match original values
		after perturbing a1-a5.

		Parameters:
		coefficients: array-like of length 7 [a1, a2, a3, a4, a5, a6, a7]

		Returns:
		Modified a6 and a7 to maintain enthalpy and entropy values at 298 K.
		"""
		T_ref = 298  # Reference temperature in Kelvin
		R = 8.314  # Gas constant in J/(mol K)
		a1, a2, a3, a4, a5, a6, a7 = coefficients

		# Perturb a1 to a5 by multiplying by 2
		#a1_mod = 2 * a1
		#a2_mod = 2 * a2
		#a3_mod = 2 * a3
		#a4_mod = 2 * a4
		#a5_mod = 2 * a5

		# Calculate original enthalpy and entropy at T_ref
		#H_original, S_original = self.calculate_enthalpy_and_entropy(T_ref, a1, a2, a3, a4, a5, a6, a7, R)

		# Calculate modified enthalpy and entropy at T_ref with perturbed coefficients
		#H_mod, S_mod = self.calculate_enthalpy_and_entropy(T_ref, a1_mod, a2_mod, a3_mod, a4_mod, a5_mod, a6, a7, R)

		# Calculate the differences needed to adjust a6 and a7
		#delta_H = H_original - H_mod
		#delta_S = S_original - S_mod

		# Adjust a6 and a7 based on the delta values
		a6_mod = a6 - T_ref*(a1 + a2 * T_ref/ 2 + a3 * T_ref**2 / 3 + a4 * T_ref**3 / 4 + a5 * T_ref**4 / 5)
		a7_mod = a7 - R*(a1 * np.log(T_ref) + a2 * T_ref + a3 * T_ref**2 / 2 + a4 * T_ref**3 / 3 + a5 * T_ref**4 / 4)

		return a6_mod, a7_mod

	def perturb_cp(self, index, beta, mechanism):
		thermo_data = mechanism["species"][index]['thermo']['data']
		a6_mod, a7_mod = self.adjust_a6_a7_for_perturbation(thermo_data[0])
		thermo_data[0][5] = float(a6_mod)
		thermo_data[0][6] = float(a7_mod)
		for i in range(2):
			thermo_data[i][:5] = [float(value) * float(beta) for value in thermo_data[i][:5]]
		
		mechanism["species"][index]['thermo']['data']	 = deepcopy(thermo_data)
		return mechanism

	def perturb_enthalpy(self, index, beta, mechanism):
		thermo_data = mechanism["species"][index]['thermo']['data']	
		for i in range(2):
			thermo_data[i][5] += float(beta)  # Perturb the 7th entry
		mechanism["species"][index]['thermo']['data']	 = deepcopy(thermo_data)
		return mechanism
		
	def perturb_entropy(self, index, beta, mechanism):
		thermo_data = mechanism["species"][index]['thermo']['data']
		for i in range(2):
			thermo_data[i][6] += float(beta)  # Perturb the 7th entry
		mechanism["species"][index]['thermo']['data'] = deepcopy(thermo_data)
		return mechanism
		
	def _extract_species_data(self, parameter_dict):
		self.species_data = {}
		for species in parameter_dict:
			for index,dict_ in enumerate(self.mechanism["species"]):
				if dict_["name"] == species:
					self.species_data[species] = index


	def getRxnDetails(self):
		rxn_dict = {}
		rxn_data = self.mechanism["reactions"]
		for rxn in self.rxn_list:
			new_rxn_data = {}
			temp = []
			index_ = []
			for index, data in enumerate(rxn_data):
				if rxn == data["equation"]:
					temp.append(data)
					index_.append(index)
			new_rxn_data["temp"] = temp
			new_rxn_data["index"] = index_
			rxn_dict[rxn] = new_rxn_datamani
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
						elif data["type"] == "falloff":
							rxn_type[data["equation"]] = "Falloff"
						elif data["type"] == "pressure-dependent-Arrhenius" and "duplicate" not in data:
							rxn_type[data["equation"]] = "PLOG"
						elif data["type"] == "pressure-dependent-Arrhenius" and "duplicate" in data:
							rxn_type[data["equation"]] = "PLOG-Duplicate"
							break
					elif "duplicate" in data:
						rxn_type[data["equation"]] = "Duplicate"
						break
					else:
						rxn_type[data["equation"]] = "Elementary"
		return rxn_type

	def ElementaryPerturbation(self, index, beta, mechanism):
		perturbation_factor = beta
		reaction_details = mechanism["reactions"][index]["rate-constant"]
		pre_exponential_factor = np.log(float(reaction_details["A"]))
		reaction_details["A"] = float(np.exp(pre_exponential_factor + perturbation_factor))
		mechanism["reactions"][index]["rate-constant"] = deepcopy(reaction_details)
		return mechanism

	def PlogPerturbation(self, index, beta, mechanism):
		if beta == 0:
			perturbation_factor = 1.0
		else:
			perturbation_factor = np.exp(beta)
		reaction_details = mechanism["reactions"][index]["rate-constants"]
		new_rxn_details = []
		for rxn in reaction_details:
			temp = {
				"P": rxn["P"],
				"A": float(float(rxn["A"]) * perturbation_factor),
				"b": rxn["b"],
				"Ea": rxn["Ea"]
			}
			new_rxn_details.append(temp)
		mechanism["reactions"][index]["rate-constants"] = deepcopy(new_rxn_details)
		return mechanism

	def DupPlogPerturbation(self, rxn_object, beta, mechanism):
		rxn_object_a = {"temp": [], "index": []}
		rxn_object_a["temp"].append(rxn_object["temp"][0])
		rxn_object_a["index"].append(int(rxn_object["index"][0]) - 1)
		rxn_object_b = {"temp": [], "index": []}
		rxn_object_b["temp"].append(rxn_object["temp"][1])
		rxn_object_b["index"].append(int(rxn_object["index"][1]) - 1)
		new_mechanism = self.PlogPerturbation(rxn_object_a, beta, mechanism)
		new_mechanism = self.PlogPerturbation(rxn_object_b, beta, new_mechanism)
		return new_mechanism

	def DupElementaryPerturbation(self, rxn_object, beta, mechanism):
		rxn_object_a = {"temp": [], "index": []}
		rxn_object_a["temp"].append(rxn_object["temp"][0])
		rxn_object_a["index"].append(int(rxn_object["index"][0]) - 1)
		rxn_object_b = {"temp": [], "index": []}
		rxn_object_b["temp"].append(rxn_object["temp"][1])
		rxn_object_b["index"].append(int(rxn_object["index"][1]) - 1)
		new_mechanism = self.ElementaryPerturbation(rxn_object_a, beta, mechanism)
		new_mechanism = self.ElementaryPerturbation(rxn_object_a, beta, new_mechanism)
		return new_mechanism

	def TroePerturbation(self, index, beta, mechanism):
		perturbation_factor = beta
		reaction_details_low = mechanism["reactions"][index]["low-P-rate-constant"]
		pre_exponential_factor_low = np.log(float(reaction_details_low["A"]))
		reaction_details_low["A"] = float(np.exp(pre_exponential_factor_low + perturbation_factor))
		mechanism["reactions"][index]["low-P-rate-constant"] = deepcopy(reaction_details_low)
		return mechanism

	def doPerturbation(self):  # purtubing data.... flag: 'reaction' for reaction perturbation, 'thermo' for thermo perturbation.
		if self.flag == "thermo":
			mechanism = self.mechanism
			self._extract_species_data(self.parameter_dict)
			perturb = ""
			count = 0
			for species in self.species_data:
				beta = self.perturbation[3 * count:3 * count + 3]  # using the purturbation array to modify cp, h and s
				index = self.species_data[species]
				if float(abs(beta[0])) > 0:
					perturb += f"{species}_cp\n"
					perturb += f"{beta[0]}"
					mechanism = self.perturb_cp(index, beta[0],mechanism)
				if float(abs(beta[1])) > 0:
					perturb += f"{species}_H\n"
					perturb += f"{beta[1]}"
					mechanism = self.perturb_enthalpy(index, beta[1],mechanism)
				if float(abs(beta[2])) > 0:
					perturb += f"{species}_S\n"
					perturb += f"{beta[2]}"
					mechanism = self.perturb_entropy(index, beta[2],mechanism)
				count+=1
			return mechanism,perturb

		elif self.flag == "reaction" :
			rxn_type = self.rxn_type
			rxn_data = self.rxn_data
			mechanism = self.mechanism
			rxn_dict = self.rxn_dict
			perturb = ""
			for i, index in enumerate(rxn_dict):
				rxn = rxn_dict[index]
				index_ = index - 1
				beta = np.asarray(self.perturbation[i])
				if float(abs(beta)) > 0:
					perturb += f"{rxn}\t{beta}"
					type_of_rxn = rxn_type[index]
					data = rxn_data[rxn]

					if type_of_rxn == "Elementary":
						new_mechanism = self.ElementaryPerturbation(index_, beta, mechanism)

					elif type_of_rxn == "PLOG":
						new_mechanism = self.PlogPerturbation(index_, beta, mechanism)

					elif type_of_rxn == "PLOG-Duplicate":
						new_mechanism = self.DupPlogPerturbation(data, beta, mechanism)

					elif type_of_rxn == "Duplicate":
						new_mechanism = self.DupElementaryPerturbation(data, beta, mechanism)

					elif type_of_rxn == "ThirdBody":
						new_mechanism = self.ElementaryPerturbation(index_, beta, mechanism)

					elif type_of_rxn == "Falloff":
						new_mechanism = self.TroePerturbation(index_, beta, mechanism)
			return new_mechanism,perturb
		else:
			raise AssertionError(f"Invalid flag: {self.flag}!!\n\t-valid flag types ['thermo','reaction']\n")

