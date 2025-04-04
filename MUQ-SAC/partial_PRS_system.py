import os,sys
import numpy as np
import DesignMatrix as DM
import simulation_manager2_0 as simulator

class PartialPRS(object):
	def __init__(self,sensitivity_dict,unsrt_data,optInputs,target_list,case_index,active_parameters,design):
		#print(len(sensitivity_dict),len(active_parameters))
		self.target_list = target_list
		self.optInputs = optInputs
		self.cut_off = float(optInputs["Stats"]["cut_off_percentage"])
		self.design = design
		self.case_index = case_index
		self.unsrt = unsrt_data
		self.s_A = []
		self.s_n = []
		self.s_Ea = []
		self.no_of_sim = None
		#print(sensitivity_dict)
		sens_SA_list = []
		sens_SA_dict = {}
		count = 0
		if self.design == "A-facto":
			for rxn in sensitivity_dict:
				for rxn_ in unsrt_data:
					if rxn == rxn_.split(":")[0]:
						#print(rxn_.split(":")[0])
						#print(rxn_,sensitivity_dict[rxn])
						sens_SA_list.append(rxn_)
						sens_SA_dict[unsrt_data[rxn_].activeParameters[0]] = float(sensitivity_dict[rxn])
						self.s_A.append(float(sensitivity_dict[rxn]))
						count+=1
			
			
			#print(len(self.s_A))
			self.active_params = active_parameters
			self.partial_active = {}
			self.partial_active_list = []
			self.selected = []
			self.coeff = []
			self.abs_coeff = []
			count = 0
			#getting the coeff according to the order of activeParameters list (unsrt_data)
			for ac in self.active_params:
				for index,sc in enumerate(sens_SA_dict):
					if ac == sc:
						self.coeff.append(sens_SA_dict[sc])
						self.abs_coeff.append(abs(sens_SA_dict[sc]))
			
			coeff_sum = sum(self.abs_coeff)
			self.normalized_coeff = np.asarray(self.abs_coeff)*100/coeff_sum
			#print(len(self.s_A),len(self.active_params))
			
			self.selected_rxn_string = ""
			#print(sens_SA_list)
			#print(self.active_params)
			#print(sens_SA_dict)
			#print(self.coeff)
			for ind,active_params in enumerate(self.active_params):
				if abs(self.coeff[ind])>=self.cut_off/100:
					#print(active_params,self.coeff[ind])
					self.partial_active[active_params] = 1
					self.partial_active_list.append(active_params)
					self.selected.append(1)
					self.selected_rxn_string+=f"{active_params}\n"
					
				else:
					self.partial_active[active_params] = 0
					self.selected.append(0)
			self.selected_rxn_count = sum(self.selected)
			#raise AssertionError("Stop!")
		else:
		
			for rxn in sensitivity_dict:
				self.s_A.append(float(sensitivity_dict[rxn][0]))
				self.s_n.append(float(sensitivity_dict[rxn][1]))
				self.s_Ea.append(float(sensitivity_dict[rxn][2]))
			self.active_params = active_parameters
			self.partial_active = {}
			self.partial_active_list = []
			self.selected = []
			self.coeff = []
			self.abs_coeff = []
			count = 0
			for index,sa in enumerate(self.s_A):
				self.coeff.append(sa)
				self.abs_coeff.append(abs(sa))
				self.coeff.append(self.s_n[index])
				self.abs_coeff.append(abs(self.s_n[index]))
				self.coeff.append(self.s_Ea[index])
				self.abs_coeff.append(abs(self.s_Ea[index]))
			
			coeff_sum = sum(self.abs_coeff)
			self.normalized_coeff = np.asarray(self.abs_coeff)*100/coeff_sum
			
			self.selected_rxn_string = ""
			for ind,active_params in enumerate(self.active_params):
				if self.normalized_coeff[ind]>self.cut_off:
					self.partial_active[active_params] = 1
					self.partial_active_list.append(active_params)
					self.selected.append(1)
					self.selected_rxn_string+=f"{active_params}\n"
				else:
					self.partial_active[active_params] = 0
					self.selected.append(0)
			self.selected_rxn_count = sum(self.selected)
	def getTotalUnknowns(self,N):
		n_ = 1 + 2*N + (N*(N-1))/2
		return int(n_)
	
	def getSim(self,n,design):
		n_ = self.getTotalUnknowns(n)
		if design == "A-facto":
			sim = 4*n_
		else:
			sim = 7*n_	
		return sim
	
	def partial_DesignMatrix(self):
		na = len(self.partial_active_list)
		n_rxn = len(self.active_params)
		self.no_of_sim = self.getSim(na,self.design)
		print("\n################################################\n###  Starting to generate Design Matrix  ###\n###  for all targets  ###\n################################################\n")
		print(f"\n[Case-{self.case_index}]\n\tNo. of Simulations required: {self.getSim(na,self.design)}\n\tNo. of selected reactions: {self.selected_rxn_count}\n")
		if "partial" not in os.listdir():
			os.mkdir("partial")
			os.mkdir(f"partial/{self.case_index}")
		else:
			if f"{self.case_index}" not in os.listdir("partial/"):
				os.mkdir(f"partial/{self.case_index}")
		
		#######################################################################
		### Goes to DesignMatrix Modules to create Matrix for Partial PRS ###
		#######################################################################
		
		if "DesignMatrix.csv" not in os.listdir(f"partial/{self.case_index}/"):
			design_matrix,p_design_matrix = DM.DesignMatrix(self.unsrt,self.design,self.getSim(na,self.design),n_rxn).getSample_partial(self.case_index,self.selected)
			g = open(f"partial/{self.case_index}/selected_parameters.csv","w").write(self.selected_rxn_string)
		else:
			design_matrix_file = open(f"partial/{self.case_index}/DesignMatrix.csv").readlines()
			p_design_matrix_file = open(f"partial/{self.case_index}/pDesignMatrix.csv").readlines()
			design_matrix = []
			for row in design_matrix_file:
				design_matrix.append([float(ele) for ele in row.strip("\n").strip(",").split(",")])
			p_design_matrix = []
			for row in p_design_matrix_file:
				p_design_matrix.append([float(ele) for ele in row.strip("\n").strip(",").split(",")])
		
		####################################################################################################
		### Goes to Simulation Manager to generate YAML mechanisms based on DesignMatrix for Partial PRS ###
		####################################################################################################
		
		SSM = simulator.SM(self.target_list,self.optInputs,self.unsrt,design_matrix)
		if f"partial_Perturbed_Mech" not in os.listdir():
			os.mkdir("partial_Perturbed_Mech")
		if f"{self.case_index}" not in os.listdir("partial_Perturbed_Mech/"):
			os.mkdir(f"partial_Perturbed_Mech/{self.case_index}")
			print("\nPerturbing the Mechanism files\n")
			
			chunk_size = 500
			params_yaml = [design_matrix[i:i+chunk_size] for i in range(0, len(design_matrix), chunk_size)]
			count = 0
			yaml_loc = []
			for params in params_yaml:
				
				yaml_list = SSM.getYAML_List(params)
				#yaml_loc = []
				location_mech = []
				index_list = []
				for i,dict_ in enumerate(yaml_list):
					index_list.append(str(count+i))
					location_mech.append(os.getcwd()+f"/partial_Perturbed_Mech/{self.case_index}/")
					yaml_loc.append(os.getcwd()+f"/partial_Perturbed_Mech/{self.case_index}/mechanism_"+str(count+i)+".yaml")
				count+=len(yaml_list)
				#gen_flag = False
				#SSM.getPerturbedMechLocation(yaml_list,location_mech,index_list)
				SSM.getPerturbedMechLocation(yaml_list,location_mech,index_list)
				print(f"\nGenerated {count} files!!\n")
			print("\nGenerated the YAML files required for simulations!!\n")
		else:
			print("\nYAML files already generated!!")
			yaml_loc = []
			location_mech = []
			index_list = []
			for i,sample in enumerate(design_matrix):
				index_list.append(i)
				location_mech.append(os.getcwd()+f"/partial_Perturbed_Mech/{self.case_index}")
				yaml_loc.append(os.getcwd()+f"/partial_Perturbed_Mech/{self.case_index}/mechanism_"+str(i)+".yaml")
		return yaml_loc,p_design_matrix,self.selected
