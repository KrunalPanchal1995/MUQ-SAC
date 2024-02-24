import os
import numpy as np
import matplotlib.pyplot as plt


class combustion_dataset():
	def __init__(self,set_name,Target_Dataset,optimized_params):
		#print(Target_Dataset)
		#self.dataSet = dataSet
		#self.dataSetTar = dtsTarget_object
		#print(optimized_params)
		self.tag = set_name
		self.opt_params = optimized_params
		self.opt_target_values = []
		self.dataType = []
		self.Fuel = []
		self.Oxidizer = []
		self.bathgas = []
		self.bathgas_x = []
		self.addendum = []
		self.BG1 = []
		self.BG2 = []
		self.BG3 = []
		self.Fuel_conc = []
		self.Oxidizer_conc = []
		self.BG1_conc = []
		self.BG2_conc = []
		self.BG3_conc = []
		self.Pressure = []
		self.observed = []
		self.Temperature = []
		self.Phi = []
		self.Unsrt_fct = []
		self.Measurnment = []
		self.Simulation = []
		self.ignition_mode = []
		self.sensitivitySimulations = []
		self.nominalSimulations = []
		self.optimizedSimulations = []
		self.sens_analysis = []
		self.responseSurfces = []
		self.blackBoxResponse = []
		self.startProfile = []
		self.flow_rate = []
		self.fuel_unit = []
		self.t_unit = []
		self.p_unit = []
		self.target_unit = []
		self.flow_unit = []
		self.unit = []		
		self.fls_id = []
		self.fls_x = []
		self.set_id = []
		self.prior_model_unsrt = []
		self.posterior_model_unsrt = []
		self.prior_model_response = []
		self.posterior_model_response = []
		for Target_object in Target_Dataset:
			#print(Target_object)
			#Target_list = Target_Dataset[str(Target_type)]
			#for Target_object in Target_list:
			#print(type(self.opt_params))
			if self.opt_params != []:
				#print(Target_object["Target_database"])
				self.opt_target_values.append(Target_object["Target_database"].calculated_target_value(self.opt_params))
			self.dataType.append(Target_object["Target_database"].target)
			self.Fuel.append(Target_object["Target_database"].fuel_id)
			self.Oxidizer.append(Target_object["Target_database"].oxidizer)
			self.set_id.append(Target_object["Target_database"].dataSet_id)
			#self.BG1.append(Target_object["Target_database"].bath_gas1)
			#self.BG2.append(Target_object["Target_database"].bath_gas2)
			#self.BG3.append(Target_object["Target_database"].bath_gas3)
			self.Fuel_conc.append(Target_object["Target_database"].fuel_x)
			self.Oxidizer_conc.append(Target_object["Target_database"].oxidizer_x)
			self.bathgas.append(Target_object["Target_database"].bath_gas)
			self.bathgas_x.append(Target_object["Target_database"].bath_gas_x)
			#self.BG1_conc.append(Target_object["Target_database"].bath_gas1_x)
			#self.BG2_conc.append(Target_object["Target_database"].bath_gas2_x)
			#self.BG3_conc.append(Target_object["Target_database"].bath_gas3_x)
			self.Pressure.append(Target_object["Target_database"].pressure)
			self.Temperature.append(Target_object["Target_database"].temperature)
			self.Phi.append(Target_object["Target_database"].phi)
			if "Fls" in self.dataType:
				if Target_object["Target_database"].phi == "":
					for i in Target_object["Target_database"].fuel_id:
						string.append("X_"+str(i))
						value.append(float(Target_object["Target_database"].fuel_x[i]))
					self.fls_id.append(("+".join(string)))
					self.fls_x.append(np.sum(np.asarray(value)))
				else:
					self.fls_id.append("Equivalence ratio")
					self.fls_x.append(Target_object["Target_database"].phi)
			self.addendum.append(Target_object["Target_database"].add)
			self.Unsrt_fct.append(np.abs(Target_object["Target_database"].std_dvtn))
			self.Measurnment.append(Target_object["Target_database"].measure_type)
			self.Simulation.append(Target_object["Target_database"].simulation)
			self.ignition_mode .append(Target_object["Target_database"].ig_mode)		
			#self.startProfile.append(Target_object["Target_database"].s_p_name)
			#self.flow_rate.append(Target_object["Target_database"].flow_unit)
			#self.fuel_unit.append(Target_object["Target_database"].f_unit)
			#self.t_unit.append(Target_object["Target_database"].t_unit)
			#self.p_unit.append(Target_object["Target_database"].p_unit)
			self.observed.append((Target_object["Target_database"].observed))
			self.target_unit.append(Target_object["Target_database"].target_unit)
			#self.flow_unit .append(Target_object["Target_database"].ig_mode)
			self.unit.append(Target_object["Target_database"].units)
			self.sensitivitySimulations.append(Target_object["SA_simulations"])
			self.nominalSimulations.append(float(Target_object["Original_simulations"].split("\t")[1].strip("\n")))
			
			#self.optimizedSimulations.append(float(Target_object["Optimized_simulations"].split("\t")[1].strip("\n")))
			if Target_object["optimization_type"] == "Direct":
				#print("Response surface was not used")
				self.optimization_type = "Direct"
			else:			
				self.optimizedSimulations.append(Target_object["Opt_simulations"])
				self.sens_analysis.append(Target_object["SA_PRS"])
				self.responseSurfces.append(Target_object["Opt_PRS"])
			self.prior_model_unsrt.append(Target_object["Prior_model_unsrt"])
			#self.posterior_model_unsrt.append(Target_object["Posterior_model_unsrt"])
			if Target_object["Prior_model_response"] != None:
				self.prior_model_response.append(Target_object["Prior_model_response"])
				self.posterior_model_response.append(Target_object["Posterior_model_response"])
			else:
				self.prior_model_response.append(0)
				self.posterior_model_response.append(0)
		#self.Unsrt_fct = sorted(self.Unsrt_fct,key=float)
		#self.Temperature = sorted(self.Temperature,key=float)
		#self.opt_target_values = sorted(self.opt_target_values,key=float)
		#self.nominalSimulations = sorted(self.nominalSimulations,key=float)
		#self.observed = sorted(self.observed,key=float)
		
#		print(self.Unsrt_fct)
#		print(self.observed)
#		print(self.nominalSimulations)
#		print(self.opt_target_values)
		#for target in Target_object:
		#	print(target)

#				
#		for target in dtaTarget_object:
#			self.dataType.append(target.target)
#			self.dtsFuel.append(target.fuel) 
#			self.dtsOxidizer.append(target.oxidizer) 
#			self.dtsBG1.append(target.bath_gas1)
#			self.dtsBG2.append(target.bath_gas2) 
#			self.dtsBG3.append(target.bath_gas3) 
#			self.dtsPressure.append(target.pressure) 
#			self.dtsTemperature.append(target.temperature)
#			self.dtsPhi.append(target.phi)
#			self.dtsUnsrt_fct.append(target.unsrt_fct)
#			self.dtsMeasurnment.append(target.measure_type)
#			self.dts_Simulation.append(target.simulation)
#	def getModelResponses(self,responseSurface):
#		
#	def getBlackBoxResponse(self,solverResponse):
#	
#	for case in case_dir:
#		if Target_object[case].target == "Tig":
#			DataSet_Tig.append(Target_object[case].d_set)
#		if Target_object[case].target == "Fsl":
#			DataSet_Fsl.append(Target_object[case].d_set)
#		if Target_object[case].target == "Scp":
#			DataSet_Scp.append(Target_object[case].d_set)

#	DataSet_Tig = list(OrderedDict.fromkeys(DataSet_Tig))
#	DataSet_Fsl = list(OrderedDict.fromkeys(DataSet_Fsl))
#	DataSet_Scp = list(OrderedDict.fromkeys(DataSet_Scp))


#	for ind in DataSet_Tig:
#		temp_T = []
#		temp_Tig = []
#		temp_P = []
#		temp_sigma_exp = []
#		temp_error = []
#		temp_phi = []
#		temp_opt_sim = []
#		temp_init_sim = []
#		prior_unsrt_res = []
#		posterior_unsrt_res = []
#		for case in Target_object:
#			if case.target == "Tig":
#				if case.d_set == ind:
#					temp_T.append(case.temperature)
#					temp_Tig.append(case.observed)
#					temp_sigma_exp.append(case.std_dvtn)
#					temp_opt_sim.append(case.calculated_target_value(opt))
#					temp_init_sim.append(case.calculated_target_value(init_guess))
#					temp_P.append(case.pressure)
#					temp_phi.append(case.phi)
#					#1/9 if specified
#					prior_cov = (1)*np.identity(len(opt))
#					prior_unsrt_res.append(case.model_response_uncertainty(init_guess,prior_cov,2))
#					posterior_unsrt_res.append(case.model_response_uncertainty(opt,Posterior_cov,2)) #Error in micro seconds
		#print(prior_unsrt_res)
		#print(posterior_unsrt_res)
		
		
