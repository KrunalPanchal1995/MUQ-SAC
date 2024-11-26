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
				#print()
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
			self.cholskyDeCorrelateMat = np.array([min(self.uncertainties)])
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
		T = np.array([self.temperatures[0],(self.temperatures[0]+self.temperatures[-1])/2,self.temperatures[-1]])
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
		T = np.array([self.temperatures[0],(self.temperatures[0]+self.temperatures[-1])/2,self.temperatures[-1]])
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

