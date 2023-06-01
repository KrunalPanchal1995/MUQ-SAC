import numpy as np
import re,os
import uncertainty
import Input_file_reader
import matplotlib.pyplot as plt
from collections import OrderedDict

class MechanismManipulator():
	def __init__(self, filePath_mech,fileType,filePath_therm,filePath_trans,reactionList,focList,thirdBodyList,thermoList,transportList, rUnsrt,fUnsrt,mUnsrt,thUnsrt,trUnsrt):
		#print("Hello")
		self.chemtag = {}
		self.chemtagCollision = {}
		self.fileType = fileType
		self.ReactionList = reactionList
		self.molecule = []
		self.mList = thirdBodyList
		self.focList = focList
		self.thermoList = thermoList
		self.transportList = transportList
		#print(thermoList)
		#parsing third body parameters
		for i in thirdBodyList:
			self.molecule.append(mUnsrt[i].branches.split(","))
		#print(self.molecule)
		
		self.rxnUnsrt = rUnsrt
		self.focUnsrt = fUnsrt
		self.tbdUnsrt = mUnsrt
		self.thermoUnsrt = thUnsrt
		self.transUnsrt	= trUnsrt

		self.Test_string = []
		self.kineticDATA = Input_file_reader.MechParsing(filePath_mech)
		self.file =   filePath_mech
		wd = filePath_mech.split("/")[:-1]
		#self.FileName = filePath.split("/")[-1]
		s = "/"
		self.directory = s.join(wd)
		self.parentFile_mech =  open(filePath_mech).read()
		self.parentFile_thermo =  open(filePath_therm).read()
		self.parentFile_transport =  open(filePath_trans).read()
		f_mech = open(filePath_mech).readlines()
		f_therm = open(filePath_therm).readlines()
		f_trans = open(filePath_trans).readlines()
		mechLine = ""
		thermLine = ""
		transLine = ""
		Line = ""
	
		temp = []
		for i in f_mech:
			if i.startswith("!"):
				continue
			else:
				temp.append(i)
		f_mech = temp
		#print(f_mech)
		self.PerturbingReactions = {}
		self.PerturbingFoC = {}
		self.PerturbingCollisionEff = {}
		self.PerturbingThermo = {}
		self.PerturbingTransport = {}
		#print(self.ReactionList)
		#print(self.transportList)
		#print(self.fileType)
		
		if "FlameMaster" in fileType:
			for r in self.ReactionList:
				for i in f_mech:
					temp = r.split("-")[0]	
					if temp+":" == i.split(" ")[0]: #logic to extract reaction index from the string
						if  "}" in i:
							mechLine += i
							self.PerturbingReactions[r] = mechLine
							mechLine = ""
						else:
							j = f_mech.index(i)
							while "}" not in mechLine:
								mechLine += f_mech[j]
								j = j+1
							#print(mechLine)
							self.PerturbingReactions[r] = mechLine
							mechLine = ""  
			#print("FoC list: \t{}\n\n".format(self.focList))
			for foc in self.focList:
				for i in f_mech:
					temp = foc.split("-")[0]
					if temp+":" == i.split(" ")[0]:
							j = f_mech.index(i)
							while "}" not in Line:
								Line += f_mech[j]
								j = j+1
							self.PerturbingFoC[foc] = Line
							Line = "" 
					
			for m in self.mList:
				for i in f_mech:
					temp = m.split("-")[0]
					if temp in i.split(" "):
						if  "." in i:
							Line += i
							self.PerturbingCollisionEff[m] = Line
							Line = ""
						else:
							j = f.index(i)
							while "." not in Line:
								Line += f[j]
								j = j+1
							self.PerturbingCollisionEff[m] = Line
							Line = "" 
			 
		
		elif "chemkin" in fileType :
			for r in self.ReactionList:
				for i in f_mech:
					if r in i: 
						#print(i)
						if "DUPLICATE" in f_mech[f_mech.index(i)+1]:
							self.chemtag[r] = 1
							start = f_mech.index(i)
							for j in range(4):
								mechLine +=f_mech[start+j]
							#print(mechLine)
							self.PerturbingReactions[r] = mechLine
							#print(mechLine)
							mechLine = "" 
							break
						
						elif "LOW" not in f_mech[f_mech.index(i)+1] and "/" in f_mech[f_mech.index(i)+2] and "LOW" not in f_mech[f_mech.index(i)+2] and "TROE" not in f_mech[f_mech.index(i)+2]:
							self.chemtag[r] = 2
							start = f_mech.index(i)  
							for j in range(2):
								mechLine +=f_mech[start+j]
							#print(mechLine)
							self.PerturbingReactions[r] = mechLine
							mechLine = "" 
						
						elif "LOW" in f_mech[f_mech.index(i)+1] and "TROE" in f_mech[f_mech.index(i)+2]and "/" not in f_mech[f_mech.index(i)+3]:
							self.chemtag[r] = 3
							start = f_mech.index(i)  
							for j in range(3):
								mechLine +=f_mech[start+j]
							#print(mechLine)
							self.PerturbingReactions[r] = mechLine
							mechLine = "" 
						
						elif "LOW" in f_mech[f_mech.index(i)+1] and "TROE" in f_mech[f_mech.index(i)+2] and "/" in f_mech[f_mech.index(i)+3]:
							self.chemtag[r] = 4
							start = f_mech.index(i)  
							for j in range(4):
								mechLine += f_mech[start+j]
							#print(mechLine)
							self.PerturbingReactions[r] = mechLine
							mechLine = "" 
						else:
							self.chemtag[r] = 0
							mechLine += i
							self.PerturbingReactions[r] = mechLine
							mechLine = ""
			for foc in self.focList:
				for i in f_mech:
					if foc in i:
						if "LOW" in f_mech[f_mech.index(i)+1] and "TROE" in f_mech[f_mech.index(i)+2]and "/" not in f_mech[f_mech.index(i)+3]:
							self.chemtag[r] = 3
							start = f_mech.index(i)  
							for j in range(3):
								mechLine +=f_mech[start+j]
							self.PerturbingFoC[foc] = mechLine
							mechLine = ""
						
						elif "LOW" in f_mech[f_mech.index(i)+1] and "TROE" in f_mech[f_mech.index(i)+2] and "/" in f_mech[f_mech.index(i)+3]:
							self.chemtag[r] = 4
							start = f_mech.index(i)  
							for j in range(4):
								mechLine += f_mech[start+j]
							self.PerturbingFoC[foc] = mechLine
							mechLine = ""
						else:
							 print("The input reaction is not a TROE type pressure dependent reaction")	
		
			for m in self.mList:
				for i in f_mech:
					if m in i:
						if "LOW" not in f_mech[f_mech.index(i)+1] and "/" in f_mech[f_mech.index(i)+2]:
							self.chemtagCollision[m] = 2
							start = f_mech.index(i)  
							for j in range(2):
								mechLine +=f_mech[start+j]
							self.PerturbingCollisionEff[m] = mechLine
							mechLine = ""
						
						elif "LOW" in f_mech[f_mech.index(i)+1] and "TROE" in f_mech[f_mech.index(i)+2] and "/" in f_mech[f_mech.index(i)+3]:
							self.chemtagCollision[m] = 4
							start = f_mech.index(i)  
							for j in range(4):
								mechLine += f_mech[start+j]
							self.PerturbingCollisionEff[m] = mechLine
							mechLine = ""
						else:
							print("The input reaction is not a third body reaction")  
							break
		#print(self.PerturbingReactions)	
		#print(self.PerturbingFoC)
		#print(self.PerturbingCollisionEff)
		for th in self.thermoList:
			for i in f_therm:
				#print(th.split("-")[0])
				#print(i.split(" "))
				if "!" not in i:
					temp = th.split("-")[0]
					
					if temp in i.split(" "):
						j =f_therm.index(i)
						count = 0
						while count<4:
							Line+=f_therm[j+count]
							Line+="\n"
							count+=1
						self.PerturbingThermo[th] = Line
						Line = ""
		for ts in self.transportList:
			for i in f_trans:
				#print(ts)
				temp = ts.split("-")[0]
				if temp in i.split(" "):
					if "!" not in i:
						Line += i
						self.PerturbingTransport[ts] = Line
						Line = ""
		#print(self.PerturbingThermo)
		#print("Perturnbing list for FoC:\t{}\n\n".format(self.PerturbingFoC))
# if doing sensitivity analysis, sens_AL = 0, else sens_AL = 1
		#print(self.chemtag)
	def GeneratePerturbedMechanism(self,target,beta,Sens_AL,isOpt):
		NewMechanismString = self.parentFile_mech
		NewThermoString = self.parentFile_thermo
		NewTransportString = self.parentFile_transport
		self.RBeta_Dict = {}
		#print(NewMechanismString)
		#Beta is an array
		beta_group = []
		
		######################
		#    REACTION        #
		######################
		
		if len(self.ReactionList) != 0:
			#print("rxn founds")
			beta_rxn = beta[0:3*len(self.ReactionList)]
			beta_per_rxn = []
			count_rxn = 0
			for i in beta_rxn:
				temp = []
				if count_rxn < 3*len(self.ReactionList):
					for i in range(3):
						temp.append(beta_rxn[count_rxn])
						count_rxn+=1
					beta_per_rxn.append(temp)
			for index, rxn in enumerate(self.ReactionList):
				self.RBeta_Dict[rxn] = beta_per_rxn[index]
			
			for index,rxn in enumerate(self.PerturbingReactions):
				#print(self.PerturbingReactions[rxn])
				NewReaction = self.ModifyReaction(self.PerturbingReactions[rxn],target,np.asarray(beta_per_rxn[index]),rxn,Sens_AL,isOpt)
				
				#print(NewReaction)
				NewMechanismString = NewMechanismString.replace(self.PerturbingReactions[rxn],NewReaction)		
		else:
			#print("rxn not founds")
			count_rxn = 0
			NewMechanismString = NewMechanismString
		
		#print(beta_rxn)
		######################
		#    FallOffCurve    #
		######################
		
		if len(self.focList) != 0:
			#print("foc founds")
			beta_foc = beta[3*len(self.ReactionList):3*len(self.ReactionList)+len(self.focList)]
			
			beta_per_foc = []
			count_foc = 0
			for i in beta_foc:
				temp = []
				if count_foc < len(self.focList):
					for i in range(1):
						temp.append(beta_foc[count_foc])
						count_foc+=1
					beta_per_foc.append(temp)
			
			for index, foc in enumerate(self.focList):
				self.RBeta_Dict[foc] = beta_per_foc[index]	
			for index,foc in enumerate(self.PerturbingFoC):
				#print("Entered the loop\n\n")
				NewFoC = self.ModifyReaction(self.PerturbingFoC[foc],target,np.asarray(beta_per_foc[index]),foc,Sens_AL,isOpt)
				NewMechanismString = NewMechanismString.replace(self.PerturbingFoC[foc],NewFoC)
			#print("Exited the loop\n\n")	
		else:
			#print("foc not founds")
			count_foc = 0 
			NewMechanismString = NewMechanismString
		
		###############################
		#    Collision Efficiencies   #
		###############################

		if len(self.molecule) != 0:
			#print("m founds")
			beta_m = beta[3*len(self.ReactionList)+len(self.focList):3*len(self.ReactionList)+len(self.focList)+len(np.asarray(self.molecule).flatten())]
			count_m = 0	
			for i in self.molecule:
				temp = {}
				if count_m < len(np.asarray(self.molecule).flatten()):
					for j in range(len(i)):
						temp[str(i[j])] = beta_m[count_m]
						count_m+=1
					beta_group.append(temp)
			#print(beta_group)
			for index, m in enumerate(self.mList):
				self.RBeta_Dict[m] = beta_group[index]
			for index, rxn in enumerate(self.PerturbingCollisionEff):
				NewLine = self.ModifyLine(self.PerturbingCollisionEff[rxn],np.asarray(beta_group[index]),rxn)
				NewMechanismString = NewMechanismString.replace(self.PerturbingCollisionEff[rxn],NewLine)
		else:
			#print("m not founds")
			count_m = 0
			NewMechanismString = NewMechanismString
		###############################
		#    Thermo Parameters        #
		###############################
		
		if 	len(self.thermoList) != 0:	
			#print("th founds")		
			beta_th = beta[3*len(self.ReactionList)+len(self.focList)+len(np.asarray(self.molecule).flatten()):3*len(self.ReactionList)+len(self.focList)+len(np.asarray(self.molecule).flatten())+7*len(self.thermoList)]
			beta_thermo = []
			#print(beta_th)
			count_th = 0
			for i in self.thermoList:
				temp = {}
				temp["Hcp"] = []
				if count_th < 7*len(self.thermoList):
					for j in range(5):
						temp["Hcp"].append(beta_th[count_th])
						count_th+=1
					temp["h"] = beta_th[count_th]
					count_th+=1
					temp["e"] = beta_th[count_th]
					count_th+=1
					beta_thermo.append(temp)
			
			
			for index, th in enumerate(self.thermoList):
				self.RBeta_Dict[th] = beta_thermo[index]
			for index,th in enumerate(self.PerturbingThermo):
				#print(self.PerturbingThermo[th])
				NewThermo = self.ModifyThermo(self.PerturbingThermo[th],np.asarray(beta_thermo[index]),th,Sens_AL,isOpt)
				
				NewThermoString = NewThermoString.replace(self.PerturbingThermo[th],NewThermo)
		else:
			#print("tn not founds")
			count_th = 0
			NewThermoString = NewThermoString
		###############################
		#    Transport parameters     #
		###############################		
		if len(self.transportList) != 0:	
			#print("ts founds")
			beta_ts = beta[3*len(self.ReactionList)+len(self.focList)+len(np.asarray(self.molecule).flatten())+7*len(self.thermoList):3*len(self.ReactionList)+len(self.focList)+len(np.asarray(self.molecule).flatten())+7*len(self.thermoList)+2*len(self.transportList)]
			#print(len(beta))
			beta_transport = []
			#print(beta_th)
			count_ts = 0
			for i in self.transportList:
				temp = {}
				if count_ts < 2*len(self.transportList):
					temp["LJe"] = beta_ts[count_ts]
					count_ts+=1
					temp["LJs"]=beta_ts[count_ts]
					count_ts+=1
					beta_transport.append(temp)
			
			for index, ts in enumerate(self.transportList):
				self.RBeta_Dict[ts] = beta_transport[index]
			for index,ts in enumerate(self.PerturbingTransport):
				NewTransport = self.ModifyTransport(self.PerturbingTransport[ts],np.asarray(beta_ts[index]),ts,Sens_AL,isOpt)
				NewTransportString = NewTransportString.replace(self.PerturbingTransport[ts],NewTransport)
		
		else:
			#print("ts not founds")
			count_ts = 0
			NewTransportString = NewTransportString
				
		if isOpt == "True":	
			mechFile = open("mechanism_opt.mech",'w')
			mechFile.write(NewMechanismString)
			mechFile.close()
		
			thermoFile = open("thermo_opt.therm",'w')
			thermoFile.write(NewThermoString)
			thermoFile.close()
		
			transportFile = open("transport_opt.trans",'w')
			transportFile.write(NewTransportString)
			transportFile.close()
		else:
			mechFile = open("mechanism.mech",'w')
			mechFile.write(NewMechanismString)
			mechFile.close()
	
			thermoFile = open("thermo.therm",'w')
			thermoFile.write(NewThermoString)
			thermoFile.close()
	
			transportFile = open("transport.trans",'w')
			transportFile.write(NewTransportString)
			transportFile.close()
		return self.RBeta_Dict
	
	def getOptBranchingData(self):
		return self.BranchingDict_Opt
	
	def getUnOptBranchingData(self):
		return self.BranchingDict_unOpt
	
	
	##Error part
	def getBranchingParameters(self,Branching_rxns,A_list,n_list,Ea_list,isOpt,R):
		#print("List after perturbation for reactions {} = \t\t{} \n\t\t\t\t{} \n\t\t\t\t{}\n".format(Branching_rxns,A_list,n_list,Ea_list))
		
		kDash = []
		A = []
		n = []
		Ea = []
		self.BranchingDict_unOpt = {}
		self.BranchingDict_Opt = {}
		T_ = self.rxnUnsrt[str(Branching_rxns[0])].temperatures
		#T_ = np.linspace(500,2500,50)
		self.ratio = []
		#print(A_list)
		
		#Find th ration for all the reactions
		for j in T_:
			k = 0
			k_temp = []
			r_temp = []
			
			#Find the rate and total rate simultaneously for all reactions
			for i,rxn in enumerate(Branching_rxns):
			    rate = A_list[i]*(j**n_list[i])*(np.exp(-Ea_list[i]/(j*float(R))))
			    #print(Ea_list[i])
			    #print(Ea_list[i]/(j*8.314E-3))
			    #print(np.exp(-Ea_list[i]/(j*8.314E-3)))
			    #total rate
			    
			    k+=rate
			    k_temp.append(rate)
			
			k_Dash = []
			for i in k_temp:
			    r = i/k
			    #print("ratio for reactions {} = {}".format(Branching_rxns,r))
			    r_temp.append(r)
			    k_Dash.append(i)
			kDash.append(k_Dash)
			self.ratio.append(r_temp)
		
		if isOpt =="True":
			self.BranchingDict_Opt[str(Branching_rxns)] = self.ratio
		
		else:
			self.BranchingDict_unOpt[str(Branching_rxns)] = self.ratio
		
		#print(kDash)
		for i in range(len(Branching_rxns)):
			A_,n_,Ea_ = Input_file_reader.Curvefit(np.asarray(kDash)[:,i],T_).getBranchingCurveFit()
			A.append(A_)
			n.append(n_)
			Ea.append(Ea_)
		return A,n,Ea
	
	def getPerturbedParameters(self,Branching_rxns,beta,isOpt,R):
		aDash = []
		cholesky_matrix = []
		
		beta = np.asarray(beta)
		#print(beta)
		#print(Branching_rxns)
		A_list = []
		n_list = []
		Ea_list = []
		A_list_ = []
		n_list_ = []
		Ea_list_ = []
		cholesky_matric = []
		LZeta = []
		Arrhenius_list = []
		ArrMax= []
		ArrMin = []
		#print("Beta si {}".format(beta))
		#print("\n\nBeta for branching rxns!!\n\n")
		for i,rxn in enumerate(Branching_rxns):
			a,n,e = Input_file_reader.MechParsing(self.file).getRxnData(rxn)
			#print("Data from Input_file_reader for {}-{}\t{}\t{}\n".format(rxn,a,n,e))
			A_list.append(float(a))
			n_list.append(float(n))
			Ea_list.append(float(e))
			alpha = np.log(a)
			#print("\n")
			#print(a)
			#print(alpha)
			#print("\n")
			epsilon = e/float(R)
			
			if rxn in self.ReactionList:
				cMatrix = self.rxnUnsrt[rxn].cholskyDeCorrelateMat
				cholesky_matrix.append(cMatrix)
				zeta = self.rxnUnsrt[rxn].zeta.x
				LZeta.append(np.asarray(np.dot(cMatrix,beta[i]*zeta).flatten()))
				#print("\nLZETA = {}\n".format(LZeta))
			
				logArrheniusMin = np.array([alpha,n,epsilon])  - np.array(np.dot(cMatrix,zeta)).flatten()
				logArrheniusMax = np.array([alpha,n,epsilon])  + np.array(np.dot(cMatrix,zeta)).flatten()
				ArrheniusMax = [np.exp(logArrheniusMax[0]), logArrheniusMax[1], float(R)*float(logArrheniusMax[2])]
				ArrheniusMin = [np.exp(logArrheniusMin[0]), logArrheniusMin[1], float(R)*float(logArrheniusMin[2])]

			
				log_Arrhenius = np.array([alpha,n,epsilon]) + np.array(np.dot(cMatrix,beta[i]*zeta)).flatten()

			
				Arrhenius = [np.exp(log_Arrhenius[0]), log_Arrhenius[1], float(R)*float(log_Arrhenius[2])]
				#print("Nominal = {}".format(np.array([alpha,n,epsilon])))
				#print("LZeta = {}".format(np.array(np.dot(cMatrix,beta[i]*zeta)).flatten()))
				#print(Arrhenius)
				A_list_.append(float(Arrhenius[0]))
				n_list_.append(float(Arrhenius[1]))
				Ea_list_.append(float(Arrhenius[2]))
				LZeta.append(np.dot(cMatrix,beta[i]*zeta).flatten())
				Arrhenius_list.append(Arrhenius)
				ArrMax.append(ArrheniusMax)
				ArrMin.append(ArrheniusMin)
				continue
				
			else: 
				cMatrix = np.identity(3,dtype=float)
				zeta = beta[i]
				cholesky_matrix.append(cMatrix)
				LZeta.append(np.asarray(np.dot(cMatrix,beta[i]*zeta).flatten()))
				#print("\nLZETA = {}\n".format(LZeta))
			
				logArrheniusMin = np.array([alpha,n,epsilon]  - np.dot(cMatrix,zeta)).flatten()
				logArrheniusMax = np.array([alpha,n,epsilon]  + np.dot(cMatrix,zeta)).flatten()
				ArrheniusMax = [np.exp(logArrheniusMax[0]), logArrheniusMax[1], float(R)*float(logArrheniusMax[2])]
				ArrheniusMin = [np.exp(logArrheniusMin[0]), logArrheniusMin[1], float(R)*float(logArrheniusMin[2])]

				log_Arrhenius = np.array([alpha,n,epsilon]) + np.array(np.dot(cMatrix,beta[i]*zeta)).flatten()
			
				Arrhenius = [np.exp(log_Arrhenius[0]), log_Arrhenius[1], float(R)*float(log_Arrhenius[2])]
				A_list_.append(float(Arrhenius[0]))
				n_list_.append(float(Arrhenius[1]))
				Ea_list_.append(float(Arrhenius[2]))
				LZeta.append(np.dot(cMatrix,beta[i]*zeta).flatten())
				Arrhenius_list.append(Arrhenius)
				ArrMax.append(ArrheniusMax)
				ArrMin.append(ArrheniusMin)
				continue
		
		
		#print("Parameters before perturbation for {} = {} \t\t\n{}\t\t\t\t{}\n".format(Branching_rxns,A_list,n_list,Ea_list))
		#print("Parameters after perturbation for {} = {} \t\t\n{}\t\t\t\t{}\n".format(Branching_rxns,A_list_,n_list_,Ea_list_))
		#print("Branching reactions :\t{}".format(Branching_rxns))
		A_o,n_o,Ea_o = self.getBranchingParameters(Branching_rxns,A_list,n_list,Ea_list,isOpt,R)
		#print("Original Branching parameters:\n\nA:\t{}\n n:\t{}\nEa:\t{}\n".format(A_o,n_o,Ea_o))
		
		A_p,n_p,Ea_p = self.getBranchingParameters(Branching_rxns,A_list_,n_list_,Ea_list_,isOpt,R)
		#print("Perturbed Branching parameters:\n\nA:\t{}\n n:\t{}\nEa:\t{}\n".format(A_p,n_p,Ea_p))
		#print("\n params = {}\t{}\t{} \n".format(A_p,n_p,Ea_p))
		return A_p,n_p,Ea_p
	
	
	def getNewBranchRxn(self,reaction,branchRxn,index,beta,isOpt):
		if self.fileType == "FlameMaster":
			w = re.compile(r'.*:?.*->.*A\s*?=\s*?(\S*E\S*)?\s*n\s*?=\s*(\S*)?\s*E\s*?=\s*(\S*)?\s*.*\}}?'.format(), re.DOTALL | re.IGNORECASE)
			match1 = re.search(w, reaction)
			aString = match1.group(1)
			a = float(aString)
			alpha = np.log(a)

			match1 = re.search(w, reaction)
			nString = match1.group(2)
			n = float(nString)


			match1 = re.search(w, reaction)
			eString = match1.group(3)
			e = float(eString)
			epsilon = e/8.314E-3
			#if target.temperature<2:
			#	T = 1000/target.temperature
			#else:
			#	T = target.temperature
			#theta = np.array([1,np.log(T),-1/T])
			rxn_list = []
			beta_list = []
			#print("Beta for branching rxns:\n")
			for i in branchRxn:
				rxn_list.append(i)
				if i in self.ReactionList:
					#print("\t{}:\t{}\n".format(i,self.RBeta_Dict[i]))
					beta_list.append(np.asarray(self.RBeta_Dict[i]))
					continue
				else:
					beta_list.append(np.array([0.0,0.0,0.0]))
					continue
					#print("\t{}:\t{}\n".format(i,np.array([0.0,0.0,0.0])))
			#print(beta_list)
			
			#print(rxn_list)
			A_perturbed,n_perturbed,Ea_perturbed = self.getPerturbedParameters(rxn_list,beta_list,isOpt,8.314E-3)
		
			if float(A_perturbed[0]) < 0:
				print("Anomaly Found!! In {}, mechanism has negative Arrheenius parameters".format(os.getcwd()))
			reaction = reaction.replace(aString,'{:.3E}'.format(A_perturbed[0]))
			reaction = reaction.replace(nString,'{:.3f}'.format(n_perturbed[0]))
			reaction = reaction.replace(eString,'{:.3f}'.format(Ea_perturbed[0]))
		elif self.fileType == "chemkin":
			
			if self.chemtag[index] == 0:
				w = re.compile(r'(.*?\<\=\>.*? )\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+)).*?\n?',re.DOTALL | re.IGNORECASE)
				match1 = re.search(w, reaction)
				aString = match1.group(2)
				a = float(aString)
				alpha = np.log(a)

				nString = match1.group(5)
				n = float(nString)

				eString = match1.group(8)
				e = float(eString)
				epsilon = e/1.987
								
			elif self.chemtag[index] == 1:
			
				w = re.compile(r'(.*?\<\=\>.*?) \s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?DUPLICATE\s*?\n?(.*?\<\=\>.*?) \s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?DUPLICATE\s*?\n?',re.DOTALL|re.IGNORECASE)

			#print("\n")
			
				match1 = re.search(w, reaction)
				#print(match1.group())
				aString = match1.group(2)
				a = float(aString)
				alpha = np.log(a)

				match1 = re.search(w, reaction)
				nString = match1.group(5)
				n = float(nString)


				match1 = re.search(w, reaction)
				eString = match1.group(8)
				e = float(eString)
				epsilon = e/1.987
				
				match1 = re.search(w, reaction)
				#print(match1.group())
				aString_dup = match1.group(12)
				a_dup = float(aString_dup)
				alpha_dup = np.log(a_dup)

				match1 = re.search(w, reaction)
				nString_dup = match1.group(15)
				n = float(nString_dup)


				match1 = re.search(w, reaction)
				eString_dup = match1.group(18)
				e_dup = float(eString_dup)
				epsilon = e_dup/1.987
			#if target.temperature<2:
			#	T = 1000/target.temperature
			#else:
			#	T = target.temperature
			#theta = np.array([1,np.log(T),-1/T])
			rxn_list = []
			beta_list = []
			#print("Beta for branching rxns:\n")
			
			for i in branchRxn:
				rxn_list.append(i)
				if i in self.ReactionList:
						
					#print("\t{}:\t{}\n".format(i,self.RBeta_Dict[i]))
					beta_list.append(np.asarray(self.RBeta_Dict[i]))
				
				elif i not in self.ReactionList:
					
					beta_list.append(np.array([0.0,0.0,0.0]))
					#print("\t{}:\t{}\n".format(i,np.array([0.0,0.0,0.0])))
			
			A_perturbed,n_perturbed,Ea_perturbed = self.getPerturbedParameters(rxn_list,beta_list,isOpt,1.987)
		
			if float(A_perturbed[0]) < 0:
				print("Anomaly Found!! In {}, mechanism has negative Arrheenius parameters".format(os.getcwd()))
			#print(reaction)
			reaction = reaction.replace(aString,'{:.3E}'.format(A_perturbed[0]))
			reaction = reaction.replace(nString,'{:.3f}'.format(n_perturbed[0]))
			reaction = reaction.replace(eString,'{:.3E}'.format(Ea_perturbed[0]))
			#print(reaction)
		else:
			print("Invalid kinetic mechanism file")
			
		return reaction
		
		
	def getNewUpperLimits(self,reaction,rxn,branchRxn,beta,isOpt,delta):
		if delta == 1:
			if self.fileType == "FlameMaster":
				w = re.compile(r'(\d*\w):\s*(.*)\s*->\s*(.*)\s*?\{\s*?Ai\s*?=\s*?(\S*E\S*)?\s*ni\s*?=\s*(\S*)?\s*Ei\s*?=\s*(\S*)?\s*A\s*?=\s*?(\S*E\S*)?\s*n\s*?=\s*(\S*)?\s*E\s*?=\s*(\S*)?\s*fca\s*?=\s*?(\S*E\S*)?\s*?fcta\s*?=\s*?(\S*E\S*)?\s*?fcb\s*?=\s*?(\S*E\S*)?\s*?fctb\s*?=\s*?(\S*E\S*)?\s*?fcc\s*?=\s*?(\S*)?\s*?fctc\s*?=\s*?(\S*E\S*).*', re.DOTALL | re.IGNORECASE)
				match1 = re.search(w, reaction)
				aiString = match1.group(4)
				ai = float(aiString)
				alpha_i = np.log(ai)

				niString = match1.group(5)
				ni = float(niString)


				eiString = match1.group(6)
				ei = float(eiString)
				epsilon_i = ei/8.314E-3
	
			
				choleskyMatrix = self.rxnUnsrt[rxn].cholskyDeCorrelateMat
				zeta = self.rxnUnsrt[rxn].zeta.x
		
				beta = self.RBeta_Dict[rxn]
				#print("\tBeta (HPL):\t{}\n".format(beta_i))
				log_Arrhenius_i = np.array([alpha_i,ni,epsilon_i]) + np.dot(choleskyMatrix,beta*zeta).flatten()
				#need to convert back to original values
				Arrhenius_i = [np.exp(log_Arrhenius_i[0]), log_Arrhenius_i[1], 8.314E-3*log_Arrhenius_i[2]]
	
				reaction = reaction.replace(aiString,'{:.3E}'.format(Arrhenius_i[0]))
				reaction = reaction.replace(niString,'{:.3f}'.format(Arrhenius_i[1]))
				reaction = reaction.replace(eiString,'{:.3f}'.format(Arrhenius_i[2]))
			elif self.fileType == "chemkin":
				if self.chemtag[rxn] == 2:
				
					w = re.compile(r'(.*?\<\=\>.*?) \s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?(.*?\/\s*?.*)\/\s*?.*?\n?',re.DOTALL|re.IGNORECASE)

				elif self.chemtag[rxn] == 3:
				
					w = re.compile(r'(.*?\<\=\>.*?) \s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?(LOW)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?\/\s*?.*?\n?\s*?(TROE)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?\/\s*?\n?',re.DOTALL|re.IGNORECASE)
			
				elif self.chemtag[rxn] == 4:					
					w = re.compile(r'(.*?\<\=\>.*?) \s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?(LOW)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?\/\s*?.*?\n?\s*?(TROE)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?\/\s*?\n(.*)',re.DOTALL|re.IGNORECASE)

			
				match1 = re.search(w, reaction)
	
				aString = match1.group(2)
				a = float(aString)
				alpha = np.log(a)

				nString = match1.group(5)
				n = float(nString)

				eString = match1.group(8)
				e = float(eString)
				epsilon = e/1.987
	
				choleskyMatrix = self.rxnUnsrt[rxn].cholskyDeCorrelateMat
				zeta = self.rxnUnsrt[rxn].zeta.x
		
				beta = self.RBeta_Dict[rxn]
				#print("\tBeta (LPL):\t{}\n".format(beta))
				self.nominal_lowPressureLimit = np.array([alpha,n,epsilon])
				log_Arrhenius = np.array([alpha,n,epsilon]) + np.dot(choleskyMatrix,beta*zeta).flatten()
				#need to convert back to original values
				self.Arrhenius_lowPressureLimit = [np.exp(log_Arrhenius[0]), log_Arrhenius[1], 1.987*log_Arrhenius[2]]
				reaction = reaction.replace(aString,'{:.3E}'.format(self.Arrhenius_lowPressureLimit[0]))
				reaction = reaction.replace(nString,'{:.3f}'.format(self.Arrhenius_lowPressureLimit[1]))
				reaction = reaction.replace(eString,'{:.3E}'.format(self.Arrhenius_lowPressureLimit[2]))
			else:
				print("Input kinetic mechansim file is invalid")
		else:
			reaction = reaction
		return reaction
	
	def getNewLowerLimits(self,reaction,rxn,branchRxn,beta,isOpt):
		if delta == 1:
			if self.fileType == "FlameMaster":
				w = re.compile(r'(\d*\w):\s*(.*)\s*->\s*(.*)\s*?\{\s*?Ai\s*?=\s*?(\S*E\S*)?\s*ni\s*?=\s*(\S*)?\s*Ei\s*?=\s*(\S*)?\s*A\s*?=\s*?(\S*E\S*)?\s*n\s*?=\s*(\S*)?\s*E\s*?=\s*(\S*)?\s*fca\s*?=\s*?(\S*E\S*)?\s*?fcta\s*?=\s*?(\S*E\S*)?\s*?fcb\s*?=\s*?(\S*E\S*)?\s*?fctb\s*?=\s*?(\S*E\S*)?\s*?fcc\s*?=\s*?(\S*)?\s*?fctc\s*?=\s*?(\S*E\S*).*', re.DOTALL | re.IGNORECASE)
				match1 = re.search(w, reaction)
	
				aString = match1.group(7)
				a = float(aString)
				alpha = np.log(a)

				nString = match1.group(8)
				n = float(nString)

				eString = match1.group(9)
				e = float(eString)
				epsilon = e/8.314E-3
	
				choleskyMatrix = self.rxnUnsrt[rxn].cholskyDeCorrelateMat
				zeta = self.rxnUnsrt[rxn].zeta.x
		
				beta = self.RBeta_Dict[rxn]
				#print("\tBeta (LPL):\t{}\n".format(beta))
				self.nominal_lowPressureLimit = np.array([alpha,n,epsilon])
				log_Arrhenius = np.array([alpha,n,epsilon]) + np.dot(choleskyMatrix,beta*zeta).flatten()
				#need to convert back to original values
				self.Arrhenius_lowPressureLimit = [np.exp(log_Arrhenius[0]), log_Arrhenius[1], 8.314E-3*log_Arrhenius[2]]
				reaction = reaction.replace(aString,'{:.3E}'.format(self.Arrhenius_lowPressureLimit[0]))
				reaction = reaction.replace(nString,'{:.3f}'.format(self.Arrhenius_lowPressureLimit[1]))
				reaction = reaction.replace(eString,'{:.3f}'.format(self.Arrhenius_lowPressureLimit[2]))
			elif self.fileType == "chemkin":
			
				if self.chemtag[rxn] == 3:
					w = re.compile(r'(.*?\<\=\>.*?) \s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?(LOW)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?\/\s*?.*?\n?\s*?(TROE)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?\/\s*?\n?',re.DOTALL|re.IGNORECASE)
			
				elif self.chemtag[rxn] == 4:					
					w = re.compile(r'(.*?\<\=\>.*?) \s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?(LOW)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?\/\s*?.*?\n?\s*?(TROE)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?\/\s*?\n(.*)',re.DOTALL|re.IGNORECASE)
			
				match1 = re.search(w, reaction)
	
				aString = match1.group(12)
				a = float(aString)
				alpha = np.log(a)

				nString = match1.group(15)
				n = float(nString)

				eString = match1.group(18)
				e = float(eString)
				epsilon = e/1.987
	
				choleskyMatrix = self.rxnUnsrt[rxn].cholskyDeCorrelateMat
				zeta = self.rxnUnsrt[rxn].zeta.x
		
				beta = self.RBeta_Dict[rxn]
				#print("\tBeta (LPL):\t{}\n".format(beta))
				self.nominal_lowPressureLimit = np.array([alpha,n,epsilon])
				log_Arrhenius = np.array([alpha,n,epsilon]) + np.dot(choleskyMatrix,beta*zeta).flatten()
				#need to convert back to original values
				self.Arrhenius_lowPressureLimit = [np.exp(log_Arrhenius[0]), log_Arrhenius[1], 1.987*log_Arrhenius[2]]
				reaction = reaction.replace(aString,'{:.3E}'.format(self.Arrhenius_lowPressureLimit[0]))
				reaction = reaction.replace(nString,'{:.3f}'.format(self.Arrhenius_lowPressureLimit[1]))
				reaction = reaction.replace(eString,'{:.3E}'.format(self.Arrhenius_lowPressureLimit[2]))
			else:
				print("Input kinetic mechansim file is invalid")
		else:
			reaction = reaction
		
		return reaction
	
	
	def getNewFallOffCurve(self,reaction,target,rxn,branchRxn,beta,isOpt):
		if delta == 1:
			if self.fileType == "FlameMaster":
				#print("Entered Foc subroutine!!")
				w = re.compile(r'(\d*\w):\s*(.*)\s*->\s*(.*)\s*?\{\s*?Ai\s*?=\s*?(\S*E\S*)?\s*ni\s*?=\s*(\S*)?\s*Ei\s*?=\s*(\S*)?\s*A\s*?=\s*?(\S*E\S*)?\s*n\s*?=\s*(\S*)?\s*E\s*?=\s*(\S*)?\s*fca\s*?=\s*?(\S*E\S*)?\s*?fcta\s*?=\s*?(\S*E\S*)?\s*?fcb\s*?=\s*?(\S*E\S*)?\s*?fctb\s*?=\s*?(\S*E\S*)?\s*?fcc\s*?=\s*?(\S*)?\s*?fctc\s*?=\s*?(\S*E\S*).*', re.DOTALL | re.IGNORECASE)
				T = target.temperature
				#print(T)
				match1 = re.search(w, reaction)
				fca_String = match1.group(10)
				fca = float(fca_String)
		
		
				fcta_String = match1.group(11)
				fcta = float(fcta_String)
		

				fcb_String = match1.group(12)
				fcb = float(fcb_String)
		
		
		
				fctb_String = match1.group(13)
				fctb = float(fctb_String)
		
		
		
				fcc_String = match1.group(14)
				fcc = float(fcc_String)
		
		
				fctc_String = match1.group(15)
				fctc = float(fctc_String)
		
		
				choleskyMatrix = self.rxnUnsrt[rxn].cholskyDeCorrelateMat
				zeta = self.rxnUnsrt[rxn].zeta
		
				beta = self.RBeta_Dict[rxn]
				#print("\tBeta (FoC):\t{}\n".format(beta))
				self.nominal_foc = np.array([fcb])
			
				self.perturbedFoC = np.array([fcb]) + np.array(np.dot(choleskyMatrix,beta*zeta)).flatten()
			
				#fcta_ = -T/np.log(self.perturbedFoC[0])
				fcb_ = self.perturbedFoC[0]
				#fca_ = self.perturbedFoC[0]
				fca_ = 1-fcb_
			
				reaction = reaction.replace(fca_String,'{:.3E}'.format(fca_))
				#reaction = reaction.replace(fcta_String,'{:.3E}'.format(fcta_))
				reaction = reaction.replace(fcb_String,'{:.3E}'.format(fcb_))
			
		
			elif self.fileType == "chemkin":
				
				if self.chemtag[rxn] == 3:
					w = re.compile(r'(.*?\<\=\>.*?) \s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?(LOW)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?\/\s*?.*?\n?\s*?(TROE)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?\/\s*?\n?',re.DOTALL|re.IGNORECASE)
								
				elif self.chemtag[rxn] == 4:					
					w = re.compile(r'(.*?\<\=\>.*?) \s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?(LOW)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?\/\s*?.*?\n?\s*?(TROE)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?\/\s*?\n(.*)',re.DOTALL|re.IGNORECASE)

				match1 = re.search(w, reaction)
				fca_String = match1.group(22)
				fca = float(fca_String)
				reaction = reaction.replace(fca_String,'{:.3E}'.format(fca_))
				choleskyMatrix = self.rxnUnsrt[rxn].cholskyDeCorrelateMat
				zeta = self.rxnUnsrt[rxn].zeta
		
				beta = self.RBeta_Dict[rxn]
				#print("\tBeta (FoC):\t{}\n".format(beta))
				self.nominal_foc = np.array([fca])
			
				self.perturbedFoC = self.nominal + np.array(np.dot(choleskyMatrix,beta*zeta)).flatten()
			
				fca_ = self.perturbedFoC[0]
				reaction = reaction.replace(fca_String,'{:.3f}'.format(fca_))
			else: 
				print("Invalid kinetic mechanism file")
		
		else:
			reaction = reaction
		return reaction
	
	
	
	def getNewRxn(self,reaction,rxn,beta):
		#print(self.fileType)
		#print(self.chemtag[rxn])
		if self.fileType == "FlameMaster":
			w = re.compile(r'.*:?.*->.*A\s*?=\s*?(\S*E\S*)?\s*n\s*?=\s*(\S*)?\s*E\s*?=\s*(\S*)?\s*.*\}}?'.format(), re.DOTALL | re.IGNORECASE)
			match1 = re.search(w, reaction)
			aString = match1.group(1)
			a = float(aString)
			alpha = np.log(a)

			match1 = re.search(w, reaction)
			nString = match1.group(2)
			n = float(nString)


			match1 = re.search(w, reaction)
			eString = match1.group(3)
			e = float(eString)
			epsilon = e/8.314E-3
		
			#print("\tBeta (p-Independent rnx):\t{}\n".format(beta))
		
			M = 3/(np.log(10.0));		
			choleskyMatrix = self.rxnUnsrt[rxn].cholskyDeCorrelateMat
			zeta = self.rxnUnsrt[rxn].zeta.x
			beta = self.RBeta_Dict[rxn]
		
			log_Arrhenius = np.array([alpha,n,epsilon]) + np.dot(choleskyMatrix,beta*zeta).flatten()
			#need to convert back to original values
			Arrhenius = [np.exp(log_Arrhenius[0]), log_Arrhenius[1], 8.314E-3*log_Arrhenius[2]]
		
			reaction = reaction.replace(aString,'{:.3E}'.format(Arrhenius[0]))
			reaction = reaction.replace(nString,'{:.3f}'.format(Arrhenius[1]))
			reaction = reaction.replace(eString,'{:.3E}'.format(Arrhenius[2]))
		
		elif self.fileType == "chemkin":
			
			if self.chemtag[rxn] == 0:
				w = re.compile(r'(.*?\<\=\>.*? )\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+)).*?\n?',re.DOTALL | re.IGNORECASE)
				match1 = re.search(w, reaction)
				aString = match1.group(2)
				a = float(aString)
				alpha = np.log(a)

				nString = match1.group(5)
				n = float(nString)

				eString = match1.group(8)
				e = float(eString)
				epsilon = e/1.987
				
			elif self.chemtag[rxn] == 1:
				
				w = re.compile(r'(.*?\<\=\>.*?) \s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?DUPLICATE\s*?\n?(.*?\<\=\>.*?) \s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?DUPLICATE\s*?\n?',re.DOTALL|re.IGNORECASE)
	
				match1 = re.search(w, reaction)
				aString = match1.group(2)
				a = float(aString)
				alpha = np.log(a)

				nString = match1.group(5)
				n = float(nString)

				eString = match1.group(8)
				e = float(eString)
				epsilon = e/1.987
		
			#print("\tBeta (p-Independent rnx):\t{}\n".format(beta))
		
			M = 3/(np.log(10.0));		
			choleskyMatrix = self.rxnUnsrt[rxn].cholskyDeCorrelateMat
			zeta = self.rxnUnsrt[rxn].zeta.x

		
			log_Arrhenius = np.array([alpha,n,epsilon]) + np.dot(choleskyMatrix,beta*zeta).flatten()
			#need to convert back to original values
			Arrhenius = [np.exp(log_Arrhenius[0]), log_Arrhenius[1], 1.987*log_Arrhenius[2]]
			
			#print(epsilon)
			#print(Arrhenius[2])
			#print('{:.3E}'.format(Arrhenius[2]))
		
			reaction = reaction.replace(aString,'{:.3E}'.format(Arrhenius[0]))
			reaction = reaction.replace(nString,'{:.3f}'.format(Arrhenius[1]))
			reaction = reaction.replace(eString,'{:.3E}'.format(Arrhenius[2]))
		else:
			print("Invalid kinetic mechanism file")
		return reaction	
	
	
	
	def getNewThermo_low(self,thermoParams,index,beta,isOpt,delta):
		if delta ==1:
			#searching the Lennard-Jonnes potential hard sphere diameter (\sigma_{ij})
			pattern = re.compile(r'\n*?([\w\d]*)\s*?(.*?G)\s*?(\d+\.?\d+)\s*?(\d+\.?\d+)\s*s*?(\d+\.?\d+)\s*?.*?1\s*?.*?\n?\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?2\n?\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?3\n?\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?4\n?',re.DOTALL|re.IGNORECASE)
			match = re.search(pattern,thermoParams)
			aString = match.group(13)
			a = float(aString)

			match = re.search(pattern,thermoParams)
			bString = match.group(14)
			b = float(bString)


			match = re.search(pattern,thermoParams)
			cString = match.group(15)
			c = float(cString)
			
			match = re.search(pattern,thermoParams)
			dString = match.group(16)
			d = float(dString)
			
			match = re.search(pattern,thermoParams)
			eString = match.group(17)
			e = float(eString)
			
			match = re.search(pattern,thermoParams)
			fString = match.group(18)
			f = float(fString)
			
			match = re.search(pattern,thermoParams)
			gString = match.group(19)
			g = float(gString)
			
			beta_cp = beta["Hcp"]
			#print(beta_cp)
			beta_h = beta["h"]
			beta_s = beta["e"]
			parameters = self.thermoUnsrt[index].branches
			if "Hcp" in parameters:
				choleskyMatrix_cp = self.thermoUnsrt[index].cholskyDeCorrelateMat["Hcp"]
				zeta_cp = self.thermoUnsrt[index].zeta["Hcp"]
				coeff_cp = np.array([a,b,c,d,e]) + np.array(np.dot(choleskyMatrix_cp,beta_cp*zeta_cp)).flatten()
			else:
				coeff_cp = np.array([a,b,c,d,e])
			if "h" in parameters:
				choleskyMatrix_H = self.thermoUnsrt[index].cholskyDeCorrelateMat["h"]
				zeta_H = self.thermoUnsrt[index].zeta["h"]
				coeff_H = np.array([f]) + np.array(np.dot(choleskyMatrix_H,beta_h*zeta_H)).flatten()		
			else:
				coeff_H = np.array([f])
			
			if "e" in parameters:
				choleskyMatrix_S = self.thermoUnsrt[index].cholskyDeCorrelateMat["e"]
				zeta_S = self.thermoUnsrt[index].zeta["e"]
				coeff_S = np.array([g]) + np.array(np.dot(choleskyMatrix_S,beta_s*zeta_S)).flatten()		
			else:
				coeff_S = np.array([g])
			
			
			#print("Before\n\n")
			#print(thermoParams)				
					
			thermoParams = thermoParams.replace(aString,'{:.8E}'.format(coeff_cp[0]))
			thermoParams = thermoParams.replace(bString,'{:.8E}'.format(coeff_cp[1]))
			thermoParams = thermoParams.replace(cString,'{:.8E}'.format(coeff_cp[2]))
			thermoParams = thermoParams.replace(dString,'{:.8E}'.format(coeff_cp[3]))
			thermoParams = thermoParams.replace(eString,'{:.8E}'.format(coeff_cp[4]))
			thermoParams = thermoParams.replace(fString,'{:.8E}'.format(coeff_H[0]))
			thermoParams = thermoParams.replace(gString,'{:.8E}'.format(coeff_S[0]))
			#print("\nAfter\n\n")
			#print(thermoParams)
		else:
			thermoParams = thermoParams
		return thermoParams	
		
	
	def getNewThermo_high(self,thermoParams,index,beta,isOpt,delta):
		if delta ==1:
			#searching the Lennard-Jonnes potential hard sphere diameter (\sigma_{ij})
			pattern = re.compile(r'\n*?([\w\d]*)\s*?(.*?G)\s*?(\d+\.?\d+)\s*?(\d+\.?\d+)\s*s*?(\d+\.?\d+)\s*?.*?1\s*?.*?\n?\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?2\n?\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?3\n?\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?4\n?',re.DOTALL|re.IGNORECASE)
			match = re.search(pattern,thermoParams)
			aString = match.group(6)
			a = float(aString)

			match = re.search(pattern,thermoParams)
			bString = match.group(7)
			b = float(bString)


			match = re.search(pattern,thermoParams)
			cString = match.group(8)
			c = float(cString)
			
			match = re.search(pattern,thermoParams)
			dString = match.group(9)
			d = float(dString)
			
			match = re.search(pattern,thermoParams)
			eString = match.group(10)
			e =float(eString)
			
			match = re.search(pattern,thermoParams)
			fString = match.group(11)
			f = float(fString)
			
			match = re.search(pattern,thermoParams)
			gString = match.group(12)
			g = float(gString)
			
			beta_cp = beta["Hcp"]
			#print(beta_cp)
			beta_h = beta["h"]
			beta_s = beta["e"]
			parameters = self.thermoUnsrt[index].branches
			if "Hcp" in parameters:
				choleskyMatrix_cp = self.thermoUnsrt[index].cholskyDeCorrelateMat["Hcp"]
				zeta_cp = self.thermoUnsrt[index].zeta["Hcp"]
				coeff_cp = np.array([a,b,c,d,e]) + np.array(np.dot(choleskyMatrix_cp,beta_cp*zeta_cp)).flatten()
			else:
				coeff_cp = np.array([a,b,c,d,e])
			if "h" in parameters:
				choleskyMatrix_H = self.thermoUnsrt[index].cholskyDeCorrelateMat["h"]
				zeta_H = self.thermoUnsrt[index].zeta["h"]
				coeff_H = np.array([f]) + np.array(np.dot(choleskyMatrix_H,beta_h*zeta_H)).flatten()		
			else:
				coeff_H = np.array([f])
			
			if "e" in parameters:
				choleskyMatrix_S = self.thermoUnsrt[index].cholskyDeCorrelateMat["e"]
				zeta_S = self.thermoUnsrt[index].zeta["e"]
				coeff_S = np.array([g]) + np.array(np.dot(choleskyMatrix_S,beta_s*zeta_S)).flatten()		
			else:
				coeff_S = np.array([g])
			
			
			thermoParams = thermoParams.replace(aString,"{:.8E}".format(coeff_cp[0]))
			thermoParams = thermoParams.replace(bString,'{:.8E}'.format(coeff_cp[1]))
			thermoParams = thermoParams.replace(cString,'{:.8E}'.format(coeff_cp[2]))
			thermoParams = thermoParams.replace(dString,'{:.8E}'.format(coeff_cp[3]))
			thermoParams = thermoParams.replace(eString,'{:.8E}'.format(coeff_cp[4]))
			thermoParams = thermoParams.replace(fString,'{:.8E}'.format(coeff_H[0]))
			thermoParams = thermoParams.replace(gString,'{:.8E}'.format(coeff_S[0]))
		else:
			thermoParams = thermoParams
		return thermoParams

	def getNewLJwellDepth(self,transParams,index,beta,isOpt,delta):
			#searching the Lennard-Jonnes potential well depth (\epsilon_{ij})
		if delta == 1:
			pattern = re.compile(r'(\w*)\s*(\d)\s*(\d+\.\d*|\d*\.\d+)\s*(\d+\.\d*|\d*\.\d+)\s*(\d+\.\d*|\d*\.\d+)\s*(\d+\.\d*|\d*\.\d+)\s*(\d+\.\d*|\d*\.\d+)',re.DOTALL|re.IGNORECASE)
			match = re.search(pattern,transParams)
			potentialWellDepth = match.group(3)
			epsilon = float(potentialWellDepth)
			parameters = self.transUnsrt[index].branches
			if "LJe" in parameters:
				choleskyMatrix = self.transUnsrt[index].cholskyDeCorrelateMat["LJe"]
				zeta = self.transUnsrt[index].zeta["LJe"]
				newEpsilon = np.array([epsilon]) + np.dot(choleskyMatrix,beta*zeta).flatten()
			else:
				newEpsilon = np.array([epsilon])
			#need to convert back to original values
			#print("\tBeta (LJe):\t{}\n".format(beta))
			transParams = transParams.replace(potentialWellDepth,'{:.3f}'.format(float(newEpsilon[0])))
		else:
			transParams = transParams
		return transParams	
	

	def getNewDiameter(self,transParams,index,beta,isOpt,delta):
		#searching the Lennard-Jonnes potential hard sphere diameter (\sigma_{ij})		
		if delta == 1:
			pattern = re.compile(r'(\w*)\s*(\d)\s*(\d+\.\d*|\d*\.\d+)\s*(\d+\.\d*|\d*\.\d+)\s*(\d+\.\d*|\d*\.\d+)\s*(\d+\.\d*|\d*\.\d+)\s*(\d+\.\d*|\d*\.\d+)',re.DOTALL|re.IGNORECASE)
			match = re.search(pattern,transParams)
			diameter = match.group(4)
			d = float(diameter)
			parameters = self.transUnsrt[index].branches
			if "LJs" in parameters:
				choleskyMatrix = self.transUnsrt[index].cholskyDeCorrelateMat["LJs"]
				zeta = self.transUnsrt[index].zeta["LJs"]
				newDia = np.array([d]) + np.dot(choleskyMatrix,beta*zeta).flatten()
			else:
				newDia = np.array([d])
			#need to convert back to original values

			transParams = transParams.replace(diameter,'{:.3f}'.format(float(newDia[0])))
			#print("\tBeta (LJs):\t{}\n".format(beta))
		else:
			transParams = transParams
		return transParams	
		
	
	def filter_list(self,List):
		temp = []
		for i in List:
			if i != "":
				temp.append(i)
		return temp
	
	
	def getNewCollisionEff(self,line,mList,beta,unsrt,rxn):
		if self.fileType == "FlameMaster":
			#molecule = index.split("-")[1]
			for index,m in enumerate(mList):
				species = m
				key = []
				param=[]		
				pattern = re.compile(r'(let?\s*\S*\s*?)=.*?(\S*?\s*?\[.*\].).*?')
				match1=re.search(pattern,line.split("\n")[0])
				#print(match1.group(0))
				header = match1.group(1)
				line_ = match1.group(2).split('+')
				#print(line_)
				for j in line_:
					list_= self.filter_list(j.split(" "))        
					if species in list_[1]:
						#print(str(float(beta)*float(list_[0])))
						choleskyMatrix = self.tbdUnsrt[rxn].cholskyDeCorrelateMat[m]
						zeta = self.tbdUnsrt[rxn].zeta[m]
						a = float(list_[0])
						#print("\tBeta ({}):\t{}\n".format(j,beta[index]))
						newCollisionEff = np.array([a])+np.asarray(np.dot(choleskyMatrix,beta[index]*zeta)).flatten()
						list_[0] = list_[0].replace(list_[0],str(newCollisionEff[0]))
						list_ = ' '.join(list_)
						param.append(list_)
				#		print(list_)
					else:
						list_ = ' '.join(list_)
						param.append(list_)
				#print(param)		
				newLine = " + ".join(param)
				key.append(header)
				key.append(newLine)
				newLine = " = ".join(key)
				line = newLine
				#print(line)
				#print(newLine)
		elif self.fileType == "chemkin":
			#molecule = index.split("-")[1]
			for index,m in enumerate(mList):
				key = []
				param=[]		
				
				if self.chemtagCollision[rxn] == 2:
					w = re.compile(r'(.*?\<\=\>.*?) \s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?(.*?\/\s*?.*)\/\s*?.*?\n?',re.DOTALL|re.IGNORECASE)
					match1 = re.search(w,line)
					String = match1.group(11)	
					line_content = match1.group(11).split('/')
				
				
				elif self.chemtagCollision[rxn] == 4:
					w = re.compile(r'(.*?\<\=\>.*?) \s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?(LOW)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?\/\s*?.*?\n?\s*?(TROE)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?\/\s*?\n(.*)',re.DOTALL|re.IGNORECASE)
					match1 = re.search(w,line)
					String = match1.group(34)	
					line_content = match1.group(34).split('/')
				#print(line_)
				index = []
				newCollisionEff = []
				for element in line_content:     
					if m in element:
						
						choleskyMatrix = self.tbdUnsrt[rxn].cholskyDeCorrelateMat[m]
						zeta = self.tbdUnsrt[rxn].zeta[m]
						index.append(line_content.index(element)+1)
						val = float(line_content[line_content.index(element)+1])
						#print("\tBeta ({}):\t{}\n".format(j,beta[index]))
						newCollisionEff.append(np.array([val])+np.asarray(np.dot(choleskyMatrix,beta[index]*zeta)).flatten())
					else:
						continue
				for ind,val in index:
					line_content[val] =  newCollisionEff[ind]
				#print(param)		
				newString = "/".join(line_content)
			 
				newLine = newLine.replace(String,newString)
				#print(line)
				#print(newLine)
		
		return newLine
		
	
	def ModifyReaction(self,reaction_str,target,beta,rxn, sens_AL,isOpt):
		#print(type(index))
		
		rxn_type = self.rxnUnsrt[rxn].type
		branch_bool = self.rxnUnsrt[rxn].branching
		
		branches = self.rxnUnsrt[rxn].branches.split(",")
		
		#print(rxn_type)
		#print(IsBranching)
		#print(branch)
		#print("\n\nbeta ({}):\t{}\n".format(index,beta))
		
		#print("\nBefore perturbation:\n\n{}".format(reaction))
		if rxn_type.strip() == "pressure_independent":
			if branch_bool.strip() == "True":
				#print("selected {} and branching is {} \n \t before perturbation {}\n".format(rxn,self.rxnUnsrt[rxn].rxn_branching,reaction_str))
				#print("Before perturbation = {}".format(reaction_str))
		#		print("\nReaction type:\n\n{}".format(rxn_type[0]))
		#		print("\nIs Branching?:\n\n{}".format(IsBranching))
				
				#print(reaction_str)
				
				reaction_str = self.getNewBranchRxn(reaction_str,branches,rxn,beta,isOpt)
				#print("\t After perturbation {} \n".format(reaction_str))
			else:
		#		print("\nReaction type:\n\n{}".format(rxn_type[0]))
		#		print("\nIs Branching?:\n\n{}".format(IsBranching))
				#print("Before perturbation = {}".format(reaction_str))
				#print("selected {} and branching is {} \n".format(rxn,self.rxnUnsrt[rxn].rxn_branching))
				#print(rxn)
				#print(reaction_str)
				#print(beta)
				reaction_str = self.getNewRxn(reaction_str,rxn,beta)			 
					
		elif rxn_type.strip() == "pressure_dependent":
			#print(rxn)
			#print(rxn[0])
			#print("Before perturbation = {}".format(reaction_str))
			kDelta_HPL = kDelta_LPL = kDelTa_FoC = 0 
			if "HPL" in branches.split(","):
				kDelta_HPL = 1
			if "LPL" in branches.split(","):
				kDelta_LPL = 1
			if "FoC" in branches.split(","):
				kDelta_FoC = 1
				#print("\nReaction type:\n\n{}".format(rxn_type[0]))
				#print("Foc FOUND IN LIST")
				reaction_str = self.getNewUpperLimits(reaction_str,rxn,branch,beta,isOpt,kDelta_HPL)
				reaction_str = self.getNewLowerLimits(reaction_str,rxn,branch,beta,isOpt,kDelta_LPL)
				reaction_str = self.getNewFallOffCurve(reaction_str,target,rxn,branch,beta,isOpt,kDelta_FoC)	
			
			
		else: 
			print("Plz give correct input\n")	
		
		#print("\nAfter perturbation\n\n{}".format(reaction))		 
		#print("\t After perturbation = {}".format(reaction_str))
		return reaction_str 
	
	def ModifyThermo(self,thermoParams,beta,index, sens_AL,isOpt):
		#Unsrt_type = self.uncert[index].getUnsrtType()
		branches = self.thermoUnsrt[index].branches
		kDelta_low = 0
		kDelta_high = 0
		thermoBeta_low = 0*np.ones(7)
		thermoBeta_high = 0*np.ones(7)
		
		if "Low" in branches.split(","):
			kDelta_low  = 1
			thermoBeta_low = self.RBeta_Dict[index] 
		if "High" in branches.split(","):
			kDelta_high = 1
			thermoBeta_high = self.RBeta_Dict[index]
		
		#print(thermoBeta_low)
		#print(thermoBeta_high)
		
		#print(knockerDelta_low)	
		thermoParams = self.getNewThermo_low(thermoParams,index,thermoBeta_low,isOpt,kDelta_low)
		
		thermoParams = self.getNewThermo_high(thermoParams,index,thermoBeta_high,isOpt,kDelta_high)
				
		return thermoParams




	def ModifyTransport(self,transportParams,beta,index, sens_AL,isOpt):
		#print("Transport parameters before perturbation:\n\n{}\n".format(transportParams))
		branches = self.transUnsrt[index].branches
		kDelta_LJs = 0
		kDelta_LJe = 0
		beta_LJe = 0*np.ones(1)
		beta_LJs = 0*np.ones(1)
		if "LJe" in branches.split(","):
			#print(self.RBeta_Dict[index])
			beta_LJe =self.RBeta_Dict[index]["LJe"]
			kDelta_LJe = 1 
		if "LJs" in branches.split(","):
			beta_LJs = self.RBeta_Dict[index]["LJs"]
			kDelta_LJs = 1 
			
			transportParams = self.getNewLJwellDepth(transportParams,index,beta_LJe,isOpt,kDelta_LJe)
			transportParams = self.getNewDiameter(transportParams,index,beta_LJs,isOpt,kDelta_LJs)
		
		return transportParams 
	
	
	##Remains to modify for chemkin format
	def ModifyLine(self,line,beta,rxn):
		#print("Third body collision efficiencies before perturbation:\n\n{}\n".format(line))
		unsrt = []	
		for mol_list in self.molecule:
			for mol in mol_list:
				unsrt.append(self.tbdUnsrt[rxn].uncertainties[mol])
			newLine = self.getNewCollisionEff(line,mol_list,beta,unsrt,rxn)
		newLine = newLine+'\n'
		#print("Third body collision efficiencies after perturbation:\n\n{}\n".format(newLine))		 
		return newLine    
	
            
        
        
