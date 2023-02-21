import numpy as np
import re,os
import Uncertainty as uncertainty
import Input_file_reader
import matplotlib.pyplot as plt
from collections import OrderedDict
import math
class MechanismManipulator():
	def __init__(self, filePath_mech,fileType,filePath_therm,filePath_trans,reactionList,focList,thirdBodyList,thermoList,transportList, rUnsrt,fUnsrt,mUnsrt,thUnsrt,trUnsrt,selectedParams,activeIndexDict,design_type=None):
		#print("Hello")
		if design_type != None:
			self.design= design_type
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
		#print(self.mList)
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
		self.reaction_search = []
		for rxn in self.ReactionList:
			if rxn.split(":")[0] not in self.reaction_search:
				self.reaction_search.append(rxn.split(":")[0])
		#print(reaction_search)
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
					if r.split(":")[0] in i: 
						#if self.rxnUnsrt[r].type == "pressure-independent" and self.rxnUnsrt[r].sub_type == "forward":
						if "DUPLICATE" in f_mech[f_mech.index(i)+1]:
							#print(r)
							self.chemtag[r] = 1
							start = f_mech.index(i)
							for j in range(4):
								mechLine +=f_mech[start+j]
							#print(mechLine)
							self.PerturbingReactions[r] = mechLine
							#print(mechLine)
							mechLine = "" 
							break
						
						elif "LOW" not in f_mech[f_mech.index(i)+1] and "/" in f_mech[f_mech.index(i)+1] and "LOW" not in f_mech[f_mech.index(i)+2] and "TROE" not in f_mech[f_mech.index(i)+2]:
							#print(r)
							self.chemtag[r] = 2
							start = f_mech.index(i)  
							mechLine += f_mech[start]
							index = 1
							while "/" in f_mech[start+index]:
								mechLine += f_mech[start+index]
								index+=1

							self.PerturbingReactions[r] = mechLine
							mechLine = "" 
						
						elif "LOW" in f_mech[f_mech.index(i)+1] and "TROE" in f_mech[f_mech.index(i)+2]and "/" not in f_mech[f_mech.index(i)+3]:
							self.chemtag[r] = 3
							start = f_mech.index(i)  
							mechLine += f_mech[start]
							index = 1
							while "LOW" in f_mech[start+index] or "TROE" in f_mech[start+index]:
								mechLine += f_mech[start+index]
								index+=1
							self.PerturbingReactions[r] = mechLine
							mechLine = "" 
						
						elif "LOW" in f_mech[f_mech.index(i)+1] and "TROE" in f_mech[f_mech.index(i)+2] and "/" in f_mech[f_mech.index(i)+3]:
							self.chemtag[r] = 4
							start = f_mech.index(i)  
							mechLine += f_mech[start]
							index = 1
							while "/" in f_mech[start+index] or "LOW" in f_mech[start+index] or "TROE" in f_mech[start+index]:
								mechLine += f_mech[start+index]
								index+=1
							self.PerturbingReactions[r] = mechLine
							mechLine = "" 
						elif "LOW" in f_mech[f_mech.index(i)+1] and "TROE" not in f_mech[f_mech.index(i)+2] and "/" in f_mech[f_mech.index(i)+2]:
							self.chemtag[r] = 5
							start = f_mech.index(i)  
							mechLine += f_mech[start]
							index = 1
							while "/" in f_mech[start+index] or "LOW" in f_mech[start+index]:
								mechLine += f_mech[start+index]
								index+=1
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
							mechLine += f_mech[start]
							index = 1
							while "LOW" in f_mech[start+index] or "TROE" in f_mech[start+index]:
								mechLine += f_mech[start+index]
								index+=1
							self.PerturbingFoC[foc] = mechLine
							mechLine = ""
						
						elif "LOW" in f_mech[f_mech.index(i)+1] and "TROE" in f_mech[f_mech.index(i)+2] and "/" in f_mech[f_mech.index(i)+3]:
							self.chemtag[r] = 4
							start = f_mech.index(i)  
							mechLine += f_mech[start]
							index = 1
							while "/" in f_mech[start+index] or "LOW" in f_mech[start+index] or "TROE" in f_mech[start+index]:
								mechLine += f_mech[start+index]
								index+=1
							self.PerturbingFoC[foc] = mechLine
							mechLine = ""
						else:
							 print("The input reaction is not a TROE type pressure dependent reaction")	
		
			for m in self.mList:
				for i in f_mech:
					if m.split(":")[0] in i:
						if "LOW" not in f_mech[f_mech.index(i)+1] and "/" in f_mech[f_mech.index(i)+1]:
							self.chemtagCollision[m] = 2
							start = f_mech.index(i)  
							mechLine += f_mech[start]
							index = 1
							while "/" in f_mech[start+index]:
								mechLine += f_mech[start+index]
								index+=1
							self.PerturbingCollisionEff[m] = mechLine
							mechLine = ""
						
						elif "LOW" in f_mech[f_mech.index(i)+1] and "TROE" in f_mech[f_mech.index(i)+2] and "/" in f_mech[f_mech.index(i)+3]:
							self.chemtagCollision[m] = 4
							start = f_mech.index(i)  
							mechLine += f_mech[start]
							index = 1
							while "/" in f_mech[start+index] or "LOW" in f_mech[start+index] or "TROE" in f_mech[start+index]:
								mechLine += f_mech[start+index]
								index+=1
							self.PerturbingCollisionEff[m] = mechLine
							mechLine = ""
						elif "LOW" in f_mech[f_mech.index(i)+1] and "TROE" not in f_mech[f_mech.index(i)+2] and "/" in f_mech[f_mech.index(i)+2]:
							self.chemtagCollision[m] = 5
							start = f_mech.index(i)  
							mechLine += f_mech[start]
							index = 1
							while "/" in f_mech[start+index] or "LOW" in f_mech[start+index]:
								mechLine += f_mech[start+index]
								index+=1
							self.PerturbingCollisionEff[m] = mechLine
							mechLine = "" 
						else:
							print("The input reaction is not a third body reaction")  
							break
		#print(self.PerturbingReactions)	
		#print(self.chemtag)
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
	def GeneratePerturbedMechanism(self,target,beta,transformation,reactionList,Sens_AL,isOpt,extra_arg = None):
		NewMechanismString = self.parentFile_mech
		NewThermoString = self.parentFile_thermo
		NewTransportString = self.parentFile_transport
		self.RBeta_Dict = {}
		self.RBeta_select_Dict = {}
		#print(f"Beta is {beta}\ntransformation is {transformation}\n")
		#print(NewMechanismString)
		#Beta is an array
		beta_group = []
		
		######################
		#    REACTION        #
		######################
		
		if len(self.ReactionList)!= 0:
			count_Arrhenius_params = 0
			temp = {}
			for i in self.ReactionList:
				if "all" in self.rxnUnsrt[i].perturbation_type:
					count_Arrhenius_params+=3
					temp[i] = 3
				elif "factor" in self.rxnUnsrt[i].perturbation_type:
					temp[i] = 1
					count_Arrhenius_params+=1
				else:
					raise AssertionError("Incorrect perturbation type")
			#print("rxn founds")
			#print(count_Arrhenius_params)
			#print(count_Arrhenius_params)
			beta_rxn = beta[0:count_Arrhenius_params]
			beta_selection = transformation[0:count_Arrhenius_params]
			beta_per_rxn = []
			beta_selection_per_rxn = []
			count_rxn = 0
			for i in self.ReactionList:
				zemp = []
				ramp = []
				if count_rxn < count_Arrhenius_params:
					for j in range(temp[i]):
						zemp.append(beta_rxn[count_rxn])
						ramp.append(beta_selection[count_rxn])
						count_rxn+=1
					beta_per_rxn.append(zemp)
					beta_selection_per_rxn.append(ramp)
			#print("\nBeta selection per reaction\n")
			#print(beta_selection_per_rxn)
			
			for index, rxn in enumerate(self.ReactionList):
				self.RBeta_Dict[rxn] = np.asarray(beta_per_rxn[index])
			for index,rxn in enumerate(self.ReactionList):
				self.RBeta_select_Dict[rxn] = beta_selection_per_rxn[index]
			
			#print("\nDictionary in Mechmanipulator\n")
			#print(self.RBeta_select_Dict)
			#raise AssertionError("Stop")
			for index,rxn in enumerate(self.PerturbingReactions):
				#print("Rxn is {}".format(rxn))
				#print("Old reaction\n")
				#print(self.PerturbingReactions[rxn])
				#print("\n")
				#print("Beta is \n")
				#print(self.RBeta_Dict[rxn])
				#print("zeta is \n")
				#print(self.rxnUnsrt[rxn].zeta.x)
				NewReaction =  self.ModifyReaction(self.PerturbingReactions[rxn],target,np.asarray(self.RBeta_Dict[rxn]),np.asarray(self.RBeta_select_Dict[rxn]).T,rxn,Sens_AL,isOpt)
				#print("New reaction after perturbation\n")
				#print(NewReaction)
				#print("\n")
				NewMechanismString = NewMechanismString.replace(self.PerturbingReactions[rxn],NewReaction)	
					
		else:
			#print("rxn not founds")
			count_rxn = 0
			NewMechanismString = NewMechanismString
		#print("Changing Mechanism")
		#print(NewMechanismString)
		#print("\n")
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
		f_mech = NewMechanismString.splitlines()
		self.PerturbingCollisionEff = {}
		mechLine = ""
		if "FlameMaster" in self.fileType:				
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
		elif "chemkin" in self.fileType :
			for m in self.mList:
				for i in f_mech:
					if m.split(":")[0] in i:
						if "LOW" not in f_mech[f_mech.index(i)+1] and "/" in f_mech[f_mech.index(i)+1]:
							self.chemtagCollision[m] = 2
							start = f_mech.index(i)  
							for j in range(2):
								mechLine +=f_mech[start+j]+"\n"
							self.PerturbingCollisionEff[m] = mechLine
							mechLine = ""
						
						elif "LOW" in f_mech[f_mech.index(i)+1] and "TROE" in f_mech[f_mech.index(i)+2] and "/" in f_mech[f_mech.index(i)+3]:
							self.chemtagCollision[m] = 4
							start = f_mech.index(i)  
							for j in range(4):
								mechLine += f_mech[start+j]+"\n"
							self.PerturbingCollisionEff[m] = mechLine
							mechLine = ""
						else:
							print("The input reaction is not a third body reaction")  
							break
		
		
		
		
		#print(self.PerturbingCollisionEff)
		
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
				#print("Old line \n")
				#print(self.PerturbingCollisionEff[rxn])
				#print("\n")
				NewLine = self.ModifyLine(self.PerturbingCollisionEff[rxn],beta_group[index],rxn)
				#print("New line")
				#print(NewLine)
				#print("\n")
				#print("Before changing MBody")
				#print(NewMechanismString)
				NewMechanismString = NewMechanismString.replace(self.PerturbingCollisionEff[rxn],NewLine)
				#print("After changing MBody")
				#print(NewMechanismString)
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
		locations = {}
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
			locations = {}
			#locations["opt_mech"] = os.gecwd()+"/mechanism_opt.mech"
			#locations["opt_thermo"] = os.gecwd()+"/thermo_opt.therm"
			#locations["opt_trans"] = os.gecwd()+"/transport_opt.trans"
			#return NewMechanismString,NewThermoString,NewTransportString,self.RBeta_Dict,locations
	
		#else:
			#mechFile = open("mechanism.mech",'w')
			#mechFile.write(NewMechanismString)
			#mechFile.close()
	
			#thermoFile = open("thermo.therm",'w')
			#thermoFile.write(NewThermoString)
			#thermoFile.close()
	
			#transportFile = open("transport.trans",'w')
			#transportFile.write(NewTransportString)
			#transportFile.close()
			#return NewMechanismString,NewThermoString,NewTransportString,self.RBeta_Dict
		return NewMechanismString,NewThermoString,NewTransportString,self.RBeta_Dict

	def ModifyReaction(self,reaction_str,target,beta,selection,rxn, sens_AL,isOpt):
		#print(rxn)
		
		#print(type(index))
		rxn_type = self.rxnUnsrt[rxn].type
		branch_bool = self.rxnUnsrt[rxn].branching
		rxn_sub_type = self.rxnUnsrt[rxn].sub_type
		#print(rxn_sub_type)
		branches = self.rxnUnsrt[rxn].branches.split(",")
		
		#print(rxn_type)
		#print(IsBranching)
		#print(branch)
		#print("\n\nbeta ({}):\t{}\n".format(index,beta))
		
		#print("\nBefore perturbation:\n\n{}".format(reaction))
		if rxn_type.strip() == "pressure_independent":
			if rxn_sub_type.strip() == "branching":
				#print("selected {} and branching is {} \n \t before perturbation {}\n".format(rxn,self.rxnUnsrt[rxn].rxn_branching,reaction_str))
				#print("Before perturbation = {}".format(reaction_str))
		#		print("\nReaction type:\n\n{}".format(rxn_type[0]))
		#		print("\nIs Branching?:\n\n{}".format(IsBranching))
				
				#print(reaction_str)
				
				reaction_str = self.getNewBranchRxn(reaction_str,branches,rxn,beta,selection,isOpt)
				#print("\t After perturbation {} \n".format(reaction_str))
			elif rxn_sub_type.strip() == "duplicate":
				#print("Duplicate")
				base_rxn = rxn.split(":")[0]
				#print(base_rxn)
				kdup_A = 0
				kdup_B = 0
				if base_rxn+":A" in self.ReactionList:
					kdup_A = 1
				if base_rxn+":B" in self.ReactionList:
					kdup_B = 1
				#print("Going to the func")
				reaction_str = self.getNewDuplicateRxn(reaction_str,rxn,beta,selection,isOpt,"A",kdup_A,kdup_B)
				#reaction_str = self.getNewDuplicateRxn(reaction_str,base_rxn+":B",beta,isOpt,"B",0,kdup_B)
			else:
				#print("Forward Rxn")
		#		print("\nReaction type:\n\n{}".format(rxn_type[0]))
		#		print("\nIs Branching?:\n\n{}".format(IsBranching))
				#print("Before perturbation = {}".format(reaction_str))
				#print("selected {} and branching is {} \n".format(rxn,self.rxnUnsrt[rxn].rxn_branching))
				#print(rxn)
				#print(reaction_str)
				#print(beta)
				reaction_str = self.getNewRxn(reaction_str,rxn,beta,selection,isOpt)			 
					
		elif rxn_type.strip() == "pressure_dependent":
			#print(rxn)
			#print(rxn[0])
			#print("Before perturbation = {}".format(reaction_str))
			base_rxn = rxn.split(":")[0]
			kDelta_HPL = 0
			kDelta_LPL = 0
			kDelta_FoC = 0 
			if base_rxn+":High" in self.ReactionList:
				kDelta_HPL = 1
			if base_rxn+":Low" in self.ReactionList:
				kDelta_LPL = 1
			if base_rxn+":FoC" in self.ReactionList:
				kDelta_FoC = 1
				#print("\nReaction type:\n\n{}".format(rxn_type[0]))
				#print("Foc FOUND IN LIST")
			#print(kDelta_HPL,kDelta_LPL,kDelta_FoC)
			
			reaction_str = self.getNewUpperLimits(reaction_str,base_rxn+":High",beta,selection,isOpt,kDelta_HPL)
			reaction_str = self.getNewLowerLimits(reaction_str,base_rxn+":Low",beta,selection,isOpt,kDelta_LPL)
			reaction_str = self.getNewFallOffCurve(reaction_str,target,base_rxn+":FoC",beta,selection,isOpt,kDelta_FoC)	
		
			
		else: 
			print("Plz give correct input\n")	
		
		#print("\nAfter perturbation\n\n{}".format(reaction))		 
		#print("\t After perturbation = {}".format(reaction_str))
		return reaction_str 
	

	
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
		
	def getNewDuplicateRxn(self,reaction,rxn,beta,selection,isOpt,whichRxn,kDa,kDb):
		#print(self.fileType)
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
			##print(f"Inside the chemkin loop\n{reaction}")
			w = re.compile(r'(.*?\<\=\>.*? \s*?)((-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))(\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?))(\n?\s*?DUPLICATE\s*?\n?(.*?\<\=\>.*?) \s*?)((-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))(\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+)))\s*?.*?)(\n?\s*?DUPLICATE\s*?\n?)',re.DOTALL|re.IGNORECASE)

		#print("\n")
			match = re.search(w,reaction)
			RxnA = match.group(1)
			ratesA = match.group(2)
			dup_flag1_rxnB = match.group(13)
			ratesB = match.group(15)
			dup_flag2 = match.group(26)
			
			pattern = re.compile(r'((-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))(\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?))',re.DOTALL|re.IGNORECASE) 
			
			match1 = re.search(pattern,ratesA)
			aString = match1.group(2)
			a = float(aString)
			alpha = np.log(a)

			nString = match1.group(6)
			n = float(nString)

			eString = match1.group(9)
			e = float(eString)
			epsilon = e/1.987
			
			match2 = re.search(pattern,ratesB)
			aString_dup = match2.group(2)
			a_dup = float(aString_dup)
			alpha_dup = np.log(a_dup)

			nString_dup = match2.group(6)
			n_dup = float(nString_dup)

			eString_dup = match2.group(9)
			e_dup = float(eString_dup)
			epsilon_dup = e_dup/1.987
			
			#print(f"kDa = {kDa} kDb = {kDb}")
			
			if kDa == 1 and kDb == 0:
				
				M = 3/(np.log(10.0));		
				choleskyMatrix = self.rxnUnsrt[rxn].cholskyDeCorrelateMat
				zeta = self.rxnUnsrt[rxn].zeta_matrix
				beta = self.RBeta_Dict[rxn]
				#print(np.dot(choleskyMatrix,beta*zeta).flatten())		
				##print(f"selected: {selection}\n")
				##print(f"Before selection rxn : {reaction}\n")	
				if self.design == "B1" or self.design == "B1-Mix" or self.design == "B1-factorial" or self.design == "B1-MC" or  self.design == "SAMAP":
					log_Arrhenius = np.array([alpha,n,epsilon]) + np.asarray(selection*np.asarray(np.dot(choleskyMatrix,beta)).flatten()).flatten()
				elif self.design == "A1":
					log_Arrhenius = np.array([alpha,n,epsilon]) + np.asarray(selection*np.asarray(np.dot(choleskyMatrix,beta)).flatten()).flatten()
				elif self.design == "Monte-Carlo-All":
					log_Arrhenius = np.array([alpha,n,epsilon]) + np.asarray(selection*np.asarray(np.dot(choleskyMatrix,zeta.T.dot(beta.T).T)).flatten()).flatten()
				#need to convert back to original values
				Arrhenius = [np.exp(log_Arrhenius[0]), log_Arrhenius[1], 1.987*log_Arrhenius[2]]			
				aString = '{:.3E}'.format(Arrhenius[0])
				ga = "   "
				nString = '{:.3f}'.format(Arrhenius[1])
				gb = "   "
				eString = '{:.3E}'.format(Arrhenius[2])
				
				ratesA = aString+ga+nString+gb+eString
				reaction = RxnA + ratesA + dup_flag1_rxnB + ratesB + dup_flag2
				#print(f"After selection rxn : {reaction}\n")	
			elif kDa == 0 and kDb == 1:
				
				M = 3/(np.log(10.0));		
				choleskyMatrix = self.rxnUnsrt[rxn].cholskyDeCorrelateMat
				zeta = self.rxnUnsrt[rxn].zeta_matrix
				beta = self.RBeta_Dict[rxn]
				#print(f"selected: {selection}\n")
				#print(f"Before selection rxn : {reaction}\n")
				#print(np.dot(choleskyMatrix,beta*zeta).flatten())
				if self.design == "B1" or self.design == "B1-Mix" or self.design == "B1-factorial" or self.design == "B1-MC" or  self.design == "SAMAP":
					log_Arrhenius = np.array([alpha_dup,n_dup,epsilon_dup]) + np.asarray(selection*np.asarray(np.dot(choleskyMatrix,beta)).flatten()).flatten()
				elif self.design == "A1":
					log_Arrhenius = np.array([alpha_dup,n_dup,epsilon_dup]) + np.asarray(selection*np.asarray(np.dot(choleskyMatrix,beta)).flatten()).flatten()
				elif self.design == "Monte-Carlo-All":
					log_Arrhenius = np.array([alpha_dup,n_dup,epsilon_dup]) + np.asarray(selection*np.asarray(np.dot(choleskyMatrix,zeta.T.dot(beta.T).T)).flatten()).flatten()
				#need to convert back to original values
				Arrhenius = [np.exp(log_Arrhenius[0]), log_Arrhenius[1], 1.987*log_Arrhenius[2]]			
			
				aString_dup = '{:.3E}'.format(Arrhenius[0])
				ga = "   "
				nString_dup = '{:.3f}'.format(Arrhenius[1])
				gb = "   "
				eString_dup = '{:.3E}'.format(Arrhenius[2])
				
				ratesB = aString_dup+ga+nString_dup+gb+eString_dup
			
				reaction = RxnA + ratesA + dup_flag1_rxnB + ratesB + dup_flag2
				#print(f"After selection rxn : {reaction}\n")		
			elif kDa == 1 and kDb == 1:
				rxn_a = rxn.split(":")[0]+":A"
				rxn_b = rxn.split(":")[0]+":B"
				
				choleskyMatrix_a = self.rxnUnsrt[rxn_a].cholskyDeCorrelateMat
				
				choleskyMatrix_b = self.rxnUnsrt[rxn_b].cholskyDeCorrelateMat
				
				M = 3/(np.log(10.0));		
				
				choleskyMatrix_a = self.rxnUnsrt[rxn_a].cholskyDeCorrelateMat
				zeta_a = self.rxnUnsrt[rxn_a].zeta_matrix
				beta_a = self.RBeta_Dict[rxn_a]
				zeta_b = self.rxnUnsrt[rxn_b].zeta_matrix
				beta_b = self.RBeta_Dict[rxn_b]
				selection_a = self.RBeta_select_Dict[rxn_a]
				selection_b = self.RBeta_select_Dict[rxn_b]
				#print(f"selected: {selection_a}\n")
				#print(f"selected: {selection_b}\n")
				#print(f"Before selection_a rxn : {reaction}\n")
				if self.design == "B1" or self.design == "B1-Mix" or self.design == "B1-factorial" or self.design == "B1-MC" or  self.design == "SAMAP":
					log_Arrhenius_a = np.array([alpha,n,epsilon]) + np.asarray(selection_a*np.asarray(np.dot(choleskyMatrix_a,beta_a)).flatten()).flatten()
				elif self.design == "A1":
					log_Arrhenius_a = np.array([alpha,n,epsilon]) + np.asarray(selection_a*np.asarray(np.dot(choleskyMatrix_a,beta_a)).flatten()).flatten()
				elif self.design == "Monte-Carlo-All":	
					log_Arrhenius_a = np.array([alpha,n,epsilon]) + np.asarray(selection_a*np.asarray(np.dot(choleskyMatrix_a,zeta_a.T.dot(beta_a.T).T)).flatten()).flatten()
				#need to convert back to original values
				Arrhenius_a = [np.exp(log_Arrhenius_a[0]), log_Arrhenius_a[1], 1.987*log_Arrhenius_a[2]]			
				#print(np.dot(choleskyMatrix,beta*zeta).flatten())		
				
				#need to convert back to original values
						
				aString = '{:.3E}'.format(Arrhenius_a[0])
				ga = "   "
				nString = '{:.3f}'.format(Arrhenius_a[1])
				gb = "   "
				eString = '{:.3E}'.format(Arrhenius_a[2])
				
				ratesA = aString+ga+nString+gb+eString
					
				#print(np.dot(choleskyMatrix,beta*zeta).flatten())
				log_Arrhenius_b = np.array([alpha_dup,n_dup,epsilon_dup]) + np.asarray(selection_b*np.asarray(np.dot(choleskyMatrix_b,beta_b)).flatten()).flatten()
				#log_Arrhenius_b = np.array([alpha_dup,n_dup,epsilon_dup]) + np.asarray(selection_b*np.asarray(np.dot(choleskyMatrix_b,zeta_b.T.dot(beta_b.T).T)).flatten()).flatten()
				#need to convert back to original values
				Arrhenius_b = [np.exp(log_Arrhenius_b[0]), log_Arrhenius_b[1], 1.987*log_Arrhenius_b[2]]			
			
				aString_dup = '{:.3E}'.format(Arrhenius_b[0])
				ga = "   "
				nString_dup = '{:.3f}'.format(Arrhenius_b[1])
				gb = "   "
				eString_dup = '{:.3E}'.format(Arrhenius_b[2])
				
				ratesB = aString_dup+ga+nString_dup+gb+eString_dup
			
				reaction = RxnA + ratesA + dup_flag1_rxnB + ratesB + dup_flag2		
				#print(f"After selection rxn : {reaction}\n")
			else:
				reaction = reaction

			#print(reaction)
		else:
			print("Invalid kinetic mechanism file")
			
		return reaction

	def getNewUpperLimits(self,reaction,rxn,beta,selection,isOpt,delta):
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
				zeta = self.rxnUnsrt[rxn].zeta_matrix
		
				beta = self.RBeta_Dict[rxn]
				#print("\tBeta (HPL):\t{}\n".format(beta_i))
				log_Arrhenius_i = np.array([alpha_i,ni,epsilon_i]) + np.asarray(selection*np.asarray(np.dot(choleskyMatrix,zeta.T.dot(beta.T).T)).flatten()).flatten()
				#need to convert back to original values
				Arrhenius_i = [np.exp(log_Arrhenius_i[0]), log_Arrhenius_i[1], 8.314E-3*log_Arrhenius_i[2]]
	
				reaction = reaction.replace(aiString,'{:.3E}'.format(Arrhenius_i[0]))
				reaction = reaction.replace(niString,'{:.3f}'.format(Arrhenius_i[1]))
				reaction = reaction.replace(eiString,'{:.3f}'.format(Arrhenius_i[2]))
			elif self.fileType == "chemkin":
				#print(rxn)
				#print(self.chemtag[rxn])
				if self.chemtag[rxn] == 2:
				
					w = re.compile(r'(.*?\<\=\>.*? \s*?)((-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?)(\s*?.*?\/\s*?.*\/?\s*?.*?\n?)',re.DOTALL|re.IGNORECASE)
					
					match1 = re.search(w,reaction)
					front = match1.group(1)
					high = match1.group(2)
					rest_of_string = match1.group(12)
														
					#groups
					#reaction = 1
					#High = 2
					#m = 12
				elif self.chemtag[rxn] == 3:
					w = re.compile(r'(.*?\<\=\>.*? \s*?)((-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?)((((LOW)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?\/\s*?.*?\n?\s*?))((TROE)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?\/?\s*?\n?))',re.DOTALL|re.IGNORECASE)
					match1 = re.search(w,reaction)
					front = match1.group(1)
					high = match1.group(2)
					rest_of_string = match1.group(12)
				#groups
				#reaction = 1
				#High = 2
				#Low + Troe = 12
				#Low = 13
				#TROE = 23
			
				elif self.chemtag[rxn] == 4:					
					w = re.compile(r'(.*?\<\=\>.*? \s*?)((-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?)(((LOW)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?\/\s*?.*?\n?\s*?)((TROE\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d*\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d*\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d*\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?)\/?\s*?\/\s*?\n?)((.*\n?\s*?.*)?))',re.DOTALL|re.IGNORECASE)
					match1 = re.search(w,reaction)
					front = match1.group(1)
					high = match1.group(2)
					rest_of_string = match1.group(12)
				#groups
				#reaction = 1
				#High = 2
				#Low+TROE+m = 12
				#Low = 13
				#Troe = 25
				#m  = 38
				elif self.chemtag[rxn] == 5:
					w = re.compile(r'(.*?\<\=\>.*? \s*?)((-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?)(((LOW)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?\/\s*?.*?\n)(\s*?\/?\s*?\n?.*))',re.DOTALL|re.IGNORECASE)
					match1 = re.search(w,reaction)
					front = match1.group(1)
					high = match1.group(2)
					rest_of_string = match1.group(12)
					#groups
					#reaction = 1
					#High = 2
					#Low+m = 12
					#Low = 13
					#m = 24
					
				if "factor" in self.rxnUnsrt[rxn].perturbation_type:
					pattern = re.compile(r'((-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))(\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?))',re.DOTALL|re.IGNORECASE)
					
					match = re.search(pattern,high)
					ne = match.group(5)
					#A = 2
					#n = 5
					#Ea = 9
					aString = match.group(2)
					a = float(aString)
					alpha = np.log(a)
					
					choleskyMatrix = self.rxnUnsrt[rxn].perturb_factor
					zeta = self.rxnUnsrt[rxn].zeta_matrix
					beta = self.RBeta_Dict[rxn]
					selection = np.asarray(self.RBeta_select_Dict[rxn])
					#print(f"selection for High is : {selection}")
					#print(f"reaction before selection: {reaction}")
					#nString = match.group(5)
					#n = float(nString)
											
					#eString = match.group(10)
					#e = float(eString)
					#epsilon = e/1.987
					
					self.nominal_highPressureLimit = np.array([alpha])
					if self.design == "B1" or self.design == "B1-Mix" or self.design == "B1-factorial" or self.design == "B1-MC" or  self.design == "SAMAP":
						log_Arrhenius = np.array([alpha])+ np.asarray(selection*np.asarray(np.dot(choleskyMatrix,beta)).flatten()).flatten()
					elif self.design == "A1":
						log_Arrhenius = np.array([alpha])+ np.asarray(selection*np.asarray(np.dot(choleskyMatrix,beta)).flatten()).flatten()
					elif self.design == "Monte-Carlo-All":
						log_Arrhenius = np.array([alpha])+ np.asarray(selection*np.asarray(np.dot(choleskyMatrix,beta*zeta)).flatten()).flatten()
					Arrhenius_highPressureLimit = [np.exp(log_Arrhenius[0]),n,1.987*epsilon]
									
					new_high = '{:.3E}'.format(Arrhenius_highPressureLimit[0])+ne
					
					reaction = front+new_high+rest_of_string
					#print(f"reaction after selection: {reaction}")
				else:
					pattern = re.compile(r'((-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))(\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?))',re.DOTALL|re.IGNORECASE)
					match = re.search(pattern,high)
					
					#print(match.group(6))
					
					aString = match.group(2)
					a = float(aString)
					alpha = np.log(a)
					
					nString = match.group(6)
					n = float(nString)
											
					eString = match.group(10)
					e = float(eString)
					epsilon = e/1.987
					
					choleskyMatrix = self.rxnUnsrt[rxn].cholskyDeCorrelateMat
					zeta = self.rxnUnsrt[rxn].zeta_matrix
					beta = self.RBeta_Dict[rxn]
					selection = np.asarray(self.RBeta_select_Dict[rxn])
					#print(f"selection for high is : {selection}")
					#print(f"reaction before selection: {reaction}")
					
					self.nominal_highPressureLimit = np.array([alpha,n,epsilon])
					if self.design == "B1" or self.design == "B1-Mix" or self.design == "B1-factorial" or self.design == "B1-MC" or  self.design == "SAMAP":
						log_Arrhenius = np.array([alpha,n,epsilon])+ np.asarray(selection*np.asarray(np.dot(choleskyMatrix,beta)).flatten()).flatten()
					elif self.design == "A1":
						log_Arrhenius = np.array([alpha,n,epsilon])+ np.asarray(selection*np.asarray(np.dot(choleskyMatrix,beta)).flatten()).flatten()
					elif self.design == "Monte-Carlo-All":
						log_Arrhenius = np.array([alpha,n,epsilon])+ np.asarray(selection*np.asarray(np.dot(choleskyMatrix,zeta.T.dot(beta.T).T)).flatten()).flatten()
					Arrhenius_highPressureLimit = [np.exp(log_Arrhenius[0]),log_Arrhenius[1],1.987*log_Arrhenius[2]]
					
					new_high = '{:.3E}'.format(Arrhenius_highPressureLimit[0])+"   "+'{:.3f}'.format(Arrhenius_highPressureLimit[1])+"    "+'{:.3E}'.format(Arrhenius_highPressureLimit[2])+"\n"
					
					reaction = front+new_high+rest_of_string
					#print(f"reaction after selection = {reaction}")
			else:
				print("Input kinetic mechansim file is invalid")
		else:
			reaction = reaction
		#print(reaction)
		return reaction
	
	def getNewLowerLimits(self,reaction,rxn,beta,selection,isOpt,delta):
		#print(reaction)
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
				zeta = self.rxnUnsrt[rxn].zeta_matrix
				beta = self.RBeta_Dict[rxn]
				selection = np.asarray(self.RBeta_select_Dict[rxn])
				#print(f"selection for low is : {selection}\n")
				#print(f"reaction before selection: {reaction}\n")
					
				#print("\tBeta (LPL):\t{}\n".format(beta))
				self.nominal_lowPressureLimit = np.array([alpha,n,epsilon])
				log_Arrhenius = np.array([alpha,n,epsilon]) + np.asarray(selection*np.asarray(np.dot(choleskyMatrix,beta)).flatten()).flatten()
				
				#log_Arrhenius = np.array([alpha,n,epsilon]) + np.asarray(selection*np.asarray(np.dot(choleskyMatrix,zeta.T.dot(beta.T).T)).flatten()).flatten()
				#need to convert back to original values
				Arrhenius_lowPressureLimit = [np.exp(log_Arrhenius[0]), log_Arrhenius[1], 8.314E-3*log_Arrhenius[2]]
				reaction = reaction.replace(aString,'{:.3E}'.format(Arrhenius_lowPressureLimit[0]))
				reaction = reaction.replace(nString,'{:.3f}'.format(Arrhenius_lowPressureLimit[1]))
				reaction = reaction.replace(eString,'{:.3f}'.format(Arrhenius_lowPressureLimit[2]))
			elif self.fileType == "chemkin":
			
				if self.chemtag[rxn] == 3:
					w = re.compile(r'(.*?\<\=\>.*? \s*?)((-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?)((((LOW)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?\/\s*?.*?\n?\s*?))((TROE)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?\/?\s*?\n?))',re.DOTALL|re.IGNORECASE)
					match1 = re.search(w,reaction)
					front = match1.group(1)
					high = match1.group(2)
					low = match1.group(13)
					rest_of_string = match1.group(23)
				#groups
				#reaction = 1
				#High = 2
				#Low + Troe = 12
				#Low = 13
				#TROE = 23
			
				elif self.chemtag[rxn] == 4:					
					w = re.compile(r'(.*?\<\=\>.*? \s*?)((-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?)(((LOW)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?\/\s*?.*?\n?\s*?)((TROE\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d*\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d*\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d*\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\/?\s*?\/\s*?\n?))((.*\n?\s*?.*)?))',re.DOTALL|re.IGNORECASE)
					match1 = re.search(w,reaction)
					front = match1.group(1)
					high = match1.group(2)
					low = match1.group(13)
					rest_of_string = match1.group(25)+match1.group(38) #TROE+m
				#groups
				#reaction = 1
				#High = 2
				#Low+TROE+m = 12
				#Low = 13
				#Troe = 25
				#m  = 38
				elif self.chemtag[rxn] == 5:
					w = re.compile(r'(.*?\<\=\>.*? \s*?)((-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?)(((LOW)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?\/\s*?.*?\n)(\s*?\/?\s*?\n?.*))',re.DOTALL|re.IGNORECASE)
					match1 = re.search(w,reaction)
					front = match1.group(1)
					high = match1.group(2)
					low = match1.group(13)
					rest_of_string = match1.group(24)
					#groups
					#reaction = 1
					#High = 2
					#Low+m = 12
					#Low = 13
					#m = 24
			
				if "factor" in self.rxnUnsrt[rxn].perturbation_type:
					pattern = re.compile(r'(LOW\s*?\/\s*?)(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))(\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))(\s*?\/\s*?.*?\n?))',re.DOTALL|re.IGNORECASE)
					
					match = re.search(pattern,low)
					x1 = match.group(1)
					A = match.group(2)
					n = match.group(6)
					Ea = match.group(9)
					x2 = match.group(5)
					#A = 2
					#n = 5
					#Ea = 9
					aString = match.group(2)
					a = float(aString)
					alpha = np.log(a)
					
					choleskyMatrix = self.rxnUnsrt[rxn].perturb_factor
					zeta = self.rxnUnsrt[rxn].zeta_matrix
					beta = self.RBeta_Dict[rxn]
					selection = np.asarray(self.RBeta_select_Dict[rxn])
					#print(f"selection for low is : {selection}\n")
					#print(f"reaction before selection: {reaction}\n")

					nString = match.group(6)
					n = float(nString)
											
					eString = match.group(9)
					e = float(eString)
					epsilon = e/1.987
					
					self.nominal_highPressureLimit = np.array([alpha,n,epsilon])
					if self.design == "B1" or self.design == "B1-Mix" or self.design == "B1-factorial" or self.design == "B1-MC" or  self.design == "SAMAP":
						log_Arrhenius = np.array([alpha])+np.asarray(selection*np.asarray(np.dot(choleskyMatrix,beta)).flatten()).flatten()
					elif self.design == "A1":
						log_Arrhenius = np.array([alpha])+np.asarray(selection*np.asarray(np.dot(choleskyMatrix,beta)).flatten()).flatten()
					elif self.design == "Monte-Carlo-All":
						log_Arrhenius = np.array([alpha])+np.asarray(selection*np.asarray(np.dot(choleskyMatrix,zeta.T.dot(beta.T).T)).flatten()).flatten()
					Arrhenius_lowPressureLimit = [np.exp(log_Arrhenius[0]),n,1.987*epsilon]
									
					new_low = x1+'{:.3E}'.format(Arrhenius_lowPressureLimit[0])+x2
					
					reaction = front+high+new_low+rest_of_string
					#print(f"reaction after selection = {reaction}\n")
				else:
					pattern = re.compile(r'(LOW\s*?\/\s*?)(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))(\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))(\s*?\/\s*?.*?\n?))',re.DOTALL|re.IGNORECASE)
					match = re.search(pattern,low)
					
					#print(reaction)
					
					#print(match.group())
					x1 = match.group(1)
					x3 = match.group(12)		
					aString = match.group(2)
					a = float(aString)
					alpha = np.log(a)
					
					nString = match.group(6)
					n = float(nString)
					#print(n)						
					eString = match.group(9)
					e = float(eString)
					epsilon = e/1.987
					
					choleskyMatrix = self.rxnUnsrt[rxn].cholskyDeCorrelateMat
					zeta = self.rxnUnsrt[rxn].zeta_matrix
					beta = self.RBeta_Dict[rxn]
					selection = np.asarray(self.RBeta_select_Dict[rxn])
					#print(f"selection for low is : {selection}")
					#print(f"reaction before selection: {reaction}")
					self.nominal_highPressureLimit = np.array([alpha,n,epsilon])
					if self.design == "B1" or self.design == "B1-Mix" or self.design == "B1-factorial" or self.design == "B1-MC" or  self.design == "SAMAP":
						log_Arrhenius = np.array([alpha,n,epsilon])+ np.asarray(selection*np.asarray(np.dot(choleskyMatrix,beta)).flatten()).flatten()
					elif self.design == "A1":
						log_Arrhenius = np.array([alpha,n,epsilon])+ np.asarray(selection*np.asarray(np.dot(choleskyMatrix,beta)).flatten()).flatten()
						
					elif self.design == "Monte-Carlo-All":
						log_Arrhenius = np.array([alpha,n,epsilon])+ np.asarray(selection*np.asarray(np.dot(choleskyMatrix,zeta.T.dot(beta.T).T)).flatten()).flatten()
					Arrhenius_lowPressureLimit = [np.exp(log_Arrhenius[0]),log_Arrhenius[1],1.987*log_Arrhenius[2]]
					
					new_low = x1+'{:.3E}'.format(Arrhenius_lowPressureLimit[0])+"   "+'{:.3f}'.format(Arrhenius_lowPressureLimit[1])+"    "+'{:.3E}'.format(Arrhenius_lowPressureLimit[2])+x3
					
					reaction = front+high+new_low+rest_of_string
					#print(f"reaction after selection: {reaction}")
			else:
				print("Input kinetic mechansim file is invalid")
		else:
			reaction = reaction
		#print(reaction)
		return reaction
	
	
	def getNewFallOffCurve(self,reaction,target,rxn,beta,selection,isOpt,delta):
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
				zeta = self.rxnUnsrt[rxn].zeta_matrix
		
				beta = self.RBeta_Dict[rxn]
				#print("\tBeta (FoC):\t{}\n".format(beta))
				self.nominal_foc = np.array([fcb])
			
				self.perturbedFoC = np.array([fcb]) + np.asarray(selection*np.asarray(np.dot(choleskyMatrix,zeta.T.dot(beta.T).T)).flatten()).flatten()
			
				#fcta_ = -T/np.log(self.perturbedFoC[0])
				fcb_ = self.perturbedFoC[0]
				#fca_ = self.perturbedFoC[0]
				fca_ = 1-fcb_
			
				reaction = reaction.replace(fca_String,'{:.3E}'.format(fca_))
				#reaction = reaction.replace(fcta_String,'{:.3E}'.format(fcta_))
				reaction = reaction.replace(fcb_String,'{:.3E}'.format(fcb_))
			
		
			elif self.fileType == "chemkin":
				
				if self.chemtag[rxn] == 3:
					w = re.compile(r'(.*?\<\=\>.*? \s*?)((-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?)((((LOW)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?\/\s*?.*?\n?\s*?))((TROE)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?\/?\s*?\n?))',re.DOTALL|re.IGNORECASE)
				#groups
				#reaction = 1
				#High = 2
				#Low + Troe = 12
				#Low = 13
				#TROE = 23
								
				elif self.chemtag[rxn] == 4:					
					w = re.compile(r'(.*?\<\=\>.*? \s*?)((-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?)(((LOW)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?\/\s*?.*?\n?\s*?)((TROE\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d*\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d*\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d*\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\/?\s*?\/\s*?\n?))((.*\n?\s*?.*)?))',re.DOTALL|re.IGNORECASE)
				#groups
				#reaction = 1
				#High = 2
				#Low+TROE+m = 12
				#Low = 13
				#Troe = 25
				#m  = 38

				match1 = re.search(w, reaction)
				fca_String = match1.group(22)
				fca = float(fca_String)
				reaction = reaction.replace(fca_String,'{:.3E}'.format(fca_))
				choleskyMatrix = self.rxnUnsrt[rxn].cholskyDeCorrelateMat
				zeta = self.rxnUnsrt[rxn].zeta_matrix
		
				beta = self.RBeta_Dict[rxn]
				#print("\tBeta (FoC):\t{}\n".format(beta))
				self.nominal_foc = np.array([fca])
			
				self.perturbedFoC = self.nominal + np.asarray(selection*np.asarray(np.dot(choleskyMatrix,zeta.T.dot(beta.T).T)).flatten()).flatten()
			
				fca_ = self.perturbedFoC[0]
				reaction = reaction.replace(fca_String,'{:.3f}'.format(fca_))
			else: 
				print("Invalid kinetic mechanism file")
		
		else:
			reaction = reaction
		return reaction
	
	
	
	def getNewRxn(self,reaction,rxn,beta,selection,isOpt):
		#print(f"----------Entered Mechanism Manipulator ------------\n\n\tReaction\n\n{rxn}\nBeta\n\n{beta}\n\nSelection\n\n\t{selection}\n")
		#print(self.fileType)
		#print(self.chemtag[rxn])
		#print(reaction)
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
			zeta = self.rxnUnsrt[rxn].zeta_matrix
			beta = self.RBeta_Dict[rxn]
			
			if self.design == "Monte-Carlo-All":
				log_Arrhenius = np.array([alpha,n,epsilon]) + np.asarray(selection*np.asarray(np.dot(choleskyMatrix,zeta.T.dot(beta.T).T)).flatten()).flatten()
			#need to convert back to original values
			Arrhenius = [np.exp(log_Arrhenius[0]), log_Arrhenius[1], 8.314E-3*log_Arrhenius[2]]
		
			reaction = reaction.replace(aString,'{:.3E}'.format(Arrhenius[0]))
			reaction = reaction.replace(nString,'{:.3f}'.format(Arrhenius[1]))
			reaction = reaction.replace(eString,'{:.3E}'.format(Arrhenius[2]))
		
		elif self.fileType == "chemkin":
			w = re.compile(r'(.*?\<\=\>.*? \s*?)((-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))(\s*?.*?\n?\s*?))',re.DOTALL|re.IGNORECASE)
			#reaction = 1
			# body = 2
			#A = 3
			#n = 6
			#Ea = 9
			#end = 12
			
			match1 = re.search(w, reaction)
			
			front = match1.group(1)
			
			aString = match1.group(3)
			a = float(aString)
			alpha = np.log(a)

			nString = match1.group(6)
			n = float(nString)

			eString = match1.group(9)
			e = float(eString)
			epsilon = e/1.987
			
			end = match1.group(12)
			
			if "factor" in self.rxnUnsrt[rxn].perturbation_type:
								
				self.nominal_highPressureLimit = np.array([alpha,n,epsilon])
				choleskyMatrix = self.rxnUnsrt[rxn].perturb_factor
				zeta = self.rxnUnsrt[rxn].zeta_matrix
				beta = self.RBeta_Dict[rxn]
				#print(f"selection is: {selection}")
				#print(f"reaction before selection: {reaction}")
				if self.design == "B1" or self.design == "B1-Mix" or self.design == "B1-factorial" or self.design == "B1-MC" or  self.design == "SAMAP":	
					log_Arrhenius = np.array([alpha])+np.asarray(selection*np.asarray(np.dot(choleskyMatrix,beta[0]).flatten())).flatten()
				elif self.design == "A1":	
					log_Arrhenius = np.array([alpha])+np.asarray(selection*np.asarray(np.dot(choleskyMatrix,beta[0]).flatten())).flatten()
				elif self.design == "Monte-Carlo-All":
					log_Arrhenius = np.array([alpha])+np.asarray(selection*np.asarray(np.dot(choleskyMatrix,beta[0]*zeta).flatten())).flatten()
				
				Arrhenius_lowPressureLimit = [np.exp(log_Arrhenius[0]),n,1.987*epsilon]
				new_A = Arrhenius_lowPressureLimit[0]
				new_n = Arrhenius_lowPressureLimit[1]
				new_Ea = Arrhenius_lowPressureLimit[2]
							 
			else:				
				choleskyMatrix = self.rxnUnsrt[rxn].cholskyDeCorrelateMat
				zeta = self.rxnUnsrt[rxn].zeta_matrix
				#print(type(zeta))
				beta = self.RBeta_Dict[rxn]
				#print(type(beta))
				#print(f"selection is: {selection}\n")
				#print(f"reaction before selection: {reaction}\n")

				p1 = np.asarray(np.dot(choleskyMatrix,zeta.T.dot(beta.T).T)).flatten()
				##print(f"p1 = {p1}")
				#print(selection*p1)
				if self.design == "B1" or self.design == "B1-Mix" or self.design == "B1-factorial" or self.design == "B1-MC" or  self.design == "SAMAP":
					log_Arrhenius = np.array([alpha,n,epsilon])+ np.asarray(selection*np.asarray(np.dot(choleskyMatrix,beta)).flatten()).flatten()
				elif self.design == "A1":
					log_Arrhenius = np.array([alpha,n,epsilon])+ np.asarray(selection*np.asarray(np.dot(choleskyMatrix,beta)).flatten()).flatten()
				elif self.design == "Monte-Carlo-All":
					log_Arrhenius = np.array([alpha,n,epsilon])+ np.asarray(selection*np.asarray(np.dot(choleskyMatrix,zeta.T.dot(beta.T).T)).flatten()).flatten()
				Arrhenius_lowPressureLimit = [np.exp(log_Arrhenius[0]),log_Arrhenius[1],1.987*log_Arrhenius[2]]
				
				new_A = Arrhenius_lowPressureLimit[0]
				new_n = Arrhenius_lowPressureLimit[1]
				new_Ea = Arrhenius_lowPressureLimit[2]
			
			reaction = front+'{:.3E}'.format(new_A)+"   "+'{:.3f}'.format(new_n)+"    "+'{:.3E}'.format(new_Ea)+"\n"
			#print(f"reaction after selection: {reaction}\n")
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
		#print(line)
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
			#for ind,m in enumerate(mList):
				
			key = []
			param=[]		
			
			if self.chemtagCollision[rxn] == 2:
				w = re.compile(r'(.*?\<\=\>.*? \s*?)((-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?)(\s*?.*?\/\s*?.*\/?\s*?.*?\n?)',re.DOTALL|re.IGNORECASE)
					#groups
					#reaction = 1
					#High = 2
					#m = 12
				match1 = re.search(w,line)
				String = match1.group(11)	
				line_content = match1.group(11).split('/')
				
			
			elif self.chemtagCollision[rxn] == 4:
				w = re.compile(r'(.*?\<\=\>.*? \s*?)((-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?)(((LOW)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?\/\s*?.*?\n?\s*?)((TROE\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d*\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d*\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d*\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\/?\s*?\/\s*?\n?))((.*\n?\s*?.*)?))',re.DOTALL|re.IGNORECASE)
				#groups
				#reaction = 1
				#High = 2
				#Low+TROE+m = 12
				#Low = 13
				#Troe = 25
				#m  = 38
				match1 = re.search(w,line)
				String = match1.group(34)	
				line_content = match1.group(34).split('/')
			elif self.chemtagCollision[rxn] == 5:
				w = re.compile(r'(.*?\<\=\>.*? \s*?)((-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?)(((LOW)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?\/\s*?.*?\n)(\s*?\/?\s*?\n?.*))',re.DOTALL|re.IGNORECASE)
					#groups
					#reaction = 1
					#High = 2
					#Low+m = 12
					#Low = 13
					#m = 24
				match1 = re.search(w,line)
				String = match1.group(24)	
				line_content = match1.group(24).split('/')
			#print(line_content)
			index = []
			newCollisionEff = []
			for ind,m in enumerate(mList):
				for element in line_content:     
					if element.strip() == m:
						
						choleskyMatrix = self.tbdUnsrt[rxn].cholskyDeCorrelateMat[m]
						zeta = self.tbdUnsrt[rxn].zeta[m]
						index.append(line_content.index(element)+1)
						val = float(line_content[line_content.index(element)+1])
						#print(type(beta))
						#print(beta[str(m)])
						sigma =abs( np.asarray(np.dot(choleskyMatrix,beta[m]*zeta)).flatten()[0])
						#print(sigma[0])
						k = math.copysign(1,beta[m])
						newCollisionEff.append(str(val*(1+sigma*k)))
					else:
						continue
			
			#print(newCollisionEff)
			#print(line_content)
			#print(mList)
			for ind,val in enumerate(index):
				line_content[val] =  newCollisionEff[ind]
			#print(param)		
			#print(line_content)
			newString = "/".join(line_content)
		 
			newLine = line.replace(String,newString)
			#print(line)
			#print(newLine)
		
		return newLine
		
	
	
	
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
		#print("Beta is {}".format(beta))
		unsrt = []	
		
		for mol_list in self.molecule:
			#print(mol_list)
			for mol in mol_list:
				unsrt.append(self.tbdUnsrt[rxn].uncertainties[mol])
		
			newLine = self.getNewCollisionEff(line,mol_list,beta,unsrt,rxn)
		newLine = newLine+'\n'
		#print("Third body collision efficiencies after perturbation:\n\n{}\n".format(newLine))		 
		return newLine    
	
            
        
        
