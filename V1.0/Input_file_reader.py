import numpy as np
import re
import os
from scipy.optimize import minimize
class MechParsing():
	def __init__ (self, Mechfile):
		self.rxn_index = [] 
		self.rxn_class = []
		self.rxn = [] 
		self.A1 = []
		self.n1 = []
		self.Ea1 = []
		self.A2 = []
		self.n2 = []
		self.Ea2 = []
		self.fca = []
		self.fcta = []
		self.fcb = []
		self.fctb = []
		self.fcc = [] 
		self.fctc = []
		self.rxn_PDept = []
		self.mechFile = Mechfile
		
		mech = open(self.mechFile).read()
		f = open(self.mechFile).readlines()
		#print(f)
		temp = []
		for i in f:
			if i.startswith("!"):
				continue
			else:
				temp.append(i)
		f = temp
		#use this for validation purpose only
		#print("----------------------------\n")
		#print(f)
		
		#print("----------------------------\n")
		#f = ("\n").join(f)
		#temp_file = open("mechanism.temp","w")
		#temp_file.write(f)
		#temp_file.close()
		#print(f)
		#print("==============================\n")
		#print(mech)
				
		if "1f:" in mech and "let allowed atoms be" in mech and "->" in mech:
			print("read the file")
			self.fileType = "FlameMaster"
		elif "REACTIONS" in mech:
			self.fileType = "chemkin"
		#print(self.fileType)
		else:
			print("File Type Reading problem")
		self.mechLine = ""
		self.tbcLine = ""
		self.Reactions = {}
		self.mLine = []
		self.mList = []
		if self.fileType == "FlameMaster":
			for i in f:
				if ":" in i:		
					if  "}" in i:
						self.mechLine += i
						pattern = re.compile(r'(\d*\w):\s*(.*)\s*->\s*(.*)\s*?\{\s*?A\s*?=\s*?(\S*E\S*)?\s*n\s*?=\s*(\S*)?\s*E\s*?=\s*(\S*)?\s*.*\}}?', re.DOTALL | re.IGNORECASE)
						self.match = re.search(pattern,self.mechLine)
						self.rxn_index.append(self.match.group(1).strip())
						reactant = (self.match.group(2))
						product = (self.match.group(3))
						self.A1.append(float(self.match.group(4)))
						self.n1.append(float(self.match.group(5)))
						self.Ea1.append(float(self.match.group(6)))
						self.A2.append(float(0.0))
						self.n2.append(float(0.0))
						self.Ea2.append(float(0.0))
						self.fca.append(float(0.0))
						self.fcta.append(float(0.0))
						self.fcb.append(float(0.0))
						self.fctb.append(float(0.0))
						self.fcc.append(float(0.0))
						self.fctc.append(float(0.0))
						self.Reactions[i.split(" ")[0]] = self.mechLine
						self.Reactions["reactant"] = reactant
						self.Reactions["product"] = product
						#print(mechLine)
						self.mechLine = ""
					
					else:
						j = f.index(i)
						while "}" not in self.mechLine:
							self.mechLine += f[j]
							j = j+1
					
						pattern = re.compile(r'(\d*\w):\s*(.*)\s*->\s*(.*)\s*?\{\s*?Ai\s*?=\s*?(\S*E\S*)?\s*ni\s*?=\s*(\S*)?\s*Ei\s*?=\s*(\S*)?\s*A\s*?=\s*?(\S*E\S*)?\s*n\s*?=\s*(\S*)?\s*E\s*?=\s*(\S*)?\s*fca\s*?=\s*?(\S*E\S*)?\s*?fcta\s*?=\s*?(\S*E\S*)?\s*?fcb\s*?=\s*?(\S*E\S*)?\s*?fctb\s*?=\s*?(\S*E\S*)?\s*?fcc\s*?=\s*?(\S*)?\s*?fctc\s*?=\s*?(\S*E\S*).*', re.DOTALL | re.IGNORECASE)
																
						self.match = re.search(pattern,self.mechLine)
						self.rxn_index.append(self.match.group(1).strip())
						reactant = (self.match.group(2))
						product = (self.match.group(3))
						self.A1.append(float(self.match.group(4)))
						self.n1.append(float(self.match.group(5)))
						self.Ea1.append(float(self.match.group(6)))
						self.A2.append(float(self.match.group(7)))
						self.n2.append(float(self.match.group(8)))
						self.Ea2.append(float(self.match.group(9)))					
						self.fca.append(float(self.match.group(10)))
						self.fcta.append(float(self.match.group(11)))
						self.fcb.append(float(self.match.group(12)))					
						self.fctb.append(float(self.match.group(13)))
						self.fcc.append(float(self.match.group(14)))
						self.fctc.append(float(self.match.group(15)))					
						self.Reactions[i.split(" ")[0]] = self.mechLine
						self.Reactions["reactant"] = reactant
						self.Reactions["product"] = product
						#print(mechLine)
						self.mechLine = "" 	
	
		
			patternThirdBody = re.compile(r'\w*?\s*?(M*?[0-9]*?)=.*?(\S*?\s*?\[.*\].).*?')
			string = open(self.mechFile).read()
			lines = string.split("\n")
			self.thirdBody = {}

			for i in lines:
				match = re.search(patternThirdBody,i)
				if match != None:
					per_line = match.group(0).split("=")
					self.thirdBodyIndex = per_line[0]
					self.mList.append(per_line[0])
					self.thirdBody[self.thirdBodyIndex] = per_line[1]
		
		
			#print(self.thirdBody)
		elif self.fileType == "chemkin":
			rxnLine = re.compile(r'(.*?<=>.*?\s)\s*?.*?')
			self.rxnList = []
			self.thirdBody = {}
			for i in f:
				if "<=>" in i:
					#print(i)
					match = re.search(rxnLine,i)
					self.rxnList.append(match.group(1).strip())
			#print(self.rxnList)
			for r in self.rxnList:
				#print(r)
				for i in f:
					if i.startswith("!"):
						continue
					else:
						#print(i)
						if r in i: 	 
							 #print(line)
							if "DUPLICATE" in f[f.index(i)+1] and "DUPLICATE" in f[f.index(i)+3]:
								self.chemtag = 1
								start = f.index(i)
								self.mechline = ""
								for j in range(4):
									self.mechLine +=f[start+j]
							
								w = re.compile(r'(.*?\<\=\>.*?) \s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?DUPLICATE\s*?\n?(.*?\<\=\>.*?) \s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?DUPLICATE\s*?\n?',re.DOTALL|re.IGNORECASE)
								
								#print("_______________________________\n")
								#print("Duplicate\n")
								#print(self.mechLine)
								#print("_______________________________\n")
								self.match = re.search(w,self.mechLine)
								#print(self.match.group())
								self.rxn_index.append(self.match.group(1).strip())
		
								self.A1.append(float(self.match.group(2)))
								self.n1.append(float(self.match.group(5)))
								self.Ea1.append(float(self.match.group(8)))
								self.A2.append(self.match.group(12))
								self.n2.append(self.match.group(15))
								self.Ea2.append(self.match.group(18))
								self.fca.append(float(0.0))
								self.fcta.append(float(0.0))
								self.fcb.append(float(0.0))
								self.fctb.append(float(0.0))
								self.fcc.append(float(0.0))
								self.fctc.append(float(0.0))
								#print(self.mechLine)
								self.Reactions[str(r)] = self.mechLine
								self.mechLine = "" 
								continue
							elif "LOW" not in f[f.index(i)+1] and "/" in f[f.index(i)+1]:
								#print("no LOW TROE\n")
								self.chemtag = 2
								start = f.index(i)  
								for j in range(2):
									self.mechLine +=f[start+j]
							
								w = re.compile(r'(.*?\<\=\>.*?) \s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?(.*?\/\s*?.*)\/\s*?.*?\n?',re.DOTALL|re.IGNORECASE)

								#print("_______________________________\n")
								#print(self.mechLine)
								#print("_______________________________\n")
								self.match = re.search(w,self.mechLine )
								#print(self.match.group())
								self.rxn_index.append(self.match.group(1).strip())
							
								self.A1.append(float(self.match.group(2)))
								self.n1.append(float(self.match.group(5)))
								self.Ea1.append(float(self.match.group(8)))
								self.A2.append(float(0.0))
								self.n2.append(float(0.0))
								self.Ea2.append(float(0.0))
								self.fca.append(float(0.0))
								self.fcta.append(float(0.0))
								self.fcb.append(float(0.0))
								self.fctb.append(float(0.0))
								self.fcc.append(float(0.0))
								self.fctc.append(float(0.0))
							
								self.mList.append(self.match.group(11))
								self.thirdBody[str(r)] = self.match.group(11)
								self.Reactions[str(r)] = self.mechLine
								self.mechLine = "" 
								continue
							
							elif "LOW" in f[f.index(i)+1] and "TROE" not in f[f.index(i)+2] and "/" in f[f.index(i)+2]:
								#print("no LOW TROE\n")
								self.chemtag = 5
								start = f.index(i)  
								for j in range(3):
									self.mechLine +=f[start+j]
							
								w = re.compile(r'(.*?\<\=\>.*?) \s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?(LOW)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?\/\s*?.*?\n?\s*?(.*?\/\s*?.*)\/\s*?.*?\n?',re.DOTALL|re.IGNORECASE)

								#print("_______________________________\n")
								#print(self.mechLine)
								#print("_______________________________\n")
								self.match = re.search(w,self.mechLine )
								#print(self.match.group())
								self.rxn_index.append(self.match.group(1).strip())
							
								self.A1.append(float(self.match.group(2)))
								self.n1.append(float(self.match.group(5)))
								self.Ea1.append(float(self.match.group(8)))
								self.A2.append(float(self.match.group(12)))
								self.n2.append(float(self.match.group(15)))
								self.Ea2.append(float(self.match.group(18)))
								self.fca.append(float(0.0))
								self.fcta.append(float(0.0))
								self.fcb.append(float(0.0))
								self.fctb.append(float(0.0))
								self.fcc.append(float(0.0))
								self.fctc.append(float(0.0))
							
								self.mList.append(self.match.group(21))
								self.thirdBody[str(r)] = self.match.group(21)
								self.Reactions[str(r)] = self.mechLine
								self.mechLine = "" 
								continue
							
							elif "LOW" in f[f.index(i)+1] and "TROE" in f[f.index(i)+2] and "/" not in f[f.index(i)+3]:
								#print("LOW TROE no m\n")
								self.chemtag = 3
								start = f.index(i)  
								for j in range(3):
									self.mechLine +=f[start+j]
							
								w = re.compile(r'(.*?\<\=\>.*?) \s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?(LOW)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?\/\s*?.*?\n?\s*?(TROE)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?\/\s*?\n?',re.DOTALL|re.IGNORECASE)
								#print("_______________________________\n")
								#print(self.mechLine)
								#print("_______________________________\n")
								self.match = re.search(w,self.mechLine )
								#print(self.match.group())
								self.rxn_index.append(self.match.group(1).strip())
							
								self.A1.append(float(self.match.group(2)))
								self.n1.append(float(self.match.group(5)))
								self.Ea1.append(float(self.match.group(8)))
							
								self.A2.append(self.match.group(12))
								self.n2.append(self.match.group(15))
								self.Ea2.append(self.match.group(18))
								self.fca.append(1.0-float(self.match.group(22)))
								self.fcta.append(self.match.group(25))
								self.fcb.append(self.match.group(22))
								self.fctb.append(self.match.group(28))
								self.fcc.append(float(1.0))
								self.fctc.append(self.match.group(31))				
							
							
								self.Reactions[str(r)] = self.mechLine
								self.mechLine = "" 
								continue
							elif "LOW" in f[f.index(i)+1] and "TROE" in f[f.index(i)+2] and "/" in f[f.index(i)+3]:
								#print("LOW TROE\n")
								self.chemtag = 4
								start = f.index(i)  
								for j in range(4):
									self.mechLine += f[start+j]
							#	print("_______________________________\n")
							#	print(self.mechLine)
							#	print("_______________________________\n")
								#print(self.mechLine)
								w = re.compile(r'(.*?\<\=\>.*?) \s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?.*?\n?\s*?(LOW)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?\/\s*?.*?\n?\s*?(TROE)\s*?\/\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?\/\s*?\n(.*)',re.DOTALL|re.IGNORECASE)
						
							
								self.match = re.search(w,self.mechLine )
								#print(self.match.group())
								self.rxn_index.append(self.match.group(1).strip())
							
								self.A1.append(float(self.match.group(2)))
								self.n1.append(float(self.match.group(5)))
								self.Ea1.append(float(self.match.group(8)))
								self.A2.append(self.match.group(12))
								self.n2.append(self.match.group(15))
								self.Ea2.append(self.match.group(18))
								self.fca.append((1.0-float(self.match.group(22))))
								self.fcta.append(float(self.match.group(25)))
								self.fcb.append(float(self.match.group(22)))
								self.fctb.append(float(self.match.group(28)))
								self.fcc.append(float(1.0))
								self.fctc.append(float(self.match.group(31)))				
							
								#print(match1)
								self.mList.append(self.match.group(34).strip())
								self.thirdBody[str(r)] = self.match.group(34).strip()
								self.Reactions[str(r)] = self.mechLine
								self.mechLine = "" 
								continue
							else:

								self.chemtag = 0
								self.mechLine += i
																
								#match A=group 2, n = group5, Ea=group 8
								w = re.compile(r'(.*?\<\=\>.*? )\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+))\s*?(-?(\d\.?\d*[Ee][+\-]?\d+|(\d+\.\d*|\d*\.\d+)|\d+)).*?\n?',re.DOTALL | re.IGNORECASE)
								#print("_______________________________\n")
								#print(self.mechLine)
								#print("_______________________________\n")
								self.match = re.search(w,self.mechLine)
								#print(self.match.group(2))
								#print(self.match.group(5))
								#print(self.match.group(8))					
								self.rxn_index.append(self.match.group(1).strip())
							
								self.A1.append(float(self.match.group(2)))
								self.n1.append(float(self.match.group(5)))
								self.Ea1.append(float(self.match.group(8)))
								self.A2.append(float(0.0))
								self.n2.append(float(0.0))
								self.Ea2.append(float(0.0))
								self.fca.append(float(0.0))
								self.fcta.append(float(0.0))
								self.fcb.append(float(0.0))
								self.fctb.append(float(0.0))
								self.fcc.append(float(0.0))
								self.fctc.append(float(0.0))
								self.mechLine = ""
								continue
		
						
	def getThirdBodyCollisionEff(self,name, species):
		if self.fileType == "FlameMaster":
			key = []
			value = []
			index = name.split("-")[0]
		
			species = name.split("-")[1]
			for i in self.mList:
				if index in i:
					line = self.thirdBody[i]
					pattern = re.compile(r'(let?\s*\S*\s*?)=.*?(\S*?\s*?\[.*\].).*?')	
					match=re.search(pattern,line.split("\n")[0])
					speciesList= line.split('+')
					for j in speciesList:
						key_and_value = self.filter_list(j.split(" ")) 
						if species in key_and_value:
							key.append(key_and_value[1])
							value.append(key_and_value[0])
			
				else:
					continue
			data = float(value[0])	
					#print(i)
					#print(index)
				
			#print("Key = {}, and Value = {}\n".format(key,value))
		elif self.fileType == "chemkin":
			 data = {}
			 #print(name)
			 #print(self.thirdBody)
			 #print(self.thirdBody)
			 #print(self.rxnList)
			 line = self.thirdBody[name]
			 listTBD = line.split("/") 
			 listTBD = self.filter_list(listTBD)
			 #print(dict(listTBD))
			 count = 0
			 while count<len(listTBD):
			 	#print(count)
			 	key = count
			 	value = count+1
			 	data[listTBD[key].strip()] = float(listTBD[value])
			 	
			 	count +=2
			
			# print(data)	
		return data			
	
	def filter_list(self,List):
		temp = []
		for i in List:
			if i != "":
				temp.append(i)
		return temp
	
	def getRxnData(self,rxn):
		#print(rxn)
		if rxn in self.rxn_index:
			k = self.rxn_index.index(rxn)
			
		else:
			print(rxn)
			print(self.rxn_index)
		return np.array([self.A1[k], self.n1[k], self.Ea1[k]])
	def getKappa(self,rxn):
		#print(rxn)
		
		if rxn in self.rxn_index:
			k = self.rxn_index.index(rxn)
			if self.fileType == "FlameMaster":
				R = 8.314E-03
			else:
				R = 1.987
		else:
			print(rxn)
			print(self.rxn_index)
		return np.array([np.log(self.A1[k]), self.n1[k], (self.Ea1[k]/R)])
	def getRxnData_HPL(self,rxn):
		if rxn in self.rxn_index:
			k = self.rxn_index.index(rxn)
			
		else:
			print(rxn)
			print(self.rxn_index)
		return np.array([self.A1[k], self.n1[k], self.Ea1[k]])
	def getRxnData_LPL(self,rxn):
		if rxn in self.rxn_index:
			k = self.rxn_index.index(rxn)
			
		else:
			print(rxn)
			print(self.rxn_index)
		return np.array([np.log(self.A2[k]), self.n2[k], self.Ea2[k]/8.314E-3])
	def getFocData(self,name):
		foc = name.split("-")[0]
		if foc in self.rxn_index:
			k = self.rxn_index.index(foc)
			data = np.array([self.fca[k], self.fcta[k], self.fcb[k], self.fctb[k], self.fcc[k], self.fctc[k]])
			#print(data)
		else:
			print(foc)
			print(self.rxn_index)
		return data
	def getRxnData_p(self):
		return self.A1,self.n1,self.Ea1,self.A2,self.n2,self.Ea2,self.fca,self.fcta,self.fcb,self.fctb,self.fcc,self.fctc
	
	def getBranchingParameters(self,Branching_rxns,A_list,n_list,Ea_list,isOpt):
		kDash = []
		A = []
		n = []
		Ea = []
		self.BranchingDict_unOpt = {}
		self.BranchingDict_Opt = {}
		T_ = np.linspace(500,2500,50)
		self.ratio = []
		for j in T_:
			k = 0
			k_temp = []
			r_temp = []
			k_Dash = []
			for i,rxn in enumerate(Branching_rxns):
			    rate = A_list[i]*j**n_list[i]*np.exp(-Ea_list[i]/(j*8.314E-3))
			    k+=rate
			    k_temp.append(rate)
			for i in k_temp:
			    r = i/k
			    r_temp.append(r)
			    k_Dash.append(r*(k))
			kDash.append(k_Dash)
			self.ratio.append(r_temp)
		
		if isOpt =="True":
			self.BranchingDict_Opt[str(Branching_rxns)] = self.ratio
		
		else:
			self.BranchingDict_unOpt[str(Branching_rxns)] = self.ratio
		
		for i in range(len(Branching_rxns)):
			A_,n_,Ea_ = Input_file_reader.Curvefit(np.asarray(kDash)[:,i],T_).getBranchingCurveFit()
			A.append(A_)
			n.append(n_)
			Ea.append(Ea_)
		return A,n,Ea
	
	
	def getRatioData(self,Branching_rxns,key_type):
		kDash = []
		A_list = []
		n_list = []
		Ea_list = []
		for i in Branching_rxns:
			a,n,e = self.getRxnData(i)
			A_list.append(a)
			n_list.append(n)
			Ea_list.append(e)
		T_ = np.linspace(500,2500,50)
		ratio = []
		for j in T_:
			k = 0
			k_temp = []
			r_temp = []
			for i,rxn in enumerate(Branching_rxns):
			    rate = A_list[i]*j**n_list[i]*np.exp(-Ea_list[i]/(j*8.314E-3))
			    k+=rate
			    k_temp.append(rate)
			for i in k_temp:
			    r = i/k
			    r_temp.append(r)
			ratio.append(r_temp)
		
		ratioFile = open("./BranchRatio/"+key_type+"_"+str(Branching_rxns),"w")
		string = ""
		for i in np.asarray(ratio):
			for j in range(len(i)):
				string+="{}\t".format(i[j])
			string+="\n"
		
		ratioFile.write(string)
		ratioFile.close()	
		return ratio

class ThermoParsing():
	def __init__ (self, thermofile): 
		
		self.T_low = []
		self.T_high = []
		self.T_mid = []
		self.others = []
		 
		self.a1 = []
		self.a2 = []
		self.a3 = []
		self.a4 = []
		self.a5 = []
		self.a6 = []
		self.a7 = []
		self.A1 = []
		self.A2 = []
		self.A3 = []
		self.A4 = []
		self.A5 = []
		self.A6 = []
		self.A7 = []
		self.thermoFile = thermofile
		f = open(self.thermoFile).read()
		self.Thermo = {}
		lines = f.split("\n")
		speciesPattern = re.compile(r'([A-Z_0-9]*).*',re.DOTALL|re.IGNORECASE)
		self.speciesList = []
		for i in lines:
			matchSpecies = re.search(speciesPattern,i)
			if matchSpecies.group(1) != "":
				self.speciesList.append(matchSpecies.group(1))
		del(self.speciesList[0],self.speciesList[0],self.speciesList[-1])
		for species in self.speciesList:
			#print(species)
			self.thermoLine = ""
			for i in lines[2:]:
				if i.startswith("!"):
					continue
				else:
					if species == i.split(" ")[0]:
						j = lines.index(i)
						count = 0
						while count<4:
							self.thermoLine+=lines[j+count]
							self.thermoLine+="\n"
							count+=1
						self.Thermo[species] = self.thermoLine
						pattern = re.compile(r'\n*?([\w\d]*)\s*?(.*?G)\s*?(\d+\.?\d+)\s*?(\d+\.?\d+)\s*s*?(\d+\.?\d+)\s*?.*?1\s*?.*?\n?\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?2\n?\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?3\n?\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?(-?\d*\.?\d+[Ee][+\-]?\d+)\s*?4\n?',re.DOTALL|re.IGNORECASE)
						#print(self.thermoLine)
						match = re.search(pattern,self.thermoLine)
						
						self.others.append(match.group(2))
						self.T_low.append(float(match.group(3)))
						self.T_mid.append(float(match.group(5)))
						self.T_high.append(float(match.group(4)))
						self.a1.append(float(match.group(6)))
						self.a2.append(float(match.group(7)))
						self.a3.append(float(match.group(8)))
						self.a4.append(float(match.group(9)))
						self.a5.append(float(match.group(10)))
						self.a6.append(float(match.group(11)))
						self.a7.append(float(match.group(12)))
						self.A1.append(float(match.group(13)))
						self.A2.append(float(match.group(14)))
						self.A3.append(float(match.group(15)))
						self.A4.append(float(match.group(16)))
						self.A5.append(float(match.group(17)))
						self.A6.append(float(match.group(18)))
						self.A7.append(float(match.group(19)))
	
		#print(self.Thermo)
	def getSpeciesList(self):
		return self.speciesList
	def getThermoData(self,species):
		species = species.split("-")[0]
		if species in self.speciesList:
			k = self.speciesList.index(species)
			data = [self.a1[k], self.a2[k], self.a3[k],self.a4[k], self.a5[k], self.a6[k],self.a7[k], self.A1[k], self.A2[k],self.A3[k], self.A4[k], self.A5[k],self.A6[k], self.A7[k]]
		else:
			print(species)
			print(self.speciesList)
		return np.asarray(data)
	
	def getThermoHigh(self,species):
		thermo_dict = {}
		species = species.split("-")[0]
		if species in self.speciesList:
			k = self.speciesList.index(species)
			hcp_data = np.array([self.a1[k], self.a2[k], self.a3[k],self.a4[k], self.a5[k]])
			h_data = np.array([self.a6[k]])
			e_data = np.array([self.a7[k]])
			thermo_dict["Hcp"] = hcp_data
			thermo_dict["h"] = h_data
			thermo_dict["e"] = e_data
		else:
			print("Species Not available")
			print(species)
			print(self.speciesList)
		
		return thermo_dict
	def getThermoLow(self,species):
		thermo_dict = {}
		species = species.split("-")[0]
		if species in self.speciesList:
			k = self.speciesList.index(species)
			data = [self.A1[k], self.A2[k], self.A3[k],self.A4[k], self.A5[k],self.A6[k],self.A7[k]]
			hcp_data = np.array([self.A1[k], self.A2[k], self.A3[k],self.A4[k], self.A5[k]])
			h_data = np.array([self.A6[k]])
			e_data = np.array([self.A7[k]])
			thermo_dict["Hcp"] = hcp_data
			thermo_dict["h"] = h_data
			thermo_dict["e"] = e_data
		else:
			print("Species Not available")
			print(species)
			print(self.speciesList)
		return thermo_dict

	
	def function(self,string,coeff,T):
		if string == "Hcp":
			func = coeff[0]*(T/T) + coeff[1]*T + coeff[2]*T**2  +coeff[3]*T**3 + coeff[4]*T**4
		if string == "h":
			func = coeff[0]*(T/T)
		if string == "e":
			func = coeff[0]*(T/T)
		return func
	def getHCP_low(self,species):
		species = species.split("-")[0]
		if species in self.speciesList:
			k = self.speciesList.index(species)
			data = [self.a1[k], self.a2[k], self.a3[k],self.a4[k], self.a5[k]]
		else:
			print("Species Not available")
			print(species)
			print(self.speciesList)
		return np.asarray(data)
	def getHCP_high(self,species):
		species = species.split("-")[0]
		if species in self.speciesList:
			k = self.speciesList.index(species)
			data = [self.A1[k], self.A2[k], self.A3[k],self.A4[k], self.A5[k]]
		else:
			print("Species Not available")
			print(species)
			print(self.speciesList)
		return np.asarray(data)
	def getEnthalpy_low(self,species):
		species = species.split("-")[0]
		if species in self.speciesList:
			k = self.speciesList.index(species)
			data = [self.a6[k]]
		else:
			print("Species Not available")
			print(species)
			print(self.speciesList)
		return np.asarray(data)
	def getEnthalpy_high(self,species):
		species = species.split("-")[0]
		if species in self.speciesList:
			k = self.speciesList.index(species)
			data = [self.A6[k]]
		else:
			print("Species Not available")
			print(species)
			print(self.speciesList)
		return np.asarray(data)
	def getEntropy_low(self,species):
		species = species.split("-")[0]
		if species in self.speciesList:
			k = self.speciesList.index(species)
			data = [self.a7[k]]
		else:
			print("Species Not available")
			print(species)
			print(self.speciesList)
		return np.asarray(data)
	def getEntropy_high(self,species):
		species = species.split("-")[0]
		if species in self.speciesList:
			k = self.speciesList.index(species)
			data = [self.A7[k]]
		else:
			print("Species Not available")
			print(species)
			print(self.speciesList)
		return np.asarray(data)
class TransportParsing():
	def __init__ (self, transportfile): 
		self.species = []
		self.molIndex = []
		self.mu = []
		self.alpha = []
		self.zRotation = []
		self.epsilon = []
		self.sigma = []
		self.transportFile = transportfile
		f = open(transportfile,'r').read()
		self.Transport = {}
		lines = f.split("\n")
		speciesPattern = re.compile(r'([A-Z_0-9]*).*',re.DOTALL|re.IGNORECASE)
		self.speciesList = []
		for i in lines:
			matchSpecies = re.search(speciesPattern,i)
			if matchSpecies.group(1) != "":
				self.speciesList.append(matchSpecies.group(1))
		for species in self.speciesList:
			self.transportLine = ""
			for i in lines:
				#print(i)
				if "!" not in i:
					if species == i.split(" ")[0]:
						
						self.transportLine+=i
						pattern = re.compile(r'(\w*)\s*(\d)\s*(\d+\.\d*|\d*\.\d+)\s*(\d+\.\d*|\d*\.\d+)\s*(\d+\.\d*|\d*\.\d+)\s*(\d+\.\d*|\d*\.\d+)\s*(\d+\.\d*|\d*\.\d+)',re.DOTALL|re.IGNORECASE)
						
						match = re.search(pattern,i)
						
						self.Transport[species] = self.transportLine
						self.species.append(species)
						self.molIndex.append(float(match.group(2)))
						self.mu.append(float(match.group(5)))
						self.alpha.append(float(match.group(6)))
						self.zRotation.append(match.group(7))
						self.epsilon.append(float(match.group(3)))
						self.sigma.append(float(match.group(4)))
	def getSpeciesList(self):
		return self.speciesList
	def getTransportData(self,species):
		data = {}
		for i in self.species:
			if str(species) in i:
				k = self.species.index(i)
				data["LJe"] = self.epsilon[k]
				data["LJs"] = self.sigma[k]
			else:

				continue
			#print(species)
			#print(self.species)
		return data
	def epsilon(self,species):
		species = species.split("-")[0]
		if species in self.species:
			k = self.species.index(species)
			data = [self.epsilon[k]]
		else:
			print("Species Not available")
			print(species)
			print(self.species)
		return np.asarray(data)
	def alpha(self,species):
		species = species.split("-")[0]
		if species in self.species:
			k = self.species.index(species)
			data = [self.aplha[k]]
		else:
			print("Species Not available")
			print(species)
			print(self.species)
		return np.asarray(data)
	def sigma(self,species):
		species = species.split("-")[0]
		if species in self.species:
			k = self.species.index(species)
			data = [self.sigma[k]]
		else:
			print("Species Not available")
			print(species)
			print(self.species)
		return np.asarray(data)
	def mu(self,species):
		species = species.split("-")[0]
		if species in self.species:
			k = self.species.index(species)
			data = [self.mu[k]]
		else:
			print("Species Not available")
			print(species)
			print(self.species)
		return np.asarray(data)
	def zRotation(self,species):
		species = species.split("-")[0]
		if species in self.species:
			k = self.species.index(species)
			data = [self.zRotation[k]]
		else:
			print("Species Not available")
			print(species)
			print(self.species)
		return np.asarray(data)
class Curvefit():
	def __init__(self,Ydata,Xdata):
		self.A = self.n = self.Ea = self.guess = self.alpha = self.func =self.reaction = self.pressure= None
		self.xdata = Xdata
		self.ydata = Ydata
		
	def branchingObj(self,guess):
		R = 1.987;
		Rdash = 8.314E-3
		a1 = guess[0];
		a2 = guess[1];
		a3 = guess[2];
		k = np.log(self.ydata) - (a1+a2*np.log(self.xdata)-a3/(Rdash*self.xdata));
		self.branchingObj = np.dot(k,k);
		return self.branchingObj

	def fallOffCurveObj(self,guess):
		#(a,fcta,fctb,fctc)
		fcc = 1.0
		fca = 1-guess[0]
		fcb = guess[0]
		fcta = guess[1]
		fctb = guess[2]
		fctc = guess[3]
		T = self.xdata
		foc_function = (1-fca)*np.exp(-T/fcta)+fcb*np.exp(-T/fctb)+fcc*np.exp(-fctc/T)	
		k = self.ydata - foc_function
		self.fallOffCurveObj = np.dot(k,k);
		return self.fallOffCurveObj
	
	
	def getBranchingCurveFit (self):
		#GUESS VALUES
		self.guess = np.array([1000,1000,5000.0])
		#BOUNDS	
		a = (-np.inf,np.inf);
		b = (-3,3);
		bnds = (a,b,a);
		
		self.sol = minimize(self.branchingObj,self.guess);
	
		alpha = self.sol.x[0];
		n = self.sol.x[1];
		Ea = self.sol.x[2];
		A = np.exp(alpha);
		return A, n, Ea
	
	def getFallOffCurveFit(self,foc):
		self.guess = foc
		a = (-np.inf,np.inf);
		bnds = (a,a,a,a)
		self.sol = minimize(self.fallOffCurveObj,self.guess);
		fca = 1-self.sol.x[0];
		fcta = self.sol.x[1];
		fcb = self.sol.x[0];
		fctb = self.sol.x[2];
		fcc  = 1.0
		fctc = self.sol.x[3]
		
		return fca,fcta,fcb,fctb,fcc,fctc
		

