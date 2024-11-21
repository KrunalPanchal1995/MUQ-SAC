try:
    import ruamel_yaml as yaml
except ImportError:
    from ruamel import yaml
import yaml
import os,re

#rxn = ["AC3H5OOH <=> C2H3CHO + H2O"]
#rxn = ["C3H6 + OH <=> C3H6OH1-2:A","C3H6 + OH <=> C3H6OH1-2:B","C3H6 + OH <=> C3H6OH1-2:A","C3H6 + OH <=> C3H6OH1-2:B","C2H + CH3 <=> C3H4-P","C4H10 (+M) <=> 2 C2H5 (+M)","C4H8OOH1-4O2 <=> C4H72-1,4OOH:A","C4H8OOH1-4O2 <=> C4H72-1,4OOH:B","C4H71-3OOH => C2H3CHO + CH3 + OH:B1","C4H71-3OOH => CH3CHO + C2H3 + OH:B2"]

#mechanism = "/home/krithika/Desktop/KineticMechanismOptimization/Mechanism/nHeptane/heptane_632.yaml"

class Parser:
	def __init__(self,MechFile):
		self.file_yaml = open(MechFile,"r").read()
		self.mech = yaml.safe_load(self.file_yaml)
		self.rxnList = self.rxn_list()
		#print(self.mech["phases"][0]["species"])
		#raise AssertionError("Stop")
	
	def rxn_list(self):
		rxnList = []
		for i,ind in enumerate(self.mech["reactions"]):
		    rxnList.append(self.mech["reactions"][i]["equation"])
		#print(f"Reactions are {len(rxnList)}")
		return rxnList
	
	#def populatePlogRxn(self,plogRxn):
		
	def PerturbingReactions(self,activeReaction):
		"""
		LOOK FOR DUPLICATE REACTIONS:
			PRESSURE INDEPENDENT
			PRESSURE DEPENDENT (PLOG)
		
		LOOK FOR BRANCHING REACTIONS
			(have same reactants but different products)
			looks for reaction next to the index and sees
			if the dict has no duplicate flag and has same
			reactant but different products
		"""
		
		TOKKENS = [" <=> "," = "," => "]
		
		PerturbingRxnDict = {}
		for rxn in activeReaction:
			#print(i)
			if ":" in rxn:
				i = rxn.split(":")[0]
				description = rxn.split(":")[1]				
				for tokken in TOKKENS:
					if tokken in i:
						
						reactants = i.split(tokken)[0].strip()
						products = i.split(tokken)[1].strip()
					else:
						continue
					#else:
					#	raise AssertionError(f"New Tokken detected in reaction {i}\n")
			else:
				i = rxn
				#print(i)
				for tokken in TOKKENS:
					if tokken in i:
						reactants = i.split(tokken)[0].strip()
						products = i.split(tokken)[1].strip()
						
					else:
						continue
					#	raise AssertionError(f"New Tokken detected in reaction {i}\n")
			
			if i in self.rxnList:
				index = self.rxnList.index(i)
				#print(index)
				r_equation = self.mech["reactions"][index]
				if "duplicate" in r_equation:
					if r_equation["duplicate"] == True:
						
						if description == "A":
							PerturbingRxnDict[rxn] = r_equation
						else:
							PerturbingRxnDict[rxn] = self.mech["reactions"][index+1]
				else:
					PerturbingRxnDict[rxn] = r_equation
				#print(r_equation)
		return PerturbingRxnDict

	#def getPerturbationData(self,selection,beta):

#trial = Parser(mechanism)
#dictionary = trial.PerturbingReactions(rxn)
#print(dictionary["C4H71-3OOH => C2H3CHO + CH3 + OH:B1"])
