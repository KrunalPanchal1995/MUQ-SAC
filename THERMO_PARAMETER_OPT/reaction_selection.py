def species_selection(species,species_data,carbon_number):
	selected_species = []
	for kind in species:
		for data in species_data:
			if kind in data["name"]:
				kind_composition = data["composition"]
				if "C" in kind_composition:
					if kind_composition["C"] >= carbon_number:
						selected_species.append(data["name"])
			
	selected_species = set(selected_species)
	return selected_species
	
	
def reaction_selection(selected_species,reactions):
	selected_reactions = []
	for species in selected_species:
		for reaction in reactions:
			if species in reaction["equation"]:
				selected_reactions.append(reaction["equation"])

	selected_reactions = set(selected_reactions)
	return selected_reactions

def reaction_index(selected_reactions,reactions):
	reaction_dict = {}
	for rxn in selected_reactions:
		#print(rxn)
		for index,reaction in enumerate(reactions):
			if rxn.strip().split(":")[0] == reaction["equation"]:
				#print("\t"+reaction["equation"])
				reaction_dict[index+1] = rxn.split(":")[0]
				#break
	
	return reaction_dict

def getRxnType(mechanism,rxn_list):
	rxn_type = {}
	rxn_data = mechanism["reactions"]
	for rxn in rxn_list:
		for index,data in enumerate(rxn_data):
			if rxn in data["equation"]:
				if "type" in data:
					if data["type"] == "three-body":
						rxn_type[index+1] = "ThirdBody"
						#print("Rxn is three-body")
					elif data["type"] == "falloff":
						rxn_type[index+1] = "Falloff"
						#print("Rxn is fall-off")
					elif data["type"] == "Chebyshev":
						rxn_type[index+1] = "Chebyshev"
						
					elif data["type"] == "pressure-dependent-Arrhenius":# and "duplicate" not in data:
						rxn_type[index+1] = "PLOG"
						#print("Rxn is Plog")
					#elif data["type"] == "pressure-dependent-Arrhenius" and "duplicate" in data:
						#rxn_type[data["equation"]] = "PLOG-Duplicate"
						#break
						#print("Rxn is Plog duplicate")
				#elif "duplicate" in data:
				#	rxn_type[data["equation"]] = "Duplicate"
				#	break
					#print("Rxn is duplicate")
				
				else:
					rxn_type[index+1] = "Elementary"
					#print("Rxn is elementary")
	return rxn_type
	
def getRxnDetails(mechanism,rxn_list):
	rxn_dict = {}
	rxn_data = mechanism["reactions"]
	for rxn in rxn_list:
		new_rxn_data = {}
		temp = []
		index_ = []
		rxn_ = []
		for index,data in enumerate(rxn_data):
			if rxn == data["equation"]:
				temp.append(data)
				index_.append(index+1)
		new_rxn_data["temp"] = temp
		new_rxn_data["index"] = index_
		rxn_dict[rxn] = new_rxn_data
	return rxn_dict

