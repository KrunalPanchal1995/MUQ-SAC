import reaction_selection as rs
import numpy as np
import sys
import os
import pickle

def dictionary_creator(flag, mechanism, carbon_number, rxn_list, species_list):
    species = mechanism['phases'][0]["species"]  # Takes all the species available in the mechanism file
    species_data = mechanism["species"]
    reactions = mechanism["reactions"]
    
    if len(rxn_list) == 0 and len(species_list) == 0:
        if carbon_number != 0:
            # Get the species list greater than or equal to the predetermined carbon number
            selected_species = rs.species_selection(species, species_data, carbon_number)

            # Get the reaction list containing the selected species
            selected_reactions = rs.reaction_selection(selected_species, reactions)

            # Get the reaction index for the selected reactions
            reaction_dict = rs.reaction_index(selected_reactions, reactions)
        else:
            selected_species = species
            selected_reactions = []
            reaction_dict = {}
            for index, rxn in enumerate(reactions):
                selected_reactions.append(rxn["equation"])
                reaction_dict[index + 1] = rxn["equation"]
    else:
        selected_reactions = rxn_list  # For sensitivity analysis of reactions 
        selected_species = species_list  # For sensitivity analysis of thermo data
        reaction_dict = rs.reaction_index(selected_reactions, reactions)

    if flag == 'reaction':
        rxn_type = rs.getRxnType(mechanism, selected_reactions)
        string_f = ""
        string_g = ""
        for index in reaction_dict:
            string_f += f"{index}\t{reaction_dict[index]}\n"
        for rxn in rxn_type:
            string_g += f"{rxn}\t{rxn_type[rxn]}\n"
        f = open("Reaction_dict.txt", "w").write(string_f)
        g = open("Reaction_type.txt", "w").write(string_g)

        # Creating reaction dictionary
        rxn_dict = {}
        rxn_dict["reaction"] = reaction_dict
        rxn_dict["type"] = rxn_type
        rxn_dict["data"] = rs.getRxnDetails(mechanism, selected_reactions)

        if "RXN_DICT.pkl" not in os.listdir():
            with open('RXN_DICT.pkl', 'wb') as file_:
                pickle.dump(rxn_dict, file_)

        string_reaction = ""
        for index in reaction_dict:
            string_reaction += f"{index}\t{reaction_dict[index]}\n"
        g = open("selected_rxn.txt", "+w").write(string_reaction)
        return rxn_dict, selected_species, selected_reactions
    
    if flag == 'thermo':  # Creating species dictionary 
        species_dict = {}
        species_data = mechanism["species"] 
        reactions = mechanism["reactions"]

        # Loop through selected species to get their thermo data
        for species in species_data:
            species_name = species['name']  # Extract species name
            if species_name in selected_species:
                thermo_data = species['thermo']  # Extract thermo data
                species_dict[species_name] = thermo_data 

        """
        # Add reaction details only for selected reactions 
        reaction_details_list = []  # Create a list to store reaction details
        for rxn in reactions:
            if rxn["equation"] in selected_reactions:
                reaction_details = rs.getRxnDetails(mechanism, [rxn["equation"]])
                reaction_details_list.append({"reaction": rxn["equation"], "details": reaction_details})

        # Add the reaction details list to the species_dict
        species_dict["reactions"] = reaction_details_list 
        species_dict["reaction"] = reaction_details_list
		"""
		
        if "SPECIES_DICT.pkl" not in os.listdir():
            with open('SPECIES_DICT.pkl', 'wb') as file_:
                pickle.dump(species_dict, file_)
        
        return species_dict, selected_species, selected_reactions

