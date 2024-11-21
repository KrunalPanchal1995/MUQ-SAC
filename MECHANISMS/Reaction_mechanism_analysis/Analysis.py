from Input_file_reader import MechParsing,ThermoParsing

#mechFile = "/home/krithika/Desktop/KineticMechanismOptimization/Mechanism/FFCM1/FFCM_H2_CO.mech"
mechFile = "/home/krithika/Desktop/KineticMechanismOptimization/Mechanism/Keromens/Keromense.mech"
#mechFile = "/home/krithika/Desktop/KineticMechanismOptimization/Mechanism/Connaire/h2_v1b_mech.txt"
#mechFile = "/home/krithika/Desktop/KineticMechanismOptimization/Mechanism/CRECK/CRECK.mech"
#mechFile = "/home/krithika/Desktop/KineticMechanismOptimization/Mechanism/NUIG-NGM_2010/NUIG.mech"
#mechFile = "/home/krithika/Desktop/KineticMechanismOptimization/Mechanism/Li_2015/Li.mech"
#mechFile = "/home/krithika/Desktop/KineticMechanismOptimization/Mechanism/SanDeiago/SanDiago_H2_CO.mech"
#mechFile = "/home/krithika/Desktop/KineticMechanismOptimization/Mechanism/GRI3/GRI_CO.mech"

mechanism = MechParsing(mechFile)

rxn_index = mechanism.rxn_index

rxn_dict = {}
string="Keromens mechanism\n"
string+=f"{len(rxn_index)}\n"
for i in rxn_index:
	rxn_dict[i] = mechanism.getKappa(i)
	string+=f"{i},{mechanism.getKappa(i)}\n"
fi = open("Keromens","w").write(string)
